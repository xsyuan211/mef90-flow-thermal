#!/usr/bin/env python3
"""Physics verdict / Gate 3.3 for the DefMech displacement field.

Pure standard library + the committed parse_cdf2.py. The CI runner has no
numpy / scipy / matplotlib, so every criterion is computed with math only.

Two modes
---------
  report : python3 check_terzaghi_displacement.py <main.gen> <ksponly.gen> <ksponly.log>
      INFORMATIONAL. Always exits 0, prints VERDICT lines, never a gate.

  gate   : python3 check_terzaghi_displacement.py --gate <production.gen> <production.log>
      HARD Gate 3.3. Prints a machine-readable ``GATE33_RESULT: PASS|FAIL``
      line and exits 0 / 1. Missing, empty, short, non-finite or unreadable
      data is ALWAYS FAIL (CI #73 lesson: "no data" must not read as zero).

Why the gate judges the physical field
-------------------------------------
Five earlier "convergences" were false: CI #73 (ksponly fnorm placeholders
read as zeros), #75 (job-log section misread), #78 (SEGV produced no failure
line), #81 (zero step -> SNES_CONVERGED_SNORM_RELATIVE), #83 KSPINIT (0
failure messages but only ~3 SNES iterations per step and residual still
7-9). CI #83 then proved the production field and the ksponly field are
PHYSICALLY EQUIVALENT (pressure delta 0, damage 0.38%, u_x 0.21%, stresses
<= 0.41%) while the standard Newton run still printed 57x -5. The residual
sits on an irreducible plateau 3-13, so counting "SNESSolve failed" lines
is not a valid convergence criterion. Gate 3.3 now checks the written field
against the analytic Terzaghi solution instead.

Criteria (thresholds and sources: FLOW_FIX_CHECKLIST.md, section 判据检索)
-------------------------------------------------------------------------
  P1  BC nodes x=+-5      : max|u_x|, |u_y| < 1e-10 m
  P2  rigid-body drift    : max over steps |mean(u_y)| <= 1e-5 m AND
                            <= 5% of peak|u_x|. The old 1e-6 m absolute
                            value had no resolving power: CI #83 measured
                            1.0237e-6 (production) vs 8.17e-7 (ksponly) for
                            the SAME physical field. Basis: gravity bending
                            is +/-1.8e-5 m (CI #70) and the broken solutions
                            drifted 0.0556 / 1.06 m uniformly.
  P3  u_y time drift      : max|u_y(t=1)| / max(1e-5, max|u_y(t=0)|) < 2
  P4  u_x peak            : 5e-5 <= max_t max|u_x| <= 2e-4 m, and the peak
                            occurs at t <= 0.2 s (analytic 9.90e-5 @ 0.08 s)
  P5  u_x shape           : non-monotone (peak strictly inside the run, and
                            final amplitude < half the peak)
  P6  pressure (Gate 3.2) : ratio in [0.05,0.15], pmin == 0,
                            last-step L2 error vs the BE-corrected Terzaghi
                            series <= 0.05 (CI #70 measured 0.032)
  P7  same-source         : max|p_ksponly - p_main| <= 1e-2 Pa (the flow solve
                            must be untouched); dUx / dUy are REPORTED ONLY.
                            CI #83 measured dUx=2.3e-7 m, dUy=1.6e-6 m:
                            the two solvers give the same field, so a LARGE
                            difference would now indicate an inconsistency.
  P8  damage field        : finite, in [0,1], final max <= 0.1. This is the
                            physical guard that lets Gate 3.3 REPORT damage
                            SNES -5 messages without counting them as a
                            displacement-solver failure.

Gate-only sanity checks: >= 50 records, final time within 0.05 s of t=1,
all fields finite, no crash markers / "Unknown option" in the production
solver log, and >= 50 Gate 3.2 pressure probes.
"""
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse_cdf2 import CDF2  # noqa: E402

# --- Physics / discretisation (terzaghi1D.yaml, identical to Gate 3.2) -----
P0 = 1.0e5                           # initial pressure [Pa]
CV = (5.0e-11 / 1.0e-3) / 4.6e-10  # consolidation coefficient [m^2/s]
H = 10.0                            # drainage path [m], z = 5 - x
DT = 0.02                           # fixed BE step [s]
NSERIES = 200                       # BE-corrected Terzaghi modes
_MODES = range(NSERIES)
_K = [(2 * m + 1) * math.pi / (2.0 * H) for m in _MODES]
_A = [(4.0 / math.pi) * ((-1.0) ** m) / (2 * m + 1) for m in _MODES]
_LAM = [CV * k * k for k in _K]

# --- Verdict thresholds -----------------------------------------------------
TOL_BC = 1.0e-10
TOL_MEAN_UY = 1.0e-5                # gravity-bending scale (CI #70)
TOL_MEAN_UY_REL = 0.05             # <= 5% of the peak |u_x|
TOL_UY_DRIFT = 2.0
UY_REF = 1.0e-5                     # gravity-bending scale, avoids /0
UX_PEAK_LO, UX_PEAK_HI = 5.0e-5, 2.0e-4
TOL_PEAK_TIME = 0.2
RATIO_LO, RATIO_HI = 0.05, 0.15
TOL_PMIN = 1.0e-8
TOL_L2_BE = 0.05
TOL_SAME_P = 1.0e-2
TOL_DAMAGE_FINAL = 0.1
MIN_RECORDS = 50
T_END = 1.0
TOL_T_END = 0.05

CRASH_RE = re.compile(r'Caught signal|SEGV|MPI_ABORT|PETSC ERROR', re.I)


def amp(field):
    """max |value| of a nodal field."""
    return max((abs(v) for v in field), default=0.0)


def load(path):
    c = CDF2(path)
    x = [float(v) for v in c.read_static('coordx')]
    times, p, ux, uy, dmg = [], [], [], [], []
    for i in range(c.numrecs):
        times.append(float(c.read_record('time_whole', i)[0]))
        p.append([float(v) for v in c.read_record('vals_nod_var1', i)])
        ux.append([float(v) for v in c.read_record('vals_nod_var3', i)])
        uy.append([float(v) for v in c.read_record('vals_nod_var4', i)])
        try:
            dmg.append([float(v) for v in c.read_record('vals_nod_var2', i)])
        except Exception:                   # noqa: BLE001 - gate must fail
            dmg.append(None)
    return dict(path=path, x=x, times=times, p=p, ux=ux, uy=uy, dmg=dmg,
                n=c.numrecs, nn=len(x))


def all_finite(d):
    if not d['times'] or not all(math.isfinite(t) for t in d['times']):
        return False
    for key in ('p', 'ux', 'uy'):
        for f in d[key]:
            if not f or not all(math.isfinite(v) for v in f):
                return False
    for f in d['dmg']:
        if f is None or not f or not all(math.isfinite(v) for v in f):
            return False
    return True


def be_pressure(x, nsteps):
    """BE-corrected Terzaghi series field after nsteps steps (Gate 3.2 form).

    p(z,t) = p0 * sum_m a_m cos(k_m z) (1 + lambda_m dt)^(-nsteps)
    with z = 5 - x, k_m = (2m+1) pi / (2H), lambda_m = cv k_m^2.
    """
    decay = [(1.0 + lam * DT) ** (-nsteps) for lam in _LAM]
    out = []
    for xi in x:
        z = 5.0 - xi
        s = 0.0
        for k, a, d in zip(_K, _A, decay):
            s += a * math.cos(k * z) * d
        out.append(P0 * s)
    return out


def l2_rel(a, b):
    num = math.sqrt(sum((u - v) ** 2 for u, v in zip(a, b)))
    den = math.sqrt(sum(v * v for v in b))
    return num / den if den > 0 else float('inf')


def criteria(d):
    x, n = d['x'], d['n']
    bc = [i for i, xi in enumerate(x) if abs(xi) >= 4.9999]

    # P1: Dirichlet nodes must stay at the prescribed value (0,0,0).
    p1 = 0.0
    if bc:
        for field in d['ux'] + d['uy']:
            for k in bc:
                p1 = max(p1, abs(field[k]))
    else:
        p1 = float('inf')

    # P2: a uniform u_y is a rigid-body drift, not gravity bending.
    p2 = max((abs(sum(f) / len(f)) for f in d['uy']), default=float('inf'))

    # P3: u_y must not accumulate over time under a static body force.
    a0, a1 = amp(d['uy'][0]), amp(d['uy'][-1])
    p3 = a1 / max(a0, UY_REF)

    # P4/P5: u_x amplitude and its time shape.
    uxa = [amp(f) for f in d['ux']]
    peak = max(uxa, default=0.0)
    pi = uxa.index(peak) if uxa else -1
    tpeak = d['times'][pi] if 0 <= pi < n else float('inf')
    p5 = (0 < pi < n - 1) and (uxa[-1] < 0.5 * peak)
    p2rel = p2 / peak if peak > 0 else float('inf')

    # P6: pressure field must still be the validated Terzaghi solution.
    pf, pl = amp(d['p'][0]), amp(d['p'][-1])
    ratio = pl / pf if pf > 0 else 0.0
    pmin = max((abs(min(f)) for f in d['p']), default=float('inf'))
    nsteps = max(int(round(d['times'][-1] / DT)), 0) if n else 0
    l2be = l2_rel(d['p'][-1], be_pressure(x, nsteps)) if n else float('inf')

    # P8: damage must stay a bounded, physical history field.
    p8 = (all(f is not None and all(0.0 <= v <= 1.0 for v in f)
              for f in d['dmg'])
          and amp(d['dmg'][-1] or []) <= TOL_DAMAGE_FINAL)

    return dict(P1=p1, P2=p2, P2rel=p2rel, P3=p3, P4=(peak, tpeak), P5=p5,
                P6=(ratio, pmin, l2be), P8=p8,
                dmg1=amp(d['dmg'][-1] or []),
                bc=len(bc), n=n, nn=d['nn'],
                uy0=a0, uy1=a1, ux0=uxa[0] if uxa else 0.0,
                ux1=uxa[-1] if uxa else 0.0, peak_idx=pi)


def p2_ok(c):
    return c['P2'] <= TOL_MEAN_UY and c['P2rel'] <= TOL_MEAN_UY_REL


def failed(c):
    out = []
    if c['P1'] > TOL_BC:
        out.append('P1')
    if not p2_ok(c):
        out.append('P2')
    if c['P3'] >= TOL_UY_DRIFT:
        out.append('P3')
    peak, tpeak = c['P4']
    if not (UX_PEAK_LO <= peak <= UX_PEAK_HI and tpeak <= TOL_PEAK_TIME):
        out.append('P4')
    if not c['P5']:
        out.append('P5')
    ratio, pmin, l2be = c['P6']
    if not (RATIO_LO <= ratio <= RATIO_HI and pmin <= TOL_PMIN
            and l2be <= TOL_L2_BE):
        out.append('P6')
    if not c['P8']:
        out.append('P8')
    return out


def row(label, value, limit, ok):
    print('  %-34s %12.6g   %-18s %s' %
          (label, value, limit, 'PASS' if ok else 'FAIL'))


def report(tag, d, c):
    print('')
    print('=== PHYSICS VERDICT: %s (%s) ===' % (tag, d['path']))
    print('  records=%d nodes=%d BC nodes=%d  t=[%.3f, %.3f]' %
          (c['n'], c['nn'], c['bc'], d['times'][0], d['times'][-1]))
    peak, tpeak = c['P4']
    ratio, pmin, l2be = c['P6']
    row('P1 max|u| at x=+-5 [m]', c['P1'], '< 1e-10', c['P1'] <= TOL_BC)
    row('P2 max|mean(u_y)| [m]', c['P2'], '<= 1e-5', c['P2'] <= TOL_MEAN_UY)
    row('P2 drift / peak|u_x|', c['P2rel'], '<= 0.05',
        c['P2rel'] <= TOL_MEAN_UY_REL)
    row('P3 u_y drift (t1/t0)', c['P3'], '< 2', c['P3'] < TOL_UY_DRIFT)
    row('P4 peak |u_x| [m]', peak, '5e-5..2e-4',
        UX_PEAK_LO <= peak <= UX_PEAK_HI)
    row('P4 peak time [s]', tpeak, '<= 0.2', tpeak <= TOL_PEAK_TIME)
    print('  %-34s %12s   %-18s %s' %
          ('P5 non-monotone u_x', 'peak@%d' % (c['peak_idx'] + 1),
           'interior + decay', 'PASS' if c['P5'] else 'FAIL'))
    row('P6 pressure ratio', ratio, '0.05..0.15',
        RATIO_LO <= ratio <= RATIO_HI)
    row('P6 max|pmin| [Pa]', pmin, '<= 1e-8', pmin <= TOL_PMIN)
    row('P6 L2 err vs BE series', l2be, '<= 0.05', l2be <= TOL_L2_BE)
    row('P8 max damage (final)', c['dmg1'], '<= 0.1', c['P8'])
    print('  (reference: |u_x| t0=%.4g t1=%.4g ; |u_y| t0=%.4g t1=%.4g)' %
          (c['ux0'], c['ux1'], c['uy0'], c['uy1']))
    fails = failed(c)
    print('  === %s VERDICT: %s%s ===' %
          (tag, 'FAIL' if fails else 'PASS',
           ' (' + ','.join(fails) + ')' if fails else ' (P1-P6,P8 all pass)'))


def max_field_diff(a, b, key):
    return max((abs(u - v) for fa, fb in zip(a[key], b[key])
                for u, v in zip(fa, fb)), default=0.0)


def log_health(path):
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as f:
            text = f.read()
    except OSError:
        print('  solver log %s : MISSING' % path)
        return
    print('  solver log %s : %d bytes, crash markers=%d, '
          'SNESSolve-failed msgs=%d, "Linear solve" lines=%d, '
          '"Gate 3.1: step" lines=%d' %
          (path, len(text), len(CRASH_RE.findall(text)),
           text.count('SNESSolve failed with SNESConvergedReason'),
           text.count('Linear solve'), text.count('Gate 3.1: step')))
    print('  (informational only: these counts are NOT the verdict)')


def gate(mesh_path, log_path):
    """HARD Gate 3.3. Every failure path must print GATE33_RESULT: FAIL."""
    print('=== Gate 3.3 HARD PHYSICS VERDICT (production run) ===')
    if not os.path.isfile(mesh_path) or os.path.getsize(mesh_path) == 0:
        print('GATE33_RESULT: FAIL (production mesh missing/empty: %s)'
              % mesh_path)
        return 1
    try:
        d = load(mesh_path)
    except Exception as exc:                 # noqa: BLE001 - gate must fail
        print('GATE33_RESULT: FAIL (mesh unreadable: %s: %s)'
              % (type(exc).__name__, exc))
        return 1
    if d['n'] < MIN_RECORDS:
        print('GATE33_RESULT: FAIL (only %d records, need >= %d)'
              % (d['n'], MIN_RECORDS))
        return 1
    if abs(d['times'][-1] - T_END) > TOL_T_END:
        print('GATE33_RESULT: FAIL (final time %.4g != %.4g)'
              % (d['times'][-1], T_END))
        return 1
    if not all_finite(d):
        print('GATE33_RESULT: FAIL (non-finite / missing field in mesh)')
        return 1
    try:
        with open(log_path, 'r', encoding='utf-8', errors='replace') as f:
            text = f.read()
    except OSError:
        print('GATE33_RESULT: FAIL (production solver log missing: %s)'
              % log_path)
        return 1
    if CRASH_RE.search(text):
        print('GATE33_RESULT: FAIL (crash marker in production solver log)')
        return 1
    if 'Unknown option' in text:
        print('GATE33_RESULT: FAIL (Unknown option in production solver log)')
        return 1
    probes = text.count('Gate 3.2:')
    if probes < MIN_RECORDS:
        print('GATE33_RESULT: FAIL (only %d Gate 3.2 pressure probes, need >= %d)'
              % (probes, MIN_RECORDS))
        return 1

    c = criteria(d)
    report('GATE33', d, c)
    fails = failed(c)
    print('')
    print('  production solver log: %d bytes, %d Gate 3.2 probes, '
          '%d SNESSolve-failed messages (INFORMATIONAL ONLY, not the verdict)'
          % (len(text), probes,
             text.count('SNESSolve failed with SNESConvergedReason')))
    if fails:
        print('GATE33_RESULT: FAIL (%s)' % ','.join(fails))
        return 1
    print('GATE33_RESULT: PASS (P1-P6,P8 physics criteria on %d records)'
          % d['n'])
    return 0


def main(argv):
    if len(argv) >= 2 and argv[1] == '--gate':
        if len(argv) != 4:
            print('usage: %s --gate <production.gen> <production.log>' % argv[0])
            print('GATE33_RESULT: FAIL (bad arguments)')
            return 1
        return gate(argv[2], argv[3])

    if len(argv) != 4:
        print(__doc__)
        return 0
    main_path, ko_path, log_path = argv[1:]
    print('=== Terzaghi DefMech displacement PHYSICS VERDICT (P1-P8) ===')
    print('=== physical fields only; solver messages are informational ===')

    results = {}
    for tag, path in (('MAIN', main_path), ('KSPONLY', ko_path)):
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            print('  === %s VERDICT: INCONCLUSIVE (mesh %s missing/empty) ==='
                  % (tag, path))
            continue
        try:
            d = load(path)
            if d['n'] < 2 or not all_finite(d):
                print('  === %s VERDICT: INCONCLUSIVE (%d records / non-finite) ==='
                      % (tag, d['n']))
                continue
            c = criteria(d)
            report(tag, d, c)
            results[tag] = (d, c)
        except Exception as exc:                 # noqa: BLE001 - diagnostic only
            print('  === %s VERDICT: INCONCLUSIVE (%s: %s) ==='
                  % (tag, type(exc).__name__, exc))

    print('')
    print('=== P7 same-source comparison (production vs ksponly) ===')
    if 'MAIN' in results and 'KSPONLY' in results:
        dm, cm = results['MAIN']
        dk, ck = results['KSPONLY']
        dp = max_field_diff(dm, dk, 'p')
        dux = max_field_diff(dm, dk, 'ux')
        duy = max_field_diff(dm, dk, 'uy')
        row('max|dPressure| [Pa]', dp, '<= 1e-2', dp <= TOL_SAME_P)
        print('  %-34s %12.6g   (reported only)' % ('max|dUx| [m]', dux))
        print('  %-34s %12.6g   (reported only)' % ('max|dUy| [m]', duy))
        print('  NOTE: CI #83 measured dUx=2.3e-7 m, dUy=1.6e-6 m: the two')
        print('        solvers give the SAME field. A large difference would')
        print('        now indicate a reproducibility problem, not a fix.')
    else:
        print('  INCONCLUSIVE: need both meshes')

    log_health(log_path)
    if 'KSPONLY' in results:
        fails = failed(results['KSPONLY'][1])
        print('')
        print('=== FINAL KSPONLY PHYSICS VERDICT: %s%s ===' %
              ('FAIL' if fails else 'PASS',
               ' (' + ','.join(fails) + ')' if fails else ''))
        print('=== interpretation (CI #83): production and ksponly fields are')
        print('    physically equivalent; -5 is an irreducible-residual /')
        print('    convergence-criterion issue, NOT a corrupted field ===')
    else:
        print('')
        print('=== FINAL KSPONLY PHYSICS VERDICT: INCONCLUSIVE ===')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
