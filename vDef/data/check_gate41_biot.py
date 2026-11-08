#!/usr/bin/env python3
"""Gate 4.1 -- Biot body-force superposition verdict.

Pure standard library + the committed ``parse_cdf2.py`` (the CI runner has
no numpy / scipy / matplotlib).

Two modes
---------
  report : python3 check_gate41_biot.py <prod.log> <prod.gen> <nobiot.gen>
      INFORMATIONAL. Always exits 0, never a gate.

  gate   : python3 check_gate41_biot.py --gate <prod.log> <prod.gen> <nobiot.gen>
      HARD Gate 4.1. Prints ``GATE41_RESULT: PASS|FAIL`` and exits 0 / 1.
      Missing / empty / short / unreadable / non-finite data is ALWAYS FAIL
      (CI #73 lesson: "no data" must never read as a pass).

What checklist item 4.1 actually asks
-------------------------------------
"开启体力缩放后结果不回滚" -- the preset body force is rescaled on EVERY
outer step by ``MEF90DefMechSetTransients`` (default ForceScaling=Linear,
Force = (0,-9.81,0) in terzaghi1D.yaml), so |baseForce| must grow
0 -> 9.81 over the 51 steps. An implementation that snapshots baseForce
only once at startup, or copies the snapshot back over the scaled force,
pins it at its t=0 value -> B3 fails.

The force written into DefMechCtx must be ``baseForce - alpha*grad(p)``.
``superpos_err`` is ||force - (base - alpha*gradp)||_inf computed from the
UNMODIFIED snapshot, so it is ~0 for a true superposition, ~|base| if the
force was overwritten by the Biot term alone, and ~|alpha*gradp| if the
Biot term was dropped. This single number distinguishes all three.

The NOBIOT control deck (terzaghi1D_nobiot.yaml, only biot_coeff 1.0->0.0)
gives an independent physical check that the coupling is live: with alpha=0
the u_x peak collapses to the gravity-only level.

Criteria
--------
  B1  >= 50 parseable ``Gate 4.1:`` rows, and row count ~= Gate 3.2 probes
  B2  alpha == 1.0 on every row
  B3  |base|inf grows linearly 0 -> >=5 (Linear scaling NOT rolled back)
  B4  superpos_err ~= 0 on every row (superposition, not overwrite/accumulate)
  B5  |gradp|inf > 0, evolves in time, and |force| exceeds |base|
  B6  peak|u_x|(main) >= 3 x peak|u_x|(nobiot), and peak in [5e-5,2e-4] m
  B7  production log has >= 50 Gate 3.2 probes, no crash / unknown option

Note on checklist 4.1b (effective stress): adding a separate
(-alpha*p*I) volume term ON TOP of the body force would DOUBLE-COUNT the
pore pressure. Integration by parts gives
    int(-alpha*grad p).v = int(alpha*p*div v) - int_boundary(alpha*p*v.n)
so the two forms share the same volume term and differ only by the natural
boundary traction. This gate verifies the existing (body-force) form; a
mutually exclusive stress form is tracked separately as 4.1b.
"""
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse_cdf2 import CDF2  # noqa: E402

# Fortran list-directed Write(IOBuffer,*) output, e.g.
# " Gate 4.1: alpha=  1.0000E+00  |base|inf=  9.8100E+00 ... superpos_err= 0.0"
ROW_RE = re.compile(
    r'Gate\s+4\.1:\s*alpha=\s*([-+0-9.eEdD]+)\s*'
    r'\|base\|inf=\s*([-+0-9.eEdD]+)\s*'
    r'\|gradp\|inf=\s*([-+0-9.eEdD]+)\s*'
    r'\|force\|inf=\s*([-+0-9.eEdD]+)\s*'
    r'superpos_err=\s*([-+0-9.eEdD]+)')
CRASH_RE = re.compile(r'Caught signal|SEGV|MPI_ABORT|PETSC ERROR', re.I)

MIN_ROWS = 50
MIN_PROBES = 50
TOL_ALPHA = 1.0e-6
BASE_MIN_END = 5.0          # 9.81 at t=1
BASE_MAX_START = 1.0        # 0 at t=0
R2_MIN = 0.98
SLOPE_MIN = 5.0
LIN_RESID_MAX = 0.5
TOL_SUPERPOS = 1.0e-8       # relative to max(1,|force|)
GRAD_EVOL_MIN = 1.5         # max/min of |gradp|
FORCE_EXCESS_MIN = 1.0      # |force|-|base| must exceed this somewhere
UX_RATIO_MIN = 3.0
UX_PEAK_LO, UX_PEAK_HI = 5.0e-5, 2.0e-4


def _f(s):
    """Fortran float (D exponent) -> python float; NaN/Inf stay non-finite."""
    return float(s.replace('D', 'E').replace('d', 'e'))


def read_rows(path):
    if not path or not os.path.isfile(path) or os.path.getsize(path) == 0:
        return None, 'missing/empty log: %s' % path
    rows, probes, bad = [], 0, []
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as fh:
            for line in fh:
                if 'Gate 3.2: step' in line:
                    probes += 1
                m = ROW_RE.search(line)
                if m:
                    rows.append(tuple(_f(v) for v in m.groups()))
                elif 'Gate 4.1' in line:
                    bad.append(line.strip()[:120])
    except Exception as exc:                       # noqa: BLE001 - gate fails
        return None, 'log unreadable: %s' % exc
    return dict(rows=rows, probes=probes, unparsed=bad), None


def read_log_text(path):
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as fh:
            return fh.read()
    except Exception:                              # noqa: BLE001 - gate fails
        return ''


def peak_ux(path):
    """max_t max_node |u_x| from Exodus slot 3; returns (peak, nrec)."""
    if not path or not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise ValueError('missing/empty mesh: %s' % path)
    c = CDF2(path)
    if c.numrecs < MIN_ROWS:
        raise ValueError('only %d records (need >= %d): %s'
                         % (c.numrecs, MIN_ROWS, path))
    peak = 0.0
    for i in range(c.numrecs):
        vals = [float(v) for v in c.read_record('vals_nod_var3', i)]
        if not vals:
            raise ValueError('empty u_x field at record %d: %s' % (i + 1, path))
        peak = max(peak, max(abs(v) for v in vals))
    if not math.isfinite(peak):
        raise ValueError('non-finite u_x peak: %s' % path)
    return peak, c.numrecs


def linfit(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    k = sxy / sxx if sxx > 0 else 0.0
    b = my - k * mx
    r2 = (sxy * sxy) / (sxx * syy) if sxx > 0 and syy > 0 else 0.0
    resid = max(abs(y - (k * x + b)) for x, y in zip(xs, ys))
    return k, b, r2, resid


def evaluate(log_path, main_mesh, nobiot_mesh):
    """Return (checks, metrics, load_error). checks: name -> (ok, detail)."""
    chk = {}
    met = {}
    data, err = read_rows(log_path)
    if err:
        return chk, met, err
    rows = data['rows']
    met['rows'] = len(rows)
    met['probes'] = data['probes']
    met['unparsed'] = data['unparsed']

    # B1 -- completeness
    ok = (len(rows) >= MIN_ROWS and not data['unparsed']
          and abs(len(rows) - data['probes']) <= 2
          and all(all(math.isfinite(v) for v in r) for r in rows))
    chk['B1'] = (ok, '%d parseable rows, %d Gate 3.2 probes, %d unparsed'
                 % (len(rows), data['probes'], len(data['unparsed'])))
    if not ok:
        return chk, met, None

    alpha = [r[0] for r in rows]
    base = [r[1] for r in rows]
    grad = [r[2] for r in rows]
    force = [r[3] for r in rows]
    serr = [r[4] for r in rows]
    met.update(alpha=alpha, base=base, grad=grad, force=force, serr=serr)

    # B2 -- alpha
    worst_a = max(abs(a - 1.0) for a in alpha)
    met['worst_alpha'] = worst_a
    chk['B2'] = (worst_a <= TOL_ALPHA,
                 'max |alpha-1| = %.3g (tol %.3g)' % (worst_a, TOL_ALPHA))

    # B3 -- Linear force scaling is NOT rolled back
    ts = [i / float(len(base) - 1) for i in range(len(base))]
    k, b, r2, resid = linfit(ts, base)
    met.update(lin_k=k, lin_b=b, lin_r2=r2, lin_resid=resid)
    ok = (max(base) >= BASE_MIN_END and min(base) <= BASE_MAX_START
          and r2 >= R2_MIN and k >= SLOPE_MIN and resid <= LIN_RESID_MAX)
    chk['B3'] = (ok, '|base| %.4g -> %.4g, slope %.4g, R2 %.5f, '
                 'max resid %.4g' % (min(base), max(base), k, r2, resid))

    # B4 -- true superposition (no overwrite, no accumulation)
    worst_s = max(s / max(1.0, abs(f)) for s, f in zip(serr, force))
    met['worst_superpos_rel'] = worst_s
    chk['B4'] = (worst_s <= TOL_SUPERPOS,
                 'max superpos_err/|force| = %.3g (tol %.3g)'
                 % (worst_s, TOL_SUPERPOS))

    # B5 -- Biot term is live and time-dependent
    gmin, gmax = min(grad), max(grad)
    evol = gmax / gmin if gmin > 0 else float('inf')
    excess = max(f - bb for f, bb in zip(force, base))
    met.update(grad_min=gmin, grad_max=gmax, grad_evol=evol,
               force_excess=excess)
    ok = (gmin > 0.0 and evol >= GRAD_EVOL_MIN
          and excess >= FORCE_EXCESS_MIN and max(force) > max(base))
    chk['B5'] = (ok, '|gradp| %.4g..%.4g (ratio %.3g), '
                 'max(|force|-|base|) %.4g' % (gmin, gmax, evol, excess))

    # B6 -- NOBIOT control: the coupling really drives the displacement
    try:
        pm, nm = peak_ux(main_mesh)
        pn, nn = peak_ux(nobiot_mesh)
        ratio = pm / pn if pn > 0 else float('inf')
        met.update(peak_main=pm, peak_nobiot=pn, ux_ratio=ratio,
                   nrec_main=nm, nrec_nobiot=nn)
        ok = (ratio >= UX_RATIO_MIN and UX_PEAK_LO <= pm <= UX_PEAK_HI)
        chk['B6'] = (ok, 'peak|u_x| main %.5g m (%d rec) vs nobiot %.5g m '
                     '(%d rec), ratio %.4g' % (pm, nm, pn, nn, ratio))
    except Exception as exc:                       # noqa: BLE001 - gate fails
        chk['B6'] = (False, 'mesh read FAILED: %s' % exc)

    # B7 -- production run completeness
    txt = read_log_text(log_path)
    crash = bool(CRASH_RE.search(txt))
    unknown = 'Unknown option' in txt
    ok = (data['probes'] >= MIN_PROBES and not crash and not unknown)
    chk['B7'] = (ok, '%d probes, crash=%s, unknown_option=%s'
                 % (data['probes'], crash, unknown))
    return chk, met, None


def print_report(chk, met, err):
    print('=== Gate 4.1 -- Biot body-force superposition verdict ===')
    if err:
        print('  DATA ERROR: %s' % err)
    for name in ('B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7'):
        if name in chk:
            ok, detail = chk[name]
            print('  %-3s %-4s %s' % (name, 'PASS' if ok else 'FAIL', detail))
        else:
            print('  %-3s ---- (not evaluated)' % name)
    if met.get('rows'):
        print('  --- time trace (first / middle / last) ---')
        n = met['rows']
        for i in (0, n // 2, n - 1):
            print('    step %2d  |base|=%10.4g  |gradp|=%10.4g  '
                  '|force|=%10.4g  err=%8.3g'
                  % (i + 1, met['base'][i], met['grad'][i],
                     met['force'][i], met['serr'][i]))
    fails = [k for k, v in chk.items() if not v[0]]
    return fails


def main(argv):
    gate = bool(argv) and argv[0] == '--gate'
    args = argv[1:] if gate else argv
    if gate and len(args) != 3:
        print('GATE41_RESULT: FAIL (usage: --gate <prod.log> <prod.gen> <nobiot.gen>)')
        return 1
    if not gate and len(args) != 3:
        print('usage: check_gate41_biot.py [--gate] <prod.log> <prod.gen> <nobiot.gen>')
        return 0
    chk, met, err = evaluate(args[0], args[1], args[2])
    fails = print_report(chk, met, err)
    if gate:
        if err or fails or len(chk) < 7:
            print('GATE41_RESULT: FAIL (%s)'
                  % (err if err else ','.join(fails) or 'incomplete criteria'))
            return 1
        print('GATE41_RESULT: PASS (B1-B7 on %d steps, '
              'Linear force scaling preserved, superposition error ~0, '
              'NOBIOT ratio %.4g)'
              % (met['rows'], met.get('ux_ratio', float('nan'))))
        return 0
    print('VERDICT: %s (report mode, informational only)'
          % ('FAIL ' + ','.join(fails) if fails else 'PASS'))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
