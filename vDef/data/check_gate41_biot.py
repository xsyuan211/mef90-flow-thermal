#!/usr/bin/env python3
"""Gate 4.1 -- Biot body-force superposition verdict.

Pure standard library + the committed ``parse_cdf2.py`` (the CI runner has
no numpy / scipy / matplotlib).

Two modes
---------
  report : python3 check_gate41_biot.py <prod.log> <prod.gen> <nobiot.gen> <bodyforce.log>
      INFORMATIONAL. Always exits 0, never a gate.

  gate   : python3 check_gate41_biot.py --gate <prod.log> <prod.gen> <nobiot.gen> <bodyforce.log>
      HARD Gate 4.1. Prints ``GATE41_RESULT: PASS|FAIL`` and exits 0 / 1.
      Missing / empty / short / unreadable / non-finite data is ALWAYS FAIL
      (CI #73 lesson: "no data" must never read as a pass).

Why a BODYFORCE variant is required (CI #85)
--------------------------------------------
The production deck is body-force free. Its old ``Force: 0, -9.81, 0`` had
spaces after the commas; PETSc 3.3's ``PetscBagRegisterRealArray`` keeps only
the first element, so the CellSet bag view reads ``Force = 0 0 0``. With
``force_scaling = linear`` the preset force is ``linear * 0 = 0`` and
``|baseForce|`` is identically zero. On such a run:

  * B3 ("Linear scaling not rolled back") cannot pass -- there is nothing
    to roll back, so CI #85 failed B3 even though the Fortran code was fine;
  * B4 ("superposition vs overwrite") has NO discriminating power, because
    ``force = base - alpha*grad(p)`` and ``force = -alpha*grad(p)`` are the
    same when ``base == 0``. CI #85's ``err = 0`` proved nothing.

``terzaghi1D_bodyforce.yaml`` is byte-identical to the production deck except
``cs0001.Force: 0,0,0 -> 0,-9.81,0`` (no spaces after commas). It is a
DIAGNOSTIC run only; Gates 3.1/3.2/3.3 keep using the frozen production
mesh. On that run ``|baseForce|`` grows 0 -> 9.81 and B3/B4 become meaningful.

The force written into DefMechCtx must be ``baseForce - alpha*grad(p)``.
``superpos_err`` is ||force - (base - alpha*gradp)||_inf computed from the
UNMODIFIED snapshot, so it is ~0 for a true superposition, ~|base| if the
force was overwritten by the Biot term alone, and ~|alpha*gradp| if the
Biot term was dropped. This distinguishes the four implementations ONLY
when |base| > 0, which B4 now checks explicitly.

The NOBIOT control deck (only biot_coeff 1.0 -> 0.0) gives an independent
physical check: with alpha=0 and no body force, the displacement field is
exactly zero (CI #85 measured ux/uy == 0 over all 51 records).

Criteria
--------
  B1  >= 50 parseable rows in BOTH logs, row count ~= Gate 3.2 probes,
      all values finite, no unparsed rows
  B2  alpha == 1.0 on every row of BOTH logs
  B3  BODYFORCE |base|inf grows linearly 0 -> >=5 (Linear scaling not
      rolled back); production |base| is expected to be 0 and NOT judged
  B4  BODYFORCE superpos_err ~= 0 AND max|base| >= 5 (so the test can
      actually tell superposition from overwrite)
  B5  production |gradp|inf > 0, evolves in time, |force| exceeds |base|
  B6  NOBIOT control: main peak|u_x| in [5e-5,2e-4] m and either
      (a) nobiot peak == 0 (exact decoupling), or
      (b) main/nobiot >= 3. CI #85's ``ratio = inf`` is no longer used.
  B7  BOTH logs have >= 50 Gate 3.2 probes, no crash / unknown option

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
BASE_MIN_END = 5.0          # 9.81 at t=1 in the BODYFORCE variant
BASE_MAX_START = 1.0        # 0 at t=0
R2_MIN = 0.98
SLOPE_MIN = 5.0
LIN_RESID_MAX = 0.5
TOL_SUPERPOS = 1.0e-8       # relative to max(1,|force|)
GRAD_EVOL_MIN = 1.5         # max/min of |gradp|
FORCE_EXCESS_MIN = 1.0      # |force|-|base| must exceed this somewhere
UX_RATIO_MIN = 3.0
UX_ZERO_TOL = 1.0e-9       # NOBIOT field counts as identically zero below this
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
                # Match the COLON: 'Gate 4.1' is also a prefix of the Gate
                # 4.12/4.13 diagnostics, and of 'Gate 4.1b:'. Counting those
                # as unparsed Gate 4.1 rows made B1 fail on CI #116 (102
                # phantom rows = 2 logs x 51 steps).
                elif 'Gate 4.1:' in line:
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


def _unpack(data):
    rows = data['rows']
    return dict(alpha=[r[0] for r in rows], base=[r[1] for r in rows],
                grad=[r[2] for r in rows], force=[r[3] for r in rows],
                serr=[r[4] for r in rows])


def evaluate(prod_log, main_mesh, nobiot_mesh, bf_log):
    """Return (checks, metrics, load_error). checks: name -> (ok, detail)."""
    chk, met = {}, {}
    dp, errp = read_rows(prod_log)
    db, errb = read_rows(bf_log)
    if errp or errb:
        return chk, met, '; '.join(x for x in (errp, errb) if x)

    p, b = _unpack(dp), _unpack(db)
    met.update(p_rows=len(dp['rows']), b_rows=len(db['rows']),
               p_probes=dp['probes'], b_probes=db['probes'],
               p_unparsed=dp['unparsed'], b_unparsed=db['unparsed'],
               p=p, b=b)

    # B1 -- completeness of BOTH logs
    ok = all((
        len(dp['rows']) >= MIN_ROWS, len(db['rows']) >= MIN_ROWS,
        not dp['unparsed'], not db['unparsed'],
        abs(len(dp['rows']) - dp['probes']) <= 2,
        abs(len(db['rows']) - db['probes']) <= 2,
        all(all(math.isfinite(v) for v in r) for r in dp['rows']),
        all(all(math.isfinite(v) for v in r) for r in db['rows']),
    ))
    chk['B1'] = (ok, 'prod %d rows/%d probes, bodyforce %d rows/%d probes, '
                 '%d unparsed' % (len(dp['rows']), dp['probes'],
                                  len(db['rows']), db['probes'],
                                  len(dp['unparsed']) + len(db['unparsed'])))
    if not ok:
        return chk, met, None

    # B2 -- alpha in BOTH runs
    worst_a = max(abs(a - 1.0) for a in p['alpha'] + b['alpha'])
    met['worst_alpha'] = worst_a
    chk['B2'] = (worst_a <= TOL_ALPHA,
                 'max |alpha-1| = %.3g (tol %.3g)' % (worst_a, TOL_ALPHA))

    # B3 -- Linear force scaling is NOT rolled back (BODYFORCE run only).
    # The production deck has no body force by design, so its |base| == 0
    # is EXPECTED and must never be judged here (CI #85 lesson).
    ts = [i / float(len(b['base']) - 1) for i in range(len(b['base']))]
    k, b0, r2, resid = linfit(ts, b['base'])
    met.update(lin_k=k, lin_b=b0, lin_r2=r2, lin_resid=resid,
               bf_base_max=max(b['base']), prod_base_max=max(p['base']))
    ok = (max(b['base']) >= BASE_MIN_END and min(b['base']) <= BASE_MAX_START
          and r2 >= R2_MIN and k >= SLOPE_MIN and resid <= LIN_RESID_MAX)
    chk['B3'] = (ok, 'BODYFORCE |base| %.4g -> %.4g, slope %.4g, R2 %.5f, '
                 'max resid %.4g (prod |base|max %.4g, expected 0)'
                 % (min(b['base']), max(b['base']), k, r2, resid,
                    max(p['base'])))

    # B4 -- true superposition, WITH a non-zero base so it can discriminate.
    worst_s = max(s / max(1.0, abs(f)) for s, f in zip(b['serr'], b['force']))
    has_power = max(b['base']) >= BASE_MIN_END
    met['worst_superpos_rel'] = worst_s
    met['b4_has_power'] = has_power
    ok = worst_s <= TOL_SUPERPOS and has_power
    power = 'with |base|>0 discriminating power' if has_power else \
            'NO POWER: |base|==0 (overwrite would also give err=0)'
    chk['B4'] = (ok, 'BODYFORCE max superpos_err/|force| = %.3g (tol %.3g), %s'
                 % (worst_s, TOL_SUPERPOS, power))

    # B5 -- Biot term is live and time-dependent (production run)
    gmin, gmax = min(p['grad']), max(p['grad'])
    evol = gmax / gmin if gmin > 0 else float('inf')
    excess = max(f - bb for f, bb in zip(p['force'], p['base']))
    met.update(grad_min=gmin, grad_max=gmax, grad_evol=evol,
               force_excess=excess)
    ok = (gmin > 0.0 and evol >= GRAD_EVOL_MIN
          and excess >= FORCE_EXCESS_MIN
          and max(p['force']) > max(p['base']))
    chk['B5'] = (ok, 'prod |gradp| %.4g..%.4g (ratio %.3g), '
                 'max(|force|-|base|) %.4g' % (gmin, gmax, evol, excess))

    # B6 -- NOBIOT control. CI #85 found an EXACTLY ZERO control field; do
    # not compute pm/0 -> inf (the inf was not caught by B7 either).
    try:
        pm, nm = peak_ux(main_mesh)
        pn, nn = peak_ux(nobiot_mesh)
        met.update(peak_main=pm, peak_nobiot=pn, nrec_main=nm,
                   nrec_nobiot=nn)
        in_band = UX_PEAK_LO <= pm <= UX_PEAK_HI
        if pn <= UX_ZERO_TOL:
            ratio = None
            ok = in_band
            detail = ('peak|u_x| main %.5g m (%d rec) vs nobiot %.5g m '
                      '(%d rec), NOBIOT identically zero (exact decoupling)'
                      % (pm, nm, pn, nn))
        else:
            ratio = pm / pn
            ok = in_band and ratio >= UX_RATIO_MIN
            detail = ('peak|u_x| main %.5g m (%d rec) vs nobiot %.5g m '
                      '(%d rec), ratio %.4g (need >= %.4g)'
                      % (pm, nm, pn, nn, ratio, UX_RATIO_MIN))
        met['ux_ratio'] = ratio
        chk['B6'] = (ok, detail)
    except Exception as exc:                       # noqa: BLE001 - gate fails
        chk['B6'] = (False, 'mesh read FAILED: %s' % exc)

    # B7 -- completeness of BOTH runs
    tp, tb = read_log_text(prod_log), read_log_text(bf_log)
    crash = bool(CRASH_RE.search(tp) or CRASH_RE.search(tb))
    unknown = 'Unknown option' in tp or 'Unknown option' in tb
    ok = (dp['probes'] >= MIN_PROBES and db['probes'] >= MIN_PROBES
          and not crash and not unknown)
    chk['B7'] = (ok, 'prod %d / bodyforce %d probes, crash=%s, '
                 'unknown_option=%s' % (dp['probes'], db['probes'],
                                        crash, unknown))
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
    if met.get('p_rows'):
        for label, key in (('production (body-force free)', 'p'),
                           ('BODYFORCE (Force=0,-9.81,0)', 'b')):
            d = met[key]
            n = len(d['base'])
            print('  --- %s: first / middle / last ---' % label)
            for i in (0, n // 2, n - 1):
                print('    step %2d  |base|=%10.4g  |gradp|=%10.4g  '
                      '|force|=%10.4g  err=%8.3g'
                      % (i + 1, d['base'][i], d['grad'][i],
                         d['force'][i], d['serr'][i]))
    fails = [k for k, v in chk.items() if not v[0]]
    return fails


def main(argv):
    gate = bool(argv) and argv[0] == '--gate'
    args = argv[1:] if gate else argv
    usage = '<prod.log> <prod.gen> <nobiot.gen> <bodyforce.log>'
    if gate and len(args) != 4:
        print('GATE41_RESULT: FAIL (usage: --gate %s)' % usage)
        return 1
    if not gate and len(args) != 4:
        print('usage: check_gate41_biot.py [--gate] %s' % usage)
        return 0
    chk, met, err = evaluate(args[0], args[1], args[2], args[3])
    fails = print_report(chk, met, err)
    if gate:
        if err or fails or len(chk) < 7:
            print('GATE41_RESULT: FAIL (%s)'
                  % (err if err else ','.join(fails) or 'incomplete criteria'))
            return 1
        ratio = met.get('ux_ratio')
        rtxt = ('NOBIOT field zero' if ratio is None
                else 'NOBIOT ratio %.4g' % ratio)
        print('GATE41_RESULT: PASS (B1-B7; prod %d steps, BODYFORCE %d steps, '
              'Linear scaling preserved with non-zero |base|, '
              'superposition error ~0, %s)'
              % (met['p_rows'], met['b_rows'], rtxt))
        return 0
    print('VERDICT: %s (report mode, informational only)'
          % ('FAIL ' + ','.join(fails) if fails else 'PASS'))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
