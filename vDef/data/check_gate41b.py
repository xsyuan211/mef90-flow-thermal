#!/usr/bin/env python3
"""Gate 4.1b -- effective-stress (total-stress) Biot form verdict.

Pure standard library + the committed ``parse_cdf2.py`` (the CI runner has
no numpy / scipy / matplotlib).

Two modes
---------
  report : python3 check_gate41b.py <stress.log> <main.log> <stress.gen> <main.gen>
      INFORMATIONAL. Always exits 0, never a gate.

  gate   : python3 check_gate41b.py --gate <stress.log> <main.log> <stress.gen> <main.gen>
      HARD Gate 4.1b. Prints ``GATE41B_RESULT: PASS|FAIL`` and exits 0 / 1.
      Missing / empty / short / unreadable / non-finite data is ALWAYS FAIL
      (CI #73 lesson: "no data" must never read as a pass).

Background (FLOW_FIX_CHECKLIST.md, item 110-117)
------------------------------------------------
Gate 4.1b makes the "Biot coupling shows up as effective stress" requirement a
form that is MUTUALLY EXCLUSIVE with the 4.1a body-force superposition:

  * body-force (4.1a, default):     residual gets  f = baseForce - alpha*grad(p)
  * effective-stress (4.1b):        residual gets  sigma_total = sigma' - alpha*p*I
                                    (the body-force -alpha*grad(p) is SUPPRESSED)

Integration by parts proves the two weak forms share the SAME volume term,
    int(-alpha*grad p).v = int(alpha*p*div v) - int_boundary(alpha*p*v.n),
so they differ ONLY by the natural boundary traction -alpha*p*n. The Terzaghi
deck pins BOTH ends (cs0010 and cs0020 carry displacementBC 1,1,1), so the
boundary term is identically zero and the two forms must produce BIT-EQUIVALENT
displacement fields (up to KSP tolerance, since both runs are -disp_snes_type
ksponly -- a single exact linear solve of an effectively-linear system).

The Fortran side (m_MEF90_FlowAssembly.F90, MEF90FlowApplyBiotForce) returns
early in the effective-stress branch and writes ONE ``Gate 4.1b: form=stress``
row per outer step; it never reaches the body-force ``Gate 4.1:`` diagnostic.
So the four criteria are all observable from the logs + result meshes:

  E1  stress deck log contains "Gate 4.1b: form=stress" (effective form active)
  E2  stress displacement field is equivalent to the superposition (production)
      field: relative L2 over u_x,u_y across all 51 records <= 1e-3, AND the
      stress deck's peak|u_x| stays in the Terzaghi band [5e-5, 2e-4] m
      (physics sanity so a zero / blown-up field cannot read as "equivalent")
  E3  production (superposition) deck log does NOT contain "Gate 4.1b:" -- the
      default build is unchanged (no unexpected stress diagnostic), AND it still
      emits its own "Gate 4.1:" body-force rows (the 4.1a path is intact)
  E4  mutual exclusion: the stress deck log has "Gate 4.1b:" but NO "Gate 4.1:"
      body-force rows, proving the two forms cannot both act on one solve

Note: E2 compares the stress deck against the PRODUCTION (superposition) mesh,
not the body-force diagnostic deck. The body-force diagnostic deck
(terzaghi1D_bodyforce.yaml) carries a non-zero preset Force 0,-9.81,0 whose
baseForce would contaminate the equivalence test; the production deck has
Force 0,0,0, exactly matching the stress deck's zero base force, so the only
difference between the two solves is -alpha*p*I (stress) vs -alpha*grad(p)
(body force) -- precisely the form equivalence under test.
"""
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse_cdf2 import CDF2  # noqa: E402

FORM_RE = re.compile(r'Gate\s+4\.1b:\s*form=stress')
BODYFORCE_ROW_RE = re.compile(r'Gate\s+4\.1:')  # 4.1a body-force diagnostic (no 'b')
CRASH_RE = re.compile(r'Caught signal|SEGV|MPI_ABORT|PETSC ERROR', re.I)

TOL_REL_L2 = 1.0e-3            # relative L2 of the displacement field
UX_PEAK_LO, UX_PEAK_HI = 5.0e-5, 2.0e-4   # Terzaghi u_x band (Gate 3.3 P4)
MIN_RECORDS = 50


def load_disp(path):
    """Load u_x (var3) / u_y (var4) + times from an Exodus mesh; return dict."""
    if not path or not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise ValueError('missing/empty mesh: %s' % path)
    c = CDF2(path)
    nn = len(c.read_static('coordx'))
    if nn <= 0:
        raise ValueError('no nodes (coordx) in: %s' % path)
    ux, uy, times = [], [], []
    for i in range(c.numrecs):
        times.append(float(c.read_record('time_whole', i)[0]))
        ux.append([float(v) for v in c.read_record('vals_nod_var3', i)])
        uy.append([float(v) for v in c.read_record('vals_nod_var4', i)])
    if not all(math.isfinite(t) for t in times):
        raise ValueError('non-finite time axis in: %s' % path)
    for key in ('ux', 'uy'):
        for f in (ux if key == 'ux' else uy):
            if not f or len(f) != nn or not all(math.isfinite(v) for v in f):
                raise ValueError('non-finite / mismatched %s field in: %s'
                                 % (key, path))
    return dict(path=path, nn=nn, nrec=c.numrecs, times=times, ux=ux, uy=uy)


def rel_l2(a, b):
    """Relative L2 of (a - b) over u_x,u_y; returns (rel, abs_b)."""
    num = 0.0
    den = 0.0
    for ra, rb in zip(a['ux'], b['ux']):
        for xa, xb in zip(ra, rb):
            num += (xa - xb) ** 2
            den += xb * xb
    for ra, rb in zip(a['uy'], b['uy']):
        for ya, yb in zip(ra, rb):
            num += (ya - yb) ** 2
            den += yb * yb
    if den <= 0.0:
        return float('inf'), 0.0
    return math.sqrt(num / den), math.sqrt(den)


def peak_ux(d):
    return max((max(abs(v) for v in rec) for rec in d['ux']), default=0.0)


def read_log_text(path):
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as fh:
            return fh.read()
    except Exception:                              # noqa: BLE001 - gate fails
        return ''


def evaluate(stress_log, main_log, stress_gen, main_gen):
    """Return (checks, metrics, load_error). checks: name -> (ok, detail)."""
    chk, met = {}, {}

    # --- E1 / E3 / E4: parse the two logs (text only, cheap) ----------------
    ts = read_log_text(stress_log)
    tm = read_log_text(main_log)
    if not ts:
        return chk, met, 'stress log unreadable/empty: %s' % stress_log
    if not tm:
        return chk, met, 'production log unreadable/empty: %s' % main_log

    stress_has_form = bool(FORM_RE.search(ts))
    stress_has_bodyforce_row = bool(BODYFORCE_ROW_RE.search(ts))
    main_has_form = bool(FORM_RE.search(tm))
    main_has_bodyforce_row = bool(BODYFORCE_ROW_RE.search(tm))
    met.update(stress_has_form=stress_has_form,
               stress_has_bodyforce_row=stress_has_bodyforce_row,
               main_has_form=main_has_form,
               main_has_bodyforce_row=main_has_bodyforce_row)

    # E1 -- effective-stress form is active in the stress deck
    chk['E1'] = (stress_has_form,
                 'stress deck emitted "Gate 4.1b: form=stress": %s'
                 % stress_has_form)

    # E3 -- default (superposition) build is unchanged: no stress diagnostic,
    # but the 4.1a body-force diagnostic is still present.
    ok = (not main_has_form) and main_has_bodyforce_row
    chk['E3'] = (ok, 'production log: Gate 4.1b present=%s (want False), '
                 'Gate 4.1 present=%s (want True)' % (main_has_form,
                                                      main_has_bodyforce_row))

    # E4 -- mutual exclusion: stress deck uses the effective form AND does NOT
    # emit the body-force diagnostic (the other branch was skipped).
    ok = stress_has_form and not stress_has_bodyforce_row
    chk['E4'] = (ok, 'stress deck: Gate 4.1b=%s, Gate 4.1 (body-force)=%s '
                 '(want True / False)' % (stress_has_form,
                                          stress_has_bodyforce_row))

    # crash / unknown-option guard on the stress deck log
    crash = bool(CRASH_RE.search(ts))
    unknown = 'Unknown option' in ts
    if crash or unknown:
        chk['E0'] = (False, 'stress deck log: crash=%s, unknown_option=%s'
                     % (crash, unknown))
        return chk, met, None

    # --- E2: displacement-field equivalence --------------------------------
    try:
        ds = load_disp(stress_gen)
        dm = load_disp(main_gen)
    except Exception as exc:                       # noqa: BLE001 - gate fails
        chk['E2'] = (False, 'mesh read FAILED: %s' % exc)
        return chk, met, None

    if ds['nrec'] < MIN_RECORDS or dm['nrec'] < MIN_RECORDS:
        chk['E2'] = (False, 'too few records (stress %d, main %d, need >= %d)'
                     % (ds['nrec'], dm['nrec'], MIN_RECORDS))
        return chk, met, None
    if ds['nn'] != dm['nn'] or ds['nrec'] != dm['nrec']:
        chk['E2'] = (False, 'mesh shape mismatch: stress nn=%d nrec=%d, '
                     'main nn=%d nrec=%d' % (ds['nn'], ds['nrec'],
                                             dm['nn'], dm['nrec']))
        return chk, met, None

    rel, abs_b = rel_l2(ds, dm)
    pk = peak_ux(ds)
    in_band = UX_PEAK_LO <= pk <= UX_PEAK_HI
    met.update(e2_rel_l2=rel, e2_abs_b=abs_b, stress_peak_ux=pk,
               stress_nrec=ds['nrec'], main_nrec=dm['nrec'],
               stress_nn=ds['nn'], main_nn=dm['nn'])
    ok = (rel <= TOL_REL_L2) and in_band
    detail = ('relL2(ux,uy)=%.4g (tol %.1g), stress peak|u_x|=%.4g m '
              '[band %.0e..%.0e]: %s' % (rel, TOL_REL_L2, pk,
                                         UX_PEAK_LO, UX_PEAK_HI,
                                         'in band' if in_band else 'OUT OF BAND'))
    chk['E2'] = (ok, detail)
    return chk, met, None


def print_report(chk, met, err):
    print('=== Gate 4.1b -- effective-stress Biot form verdict ===')
    if err:
        print('  DATA ERROR: %s' % err)
    order = ['E1', 'E2', 'E3', 'E4', 'E0']
    for name in order:
        if name in chk:
            ok, detail = chk[name]
            print('  %-3s %-4s %s' % (name, 'PASS' if ok else 'FAIL', detail))
    if met.get('e2_rel_l2') is not None:
        print('  --- E2 metrics ---')
        print('    relL2(stress vs main) = %.4g' % met['e2_rel_l2'])
        print('    stress peak|u_x|      = %.5g m' % met.get('stress_peak_ux', 0.0))
        print('    records/nodes         = %d / %d' % (met.get('stress_nrec', 0),
                                                        met.get('stress_nn', 0)))
    fails = [k for k, v in chk.items() if not v[0]]
    return fails


def main(argv):
    gate = bool(argv) and argv[0] == '--gate'
    args = argv[1:] if gate else argv
    usage = '<stress.log> <main.log> <stress.gen> <main.gen>'
    if len(args) != 4:
        if gate:
            print('GATE41B_RESULT: FAIL (usage: --gate %s)' % usage)
            return 1
        print('usage: check_gate41b.py [--gate] %s' % usage)
        return 0
    chk, met, err = evaluate(args[0], args[1], args[2], args[3])
    fails = print_report(chk, met, err)
    if gate:
        if err or fails or len(chk) < 4:
            print('GATE41B_RESULT: FAIL (%s)'
                  % (err if err else ','.join(fails) or 'incomplete criteria'))
            return 1
        print('GATE41B_RESULT: PASS (E1 effective form active; E2 relL2=%.4g '
              '<= %.1g with peak|u_x|=%.4g m in band; E3 production unchanged; '
              'E4 forms mutually exclusive)'
              % (met['e2_rel_l2'], TOL_REL_L2, met['stress_peak_ux']))
        return 0
    print('VERDICT: %s (report mode, informational only)'
          % ('FAIL ' + ','.join(fails) if fails else 'PASS'))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
