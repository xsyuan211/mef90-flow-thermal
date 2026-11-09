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
Gate 4.1b makes "Biot coupling shows up as effective stress" a form that is
MUTUALLY EXCLUSIVE with the 4.1a body-force superposition:

  * body-force (4.1a, default):     residual gets  f = baseForce - alpha*grad(p)
  * effective-stress (4.1b):        residual gets  sigma_total = sigma' - alpha*p*I
                                    (the body-force -alpha*grad(p) is SUPPRESSED)

IMPORTANT CORRECTION (CI #97 offline diagnosis, 2026-11-09)
------------------------------------------------------------
The original version of this script asserted that the two forms must be
BIT-EQUIVALENT because "the Terzaghi deck pins BOTH ends, so the boundary
term is identically zero".  That premise is WRONG.

Integration by parts gives

    int(-alpha*grad p).v = int(alpha*p*div v) - int_boundary(alpha*p*v.n)

The deck pins only the two END facets x=-5 and x=+5 (n = +/-e_x).  The LATERAL
facets y=-0.5 and y=+0.5 are free Neumann surfaces: v_y is unconstrained and
p is NOT zero there (1.0e5 Pa at step 1, 3.58e4 Pa at step 25).  Hence the
boundary term is NON-ZERO.  The two forms are two DIFFERENT boundary value
problems, and no amount of assembly fixing can make the displacement fields
bit-equivalent:

  * 4.1a natural BC:  sigma'.n = 0
  * 4.1b natural BC:  sigma_total.n = 0   <=>   sigma'.n = +alpha*p*n

The second one produces an OUTWARD lateral swelling driven by the local pore
pressure.  That is a real physical effect, not a bug.  CI #97 measured (356
nodes, 51 records, out_stress.gen vs out_main.gen):

    u_y^stress ~= k*(y*p),  k ~= 7.5e-10 1/Pa, constant over all 51 steps
    cos(u_y^stress, y*p)        ~= +0.97 .. +0.98
    |u_y^stress| ~= 4.5*|u_y^main|, with NEGATIVE correlation (-0.23..-0.34)
    u_x^stress   ~= 0.69*u_x^main, cos ~= 1.000
    pressure fields IDENTICAL (max diff 0 Pa over all steps)

The "0.69" is therefore a CONSEQUENCE (Poisson-coupled axial contraction from
the lateral swelling), not an assembly scaling bug.  The Fortran assembly is
correct: for TRI3/P1 ``GradS_BF`` is element-wise constant and ``pCell`` is the
Gauss-weighted mean, so ``int(p*div phi_i) = pCell*Vol*tr(GradS_BF_i)`` holds
exactly.  See m_MEF90_DefMechAssembly.F90:322-326, :368-374 and
m_MEF90_FlowAssembly.F90:1076-1153, :1156-1249, :1269-1330.

E2 therefore no longer tests bit-equivalence.  It tests the PHYSICAL SIGNATURE
of the effective-stress form, with hard bounds:

  a. pressure drive identical            : max|pS-pM| / max|pM| <= 1e-6
  b. axial shape preserved               : cos(u_x^S, u_x^M) >= 0.99
  c. axial response shortened            : 0.40 <= ||u_x^S||/||u_x^M|| <= 0.95
                                          (upper bound < 1 catches a 4.1b that
                                          silently collapsed back into 4.1a)
  d. peak|u_x| in the Terzaghi band       : [5e-5, 2e-4] m (no zero / blow-up)
  e. excess u_y is a real field           : ||u_y^S-u_y^M|| >= 0.5*||u_y^M||
  f. excess u_y is the boundary swelling  : cos(u_y^S-u_y^M, y*p) >= 0.90
  g. swelling compliance in an envelope   : k = dot(excess,y*p)/dot(y*p,y*p)
                                          within [0.4*k_min, 2.0*k_max], where
                                          k_min/k_max are the plane-strain
                                          analytic bounds for alpha=1, E=1e9,
                                          nu=0.25 (see constants below).

E1 / E3 / E4 are unchanged (log-only); E0 crash guard is unchanged.
The old relL2 is still computed and printed, but strictly INFORMATIONAL.

Note: E2 compares the stress deck against the PRODUCTION (superposition) mesh,
not the body-force diagnostic deck.  terzaghi1D_bodyforce.yaml carries a
non-zero preset Force 0,-9.81,0 whose baseForce would contaminate the
comparison; the production deck has Force 0,0,0, exactly matching the stress
deck's zero base force, so -alpha*p*I (stress) vs -alpha*grad(p) (body force)
plus the natural-boundary traction is the ONLY difference between the solves.
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

TOL_REL_L2 = 1.0e-3            # INFORMATIONAL ONLY -- not a criterion any more
UX_PEAK_LO, UX_PEAK_HI = 5.0e-5, 2.0e-4   # Terzaghi u_x band (Gate 3.3 P4)
MIN_RECORDS = 50

# --- E2 physical-signature thresholds --------------------------------------
TOL_PRESSURE_REL = 1.0e-6      # same pressure drive in both decks
COS_UX_MIN = 0.99              # axial shape must be preserved
UX_RATIO_LO, UX_RATIO_HI = 0.40, 0.95   # swelling shortens the axial response
COS_BND_MIN = 0.90             # excess u_y must align with y*p
EXCESS_UY_MIN_REL = 0.50       # excess u_y must be a real field, not noise
TOL_GEOM = 1.0e-9              # node coordinates must agree between meshes

# Plane-strain analytic envelope for u_y ~ k*(y*p).  With alpha=1, E=1e9 Pa,
# nu=0.25 (terzaghi1D.yaml / terzaghi1D_stress.yaml cs0001):
#   axial strain fully constrained : k = (1+nu)(1-2nu)/(E(1-nu)) = 8.33e-10
#   axial stress free              : k = (1-nu^2)/E             = 9.38e-10
# CI #97 measures k ~ 7.5e-10, slightly BELOW the ideal lower bound because
# both ends are fully pinned and p is non-uniform along x.  The gate therefore
# uses a generous envelope (0.4x .. 2.0x) that locks the order of magnitude
# and the sign, without hard-coding the measured answer.
E_YOUNG, POISSON, ALPHA_BIOT = 1.0e9, 0.25, 1.0
K_AXIAL_CONSTRAINED = (ALPHA_BIOT * (1.0 + POISSON) * (1.0 - 2.0 * POISSON)
                       / (E_YOUNG * (1.0 - POISSON)))
K_AXIAL_FREE = ALPHA_BIOT * (1.0 - POISSON ** 2) / E_YOUNG
K_LO, K_HI = 0.4 * K_AXIAL_CONSTRAINED, 2.0 * K_AXIAL_FREE


def load_disp(path):
    """Load p (var1), u_x (var3), u_y (var4), coordy + times from an Exodus mesh."""
    if not path or not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise ValueError('missing/empty mesh: %s' % path)
    c = CDF2(path)
    nn = len(c.read_static('coordx'))
    if nn <= 0:
        raise ValueError('no nodes (coordx) in: %s' % path)
    try:
        ycoord = [float(v) for v in c.read_static('coordy')]
    except Exception as exc:                       # noqa: BLE001 - gate fails
        raise ValueError('coordy unreadable in %s: %s' % (path, exc))
    if len(ycoord) != nn or not all(math.isfinite(v) for v in ycoord):
        raise ValueError('bad coordy in: %s' % path)
    p, ux, uy, times = [], [], [], []
    for i in range(c.numrecs):
        times.append(float(c.read_record('time_whole', i)[0]))
        p.append([float(v) for v in c.read_record('vals_nod_var1', i)])
        ux.append([float(v) for v in c.read_record('vals_nod_var3', i)])
        uy.append([float(v) for v in c.read_record('vals_nod_var4', i)])
    if not all(math.isfinite(t) for t in times):
        raise ValueError('non-finite time axis in: %s' % path)
    for key, rows in (('p', p), ('ux', ux), ('uy', uy)):
        for f in rows:
            if not f or len(f) != nn or not all(math.isfinite(v) for v in f):
                raise ValueError('non-finite / mismatched %s field in: %s'
                                 % (key, path))
    return dict(path=path, nn=nn, nrec=c.numrecs, y=ycoord, times=times,
                p=p, ux=ux, uy=uy)


def rel_l2(a, b):
    """INFORMATIONAL ONLY: relative L2 of (a - b) over u_x,u_y.

    Kept so the verdict log still shows that the two forms ARE different;
    a non-zero value here is EXPECTED, not a failure (see module docstring).
    """
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


def _flat(rows):
    return [v for row in rows for v in row]


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def _nrm(a):
    return math.sqrt(max(_dot(a, a), 0.0))


def _cosine(a, b):
    na, nb = _nrm(a), _nrm(b)
    if na <= 0.0 or nb <= 0.0:
        return None
    return _dot(a, b) / (na * nb)


def physics_signature(ds, dm):
    """Hard E2: verify the effective-stress form's PHYSICAL SIGNATURE.

    The two forms are NOT bit-equivalent on this mesh (free lateral Neumann
    facets carry a non-zero ``int(alpha*p*v.n)`` boundary term).  Instead we
    require that the difference be exactly what the missing natural traction
    predicts: identical pressure drive, preserved u_x shape shortened by the
    Poisson-coupled axial contraction, and an extra u_y field proportional to
    ``y*p`` with an elastic compliance in the plane-strain envelope.
    """
    met = {}

    # (a) identical geometry and pressure drive
    ydiff = max((abs(u - v) for u, v in zip(ds['y'], dm['y'])), default=0.0)
    max_dp, max_p = 0.0, 0.0
    for rs, rm in zip(ds['p'], dm['p']):
        for u, v in zip(rs, rm):
            max_dp = max(max_dp, abs(u - v))
            max_p = max(max_p, abs(v))
    p_rel = (max_dp / max_p) if max_p > 0.0 else float('inf')
    met.update(y_max_diff=ydiff, p_max_diff=max_dp, p_rel=p_rel)

    # (b,c) axial shape and amplitude
    fs, fm = _flat(ds['ux']), _flat(dm['ux'])
    cos_ux = _cosine(fs, fm)
    nfm = _nrm(fm)
    ux_ratio = (_nrm(fs) / nfm) if nfm > 0.0 else float('inf')
    met.update(cos_ux=cos_ux, ux_ratio=ux_ratio)

    # (d) zero / blow-up guard on the stress deck
    pk = peak_ux(ds)
    in_band = UX_PEAK_LO <= pk <= UX_PEAK_HI
    met['stress_peak_ux'] = pk

    # (e,f,g) the extra lateral swelling u_y^S - u_y^M ~ k*(y*p)
    excess, yp = [], []
    for i in range(min(ds['nrec'], dm['nrec'])):
        for j in range(ds['nn']):
            excess.append(ds['uy'][i][j] - dm['uy'][i][j])
            yp.append(ds['y'][j] * ds['p'][i][j])
    cos_bnd = _cosine(excess, yp)
    den_yp = _dot(yp, yp)
    k = (_dot(excess, yp) / den_yp) if den_yp > 0.0 else 0.0
    den_uy = _nrm(_flat(dm['uy'])) or _nrm(_flat(ds['uy']))
    excess_rel = (_nrm(excess) / den_uy) if den_uy > 0.0 else 0.0
    met.update(cos_bnd=cos_bnd, k_compliance=k, excess_uy_rel=excess_rel,
               k_lo=K_LO, k_hi=K_HI)

    checks = [
        ('pressure-drive-identical', p_rel <= TOL_PRESSURE_REL),
        ('geometry-identical', ydiff <= TOL_GEOM),
        ('ux-shape-preserved', cos_ux is not None and cos_ux >= COS_UX_MIN),
        ('ux-shortened-not-collapsed', UX_RATIO_LO <= ux_ratio <= UX_RATIO_HI),
        ('peak-ux-in-band', in_band),
        ('excess-uy-is-real', excess_rel >= EXCESS_UY_MIN_REL),
        ('excess-uy-is-y-times-p', cos_bnd is not None and cos_bnd >= COS_BND_MIN),
        ('compliance-in-envelope', K_LO <= k <= K_HI),
    ]
    failed = [name for name, ok in checks if not ok]
    summary = ('p_rel=%.3g cos_ux=%s ux_ratio=%.4f peak|ux|=%.4g '
               'excess_rel=%.3f cos_bnd=%s k=%.3e 1/Pa [%.3e..%.3e]'
               % (p_rel, 'n/a' if cos_ux is None else '%.4f' % cos_ux,
                  ux_ratio, pk, excess_rel,
                  'n/a' if cos_bnd is None else '%.4f' % cos_bnd,
                  k, K_LO, K_HI))
    if failed:
        detail = 'FAILED [%s] :: %s' % (', '.join(failed), summary)
    else:
        detail = '%d/%d physics subchecks PASS :: %s' % (
            len(checks), len(checks), summary)
    return dict(ok=not failed, detail=detail, metrics=met,
                rel_l2=rel_l2(ds, dm)[0])


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

    chk['E1'] = (stress_has_form,
                 'stress deck emitted "Gate 4.1b: form=stress": %s'
                 % stress_has_form)
    ok = (not main_has_form) and main_has_bodyforce_row
    chk['E3'] = (ok, 'production log: Gate 4.1b present=%s (want False), '
                 'Gate 4.1 present=%s (want True)' % (main_has_form,
                                                      main_has_bodyforce_row))
    ok = stress_has_form and not stress_has_bodyforce_row
    chk['E4'] = (ok, 'stress deck: Gate 4.1b=%s, Gate 4.1 (body-force)=%s '
                 '(want True / False)' % (stress_has_form,
                                          stress_has_bodyforce_row))

    crash = bool(CRASH_RE.search(ts))
    unknown = 'Unknown option' in ts
    if crash or unknown:
        chk['E0'] = (False, 'stress deck log: crash=%s, unknown_option=%s'
                     % (crash, unknown))
        return chk, met, None

    # --- E2: PHYSICAL SIGNATURE of the effective-stress form ----------------
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

    sig = physics_signature(ds, dm)
    met.update(sig['metrics'])
    met.update(e2_rel_l2=sig['rel_l2'],
               stress_nrec=ds['nrec'], main_nrec=dm['nrec'],
               stress_nn=ds['nn'], main_nn=dm['nn'])
    chk['E2'] = (sig['ok'], sig['detail'])
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
        print('  --- E2 physics-signature metrics ---')
        print('    relL2(stress vs main) = %.4g   <-- INFORMATIONAL, NOT a '
              'criterion (non-zero boundary term is expected)' % met['e2_rel_l2'])
        rows = [
            ('p_max_diff [Pa]', 'p_max_diff'),
            ('p_rel', 'p_rel'),
            ('y_max_diff', 'y_max_diff'),
            ('cos_ux', 'cos_ux'),
            ('ux_ratio', 'ux_ratio'),
            ('stress_peak_ux [m]', 'stress_peak_ux'),
            ('excess_uy_rel', 'excess_uy_rel'),
            ('cos_bnd (excess vs y*p)', 'cos_bnd'),
            ('k_compliance [1/Pa]', 'k_compliance'),
            ('k envelope lo', 'k_lo'),
            ('k envelope hi', 'k_hi'),
            ('records / nodes', None),
        ]
        for label, key in rows:
            if key is None:
                print('    %-26s = %d / %d' % (label, met.get('stress_nrec', 0),
                                               met.get('stress_nn', 0)))
            elif key in met and met[key] is not None:
                print('    %-26s = %s' % (label, met[key]))
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
        print('GATE41B_RESULT: PASS (E1 effective form active; E2 physical '
              'signature 8/8 [identical p, cos_ux=%s, ux_ratio=%.4f, '
              'cos_bnd=%s, k=%.3e 1/Pa]; E3 production unchanged; '
              'E4 forms mutually exclusive)'
              % ('%.4f' % met['cos_ux'] if met.get('cos_ux') is not None else 'n/a',
                 met.get('ux_ratio', 0.0),
                 '%.4f' % met['cos_bnd'] if met.get('cos_bnd') is not None else 'n/a',
                 met.get('k_compliance', 0.0)))
        return 0
    print('VERDICT: %s (report mode, informational only)'
          % ('FAIL ' + ','.join(fails) if fails else 'PASS'))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
