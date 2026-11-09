#!/usr/bin/env python3
"""Gate 4.4 -- Darcy SOURCE SIGN verdict (Bourdin G-F contract).

Pure standard library + the committed ``parse_cdf2.py`` (the CI runner has
no numpy / scipy / PyYAML).

Two modes
---------
  report : python3 check_gate44_source.py <source.log> <source.gen> <main.gen>
      INFORMATIONAL. Always exits 0, never a gate.

  gate   : python3 check_gate44_source.py --gate <source.log> <source.gen> <main.gen>
      HARD Gate 4.4. Prints ``GATE44_RESULT: PASS|FAIL`` and exits 0 / 1.
      Missing / empty / short / unreadable / non-finite data is ALWAYS FAIL
      (CI #73 lesson: "no data" must never read as a pass).

  selftest : python3 check_gate44_source.py --selftest
      Offline analytic replay + negative controls. No PETSc, no solver, no CI.

Why this gate exists
--------------------
MEF90/m_MEF90_DiffusionImplementation.F90:27-60 defines the diffusion
contract

    -div(A grad p) + lambda p = f
    Operator G = +int(A grad p . grad w + lambda p w)      [ADD_VALUES]
    RHS      F = +int(f w)                                 [ADD_VALUES]
    SNES Function MUST be  G - F , Jacobian is K = dG/dp.

MEF90HeatXfer honours it by negating the cell flux before calling
MEF90DiffusionRHSSetCell (m_MEF90_HeatXferAssembly.F90:113 and :476:
``modifiedFluxPtr = -fluxPtr``).  Before Gate 4.4, m_MEF90_FlowAssembly.F90
passed the PHYSICAL +f straight through, so MEF90FlowOperatorCore assembled
G + F: a positive deck source acted as a SINK.  Every production and
diagnostic deck had ``cs0001.source: 0.0``, so none of the Gate 3.3 / 4.1 /
4.2 / 4.3 / 4.1b regressions could see it.  Gate 3.3 P6 (Terzaghi decay)
is still correct: with f=0 the sign of f is unobservable.

Isolation strategy (the project's "two decks differ by ONE line" method)
------------------------------------------------------------------------
``terzaghi1D_source.yaml`` is byte-identical to ``terzaghi1D.yaml`` except
``cs0001.source``: 0.0 -> 1.0e-6.  p0, BCs, geometry, material and time grid
are unchanged.  The Flow pressure solve is linear in the RHS (k(v)=k_min is
constant for this pure-poroelastic benchmark, see Gate 3.3 P6), therefore

    delta_p(n,j) = p_source(n,j) - p_main(n,j)

is the PURE source response, with the whole 1e5 Pa Terzaghi decay cancelled.
S2 makes the one-line isolation itself a hard criterion, so a deck that
silently changes p0 or a BC cannot sneak through.

Analytic reference (1D, H=10 m, s=5-x, drained at s=H, closed at s=0)
----------------------------------------------------------------------
    S dp/dt - A p'' = f,  A = k_min/mu = 5e-8,  S = 4.6e-10
    p(s,0)=0, p(H,t)=0, p'(0,t)=0

Steady part           p_ss(s) = f/(2A) (H^2 - s^2)
Modes                 lambda_m = (2m+1) pi/(2H)
Backward-Euler decay  (1 + c_v lambda_m^2 dt)^(-n),  c_v = A/S
Source response       p_f(s,n) = p_ss(s) + sum_m a_m cos(lambda_m s)
                                          (1 + c_v lambda_m^2 dt)^(-n)
with a_m = -2 (-1)^m f / (A H lambda_m^3), i.e. the cosine series of -p_ss,
so p_f(s,0) = 0 exactly.  The production baseline is the usual Terzaghi BE
series p_T(s,n) = P0 sum_m (4/pi)(-1)^m/(2m+1) cos(lambda_m s) decay_m(n).

With the OLD (unfixed) code the numerical field would be p_T - p_f, so
cos(delta_p, p_f) = -1 and S4/S5/S6/S7 all FAIL.  See --selftest.
"""
import math
import os
import re
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
try:
    from parse_cdf2 import CDF2
except Exception as exc:                                # pragma: no cover
    print('GATE44_RESULT: FAIL (cannot import parse_cdf2: %s)' % exc)
    sys.exit(1)

# --- Physics / discretisation (terzaghi1D.yaml, CI #66) -------------------
P0 = 1.0e5
K_MIN = 5.0e-11
MU = 1.0e-3
A = K_MIN / MU                 # k_min / viscosity [m^2/(Pa s)] = 5e-8
STOR = 4.6e-10                 # storativity [1/Pa]
CV = A / STOR                  # consolidation coefficient [m^2/s]
H = 10.0                       # drainage path [m], x in [-5,5]
DT = 0.02                      # fixed BE step [s]
SOURCE = 1.0e-6                # cs0001.source in terzaghi1D_source.yaml
NSERIES = 300                  # cosine modes (CI awk uses 200; 300 is safe)
MIN_RECORDS = 50

# --- Hard thresholds -----------------------------------------------------
COS_MIN = 0.95                 # delta_p must ALIGN with +p_f (sign + shape)
RATIO_LO, RATIO_HI = 0.70, 1.30   # ||delta_p|| / ||p_f||
L2_MAX = 0.15                  # FE vs continuum BE series (Gate 3.3 P6 ~0.032)
BASE_L2_MAX = 0.05             # production baseline vs Terzaghi BE series
PEAK_LO, PEAK_HI = 0.70, 1.30 # final closed-end delta_p / analytic p_f
TOL_GEOM = 1.0e-9

DECK_MAIN = os.path.join(_HERE, 'terzaghi1D.yaml')
DECK_SOURCE = os.path.join(_HERE, 'terzaghi1D_source.yaml')

FATAL_RE = re.compile(
    r'Caught signal|SEGV|MPI_ABORT|PETSC ERROR|Unknown option|'
    r'attempt to deallocate|Error opening EXO file|TSSolve failed', re.I)
BAG_SOURCE_RE = re.compile(r'^\s{2}source\s*=\s*(.*?)\s*;', re.M)


# --------------------------------------------------------------------------
# Analytic reference
# --------------------------------------------------------------------------
_LAM = [(2.0 * m + 1.0) * math.pi / (2.0 * H) for m in range(NSERIES)]
_DECAY_BASE = [1.0 + CV * l * l * DT for l in _LAM]
_AMP_T = [(4.0 / math.pi) * ((-1.0) ** m) / (2.0 * m + 1.0)
          for m in range(NSERIES)]
_AMP_F = [-2.0 * ((-1.0) ** m) * SOURCE / (A * H * l ** 3)
          for m, l in enumerate(_LAM)]


def _decay(n):
    return [b ** (-n) for b in _DECAY_BASE]


def build_analytic(x):
    """Return (p_T, p_f) as [record][node] for a given x-coordinate list."""
    z = [5.0 - float(xi) for xi in x]
    p_ss = [0.0 if zi >= H - 1e-12 else SOURCE / (2.0 * A) * (H * H - zi * zi)
            for zi in z]
    cos_tab = [[math.cos(l * zi) for l in _LAM] for zi in z]
    return z, p_ss, cos_tab


def analytic_fields(z, p_ss, cos_tab, times):
    """p_T (production Terzaghi decay) and p_f (pure source response)."""
    p_t, p_f = [], []
    for t in times:
        n = int(round(t / DT))
        dec = _decay(n)
        row_t, row_f = [], []
        for j, zi in enumerate(z):
            ct = cos_tab[j]
            if zi >= H - 1e-12:
                row_t.append(0.0)
                row_f.append(0.0)
                continue
            row_t.append(P0 * sum(a * c * d for a, c, d in
                                  zip(_AMP_T, ct, dec)))
            if n <= 0:
                row_f.append(0.0)
            else:
                row_f.append(p_ss[j] + sum(a * c * d for a, c, d in
                                           zip(_AMP_F, ct, dec)))
        p_t.append(row_t)
        p_f.append(row_f)
    return p_t, p_f


# --------------------------------------------------------------------------
# IO
# --------------------------------------------------------------------------
def load_mesh(path):
    """Read coordx, time_whole and nodal pressure (var1) from an Exodus mesh."""
    if not path or not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise ValueError('missing/empty mesh: %s' % path)
    c = CDF2(path)
    x = [float(v) for v in c.read_static('coordx')]
    nn = len(x)
    if nn <= 0 or not all(math.isfinite(v) for v in x):
        raise ValueError('bad coordx in: %s' % path)
    times, p = [], []
    for i in range(c.numrecs):
        times.append(float(c.read_record('time_whole', i)[0]))
        p.append([float(v) for v in c.read_record('vals_nod_var1', i)])
    if not all(math.isfinite(t) for t in times):
        raise ValueError('non-finite time axis in: %s' % path)
    for row in p:
        if len(row) != nn or not all(math.isfinite(v) for v in row):
            raise ValueError('non-finite / mismatched pressure field in: %s'
                             % path)
    return dict(path=path, nn=nn, nrec=c.numrecs, x=x, times=times, p=p)


def read_log_text(path):
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as fh:
            return fh.read()
    except Exception:                                  # noqa: BLE001 - gate fails
        return ''


def strip_yaml(path):
    """Non-comment, non-blank lines of a vDef option deck."""
    out = []
    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        for raw in fh:
            line = raw.split('#', 1)[0].rstrip()
            if line.strip():
                out.append(line.rstrip())
    return out


def deck_isolation():
    """S2a: the two decks must differ in EXACTLY the cs0001.source line."""
    try:
        a, b = strip_yaml(DECK_MAIN), strip_yaml(DECK_SOURCE)
    except Exception as exc:                           # noqa: BLE001
        return False, 'deck read failed: %s' % exc
    if len(a) != len(b):
        return False, 'deck line count differs (%d vs %d)' % (len(a), len(b))
    diffs = [(i, la.strip(), lb.strip()) for i, (la, lb) in enumerate(zip(a, b))
             if la.strip() != lb.strip()]
    if len(diffs) != 1:
        return False, ('expected exactly 1 differing option line, got %d: %s'
                       % (len(diffs), [d[1:3] for d in diffs[:4]]))
    _, old, new = diffs[0]
    if not (old.startswith('source:') and new.startswith('source:')):
        return False, 'the differing line is not cs0001.source: %r -> %r' % (old, new)
    try:
        v_old, v_new = float(old.split(':', 1)[1]), float(new.split(':', 1)[1])
    except ValueError:
        return False, 'unparseable source values: %r -> %r' % (old, new)
    ok = abs(v_old) < 1e-30 and abs(v_new - SOURCE) <= 1e-12
    return ok, ('only difference: cs0001.source %g -> %g (want 0 -> %g)'
                % (v_old, v_new, SOURCE))


def bag_source_ok(log_text):
    """S2b: the solver's Flow CellSet bag must really carry the non-zero source."""
    vals = []
    for raw in BAG_SOURCE_RE.findall(log_text):
        try:
            vals.append(float(raw.strip()))
        except ValueError:
            continue
    nz = [v for v in vals if abs(v) > 1e-30]
    ok = any(abs(v - SOURCE) <= SOURCE * 1e-3 for v in nz)
    return ok, 'bag "source" values=%s, non-zero matching %g: %s' % (
        vals[:8], SOURCE, ok)


# --------------------------------------------------------------------------
# Physics
# --------------------------------------------------------------------------
def _flat(rows):
    return [v for row in rows for v in row]


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def _nrm(a):
    return math.sqrt(max(_dot(a, a), 0.0))


def _cos(a, b):
    na, nb = _nrm(a), _nrm(b)
    if na <= 0.0 or nb <= 0.0:
        return None
    return _dot(a, b) / (na * nb)


def _rel_l2(a, b):
    nb = _nrm(b)
    if nb <= 0.0:
        return float('inf')
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b))) / nb


def physics_checks(d_src, d_main):
    """S3-S8: compare the measured source response with the analytic p_f."""
    chk, met = {}, {}
    x, times = d_src['x'], d_src['times']
    z, p_ss, cos_tab = build_analytic(x)
    p_t, p_f = analytic_fields(z, p_ss, cos_tab, times)

    # S3 -- production baseline really is the Terzaghi BE series
    base_l2 = _rel_l2(_flat(d_main['p']), _flat(p_t))
    met['base_l2'] = base_l2
    chk['S3'] = (base_l2 <= BASE_L2_MAX,
                 'production p vs Terzaghi BE series relL2=%.5f (<= %.2f)'
                 % (base_l2, BASE_L2_MAX))

    delta = [[a - b for a, b in zip(rs, rm)]
             for rs, rm in zip(d_src['p'], d_main['p'])]
    fd, fp = _flat(delta), _flat(p_f)
    cos_d = _cos(fd, fp)
    ratio = (_nrm(fd) / _nrm(fp)) if _nrm(fp) > 0.0 else float('inf')
    l2 = _rel_l2(fd, fp)
    met.update(cos_delta=cos_d, amp_ratio=ratio, l2_delta=l2,
               delta_norm=_nrm(fd), analytic_norm=_nrm(fp))
    chk['S4'] = (cos_d is not None and cos_d >= COS_MIN,
                 'cos(delta_p, +p_f)=%s (want >= %.2f; -1 means OLD G+F bug)'
                 % ('n/a' if cos_d is None else '%.5f' % cos_d, COS_MIN))
    chk['S5'] = (RATIO_LO <= ratio <= RATIO_HI,
                 '||delta_p||/||p_f||=%.4f (want [%.2f,%.2f])'
                 % (ratio, RATIO_LO, RATIO_HI))
    chk['S6'] = (l2 <= L2_MAX,
                 'relL2(delta_p,p_f)=%.5f (<= %.2f)' % (l2, L2_MAX))

    # S7 -- closed end (z=0, x=+5) final increment: positive, analytic band
    closed = [j for j, zi in enumerate(z) if zi <= TOL_GEOM]
    if not closed:
        chk['S7'] = (False, 'no closed-end (x=+5) nodes found')
    else:
        meas = sum(delta[-1][j] for j in closed) / len(closed)
        ref = sum(p_f[-1][j] for j in closed) / len(closed)
        rel = (meas / ref) if abs(ref) > 0.0 else float('inf')
        met.update(closed_delta=meas, closed_analytic=ref, closed_rel=rel)
        ok = meas > 0.0 and PEAK_LO <= rel <= PEAK_HI
        chk['S7'] = (ok, 'final closed-end delta_p=%.4f Pa, analytic=%.4f Pa, '
                     'ratio=%.4f (want >0 and in [%.2f,%.2f])'
                     % (meas, ref, rel, PEAK_LO, PEAK_HI))

    # S8 -- the source response must GROW monotonically from 0
    if not closed:
        chk['S8'] = (False, 'no closed-end nodes for time monotonicity')
    else:
        probes = [0, 10, 25, min(50, len(delta) - 1)]
        seq = [sum(delta[n][j] for j in closed) / len(closed) for n in probes]
        met['closed_series'] = list(zip(probes, seq))
        ok = (abs(seq[0]) <= max(1.0, 0.02 * abs(seq[-1]))
              and all(b > a for a, b in zip(seq, seq[1:])))
        chk['S8'] = (ok, 'closed-end delta_p at n=%s = %s (want 0 then strictly '
                     'increasing)' % (probes, ['%.3f' % v for v in seq]))
    return chk, met


def evaluate(src_log, src_gen, main_gen):
    """Return (checks, metrics, load_error). checks: name -> (ok, detail)."""
    chk, met = {}, {}
    log = read_log_text(src_log)
    if not log:
        return chk, met, 'source log unreadable/empty: %s' % src_log
    chk['S1a'] = (not bool(FATAL_RE.search(log)),
                  'source log: crash/unknown-option signature=%s'
                  % bool(FATAL_RE.search(log)))

    ok, detail = deck_isolation()
    chk['S2a'] = (ok, 'deck isolation -- ' + detail)
    ok, detail = bag_source_ok(log)
    chk['S2b'] = (ok, 'bag wiring -- ' + detail)

    try:
        ds, dm = load_mesh(src_gen), load_mesh(main_gen)
    except Exception as exc:                           # noqa: BLE001 - gate fails
        chk['S1b'] = (False, 'mesh read FAILED: %s' % exc)
        return chk, met, None
    if ds['nrec'] < MIN_RECORDS or dm['nrec'] < MIN_RECORDS:
        chk['S1b'] = (False, 'too few records (source %d, main %d, need >= %d)'
                      % (ds['nrec'], dm['nrec'], MIN_RECORDS))
        return chk, met, None
    if (ds['nn'] != dm['nn'] or ds['nrec'] != dm['nrec']
            or len(ds['times']) != len(dm['times'])):
        chk['S1b'] = (False, 'mesh shape mismatch: src nn=%d nrec=%d, '
                      'main nn=%d nrec=%d' % (ds['nn'], ds['nrec'],
                                              dm['nn'], dm['nrec']))
        return chk, met, None
    tdiff = max((abs(a - b) for a, b in zip(ds['times'], dm['times'])),
                default=float('inf'))
    xdiff = max((abs(a - b) for a, b in zip(ds['x'], dm['x'])),
                default=float('inf'))
    ok = tdiff <= 1e-9 and xdiff <= TOL_GEOM
    chk['S1b'] = (ok, 'mesh integrity: nn=%d nrec=%d, max|dt|=%.3g, '
                  'max|dx|=%.3g' % (ds['nn'], ds['nrec'], tdiff, xdiff))
    if not ok:
        return chk, met, None

    more, mmet = physics_checks(ds, dm)
    chk.update(more)
    met.update(mmet)
    return chk, met, None


def print_report(chk, met, err):
    print('=== Gate 4.4 -- Darcy source sign (Bourdin G-F) verdict ===')
    if err:
        print('  DATA ERROR: %s' % err)
    for name in sorted(chk):
        ok, detail = chk[name]
        print('  %-4s %-4s %s' % (name, 'PASS' if ok else 'FAIL', detail))
    if met:
        print('  --- metrics ---')
        for key in ('base_l2', 'cos_delta', 'amp_ratio', 'l2_delta',
                    'delta_norm', 'analytic_norm', 'closed_delta',
                    'closed_analytic', 'closed_rel', 'closed_series'):
            if key in met:
                print('    %-20s = %s' % (key, met[key]))
    return [k for k, v in chk.items() if not v[0]]


# --------------------------------------------------------------------------
# Offline self-test: analytic replay + negative controls (no PETSc, no CI)
# --------------------------------------------------------------------------
def selftest():
    x = [5.0 - float(i) for i in range(11)]     # z = 0..10
    times = [n * DT for n in range(51)]
    z, p_ss, cos_tab = build_analytic(x)
    p_t, p_f = analytic_fields(z, p_ss, cos_tab, times)
    dm = dict(nn=len(x), nrec=len(times), x=x, times=times, p=p_t)

    def mesh(fn):
        return dict(nn=len(x), nrec=len(times), x=x, times=times,
                    p=[[fn(p_t[n][j], p_f[n][j]) for j in range(len(x))]
                       for n in range(len(times))])

    cases = [
        ('positive  (fixed G-F)      ', lambda t, f: t + f, True),
        ('negative  (OLD G+F bug)    ', lambda t, f: t - f, False),
        ('zero      (source not wired)', lambda t, f: t, False),
        ('doubled   (wrong magnitude)', lambda t, f: t + 2.0 * f, False),
        ('reversed  (wrong time dir) ', lambda t, f: t + f, False),
    ]
    print('=== Gate 4.4 offline selftest (analytic replay, no solver) ===')
    print('  analytic p_f(0,t=1s) = %.4f Pa; p_f(H,t) = %.4f Pa'
          % (p_f[-1][0], p_f[-1][-1]))
    all_ok = True
    for idx, (label, fn, want) in enumerate(cases):
        if idx == 4:                       # time-reversed source response
            src = mesh(lambda t, f: t + f)
            src['p'] = list(reversed(src['p']))
        else:
            src = mesh(fn)
        chk, _ = physics_checks(src, dm)
        fails = [k for k, v in chk.items() if not v[0]]
        got = not fails
        ok = (got == want)
        all_ok &= ok
        print('  %s -> %-4s (want %-4s) %s  failing=%s'
              % (label, 'PASS' if got else 'FAIL', 'PASS' if want else 'FAIL',
                 'OK' if ok else 'SELFTEST ERROR', ','.join(fails) or 'none'))
    iso, detail = deck_isolation()
    print('  deck isolation             -> %s  %s'
          % ('PASS' if iso else 'FAIL', detail))
    all_ok &= iso
    print('SELFTEST_RESULT: %s' % ('PASS' if all_ok else 'FAIL'))
    return 0 if all_ok else 1


def main(argv):
    if argv and argv[0] == '--selftest':
        return selftest()
    gate = bool(argv) and argv[0] == '--gate'
    args = argv[1:] if gate else argv
    usage = '<source.log> <source.gen> <main.gen>'
    if len(args) != 3:
        if gate:
            print('GATE44_RESULT: FAIL (usage: --gate %s)' % usage)
            return 1
        print('usage: check_gate44_source.py [--gate|--selftest] %s' % usage)
        return 0
    chk, met, err = evaluate(args[0], args[1], args[2])
    fails = print_report(chk, met, err)
    if gate:
        if err or fails or len(chk) < 8:
            print('GATE44_RESULT: FAIL (%s)'
                  % (err if err else ','.join(fails) or 'incomplete criteria'))
            return 1
        print('GATE44_RESULT: PASS (S1 data, S2 one-line deck isolation + bag, '
              'S3 Terzaghi baseline, S4-S6 delta_p aligns with +p_f '
              '[cos=%s, amp=%.4f, L2=%.4f], S7 positive closed-end peak, '
              'S8 monotone growth)'
              % ('%.4f' % met['cos_delta'] if met.get('cos_delta') is not None
                 else 'n/a',
                 met.get('amp_ratio', 0.0), met.get('l2_delta', 0.0)))
        return 0
    print('VERDICT: %s (report mode, informational only)'
          % ('FAIL ' + ','.join(fails) if fails else 'PASS'))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
