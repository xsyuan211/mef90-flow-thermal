#!/usr/bin/env python3
"""Gate 4.5 -- Darcy ENERGY diagnostics verdict (Bourdin energy/work contract).

Pure standard library (the CI runner has no numpy / scipy / PyYAML).

Modes
-----
  report : python3 check_gate45_energy.py <main.log> <source.log>
      INFORMATIONAL. Always exits 0, never a gate.

  gate   : python3 check_gate45_energy.py --gate <main.log> <source.log>
      HARD Gate 4.5. Prints ``GATE45_RESULT: PASS|FAIL`` and exits 0 / 1.
      Missing / empty / short / unreadable / non-finite data is ALWAYS FAIL
      (CI #73 lesson: "no data" must never read as a pass).

  selftest : python3 check_gate45_energy.py --selftest
      Offline analytic replay + negative controls. No PETSc, no solver, no CI.

What Gate 4.5 fixes
-------------------
``MEF90FlowEnergy`` was exported and dispatched but NEVER called, and was
wrong even if it had been called:

  * ``permFieldSec`` / ``sourceSec`` / ``p0Sec`` were ``SectionRealDuplicate``d
    but never filled from their Vecs -> E and W read uninitialised storage.
  * ``MEF90DiffusionEnergySet`` / ``MEF90DiffusionWorkSetCell`` declare
    energy/work as INTENT(OUT) but only ACCUMULATE (``+=``).  HeatXfer zeroes
    them before the call (m_MEF90_HeatXferAssembly.F90:354-355); Flow did not.
  * the mesh argument was FlowCtx%DM (2 dof/vertex default section) while
    pSec lives on DMScalSec (1 dof/vertex).  Now DMScal, like HeatXfer.
  * per-rank partial sums were returned without an MPI reduction.

vDef.F90 now calls it once per outer time step and prints

    Gate 4.5: step N t T E <E> W <W>

where E = 1/2 int (k/mu)|grad p|^2 dOmega and W = int f*p dOmega.  It is PURE
post-processing: nothing feeds back into the Flow solve, so every pre-existing
physics gate must be bit-identical (see Gate 3.3 / 4.1 / 4.2 / 4.3 / 4.1b /
4.4 regressions).

Isolation
---------
terzaghi1D_source.yaml differs from terzaghi1D.yaml in EXACTLY ONE option
(cs0001.source 0.0 -> 1.0e-6, Gate 4.4).  Therefore:

  * the production deck has f = 0  =>  W_main == 0 for ALL 51 steps;
  * the source deck has f > 0     =>  W_src > 0 and tracks int p dOmega.

An all-zero W_src is the signature of an unwired ``sourceSec``; an all-zero E
is the signature of the old uninitialised sections.  Both are hard FAILs.

Analytic reference (1D, H=10 m, s=5-x, drained at s=H, closed at s=0)
----------------------------------------------------------------------
    A = k_min/mu = 5e-8,  S = 4.6e-10,  c_v = A/S
    lambda_m = (2m+1) pi/(2H),  d_m(n) = (1 + c_v lambda_m^2 dt)^(-n)

Terzaghi pressure  p_T(s,n) = P0 sum_m a_m cos(lambda_m s) d_m(n),
with a_m = (4/pi)(-1)^m/(2m+1).  Because a_m^2 lambda_m^2 == 4/H^2,

    E_an(n) = (A W P0^2 / H) sum_m d_m(n)^2

The n=0 field is a step (Dirichlet applied at t=0), so the continuum series
diverges there; only n>=1 is compared.  At n=50 E_an = 0.269 J.

Source work (f = 1e-6, domain width W=1 m):

    I_T(n) = int p_T ds = (8 P0 H/pi^2) sum_m d_m(n)/(2m+1)^2
    I_f(n) = int p_f  ds = f H^3/(3A) - (2f/(A H)) sum_m d_m(n)/lambda_m^4
    W_an(n) = f W (I_T + I_f)

I_f(0)=0 (the source response starts at zero).  W_an decays 0.999 -> 0.066.
"""
import math
import os
import re
import sys

# --- Physics / discretisation (terzaghi1D.yaml, CI #66) -------------------
P0 = 1.0e5
K_MIN = 5.0e-11
MU = 1.0e-3
A = K_MIN / MU                 # k_min/mu = 5e-8 [m^2/(Pa s)]
STOR = 4.6e-10
CV = A / STOR
H = 10.0                       # drainage path [m]
W_DOMAIN = 1.0                 # domain width [m] (y in [-0.5,0.5])
DT = 0.02
F_SOURCE = 1.0e-6              # cs0001.source in terzaghi1D_source.yaml
NSERIES = 400
MIN_RECORDS = 50

# --- Hard thresholds -----------------------------------------------------
TOL_W_ZERO = 1.0e-10           # |W_main| must be ~0 (f=0)
E_PEAK_MIN = 1.0e-3            # E must be a real, non-zero field
DECAY_RATIO_MIN = 5.0          # E(1)/E(50); analytic ~315
E_LAST_LO, E_LAST_HI = 0.05, 20.0   # E(50)/E_an(50), wide FE-vs-continuum band
W_DECAY_RATIO_MIN = 3.0        # W_src(0)/W_src(50); analytic ~15
COS_W_MIN = 0.90               # W_src shape vs W_an (n>=1)
W_LAST_LO, W_LAST_HI = 0.50, 1.50   # W_src(50)/W_an(50)
# The source steady state p_ss = f/(2A)(H^2-s^2) has grad = -f s/A, which
# points the SAME way as grad(p_T) (high at the closed end, 0 at the drained
# end).  The cross term 2*grad(p_T).grad(p_f) therefore RAISES E noticeably;
# E_src is NOT a small perturbation of E_main.  Keep a wide magnitude band
# that still catches a zero / blown-up E_src.
E_SRC_RATIO_LO, E_SRC_RATIO_HI = 0.30, 3.00  # E_src(50)/E_main(50)
TOL_T_END = 0.05

_LAM = [(2.0 * m + 1.0) * math.pi / (2.0 * H) for m in range(NSERIES)]
_DECAY_BASE = [1.0 + CV * l * l * DT for l in _LAM]

G45_RE = re.compile(
    r'Gate\s+4\.5:\s*step\s+(\d+)\s+t\s+([-+0-9.eEdD]+)\s+'
    r'E\s+([-+0-9.eEdD]+)\s+W\s+([-+0-9.eEdD]+)', re.I)
FATAL_RE = re.compile(
    r'Caught signal|SEGV|MPI_ABORT|PETSC ERROR|Unknown option|'
    r'attempt to deallocate|TSSolve failed', re.I)


def _num(s):
    """Parse a Fortran real (E or legacy D exponent)."""
    return float(s.replace('D', 'E').replace('d', 'e'))


def analytic(n):
    """(E_an, W_an) after exactly n backward-Euler steps."""
    dec = [b ** (-n) for b in _DECAY_BASE]
    e = (A * W_DOMAIN * P0 * P0 / H) * sum(d * d for d in dec)
    i_t = (8.0 * P0 * H / math.pi ** 2) * sum(
        d / (2.0 * m + 1.0) ** 2 for m, d in enumerate(dec))
    i_f = (F_SOURCE * H ** 3 / (3.0 * A)
           - (2.0 * F_SOURCE / (A * H)) * sum(d / l ** 4
                                              for d, l in zip(dec, _LAM)))
    return e, F_SOURCE * W_DOMAIN * (i_t + i_f)


def parse_rows(text):
    rows = []
    for m in G45_RE.finditer(text):
        try:
            rows.append((int(m.group(1)), _num(m.group(2)),
                         _num(m.group(3)), _num(m.group(4))))
        except ValueError:
            continue
    return rows


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def _nrm(a):
    return math.sqrt(max(_dot(a, a), 0.0))


def _cos(a, b):
    na, nb = _nrm(a), _nrm(b)
    if na <= 0.0 or nb <= 0.0:
        return None
    return _dot(a, b) / (na * nb)


def read_log_text(path):
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as fh:
            return fh.read()
    except Exception:                                  # noqa: BLE001 - gate fails
        return ''


def evaluate_text(main_text, src_text):
    """Return (checks, metrics, load_error); checks: name -> (ok, detail)."""
    chk, met = {}, {}
    if not main_text:
        return chk, met, 'production log unreadable/empty'
    if not src_text:
        return chk, met, 'source log unreadable/empty'

    pm, ps = parse_rows(main_text), parse_rows(src_text)
    met.update(main_rows=len(pm), source_rows=len(ps))

    # --- T1: integrity of both diagnostics --------------------------------
    ok = (len(pm) >= MIN_RECORDS and len(ps) >= MIN_RECORDS
          and len(pm) == len(ps))
    if ok:
        for rows in (pm, ps):
            steps = [r[0] for r in rows]
            ts = [r[1] for r in rows]
            vals = [v for r in rows for v in (r[2], r[3])]
            ok &= (steps == list(range(1, len(rows) + 1))
                   and all(b > a for a, b in zip(ts, ts[1:]))
                   and abs(ts[-1] - 1.0) <= TOL_T_END
                   and all(math.isfinite(v) and v >= 0.0 for v in vals))
    chk['T1'] = (ok, 'rows main=%d source=%d (need >=%d, equal, steps 1..N, '
                 'monotone t to ~1.0, finite non-negative E/W)'
                 % (len(pm), len(ps), MIN_RECORDS))
    fatal = bool(FATAL_RE.search(main_text) or FATAL_RE.search(src_text))
    chk['T1b'] = (not fatal, 'crash/unknown-option signature=%s' % fatal)
    if not ok or fatal:
        return chk, met, None

    nrec = len(pm)
    e_main = [r[2] for r in pm]
    w_main = [r[3] for r in pm]
    e_src = [r[2] for r in ps]
    w_src = [r[3] for r in ps]
    e_an = [analytic(n)[0] for n in range(nrec)]
    w_an = [analytic(n)[1] for n in range(nrec)]

    # --- T2: production W must be identically zero (f=0) ------------------
    wmax = max(abs(v) for v in w_main)
    met['max_abs_W_main'] = wmax
    chk['T2'] = (wmax <= TOL_W_ZERO,
                 'max|W_main|=%.4g (f=0, want <= %.1e)' % (wmax, TOL_W_ZERO))

    # --- T3: E is a real, decaying dissipation ----------------------------
    e_peak = max(e_main)
    tail = e_main[1:]
    monotone = all(b < a for a, b in zip(tail, tail[1:]))
    ratio = (e_main[1] / e_main[-1]) if e_main[-1] > 0.0 else float('inf')
    met.update(E_peak=e_peak, E_first=e_main[1], E_last=e_main[-1],
               E_decay_ratio=ratio)
    ok = e_peak > E_PEAK_MIN and monotone and ratio >= DECAY_RATIO_MIN
    chk['T3'] = (ok, 'peak E=%.5g, E(1)/E(50)=%.3f (want >=%.1f), '
                 'strictly decreasing from n=1: %s'
                 % (e_peak, ratio, DECAY_RATIO_MIN, monotone))

    # --- T4: final E magnitude vs the analytic BE series ------------------
    er = (e_main[-1] / e_an[-1]) if e_an[-1] > 0.0 else float('inf')
    met.update(E_an_last=e_an[-1], E_over_analytic=er)
    chk['T4'] = (E_LAST_LO <= er <= E_LAST_HI,
                 'E(50)=%.6g, E_an(50)=%.6g, ratio=%.3f (want [%.2f,%.2f])'
                 % (e_main[-1], e_an[-1], er, E_LAST_LO, E_LAST_HI))

    # --- T5: source work is alive and decays with p -----------------------
    wmin = min(w_src)
    src_monotone = all(b < a for a, b in zip(w_src[1:], w_src[2:]))
    wr = (w_src[0] / w_src[-1]) if w_src[-1] > 0.0 else float('inf')
    met.update(W_src_first=w_src[0], W_src_last=w_src[-1],
               W_decay_ratio=wr)
    ok = wmin > 0.0 and src_monotone and wr >= W_DECAY_RATIO_MIN
    chk['T5'] = (ok, 'min W_src=%.6g (>0), W_src(0)/W_src(50)=%.3f '
                 '(want >=%.1f), strictly decreasing: %s'
                 % (wmin, wr, W_DECAY_RATIO_MIN, src_monotone))

    # --- T6/T7: source work follows the analytic int(f*p) ----------------
    cos_w = _cos(w_src[1:], w_an[1:])
    met['cos_W'] = cos_w
    chk['T6'] = (cos_w is not None and cos_w >= COS_W_MIN,
                 'cos(W_src, W_an) [n>=1]=%s (want >=%.2f)'
                 % ('n/a' if cos_w is None else '%.5f' % cos_w, COS_W_MIN))
    wl = (w_src[-1] / w_an[-1]) if w_an[-1] > 0.0 else float('inf')
    met.update(W_an_last=w_an[-1], W_over_analytic=wl)
    chk['T7'] = (W_LAST_LO <= wl <= W_LAST_HI,
                 'W_src(50)=%.6g, W_an(50)=%.6g, ratio=%.3f '
                 '(want [%.2f,%.2f])'
                 % (w_src[-1], w_an[-1], wl, W_LAST_LO, W_LAST_HI))

    # --- T8: source deck E stays alive / decaying / same order of magnitude
    esr = (e_src[-1] / e_main[-1]) if e_main[-1] > 0.0 else float('inf')
    src_e_ok = max(e_src) > E_PEAK_MIN and all(
        b < a for a, b in zip(e_src[1:], e_src[2:]))
    met.update(E_src_last=e_src[-1], E_src_over_main=esr)
    chk['T8'] = (E_SRC_RATIO_LO <= esr <= E_SRC_RATIO_HI and src_e_ok,
                 'E_src(50)/E_main(50)=%.3f (want [%.2f,%.2f]); '
                 'E_src alive+decaying=%s'
                 % (esr, E_SRC_RATIO_LO, E_SRC_RATIO_HI, src_e_ok))
    return chk, met, None


def evaluate(main_log, src_log):
    return evaluate_text(read_log_text(main_log), read_log_text(src_log))


def print_report(chk, met, err):
    print('=== Gate 4.5 -- Darcy energy diagnostics verdict ===')
    if err:
        print('  DATA ERROR: %s' % err)
    for name in sorted(chk):
        ok, detail = chk[name]
        print('  %-4s %-4s %s' % (name, 'PASS' if ok else 'FAIL', detail))
    if met:
        print('  --- metrics ---')
        for key in ('main_rows', 'source_rows', 'E_peak', 'E_first', 'E_last',
                    'E_an_last', 'E_over_analytic', 'E_decay_ratio',
                    'max_abs_W_main', 'W_src_first', 'W_src_last',
                    'W_an_last', 'W_over_analytic', 'W_decay_ratio',
                    'cos_W', 'E_src_last', 'E_src_over_main'):
            if key in met:
                print('    %-20s = %s' % (key, met[key]))
    return [k for k, v in chk.items() if not v[0]]


# --------------------------------------------------------------------------
# Offline selftest: analytic replay + negative controls (no PETSc, no CI)
# --------------------------------------------------------------------------
def _fmt_log(e_rows, w_rows):
    lines = []
    for n, (e, w) in enumerate(zip(e_rows, w_rows)):
        lines.append('Gate 4.5: step %d t %.10E  E %.10E  W %.10E'
                     % (n + 1, n * DT, e, w))
    return '\n'.join(lines) + '\n'


def selftest():
    nrec = 51
    e_an = [analytic(n)[0] for n in range(nrec)]
    w_an = [analytic(n)[1] for n in range(nrec)]
    # n=0 continuum E diverges (step initial condition); the FE field is
    # finite. Use a representative FE-like value; T3 only checks n>=1.
    e_fe = [250.0] + [1.004 * v for v in e_an[1:]]
    e_src = [250.0] + [1.012 * v for v in e_an[1:]]

    main_ok = _fmt_log(e_fe, [0.0] * nrec)
    src_ok = _fmt_log(e_src, [1.004 * v for v in w_an])
    cases = [
        ('positive (wired G-F energy)', main_ok, src_ok, True),
        ('E all zero (old empty sections)',
         _fmt_log([0.0] * nrec, [0.0] * nrec), src_ok, False),
        ('W_src all zero (sourceSec unwired)',
         main_ok, _fmt_log(e_src, [0.0] * nrec), False),
        ('W_src reversed sign', main_ok,
         _fmt_log(e_src, [-1.004 * v for v in w_an]), False),
        ('E increasing (wrong sign / blow-up)',
         _fmt_log(list(reversed(e_fe)), [0.0] * nrec), src_ok, False),
        ('W_main non-zero (source leaked)',
         _fmt_log(e_fe, [0.5] * nrec), src_ok, False),
        ('missing rows', '\n'.join(main_ok.splitlines()[:10]) + '\n',
         src_ok, False),
    ]
    print('=== Gate 4.5 offline selftest (analytic replay, no solver) ===')
    print('  E_an(50)=%.6g  W_an(0)=%.6g  W_an(50)=%.6g'
          % (e_an[-1], w_an[0], w_an[-1]))
    all_ok = True
    for label, mt, st, want in cases:
        chk, met, err = evaluate_text(mt, st)
        fails = [k for k, v in chk.items() if not v[0]]
        got = (not fails) and not err
        ok = got == want
        all_ok &= ok
        print('  %-38s -> %-4s (want %-4s) %s  failing=%s'
              % (label, 'PASS' if got else 'FAIL', 'PASS' if want else 'FAIL',
                 'OK' if ok else 'SELFTEST ERROR', ','.join(fails) or 'none'))
    print('SELFTEST_RESULT: %s' % ('PASS' if all_ok else 'FAIL'))
    return 0 if all_ok else 1


def main(argv):
    if argv and argv[0] == '--selftest':
        return selftest()
    gate = bool(argv) and argv[0] == '--gate'
    args = argv[1:] if gate else argv
    usage = '<main.log> <source.log>'
    if len(args) != 2:
        if gate:
            print('GATE45_RESULT: FAIL (usage: --gate %s)' % usage)
            return 1
        print('usage: check_gate45_energy.py [--gate|--selftest] %s' % usage)
        return 0
    chk, met, err = evaluate(args[0], args[1])
    fails = print_report(chk, met, err)
    if gate:
        if err or fails or len(chk) < 8:
            print('GATE45_RESULT: FAIL (%s)'
                  % (err if err else ','.join(fails) or 'incomplete criteria'))
            return 1
        print('GATE45_RESULT: PASS (T1 data, T2 W_main=0, T3 E nonzero+decaying '
              '[E(1)/E(50)=%.1f], T4 E(50)/E_an=%.3f, T5 W_src>0+decaying, '
              'T6 cos_W=%s, T7 W_src/W_an=%.3f, T8 E_src/E_main=%.3f)'
              % (met.get('E_decay_ratio', 0.0), met.get('E_over_analytic', 0.0),
                 '%.4f' % met['cos_W'] if met.get('cos_W') is not None else 'n/a',
                 met.get('W_over_analytic', 0.0),
                 met.get('E_src_over_main', 0.0)))
        return 0
    print('VERDICT: %s (report mode, informational only)'
          % ('FAIL ' + ','.join(fails) if fails else 'PASS'))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
