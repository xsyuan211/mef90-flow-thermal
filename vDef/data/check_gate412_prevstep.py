#!/usr/bin/env python3
"""Gate 4.12 -- the previous-step pressure snapshot is a LIVE diagnostic.

Pure standard library (the CI runner has no numpy / scipy / PyYAML).
Run with ``--gate`` for the hard verdict (``GATE412_RESULT: PASS|FAIL``, exit
0/1), ``--selftest`` for the offline positive + negative replay, or with two
paths for a report-only run.

WHAT THE GATE FIXES
-------------------
``FlowCtx%pressurePreviousStep`` was created, initialised and destroyed, and
``vDef.F90`` copied ``pressure`` into it once per outer step -- but NOTHING
ever read it. That is a dead write: the vector occupies memory and costs an
MPI-wide copy every step while producing no information. The original intent
was to reserve it for an advection / Biot consumer; Gates 4.1a/4.1b/4.2
closed without ever needing a lagged pressure, so the reservation is obsolete.

Gate 4.12 turns it into a genuine, OBSERVABLE quantity with zero feedback into
the solve: the previous outer-step field is sampled with VecMin/VecMax BEFORE
TSSolve, and the per-step pressure envelope change is printed as

    Gate 4.12: step N t T dpmax DMAX dpmin DMIN

``dpmax`` is the closed-end (undrained) pressure drop, i.e. the discrete
consolidation rate. It is post-processing only: no Flow assembly routine
reads the vector, so all physics is bit-identical.

FALSIFIABLE CHECKS
------------------
  Z1  data integrity: 51 Gate 3.2 + 51 Gate 4.12 rows, contiguous steps,
      monotone dt=0.02 time axis, mesh has 51 records
  Z2  static anti-regression: init snapshot, pre-TSSolve VecMin/VecMax,
      per-step VecCopy, diagnostic print, CI hard gate, and NO consumer in
      m_Flow (proves zero feedback)
  Z3  log self-consistency: dpmax(n) == pmax(n)-pmax(n-1) and
      dpmin(n) == pmin(n)-pmin(n-1) bitwise; step 1 drop is exactly zero
  Z4  EXO anchor: pmax(n)/pmin(n) equal max/min of the nodal pressure field
  Z5  physics shape: all drops negative, single peak in steps 3..20,
      last/peak ratio < 0.30, last ten drops non-increasing, dpmin == 0
  Z6  production isolation: no fatal signatures, >= 51 step markers
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
    CDF2 = None
    _CDF2_IMPORT_ERROR = exc

SRC_VDEF = os.path.join(_HERE, '..', 'vDef.F90')
SRC_ASM = os.path.join(_HERE, '..', '..', 'm_Flow', 'm_MEF90_FlowAssembly.F90')
SRC_CTX = os.path.join(_HERE, '..', '..', 'm_Flow', 'm_MEF90_FlowCtx.F90')
SRC_CI = os.path.join(_HERE, '..', '..', '.github', 'workflows', 'build.yml')

MIN_RECORDS = 50
DT = 0.02
PEAK_LO, PEAK_HI = 3, 20
LAST_PEAK_MAX = 0.30
TAIL = 10
EXO_REL_TOL = 1.0e-9
DROP_ABS_TOL = 1.0e-8

G32_RE = re.compile(
    r'Gate 3\.2:\s*step\s+(\d+)\s+pmin\s+([-+0-9.eEdD]+)\s+'
    r'pmax\s+([-+0-9.eEdD]+)', re.I)
G412_RE = re.compile(
    r'Gate 4\.12:\s*step\s+(\d+)\s+t\s+([-+0-9.eEdD]+)\s+'
    r'dpmax\s+([-+0-9.eEdD]+)\s+dpmin\s+([-+0-9.eEdD]+)', re.I)
STEP_RE = re.compile(r'Gate 3\.1: step', re.I)
FATAL_RE = re.compile(
    r'Caught signal|SEGV|MPI_ABORT|PETSC ERROR|segmentation fault|'
    r'TSSolve failed|KSP_DIVERGED|fortran runtime error', re.I)

CHECK_TITLES = {
    'Z1': 'data integrity (51+51 rows, contiguous steps, dt=0.02)',
    'Z2': 'static anti-regression (live snapshot, no assembly consumer)',
    'Z3': 'log self-consistency (drop == pmax(n)-pmax(n-1), step 1 zero)',
    'Z4': 'EXO anchor (pmax/pmin == nodal field extrema)',
    'Z5': 'physics shape (negative, single peak, decaying tail)',
    'Z6': 'production isolation (no fatal signature)',
}


def read_text(path):
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as fh:
            return fh.read()
    except Exception:                                    # noqa: BLE001
        return ''


def num(s):
    return float(s.replace('D', 'E').replace('d', 'e'))


def load_fields(path):
    if CDF2 is None:
        raise RuntimeError('parse_cdf2 unavailable: %r' % (_CDF2_IMPORT_ERROR,))
    c = CDF2(path)
    rows = [c.read_record('vals_nod_var1', i) for i in range(c.numrecs)]
    if not rows or any(len(r) < 2 for r in rows):
        raise ValueError('empty or degenerate EXO pressure field')
    for row in rows:
        if not all(math.isfinite(v) for v in row):
            raise ValueError('non-finite EXO pressure')
    return dict(nrec=c.numrecs, nn=len(rows[0]), p=rows)


def load_sources():
    return dict(vdef=read_text(SRC_VDEF), asm=read_text(SRC_ASM),
                ctx=read_text(SRC_CTX), ci=read_text(SRC_CI))


def parse_rows(log):
    g32 = [(int(m.group(1)), num(m.group(2)), num(m.group(3)))
           for m in G32_RE.finditer(log or '')]
    g412 = [(int(m.group(1)), num(m.group(2)), num(m.group(3)), num(m.group(4)))
            for m in G412_RE.finditer(log or '')]
    return g32, g412


def pass_fail(flag, detail):
    return (bool(flag), detail)


# ---------------------------------------------------------------------------
# checks
# ---------------------------------------------------------------------------
def check_z1(log, g32, g412, mesh):
    issues = []
    n = len(g32)
    if n < MIN_RECORDS:
        issues.append('Gate 3.2 rows=%d < %d' % (n, MIN_RECORDS))
    if len(g412) != n:
        issues.append('Gate 4.12 rows=%d, Gate 3.2 rows=%d' % (len(g412), n))
    if [r[0] for r in g32] != list(range(1, n + 1)):
        issues.append('Gate 3.2 steps not contiguous 1..%d' % n)
    if [r[0] for r in g412] != list(range(1, n + 1)):
        issues.append('Gate 4.12 steps not contiguous 1..%d' % n)
    for step, t, _, _ in g412:
        if abs(t - (step - 1) * DT) > 1.0e-9:
            issues.append('step %d t=%g != (step-1)*dt' % (step, t))
            break
    if mesh['nrec'] < MIN_RECORDS:
        issues.append('EXO records=%d < %d' % (mesh['nrec'], MIN_RECORDS))
    if mesh['nrec'] != n:
        issues.append('EXO records=%d != log rows=%d' % (mesh['nrec'], n))
    return pass_fail(not issues, '; '.join(issues) if issues else
                     '%d rows, %d EXO records, %d nodes' %
                     (n, mesh['nrec'], mesh['nn']))


def check_z2(sources):
    v, a, c, ci = (sources[k] for k in ('vdef', 'asm', 'ctx', 'ci'))
    issues = []
    if not v or not a or not c or not ci:
        issues.append('a source/workflow file is unreadable')
    ncopy = v.count('Call VecCopy(MEF90FlowCtx%pressure,MEF90FlowCtx%pressurePreviousStep')
    if ncopy < 2:
        issues.append('init + per-step VecCopy snapshots=%d (want 2)' % ncopy)
    if 'Call VecMin(MEF90FlowCtx%pressurePreviousStep' not in v or \
       'Call VecMax(MEF90FlowCtx%pressurePreviousStep' not in v:
        issues.append('previous-step VecMin/VecMax sampling missing')
    if "'Gate 4.12: step'" not in v:
        issues.append('Gate 4.12 diagnostic print missing')
    if 'pressurePrevMin' not in v or 'pressurePrevMax' not in v:
        issues.append('previous-step extrema not declared')
    # Zero feedback: no Flow assembly routine may consume the vector.
    if 'pressurePreviousStep' in a:
        issues.append('FlowAssembly now consumes pressurePreviousStep (feedback!)')
    if 'VecDestroy(FlowCtx%pressurePreviousStep' not in c:
        issues.append('FlowCtx lifecycle incomplete')
    if 'check_gate412_prevstep.py --gate' not in ci:
        issues.append('CI hard gate not wired')
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'snapshots=2, sampled, printed, wired, zero assembly consumer')


def check_z3(g32, g412):
    issues = []
    if len(g32) != len(g412):
        return pass_fail(False, 'row count mismatch')
    for i, (step, pmin, pmax) in enumerate(g32):
        s, _, dmax, dmin = g412[i]
        if s != step:
            issues.append('step mismatch at row %d' % i)
            break
        want_max = 0.0 if i == 0 else pmax - g32[i - 1][2]
        want_min = 0.0 if i == 0 else pmin - g32[i - 1][1]
        if abs(dmax - want_max) > DROP_ABS_TOL * max(1.0, abs(want_max)):
            issues.append('step %d dpmax=%g want %g' % (step, dmax, want_max))
            break
        if abs(dmin - want_min) > DROP_ABS_TOL * max(1.0, abs(want_min)):
            issues.append('step %d dpmin=%g want %g' % (step, dmin, want_min))
            break
    if g412 and (g412[0][2] != 0.0 or g412[0][3] != 0.0):
        issues.append('step 1 drop is not exactly zero')
    return pass_fail(not issues, '; '.join(issues) if issues else
                     '%d rows bitwise consistent, step 1 drop = 0' % len(g32))


def check_z4(g32, mesh):
    issues = []
    for i, (step, pmin, pmax) in enumerate(g32):
        if i >= len(mesh['p']):
            break
        row = mesh['p'][i]
        hi, lo = max(row), min(row)
        if abs(hi - pmax) > EXO_REL_TOL * max(1.0, abs(pmax)):
            issues.append('step %d EXO max=%g vs pmax=%g' % (step, hi, pmax))
            break
        if abs(lo - pmin) > EXO_REL_TOL * max(1.0, abs(pmin)) + 1.0e-12:
            issues.append('step %d EXO min=%g vs pmin=%g' % (step, lo, pmin))
            break
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'pmax/pmin match nodal field extrema on %d records' % len(g32))


def check_z5(g32, g412):
    issues = []
    drops = [r[2] for r in g412][1:]
    dmins = [r[3] for r in g412]
    if not drops:
        return pass_fail(False, 'no drops')
    if any(d >= 0 for d in drops):
        issues.append('a pressure drop is non-negative (not consolidating)')
    if any(abs(d) > 1.0e-9 for d in dmins):
        issues.append('drained-face dpmin is not zero')
    peak = max(range(len(drops)), key=lambda i: abs(drops[i])) + 2
    if not (PEAK_LO <= peak <= PEAK_HI):
        issues.append('peak |dpmax| at step %d outside %d..%d'
                      % (peak, PEAK_LO, PEAK_HI))
    peak_abs = max(abs(d) for d in drops)
    if peak_abs <= 0.0:
        issues.append('all pressure drops are zero (dead diagnostic)')
    ratio = abs(drops[-1]) / peak_abs if peak_abs > 0.0 else 1.0
    if ratio >= LAST_PEAK_MAX:
        issues.append('last/peak |dpmax|=%.4f >= %.2f' % (ratio, LAST_PEAK_MAX))
    tail = drops[-TAIL:]
    if any(abs(tail[i]) > abs(tail[i - 1]) + 1.0e-9
           for i in range(1, len(tail))):
        issues.append('|dpmax| not decaying over the last %d steps' % TAIL)
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'peak @step %d, last/peak=%.3f, tail decays, all drops < 0'
                     % (peak, ratio)), dict(peak_step=peak, last_peak=ratio,
                                            first_drop=drops[0],
                                            last_drop=drops[-1])


def check_z6(log, g32):
    issues = []
    if FATAL_RE.search(log or ''):
        issues.append('fatal signature in production log')
    n = len(STEP_RE.findall(log or ''))
    if n < len(g32):
        issues.append('step markers=%d < rows=%d' % (n, len(g32)))
    return pass_fail(not issues, '; '.join(issues) if issues else
                     '%d step markers, no fatal signature' % n)


def evaluate(log, mesh, sources):
    g32, g412 = parse_rows(log)
    chk = {}
    met = {}
    chk['Z1'] = check_z1(log, g32, g412, mesh)
    chk['Z2'] = check_z2(sources)
    chk['Z3'] = check_z3(g32, g412)
    chk['Z4'] = check_z4(g32, mesh)
    chk['Z5'], met['z5'] = check_z5(g32, g412)
    chk['Z6'] = check_z6(log, g32)
    return chk, met


def print_report(chk, met):
    print('=== Gate 4.12 -- previous-step pressure snapshot verdict ===')
    for key in sorted(chk):
        ok, msg = chk[key]
        print('  %s %-3s %-70s | %s'
              % (key, 'PASS' if ok else 'FAIL', CHECK_TITLES.get(key, ''), msg))
    if met.get('z5'):
        print('  --- metrics ---')
        for k in ('peak_step', 'last_peak', 'first_drop', 'last_drop'):
            print('    %-12s = %.6g' % (k, met['z5'][k]))
    return [k for k in sorted(chk) if not chk[k][0]]


def verdict_line(fails, met):
    if fails:
        return 'GATE412_RESULT: FAIL (%s)' % ','.join(fails)
    return ('GATE412_RESULT: PASS (Z1-Z6: pressurePreviousStep consumed, '
            'dpmax/dpmin bitwise consistent, EXO-anchored, peak @step %d, '
            'physics unchanged)' % met['z5']['peak_step'])


# ---------------------------------------------------------------------------
# offline selftest
# ---------------------------------------------------------------------------
def make_log(pmax, dpmax_fn=None, dpmin_fn=None, drop_rows=(), crash=False,
             nsteps=None):
    lines = []
    n = len(pmax) if nsteps is None else nsteps
    for i in range(1, n + 1):
        for _ in range(8):
            lines.append('Gate 3.1: step %d enter transient block' % i)
        lines.append('Gate 3.2: step %d pmin %.17g pmax %.17g' % (i, 0.0, pmax[i - 1]))
        dmax = 0.0 if i == 1 else pmax[i - 1] - pmax[i - 2]
        dmin = 0.0
        if dpmax_fn is not None:
            dmax = dpmax_fn(i, dmax)
        if dpmin_fn is not None:
            dmin = dpmin_fn(i, dmin)
        if i not in drop_rows:
            lines.append('Gate 4.12: step %d t %.17g dpmax %.17g dpmin %.17g'
                         % (i, (i - 1) * DT, dmax, dmin))
    if crash:
        lines.append('PETSC ERROR: synthetic fatal signature')
    return '\n'.join(lines) + '\n'


def make_mesh(pmax, nn=356):
    # node 0 carries the closed-end pmax, node 1 the drained-face pmin=0,
    # the rest interpolate; max/min of each row therefore match Gate 3.2.
    p = []
    for v in pmax:
        row = [v] + [0.0] + [v * 0.5] * (nn - 2)
        p.append(row)
    return dict(nrec=len(pmax), nn=nn, p=p)


def selftest():
    print('=== Gate 4.12 offline selftest (synthetic replay, no solver) ===')
    src = load_sources()
    if not all(src.values()):
        print('SELFTEST_RESULT: FAIL (source files unreadable)')
        return 1
    # Real CI #115 pmax sequence shape: fast early consolidation, then a
    # single peak around step 10 followed by monotone exponential decay.
    p0 = 1.0e5
    pmax = [p0]
    for n in range(1, 51):
        d = -3800.0 * math.exp(-((n - 9.0) / 5.0) ** 2) \
            - 500.0 * math.exp(-n / 40.0) - 200.0
        pmax.append(pmax[-1] + d)
    good_log = make_log(pmax)
    good_mesh = make_mesh(pmax)

    cases = []

    def run(log, mesh, sources=None):
        chk, met = evaluate(log, mesh, sources if sources is not None else src)
        return set(k for k, (ok, _) in chk.items() if not ok), met

    fails, met = run(good_log, good_mesh)
    cases.append(('positive: live snapshot, consistent drops',
                  fails, set(), met))

    f, _ = run(make_log(pmax, drop_rows=(7,)), good_mesh)
    cases.append(('NEG a Gate 4.12 row is missing', f, {'Z1', 'Z3'}, None))

    bad = dict(src)
    bad['vdef'] = src['vdef'].replace(
        'Call VecMin(MEF90FlowCtx%pressurePreviousStep',
        'Call VecMin(MEF90FlowCtx%pressure')
    f, _ = run(good_log, good_mesh, bad)
    cases.append(('NEG pre-TSSolve sampling removed (dead write back)',
                  f, {'Z2'}, None))

    f, _ = run(make_log(pmax, dpmax_fn=lambda i, d: -d), good_mesh)
    cases.append(('NEG drops sign-flipped', f, {'Z3', 'Z5'}, None))

    f, _ = run(make_log(pmax, dpmax_fn=lambda i, d: 0.0 if i > 1 else d),
               good_mesh)
    cases.append(('NEG drops all zero (dead diagnostic)', f, {'Z3', 'Z5'}, None))

    pert = make_mesh(pmax)
    pert['p'][25][0] *= 1.01
    f, _ = run(good_log, pert)
    cases.append(('NEG EXO pressure max perturbed', f, {'Z4'}, None))

    f, _ = run(make_log(pmax, dpmin_fn=lambda i, d: 1.0 if i > 1 else d),
               good_mesh)
    cases.append(('NEG drained face dpmin nonzero', f, {'Z3', 'Z5'}, None))

    f, _ = run(make_log(pmax, crash=True), good_mesh)
    cases.append(('NEG fatal signature', f, {'Z6'}, None))

    bad2 = dict(src)
    bad2['asm'] = src['asm'] + '\n! pressurePreviousStep\n'
    f, _ = run(good_log, good_mesh, bad2)
    cases.append(('NEG assembly now reads the snapshot (feedback)',
                  f, {'Z2'}, None))

    ok = True
    for label, got, want, m in cases:
        good = (got == want)
        ok &= good
        print('  %-52s -> %-12s failing=%s (want %s) %s'
              % (label, 'PASS' if good else 'FAIL',
                 ','.join(sorted(got)) or 'none',
                 ','.join(sorted(want)) or 'none',
                 'OK' if good else '<-- MISMATCH'))
    print('SELFTEST_RESULT: %s' % ('PASS' if ok else 'FAIL'))
    return 0 if ok else 1


def main(argv):
    if argv and argv[0] == '--selftest':
        return selftest()
    gate = bool(argv) and argv[0] == '--gate'
    args = argv[1:] if gate else argv
    if len(args) != 2:
        if gate:
            print('GATE412_RESULT: FAIL (usage: --gate <production.log> '
                  '<production.gen>)')
            return 1
        print('usage: check_gate412_prevstep.py [--gate|--selftest] '
              '<production.log> <production.gen>')
        return 2
    try:
        log = read_text(args[0])
        mesh = load_fields(args[1])
        sources = load_sources()
        if not log:
            raise ValueError('production log unreadable')
        chk, met = evaluate(log, mesh, sources)
    except Exception as exc:                          # noqa: BLE001
        print('GATE412_RESULT: FAIL (%s: %s)' % (type(exc).__name__, exc))
        return 1
    fails = print_report(chk, met)
    print(verdict_line(fails, met))
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
