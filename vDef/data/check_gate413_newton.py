#!/usr/bin/env python3
"""Gate 4.13 -- standard-Newton residual plateau: root cause frozen in CI.

Pure standard library. Run with ``--gate`` for the hard verdict
(``GATE413_RESULT: PASS|FAIL``, exit 0/1), ``--selftest`` for the offline
positive + negative replay, or with two paths for a report-only run.

ROOT CAUSE (static, backed by the CI #83 replay recorded in the checklist)
--------------------------------------------------------------------------
On CI #83 the displacement SNES ran standard Newton: the first iteration cut
|F| from 2.082e4 to ~9, then 51 iterations oscillated in an irreducible
3..13 plateau while ksponly produced a PHYSICALLY EQUIVALENT field in 6 s
(dPressure=0, dDamage=0.38%, dUx=0.21%, stress<=0.41%). So this is NOT an
assembly bug: the displacement system is effectively LINEAR for this deck
(SplitNone, damage~3.7e-5, zero damping/cohesion) and one linear solve is the
exact solve.

The plateau has a concrete scale mechanism in the source:

  residual, free dofs : integral(C:eps(u)) : eps(v) dOmega   -> scale E ~ 1e9
  residual, BC dofs   : residualPtr = u - u_bc               -> scale 1
  Jacobian,  BC dofs  : MatZeroRowsColumnsIS(A, IS, 1.0_Kr)  -> diagonal 1

The Jacobian therefore mixes rows of magnitude 1 and ~1e9 (condition number
~E). The unpreconditioned CG solve reaches its RELATIVE tolerance while the
physical force residual still sits at a few Newton, and SNES keeps iterating
on a fnorm that cannot decrease. CI #84 made ksponly the production solver,
which removes the futile loop without changing the solution.

This gate does not "prove" the mechanism dynamically (standard Newton is no
longer run in production); it FALSIFIABLY freezes the diagnosis: the scale
anchors must still be in the source, production must stay on ksponly, and the
ksponly field must remain bitwise identical to the production field.

  Q1  production solver is ksponly (never reverted to standard Newton)
  Q2  static scale anchors: unit-diagonal BC elimination + unit-scale BC
      residual + elastic-scale interior residual
  Q3  production log has NO Newton fnorm trace / displacement -5
  Q4  ksponly log has NO Newton fnorm trace, 51 Gate 3.2 rows
  Q5  ksponly and production pressure envelopes are bitwise identical
"""

import os
import re
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
SRC_ASM = os.path.join(_HERE, '..', '..', 'm_DefMech',
                       'm_MEF90_DefMechAssembly.F90')
SRC_VDEF = os.path.join(_HERE, '..', 'vDef.F90')
SRC_CI = os.path.join(_HERE, '..', '..', '.github', 'workflows', 'build.yml')

MIN_RECORDS = 50
G32_RE = re.compile(
    r'Gate 3\.2:\s*step\s+(\d+)\s+pmin\s+([-+0-9.eEdD]+)\s+'
    r'pmax\s+([-+0-9.eEdD]+)', re.I)
STEP_RE = re.compile(r'Gate 3\.1: step', re.I)
FNORM_RE = re.compile(r'SNES Function norm', re.I)
DISP_FAIL_RE = re.compile(r'displacement SNESSolve failed', re.I)
DAMAGE_FAIL_RE = re.compile(r'damage field SNESSolve failed', re.I)
FATAL_RE = re.compile(
    r'Caught signal|SEGV|MPI_ABORT|PETSC ERROR|segmentation fault|'
    r'TSSolve failed|KSP_DIVERGED|fortran runtime error', re.I)

CHECK_TITLES = {
    'Q1': 'production solver is ksponly (CI #84 decision not reverted)',
    'Q2': 'static scale anchors (unit BC elimination vs elastic interior)',
    'Q3': 'production log has no Newton fnorm trace / displacement -5',
    'Q4': 'ksponly log clean, 51 Gate 3.2 rows',
    'Q5': 'ksponly and production pressure envelopes bitwise identical',
}


def read_text(path):
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as fh:
            return fh.read()
    except Exception:                                    # noqa: BLE001
        return ''


def num(s):
    return float(s.replace('D', 'E').replace('d', 'e'))


def load_sources():
    return dict(asm=read_text(SRC_ASM), vdef=read_text(SRC_VDEF),
                ci=read_text(SRC_CI))


def envelopes(log):
    return [(int(m.group(1)), num(m.group(2)), num(m.group(3)))
            for m in G32_RE.finditer(log or '')]


def pass_fail(flag, detail):
    return (bool(flag), detail)


def check_q1(sources):
    ci = sources['ci']
    issues = []
    if not ci:
        issues.append('build.yml unreadable')
    else:
        i = ci.find('-options_file_yaml data/terzaghi1D.yaml')
        if i < 0:
            issues.append('production mpirun not found')
        else:
            block = ci[i:i + 900]
            if '-disp_snes_type ksponly' not in block:
                issues.append('production run is not ksponly')
            if re.search(r'-disp_snes_type\s+(newton|newtonls)', block):
                issues.append('production run reverted to standard Newton')
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'production mpirun uses -disp_snes_type ksponly')


def check_q2(sources):
    a = sources['asm']
    issues = []
    if not a:
        issues.append('DefMechAssembly unreadable')
    if a.count('MatZeroRowsColumnsIS(A,setISdof,1.0_Kr') < 2:
        issues.append('unit-diagonal Dirichlet elimination anchors missing')
    if 'residualPtr = displacementPtr - boundaryDisplacementPtr' not in a:
        issues.append('unit-scale BC residual anchor missing')
    if 'stressGauss .DotP.' not in a:
        issues.append('elastic-scale interior residual anchor missing')
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'BC diag=1, BC residual=u-u_bc, interior=integral(sigma:eps)')


def check_q3(log):
    issues = []
    n = len(FNORM_RE.findall(log or ''))
    d = len(DISP_FAIL_RE.findall(log or ''))
    g = len(DAMAGE_FAIL_RE.findall(log or ''))
    if n:
        issues.append('%d Newton fnorm rows in production' % n)
    if d:
        issues.append('%d displacement -5 in production' % d)
    if g > 10:
        issues.append('%d damage -5 (informational cap 10)' % g)
    if FATAL_RE.search(log or ''):
        issues.append('fatal signature in production log')
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'fnorm=0, displacement -5=0, damage -5=%d (informational), no fatal' % g)


def check_q4(log, env):
    issues = []
    if len(FNORM_RE.findall(log or '')):
        issues.append('ksponly log has a Newton fnorm trace')
    if len(DISP_FAIL_RE.findall(log or '')):
        issues.append('ksponly log has displacement -5')
    if FATAL_RE.search(log or ''):
        issues.append('fatal signature in ksponly log')
    if len(env) < MIN_RECORDS:
        issues.append('Gate 3.2 rows=%d < %d' % (len(env), MIN_RECORDS))
    return pass_fail(not issues, '; '.join(issues) if issues else
                     '%d Gate 3.2 rows, no fnorm, no fatal' % len(env))


def check_q5(prod_log, ko_log, env_p, env_k):
    issues = []
    if len(env_p) != len(env_k):
        issues.append('row counts %d vs %d' % (len(env_p), len(env_k)))
    else:
        for a, b in zip(env_p, env_k):
            if a != b:
                issues.append('step %d envelope differs: %r vs %r' % (a[0], a, b))
                break
    n1 = len(STEP_RE.findall(prod_log or ''))
    n2 = len(STEP_RE.findall(ko_log or ''))
    if n1 != n2:
        issues.append('step markers %d vs %d' % (n1, n2))
    return pass_fail(not issues, '; '.join(issues) if issues else
                     '%d envelope rows bitwise identical, %d step markers'
                     % (len(env_p), n1))


def evaluate(sources, prod_log, ko_log):
    ep, ek = envelopes(prod_log), envelopes(ko_log)
    return {
        'Q1': check_q1(sources),
        'Q2': check_q2(sources),
        'Q3': check_q3(prod_log),
        'Q4': check_q4(ko_log, ek),
        'Q5': check_q5(prod_log, ko_log, ep, ek),
    }, dict(rows=len(ep), damage_info=len(DAMAGE_FAIL_RE.findall(prod_log or '')))


def print_report(chk, met):
    print('=== Gate 4.13 -- standard-Newton plateau root cause verdict ===')
    for key in sorted(chk):
        ok, msg = chk[key]
        print('  %s %-3s %-70s | %s'
              % (key, 'PASS' if ok else 'FAIL', CHECK_TITLES.get(key, ''), msg))
    print('  --- archived CI #83 evidence (not a live run) ---')
    print('    Newton f0=2.082e4 -> it1~9.3 -> plateau flast 2.667/6.583/12.57')
    print('    ksponly field equivalent: dP=0, dUx=0.21%, dDamage=0.38%')
    return [k for k in sorted(chk) if not chk[k][0]]


def verdict_line(fails, met):
    if fails:
        return 'GATE413_RESULT: FAIL (%s)' % ','.join(fails)
    return ('GATE413_RESULT: PASS (Q1-Q5: scale-induced irreducible fnorm '
            'plateau diagnosed, production stays ksponly, %d envelopes '
            'bitwise identical)' % met['rows'])


def make_log(env, fnorm=False, disp_fail=False, fatal=False, steps=408):
    out = []
    for step, pmin, pmax in env:
        out.append('Gate 3.2: step %d pmin %.17g pmax %.17g' % (step, pmin, pmax))
    out.extend(['Gate 3.1: step %d enter' % (i % 51 + 1) for i in range(steps)])
    for _ in range(7):
        out.append('damage field SNESSolve failed with SNESConvergedReason -5')
    if fnorm:
        out.append('  0 SNES Function norm 9.345')
    if disp_fail:
        out.append('displacement SNESSolve failed with SNESConvergedReason -5')
    if fatal:
        out.append('PETSC ERROR: synthetic fatal signature')
    return '\n'.join(out) + '\n'


def selftest():
    print('=== Gate 4.13 offline selftest (synthetic replay, no solver) ===')
    src = load_sources()
    if not all(src.values()):
        print('SELFTEST_RESULT: FAIL (source files unreadable)')
        return 1
    env = [(i + 1, 0.0, 1.0e5 - 1800.0 * i) for i in range(51)]
    good_p, good_k = make_log(env), make_log(env)
    cases = []

    def run(p, k, s=None):
        chk, met = evaluate(s if s is not None else src, p, k)
        return set(x for x, (ok, _) in chk.items() if not ok), met

    fails, met = run(good_p, good_k)
    cases.append(('positive: ksponly production, scale anchors, identical fields',
                  fails, set(), met))

    f, _ = run(make_log(env, fnorm=True, disp_fail=True), good_k)
    cases.append(('NEG production reverted to standard Newton', f, {'Q3'}, None))

    f, _ = run(good_p, make_log(env, fnorm=True))
    cases.append(('NEG ksponly log has a Newton trace', f, {'Q4'}, None))

    bad = dict(src)
    bad['asm'] = src['asm'].replace('MatZeroRowsColumnsIS(A,setISdof,1.0_Kr',
                                    'MatZeroRowsColumnsIS(A,setISdof,2.0_Kr')
    f, _ = run(good_p, good_k, bad)
    cases.append(('NEG unit-diagonal BC elimination removed',
                  f, {'Q2'}, None))

    bad2 = dict(src)
    i = bad2['ci'].find('-options_file_yaml data/terzaghi1D.yaml')
    bad2['ci'] = (bad2['ci'][:i] + '-options_file_yaml data/terzaghi1D.yaml '
                  '-disp_snes_type newtonls ' + bad2['ci'][i + len('-options_file_yaml data/terzaghi1D.yaml'):])
    f, _ = run(good_p, good_k, bad2)
    cases.append(('NEG build.yml reverts production to Newton', f, {'Q1'}, None))

    pert = [(s, pmin, pmax - 1.0) if s == 25 else (s, pmin, pmax)
            for s, pmin, pmax in env]
    f, _ = run(good_p, make_log(pert))
    cases.append(('NEG ksponly field differs from production', f, {'Q5'}, None))

    f, _ = run(make_log(env, fatal=True), good_k)
    cases.append(('NEG fatal signature in production', f, {'Q3'}, None))

    ok = True
    for label, got, want, m in cases:
        good = (got == want)
        ok &= good
        print('  %-58s -> %-12s failing=%s (want %s) %s'
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
            print('GATE413_RESULT: FAIL (usage: --gate <production.log> '
                  '<ksponly.log>)')
            return 1
        print('usage: check_gate413_newton.py [--gate|--selftest] '
              '<production.log> <ksponly.log>')
        return 2
    try:
        p, k = read_text(args[0]), read_text(args[1])
        if not p or not k:
            raise ValueError('a solver log is unreadable')
        chk, met = evaluate(load_sources(), p, k)
    except Exception as exc:                          # noqa: BLE001
        print('GATE413_RESULT: FAIL (%s: %s)' % (type(exc).__name__, exc))
        return 1
    fails = print_report(chk, met)
    print(verdict_line(fails, met))
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
