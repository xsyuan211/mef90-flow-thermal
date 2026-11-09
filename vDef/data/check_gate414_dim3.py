#!/usr/bin/env python3
"""Gate 4.14 -- hard verdict for the 3D dimension-dispatch diagnostic.

Why this gate exists
--------------------
Every one of the 15 existing verdicts runs on ``Beam1x10-tri3`` (num_dim=2).
MEF90 nevertheless compiles BOTH dimension twins: m_Flow/Makefile and
m_DefMech/Makefile build -DMEF90_DIM=2 and -DMEF90_DIM=3 objects, and
m_MEF90_Flow.F90 dispatches at run time on DMMeshGetDimension. In particular
MEF90/mef90.inc maps MEF90_MATS_IDENTITY to MEF90MatS2DIdentity (3-component
symmetric tensor) in 2D and MEF90MatS3DIdentity (6-component) in 3D -- the
same macro used by Gate 4.1b's effective-stress term and Gate 4.9's batched
diffusion operator. A wrong 3D identity, a missing dispatch or a 2D-only
assembly bug would never have been caught by any previous gate.

The 3D mesh is already committed: Beam1x10-tetra4.gen, num_dim=3, 2635 nodes,
11761 TETRA4 cells, cell set labels [1,10,20] identical to the 2D mesh, same
10 x 1 x 1 m beam. terzaghi1D_dim3.yaml differs from the control deck by ONE
code line (time.numstep 51 -> 11, dt=0.1 s, same t=0..1 s span); the dimension
is selected by -prefix, not by a YAML key.

The 3D diagnostic is the only deck that runs on one MPI rank with a direct
block-LU preconditioner (CI #120): on two ranks the CG + bjacobi/ILU(0) chain
never finished the FIRST TSSolve on the 2635-node tetra4 mesh within 900 s,
leaving a truncated EXO. Production 2D decks are unchanged. W6 now also
requires one "Linear solve converged" row per advanced time step, so a
3D run that only initialises can never pass on a pre-allocated EXO file.

W1 run/mesh data integrity
W2 3D dispatch evidence (Displacement_Z, num_dim=3, coordz, same cell sets)
W3 static anti-regression (3D identity macro, 3D twins, runtime dispatch)
W4 consolidation physics on the 11 shared physical times (3D vs 2D pmax)
W5 3D displacement field is live, finite and z-bounded
W6 3D solver chain completed (incl. per-step KSP convergence), no fatal signature
W7 one-line deck + CI wiring (incl. solver config) + production isolation

Run with --gate for the hard verdict (exit 0/1) or --selftest for an offline
positive + negative replay. Pure standard library; no local compile, no solver.
"""

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse_cdf2 import CDF2  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, '..', '..'))
DECK_3D = os.path.join(_HERE, 'terzaghi1D_dim3.yaml')
DECK_2D = os.path.join(_HERE, 'terzaghi1D.yaml')
SRC = {
    'inc': os.path.join(_ROOT, 'MEF90', 'mef90.inc'),
    'flow_facade': os.path.join(_ROOT, 'm_Flow', 'm_MEF90_Flow.F90'),
    'flow_mk': os.path.join(_ROOT, 'm_Flow', 'Makefile'),
    'dm_mk': os.path.join(_ROOT, 'm_DefMech', 'Makefile'),
    'flow_asm': os.path.join(_ROOT, 'm_Flow', 'm_MEF90_FlowAssembly.F90'),
    'dm_asm': os.path.join(_ROOT, 'm_DefMech', 'm_MEF90_DefMechAssembly.F90'),
    'dm_facade': os.path.join(_ROOT, 'm_DefMech', 'm_MEF90_DefMech.F90'),
    'workflow': os.path.join(_ROOT, '.github', 'workflows', 'build.yml'),
}
NAME_WIDTH = 33
N_REC_3D = 11
N_REC_2D = 51
FATAL = ('PETSC ERROR', 'Caught signal', 'SEGV', 'fortran runtime error',
         'Fatal Error:', 'NaN detected')
# The 3D deck uses dt=0.1 s while the 2D control uses dt=0.02 s. Backward
# Euler on a decaying mode with lambda~2.4 gives ~0.093 (51 steps) vs ~0.113
# (11 steps), i.e. a legitimate ~25% discretisation difference. Keep the
# tolerance above that, and rely on W4's monotonicity / p0 / pend checks plus
# W2/W3/W5 for the discriminating power.
MAX_REL_P = 0.40
MAX_UZ_RATIO = 50.0
MIN_UX = 1.0e-9


def pass_fail(ok, detail):
    return ('PASS' if ok else 'FAIL'), detail


def var_names(c):
    raw = b''.join(c.read_static('name_nod_var'))
    return [raw[i:i + NAME_WIDTH].rstrip(b'\x00').decode('ascii', 'replace')
            for i in range(0, len(raw), NAME_WIDTH)]


def mesh_meta(path):
    c = CDF2(path)
    d = dict(zip(c.dims, c.sizes))
    out = {'num_dim': d.get('num_dim'), 'num_nodes': d.get('num_nodes'),
           'num_elem': d.get('num_elem'), 'num_el_blk': d.get('num_el_blk'),
           'sets': list(c.read_static('eb_prop1')),
           'has_z': any(v['name'] == 'coordz' for v in c.vars)}
    for ax in ('x', 'y', 'z'):
        try:
            v = c.read_static('coord' + ax)
            out['r' + ax] = (min(v), max(v))
        except KeyError:
            out['r' + ax] = None
    return out


def exo_fields(path):
    c = CDF2(path)
    names = var_names(c)
    p, ux, uy, uz = [], [], [], []
    for i in range(c.numrecs):
        p.append(max(c.read_record('vals_nod_var1', i)))
        ux.append(max(abs(x) for x in c.read_record('vals_nod_var%d' % (names.index('Displacement_X') + 1), i)))
        uy.append(max(abs(x) for x in c.read_record('vals_nod_var%d' % (names.index('Displacement_Y') + 1), i)))
        if 'Displacement_Z' in names:
            uz.append(max(abs(x) for x in c.read_record('vals_nod_var%d' % (names.index('Displacement_Z') + 1), i)))
    return {'numrecs': c.numrecs, 'names': names, 'pmax': p, 'ux': ux,
            'uy': uy, 'uz': uz, 'finite': all(
                math.isfinite(x) for seq in (p, ux, uy, uz) for x in seq)}


def code_lines(text):
    out = []
    for line in text.splitlines():
        s = line.split('#', 1)[0].rstrip()
        if s.strip():
            out.append(s.strip())
    return out


def read(path):
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        return f.read()


def check_w1(d):
    issues = []
    if len(d['log3']) < 1000:
        issues.append('3D solver log too short')
    if d['exo3']['numrecs'] != N_REC_3D:
        issues.append('3D records %d != %d' % (d['exo3']['numrecs'], N_REC_3D))
    if d['exo2']['numrecs'] != N_REC_2D:
        issues.append('2D records %d != %d' % (d['exo2']['numrecs'], N_REC_2D))
    if not d['exo3']['finite'] or not d['exo2']['finite']:
        issues.append('non-finite EXO field')
    if not d['m3']['num_nodes'] or d['m3']['num_nodes'] <= d['m2']['num_nodes']:
        issues.append('3D mesh is not larger than the 2D mesh')
    for sig in FATAL:
        if sig in d['log3']:
            issues.append('fatal signature in 3D log: %s' % sig)
    return pass_fail(not issues, '; '.join(issues) if issues else
                     '%d 3D / %d 2D records, %d 3D nodes, all fields finite'
                     % (d['exo3']['numrecs'], d['exo2']['numrecs'],
                        d['m3']['num_nodes']))


def check_w2(d):
    issues = []
    n3, n2 = d['exo3']['names'], d['exo2']['names']
    if 'Displacement_Z' not in n3:
        issues.append('3D EXO has no Displacement_Z (2D-only dispatch)')
    if len(n3) < 5:
        issues.append('3D nodal table has only %d vars, expected >=5' % len(n3))
    if 'Displacement_Z' in n2 or len(n2) != 4:
        issues.append('2D nodal table is not the 4-var control')
    if d['m3']['num_dim'] != 3 or d['m2']['num_dim'] != 2:
        issues.append('mesh dimensions %s/%s != 3/2' % (d['m3']['num_dim'], d['m2']['num_dim']))
    if not d['m3']['has_z']:
        issues.append('3D mesh has no coordz')
    if d['m3']['sets'] != d['m2']['sets'] or d['m3']['sets'] != [1, 10, 20]:
        issues.append('cell set labels differ: %s vs %s' % (d['m3']['sets'], d['m2']['sets']))
    for ax in ('x', 'y'):
        a, b = d['m3']['r' + ax], d['m2']['r' + ax]
        if not a or not b or max(abs(a[i] - b[i]) for i in (0, 1)) > 1e-9:
            issues.append('%s geometry differs between 2D and 3D meshes' % ax)
    if not d['m3']['rz'] or abs(d['m3']['rz'][0] + 0.5) > 1e-9 or abs(d['m3']['rz'][1] - 0.5) > 1e-9:
        issues.append('3D z extent is not [-0.5,0.5]')
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'num_dim=3, %d nodal vars incl. Displacement_Z, cell sets [1,10,20], geometry 10x1x1'
                     % len(n3))


def check_w3(d):
    issues = []
    s = d['src']
    if 'MEF90MatS3DIdentity' not in s['inc'] or 'MEF90MatS2DIdentity' not in s['inc']:
        issues.append('mef90.inc lost a dimension identity macro')
    if 'MEF90_DIM' not in s['inc'] or 'MEF90_MATS_IDENTITY' not in s['inc']:
        issues.append('mef90.inc lost the MEF90_DIM conditional')
    if 'DMMeshGetDimension' not in s['flow_facade'] or 'dim == 3' not in s['flow_facade']:
        issues.append('runtime 3D dispatch missing from m_MEF90_Flow.F90')
    if '-DMEF90_DIM=3' not in s['flow_mk'] or '-DMEF90_DIM=3' not in s['dm_mk']:
        issues.append('a 3D object twin is no longer compiled')
    if 'MEF90_MATS_IDENTITY' not in s['flow_asm'] or 'MEF90_MATS_IDENTITY' not in s['dm_asm']:
        issues.append('assembly no longer uses the dimension-aware identity')
    if 'Displacement_Z' not in s['dm_facade']:
        issues.append('DefMech EXO naming lost its 3D branch')
    # CI #119: the Flow facade is compiled WITHOUT -DMEF90_DIM, so the
    # compile-time MEF90_ELEMENT_SCAL macro silently bound to the 2D element
    # type. MEF90Element_Create is a GENERIC resolved by the DECLARED type of
    # elem, so a 3D mesh hit "Element2D_Scal_Init: not implemented yet" in
    # MEF90FlowUpdatePerm. Freeze the 2D/3D implementation-pair fix here.
    ff = s['flow_facade']
    for token in ('MEF90FlowUpdatePerm2D', 'MEF90FlowUpdatePerm3D',
                  'MEF90FlowUpdateMuFromTemperature2D',
                  'MEF90FlowUpdateMuFromTemperature3D',
                  'Type(MEF90Element2D_Scal)', 'Type(MEF90Element3D_Scal)',
                  'Call MEF90FlowUpdatePerm2D(', 'Call MEF90FlowUpdatePerm3D(',
                  'Call MEF90FlowUpdateMuFromTemperature2D(',
                  'Call MEF90FlowUpdateMuFromTemperature3D('):
        if token not in ff:
            issues.append('flow facade lost %s' % token)
    if 'Type(MEF90_ELEMENT_SCAL)' in ff:
        issues.append('flow facade bound elem to the compile-time 2D macro again')
    return pass_fail(not issues, '; '.join(issues) if issues else
                     '3D/2D identity macros, -DMEF90_DIM=3 twins, runtime dispatch, '
                     '2D/3D facade update pair, 3D EXO naming')


def check_w4(d):
    issues = []
    p3, p2 = d['exo3']['pmax'], d['exo2']['pmax']
    if len(p3) < N_REC_3D or len(p2) < N_REC_2D:
        return pass_fail(False, 'not enough records for the shared-time comparison')
    if abs(p3[0] - 1.0e5) > 1.0e3:
        issues.append('3D initial pmax %.6g != 1e5' % p3[0])
    if p3[-1] >= p3[0] * 0.5:
        issues.append('3D pressure did not consolidate (end %.6g)' % p3[-1])
    if any(p3[i] >= p3[i - 1] + 1.0 for i in range(1, len(p3))):
        issues.append('3D pmax is not monotonically decaying')
    worst, at = 0.0, -1
    for i in range(N_REC_3D):
        a, b = p3[i], p2[5 * i]
        rel = abs(a - b) / max(abs(b), 1.0)
        if rel > worst:
            worst, at = rel, i
    if worst > MAX_REL_P:
        issues.append('3D vs 2D pmax rel-diff %.3f at shared time %d > %.2f'
                      % (worst, at, MAX_REL_P))
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'p0=%.4g, pend=%.4g, monotonic decay, max shared-time rel-diff %.3f'
                     % (p3[0], p3[-1], worst))


def check_w5(d):
    issues = []
    uz = d['exo3']['uz']
    ux = d['exo3']['ux']
    if not uz:
        issues.append('no Displacement_Z field in 3D output')
    elif max(ux) < MIN_UX:
        issues.append('3D displacement field is identically zero (DefMech 3D never ran)')
    else:
        ratio = max(uz) / max(max(ux), MIN_UX)
        if ratio > MAX_UZ_RATIO:
            issues.append('|uz|/|ux| = %.3g > %.0f' % (ratio, MAX_UZ_RATIO))
    if not d['exo3']['finite']:
        issues.append('non-finite 3D displacement')
    return pass_fail(not issues, '; '.join(issues) if issues else
                     '|ux|max=%.4g, |uz|max=%.4g, ratio=%.3g, all finite'
                     % (max(ux), max(uz) if uz else 0.0,
                        (max(uz) / max(ux)) if uz and max(ux) else 0.0))


def check_w6(d):
    issues = []
    log = d['log3']
    for tag in ('Gate 3.1: step', 'Gate 3.2: step', 'Gate 4.12: step', 'Gate 4.1:'):
        n = log.count(tag)
        if n < 5:
            issues.append('%s appears only %d times' % (tag, n))
    # CI #120: the 3D run reached step 2 but its FIRST flow TSSolve never
    # printed a convergence line; the timeout killed MPI and left a
    # truncated EXO. Diagnostic rows alone are not proof that the 3D
    # assembly/solve path completed. Step 1 skips TSSolve ("no time
    # advance"), so N_REC_3D-1 completed linear solves are required.
    n_solves = log.count('Linear solve converged')
    if n_solves < N_REC_3D - 1:
        issues.append('only %d completed flow linear solves (expected >= %d)'
                      % (n_solves, N_REC_3D - 1))
    if 'Gate 4.9:' in log:
        issues.append('3D deck unexpectedly enabled the batched fast path')
    if 'Gate 4.1b:' in log:
        issues.append('3D deck unexpectedly enabled the effective-stress mode')
    for sig in FATAL:
        if sig in log:
            issues.append('fatal signature: %s' % sig)
    return pass_fail(not issues, '; '.join(issues) if issues else
                     '3.1/3.2/4.12/4.1 chains + %d converged flow solves, '
                     'no fatal signature, default physics mode' % n_solves)


def check_w7(d):
    issues = []
    a, b = code_lines(d['deck2']), code_lines(d['deck3'])
    diff = [(x, y) for x, y in zip(a, b) if x != y]
    if len(a) != len(b) or len(diff) != 1:
        issues.append('deck code-line diff is %d (expected 1)' %
                      (len(diff) if len(a) == len(b) else -1))
    else:
        x, y = diff[0]
        if not (x.startswith('numstep:') and y.startswith('numstep:')):
            issues.append('the single deck diff is not time.numstep')
        elif x.split(':', 1)[1].strip() != '51' or y.split(':', 1)[1].strip() != '11':
            issues.append('unexpected numstep values %s -> %s' % (x, y))
    for key in ('batched_assembly', 'biot_effective_stress', 'temp_dep_viscosity'):
        if key in d['deck3']:
            issues.append('3D deck enables non-default switch %s' % key)
    wf = d['src']['workflow']
    for token in ('check_gate414_dim3.py --gate', 'Beam1x10-tetra4',
                  'gate414_verdict.log', 'out_dim3.gen', 'g414_rc'):
        if token not in wf:
            issues.append('build.yml missing %s' % token)
    # CI #120: freeze the diagnostic-only solver settings. Reverting to the
    # 2-rank ILU chain silently brings back the 900 s first-solve stall and
    # a truncated EXO that no W* check can read. -flow_ksp_monitor is what
    # supplies the "Linear solve converged" rows audited by W6.
    for token in ('-np 1 -x LD_LIBRARY_PATH', '-flow_sub_pc_type lu',
                  '-flow_ksp_monitor', 'timeout -k 10s 1800',
                  'if [ "$s414_rc" -ne 0 ]'):
        if token not in wf:
            issues.append('build.yml lost the CI #120 solver fix: %s' % token)
    if 'Beam1x10-tetra4' in d['log2'] or 'Gate 4.14:' in d['log2']:
        issues.append('3D diagnostic leaked into the 2D production log')
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'one-line deck (numstep 51->11), default physics, CI wired, 2D production isolated')


CHECKS = (('W1', 'run / mesh data integrity', check_w1),
          ('W2', '3D dimension dispatch evidence', check_w2),
          ('W3', 'static anti-regression', check_w3),
          ('W4', 'consolidation physics on shared times', check_w4),
          ('W5', '3D displacement field', check_w5),
          ('W6', '3D solver chain + fatal signatures', check_w6),
          ('W7', 'one-line deck + CI wiring + isolation', check_w7))


def run(d):
    rows = []
    for tag, label, fn in CHECKS:
        verdict, detail = fn(d)
        rows.append((tag, label, verdict, detail))
    fails = [r[0] for r in rows if r[2] == 'FAIL']
    return rows, fails


def report(rows, fails, header):
    print(header)
    for tag, label, verdict, detail in rows:
        print('%s %s %-42s | %s' % (verdict, tag, label, detail))
    print('GATE414_RESULT: %s%s' % ('PASS' if not fails else 'FAIL',
                                    '' if not fails else ' (' + ','.join(fails) + ')'))
    for tag, label, verdict, detail in rows:
        if verdict == 'FAIL':
            print('  %s detail: %s' % (tag, detail))


def synth_sources():
    return {k: read(v) for k, v in SRC.items()}


def synth():
    """Offline positive case: real 2D/3D mesh + source files, synthetic fields."""
    p2 = [1.0e5 * (0.0906 ** (j / 50.0)) for j in range(N_REC_2D)]
    # Identical initial pressure (Dirichlet-corrected p0), then a slightly
    # faster decay from the coarser dt=0.1 s BE step (3% at the final time).
    p3 = [p2[5 * i] * (1.0 - 0.03 * i / 10.0) for i in range(N_REC_3D)]
    ux = [1.13e-5 * (i / 10.0) for i in range(N_REC_3D)]
    uy = [7.6e-7 * (i / 10.0) for i in range(N_REC_3D)]
    uz = [1.2e-7 * (i / 10.0) for i in range(N_REC_3D)]
    exo3 = {'numrecs': N_REC_3D,
            'names': ['Pressure', 'Damage', 'Displacement_X', 'Displacement_Y', 'Displacement_Z'],
            'pmax': p3, 'ux': ux, 'uy': uy, 'uz': uz, 'finite': True}
    exo2 = {'numrecs': N_REC_2D,
            'names': ['Pressure', 'Damage', 'Displacement_X', 'Displacement_Y'],
            'pmax': p2, 'ux': ux[:], 'uy': uy[:], 'uz': [], 'finite': True}
    m3 = {'num_dim': 3, 'num_nodes': 2635, 'num_elem': 11873, 'num_el_blk': 3,
          'sets': [1, 10, 20], 'has_z': True,
          'rx': (-5.0, 5.0), 'ry': (-0.5, 0.5), 'rz': (-0.5, 0.5)}
    m2 = {'num_dim': 2, 'num_nodes': 356, 'num_elem': 610, 'num_el_blk': 3,
          'sets': [1, 10, 20], 'has_z': False,
          'rx': (-5.0, 5.0), 'ry': (-0.5, 0.5), 'rz': None}
    row = 'Gate 3.2: step %d pmin 0.0 pmax %.17g\n'
    log3 = 'x' * 1200 + '\n'
    for i in range(1, N_REC_3D + 1):
        log3 += 'Gate 3.1: step %d\n' % i
        log3 += row % (i, p3[i - 1])
        log3 += 'Gate 4.12: step %d t %.4g dpmax %.4g dpmin 0.0\n' % (i, (i - 1) * 0.1, -1e3)
        log3 += 'Gate 4.1: alpha 1.0\n'
        # Step 1 has no time advance; every later step must complete one
        # flow linear solve (CI #120 / W6).
        if i > 1:
            log3 += '    Linear solve converged due to CONVERGED_RTOL iterations 1\n'
    log2 = 'x' * 1200 + '\n' + 'Gate 3.1: step 1\n' * 51
    return {'log3': log3, 'log2': log2, 'exo3': exo3, 'exo2': exo2,
            'm3': m3, 'm2': m2, 'src': synth_sources(),
            'deck3': read(DECK_3D), 'deck2': read(DECK_2D)}


def selftest():
    import copy
    base = synth()
    cases = []
    rows, fails = run(base)
    cases.append(('POS full 3D verdict', rows, set(), None))
    print('POS  %-46s %s' % ('full 3D verdict', 'PASS' if not fails else 'FAIL ' + ','.join(fails)))

    def neg(name, mutate, expect):
        d = copy.deepcopy(base)
        mutate(d)
        r, f = run(d)
        ok = set(f) == set(expect)
        cases.append((name, r, set(expect), f))
        print('NEG  %-46s %s (failed %s, expected %s)'
              % (name, 'OK' if ok else 'MISMATCH', sorted(f), sorted(expect)))
        return ok

    allok = not fails
    allok &= neg('2D nodal table substituted for 3D',
                 lambda d: d['exo3'].update(names=d['exo2']['names']), {'W2'})
    allok &= neg('3D mesh reports num_dim=2',
                 lambda d: d['m3'].update(num_dim=2), {'W2'})
    allok &= neg('3D identity macro deleted',
                 lambda d: d['src'].update(inc=d['src']['inc'].replace('MEF90MatS3DIdentity', 'MEF90MatS2DIdentity')), {'W3'})
    allok &= neg('3D object twin no longer compiled',
                 lambda d: d['src'].update(flow_mk=d['src']['flow_mk'].replace('-DMEF90_DIM=3', '')), {'W3'})
    allok &= neg('runtime 3D dispatch removed',
                 lambda d: d['src'].update(flow_facade=d['src']['flow_facade'].replace('dim == 3', 'dim == 9')), {'W3'})
    allok &= neg('3D facade update pair removed (CI #119 regression)',
                 lambda d: d['src'].update(flow_facade=d['src']['flow_facade'].replace('Call MEF90FlowUpdatePerm3D(', 'Call MEF90FlowUpdatePerm2D(')), {'W3'})
    allok &= neg('2D facade update pair removed',
                 lambda d: d['src'].update(flow_facade=d['src']['flow_facade'].replace('Call MEF90FlowUpdatePerm2D(', 'Call MEF90FlowUpdatePerm3D(')), {'W3'})
    allok &= neg('3D mu(T) twin removed',
                 lambda d: d['src'].update(flow_facade=d['src']['flow_facade'].replace('Call MEF90FlowUpdateMuFromTemperature3D(', 'Call MEF90FlowUpdateMuFromTemperature2D(')), {'W3'})
    allok &= neg('compile-time 2D element macro reintroduced',
                 lambda d: d['src'].update(flow_facade=d['src']['flow_facade'].replace('Type(MEF90Element3D_Scal)', 'Type(MEF90_ELEMENT_SCAL)')), {'W3'})
    allok &= neg('3D pressure does not consolidate',
                 lambda d: d['exo3'].update(pmax=[1.0e5] * N_REC_3D), {'W4'})
    allok &= neg('3D displacement identically zero',
                 lambda d: d['exo3'].update(ux=[0.0] * N_REC_3D), {'W5'})
    allok &= neg('fatal signature in 3D log',
                 lambda d: d.update(log3=d['log3'] + 'PETSC ERROR: boom\n'), {'W1', 'W6'})
    allok &= neg('3D flow solve never converged (CI #120 timeout)',
                 lambda d: d.update(log3=d['log3'].replace('Linear solve converged', 'Linear solve stalled')), {'W6'})
    allok &= neg('CI #120 direct-PC fix reverted',
                 lambda d: d['src'].update(workflow=d['src']['workflow'].replace('-flow_sub_pc_type lu', '')), {'W7'})
    allok &= neg('CI #120 solver exit code no longer judged',
                 lambda d: d['src'].update(workflow=d['src']['workflow'].replace('if [ "$s414_rc" -ne 0 ]', 'if false')), {'W7'})
    allok &= neg('deck carries a second code change',
                 lambda d: d.update(deck3=d['deck3'].replace('p0: 1.0e5', 'p0: 2.0e5')), {'W7'})
    allok &= neg('CI wiring removed',
                 lambda d: d['src'].update(workflow=d['src']['workflow'].replace('out_dim3.gen', '')), {'W7'})
    allok &= neg('3D prefix leaked into production log',
                 lambda d: d.update(log2=d['log2'] + '-prefix data/Beam1x10-tetra4\n'), {'W7'})
    print('SELFTEST_RESULT: %s (%d positive + %d negatives)'
          % ('PASS' if allok else 'FAIL', 1, len(cases) - 1))
    return 0 if allok else 1


def main(argv):
    if argv and argv[0] == '--selftest':
        return selftest()
    if argv and argv[0] == '--gate':
        if len(argv) != 7:
            print('GATE414_RESULT: FAIL (usage: --gate <log3d> <exo3d> <log2d> '
                  '<exo2d> <mesh3d> <mesh2d>)')
            return 1
        try:
            d = {'log3': read(argv[1]), 'log2': read(argv[3]),
                 'exo3': exo_fields(argv[2]), 'exo2': exo_fields(argv[4]),
                 'm3': mesh_meta(argv[5]), 'm2': mesh_meta(argv[6]),
                 'src': synth_sources(),
                 'deck3': read(DECK_3D), 'deck2': read(DECK_2D)}
        except Exception as exc:
            # CI #120: a truncated out_dim3.gen (timeout killed MPI mid-write)
            # raised KeyError('name_nod_var') / struct.error deep inside CDF2
            # and produced a one-line verdict with no attributable check.
            # Report it as an explicit W1 failure with a diagnostic hint.
            print('=== Gate 4.14: 3D dimension dispatch (real CI logs/EXO) ===')
            print('FAIL W1 run / mesh data integrity | cannot parse 3D EXO/log: %r' % exc)
            print('GATE414_RESULT: FAIL (W1 data error: %r)' % exc)
            print('  hint: a truncated out_dim3.gen usually means the 3D run was '
                  'killed by the timeout before its first flow TSSolve converged')
            return 1
        rows, fails = run(d)
        report(rows, fails, '=== Gate 4.14: 3D dimension dispatch (real CI logs/EXO) ===')
        print('g414_records_3d = %d' % d['exo3']['numrecs'])
        print('g414_nodes_3d   = %d' % d['m3']['num_nodes'])
        print('g414_pend_3d    = %.6g' % (d['exo3']['pmax'][-1] if d['exo3']['pmax'] else -1))
        return 0 if not fails else 1
    print('usage: check_gate414_dim3.py [--gate <log3d> <exo3d> <log2d> '
          '<exo2d> <mesh3d> <mesh2d> | --selftest]')
    return 2


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
