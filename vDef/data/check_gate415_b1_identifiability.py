#!/usr/bin/env python3
"""Gate 4.15 / Stage B1 -- rate x fracability IDENTIFIABILITY hard verdict.

Why this gate exists
--------------------
Stage B0 (CI #125) calibrated the pre-initiation, storage-dominated window:
T = c_v*t/H^2 = 1.5667e-7, beta = dp_bulk/dt = f/S = 1.0e5 Pa/s,
S_eff = f/beta = 1.5e-10, Damage_max = 1.48e-4 << 1.

The ORIGINAL roadmap (RATE_FRACABILITY_ROADMAP.md sec. 5.2, criterion R4)
proposed to verify propagation-regime scalings on this window, e.g.
"viscosity-dominated p_net is independent of Q" and "d ln p_net / d ln Q ->
-1/3".  That premise is WRONG and the literature survey for B1 corrected it:

  Lecampion, Desroches, Jeffrey & Bunger (2017),
  "Experiments versus theory for the initiation and propagation of radial
  hydraulic fractures in low-permeability materials",
  JGR: Solid Earth, doi:10.1002/2016JB013183:
    "Prior to fracture initiation, the flow rate entering the fracture Qin
     is negligible, which for a constant injection rate Q0 results in a
     linear pressurization rate. The injection system compliance U can thus
     be directly obtained from the slope of the pressurization phase ...
     knowing Q0."
  They also distinguish INITIATION pressure from BREAKDOWN (peak) pressure.
  Detournay (2016) / Savitski & Detournay (2002) M/K vertices, and the
  Santillan, Juanes & Cueto-Felgueroso (2018) rate-mode thresholds, all
  describe the PROPAGATION stage and therefore do not constrain this deck.

In vDef the injection-system compliance is the storativity S, so the only
parameter identifiable from the pressure curve before initiation is U = S/f:

    beta = f/S,      S_eff = f/beta.

B1 therefore tests IDENTIFIABILITY LAYERING, which is falsifiable:

  Level 1 (pressure curve)  beta is linear in f, inversely linear in S, and
                            INVARIANT to Gc, k, mu as long as T << 1.
  Level 2 (phase field)     surface energy and the damage field DO respond
                            to Q and to Gc.  The surface-energy integrand is
                            C1 = Gc/(ell*cw)/4, C2 = Gc*ell/(4*cw)*A with
                            A = toughnessAnisotropyMatrix
                            (m_MEF90_DefMechAssembly.F90:2340-2341,2355),
                            so Gc is identifiable there even though beta
                            cannot see it.

Six new 2D decks, each with EXACTLY ONE code-line delta vs
terzaghi415_inject.yaml (comments excluded), plus the B0 control/injection
runs as the Q=0 / Q=1 anchors:

    b1_qhalf    source           1.5e-5 -> 7.5e-6     (Q = 0.5 x)
    b1_qdouble  source           1.5e-5 -> 3.0e-5     (Q = 2   x)
    b1_stor     storativity      1.5e-10 -> 3.0e-10   (S = 2   x)
    b1_gc       fracturetoughness 27000 -> 13500      (Gc = 0.5 x)
    b1_perm     perm_min         2.35e-18 -> 2.35e-17 (k = 10  x, T=1.57e-6)
    b1_visc     viscosity        1e-3 -> 2e-3         (mu = 2  x)

NO Fortran was changed.  NO propagation-regime claim is made.  The damage
SNES can still report reason -5 (known residual plateau, CI #84); the
"======= Total: ... surface: ..." rows are deterministic integrals of the
current field and are used ONLY as relative, falsifiable observables.

Verdicts
--------
A1 run / EXO / log / energy-row integrity (+ EXO vs Gate 3.2 cross-check)
A2 rate linearity: beta ~ f, S_eff = f/beta constant over 0.5x/1x/2x Q
A3 compliance dominance: beta ~ 1/S, S_eff recovers 2x S
A4 pressure-curve NON-identifiability of Gc / k / mu (Lecampion 2017)
A5 phase-field IDENTIFIABILITY of Gc through surface energy
A6 rate monotonicity of injection work, surface energy and damage
A7 static anti-regression + one-line deck contracts + CI wiring
A8 B0 verdict unchanged (CI #125) + production isolation

Pure standard library (CI runner has no numpy / scipy / PyYAML). --gate reads
real CI artefacts from one directory; --selftest replays one synthetic
positive case plus >=14 negative mutations and an invariance case.

Usage:
  check_gate415_b1_identifiability.py --gate <dir>
  check_gate415_b1_identifiability.py --selftest
"""

import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse_cdf2 import CDF2  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, '..', '..'))

N_REC = 11
DT = 0.1
P0 = 1.0e5
F0 = 1.5e-5
S_NUM = 1.50e-10
BETA0 = F0 / S_NUM                 # 1.0e5 Pa/s
AREA = 10.0
NAME_WIDTH = 33

TOL_BETA = 0.06                    # rate/compliance linearity
TOL_S_EFF = 0.10                   # recovered storativity
TOL_BETA_INV = 0.06                # Gc/k/mu must NOT move the pressure ramp
TOL_SURF_DECOUPLE = 0.05           # Gc must move surface energy by >5%
TOL_SURF_INV = 0.60                # k/mu must NOT move surface energy by >60%
TOL_CROSSCHECK = 1.0e-6
FATAL = ('PETSC ERROR', 'Caught signal', 'SEGV', 'fortran runtime error',
         'Fatal Error:', 'NaN detected')

B1_TAGS = ('qhalf', 'qdouble', 'stor', 'gc', 'perm', 'visc')
ALL_TAGS = B1_TAGS + ('b0_inject', 'b0_scale')
DECK_KEY = {
    'qhalf': 'source', 'qdouble': 'source', 'stor': 'storativity',
    'gc': 'fracturetoughness', 'perm': 'perm_min', 'visc': 'viscosity',
}
Q_FACTOR = {'b0_scale': 0.0, 'qhalf': 0.5, 'b0_inject': 1.0, 'qdouble': 2.0}

SRC = {
    'flowctx': os.path.join(_ROOT, 'm_Flow', 'm_MEF90_FlowCtx.F90'),
    'flowasm': os.path.join(_ROOT, 'm_Flow', 'm_MEF90_FlowAssembly.F90'),
    'defasm': os.path.join(_ROOT, 'm_DefMech', 'm_MEF90_DefMechAssembly.F90'),
    'workflow': os.path.join(_ROOT, '.github', 'workflows', 'build.yml'),
}
DECKS = {
    'inject': os.path.join(_HERE, 'terzaghi415_inject.yaml'),
    'scale': os.path.join(_HERE, 'terzaghi415_scale.yaml'),
}
for _t in B1_TAGS:
    DECKS[_t] = os.path.join(_HERE, 'terzaghi415_b1_%s.yaml' % _t)

FILEMAP = {}
for _t in ALL_TAGS:
    FILEMAP[_t] = ('%s.log' % _t, '%s.gen' % _t)
FILEMAP['b0_verdict'] = ('b0_verdict.log', None)
FILEMAP['main'] = ('main.log', None)
FILEMAP['mesh'] = (None, 'mesh.gen')

G32 = re.compile(
    r'Gate 3\.2:\s*step\s+(\d+)\s+pmin\s+([-+0-9.eEdD]+)\s+'
    r'pmax\s+([-+0-9.eEdD]+)', re.I)
G45 = re.compile(
    r'Gate 4\.5:\s*step\s+(\d+)\s+t\s+([-+0-9.eEdD]+)\s+'
    r'E\s+([-+0-9.eEdD]+)\s+W\s+([-+0-9.eEdD]+)', re.I)
NUM = r'([-+0-9.eEdD]+)'
ENERGY = re.compile(
    r'=======\s*Total:\s*elastic energy:\s*' + NUM +
    r'\s+work:\s*' + NUM + r'\s+cohesive:\s*' + NUM +
    r'\s+surface:\s*' + NUM + r'\s+total:\s*' + NUM, re.I)


def pass_fail(ok, detail):
    return ('PASS' if ok else 'FAIL'), detail


def num(s):
    return float(s.replace('D', 'E').replace('d', 'e'))


def read(path):
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        return f.read()


def var_names(c):
    raw = b''.join(c.read_static('name_nod_var'))
    return [raw[i:i + NAME_WIDTH].rstrip(b'\x00').decode('ascii', 'replace')
            for i in range(0, len(raw), NAME_WIDTH)]


def quantile(v, q):
    s = sorted(v)
    return s[min(len(s) - 1, int(q * len(s)))]


def exo_fields(path):
    c = CDF2(path)
    names = var_names(c)
    pbulk, pmax, dmax = [], [], []
    for i in range(c.numrecs):
        p = list(c.read_record('vals_nod_var1', i))
        pbulk.append(quantile(p, 0.50))
        pmax.append(max(p))
        d = list(c.read_record('vals_nod_var2', i))
        dmax.append(max(d))
    vals = pbulk + pmax + dmax
    return {'numrecs': c.numrecs, 'names': names,
            'pbulk': pbulk, 'pmax': pmax, 'dmax': dmax,
            'finite': all(math.isfinite(x) for x in vals)}


def mesh_meta(path):
    c = CDF2(path)
    d = dict(zip(c.dims, c.sizes))
    return {'num_dim': d.get('num_dim'), 'num_nodes': d.get('num_nodes'),
            'sets': list(c.read_static('eb_prop1'))}


def energy_rows(log):
    out = []
    for m in ENERGY.finditer(log):
        out.append(tuple(num(m.group(i)) for i in range(1, 6)))
    return out


def g32_rows(log):
    return [(int(m.group(1)), num(m.group(2)), num(m.group(3)))
            for m in G32.finditer(log)]


def g45_rows(log):
    return [(int(m.group(1)), num(m.group(2)), num(m.group(3)), num(m.group(4)))
            for m in G45.finditer(log)]


def linfit(ts, ys):
    tb = sum(ts) / len(ts)
    yb = sum(ys) / len(ys)
    den = sum((t - tb) ** 2 for t in ts)
    if den <= 0.0:
        return 0.0, yb
    m = sum((t - tb) * (y - yb) for t, y in zip(ts, ys)) / den
    return m, yb - m * tb


def beta_of(d, tag):
    """Bulk pressure ramp [Pa/s], always measured from the real EXO."""
    p = d['exo'][tag]['pbulk'][:N_REC]
    return linfit([i * DT for i in range(N_REC)], [x - P0 for x in p])[0]


def code_lines(text):
    out = []
    for line in text.splitlines():
        s = line.split('#', 1)[0].rstrip()
        if s.strip():
            out.append(s.strip())
    return out


def code_diff(a, b):
    la, lb = code_lines(a), code_lines(b)
    if len(la) != len(lb):
        return None
    return [(x, y) for x, y in zip(la, lb) if x != y]


def key_of(line):
    return line.split(':', 1)[0].strip()


def fatal_in(log):
    return [s for s in FATAL if s in log]


# ---------------------------------------------------------------------------
# Verdicts
# ---------------------------------------------------------------------------

def check_a1(d):
    issues = []
    for tag in ALL_TAGS:
        e, log = d['exo'][tag], d['log'][tag]
        if len(log) < 1000:
            issues.append('%s log too short' % tag)
        if e['numrecs'] != N_REC:
            issues.append('%s EXO records %d != %d' % (tag, e['numrecs'], N_REC))
        if not e['finite']:
            issues.append('%s EXO has non-finite pressure/damage' % tag)
        for nm in ('Pressure', 'Damage', 'Displacement_X', 'Displacement_Y'):
            if nm not in e['names']:
                issues.append('%s EXO lost %s' % (tag, nm))
        # CI #124 lesson: a shortened run can silently APPEND to an old mesh.
        logp = [r[2] for r in g32_rows(log)]
        if len(logp) < N_REC:
            issues.append('%s has only %d Gate 3.2 rows' % (tag, len(logp)))
        elif len(e['pmax']) >= N_REC:
            worst = max(abs(a - b) / max(P0, abs(b))
                        for a, b in zip(e['pmax'][:N_REC], logp[:N_REC]))
            if worst > TOL_CROSSCHECK:
                issues.append('%s EXO pmax != Gate 3.2 log (rel %.3g): stale mesh'
                              % (tag, worst))
        for label, n, need in (('Gate 4.5', len(g45_rows(log)), N_REC),
                               ('energy row', len(energy_rows(log)), N_REC),
                               ('KSP solve', log.count('Linear solve converged'),
                                N_REC - 1)):
            if n < need:
                issues.append('%s has %d %s (need >=%d)' % (tag, n, label, need))
        for sig in fatal_in(log):
            issues.append('fatal in %s: %s' % (tag, sig))
    m = d['mesh']
    if m['num_dim'] != 2 or m['num_nodes'] != 356 or m['sets'] != [1, 10, 20]:
        issues.append('2D mesh metadata changed: %r' % m)
    return pass_fail(not issues, '; '.join(issues) if issues else
                     '%d runs x %d clean records, 11 Gate 3.2/4.5/energy rows '
                     'each, EXO pmax == solver, 356-node tri3, no fatal'
                     % (len(ALL_TAGS), N_REC))


def check_a2(d):
    issues = []
    betas = {}
    for tag, q in (('qhalf', 0.5), ('b0_inject', 1.0), ('qdouble', 2.0)):
        b = beta_of(d, tag)
        betas[tag] = b
        d['metrics']['beta_' + tag] = b
        f = F0 * q
        if b <= 0.0:
            issues.append('%s beta=%.6g <= 0' % (tag, b))
            continue
        if abs(b / (f / S_NUM) - 1.0) > TOL_BETA:
            issues.append('%s beta=%.6g vs f/S=%.6g (%.1f%% off)'
                          % (tag, b, f / S_NUM, 100 * abs(b / (f / S_NUM) - 1)))
        s_eff = f / b
        d['metrics']['S_eff_' + tag] = s_eff
        if abs(s_eff / S_NUM - 1.0) > TOL_S_EFF:
            issues.append('%s S_eff=%.4g != S_num=%.4g' % (tag, s_eff, S_NUM))
    if betas.get('b0_inject', 0) > 0:
        for tag, q in (('qhalf', 0.5), ('qdouble', 2.0)):
            r = betas.get(tag, 0) / betas['b0_inject']
            d['metrics']['beta_ratio_' + tag] = r
            if abs(r / q - 1.0) > TOL_BETA:
                issues.append('beta(%s)/beta0=%.4f != %.1f' % (tag, r, q))
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'beta(0.5x/1x/2x Q) = %.4g / %.4g / %.4g Pa/s, '
                     'ratios %.3f / %.3f, S_eff all = %.4g'
                     % (betas.get('qhalf', -1), betas.get('b0_inject', -1),
                        betas.get('qdouble', -1),
                        d['metrics'].get('beta_ratio_qhalf', -1),
                        d['metrics'].get('beta_ratio_qdouble', -1), S_NUM))


def check_a3(d):
    issues = []
    b0 = beta_of(d, 'b0_inject')
    bs = beta_of(d, 'stor')
    d['metrics']['beta_stor'] = bs
    d['metrics']['S_eff_stor'] = (F0 / bs) if bs > 0 else -1.0
    if b0 <= 0 or bs <= 0:
        issues.append('non-positive beta in the compliance arm')
    else:
        r = bs / b0
        d['metrics']['beta_stor_ratio'] = r
        if abs(r / 0.5 - 1.0) > TOL_BETA:
            issues.append('2x S gave beta ratio %.4f != 0.5' % r)
        # Falsification: a model where beta depends on f only is rejected.
        if r > 0.9:
            issues.append('beta barely responds to S (%.3f): compliance is '
                          'not the controlling storage' % r)
        if abs((F0 / bs) / (2.0 * S_NUM) - 1.0) > TOL_S_EFF:
            issues.append('S_eff(stor)=%.4g != 2x S_num' % (F0 / bs))
    return pass_fail(not issues, '; '.join(issues) if issues else
                     '2x storativity -> beta ratio %.4f, recovered S_eff=%.4g '
                     '(2 x S_num), U = S/f is the identifiable compliance'
                     % (d['metrics'].get('beta_stor_ratio', -1),
                        d['metrics'].get('S_eff_stor', -1)))


def check_a4(d):
    """Lecampion et al. (2017): rock parameters are invisible BEFORE initiation."""
    issues = []
    b0 = beta_of(d, 'b0_inject')
    if b0 <= 0:
        return pass_fail(False, 'B0 reference beta is non-positive')
    for tag in ('gc', 'perm', 'visc'):
        b = beta_of(d, tag)
        dev = abs(b / b0 - 1.0)
        d['metrics']['beta_' + tag] = b
        d['metrics']['beta_dev_' + tag] = dev
        if dev > TOL_BETA_INV:
            issues.append('%s moved beta by %.2f%% (>%.0f%%): the pressure '
                          'curve is no longer compliance-dominated'
                          % (tag, 100 * dev, 100 * TOL_BETA_INV))
        if b > 0 and abs((F0 / b) / S_NUM - 1.0) > TOL_S_EFF:
            issues.append('%s no longer recovers S_eff' % tag)
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'Gc / k / mu leave beta within %.2f%% / %.2f%% / %.2f%%: '
                     'pre-initiation rock parameters are NOT identifiable on '
                     'the pressure curve (Lecampion 2017)'
                     % (100 * d['metrics'].get('beta_dev_gc', -1),
                        100 * d['metrics'].get('beta_dev_perm', -1),
                        100 * d['metrics'].get('beta_dev_visc', -1)))


def check_a5(d):
    issues = []
    surf = {}
    for tag in ALL_TAGS:
        er = energy_rows(d['log'][tag])
        if len(er) < N_REC:
            return pass_fail(False, '%s has only %d energy rows' % (tag, len(er)))
        surf[tag] = er[N_REC - 1][3]
        if any(not all(math.isfinite(x) for x in r) for r in er[:N_REC]):
            issues.append('%s energy row is non-finite' % tag)
        if er[0][0] <= 0.0:
            issues.append('%s elastic energy at step 1 is not positive' % tag)
    ref = surf['b0_inject']
    d['metrics']['surface_ref'] = ref
    if ref <= 0.0:
        issues.append('reference surface energy is not positive')
    else:
        for tag in ('gc', 'perm', 'visc'):
            r = surf[tag] / ref
            d['metrics']['surface_ratio_' + tag] = r
            if tag == 'gc':
                if abs(r - 1.0) < TOL_SURF_DECOUPLE:
                    issues.append('Gc halved but surface energy moved only '
                                  '%.2f%%: the C1/C2 ~ Gc channel is unwired'
                                  % (100 * abs(r - 1.0)))
                if not (0.02 < r < 5.0):
                    issues.append('Gc arm surface ratio %.4g is unphysical' % r)
            elif abs(r - 1.0) > TOL_SURF_INV:
                issues.append('%s moved surface energy by %.1f%% (>%.0f%%): '
                              'k/mu should be invisible here'
                              % (tag, 100 * abs(r - 1.0), 100 * TOL_SURF_INV))
    # The DECOUPLING is the actual falsifiable claim: beta unchanged (A4),
    # surface energy changed (A5).  Reject a deck where both are frozen.
    if d['metrics'].get('beta_dev_gc', 1.0) > TOL_BETA_INV:
        issues.append('A4/A5 decoupling broken: Gc also moved the pressure ramp')
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'surface_ref=%.4g; Gc x0.5 -> x%.3f, k x10 -> x%.3f, '
                     'mu x2 -> x%.3f; beta invariant while Gc is visible in '
                     'the phase field (identifiability layering)'
                     % (ref, d['metrics'].get('surface_ratio_gc', -1),
                        d['metrics'].get('surface_ratio_perm', -1),
                        d['metrics'].get('surface_ratio_visc', -1)))


def check_a6(d):
    issues = []
    order = ['b0_scale', 'qhalf', 'b0_inject', 'qdouble']
    w, surf, dmg = [], [], []
    for tag in order:
        g45 = g45_rows(d['log'][tag])
        er = energy_rows(d['log'][tag])
        if len(g45) < N_REC or len(er) < N_REC:
            return pass_fail(False, '%s lacks the 11-step work/energy series' % tag)
        w.append(g45[N_REC - 1][3])
        surf.append(er[N_REC - 1][3])
        dmg.append(max(d['exo'][tag]['dmax'][:N_REC]))
    d['metrics']['w_series'] = w
    d['metrics']['surface_series'] = surf
    d['metrics']['damage_series'] = dmg
    if abs(w[0]) > 1.0e-12:
        issues.append('zero-source control has non-zero injection work')
    for label, series, tol in (('W', w, 1.0e-9),
                               ('surface energy', surf, 0.0),
                               ('damage max', dmg, 0.0)):
        for i in range(1, len(series)):
            if series[i] < series[i - 1] - tol:
                issues.append('%s is not monotone from Q=0 to Q=2x '
                              '(%.4g -> %.4g at %s)' %
                              (label, series[i - 1], series[i], order[i]))
    if w[-1] <= w[1]:
        issues.append('Q=2x injection work is not above Q=0.5x')
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'Q=0/0.5/1/2x: W_last=%.4g/%.4g/%.4g/%.4g, '
                     'surface_last=%.3g/%.3g/%.3g/%.3g, '
                     'damage_max=%.3g/%.3g/%.3g/%.3g (all monotone)'
                     % (w[0], w[1], w[2], w[3],
                        surf[0], surf[1], surf[2], surf[3],
                        dmg[0], dmg[1], dmg[2], dmg[3]))


def check_a7(d):
    issues = []
    s = d['src']
    if 'Volumetric source/sink rate (1/s)' not in s['flowctx']:
        issues.append('flow.source [1/s] registration lost')
    for token in ('MEF90Scaling_Linear', 'MEF90Scaling_CST',
                  'MEF90FlowSetSourceCst'):
        if token not in s['flowasm']:
            issues.append('m_MEF90_FlowAssembly.F90 lost %s' % token)
    for token in ('matpropSet%fractureToughness',
                  'matpropSet%toughnessAnisotropyMatrix',
                  'Subroutine MEF90DefMechSurfaceEnergy'):
        if token not in s['defasm']:
            issues.append('m_MEF90_DefMechAssembly.F90 lost %s' % token)
    base = d['deck']['inject']
    for tag in B1_TAGS:
        diff = code_diff(base, d['deck'][tag])
        if diff is None or len(diff) != 1:
            issues.append('%s deck has %s code deltas (expected 1)'
                          % (tag, len(diff) if diff is not None else 'length mismatch'))
        elif key_of(diff[0][0]) != DECK_KEY[tag]:
            issues.append('%s deck changed %s, expected %s'
                          % (tag, key_of(diff[0][0]), DECK_KEY[tag]))
        if 'Lecampion' not in d['deck'][tag]:
            issues.append('%s deck lost the literature anchor' % tag)
        for bad in ('batched_assembly', 'biot_effective_stress',
                    'temp_dep_viscosity'):
            if bad in d['deck'][tag]:
                issues.append('%s deck enables non-default switch %s' % (tag, bad))
    wf = s['workflow']
    for tag in B1_TAGS:
        if 'terzaghi415_b1_%s.yaml' % tag not in wf:
            issues.append('build.yml missing the %s run' % tag)
    for token in ('check_gate415_b1_identifiability.py --gate',
                  'gate415_b1_verdict.log', 'gb1_rc',
                  'if [ "$gb1_rc" -ne 0 ]',
                  'if [ "$sb1_qdouble_rc" -ne 0 ]'):
        if token not in wf:
            issues.append('build.yml missing %s' % token)
    # CI #124 lesson: every 11-step 2D diagnostic must rebuild its EXO.
    n_rebuild = wf.count('rm -f data/Beam1x10-tri3_out.gen')
    d['metrics']['n_rebuild'] = n_rebuild
    if n_rebuild < 8:
        issues.append('only %d EXO rebuilds (need >=8: B0 x2 + B1 x6)'
                      % n_rebuild)
    return pass_fail(not issues, '; '.join(issues) if issues else
                     '6 one-line decks (%s), source/energy/anisotropy wiring, '
                     'CI runs + verdict + deferred fail + %d EXO rebuilds'
                     % (','.join(sorted(set(DECK_KEY.values()))), n_rebuild))


def check_a8(d):
    issues = []
    v = d['log']['b0_verdict']
    if 'GATE415_RESULT: PASS' not in v:
        issues.append('B0 verdict is no longer PASS')
    for token in ('g415_slope_measured        = 1.000000e+05',
                  'g415_S_eff                 = 1.500000e-10'):
        if token not in v:
            issues.append('B0 CI #125 number drifted: %s' % token)
    main = d['log']['main']
    rows = g32_rows(main)
    if len(rows) < 50:
        issues.append('production log has only %d Gate 3.2 rows' % len(rows))
    elif rows[-1][2] >= 0.5 * P0:
        issues.append('production pressure did not consolidate')
    wm = [r[3] for r in g45_rows(main)]
    if wm and any(abs(w) > 1.0e-12 for w in wm):
        issues.append('production injection work is non-zero')
    for bad in ('terzaghi415_b1_', 'check_gate415_b1_'):
        if bad in main:
            issues.append('a B1 diagnostic leaked into the production log: %s' % bad)
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'B0 verdict still CI #125 PASS with beta=1e5 / S_eff=1.5e-10; '
                     'production consolidates over %d steps, W=0, no B1 leakage'
                     % len(rows))


CHECKS = (('A1', 'run / EXO / log / energy-row integrity', check_a1),
          ('A2', 'rate linearity beta ~ f, S_eff constant', check_a2),
          ('A3', 'compliance dominance beta ~ 1/S', check_a3),
          ('A4', 'Gc/k/mu invisible on the pressure curve', check_a4),
          ('A5', 'Gc visible in phase-field surface energy', check_a5),
          ('A6', 'rate monotonicity of W / surface / damage', check_a6),
          ('A7', 'static anti-regression + deck/CI contracts', check_a7),
          ('A8', 'B0 unchanged (CI #125) + production isolation', check_a8))


def run(d):
    d['metrics'] = {}
    rows = []
    for tag, label, fn in CHECKS:
        verdict, detail = fn(d)
        rows.append((tag, label, verdict, detail))
    fails = [r[0] for r in rows if r[2] == 'FAIL']
    return rows, fails


def report(rows, fails, header):
    print(header)
    for tag, label, verdict, detail in rows:
        print('%s %s %-46s | %s' % (verdict, tag, label, detail))
    print('GATE415_B1_RESULT: %s%s' % ('PASS' if not fails else 'FAIL',
                                       '' if not fails else ' (' + ','.join(fails) + ')'))
    for tag, label, verdict, detail in rows:
        if verdict == 'FAIL':
            print('  %s detail: %s' % (tag, detail))
    m = ROWS_METRICS.get('metrics', {})
    for key in ('beta_qhalf', 'beta_b0_inject', 'beta_qdouble', 'beta_stor',
                'beta_gc', 'beta_perm', 'beta_visc',
                'beta_dev_gc', 'beta_dev_perm', 'beta_dev_visc',
                'surface_ref', 'surface_ratio_gc', 'surface_ratio_perm',
                'surface_ratio_visc', 'S_eff_stor'):
        print('gb1_%-26s = %.6e' % (key, m.get(key, -1)))


ROWS_METRICS = {}


def synth_sources():
    return {k: read(v) for k, v in SRC.items()}


def synth_decks():
    return {k: read(v) for k, v in DECKS.items()}


B0_VERDICT = """=== Gate 4.15 / Stage B0: rate/fracability similitude ===
PASS V1 run / EXO / mesh data integrity | 11 + 11 clean records
PASS V4 injection slope f/S + source-unit calibration | bulk slope=100000 Pa/s
GATE415_RESULT: PASS
g415_T                     = 1.566667e-07
g415_slope_measured        = 1.000000e+05
g415_S_eff                 = 1.500000e-10
g415_peak_ratio            = 2.918774e+00
"""


def make_log(tag, pmax, work, surface, damage, n=N_REC):
    """Synthetic solver log with the exact real row formats."""
    log = '=== Gate 4.15 B1 synthetic run %s ===\n' % tag + 'x' * 1200 + '\n'
    for i in range(1, n + 1):
        log += 'Gate 3.1: step %d enter transient block\n' % i
        log += 'Gate 3.2: step %d pmin 0.0 pmax %.10E\n' % (i, pmax[i - 1])
        log += ('Gate 4.12: step %d t %.4g dpmax %.4g dpmin 0.0\n'
                % (i, (i - 1) * DT, pmax[i - 1] - pmax[i - 2] if i > 1 else 0.0))
        log += ('Gate 4.5: step %d t %.10E E %.10E W %.10E\n'
                % (i, (i - 1) * DT, 1.0e-3, work[i - 1]))
        log += ('======= Total: elastic energy:  %.5E work:  %.5E '
                'cohesive:  0.00000E+00 surface:  %.5E total:  %.5E\n'
                % (0.5 * work[i - 1] or 1.0e-6, work[i - 1] or 1.0e-6,
                   surface[i - 1], 1.0))
        if i > 1:
            log += '    Linear solve converged due to CONVERGED_RTOL iterations 1\n'
    return log


def synth():
    """Offline positive case: real decks/sources, analytic solver artefacts.

    Pressures use the exact backward-Euler anchor beta = f/S (CI #125 gave
    1.000000e5 Pa/s exactly).  Surface energies use the B0 measured
    surface_last = 6.49012e-4 for Q=1x and the C1/C2 ~ Gc scaling for the
    Gc arm; k/mu arms are held invariant as the Level-1 hypothesis predicts.
    """
    d = {'log': {}, 'exo': {}, 'metrics': {}}
    names = ['Pressure', 'Damage', 'Displacement_X', 'Displacement_Y']
    for tag in ALL_TAGS:
        q = Q_FACTOR.get(tag, 1.0)
        f = F0 * q
        S = S_NUM * (2.0 if tag == 'stor' else 1.0)
        beta = f / S
        pbulk = [P0 + beta * i * DT for i in range(N_REC)]
        pmax = [P0 + 2.9184 * beta * i * DT for i in range(N_REC)]
        work = [f * AREA * p for p in pbulk]
        # Gc halves -> surface energy drops, but stays finite and non-zero.
        sf = 0.55 if tag == 'gc' else 1.0
        surf_last = (6.49012e-4 * sf *
                     (max(q, 1.0e-3) ** 1.6) if q > 0 else 1.20755e-6)
        surface = [surf_last * ((i + 1) / N_REC) ** 2 for i in range(N_REC)]
        dmg_last = (1.4803e-4 * (max(q, 1.0e-3) ** 1.5)
                    if q > 0 else 8.6005e-6)
        dmax = [dmg_last * (i + 1) / N_REC for i in range(N_REC)]
        d['exo'][tag] = {'numrecs': N_REC, 'names': names, 'finite': True,
                         'pbulk': pbulk, 'pmax': pmax, 'dmax': dmax}
        d['log'][tag] = make_log(tag, pmax, work, surface, dmax)
    p_main = [P0 * (0.0906 ** (j / 50.0)) for j in range(51)]
    d['log']['main'] = make_log('main', p_main, [0.0] * 51,
                                [1.0e-6] * 51, [1.0e-6] * 51, n=51)
    d['log']['b0_verdict'] = B0_VERDICT
    d['mesh'] = {'num_dim': 2, 'num_nodes': 356, 'sets': [1, 10, 20]}
    d['src'] = synth_sources()
    d['deck'] = synth_decks()
    return d


def set_last_work(log, value):
    """Rewrite the W entry of the LAST Gate 4.5 row (used by A6)."""
    pat = r'(Gate 4\.5:\s*step\s+%d\s+.*?\s+W\s+)[-+0-9.eEdD]+' % N_REC
    return re.sub(pat, lambda m: m.group(1) + '%.10E' % value, log, count=1)


def selftest():
    import copy
    base = synth()
    cases = []
    rows, fails = run(base)
    ROWS_METRICS['metrics'] = base['metrics']
    cases.append(('POS full B1 verdict', rows, set(), None))
    print('POS  %-54s %s' % ('full B1 verdict', 'PASS' if not fails else 'FAIL ' + ','.join(fails)))

    def neg(name, mutate, expect, label='NEG'):
        d = copy.deepcopy(base)
        mutate(d)
        r, f = run(d)
        ok = set(f) == set(expect)
        cases.append((name, r, set(expect), f))
        print('%s  %-54s %s (failed %s, expected %s)'
              % (label, name, 'OK' if ok else 'MISMATCH',
                 sorted(f), sorted(expect)))
        return ok

    def inv(name, mutate):
        return neg(name, mutate, set(), 'INV')

    allok = not fails
    allok &= neg('qhalf EXO truncated to 10 records',
                 lambda d: d['exo']['qhalf'].update(numrecs=10), {'A1'})
    allok &= neg('fatal signature in the gc log',
                 lambda d: d['log'].update(
                     gc=d['log']['gc'] + 'PETSC ERROR: boom\n'), {'A1'})
    allok &= neg('EXO pmax disagrees with Gate 3.2 (stale mesh)',
                 lambda d: d['exo']['qdouble']['pmax'].__setitem__(
                     3, d['exo']['qdouble']['pmax'][3] * 0.5), {'A1'})
    allok &= neg('energy rows missing from the perm log',
                 lambda d: d['log'].update(
                     perm=d['log']['perm'].replace('======= Total:', '==== NOPE:')),
                 {'A1', 'A5'})
    allok &= neg('Q=2x deck secretly runs at Q=1x',
                 lambda d: d['exo']['qdouble'].update(
                     pbulk=d['exo']['b0_inject']['pbulk'][:]), {'A2'})
    allok &= neg('Q=0.5x deck secretly runs at Q=1x',
                 lambda d: d['exo']['qhalf'].update(
                     pbulk=d['exo']['b0_inject']['pbulk'][:]), {'A2'})
    allok &= neg('storativity arm still ramps at beta0',
                 lambda d: d['exo']['stor'].update(
                     pbulk=d['exo']['b0_inject']['pbulk'][:]), {'A3'})
    allok &= neg('Gc arm breaks the pressure ramp (propagation-like)',
                 lambda d: d['exo']['gc'].update(
                     pbulk=[P0 + 0.5 * BETA0 * i * DT for i in range(N_REC)]),
                 {'A4', 'A5'})
    allok &= neg('k arm breaks the pressure ramp (diffusion leaks in)',
                 lambda d: d['exo']['perm'].update(
                     pbulk=[P0 + 1.5 * BETA0 * i * DT for i in range(N_REC)]),
                 {'A4'})
    allok &= neg('mu arm breaks the pressure ramp',
                 lambda d: d['exo']['visc'].update(
                     pbulk=[P0 + 0.7 * BETA0 * i * DT for i in range(N_REC)]),
                 {'A4'})
    allok &= neg('Gc halved but surface energy frozen (C1/C2 unwired)',
                 lambda d: d['log'].update(
                     gc=d['log']['gc'].replace(
                         'surface:  3.56957E-04', 'surface:  6.49012E-04')),
                 {'A5'})
    allok &= neg('k arm corrupts the surface energy',
                 lambda d: d['log'].update(
                     perm=d['log']['perm'].replace(
                         'surface:  6.49012E-04', 'surface:  6.49012E-02', 1)),
                 {'A5'})
    allok &= neg('zero-source control reports non-zero work',
                 lambda d: d['log'].update(
                     b0_scale=d['log']['b0_scale'].replace(
                         'W 0.0000000000E+00', 'W 3.5000000000E+00')), {'A6'})
    allok &= neg('Q=2x injection work not monotone',
                 lambda d: d['log'].update(
                     qdouble=set_last_work(d['log']['qdouble'], 1.0e-3)),
                 {'A6'})
    allok &= neg('one-line deck contract broken (second code delta)',
                 lambda d: d['deck'].update(
                     gc=d['deck']['gc'].replace('p0: 1.0e5', 'p0: 2.0e5')),
                 {'A7'})
    allok &= neg('all EXO rebuilds removed from the B1/B2 runs',
                 lambda d: d['src'].update(workflow=d['src']['workflow'].replace(
                     'rm -f data/Beam1x10-tri3_out.gen', 'true # no rebuild')),
                 {'A7'})
    allok &= neg('deferred gb1_rc hard fail removed',
                 lambda d: d['src'].update(workflow=d['src']['workflow'].replace(
                     'if [ "$gb1_rc" -ne 0 ]', 'if false')), {'A7'})
    allok &= neg('B0 verdict silently became FAIL',
                 lambda d: d['log'].update(
                     b0_verdict=B0_VERDICT.replace('GATE415_RESULT: PASS',
                                                   'GATE415_RESULT: FAIL')), {'A8'})
    allok &= neg('B1 deck name leaked into the production log',
                 lambda d: d['log'].update(
                     main=d['log']['main'] + '-options_file_yaml data/terzaghi415_b1_gc.yaml\n'),
                 {'A8'})
    allok &= inv('damage SNES residual plateau (-5) appears, physics unchanged',
                 lambda d: d['log'].update(
                     gc=d['log']['gc'] + '[ERROR]: damage field SNESSolve failed with SNESConvergedReason -5.\n'))
    print('SELFTEST_RESULT: %s (%d positive + %d mutations, incl. 1 invariance)'
          % ('PASS' if allok else 'FAIL', 1, len(cases) - 1))
    return 0 if allok else 1


def main(argv):
    if argv and argv[0] == '--selftest':
        return selftest()
    if argv and argv[0] == '--gate':
        if len(argv) != 2:
            print('GATE415_B1_RESULT: FAIL (usage: --gate <artifact_dir>)')
            return 1
        root = argv[1]
        try:
            d = {'log': {}, 'exo': {}, 'metrics': {}}
            for tag in ALL_TAGS + ('b0_verdict', 'main'):
                lp, ep = FILEMAP[tag]
                d['log'][tag] = read(os.path.join(root, lp))
            for tag in ALL_TAGS:
                d['exo'][tag] = exo_fields(os.path.join(root, FILEMAP[tag][1]))
            d['mesh'] = mesh_meta(os.path.join(root, FILEMAP['mesh'][1]))
            d['src'] = synth_sources()
            d['deck'] = synth_decks()
        except Exception as exc:
            print('=== Gate 4.15 / Stage B1: rate x fracability identifiability ===')
            print('FAIL A1 run / EXO / log / energy-row integrity | cannot parse artefacts: %r' % exc)
            print('GATE415_B1_RESULT: FAIL (A1 data error: %r)' % exc)
            return 1
        rows, fails = run(d)
        ROWS_METRICS['metrics'] = d['metrics']
        report(rows, fails, '=== Gate 4.15 / Stage B1: rate x fracability '
               'identifiability (real CI logs/EXO) ===')
        return 0 if not fails else 1
    print('usage: check_gate415_b1_identifiability.py '
          '[--gate <artifact_dir> | --selftest]')
    return 2


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
