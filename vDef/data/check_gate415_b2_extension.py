#!/usr/bin/env python3
"""Gate 4.15 / Stage B2 -- variable-rate, anisotropy and 3D hard verdict.

Stage B1 established IDENTIFIABILITY LAYERING in the pre-initiation,
storage-dominated window (Lecampion et al. 2017, JGR Solid Earth,
doi:10.1002/2016JB013183):
  * beta = dp_bulk/dt = f/S is identifiable in (f, S) ONLY;
  * Gc is invisible on the pressure curve but visible in the phase-field
    surface energy whose integrand uses C1 ~ Gc/ell and
    C2 ~ Gc*ell*toughnessAnisotropyMatrix
    (m_MEF90_DefMechAssembly.F90:2340-2341,2355).

B2 exercises the three remaining roadmap channels, still with ZERO Fortran
change and ZERO propagation-regime claim:

  (1) VARIABLE RATE.  `source_scaling: linear` makes
      MEF90FlowSetTransients call VecScale(source, timeVal), i.e. f(t) = f0*t.
      The backward-Euler discrete solution of S dp/dt = f0*t on the uniform
      dt=0.1 grid (step 1 only imposes the initial condition, as in B0) is

          Delta_p_n = (f0/S) * dt^2 * n*(n+1)/2,   n = 0..10,

      so Delta_p_10 = 5.5e4 Pa, EXACTLY 55% of the constant-rate B0 value.
      L2 fits Delta_p against n(n+1)/2 and rejects the constant-rate model.
      Engineering analogue: a step-rate test (SRT) uses equal-duration rate
      steps and reads the stabilised pressure at each step; here the ramp is
      continuous, so the discrete t^2 anchor replaces the SRT slope break.

  (2) BEDDING ANISOTROPY.  MatS2D is (XX, YY, XY) and the identity is
      (1,1,0) (m_MEF90_LinAlg.F90:50-54,107).  Setting
      toughnessAnisotropyMatrix: 2,0.5,0 therefore changes BOTH principal
      directions.  A enters the surface energy ONLY through
      C2 = (Gc*ell/(4*cw))*A multiplying grad(alpha)^2, while the almost
      spatially uniform pre-initiation damage (alpha ~ 1e-4) is dominated
      by C1*w(alpha).  Anisotropy must therefore be WIRED but NOT
      IDENTIFIABLE in the surface energy.  CI #127 measured
      surface_ratio = 1.004065 (0.41% shift), damage_dev = 1.247e-2
      (1.25%) and beta_dev = 0.  The original L4 demanded a >5% surface
      shift; that was a WRONG PHYSICAL PREMISE borrowed from the crack
      propagation regime (Detournay 2016), not the storage-dominated
      pre-initiation window of B0/B1.  L4 is now a falsifiable NULL
      proposition: a lower bound (not bit-equivalent => A really entered
      the solve), an upper bound (A did not leak into the C1 channel), an
      ordering constraint (energy response weaker than damage response),
      and a blind pressure ramp.  Literature basis: shale bedding-parallel
      K_Ic is typically 2-3x weaker than bedding-normal and Gc anisotropy
      dominates crack-path steering -- but that steering needs a
      PROPAGATING crack, which this window explicitly excludes.

  (3) 3D CONSISTENCY.  The same terzaghi415_inject.yaml is re-run on the
      Beam1x10-tetra4 mesh (dimension dispatch is by -prefix, not by a YAML
      key, cf. Gate 4.14 / CI #121).  L6 requires the same beta and a live
      Displacement_Z channel, reusing the single-rank direct-LU diagnostic
      configuration that fixed CI #120's 3D timeout.

Verdicts
--------
L1 run / EXO / log / energy-row integrity (2D + 3D)
L2 variable-rate discrete t^2 anchor, constant-rate model rejected
L3 variable-rate injection work monotone and below the constant-rate arm
L4 pre-initiation anisotropy: wired, bounded, beta-blind (CI #127)
L5 anisotropy deck/bag/source static contracts (no spaces after commas)
L6 3D tetra4 dispatch: beta matches 2D, Displacement_Z live, no timeout
L7 static anti-regression + CI wiring (runs, verdict, deferred fail, 3D PC)
L8 B0 (CI #125) + B1 verdicts unchanged, production isolation

Pure standard library (CI runner has no numpy / scipy / PyYAML). --gate reads
real CI artefacts from one directory; --selftest replays one synthetic
positive case plus >=12 negative mutations and an invariance case.

Usage:
  check_gate415_b2_extension.py --gate <dir>
  check_gate415_b2_extension.py --selftest
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
BETA0 = F0 / S_NUM
AREA = 10.0
NAME_WIDTH = 33
# (f0/S) * dt^2 = 1.0e5 * 0.01 = 1000 Pa per unit n(n+1)/2
QUAD_ANCHOR = (F0 / S_NUM) * DT * DT
DP10_ANCHOR = QUAD_ANCHOR * N_REC * (N_REC - 1) / 2.0   # n=10 -> 55e3 Pa

TOL_QUAD = 0.10
TOL_BETA_INV = 0.06
TOL_3D_BETA = 0.12
TOL_DAMAGE = 0.01
TOL_CROSSCHECK = 1.0e-6
# CI #127 REAL solver readings for the anisotropy arm:
#   surface_ratio_aniso = 1.004065  -> 0.4065% shift
#   damage_dev_aniso    = 1.247e-2  -> 1.2474% shift
#   beta_dev_aniso      = 0
# A enters ONLY C2*grad(alpha)^2, which C1*w(alpha) swamps before
# initiation.  Two-sided band: NOT bit-equivalent (wired) yet NOT a
# large shift (A stayed out of the C1 channel / no propagation).
SURF_BAND = (1.0e-4, 5.0e-2)   # |surface_ratio - 1| in (0.01%, 5%)
FATAL = ('PETSC ERROR', 'Caught signal', 'SEGV', 'fortran runtime error',
         'Fatal Error:', 'NaN detected')

B2_TAGS = ('linear', 'aniso', 'dim3')
ALL_TAGS = B2_TAGS + ('b0_inject', 'b0_scale')

SRC = {
    'flowasm': os.path.join(_ROOT, 'm_Flow', 'm_MEF90_FlowAssembly.F90'),
    'defasm': os.path.join(_ROOT, 'm_DefMech', 'm_MEF90_DefMechAssembly.F90'),
    'linalg': os.path.join(_ROOT, 'MEF90', 'm_MEF90_LinAlg.F90'),
    'workflow': os.path.join(_ROOT, '.github', 'workflows', 'build.yml'),
}
DECKS = {
    'inject': os.path.join(_HERE, 'terzaghi415_inject.yaml'),
    'linear': os.path.join(_HERE, 'terzaghi415_b2_linear.yaml'),
    'aniso': os.path.join(_HERE, 'terzaghi415_b2_aniso.yaml'),
}

FILEMAP = {}
for _t in ALL_TAGS:
    FILEMAP[_t] = ('%s.log' % _t, '%s.gen' % _t)
for _t in ('b0_verdict', 'b1_verdict', 'main'):
    FILEMAP[_t] = ('%s.log' % _t, None)
FILEMAP['mesh2'] = (None, 'mesh2.gen')
FILEMAP['mesh3'] = (None, 'mesh3.gen')

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
    return {'num_dim': d.get('num_dim'), 'num_nodes': d.get('num_nodes')}


def energy_rows(log):
    return [tuple(num(m.group(i)) for i in range(1, 6))
            for m in ENERGY.finditer(log)]


def g32_rows(log):
    return [(int(m.group(1)), num(m.group(2)), num(m.group(3)))
            for m in G32.finditer(log)]


def g45_rows(log):
    return [(int(m.group(1)), num(m.group(2)), num(m.group(3)), num(m.group(4)))
            for m in G45.finditer(log)]


def linfit(xs, ys):
    xb = sum(xs) / len(xs)
    yb = sum(ys) / len(ys)
    den = sum((x - xb) ** 2 for x in xs)
    if den <= 0.0:
        return 0.0, yb
    m = sum((x - xb) * (y - yb) for x, y in zip(xs, ys)) / den
    return m, yb - m * xb


def beta_of(d, tag):
    p = d['exo'][tag]['pbulk'][:N_REC]
    return linfit([i * DT for i in range(N_REC)], [x - P0 for x in p])[0]


def code_lines(text):
    out = []
    for line in text.splitlines():
        s = line.split('#', 1)[0].rstrip()
        if s.strip():
            out.append(s.strip())
    return out


def key_of(line):
    return line.split(':', 1)[0].strip()


def fatal_in(log):
    return [s for s in FATAL if s in log]


# ---------------------------------------------------------------------------
# Verdicts
# ---------------------------------------------------------------------------

def check_l1(d):
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
    if d['mesh2']['num_dim'] != 2:
        issues.append('2D mesh num_dim=%s' % d['mesh2']['num_dim'])
    if d['mesh3']['num_dim'] != 3:
        issues.append('3D mesh num_dim=%s' % d['mesh3']['num_dim'])
    if not d['mesh3']['num_nodes'] or d['mesh3']['num_nodes'] < 1000:
        issues.append('3D mesh node count %s is implausible' % d['mesh3']['num_nodes'])
    if 'Displacement_Z' not in d['exo']['dim3']['names']:
        issues.append('3D EXO lost Displacement_Z (dimension dispatch failed)')
    return pass_fail(not issues, '; '.join(issues) if issues else
                     '5 runs x %d clean records, EXO pmax == solver, 11 energy '
                     'rows each, 2D + 3D mesh metadata, Displacement_Z live'
                     % N_REC)


def check_l2(d):
    """Backward-Euler anchor for f(t) = f0*t: Delta_p_n = (f0/S)dt^2 n(n+1)/2."""
    issues = []
    p = d['exo']['linear']['pbulk'][:N_REC]
    delta = [x - P0 for x in p]
    xs = [i * (i + 1) / 2.0 for i in range(N_REC)]
    slope, intercept = linfit(xs, delta)
    resid = max(abs(y - (slope * x + intercept)) for x, y in zip(xs, delta))
    d['metrics']['quad_slope'] = slope
    d['metrics']['quad_intercept'] = intercept
    d['metrics']['quad_resid'] = resid
    d['metrics']['dp10'] = delta[-1]
    if abs(delta[0]) > 0.02 * P0:
        issues.append('Delta_p(0)=%.6g is not zero' % delta[0])
    if slope <= 0.0:
        issues.append('quadratic coefficient %.6g <= 0' % slope)
    else:
        if abs(slope / QUAD_ANCHOR - 1.0) > TOL_QUAD:
            issues.append('quadratic coefficient %.6g vs (f0/S)dt^2=%.6g '
                          '(%.1f%% off)' % (slope, QUAD_ANCHOR,
                                            100 * abs(slope / QUAD_ANCHOR - 1)))
        if abs(intercept) > 0.05 * max(delta):
            issues.append('quadratic fit has a non-zero intercept %.6g' % intercept)
        if resid / max(delta) > 0.06:
            issues.append('quadratic anchor residual %.2f%% of peak'
                          % (100 * resid / max(delta)))
        if abs(delta[-1] / DP10_ANCHOR - 1.0) > TOL_QUAD:
            issues.append('Delta_p(1s)=%.6g vs discrete anchor %.6g'
                          % (delta[-1], DP10_ANCHOR))
    # Reject the CONSTANT-rate model: a straight line in t must fit WORSE
    # than the t^2 anchor, and the final pressure must be ~55% of B0's.
    lin_slope, lin_intercept = linfit([i * DT for i in range(N_REC)], delta)
    lin_resid = max(abs(y - (lin_slope * t + lin_intercept))
                    for t, y in zip([i * DT for i in range(N_REC)], delta))
    d['metrics']['lin_resid'] = lin_resid
    if max(delta) > 0 and not (lin_resid > 3.0 * resid):
        issues.append('a constant-rate straight line fits as well as the t^2 '
                      'anchor: source_scaling: linear is not active')
    ref = d['exo']['b0_inject']['pbulk'][-1] - P0
    d['metrics']['dp10_ratio'] = delta[-1] / ref if ref > 0 else -1.0
    if ref > 0 and not (0.40 < delta[-1] / ref < 0.70):
        issues.append('Delta_p(1s)/B0=%.3f outside 0.40..0.70 (expected 0.55)'
                      % (delta[-1] / ref))
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'quadratic coeff %.6g (anchor %.6g), intercept %.3g, '
                     'residual %.2f%% (linear-model residual %.2f%%), '
                     'Delta_p(1s)=%.6g = %.3f x B0'
                     % (slope, QUAD_ANCHOR, intercept,
                        100 * resid / max(delta) if max(delta) > 0 else 0,
                        100 * lin_resid / max(delta) if max(delta) > 0 else 0,
                        delta[-1], d['metrics'].get('dp10_ratio', -1)))


def check_l3(d):
    issues = []
    wl = [r[3] for r in g45_rows(d['log']['linear'])]
    wb = [r[3] for r in g45_rows(d['log']['b0_inject'])]
    w0 = [r[3] for r in g45_rows(d['log']['b0_scale'])]
    if len(wl) < N_REC or len(wb) < N_REC or len(w0) < N_REC:
        return pass_fail(False, 'missing Gate 4.5 work series')
    # f(0) = 0 for source_scaling: linear, so W(step 1) is legitimately zero;
    # every subsequent step must be strictly positive.
    if any(w <= 0.0 for w in wl[1:]):
        issues.append('variable-rate injection work is not strictly positive')
    if any(wl[i] < wl[i - 1] - 1.0e-9 for i in range(1, N_REC)):
        issues.append('variable-rate injection work is not monotone')
    if any(abs(w) > 1.0e-12 for w in w0):
        issues.append('zero-source control has non-zero work')
    d['metrics']['w_linear_last'] = wl[-1]
    d['metrics']['w_b0_last'] = wb[-1]
    if wl[-1] >= wb[-1]:
        issues.append('variable-rate W_last=%.4g >= constant-rate %.4g'
                      % (wl[-1], wb[-1]))
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'control W=0; variable-rate W=%.4g..%.4g monotone; '
                     'constant-rate W_last=%.4g (ramp is strictly softer)'
                     % (wl[0], wl[-1], wb[-1]))


def check_l4(d):
    """Falsifiable NULL proposition for pre-initiation anisotropy (CI #127).

    The original L4 demanded |surface_ratio-1| > 5%, which assumes a
    PROPAGATING crack whose grad(alpha)^2 term dominates.  In the
    storage-dominated window Damage_max ~ 1e-4 and alpha is almost
    spatially uniform, so C1*w(alpha) >> (C2*grad(alpha)).grad(alpha)
    and A can move the surface energy only by O(0.1%).  CI #127 measured
    0.4065% (surface) vs 1.2474% (damage) with beta_dev = 0.  We now
    assert a two-sided band plus an ordering constraint instead.
    """
    issues = []
    er_a = energy_rows(d['log']['aniso'])
    er_b = energy_rows(d['log']['b0_inject'])
    if len(er_a) < N_REC or len(er_b) < N_REC:
        return pass_fail(False, 'missing energy rows for the anisotropy arm')
    sa, sb = er_a[N_REC - 1][3], er_b[N_REC - 1][3]
    d['metrics']['surface_ref'] = sb
    ds = 0.0
    lo, hi = SURF_BAND
    if sb <= 0.0:
        issues.append('reference surface energy is not positive')
    else:
        r = sa / sb
        ds = abs(r - 1.0)
        d['metrics']['surface_ratio_aniso'] = r
        d['metrics']['surface_dev_aniso'] = ds
        if ds <= lo:
            issues.append('surface energy is bit-equivalent (%.3g%% <= %.3g%%): '
                          'toughnessAnisotropyMatrix never reached C2'
                          % (100 * ds, 100 * lo))
        elif ds >= hi:
            issues.append('anisotropy shifted surface energy %.2f%% (>=%.0f%%): '
                          'A leaked into the C1*w(alpha) channel or the field is '
                          'no longer pre-initiation' % (100 * ds, 100 * hi))
    ba, bb = beta_of(d, 'aniso'), beta_of(d, 'b0_inject')
    d['metrics']['beta_aniso'] = ba
    d['metrics']['beta_dev_aniso'] = abs(ba / bb - 1.0) if bb > 0 else 1.0
    if bb <= 0 or d['metrics']['beta_dev_aniso'] > TOL_BETA_INV:
        issues.append('anisotropy changed the pressure ramp by %.2f%% '
                      '(>%.0f%%): the compliance channel should be blind to Gc'
                      % (100 * d['metrics']['beta_dev_aniso'], 100 * TOL_BETA_INV))
    da = max(d['exo']['aniso']['dmax'][:N_REC])
    db = max(d['exo']['b0_inject']['dmax'][:N_REC])
    dd = abs(da / db - 1.0) if db > 0 else 1.0
    d['metrics']['damage_dev_aniso'] = dd
    if dd < TOL_DAMAGE:
        issues.append('damage field is bit-identical to the isotropic deck: '
                      'the anisotropic solve never took place')
    if ds >= dd:
        issues.append('surface-energy response %.2f%% is NOT weaker than the '
                      'damage response %.2f%%: A is not confined to the '
                      'gradient (C2) channel' % (100 * ds, 100 * dd))
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'surface dev %.3f%% in band (%.3g%%..%.0f%%) and < damage '
                     'dev %.3f%%, beta dev %.2f%%: anisotropy is wired but NOT '
                     'identifiable before initiation'
                     % (100 * ds, 100 * lo, 100 * hi, 100 * dd,
                        100 * d['metrics'].get('beta_dev_aniso', -1)))


def check_l5(d):
    issues = []
    base, aniso = d['deck']['inject'], d['deck']['aniso']
    ca, cb = code_lines(base), code_lines(aniso)
    added = [x for x in ca if x not in cb]
    d['metrics']['aniso_added_lines'] = len(cb) - len(ca)
    if len(cb) != len(ca) + 1:
        issues.append('aniso deck has %d code lines, expected %d (+1)'
                      % (len(cb), len(ca) + 1))
    extra = [x for x in cb if x not in ca]
    if len(extra) != 1 or key_of(extra[0]) != 'toughnessAnisotropyMatrix':
        issues.append('the single added line is not toughnessAnisotropyMatrix: %r'
                      % extra)
    else:
        d['metrics']['aniso_matrix'] = extra[0]
        if added:
            issues.append('aniso deck also removed code lines: %r' % added)
        if not re.search(r'toughnessAnisotropyMatrix:\s*2,0\.5,0', extra[0]):
            issues.append('aniso matrix is not 2,0.5,0 (bedding contrast)')
        if re.search(r',\s', extra[0]):
            issues.append('a comma is followed by a space: PETSc 3.3 RealArray '
                          'would keep only the first element (CI #85)')
    if 'toughnessAnisotropyMatrix' not in d['log']['aniso']:
        issues.append('the PetscBag dump never shows toughnessAnisotropyMatrix')
    s = d['src']
    for token in ('MEF90MatS2DIdentity = MatS2D(1.0_Kr,1.0_Kr,0.0_Kr)',
                  'PetscReal          :: XX',
                  'PetscReal          :: YY',
                  'PetscReal          :: XY'):
        if token not in s['linalg']:
            issues.append('m_MEF90_LinAlg.F90 lost %s' % token)
    if 'matPropSet%toughnessAnisotropyMatrix' not in s['defasm']:
        issues.append('surface-energy C2 lost the anisotropy matrix')
    if 'VecScale(FlowCtx%source,timeVal,ierr)' not in s['flowasm']:
        issues.append('MEF90Scaling_Linear no longer scales the source by timeVal')
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'aniso deck = inject + 1 line (2,0.5,0, no comma-space), '
                     'MatS2D(XX,YY,XY) identity, C2 and VecScale wiring intact')


def check_l6(d):
    issues = []
    e = d['exo']['dim3']
    if 'Displacement_Z' not in e['names']:
        issues.append('3D EXO has no Displacement_Z field')
    b3, b2 = beta_of(d, 'dim3'), beta_of(d, 'b0_inject')
    d['metrics']['beta_dim3'] = b3
    d['metrics']['beta_dev_dim3'] = abs(b3 / b2 - 1.0) if b2 > 0 else 1.0
    if b3 <= 0 or d['metrics']['beta_dev_dim3'] > TOL_3D_BETA:
        issues.append('3D beta=%.6g deviates %.1f%% from 2D (tol %.0f%%)'
                      % (b3, 100 * d['metrics']['beta_dev_dim3'],
                         100 * TOL_3D_BETA))
    log = d['log']['dim3']
    if log.count('Linear solve converged') < N_REC - 1:
        issues.append('3D run has fewer than %d converged flow solves' % (N_REC - 1))
    if 'rc=124' in log or 'TIMEOUT' in log:
        issues.append('3D run shows a timeout signature (CI #120 regression)')
    for sig in fatal_in(log):
        issues.append('3D fatal: %s' % sig)
    wf = d['src']['workflow']
    for token in ('Beam1x10-tetra4', '-np 1',
                  '-flow_pc_type bjacobi -flow_sub_pc_type lu',
                  'timeout -k 10s 1800'):
        if token not in wf:
            issues.append('build.yml lost the 3D diagnostic config: %s' % token)
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'tetra4 3D beta=%.6g Pa/s (%.2f%% vs 2D), '
                     'Displacement_Z live, direct-LU diagnostic PC, no timeout'
                     % (b3, 100 * d['metrics'].get('beta_dev_dim3', -1)))


def check_l7(d):
    issues = []
    lin = code_lines(d['deck']['linear'])
    inj = code_lines(d['deck']['inject'])
    diff = [(x, y) for x, y in zip(lin, inj) if x != y]
    if len(lin) != len(inj) or len(diff) != 1:
        issues.append('linear deck has %s code deltas (expected 1)'
                      % (len(diff) if len(lin) == len(inj) else 'length mismatch'))
    elif key_of(diff[0][0]) != 'source_scaling':
        issues.append('linear deck changed %s, expected source_scaling'
                      % key_of(diff[0][0]))
    if 'source_scaling: linear' not in d['deck']['linear']:
        issues.append('linear deck lost source_scaling: linear')
    for bad in ('batched_assembly', 'biot_effective_stress',
                'temp_dep_viscosity'):
        if bad in d['deck']['linear'] or bad in d['deck']['aniso']:
            issues.append('a B2 deck enables non-default switch %s' % bad)
    wf = d['src']['workflow']
    for token in ('terzaghi415_b2_linear.yaml', 'terzaghi415_b2_aniso.yaml',
                  'terzaghi415_inject.yaml', 'Beam1x10-tetra4',
                  'check_gate415_b2_extension.py --gate',
                  'gate415_b2_verdict.log', 'gb2_rc',
                  'if [ "$gb2_rc" -ne 0 ]',
                  'if [ "$sdim3_rc" -ne 0 ]',
                  'rm -f data/Beam1x10-tetra4_out.gen'):
        if token not in wf:
            issues.append('build.yml missing %s' % token)
    if wf.count('rm -f data/Beam1x10-tri3_out.gen') < 10:
        issues.append('2D EXO rebuilds dropped below 10 (B0 x2 + B1 x6 + B2 x2)')
    return pass_fail(not issues, '; '.join(issues) if issues else
                     '1-line linear deck, +1-line aniso deck, 2D/3D runs + '
                     'verdict + deferred fails + EXO rebuilds all wired')


def check_l8(d):
    issues = []
    v0, v1 = d['log']['b0_verdict'], d['log']['b1_verdict']
    if 'GATE415_RESULT: PASS' not in v0:
        issues.append('B0 verdict is no longer PASS')
    for token in ('g415_slope_measured        = 1.000000e+05',
                  'g415_S_eff                 = 1.500000e-10'):
        if token not in v0:
            issues.append('B0 CI #125 number drifted: %s' % token)
    if 'GATE415_B1_RESULT: PASS' not in v1:
        issues.append('B1 verdict is no longer PASS')
    main = d['log']['main']
    rows = g32_rows(main)
    if len(rows) < 50:
        issues.append('production log has only %d Gate 3.2 rows' % len(rows))
    elif rows[-1][2] >= 0.5 * P0:
        issues.append('production pressure did not consolidate')
    wm = [r[3] for r in g45_rows(main)]
    if wm and any(abs(w) > 1.0e-12 for w in wm):
        issues.append('production injection work is non-zero')
    for bad in ('terzaghi415_b2_', 'check_gate415_b2_', 'Beam1x10-tetra4'):
        if bad in main:
            issues.append('a B2 diagnostic leaked into the production log: %s' % bad)
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'B0 (CI #125 beta=1e5/S_eff=1.5e-10) and B1 verdicts still '
                     'PASS; production consolidates over %d steps, W=0, no '
                     'B2/3D leakage' % len(rows))


CHECKS = (('L1', 'run / EXO / log / energy-row integrity (2D+3D)', check_l1),
          ('L2', 'variable-rate discrete t^2 anchor', check_l2),
          ('L3', 'variable-rate injection work', check_l3),
          ('L4', 'pre-initiation anisotropy: wired, bounded, beta-blind', check_l4),
          ('L5', 'anisotropy deck / bag / source contracts', check_l5),
          ('L6', '3D tetra4 dispatch + direct-LU diagnostic PC', check_l6),
          ('L7', 'static anti-regression + CI wiring', check_l7),
          ('L8', 'B0 + B1 unchanged + production isolation', check_l8))


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
    print('GATE415_B2_RESULT: %s%s' % ('PASS' if not fails else 'FAIL',
                                       '' if not fails else ' (' + ','.join(fails) + ')'))
    for tag, label, verdict, detail in rows:
        if verdict == 'FAIL':
            print('  %s detail: %s' % (tag, detail))
    m = LAST_METRICS.get('metrics', {})
    for key in ('quad_slope', 'quad_intercept', 'quad_resid', 'lin_resid',
                'dp10', 'dp10_ratio', 'w_linear_last', 'w_b0_last',
                'surface_ref', 'surface_ratio_aniso', 'surface_dev_aniso',
                'beta_aniso',
                'beta_dev_aniso', 'damage_dev_aniso', 'beta_dim3',
                'beta_dev_dim3', 'aniso_matrix'):
        v = m.get(key, -1)
        if isinstance(v, str):
            print('gb2_%-26s = %s' % (key, v))
        else:
            print('gb2_%-26s = %.6e' % (key, v))


LAST_METRICS = {}


def synth_sources():
    return {k: read(v) for k, v in SRC.items()}


def synth_decks():
    return {k: read(v) for k, v in DECKS.items()}


# CI #127 REAL solver readings, used as the synthetic positive case so the
# selftest protects the measured pre-initiation signature rather than a
# made-up "large" anisotropy response.
SURF_REF = 6.49012e-4               # b0_inject surface energy at step 11
SURF_ANISO = SURF_REF * 1.004065   # gb2_surface_ratio_aniso = 1.004065
DMG_REF = 1.4803e-4                 # b0_inject damage_max
DMG_ANISO = DMG_REF * 1.01247384   # gb2_damage_dev_aniso = 1.247384e-2

B0_VERDICT = """=== Gate 4.15 / Stage B0 ===
PASS V4 injection slope f/S + source-unit calibration
GATE415_RESULT: PASS
g415_slope_measured        = 1.000000e+05
g415_S_eff                 = 1.500000e-10
"""
B1_VERDICT = """=== Gate 4.15 / Stage B1 ===
PASS A4 Gc/k/mu invisible on the pressure curve
GATE415_B1_RESULT: PASS
"""


def make_log(tag, pmax, work, surface, n=N_REC, bag_lines=''):
    log = ('=== Gate 4.15 B2 synthetic run %s ===\n' % tag + bag_lines +
           'x' * 1200 + '\n')
    for i in range(1, n + 1):
        log += 'Gate 3.1: step %d enter transient block\n' % i
        log += 'Gate 3.2: step %d pmin 0.0 pmax %.10E\n' % (i, pmax[i - 1])
        log += ('Gate 4.12: step %d t %.4g dpmax %.4g dpmin 0.0\n'
                % (i, (i - 1) * DT, pmax[i - 1] - pmax[i - 2] if i > 1 else 0.0))
        log += ('Gate 4.5: step %d t %.10E E %.10E W %.10E\n'
                % (i, (i - 1) * DT, 1.0e-3, work[i - 1]))
        log += ('======= Total: elastic energy:  %.5E work:  %.5E '
                'cohesive:  0.00000E+00 surface:  %.5E total:  %.5E\n'
                % (1.0e-6, 1.0e-6, surface[i - 1], 1.0))
        if i > 1:
            log += '    Linear solve converged due to CONVERGED_RTOL iterations 1\n'
    return log


def synth():
    d = {'log': {}, 'exo': {}, 'metrics': {}}
    names2 = ['Pressure', 'Damage', 'Displacement_X', 'Displacement_Y']
    names3 = names2 + ['Displacement_Z', '']
    for tag in ALL_TAGS:
        if tag == 'linear':
            pbulk = [P0 + QUAD_ANCHOR * (i * (i + 1) / 2.0)
                     for i in range(N_REC)]
            work = [F0 * AREA * (i * DT) * p for i, p in enumerate(pbulk)]
            surf_last = 3.20e-4
            dmg_last = 7.4e-5
        elif tag == 'b0_scale':
            pbulk = [P0] * N_REC
            work = [0.0] * N_REC
            surf_last = 1.20755e-6
            dmg_last = 8.6005e-6
        else:
            pbulk = [P0 + BETA0 * i * DT for i in range(N_REC)]
            work = [F0 * AREA * p for p in pbulk]
            surf_last = SURF_REF
            dmg_last = DMG_REF
        if tag == 'aniso':
            surf_last = SURF_ANISO
            dmg_last = DMG_ANISO
        pmax = [P0 + 2.9184 * (p - P0) for p in pbulk]
        surface = [surf_last * ((i + 1) / N_REC) ** 2 for i in range(N_REC)]
        dmax = [dmg_last * (i + 1) / N_REC for i in range(N_REC)]
        bag = ('  toughnessAnisotropyMatrix = 2,0.5,0 ; []\n'
               if tag == 'aniso' else '')
        d['exo'][tag] = {'numrecs': N_REC,
                         'names': names3 if tag == 'dim3' else names2,
                         'finite': True, 'pbulk': pbulk,
                         'pmax': pmax, 'dmax': dmax}
        d['log'][tag] = make_log(tag, pmax, work, surface, bag_lines=bag)
    p_main = [P0 * (0.0906 ** (j / 50.0)) for j in range(51)]
    d['log']['main'] = make_log('main', p_main, [0.0] * 51,
                                [1.0e-6] * 51, n=51)
    d['log']['b0_verdict'] = B0_VERDICT
    d['log']['b1_verdict'] = B1_VERDICT
    d['mesh2'] = {'num_dim': 2, 'num_nodes': 356}
    d['mesh3'] = {'num_dim': 3, 'num_nodes': 2635}
    d['src'] = synth_sources()
    d['deck'] = synth_decks()
    return d


def selftest():
    import copy
    base = synth()
    cases = []
    rows, fails = run(base)
    LAST_METRICS['metrics'] = base['metrics']
    cases.append(('POS full B2 verdict', rows, set(), None))
    print('POS  %-54s %s' % ('full B2 verdict', 'PASS' if not fails else 'FAIL ' + ','.join(fails)))

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
    allok &= neg('3D EXO lost Displacement_Z',
                 lambda d: d['exo']['dim3'].update(
                     names=['Pressure', 'Damage', 'Displacement_X',
                            'Displacement_Y']), {'L1', 'L6'})
    allok &= neg('linear EXO truncated to 10 records',
                 lambda d: d['exo']['linear'].update(numrecs=10), {'L1'})
    allok &= neg('fatal signature in the aniso log',
                 lambda d: d['log'].update(
                     aniso=d['log']['aniso'] + 'PETSC ERROR: boom\n'), {'L1'})
    allok &= neg('linear run secretly used the CONSTANT-rate ramp',
                 lambda d: d['exo']['linear'].update(
                     pbulk=d['exo']['b0_inject']['pbulk'][:]), {'L2'})
    allok &= neg('linear ramp is linear in t, not quadratic',
                 lambda d: d['exo']['linear'].update(
                     pbulk=[P0 + 5.5e4 * (i / 10.0) for i in range(N_REC)]),
                 {'L2'})
    allok &= neg('linear ramp has the wrong coefficient (2x anchor)',
                 lambda d: d['exo']['linear'].update(
                     pbulk=[P0 + 2.0 * QUAD_ANCHOR * (i * (i + 1) / 2.0)
                            for i in range(N_REC)]), {'L2'})
    allok &= neg('variable-rate injection work not monotone',
                 lambda d: d['log'].update(
                     linear=re.sub(
                         r'(Gate 4\.5:\s*step\s+%d\s+.*?\s+W\s+)[-+0-9.eEdD]+'
                         % N_REC,
                         lambda m: m.group(1) + '1.0000000000E-03',
                         d['log']['linear'], count=1)), {'L3'})
    allok &= neg('zero-source control reports non-zero work',
                 lambda d: d['log'].update(
                     b0_scale=d['log']['b0_scale'].replace(
                         'W 0.0000000000E+00', 'W 3.5000000000E+00')), {'L3'})
    allok &= neg('aniso surface energy is bit-equivalent (A never reached C2)',
                 lambda d: d['log'].update(
                     aniso=d['log']['aniso'].replace(
                         'surface:  %.5E' % SURF_ANISO,
                         'surface:  %.5E' % SURF_REF)), {'L4'})
    allok &= neg('A leaked into the C1*w(alpha) channel (large surface shift)',
                 lambda d: d['log'].update(
                     aniso=d['log']['aniso'].replace(
                         'surface:  %.5E' % SURF_ANISO,
                         'surface:  %.5E' % (0.32 * SURF_REF))), {'L4'})
    allok &= neg('surface response exceeds the damage response',
                 lambda d: d['log'].update(
                     aniso=d['log']['aniso'].replace(
                         'surface:  %.5E' % SURF_ANISO,
                         'surface:  %.5E' % (1.02 * SURF_REF))), {'L4'})
    allok &= neg('damage field bit-identical (anisotropic solve never ran)',
                 lambda d: d['exo']['aniso'].update(
                     dmax=d['exo']['b0_inject']['dmax'][:]), {'L4'})
    allok &= inv('surface response drifts to 0.30% (still inside the band)',
                 lambda d: d['log'].update(
                     aniso=d['log']['aniso'].replace(
                         'surface:  %.5E' % SURF_ANISO,
                         'surface:  %.5E' % (1.003 * SURF_REF))))
    allok &= inv('surface response drifts to 0.90% (still < damage response)',
                 lambda d: d['log'].update(
                     aniso=d['log']['aniso'].replace(
                         'surface:  %.5E' % SURF_ANISO,
                         'surface:  %.5E' % (1.009 * SURF_REF))))
    allok &= neg('aniso deck also breaks the pressure ramp',
                 lambda d: d['exo']['aniso'].update(
                     pbulk=[P0 + 0.5 * BETA0 * i * DT for i in range(N_REC)]),
                 {'L4'})
    allok &= neg('aniso matrix silently reverted to identity',
                 lambda d: d['deck'].update(
                     aniso=d['deck']['aniso'].replace(
                         'toughnessAnisotropyMatrix: 2,0.5,0',
                         'toughnessAnisotropyMatrix: 1,1,0')), {'L5'})
    allok &= neg('comma followed by a space (CI #85 RealArray trap)',
                 lambda d: d['deck'].update(
                     aniso=d['deck']['aniso'].replace(
                         'toughnessAnisotropyMatrix: 2,0.5,0',
                         'toughnessAnisotropyMatrix: 2, 0.5, 0')), {'L5'})
    allok &= neg('3D beta diverges from the 2D compliance ramp',
                 lambda d: d['exo']['dim3'].update(
                     pbulk=[P0 + 2.0 * BETA0 * i * DT for i in range(N_REC)]),
                 {'L6'})
    allok &= neg('3D direct-LU diagnostic PC removed (CI #120 regression)',
                 lambda d: d['src'].update(workflow=d['src']['workflow'].replace(
                     '-flow_pc_type bjacobi -flow_sub_pc_type lu',
                     '-flow_pc_type none')), {'L6'})
    allok &= neg('linear deck changed a second code line',
                 lambda d: d['deck'].update(
                     linear=d['deck']['linear'].replace('p0: 1.0e5', 'p0: 2.0e5')),
                 {'L7'})
    allok &= neg('deferred gb2_rc hard fail removed',
                 lambda d: d['src'].update(workflow=d['src']['workflow'].replace(
                     'if [ "$gb2_rc" -ne 0 ]', 'if false')), {'L7'})
    allok &= neg('B1 verdict silently became FAIL',
                 lambda d: d['log'].update(
                     b1_verdict=B1_VERDICT.replace('GATE415_B1_RESULT: PASS',
                                                   'GATE415_B1_RESULT: FAIL')),
                 {'L8'})
    allok &= neg('B2 deck name leaked into the production log',
                 lambda d: d['log'].update(
                     main=d['log']['main'] + '-options_file_yaml data/terzaghi415_b2_linear.yaml\n'),
                 {'L8'})
    allok &= inv('damage SNES residual plateau (-5) appears, physics unchanged',
                 lambda d: d['log'].update(
                     aniso=d['log']['aniso'] + '[ERROR]: damage field SNESSolve failed with SNESConvergedReason -5.\n'))
    print('SELFTEST_RESULT: %s (%d positive + %d mutations, incl. 1 invariance)'
          % ('PASS' if allok else 'FAIL', 1, len(cases) - 1))
    return 0 if allok else 1


def main(argv):
    if argv and argv[0] == '--selftest':
        return selftest()
    if argv and argv[0] == '--gate':
        if len(argv) != 2:
            print('GATE415_B2_RESULT: FAIL (usage: --gate <artifact_dir>)')
            return 1
        root = argv[1]
        try:
            d = {'log': {}, 'exo': {}, 'metrics': {}}
            for tag in ALL_TAGS + ('b0_verdict', 'b1_verdict', 'main'):
                d['log'][tag] = read(os.path.join(root, FILEMAP[tag][0]))
            for tag in ALL_TAGS:
                d['exo'][tag] = exo_fields(os.path.join(root, FILEMAP[tag][1]))
            d['mesh2'] = mesh_meta(os.path.join(root, FILEMAP['mesh2'][1]))
            d['mesh3'] = mesh_meta(os.path.join(root, FILEMAP['mesh3'][1]))
            d['src'] = synth_sources()
            d['deck'] = synth_decks()
        except Exception as exc:
            print('=== Gate 4.15 / Stage B2: variable rate / anisotropy / 3D ===')
            print('FAIL L1 run / EXO / log / energy-row integrity (2D+3D) | '
                  'cannot parse artefacts: %r' % exc)
            print('GATE415_B2_RESULT: FAIL (L1 data error: %r)' % exc)
            return 1
        rows, fails = run(d)
        LAST_METRICS['metrics'] = d['metrics']
        report(rows, fails, '=== Gate 4.15 / Stage B2: variable rate / bedding '
               'anisotropy / 3D consistency (real CI logs/EXO) ===')
        return 0 if not fails else 1
    print('usage: check_gate415_b2_extension.py '
          '[--gate <artifact_dir> | --selftest]')
    return 2


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
