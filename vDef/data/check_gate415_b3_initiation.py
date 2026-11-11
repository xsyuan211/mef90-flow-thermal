#!/usr/bin/env python3
"""Gate 4.15 / Stage B3 -- FRACTURE-INITIATION ACCESS hard verdict.

Why this gate exists
--------------------
B0 (CI #125), B1 (CI #127) and B2 (CI #128) are strictly PRE-initiation:
the real solver gives Damage_max = 1.480310e-4 << 1, T = 1.5667e-7 and
beta = f/S = 1.0e5 Pa/s.  Every published conclusion of stages A/B0/B1/B2
is therefore limited to the storage-dominated window, where the Lecampion
et al. (2017) compliance law holds and rock parameters are invisible on the
pressure curve.  Santillan (2018) rate thresholds, Detournay (2016) M/K
vertices, K/C regime maps and any "rate optimization table" have NO
definition in that window (roadmap sec. 5.4 decision 3, sec. 14).

B3 changes the PHYSICAL PREMISE: it asks whether this numerical platform can
reach an honest initiation window at all, and where the threshold is.  It
does so with ZERO Fortran changes by exploiting the AT2 damage residual.

AT2 contract (verified in the live sources)
    m_MEF90_DefMechAT2.F90:33        cw = 0.5
    m_MEF90_DefMechAT2.F90:53,72     a(alpha)  = (1-alpha)^2,  a' = -2(1-alpha)
    m_MEF90_DefMechAT2.F90:105,124   w(alpha)  = alpha^2,      w' = 2 alpha
    m_MEF90_DefMechAssembly.F90:1687 C1 = Gc/(4*cw*ell) = Gc/(2*ell)
    m_MEF90_DefMechAssembly.F90:1688 C2 = Gc*ell/(2*cw)*A
    m_MEF90_DefMechAssembly.F90:1786 C3 = a'(alpha)*Psi + C1*w'(alpha)

Homogeneous (gradient-free) equilibrium C3 = 0:

    -2(1-alpha)*Psi + (Gc/ell)*alpha = 0
    => alpha_eq(Gc) = 2 Psi / (Gc/ell + 2 Psi).

Inverting the REAL CI #125/127 B0 reading (EXO vals_nod_var2, record 10):

    alpha(B0) = 1.480310e-4,  Gc0 = 27000 N/m, ell = 0.2 m
    2 Psi     = 1.998714e+01 Pa,  Gc0/ell = 1.35e+05 Pa,
    Gc*(alpha = 1 - 1/sqrt(2)) = ell*2Psi*(1+sqrt(2)) = 9.650647 N/m,
    Gc*(alpha = 0.5)           = ell*2Psi             = 3.997429 N/m.

This is NOT an unchecked extrapolation.  CI #127 measured
    gb1_surface_ratio_gc = 2.000410   (Gc halved -> surface energy x2.000410)
while the same analytic anchor predicts 1.999408 (0.05% error), and
    gb1_beta_dev_gc = 0               (Gc invisible on the pressure curve).
So alpha ~ 1/Gc and surface ~ 1/Gc are already verified on the real solver
in the low-damage branch.  B3 extends that ladder by 4 decades of Gc.

Ladder (each deck is EXACTLY ONE code-line delta of terzaghi415_inject.yaml,
only cs0001.fracturetoughness; NO Fortran, NO Flow key, NO time-grid change):

    b0_inject  Gc = 27000   alpha_eq = 1.480310e-04   (existing B0 run)
    gc_half    Gc = 13500   alpha_eq = 2.960182e-04   (existing B1 gc run)
    gc2700     Gc =  2700   alpha_eq = 1.478340e-03
    gc270      Gc =   270   alpha_eq = 1.458929e-02
    gc27       Gc =    27   alpha_eq = 1.289600e-01
    gc2p7      Gc =   2.7   alpha_eq = 5.968602e-01   (crosses initiation)

INITIATION is defined without an arbitrary number: a(alpha)=(1-alpha)^2 is
the stiffness degradation, so a = 1/2 at

    alpha_stiff_half = 1 - 1/sqrt(2) = 0.2928932188.

Crossing it means the material has lost half its load-bearing stiffness.
N4 brackets the measured threshold between the last non-initiated and the
first initiated rung and compares it with the analytic Gc*.

Verdicts
--------
N1 run / EXO / log / energy-row integrity (+ EXO vs Gate 3.2 cross-check)
N2 Flow invariance: beta = f/S survives the whole Gc ladder
N3 AT2 analytic anchor: THREE independent implementations (closed form,
   algebraic rearrangement, numerical root of C3) + low-damage 1/Gc power
   law + independent CI #127 surface-ratio calibration
N4 initiation access: >=1 rung crosses alpha_stiff_half, damage monotone in
   1/Gc, >=3 decades of growth, measured Gc* bracket vs analytic Gc*
N5 localisation & solver hygiene: damage spatially non-uniform (not a
   uniformly softened blob), bounded in [0,1], irreversible, -5 reported
N6 energy consistency: surface ~ 1/Gc in the low-damage branch, W invariant
N7 static anti-regression: one-line deck contract, source C1/C2/C3/AT2
   contract, EXO rebuilds, gb3_rc deferred hard fail
N8 state discipline + zero production perturbation: B0/B1/B2 verdicts PASS,
   no B3 deck name in the production log, no propagation-regime claim

STATE DISCIPLINE: B3 establishes an initiation threshold, NOT a propagating
fracture.  The damage SNES can still report reason -5 (CI #84 residual
plateau); energy rows are deterministic integrals of the CURRENT field and
are used only as relative, falsifiable observables.  NO claim is made about
fracture length, fracture-network area, an absolute breakdown pressure, a
propagation-regime scaling or a field rate optimization table.

Pure standard library (CI runner has no numpy / scipy / PyYAML).
  --gate <dir>  read real CI artefacts from one directory
  --selftest    1 positive case + >=24 mutations (incl. 3 invariances)
  --predict     print the offline ladder table (no artefacts needed)
"""

import copy
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse_cdf2 import CDF2  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, '..', '..'))

# ---------------------------------------------------------------------------
# Frozen constants.  Everything with a "REAL" comment comes from a CI solver
# run, not from an offline guess.
# ---------------------------------------------------------------------------
N_REC = 11
DT = 0.1
P0 = 1.0e5
F0 = 1.5e-5
S_NUM = 1.50e-10
BETA0 = F0 / S_NUM
AREA = 10.0
NAME_WIDTH = 33

ELL = 0.2                     # internalLength [m], all B0/B1/B2/B3 decks
CW = 0.5                      # m_MEF90_DefMechAT2.F90:33
GC0 = 27000.0                 # cs0001.fracturetoughness of B0 [N/m]
ALPHA0_B0 = 1.480310e-4       # REAL CI #125 EXO vals_nod_var2 rec10 max
TWO_PSI_FROZEN = 1.998714e+01  # REAL inversion of ALPHA0_B0 at Gc0/ell
SURF0 = 6.490120e-4           # REAL CI #125 last "======= Total: surface:"
W_LAST_REAL = 30.1605         # REAL CI #125 last Gate 4.5 W (integral of f*p)
SURF_RATIO_GC_REAL = 2.000410  # REAL CI #127 gb1_surface_ratio_gc
BETA_DEV_GC_REAL = 0.0         # REAL CI #127 gb1_beta_dev_gc
# The synthetic W must use the solver's integral calculus, not the naive
# F0*Area*pbulk estimate (which is 0.53% low because W integrates the full
# pressure field, including the local overshoot near the closed end).
W_CALIB = W_LAST_REAL / (F0 * AREA * (P0 + BETA0 * (N_REC - 1) * DT))

TWO_PSI = ALPHA0_B0 * (GC0 / ELL) / (1.0 - ALPHA0_B0)
ALPHA_STIFF_HALF = 1.0 - 1.0 / math.sqrt(2.0)
# alpha_eq(Gc) = 2Psi/(Gc/ell + 2Psi)  =>  Gc(alpha) = ell*2Psi*(1-alpha)/alpha
GC_STAR_HALF = TWO_PSI * ELL                                    # alpha = 0.5
GC_STAR_STIFF = TWO_PSI * ELL * \
    (1.0 - ALPHA_STIFF_HALF) / ALPHA_STIFF_HALF                # a(alpha) = 1/2

TOL_BETA_B3 = 5.0e-3          # k(v) feedback near alpha~1 may move beta
TOL_BETA_STRICT = 1.0e-6      # low-damage rungs must be bitwise invariant
TOL_CROSSCHECK = 1.0e-6
TOL_ANCHOR = 1.0e-10          # three independent alpha_eq implementations
TOL_RATIO_DECADES = 0.30      # low-damage 1/Gc power law: +/- factor 2
TOL_SURF_CALIB = 1.0e-2       # analytic surface ratio vs CI #127 real
TOL_W = 5.0e-3
TOL_GCSTAR_DECADES = 1.0      # measured bracket within 1 decade of Gc*
POWER_BAND_DMG = (-1.30, -0.70)
POWER_BAND_SURF = (-1.45, -0.55)
LOCAL_MIN = 1.5               # dmax/dmean must exceed this (REAL B0 = 7.00)

LADDER = (('b0_inject', GC0),
          ('gc_half', 13500.0),
          ('gc2700', 2700.0),
          ('gc270', 270.0),
          ('gc27', 27.0),
          ('gc2p7', 2.7))
GC_OF = dict(LADDER)
ALL_TAGS = tuple(t for t, _ in LADDER)
NEW_TAGS = ('gc2700', 'gc270', 'gc27', 'gc2p7')
LOW_TAGS = ('b0_inject', 'gc_half', 'gc2700', 'gc270')
STRICT_TAGS = ('gc_half', 'gc2700', 'gc270')

SRC = {
    'defasm': os.path.join(_ROOT, 'm_DefMech', 'm_MEF90_DefMechAssembly.F90'),
    'at2': os.path.join(_ROOT, 'm_DefMech', 'm_MEF90_DefMechAT2.F90'),
    'workflow': os.path.join(_ROOT, '.github', 'workflows', 'build.yml'),
}
DECKS = {'inject': os.path.join(_HERE, 'terzaghi415_inject.yaml'),
         'b1_gc': os.path.join(_HERE, 'terzaghi415_b1_gc.yaml')}
for _t in NEW_TAGS:
    DECKS[_t] = os.path.join(_HERE, 'terzaghi415_b3_%s.yaml' % _t)

FILEMAP = {}
for _t in ALL_TAGS:
    FILEMAP[_t] = ('%s.log' % _t, '%s.gen' % _t)
for _t in ('b0_verdict', 'b1_verdict', 'b2_verdict', 'main'):
    FILEMAP[_t] = ('%s.log' % _t, None)
FILEMAP['mesh'] = (None, 'mesh.gen')

NUM = r'([-+0-9.eEdD]+)'
G32 = re.compile(
    r'Gate 3\.2:\s*step\s+(\d+)\s+pmin\s+' + NUM +
    r'\s+pmax\s+' + NUM, re.I)
G45 = re.compile(
    r'Gate 4\.5:\s*step\s+(\d+)\s+t\s+' + NUM +
    r'\s+E\s+' + NUM + r'\s+W\s+' + NUM, re.I)
ENERGY = re.compile(
    r'=======\s*Total:\s*elastic energy:\s*' + NUM +
    r'\s+work:\s*' + NUM + r'\s+cohesive:\s*' + NUM +
    r'\s+surface:\s*' + NUM + r'\s+total:\s*' + NUM, re.I)

FATAL = ('PETSC ERROR', 'Caught signal', 'SEGV', 'fortran runtime error',
         'Fatal Error:', 'NaN detected')


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
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


def linfit(xs, ys):
    xb = sum(xs) / len(xs)
    yb = sum(ys) / len(ys)
    den = sum((x - xb) ** 2 for x in xs)
    if den <= 0.0:
        return 0.0, yb
    m = sum((x - xb) * (y - yb) for x, y in zip(xs, ys)) / den
    return m, yb - m * xb


def code_lines(text):
    return [l.split('#', 1)[0].strip() for l in text.splitlines()
            if l.split('#', 1)[0].strip()]


def code_diff(a, b):
    la, lb = code_lines(a), code_lines(b)
    if len(la) != len(lb):
        return None
    return [(x, y) for x, y in zip(la, lb) if x != y]


def fatal_in(log):
    return [s for s in FATAL if s in log]


# ---------------------------------------------------------------------------
# THREE independent implementations of alpha_eq (N3)
# ---------------------------------------------------------------------------
def alpha_ref(gc):
    """Closed form derived from C3 = a'(alpha)*Psi + C1*w'(alpha) = 0."""
    return TWO_PSI / (gc / ELL + TWO_PSI)


def alpha_lit(gc):
    """Algebraic rearrangement: alpha = 1/(1 + Gc/(ell*2Psi))."""
    return 1.0 / (1.0 + gc / (ELL * TWO_PSI))


def alpha_root(gc):
    """Numerical root of the AT2 residual itself, independent of the
    closed forms.  C1 uses the cw-1/4 contract from the Fortran source."""
    C1 = gc * 0.25 / CW / ELL
    psi = TWO_PSI / 2.0

    def C3(a):
        return (-2.0 * (1.0 - a)) * psi + C1 * (2.0 * a)

    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if C3(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def surface_ratio_pred(gc):
    """surface ~ C1*<alpha^2> with C1 ~ Gc and alpha ~ 1/Gc => 1/Gc."""
    return (gc / GC0) * (alpha_ref(gc) / alpha_ref(GC0)) ** 2


# ---------------------------------------------------------------------------
# Artefact readers
# ---------------------------------------------------------------------------
def exo_fields(path):
    c = CDF2(path)
    names = var_names(c)
    pbulk, pmax, dmax, dmean, dmin = [], [], [], [], []
    finite = True
    for i in range(c.numrecs):
        p = list(c.read_record('vals_nod_var1', i))
        d = list(c.read_record('vals_nod_var2', i))
        pbulk.append(quantile(p, 0.50))
        pmax.append(max(p))
        dmax.append(max(d))
        dmin.append(min(d))
        dmean.append(sum(d) / len(d))
        finite = finite and all(math.isfinite(x) for x in p + d)
    return {'numrecs': c.numrecs, 'names': names,
            'pbulk': pbulk, 'pmax': pmax,
            'dmax': dmax, 'dmean': dmean, 'dmin': dmin,
            'finite': finite}


def mesh_meta(path):
    c = CDF2(path)
    d = dict(zip(c.dims, c.sizes))
    return {'num_dim': d.get('num_dim'), 'num_nodes': d.get('num_nodes'),
            'sets': list(c.read_static('eb_prop1'))}


def energy_rows(log):
    return [tuple(num(m.group(i)) for i in range(1, 6))
            for m in ENERGY.finditer(log)]


def g32_rows(log):
    return [(int(m.group(1)), num(m.group(2)), num(m.group(3)))
            for m in G32.finditer(log)]


def g45_rows(log):
    return [(int(m.group(1)), num(m.group(2)), num(m.group(3)), num(m.group(4)))
            for m in G45.finditer(log)]


def beta_of(d, tag):
    p = d['exo'][tag]['pbulk'][:N_REC]
    return linfit([i * DT for i in range(N_REC)], [x - P0 for x in p])[0]


def dmax_last(d, tag):
    return d['exo'][tag]['dmax'][N_REC - 1]


def surface_last(d, tag):
    er = energy_rows(d['log'][tag])
    return er[N_REC - 1][3] if len(er) >= N_REC else -1.0


def w_last(d, tag):
    r = g45_rows(d['log'][tag])
    return r[N_REC - 1][3] if len(r) >= N_REC else -1.0


# ---------------------------------------------------------------------------
# Verdicts
# ---------------------------------------------------------------------------
def check_n1(d):
    """Run / EXO / log integrity, including EXO-vs-solver cross-check."""
    issues = []
    for tag in ALL_TAGS:
        e, log = d['exo'][tag], d['log'][tag]
        if e['numrecs'] != N_REC:
            issues.append('%s EXO has %d records (want %d)'
                          % (tag, e['numrecs'], N_REC))
        if not e['finite']:
            issues.append('%s has non-finite pressure/damage' % tag)
        if not all(len(e[k]) >= N_REC for k in
                   ('pbulk', 'pmax', 'dmax', 'dmean', 'dmin')):
            issues.append('%s EXO series truncated' % tag)
            continue
        if any(x < -1e-9 or x > 1.0 + 1e-9 for x in
               e['dmax'][:N_REC] + e['dmin'][:N_REC]):
            issues.append('%s damage left [0,1]' % tag)
        if len(g32_rows(log)) < N_REC or len(g45_rows(log)) < N_REC:
            issues.append('%s missing Gate 3.2 / Gate 4.5 rows' % tag)
        if len(energy_rows(log)) < N_REC:
            issues.append('%s has only %d energy rows'
                          % (tag, len(energy_rows(log))))
        f = fatal_in(log)
        if f:
            issues.append('%s fatal signature %s' % (tag, ','.join(f)))
        rows = g32_rows(log)
        for i in range(min(N_REC, len(rows))):
            solver_max = rows[i][2]
            exo_max = e['pmax'][i]
            if abs(exo_max - solver_max) > TOL_CROSSCHECK * max(1.0, abs(solver_max)):
                issues.append('%s rec %d EXO pmax != Gate 3.2 (stale mesh?)'
                              % (tag, i))
                break
    m = d.get('mesh', {})
    if m.get('num_dim') != 2 or m.get('num_nodes') != 356:
        issues.append('mesh is not the 356-node tri3 B0 mesh: %r' % m)
    return pass_fail(not issues, '; '.join(issues) if issues else
                     '%d runs x %d clean records, %d Gate 3.2/4.5/energy '
                     'rows each, EXO pmax == solver, 356-node tri3, no fatal'
                     % (len(ALL_TAGS), N_REC, N_REC))


def check_n2(d):
    """Gc must remain invisible on p(t): beta = f/S for the WHOLE ladder."""
    issues = []
    b0 = beta_of(d, 'b0_inject')
    d['metrics']['beta_b0_inject'] = b0
    if b0 <= 0 or abs(b0 / BETA0 - 1.0) > TOL_BETA_B3:
        issues.append('B0 reference beta=%.6g is not f/S=%.6g'
                      % (b0, BETA0))
    for tag in ALL_TAGS[1:]:
        b = beta_of(d, tag)
        dev = abs(b / b0 - 1.0) if b0 else 9.9
        d['metrics']['beta_' + tag] = b
        d['metrics']['beta_dev_' + tag] = dev
        tol = TOL_BETA_STRICT if tag in STRICT_TAGS else TOL_BETA_B3
        if dev > tol:
            issues.append('%s moved beta by %.4g (>%.4g): the initiating '
                          'Gc ladder leaked into the storage-dominated flow'
                          % (tag, dev, tol))
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'beta=%.6g Pa/s; low-damage rungs bitwise invariant '
                     '(<1e-6), max deviation %.3g across 4 decades of Gc'
                     % (b0, max(d['metrics'].get('beta_dev_' + t, 0.0)
                                for t in ALL_TAGS[1:])))


def check_n3(d):
    """AT2 analytic anchor: three implementations + real low-damage law."""
    issues = []
    s = d['src']
    if abs(TWO_PSI / TWO_PSI_FROZEN - 1.0) > 1e-6:
        issues.append('frozen 2Psi inversion drifted: %.8g vs %.8g'
                      % (TWO_PSI, TWO_PSI_FROZEN))
    for gc in (GC0, 13500.0, 2700.0, 270.0, 27.0, 2.7):
        a, b, c = alpha_ref(gc), alpha_lit(gc), alpha_root(gc)
        for label, v in (('literal', b), ('root', c)):
            if abs(v - a) > TOL_ANCHOR * max(1e-12, a):
                issues.append('alpha_eq %s != reference at Gc=%.4g '
                              '(%.12g vs %.12g)' % (label, gc, v, a))
    # Static AT2 / assembly contract (a prime source error would break the anchor).
    for key, needle, mincount in (
        ('at2', 'Constructor%cw                = 0.5_Kr', 1),
        ('at2', 'aAT2 = (1.0_kr - alpha)**2', 1),
        ('at2', 'wAT2 = alpha**2', 1),
        # C1 must be present in BOTH the operator and the Jacobian assembly;
        # a one-sided prefactor drift would otherwise go undetected.
        ('defasm', 'C1 = matpropSet%fractureToughness / ATModel%cw * 0.25_Kr / matpropSet%internalLength', 2),
        ('defasm', 'C3 = ATModel%Da(damageGauss) * EEDGaussPlus + C1 * ATModel%Dw(damageGauss)', 1),
    ):
        got = s.get(key, '').count(needle)
        if got < mincount:
            issues.append('AT2 / C3 assembly contract %s: found %d, want %d (%s)'
                          % (key, got, mincount, needle[:48]))
    # Independent REAL calibration from CI #127 (B1 gc arm).
    pred = surface_ratio_pred(13500.0)
    d['metrics']['surface_ratio_gc_pred'] = pred
    d['metrics']['surface_ratio_gc_real'] = SURF_RATIO_GC_REAL
    if abs(pred / SURF_RATIO_GC_REAL - 1.0) > TOL_SURF_CALIB:
        issues.append('analytic surface ratio %.6f disagrees with CI #127 '
                      'real %.6f' % (pred, SURF_RATIO_GC_REAL))
    if abs(BETA_DEV_GC_REAL) > TOL_BETA_STRICT:
        issues.append('CI #127 beta_dev_gc is not zero: %.6g'
                      % BETA_DEV_GC_REAL)
    # Low-damage 1/Gc law measured from the real B3 EXO.
    a0 = dmax_last(d, 'b0_inject')
    if a0 <= 0:
        issues.append('B0 damage_max is non-positive')
    else:
        for tag in STRICT_TAGS:
            r = dmax_last(d, tag) / a0
            rp = alpha_ref(GC_OF[tag]) / alpha_ref(GC0)
            d['metrics']['dmg_ratio_' + tag] = r
            d['metrics']['dmg_ratio_pred_' + tag] = rp
            if rp <= 0 or r <= 0 or abs(math.log10(r / rp)) > TOL_RATIO_DECADES:
                issues.append('%s dmax ratio %.4g vs analytic %.4g '
                              '(>%.2f decades): AT2 1/Gc branch broken'
                              % (tag, r, rp, TOL_RATIO_DECADES))
    xs = [math.log(GC_OF[t]) for t in LOW_TAGS]
    ys = [math.log(max(dmax_last(d, t), 1e-300)) for t in LOW_TAGS]
    slope, _ = linfit(xs, ys)
    d['metrics']['dmg_power_slope'] = slope
    if not (POWER_BAND_DMG[0] < slope < POWER_BAND_DMG[1]):
        issues.append('low-damage dmax ~ Gc^%.4f outside %s'
                      % (slope, POWER_BAND_DMG))
    return pass_fail(not issues, '; '.join(issues) if issues else
                     '2Psi=%.6e Pa, Gc*(a=1/2)=%.6f N/m; 3 independent '
                     'alpha_eq implementations agree <1e-10; CI #127 real '
                     'surface x%.6f vs analytic x%.6f; dmg ~ Gc^%.4f'
                     % (TWO_PSI, GC_STAR_HALF, SURF_RATIO_GC_REAL,
                        pred, slope))


def check_n4(d):
    """Initiation access: at least one rung crosses stiffness half-life."""
    issues = []
    al = [dmax_last(d, t) for t, _ in LADDER]
    for t, a in zip([x[0] for x in LADDER], al):
        d['metrics']['alpha_' + t] = a
    if any(a < -1e-9 or a > 1.0 + 1e-9 for a in al):
        issues.append('damage outside [0,1]: %s'
                      % ','.join('%.4g' % a for a in al))
    for i in range(1, len(al)):
        if al[i] < al[i - 1] * (1.0 - 1e-9) - 1e-12:
            issues.append('damage NOT monotone in 1/Gc at %s (%.4g < %.4g)'
                          % (LADDER[i][0], al[i], al[i - 1]))
    first = next((i for i, a in enumerate(al) if a >= ALPHA_STIFF_HALF), None)
    if first is None:
        issues.append('NO rung reached alpha_stiff_half=%.6f (max %.4g): '
                      'the platform did NOT reach initiation'
                      % (ALPHA_STIFF_HALF, max(al)))
        return pass_fail(False, '; '.join(issues))
    if first == 0:
        issues.append('B0 itself is already initiated: the ladder has no '
                      'pre-initiation anchor')
    g_lo = GC_OF[LADDER[first][0]]
    g_hi = GC_OF[LADDER[first - 1][0]] if first > 0 else GC0 * 10.0
    d['metrics']['initiation_first_tag'] = LADDER[first][0]
    d['metrics']['gc_bracket_lo'] = g_lo
    d['metrics']['gc_bracket_hi'] = g_hi
    d['metrics']['alpha_stiff_half'] = ALPHA_STIFF_HALF
    d['metrics']['gc_star_stiff'] = GC_STAR_STIFF
    # The analytic Gc* must sit within one decade of the measured bracket.
    if GC_STAR_STIFF > g_hi * (10.0 ** TOL_GCSTAR_DECADES) or \
       GC_STAR_STIFF < g_lo / (10.0 ** TOL_GCSTAR_DECADES):
        issues.append('analytic Gc*=%.4g N/m is outside the measured '
                      'bracket [%.4g, %.4g] (+/-1 decade)'
                      % (GC_STAR_STIFF, g_lo, g_hi))
    growth = al[first] / max(al[0], 1e-300)
    d['metrics']['initiation_growth'] = growth
    if growth < 1.0e3:
        issues.append('initiation grew damage only x%.4g (<3 decades)'
                      % growth)
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'initiation at %s (Gc=%.4g, alpha=%.4g >= %.6f); '
                     'measured bracket [%.4g, %.4g] N/m brackets analytic '
                     'Gc*=%.4f N/m; damage grew x%.4g over 4 decades of Gc'
                     % (LADDER[first][0], g_lo, al[first], ALPHA_STIFF_HALF,
                        g_lo, g_hi, GC_STAR_STIFF, growth))


def check_n5(d):
    """Localisation (a crack, not a uniformly softened blob) + solver hygiene."""
    issues = []
    ratios = {}
    for tag in ALL_TAGS:
        e = d['exo'][tag]
        dm = e['dmean'][N_REC - 1]
        ratio = e['dmax'][N_REC - 1] / dm if dm > 0 else 0.0
        ratios[tag] = ratio
        d['metrics']['local_' + tag] = ratio
        if ratio <= LOCAL_MIN:
            issues.append('%s damage is spatially uniform (dmax/dmean=%.4g '
                          '<= %.2f): uniform softening is NOT a crack'
                          % (tag, ratio, LOCAL_MIN))
        seq = e['dmax'][:N_REC]
        for i in range(1, len(seq)):
            if seq[i] < seq[i - 1] * (1.0 - 1e-6) - 1e-12:
                issues.append('%s damage decreased at rec %d: irreversibility '
                              'bound broken' % (tag, i))
                break
    n5 = sum(len(re.findall(r'damage field SNESSolve failed', d['log'][t]))
             for t in ALL_TAGS)
    d['metrics']['damage_sneg5_msgs'] = n5
    d['metrics']['evidence'] = ('relative-only (damage SNES -5 present)'
                                if n5 else 'current-field integrals')
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'damage localised on every rung (min dmax/dmean=%.2f, '
                     'REAL B0=7.00), bounded in [0,1], irreversible; '
                     '%d damage-SNES -5 messages REPORTED, not hidden'
                     % (min(ratios.values()), n5))


def check_n6(d):
    """Energy consistency: surface ~ 1/Gc low branch, W invariant."""
    issues = []
    surf = {}
    for tag in ALL_TAGS:
        er = energy_rows(d['log'][tag])
        if len(er) < N_REC:
            issues.append('%s has only %d energy rows' % (tag, len(er)))
            continue
        if any(not all(math.isfinite(x) for x in r) for r in er[:N_REC]):
            issues.append('%s energy row is non-finite' % tag)
        surf[tag] = er[N_REC - 1][3]
        d['metrics']['surface_' + tag] = surf[tag]
        if surf[tag] <= 0.0:
            issues.append('%s surface energy is non-positive' % tag)
        w = w_last(d, tag)
        d['metrics']['W_' + tag] = w
        if w < 0.0:
            issues.append('%s injection work is negative' % tag)
    if len(surf) == len(ALL_TAGS):
        wref = w_last(d, 'b0_inject')
        for tag in ALL_TAGS[1:]:
            w = w_last(d, tag)
            if wref > 0 and w > 0 and abs(w / wref - 1.0) > TOL_W:
                issues.append('%s W moved by %.4g (>%.4g): flow field '
                              'changed despite N2' % (tag, abs(w / wref - 1.0), TOL_W))
        for i in range(1, len(LOW_TAGS)):
            if surf[LOW_TAGS[i]] <= surf[LOW_TAGS[i - 1]]:
                issues.append('surface energy NOT increasing as Gc falls '
                              'at %s (%.4g <= %.4g)'
                              % (LOW_TAGS[i], surf[LOW_TAGS[i]],
                                 surf[LOW_TAGS[i - 1]]))
        xs = [math.log(GC_OF[t]) for t in LOW_TAGS]
        ys = [math.log(max(surf[t], 1e-300)) for t in LOW_TAGS]
        slope, _ = linfit(xs, ys)
        d['metrics']['surface_power_slope'] = slope
        if not (POWER_BAND_SURF[0] < slope < POWER_BAND_SURF[1]):
            issues.append('low-damage surface ~ Gc^%.4f outside %s'
                          % (slope, POWER_BAND_SURF))
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'surface_ref=%.4g; W invariant to %.2g; low-damage '
                     'surface ~ Gc^%.4f (analytic -1); all energy rows finite'
                     % (d['metrics'].get('surface_b0_inject', -1), TOL_W,
                        d['metrics'].get('surface_power_slope', -1)))


def check_n7(d):
    """Static anti-regression: deck contract, source contract, CI wiring."""
    issues = []
    base = d['deck'].get('inject', '')
    if not base:
        issues.append('B0 injection deck missing')
    for tag, fname in (('b1_gc', 'b1_gc'),) + tuple((t, t) for t in NEW_TAGS):
        txt = d['deck'].get(fname, '')
        diff = code_diff(base, txt)
        if diff is None or len(diff) != 1:
            issues.append('%s deck is not a ONE-code-line delta (%s)'
                          % (tag, 'length changed' if diff is None
                             else len(diff)))
            continue
        x, y = diff[0]
        if not (x.startswith('fracturetoughness:') and
                y.startswith('fracturetoughness:')):
            issues.append('%s deck changed the wrong key: %s -> %s'
                          % (tag, x, y))
    s = d['src']
    wf = s.get('workflow', '')
    for needle in ('check_gate415_b3_initiation.py --gate /tmp/gate415b3',
                   'gb3_rc=${PIPESTATUS[0]}',
                   'if [ "$gb3_rc" -ne 0 ]'):
        if needle not in wf:
            issues.append('CI wiring missing: %s' % needle)
    for t in NEW_TAGS:
        if 'terzaghi415_b3_%s.yaml' % t not in wf:
            issues.append('CI has no run for B3 deck %s' % t)
    n_rebuild = wf.count('rm -f data/Beam1x10-tri3_out.gen')
    d['metrics']['exo_rebuild_count'] = n_rebuild
    if n_rebuild < 14:
        issues.append('only %d EXO rebuilds (B0/B1/B2=10 + B3=4 required)'
                      % n_rebuild)
    # The anchor formula must stay in the deck comments (falsifiable docs).
    if 'alpha_eq(Gc) = 2*Psi / (Gc/ell + 2*Psi)' not in d['deck'].get('gc2p7', ''):
        issues.append('B3 deck lost the analytic-anchor comment contract')
    return pass_fail(not issues, '; '.join(issues) if issues else
                     '5 one-line Gc decks; AT2 C1/C2/C3/cw source contract '
                     'in place; %d EXO rebuilds; gb3_rc deferred hard fail'
                     % n_rebuild)


def check_n8(d):
    """State discipline + zero production perturbation."""
    issues = []
    want = {'b0_verdict': 'GATE415_RESULT: PASS',
            'b1_verdict': 'GATE415_B1_RESULT: PASS',
            'b2_verdict': 'GATE415_B2_RESULT: PASS'}
    for tag, needle in want.items():
        if needle not in d['log'].get(tag, ''):
            issues.append('%s lost "%s"' % (tag, needle))
    main = d['log'].get('main', '')
    if 'terzaghi415_b3_' in main:
        issues.append('a B3 deck leaked into the 51-step production log')
    # B3 must NOT claim any propagation-regime deliverable.
    doc = read(os.path.join(_HERE, 'check_gate415_b3_initiation.py'))
    for needle in ('STATE DISCIPLINE: B3 establishes an initiation threshold, '
                   'NOT a propagating',
                   'NO claim is made about'):
        if needle not in doc:
            issues.append('state-discipline declaration missing: %s'
                          % needle[:40])
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'B0/B1/B2 verdicts PASS and production log untouched; '
                     'B3 claims initiation access ONLY, no propagation '
                     'scaling / fracture length / rate-optimization table')


CHECKS = (('N1', 'run / EXO / log / energy-row integrity', check_n1),
          ('N2', 'Flow invariance: beta=f/S over 4 decades of Gc', check_n2),
          ('N3', 'AT2 analytic anchor (3 implementations + CI #127)', check_n3),
          ('N4', 'initiation access + measured Gc* bracket', check_n4),
          ('N5', 'localisation, bounds, irreversibility, -5 reporting', check_n5),
          ('N6', 'surface ~ 1/Gc + injection-work invariance', check_n6),
          ('N7', 'static anti-regression + deck/CI contracts', check_n7),
          ('N8', 'state discipline + zero production perturbation', check_n8))


def run(d):
    d['metrics'] = {}
    rows = []
    for tag, label, fn in CHECKS:
        verdict, detail = fn(d)
        rows.append((tag, label, verdict, detail))
    fails = [r[0] for r in rows if r[2] == 'FAIL']
    return rows, fails


METRICS_ROW = {}


def report(rows, fails, header):
    print(header)
    for tag, label, verdict, detail in rows:
        print('%s %s %-52s | %s' % (verdict, tag, label, detail))
    print('GATE415_B3_RESULT: %s%s' % ('PASS' if not fails else 'FAIL',
                                       '' if not fails else ' (' + ','.join(fails) + ')'))
    for tag, label, verdict, detail in rows:
        if verdict == 'FAIL':
            print('  %s detail: %s' % (tag, detail))
    m = METRICS_ROW.get('metrics', {})
    for key in ('alpha_b0_inject', 'alpha_gc_half', 'alpha_gc2700',
                'alpha_gc270', 'alpha_gc27', 'alpha_gc2p7',
                'alpha_stiff_half', 'gc_star_stiff', 'gc_bracket_lo',
                'gc_bracket_hi', 'initiation_growth', 'dmg_power_slope',
                'surface_power_slope', 'surface_ratio_gc_pred',
                'surface_ratio_gc_real', 'beta_dev_gc2p7',
                'local_b0_inject', 'local_gc2p7', 'damage_sneg5_msgs',
                'exo_rebuild_count'):
        v = m.get(key, -1)
        if isinstance(v, str):
            print('gb3_%-28s = %s' % (key, v))
        else:
            print('gb3_%-28s = %.6e' % (key, v))


# ---------------------------------------------------------------------------
# Synthetic offline case (real constants + analytic solver artefacts)
# ---------------------------------------------------------------------------
def make_log(tag, pmax, work, surface, damage, n=N_REC):
    log = '=== Gate 4.15 B3 synthetic run %s ===\n' % tag + 'x' * 1200 + '\n'
    for i in range(1, n + 1):
        log += 'Gate 3.1: step %d enter transient block\n' % i
        log += 'Gate 3.2: step %d pmin 0.0 pmax %.10E\n' % (i, pmax[i - 1])
        log += ('Gate 4.12: step %d t %.4g dpmax %.4g dpmin 0.0\n'
                % (i, (i - 1) * DT, pmax[i - 1] - pmax[i - 2] if i > 1 else 0.0))
        log += ('Gate 4.5: step %d t %.10E E %.10E W %.10E\n'
                % (i, (i - 1) * DT, 1.0e-3, work[i - 1]))
        log += ('======= Total: elastic energy:  %.5E work:  %.5E '
                'cohesive:  0.00000E+00 surface:  %.5E total:  %.5E\n'
                % (max(0.5 * work[i - 1], 1.0e-6), work[i - 1] or 1.0e-6,
                   surface[i - 1], 1.0))
        if i > 1:
            log += '    Linear solve converged due to CONVERGED_RTOL iterations 1\n'
    return log


B0_VERDICT = """=== Gate 4.15 / Stage B0: rate/fracability similitude ===
PASS V1 run / EXO / mesh data integrity | 11 + 11 clean records
PASS V4 injection slope f/S + source-unit calibration | bulk slope=100000 Pa/s
GATE415_RESULT: PASS
"""
B1_VERDICT = """=== Gate 4.15 / Stage B1 ===
PASS A4 Gc/k/mu invisible on the pressure curve
PASS A5 Gc visible in phase-field surface energy
GATE415_B1_RESULT: PASS
"""
B2_VERDICT = """=== Gate 4.15 / Stage B2 ===
PASS L4 pre-initiation anisotropy null proposition
GATE415_B2_RESULT: PASS
"""


def synth_sources():
    return {k: read(v) for k, v in SRC.items()}


def synth_decks():
    return {k: read(v) for k, v in DECKS.items()}


def synth():
    """Analytic positive case anchored on REAL CI #125/#127 readings."""
    d = {'log': {}, 'exo': {}, 'metrics': {}}
    names = ['Pressure', 'Damage', 'Displacement_X', 'Displacement_Y']
    pbulk = [P0 + BETA0 * i * DT for i in range(N_REC)]
    pmax = [P0 + 2.9184 * BETA0 * i * DT for i in range(N_REC)]
    work = [W_CALIB * F0 * AREA * p for p in pbulk]
    # Offline damage ladder.  The high-damage rungs keep the analytic branch
    # (the solver may snap higher; N4 only needs the threshold crossing).
    for tag, gc in LADDER:
        a_last = alpha_ref(gc)
        dmax = [a_last * ((i + 1.0) / N_REC) ** 2 for i in range(N_REC)]
        dmean = [x / 7.0 for x in dmax]           # REAL B0 dmax/dmean = 7.00
        dmin = [x / 26.0 for x in dmax]           # REAL B0 dmax/dmin = 26.28
        s_last = SURF0 * surface_ratio_pred(gc)
        surface = [s_last * ((i + 1.0) / N_REC) ** 2 for i in range(N_REC)]
        d['exo'][tag] = {'numrecs': N_REC, 'names': names, 'finite': True,
                         'pbulk': list(pbulk), 'pmax': list(pmax),
                         'dmax': dmax, 'dmean': dmean, 'dmin': dmin}
        d['log'][tag] = make_log(tag, pmax, work, surface, dmax)
    d['log']['b0_verdict'] = B0_VERDICT
    d['log']['b1_verdict'] = B1_VERDICT
    d['log']['b2_verdict'] = B2_VERDICT
    d['log']['main'] = make_log('main', [P0 * (0.0906 ** (j / 50.0))
                                        for j in range(51)],
                                [0.0] * 51, [1.0e-6] * 51,
                                [1.0e-6] * 51, n=51)
    d['mesh'] = {'num_dim': 2, 'num_nodes': 356, 'sets': [1, 10, 20]}
    d['src'] = synth_sources()
    d['deck'] = synth_decks()
    return d


def selftest():
    base = synth()
    cases = []
    rows, fails = run(base)
    METRICS_ROW['metrics'] = base['metrics']
    cases.append(('POS full B3 verdict', rows, set(), None))
    print('POS  %-58s %s' % ('full B3 verdict',
                            'PASS' if not fails else 'FAIL ' + ','.join(fails)))

    def neg(name, mutate, expect, label='NEG'):
        d = copy.deepcopy(base)
        mutate(d)
        r, f = run(d)
        ok = set(f) == set(expect)
        cases.append((name, r, set(expect), f))
        print('%s  %-58s %s (failed %s, expected %s)'
              % (label, name, 'OK' if ok else 'MISMATCH',
                 sorted(f), sorted(expect)))
        return ok

    def inv(name, mutate):
        return neg(name, mutate, set(), 'INV')

    allok = not fails
    # --- N1 integrity -----------------------------------------------------
    allok &= neg('gc2700 EXO truncated to 10 records',
                 lambda d: d['exo']['gc2700'].update(numrecs=10), {'N1'})
    allok &= neg('fatal signature in the gc2p7 log',
                 lambda d: d['log'].update(
                     gc2p7=d['log']['gc2p7'] + 'PETSC ERROR: boom\n'), {'N1'})
    allok &= neg('EXO pmax disagrees with Gate 3.2 (stale mesh)',
                 lambda d: d['exo']['gc27']['pmax'].__setitem__(
                     4, d['exo']['gc27']['pmax'][4] * 0.5), {'N1'})
    allok &= neg('damage exceeds the physical upper bound',
                 lambda d: d['exo']['gc2p7']['dmax'].__setitem__(
                     N_REC - 1, 1.7), {'N1', 'N4'})
    allok &= neg('energy rows missing from gc270 log',
                 lambda d: d['log'].update(
                     gc270=d['log']['gc270'].replace('======= Total:', '==== NOPE:')),
                 {'N1', 'N6'})
    # --- N2 flow invariance ------------------------------------------------
    allok &= neg('Gc ladder breaks the pressure ramp',
                 lambda d: d['exo']['gc2p7'].update(
                     pbulk=[P0 + 1.05 * BETA0 * i * DT for i in range(N_REC)]),
                 {'N2'})
    allok &= neg('low-damage rung perturbs beta by 1e-4 (not bitwise)',
                 lambda d: d['exo']['gc2700'].update(
                     pbulk=[P0 + 1.0001 * BETA0 * i * DT for i in range(N_REC)]),
                 {'N2'})
    # --- N3 analytic anchor ------------------------------------------------
    allok &= neg('AT2 cw constant changed in the source',
                 lambda d: d['src'].update(
                     at2=d['src']['at2'].replace(
                         'Constructor%cw                = 0.5_Kr',
                         'Constructor%cw                = 2.0_Kr')), {'N3'})
    allok &= neg('C3 assembly no longer uses Da*Psi + C1*Dw',
                 lambda d: d['src'].update(
                     defasm=d['src']['defasm'].replace(
                         'C3 = ATModel%Da(damageGauss) * EEDGaussPlus + C1 * ATModel%Dw(damageGauss)',
                         'C3 = 0.0_Kr')), {'N3'})
    allok &= neg('C1 prefactor drifted from 0.25/cw',
                 lambda d: d['src'].update(
                     defasm=d['src']['defasm'].replace(
                         'C1 = matpropSet%fractureToughness / ATModel%cw * 0.25_Kr / matpropSet%internalLength',
                         'C1 = matpropSet%fractureToughness / ATModel%cw * 1.25_Kr / matpropSet%internalLength',
                         1)), {'N3'})
    allok &= neg('damage frozen across the whole Gc ladder (C1 unwired)',
                 lambda d: [d['exo'][t].update(
                     dmax=[ALPHA0_B0] * N_REC,
                     dmean=[ALPHA0_B0 / 7.0] * N_REC,
                     dmin=[ALPHA0_B0 / 26.0] * N_REC) for t in ALL_TAGS[1:]],
                 {'N3', 'N4'})
    allok &= neg('low-damage 1/Gc power law inverted',
                 lambda d: [d['exo'][t].update(
                     dmax=list(reversed(d['exo'][t]['dmax']))) for t in STRICT_TAGS],
                 {'N3', 'N4', 'N5'})
    # --- N4 initiation -----------------------------------------------------
    allok &= neg('no rung reaches initiation (Gc ladder ineffective)',
                 lambda d: [d['exo'][t].update(
                     dmax=[min(x, 1.0e-6) for x in d['exo'][t]['dmax']],
                     dmean=[min(x, 1.0e-7) for x in d['exo'][t]['dmean']],
                     dmin=[min(x, 1.0e-8) for x in d['exo'][t]['dmin']])
                     for t in ('gc27', 'gc2p7')], {'N4'})
    allok &= neg('initiation growth is only 2 decades',
                 lambda d: d['exo']['gc2p7'].update(
                     dmax=[min(x, 0.01) for x in d['exo']['gc2p7']['dmax']],
                     dmean=[min(x, 0.01 / 7.0) for x in d['exo']['gc2p7']['dmean']],
                     dmin=[min(x, 0.01 / 26.0) for x in d['exo']['gc2p7']['dmin']]),
                 {'N4'})
    allok &= neg('damage decreases as Gc falls (monotonicity broken)',
                 lambda d: d['exo']['gc270'].update(
                     dmax=[x * 0.01 for x in d['exo']['gc270']['dmax']]),
                 {'N3', 'N4', 'N5'})
    # --- N5 localisation ---------------------------------------------------
    allok &= neg('uniform softening instead of a localised crack',
                 lambda d: d['exo']['gc2p7'].update(
                     dmean=d['exo']['gc2p7']['dmax'][:]), {'N5'})
    allok &= neg('irreversibility broken (damage shrinks in time)',
                 lambda d: d['exo']['gc27']['dmax'].__setitem__(N_REC - 1, 1.0e-8),
                 {'N4', 'N5'})
    # --- N6 energy ---------------------------------------------------------
    allok &= neg('surface energy frozen while Gc falls 4 decades',
                 lambda d: [d['log'].__setitem__(
                     t, re.sub(r'surface:  [-+0-9.E]+',
                               'surface:  6.49012E-04', d['log'][t]))
                     for t in NEW_TAGS], {'N6'})
    allok &= neg('injection work changes despite invariant pressure',
                 lambda d: d['log'].__setitem__(
                     'gc2p7', re.sub(
                         r'(Gate 4\.5:\s*step\s+11\s+.*?\s+W\s+)[-+0-9.eEdD]+',
                         lambda m: m.group(1) + '9.99999E+09',
                         d['log']['gc2p7'], count=1)), {'N6'})
    # --- N7 static contracts ----------------------------------------------
    allok &= neg('one-line deck contract broken (second code delta)',
                 lambda d: d['deck'].update(
                     gc27=d['deck']['gc27'].replace('p0: 1.0e5', 'p0: 2.0e5')),
                 {'N7'})
    allok &= neg('B3 deck changed the wrong key',
                 lambda d: d['deck'].update(
                     gc270=d['deck']['gc270'].replace(
                         'fracturetoughness: 270.0', 'fracturetoughness: 27000.0')),
                 {'N7'})
    allok &= neg('all EXO rebuilds removed from the workflow',
                 lambda d: d['src'].update(workflow=d['src']['workflow'].replace(
                     'rm -f data/Beam1x10-tri3_out.gen', 'true # no rebuild')),
                 {'N7'})
    allok &= neg('deferred gb3_rc hard fail removed',
                 lambda d: d['src'].update(workflow=d['src']['workflow'].replace(
                     'if [ "$gb3_rc" -ne 0 ]', 'if false')), {'N7'})
    # --- N8 state discipline / isolation ----------------------------------
    allok &= neg('B0 verdict silently became FAIL',
                 lambda d: d['log'].update(
                     b0_verdict=B0_VERDICT.replace('GATE415_RESULT: PASS',
                                                   'GATE415_RESULT: FAIL')), {'N8'})
    allok &= neg('B3 deck name leaked into the production log',
                 lambda d: d['log'].update(
                     main=d['log']['main'] +
                     '-options_file_yaml data/terzaghi415_b3_gc2p7.yaml\n'), {'N8'})
    # --- invariances -------------------------------------------------------
    allok &= inv('damage SNES residual plateau (-5) appears, physics unchanged',
                 lambda d: d['log'].update(
                     gc2p7=d['log']['gc2p7'] +
                     '[ERROR]: damage field SNESSolve failed with SNESConvergedReason -5.\n'))
    allok &= inv('deck comments rewritten (code lines untouched)',
                 lambda d: d['deck'].update(
                     gc2p7=d['deck']['gc2p7'].replace(
                         '# STATE DISCIPLINE:', '# REWORDED DISCIPLINE:')))
    allok &= inv('initiation rung snaps to alpha=1 (still a valid crossing)',
                 lambda d: d['exo']['gc2p7'].update(
                     dmax=[min(1.0, x * 1.7) for x in d['exo']['gc2p7']['dmax']],
                     dmean=[min(1.0 / 7.0, x * 1.7 / 7.0)
                            for x in d['exo']['gc2p7']['dmean']],
                     dmin=[min(1.0 / 26.0, x * 1.7 / 26.0)
                           for x in d['exo']['gc2p7']['dmin']]))
    print('SELFTEST_RESULT: %s (%d positive + %d mutations, incl. 3 invariances)'
          % ('PASS' if allok else 'FAIL', 1, len(cases) - 1))
    return 0 if allok else 1


def predict():
    print('=== Gate 4.15 / Stage B3 OFFLINE initiation ladder ===')
    print('2*Psi = %.6e Pa   Gc0/ell = %.6e Pa   alpha_stiff_half = %.10f'
          % (TWO_PSI, GC0 / ELL, ALPHA_STIFF_HALF))
    print('Gc*(alpha=0.5) = %.6f N/m   Gc*(a=1/2) = %.6f N/m'
          % (GC_STAR_HALF, GC_STAR_STIFF))
    print('%-12s %-12s %-14s %-14s %-12s' %
          ('tag', 'Gc [N/m]', 'alpha_eq', 'dmax/a0', 'surface_pred'))
    for tag, gc in LADDER:
        a = alpha_ref(gc)
        print('%-12s %-12.4g %-14.6e %-14.6g %-12.6e' %
              (tag, gc, a, a / ALPHA0_B0, SURF0 * surface_ratio_pred(gc)))
    print('CI #127 REAL surface ratio (Gc x0.5) = %.6f; analytic = %.6f'
          % (SURF_RATIO_GC_REAL, surface_ratio_pred(13500.0)))
    return 0


def main(argv):
    if argv and argv[0] == '--predict':
        return predict()
    if argv and argv[0] == '--selftest':
        return selftest()
    if argv and argv[0] == '--gate':
        if len(argv) != 2:
            print('GATE415_B3_RESULT: FAIL (usage: --gate <artifact_dir>)')
            return 1
        root = argv[1]
        try:
            d = {'log': {}, 'exo': {}, 'metrics': {}}
            for tag in ALL_TAGS + ('b0_verdict', 'b1_verdict',
                                   'b2_verdict', 'main'):
                lp = FILEMAP[tag][0]
                d['log'][tag] = read(os.path.join(root, lp))
            for tag in ALL_TAGS:
                d['exo'][tag] = exo_fields(os.path.join(root, FILEMAP[tag][1]))
            d['mesh'] = mesh_meta(os.path.join(root, FILEMAP['mesh'][1]))
            d['src'] = synth_sources()
            d['deck'] = synth_decks()
        except Exception as exc:
            print('=== Gate 4.15 / Stage B3: fracture-initiation access ===')
            print('FAIL N1 run / EXO / log / energy-row integrity | '
                  'cannot parse artefacts: %r' % exc)
            print('GATE415_B3_RESULT: FAIL (N1 data error: %r)' % exc)
            return 1
        rows, fails = run(d)
        METRICS_ROW['metrics'] = d['metrics']
        report(rows, fails,
               '=== Gate 4.15 / Stage B3: fracture-initiation access ===')
        return 1 if fails else 0
    print('usage: check_gate415_b3_initiation.py '
          '(--gate <dir> | --selftest | --predict)')
    return 2


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
