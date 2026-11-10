#!/usr/bin/env python3
"""Gate 4.15 / Stage B0 -- hard verdict for the rate/fracability similitude
calibration and source-unit calibration smoke run.

Why this gate exists
--------------------
Every earlier Gate ran the Terzaghi benchmark at k_min=5e-11 m2, S=4.6e-10
1/Pa, so c_v = 108.7 m2/s and T = c_v*t/H^2 ~ 1 at t=1 s: the whole 1 s
experiment is ONE consolidation event.  Field shale at t_f=3600 s, H_f=50 m
has T ~ 1.6e-7 (RATE_FRACABILITY_ROADMAP.md, Table A/B), i.e. the injected
fluid is stored, not diffused away.  A rate/fracability study that keeps the
production permeability therefore cannot raise any net pressure at all.

Stage B0 introduces TWO new 2D decks (no Fortran change):
  terzaghi415_scale.yaml   vs terzaghi1D.yaml        -> 4 code lines changed
      numstep 51->11, perm_min 5e-11->2.35e-18,
      perm_max 1e-9->4.70e-17, storativity 4.6e-10->1.50e-10
  terzaghi415_inject.yaml  vs terzaghi415_scale.yaml -> 1 code line changed
      cs0001.source 0.0->1.5e-5 [1/s]

With c_v = k/(mu*S) = 1.5667e-5 m2/s and T(1 s) = 1.5667e-7, the Darcy term
is negligible over 1 s except in a ~4 mm boundary layer, so the discrete
balance is exactly

    S dp/dt = f    =>    p(t) = p0 + (f/S) t,   f/S = 1.0e5 Pa/s

and the zero-source control stays frozen at p0.  This is a FALSIFIABLE anchor:
it fixes the source unit.  PetscBag registers "Volumetric source/sink rate
(1/s)" and MEF90FlowEnergy computes W = int f*p dOmega, whereas an old deck
comment claimed [Pa/s].  If f were [Pa/s] the measured slope would be f, not
f/S (a factor 1/S = 6.67e9).  V4 recovers S_eff = f/slope from the real EXO
and must give back the deck storativity.

V1 run / EXO / mesh data integrity
V2 similitude parameters frozen from the decks (T, c_v, k/S ratios, Q_total)
V3 zero-source control is storage-dominated / frozen at p0
V4 injection is linear with slope f/S and rejects the [Pa/s] interpretation
V5 Gate 4.5 injection work (control W=0, injection W>0 and monotone)
V6 solver chain completed, per-step KSP convergence, no fatal signature
V7 static anti-regression + code-line deck contracts + CI wiring
V8 production isolation (main log still consolidates, no 415 leakage)

CI #124 post-mortem (two INDEPENDENT defects, both fixed here)
--------------------------------------------------------------
(1) EXO reuse.  MEF90CtxOpenEXO appends to an existing result mesh, so the two
    11-step B0 runs inherited the 51-record production file and left 40 stale
    records after their own 11.  V1 rejected 51 != 11, and V3 was reading ALL
    51 records (it did not slice to N_REC), so it reported a 90.94% "drift"
    that was really the old production decay tail.  The workflow now deletes
    the result mesh before EACH B0 run so vDef takes the EXCRE(EXCLOB) branch;
    V1 additionally cross-checks EXO pmax against the solver's Gate 3.2 rows.

(2) Wrong observable.  The analytic anchor p = p0 + (f/S)t describes the BULK
    pressure.  In a storage-dominated run (A=k/mu=2.35e-15, a*S=1.5e-9) the
    consistent mass matrix plus the strong Dirichlet inlet produces a LOCAL
    numerical overshoot at a handful of nodes: CI #124 measured pmax=3.92e5 Pa
    at step 10 while the median was exactly p50=2.00e5 Pa.  Gate 4.5's work
    W=f*Area*p_bulk independently confirmed the 1.0e5 Pa/s bulk ramp.  Using
    pmax therefore measured the overshoot (2.918e5 Pa/s), not f/S.  V3/V4 now
    use the median pbulk as the physical observable, keep pmax only as a
    LOCALITY diagnostic (p90 must stay near pbulk), and still hard-fail if the
    bulk ramp is wrong or the source is read as [Pa/s].

Pure standard library (the CI runner has no numpy / scipy / PyYAML), no local
compile and no local solver.  --gate reads real CI artefacts; --selftest
replays one synthetic positive case (including the real overshoot) plus >=12
negative mutations and an invariance case.
"""

import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse_cdf2 import CDF2  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, '..', '..'))
DECKS = {
    'main': os.path.join(_HERE, 'terzaghi1D.yaml'),
    'scale': os.path.join(_HERE, 'terzaghi415_scale.yaml'),
    'inject': os.path.join(_HERE, 'terzaghi415_inject.yaml'),
}
SRC = {
    'flowctx': os.path.join(_ROOT, 'm_Flow', 'm_MEF90_FlowCtx.F90'),
    'flowasm': os.path.join(_ROOT, 'm_Flow', 'm_MEF90_FlowAssembly.F90'),
    'vdef': os.path.join(_ROOT, 'vDef', 'vDef.F90'),
    'workflow': os.path.join(_ROOT, '.github', 'workflows', 'build.yml'),
}

NAME_WIDTH = 33
N_REC = 11                 # 11 EXO records, dt = 0.1 s, t = 0..1 s
DT = 0.1
P0 = 1.0e5                 # initial pore pressure [Pa]
F_SOURCE = 1.5e-5          # cs0001.source [1/s]
K_DECK = 5.0e-11           # production perm_min [m2]
S_DECK = 4.6e-10           # production storativity [1/Pa]
MU = 1.0e-3                # viscosity [Pa.s]
S_NUM = 1.50e-10           # calibrated storativity [1/Pa] (frozen by V4)
H = 10.0                   # drainage path [m]
AREA = 10.0                # 10 m x 1 m, unit thickness [m2]
T_TARGET = (2.35e-18 / (MU * S_NUM)) * 1.0 / (H * H)   # 1.5667e-7
SLOPE_TARGET = F_SOURCE / S_NUM                         # 1.0e5 Pa/s
Q_TOTAL = F_SOURCE * AREA                                # 1.5e-4 m2/s
S_RATIO_TARGET = S_NUM / S_DECK                          # 0.326
TOL_SLOPE = 0.10
TOL_FREEZE = 0.02
TOL_LINEAR = 0.06
TOL_PEAK_LOCAL = 0.05        # p90 may not depart from the bulk ramp by >5%
TOL_CROSSCHECK = 1.0e-6      # EXO pmax vs solver Gate 3.2 pmax (relative)
PEAK_RATIO_BAND = (1.0, 6.0)  # local overshoot is a fixed spatial mode
FATAL = ('PETSC ERROR', 'Caught signal', 'SEGV', 'fortran runtime error',
         'Fatal Error:', 'NaN detected')

G45 = re.compile(
    r'Gate 4\.5:\s*step\s+(\d+)\s+t\s+([-+0-9.eEdD]+)\s+'
    r'E\s+([-+0-9.eEdD]+)\s+W\s+([-+0-9.eEdD]+)', re.I)
G32 = re.compile(
    r'Gate 3\.2:\s*step\s+(\d+)\s+pmin\s+([-+0-9.eEdD]+)\s+'
    r'pmax\s+([-+0-9.eEdD]+)', re.I)


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


def mesh_meta(path):
    c = CDF2(path)
    d = dict(zip(c.dims, c.sizes))
    return {'num_dim': d.get('num_dim'), 'num_nodes': d.get('num_nodes'),
            'sets': list(c.read_static('eb_prop1'))}


def quantile(v, q):
    """A robust bulk observable: pmax alone is a LOCAL mass-matrix overshoot."""
    s = sorted(v)
    return s[min(len(s) - 1, int(q * len(s)))]


def exo_fields(path):
    c = CDF2(path)
    names = var_names(c)
    pmax, pbulk, p90, pmean = [], [], [], []
    for i in range(c.numrecs):
        v = list(c.read_record('vals_nod_var1', i))
        pmax.append(max(v))
        pbulk.append(quantile(v, 0.50))
        p90.append(quantile(v, 0.90))
        pmean.append(sum(v) / len(v))
    vals = pmax + pbulk + p90
    return {'numrecs': c.numrecs, 'names': names,
            'pmax': pmax, 'pbulk': pbulk, 'p90': p90, 'pmean': pmean,
            'finite': all(math.isfinite(x) for x in vals)}


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


def block_values(text, block):
    """Minimal YAML-subset reader for a top-level block's numeric scalars."""
    vals = {}
    active = False
    for line in text.splitlines():
        code = line.split('#', 1)[0].rstrip()
        if not code.strip():
            continue
        if not line[0].isspace():
            active = (code.strip().split(':', 1)[0].strip() == block)
            continue
        if active and ':' in code:
            k, v = code.strip().split(':', 1)
            try:
                vals[k.strip()] = float(v.strip())
            except ValueError:
                pass
    return vals


def g45_rows(log):
    out = []
    for m in G45.finditer(log):
        out.append((int(m.group(1)), num(m.group(2)),
                    num(m.group(3)), num(m.group(4))))
    return out


def g32_rows(log):
    return [(int(m.group(1)), num(m.group(2)), num(m.group(3)))
            for m in G32.finditer(log)]


def linfit(ts, ys):
    tb = sum(ts) / len(ts)
    yb = sum(ys) / len(ys)
    den = sum((t - tb) ** 2 for t in ts)
    if den <= 0.0:
        return 0.0, yb
    m = sum((t - tb) * (y - yb) for t, y in zip(ts, ys)) / den
    return m, yb - m * tb


def fatal_in(log):
    return [s for s in FATAL if s in log]


# ---------------------------------------------------------------------------
# Verdicts
# ---------------------------------------------------------------------------

def check_v1(d):
    issues = []
    for tag in ('scale', 'inject'):
        if len(d['log_' + tag]) < 1000:
            issues.append('%s solver log too short' % tag)
        e = d['exo_' + tag]
        if e['numrecs'] != N_REC:
            issues.append('%s records %d != %d' % (tag, e['numrecs'], N_REC))
        if not e['finite']:
            issues.append('%s EXO pressure is non-finite' % tag)
        for key in ('pmax', 'pbulk', 'p90'):
            if key not in e or len(e[key]) < N_REC:
                issues.append('%s EXO lost the %s field' % (tag, key))
        for nm in ('Pressure', 'Damage', 'Displacement_X', 'Displacement_Y'):
            if nm not in e['names']:
                issues.append('%s EXO lost %s' % (tag, nm))
        # CI #124: an 11-step run APPENDED to the 51-record production mesh.
        # Even if numrecs happened to match, a stale/wrong result mesh is
        # caught by comparing EXO pmax with the solver's own Gate 3.2 rows.
        logp = [r[2] for r in g32_rows(d['log_' + tag])]
        if len(logp) < N_REC:
            issues.append('%s log has only %d Gate 3.2 rows' % (tag, len(logp)))
        elif 'pmax' in e and len(e['pmax']) >= N_REC:
            worst = max(abs(a - b) / max(P0, abs(b))
                        for a, b in zip(e['pmax'][:N_REC], logp[:N_REC]))
            if worst > TOL_CROSSCHECK:
                issues.append('%s EXO pmax vs Gate 3.2 log rel=%.3g: stale or '
                              'wrong result mesh' % (tag, worst))
        for sig in fatal_in(d['log_' + tag]):
            issues.append('fatal in %s log: %s' % (tag, sig))
    m = d['mesh']
    if m['num_dim'] != 2:
        issues.append('mesh num_dim=%s != 2' % m['num_dim'])
    if m['sets'] != [1, 10, 20]:
        issues.append('mesh cell sets %s != [1,10,20]' % m['sets'])
    if not m['num_nodes'] or m['num_nodes'] != 356:
        issues.append('mesh node count %s != 356' % m['num_nodes'])
    return pass_fail(not issues, '; '.join(issues) if issues else
                     '%d + %d clean records, EXO pmax == Gate 3.2 log, 4 nodal '
                     'vars, 356-node tri3 mesh, no fatal'
                     % (d['exo_scale']['numrecs'], d['exo_inject']['numrecs']))


def check_v2(d):
    issues = []
    t = block_values(d['deck_scale'], 'time')
    fl = block_values(d['deck_scale'], 'flow')
    cs = block_values(d['deck_inject'], 'cs0001')
    cs0 = block_values(d['deck_scale'], 'cs0001')
    k = fl.get('perm_min')
    kmax = fl.get('perm_max')
    S = fl.get('storativity')
    mu = fl.get('viscosity')
    f = cs.get('source')
    if not all(v is not None for v in (k, kmax, S, mu, f)):
        return pass_fail(False, 'cannot read k/perm_max/S/mu/source from the decks')
    if t.get('numstep') != N_REC:
        issues.append('numstep %s != %d' % (t.get('numstep'), N_REC))
    cv = k / (mu * S)
    T = cv * 1.0 / (H * H)
    d['metrics']['T'] = T
    d['metrics']['cv'] = cv
    d['metrics']['S'] = S
    d['metrics']['f'] = f
    d['metrics']['q_total'] = f * AREA
    if abs(T / T_TARGET - 1.0) > 0.05:
        issues.append('T=%.4g deviates %.1f%% from the field value %.4g'
                      % (T, 100 * abs(T / T_TARGET - 1.0), T_TARGET))
    if not (1.0e-8 <= T <= 1.0e-5):
        issues.append('T=%.4g outside the storage-dominated window' % T)
    if k / K_DECK >= 1.0e-6:
        issues.append('k is only %.2g x the production deck (need <1e-6)' % (k / K_DECK))
    if abs(S / S_DECK - S_RATIO_TARGET) > 0.02:
        issues.append('S/S_deck=%.4f != %.4f' % (S / S_DECK, S_RATIO_TARGET))
    if abs(kmax / k - 20.0) > 0.01:
        issues.append('perm_max/perm_min=%.3f != 20' % (kmax / k))
    if abs(f / F_SOURCE - 1.0) > 0.01:
        issues.append('injected source %.4g != %.4g' % (f, F_SOURCE))
    if cs0.get('source') != 0.0:
        issues.append('control deck source is not zero')
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'T=%.4g, c_v=%.4g m2/s, k/k_deck=%.2g, S/S_deck=%.3f, '
                     'Q_total=%.4g m2/s per unit thickness'
                     % (T, cv, k / K_DECK, S / S_DECK, f * AREA))


def check_v3(d):
    """BULK pressure must stay frozen; pmax is a local overshoot, not physics."""
    issues = []
    p = d['exo_scale']['pbulk'][:N_REC]
    if len(p) < N_REC:
        return pass_fail(False, 'not enough control records')
    if abs(p[0] - P0) / P0 > TOL_FREEZE:
        issues.append('control p_bulk(0)=%.6g != p0=%.6g' % (p[0], P0))
    worst = max(abs(x - P0) / P0 for x in p)
    if worst > TOL_FREEZE:
        issues.append('control bulk pressure drifted %.2f%% (>%.0f%%): diffusion was not frozen'
                      % (100 * worst, 100 * TOL_FREEZE))
    if max(p) - min(p) > TOL_FREEZE * P0:
        issues.append('control bulk pressure is not flat over the 11 steps')
    peak = max(d['exo_scale']['pmax'][:N_REC])
    d['metrics']['ctrl_peak'] = peak
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'p_bulk=%.6g..%.6g Pa, max drift %.3f%% (storage-dominated, '
                     'no leak-off); local pmax=%.6g Pa (mass-matrix overshoot)'
                     % (min(p), max(p), 100 * worst, peak))


def check_v4(d):
    """BULK ramp p_bulk = p0 + (f/S)t; pmax is only a LOCALITY diagnostic.

    CI #124 taught us that pmax sits on a handful of overshooting nodes when
    the consistent mass matrix dominates the (nearly frozen) Darcy operator.
    The physical observable is the median pbulk, corroborated independently by
    Gate 4.5's W = int f*p dOmega.
    """
    issues = []
    pi = d['exo_inject']['pbulk'][:N_REC]
    pc = d['exo_scale']['pbulk'][:N_REC]
    if len(pi) < N_REC or len(pc) < N_REC:
        return pass_fail(False, 'not enough records for the source response')
    delta = [a - b for a, b in zip(pi, pc)]
    # V4 deliberately uses the FROZEN constants, not the values V2 read back
    # from the decks: the physical calibration must be recoverable from the
    # solver artefacts alone.  A deck edit therefore fails V2/V7, not V4.
    f, S = F_SOURCE, S_NUM
    if abs(delta[0]) > TOL_FREEZE * P0:
        issues.append('delta_p_bulk(0)=%.6g is not zero' % delta[0])
    if any(delta[i] < delta[i - 1] - 1.0 for i in range(1, N_REC)):
        issues.append('delta_p_bulk is not monotonically increasing')
    ts = [i * DT for i in range(N_REC)]
    slope, intercept = linfit(ts, delta)
    resid = max(abs(y - (slope * t + intercept)) for t, y in zip(ts, delta))
    d['metrics']['slope'] = slope
    d['metrics']['S_eff'] = (f / slope) if slope > 0 else -1.0
    d['metrics']['resid'] = resid
    if slope <= 0.0:
        issues.append('bulk pressure ramp slope %.6g <= 0' % slope)
    else:
        S_eff = f / slope
        if abs(S_eff / S - 1.0) > TOL_SLOPE:
            issues.append('S_eff=f/slope=%.4g vs deck S=%.4g (%.1f%% off)'
                          % (S_eff, S, 100 * abs(S_eff / S - 1.0)))
        if abs(slope / SLOPE_TARGET - 1.0) > TOL_SLOPE:
            issues.append('bulk slope %.6g vs f/S=%.6g' % (slope, SLOPE_TARGET))
        # The [Pa/s] misreading would give slope ~ f, i.e. slope/f ~ 1.
        if slope / f < 1.0e6:
            issues.append('slope/f=%.3g: source behaves like [Pa/s], not [1/s]'
                          % (slope / f))
    if max(delta) > 0 and resid / max(delta) > TOL_LINEAR:
        issues.append('bulk response is not linear (max residual %.2f%% of peak)'
                      % (100 * resid / max(delta)))
    # The overshoot must stay LOCAL (p90 ~ pbulk) and be a FIXED spatial mode
    # (its own ramp is proportional to the bulk ramp), not a diverging solve.
    p90 = d['exo_inject']['p90'][:N_REC]
    pmax = d['exo_inject']['pmax'][:N_REC]
    dev90 = max(abs(q - b) / max(P0, abs(b)) for q, b in zip(p90, pi))
    d['metrics']['p90_dev'] = dev90
    if dev90 > TOL_PEAK_LOCAL:
        issues.append('overshoot is GLOBAL: p90 departs %.2f%% from the bulk '
                      '(>%.0f%%)' % (100 * dev90, 100 * TOL_PEAK_LOCAL))
    peak_slope, _ = linfit(ts, [x - P0 for x in pmax])
    ratio_pk = (peak_slope / slope) if slope > 0 else -1.0
    d['metrics']['peak_ratio'] = ratio_pk
    if slope > 0 and not (PEAK_RATIO_BAND[0] <= ratio_pk <= PEAK_RATIO_BAND[1]):
        issues.append('pmax/bulk ramp ratio %.3g outside %s: peak is not a '
                      'fixed local mode' % (ratio_pk, PEAK_RATIO_BAND))
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'bulk slope=%.6g Pa/s (f/S=%.6g), S_eff=%.4g, linear '
                     'residual %.2f%%, slope/f=%.3g; local pmax ramp x%.3f '
                     '(p90 dev %.3f%%)'
                     % (slope, SLOPE_TARGET, d['metrics']['S_eff'],
                        100 * resid / max(delta) if max(delta) > 0 else 0.0,
                        slope / f, ratio_pk, 100 * dev90))


def check_v5(d):
    issues = []
    rs, ri = g45_rows(d['log_scale']), g45_rows(d['log_inject'])
    if len(rs) != N_REC or len(ri) != N_REC:
        issues.append('Gate 4.5 rows %d/%d != %d' % (len(rs), len(ri), N_REC))
        return pass_fail(False, '; '.join(issues))
    if [r[0] for r in ri] != list(range(1, N_REC + 1)):
        issues.append('injection step numbers are not 1..%d' % N_REC)
    ws = [r[3] for r in rs]
    wi = [r[3] for r in ri]
    if any(abs(w) > 1.0e-12 for w in ws):
        issues.append('zero-source control reported non-zero W')
    if any(w <= 0.0 for w in wi):
        issues.append('injection work W is not strictly positive')
    if any(wi[i] < wi[i - 1] - 1.0e-9 for i in range(1, len(wi))):
        issues.append('injection work W is not monotone')
    if any(not math.isfinite(r[2]) or r[2] < 0.0 for r in ri):
        issues.append('dissipation E is non-finite or negative')
    # Loose magnitude band against the frozen analytical expectation; the
    # exact cell-wise quadrature is left to Gate 4.5 / 4.8.
    w_an = F_SOURCE * AREA * (P0 + SLOPE_TARGET * 1.0)
    ratio = wi[-1] / w_an
    d['metrics']['w_last'] = wi[-1]
    d['metrics']['w_ratio'] = ratio
    if not (0.05 <= ratio <= 20.0):
        issues.append('W_last/W_an=%.4g outside the loose magnitude band' % ratio)
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'control W=0, injection W=%.4g..%.4g (monotone), W/W_an=%.3f'
                     % (wi[0], wi[-1], ratio))


def check_v6(d):
    issues = []
    for tag in ('scale', 'inject'):
        log = d['log_' + tag]
        for label, needle, need in (
                ('Gate 3.1:', 'Gate 3.1: step', N_REC),
                ('Gate 3.2:', 'Gate 3.2: step', N_REC),
                ('Gate 4.12:', 'Gate 4.12: step', N_REC),
                ('Gate 4.5:', 'Gate 4.5: step', N_REC)):
            n = log.count(needle)
            if n < need:
                issues.append('%s %s appears only %d times' % (tag, label, n))
        n_solves = log.count('Linear solve converged')
        if n_solves < N_REC - 1:
            issues.append('%s has only %d converged flow solves (>=%d required)'
                          % (tag, n_solves, N_REC - 1))
        for sig in fatal_in(log):
            issues.append('%s fatal signature: %s' % (tag, sig))
    wf = d['src']['workflow']
    for token in ('timeout -k 10s 600', '-flow_ksp_monitor',
                  '-flow_ksp_converged_reason'):
        if token not in wf:
            issues.append('build.yml lost solver observability: %s' % token)
    return pass_fail(not issues, '; '.join(issues) if issues else
                     '11 diagnostic chains + >=10 converged flow solves per deck, '
                     'timeout/KSP monitor wired, no fatal')


def check_v7(d):
    issues = []
    s = d['src']
    if 'Volumetric source/sink rate (1/s)' not in s['flowctx']:
        issues.append('flow.source registration lost its [1/s] unit')
    for token in ('MEF90Scaling_CST', 'MEF90FlowSetSourceCst',
                  'MEF90DiffusionWorkSetCell'):
        if token not in s['flowasm']:
            issues.append('m_MEF90_FlowAssembly.F90 lost %s' % token)
    if 'Call MEF90FlowEnergy(' not in s['vdef']:
        issues.append('vDef.F90 no longer emits the Gate 4.5 energy row')
    if "source_scaling: constant" not in d['deck_inject']:
        issues.append('injection deck lost source_scaling: constant')

    d1 = code_diff(d['deck_main'], d['deck_scale'])
    if d1 is None or len(d1) != 4:
        issues.append('main->scale code diff is %s (expected 4)' %
                      (len(d1) if d1 is not None else 'length mismatch'))
    elif sorted(key_of(x) for x, y in d1) != ['numstep', 'perm_max',
                                             'perm_min', 'storativity']:
        issues.append('main->scale changed non-similitude keys')
    d2 = code_diff(d['deck_scale'], d['deck_inject'])
    if d2 is None or len(d2) != 1:
        issues.append('scale->inject code diff is %s (expected 1)' %
                      (len(d2) if d2 is not None else 'length mismatch'))
    elif key_of(d2[0][0]) != 'source':
        issues.append('the single scale->inject diff is not cs0001.source')

    wf = s['workflow']
    for token in ('check_gate415_rate_fracability.py --gate',
                  'gate415_verdict.log', 'terzaghi415_scale.yaml',
                  'terzaghi415_inject.yaml', 'out_415_scale.gen',
                  'out_415_inject.gen', 's415_scale_rc', 's415_inject_rc',
                  # CI #124: without these the 11-step EXO appends to the
                  # 51-record production mesh and V1/V3 read stale records.
                  'rm -f data/Beam1x10-tri3_out.gen',
                  'g415_rc',
                  'if [ "$s415_scale_rc" -ne 0 ]',
                  'if [ "$s415_inject_rc" -ne 0 ]',
                  'if [ "$g415_rc" -ne 0 ]'):
        if token not in wf:
            issues.append('build.yml missing %s' % token)
    if wf.count('rm -f data/Beam1x10-tri3_out.gen') < 2:
        issues.append('build.yml must rebuild a clean EXO before BOTH 415 runs')
    for bad in ('batched_assembly', 'biot_effective_stress',
                'temp_dep_viscosity'):
        if bad in d['deck_inject'] or bad in d['deck_scale']:
            issues.append('a 415 deck enables non-default switch %s' % bad)
    return pass_fail(not issues, '; '.join(issues) if issues else
                     '[1/s] unit + source/energy wiring, 4-line + 1-line deck '
                     'contracts, CI run/verdict/deferred-fail wiring, default physics')


def check_v8(d):
    issues = []
    log = d['log_main']
    rows = g32_rows(log)
    if len(rows) < 50:
        issues.append('production log has only %d Gate 3.2 rows' % len(rows))
    else:
        pm = [r[2] for r in rows]
        if abs(pm[0] - P0) / P0 > TOL_FREEZE:
            issues.append('production p(0)=%.6g != p0' % pm[0])
        if pm[-1] >= 0.5 * P0:
            issues.append('production pressure did not consolidate (end %.6g)' % pm[-1])
    wm = [r[3] for r in g45_rows(log)]
    if wm and any(abs(w) > 1.0e-12 for w in wm):
        issues.append('production injection work is non-zero (source leaked in)')
    if 'terzaghi415' in log or 'Gate 4.15:' in log:
        issues.append('a 4.15 diagnostic leaked into the production log')
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'production still consolidates from p0 over %d steps, W=0, '
                     'no 415 leakage' % len(rows))


CHECKS = (('V1', 'run / EXO / mesh data integrity', check_v1),
          ('V2', 'similitude parameters frozen', check_v2),
          ('V3', 'zero-source control frozen (storage-dominated)', check_v3),
          ('V4', 'injection slope f/S + source-unit calibration', check_v4),
          ('V5', 'Gate 4.5 injection work', check_v5),
          ('V6', 'solver chain + KSP convergence', check_v6),
          ('V7', 'static anti-regression + deck/CI contracts', check_v7),
          ('V8', 'production isolation', check_v8))


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
    print('GATE415_RESULT: %s%s' % ('PASS' if not fails else 'FAIL',
                                    '' if not fails else ' (' + ','.join(fails) + ')'))
    for tag, label, verdict, detail in rows:
        if verdict == 'FAIL':
            print('  %s detail: %s' % (tag, detail))


def synth_sources():
    return {k: read(v) for k, v in SRC.items()}


def synth_decks():
    return {k: read(v) for k, v in DECKS.items()}


def make_log(pmax, work, n=11, prefix=''):
    log = prefix + 'x' * 1200 + '\n'
    for i in range(1, n + 1):
        log += 'Gate 3.1: step %d\n' % i
        log += 'Gate 3.2: step %d pmin 0.0 pmax %.10E\n' % (i, pmax[i - 1])
        log += 'Gate 4.12: step %d t %.4g dpmax %.4g dpmin 0.0\n' % (
            i, (i - 1) * 0.1, pmax[i - 1] - (pmax[i - 2] if i > 1 else pmax[0]))
        log += 'Gate 4.5: step %d t %.10E E %.10E W %.10E\n' % (
            i, (i - 1) * 0.1, 0.125, work[i - 1])
        if i > 1:
            log += '    Linear solve converged due to CONVERGED_RTOL iterations 1\n'
    return log


def synth():
    """Offline positive case: real decks/sources, synthetic solver artefacts.

    The positive case deliberately carries the CI #124 LOCAL pmax overshoot
    (2.9184x the bulk ramp) so that a PASS proves the verdict is anchored on
    pbulk and is NOT fooled by the mass-matrix spike.
    """
    p_ctrl = [P0] * N_REC
    p_inj = [P0 + SLOPE_TARGET * i * DT for i in range(N_REC)]          # bulk
    peak_inj = [P0 + 2.9184 * SLOPE_TARGET * i * DT for i in range(N_REC)]
    p90_inj = [P0 + 1.0038 * SLOPE_TARGET * i * DT for i in range(N_REC)]
    peak_ctrl = [P0 + 3.67 * i for i in range(N_REC)]
    w_ctrl = [0.0] * N_REC
    w_inj = [F_SOURCE * AREA * p for p in p_inj]
    e = {'numrecs': N_REC, 'names': ['Pressure', 'Damage',
                                     'Displacement_X', 'Displacement_Y'],
         'finite': True}
    exo_c = dict(e, pmax=peak_ctrl, pbulk=p_ctrl, p90=p_ctrl, pmean=p_ctrl)
    exo_i = dict(e, pmax=peak_inj, pbulk=p_inj, p90=p90_inj, pmean=p_inj)
    p_main = [P0 * (0.0906 ** (j / 50.0)) for j in range(51)]
    log_main = make_log(p_main, [0.0] * 51, n=51,
                        prefix='Gate 3.1: transient Darcy TS created\n')
    d = {
        'log_scale': make_log(peak_ctrl, w_ctrl),
        'log_inject': make_log(peak_inj, w_inj),
        'log_main': log_main,
        'exo_scale': exo_c, 'exo_inject': exo_i,
        'mesh': {'num_dim': 2, 'num_nodes': 356, 'sets': [1, 10, 20]},
        'src': synth_sources(),
    }
    d.update({'deck_' + k: v for k, v in synth_decks().items()})
    return d


def selftest():
    import copy
    base = synth()
    cases = []
    rows, fails = run(base)
    cases.append(('POS full B0 verdict', rows, set(), None))
    print('POS  %-52s %s' % ('full B0 verdict', 'PASS' if not fails else 'FAIL ' + ','.join(fails)))

    def neg(name, mutate, expect, label='NEG'):
        d = copy.deepcopy(base)
        mutate(d)
        r, f = run(d)
        ok = set(f) == set(expect)
        cases.append((name, r, set(expect), f))
        print('%s  %-52s %s (failed %s, expected %s)'
              % (label, name, 'OK' if ok else 'MISMATCH',
                 sorted(f), sorted(expect)))
        return ok

    def inv(name, mutate):
        """A permitted perturbation: the verdict must stay fully green."""
        return neg(name, mutate, set(), 'INV')

    allok = not fails
    allok &= neg('truncated injection EXO (10 records)',
                 lambda d: d['exo_inject'].update(numrecs=10), {'V1'})
    allok &= neg('fatal signature in scale log',
                 lambda d: d.update(log_scale=d['log_scale'] + 'PETSC ERROR: boom\n'),
                 {'V1', 'V6'})
    allok &= neg('k_num reverted to the production value',
                 lambda d: d.update(deck_scale=d['deck_scale'].replace(
                     'perm_min: 2.35e-18', 'perm_min: 5.0e-11')), {'V2', 'V7'})
    allok &= neg('S_num reverted to the production value',
                 lambda d: d.update(deck_scale=d['deck_scale'].replace(
                     'storativity: 1.50e-10', 'storativity: 4.6e-10')), {'V2', 'V7'})
    allok &= neg('numstep reverted to 51',
                 lambda d: d.update(deck_scale=d['deck_scale'].replace(
                     'numstep: 11', 'numstep: 51')), {'V2', 'V7'})
    allok &= neg('injection source silently zeroed',
                 lambda d: d.update(deck_inject=d['deck_inject'].replace(
                     'source: 1.5e-5', 'source: 0.0')), {'V2', 'V7'})
    allok &= neg('control deck pressure also ramps (diffusion not frozen)',
                 lambda d: d['exo_scale'].update(pbulk=d['exo_inject']['pbulk'][:]),
                 {'V3', 'V4'})
    allok &= neg('injection pressure frozen at p0 (source unwired)',
                 lambda d: d['exo_inject'].update(pbulk=[P0] * N_REC), {'V4'})
    allok &= neg('ramp slope halved (wrong storativity)',
                 lambda d: d['exo_inject'].update(
                     pbulk=[P0 + 0.5 * SLOPE_TARGET * i * DT for i in range(N_REC)]),
                 {'V4'})
    allok &= neg('slope equals f, i.e. the [Pa/s] misreading',
                 lambda d: d['exo_inject'].update(
                     pbulk=[P0 + F_SOURCE * i * DT for i in range(N_REC)]), {'V4'})
    allok &= neg('local overshoot becomes GLOBAL (p90 follows pmax)',
                 lambda d: d['exo_inject'].update(p90=d['exo_inject']['pmax'][:]),
                 {'V4'})

    def stale_exo(d):
        """Replay CI #124: 11 real records + 40 stale production records."""
        tail = [P0 * (0.0906 ** (j / 50.0)) for j in range(11, 51)]
        for tag in ('scale', 'inject'):
            e = d['exo_' + tag]
            e['numrecs'] = 51
            for key in ('pmax', 'pbulk', 'p90', 'pmean'):
                e[key] = list(e[key]) + tail

    allok &= neg('CI #124 stale 51-record EXO not rebuilt', stale_exo, {'V1'})
    allok &= neg('EXO pmax disagrees with the Gate 3.2 solver log',
                 lambda d: d['exo_inject']['pmax'].__setitem__(
                     2, d['exo_inject']['pmax'][2] * 0.5), {'V1'})
    allok &= inv('local pmax spike removed, bulk physics unchanged',
                 lambda d: (d['exo_inject'].update(
                     pmax=d['exo_inject']['pbulk'][:]),
                     d.update(log_inject=make_log(
                         d['exo_inject']['pbulk'],
                         [F_SOURCE * AREA * p
                          for p in d['exo_inject']['pbulk']]))))
    allok &= neg('build.yml EXO rebuild before both 415 runs removed',
                 lambda d: d['src'].update(workflow=d['src']['workflow'].replace(
                     'rm -f data/Beam1x10-tri3_out.gen', 'true # no rebuild')),
                 {'V7'})
    allok &= neg('zero-source control reports non-zero injection work',
                 lambda d: d.update(log_scale=d['log_scale'].replace('W 0.0000000000E+00',
                                                                    'W 3.5000000000E+00')),
                 {'V5'})
    allok &= neg('injection work identically zero',
                 lambda d: d.update(log_inject=d['log_inject'].replace(
                     re.search(r'W [0-9.]+E[+-][0-9]+', d['log_inject']).group(0),
                     'W 0.0000000000E+00', 1)), {'V5'})
    allok &= neg('injection work not monotone',
                 lambda d: d.update(log_inject=d['log_inject'].replace(
                     'W 3.0000000000E+01', 'W 1.0000000000E-03', 1)), {'V5'})
    allok &= neg('flow linear solves never converged (timeout signature)',
                 lambda d: d.update(log_inject=d['log_inject'].replace(
                     'Linear solve converged', 'Linear solve stalled')), {'V6'})
    allok &= neg('source unit registration changed back to [Pa/s]',
                 lambda d: d['src'].update(flowctx=d['src']['flowctx'].replace(
                     'Volumetric source/sink rate (1/s)',
                     'Volumetric source/sink rate (Pa/s)')), {'V7'})
    allok &= neg('MEF90FlowSetSourceCst removed',
                 lambda d: d['src'].update(flowasm=d['src']['flowasm'].replace(
                     'MEF90FlowSetSourceCst', 'MEF90FlowSetSourceDisabled')), {'V7'})
    allok &= neg('deferred s415_inject_rc hard fail removed',
                 lambda d: d['src'].update(workflow=d['src']['workflow'].replace(
                     'if [ "$s415_inject_rc" -ne 0 ]', 'if false')), {'V7'})
    allok &= neg('injection deck carries a second code change',
                 lambda d: d.update(deck_inject=d['deck_inject'].replace(
                     'p0: 1.0e5', 'p0: 2.0e5')), {'V7'})
    allok &= neg('415 prefix leaked into the production log',
                 lambda d: d.update(log_main=d['log_main'] +
                                    '-options_file_yaml data/terzaghi415_inject.yaml\n'), {'V8'})
    allok &= neg('production consolidation disabled',
                 lambda d: d.update(log_main=d['log_main'].replace(
                     'pmax 9.0600000000E+03', 'pmax 1.0000000000E+05')), {'V8'})
    print('SELFTEST_RESULT: %s (%d positive + %d mutations, incl. 1 invariance)'
          % ('PASS' if allok else 'FAIL', 1, len(cases) - 1))
    return 0 if allok else 1


def main(argv):
    if argv and argv[0] == '--selftest':
        return selftest()
    if argv and argv[0] == '--gate':
        if len(argv) != 7:
            print('GATE415_RESULT: FAIL (usage: --gate <log_scale> <exo_scale> '
                  '<log_inject> <exo_inject> <mesh2d> <log_main>)')
            return 1
        try:
            d = {'log_scale': read(argv[1]), 'log_inject': read(argv[3]),
                 'log_main': read(argv[6]),
                 'exo_scale': exo_fields(argv[2]),
                 'exo_inject': exo_fields(argv[4]),
                 'mesh': mesh_meta(argv[5]),
                 'src': synth_sources(),
                 'deck_main': read(DECKS['main']),
                 'deck_scale': read(DECKS['scale']),
                 'deck_inject': read(DECKS['inject'])}
        except Exception as exc:
            print('=== Gate 4.15 / Stage B0 (real CI logs/EXO) ===')
            print('FAIL V1 run / EXO / mesh data integrity | cannot parse artefacts: %r' % exc)
            print('GATE415_RESULT: FAIL (V1 data error: %r)' % exc)
            print('  hint: a truncated out_415_*.gen usually means the run was '
                  'killed by the timeout before its first flow solve converged')
            return 1
        rows, fails = run(d)
        report(rows, fails, '=== Gate 4.15 / Stage B0: rate/fracability '
               'similitude + source-unit calibration (real CI logs/EXO) ===')
        m = d['metrics']
        print('g415_T                     = %.6e' % m.get('T', -1))
        print('g415_cv                    = %.6e' % m.get('cv', -1))
        # g415_slope_measured is now the BULK slope (key kept for continuity);
        # the local pmax overshoot is reported separately, never used as f/S.
        print('g415_slope_measured        = %.6e' % m.get('slope', -1))
        print('g415_S_eff                 = %.6e' % m.get('S_eff', -1))
        print('g415_peak_ratio            = %.6e' % m.get('peak_ratio', -1))
        print('g415_p90_dev               = %.6e' % m.get('p90_dev', -1))
        print('g415_q_total_per_thickness = %.6e' % m.get('q_total', -1))
        print('g415_w_last                = %.6e' % m.get('w_last', -1))
        return 0 if not fails else 1
    print('usage: check_gate415_rate_fracability.py [--gate <log_scale> '
          '<exo_scale> <log_inject> <exo_inject> <mesh2d> <log_main> | --selftest]')
    return 2


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
