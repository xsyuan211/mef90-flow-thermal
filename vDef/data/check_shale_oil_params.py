#!/usr/bin/env python3
"""Shale-oil literature parameter recomputation + falsifiable M1-M8 verdict.

Zero Fortran / zero deck / zero solver changes.  This is a pure post-processing
checker over `shale_oil_params.json`, the first-edition tables A-C in
RATE_FRACABILITY_ROADMAP.md 4.2, and the frozen Gate 4.15 B0 constants.

Three independent formula implementations are compared:

  REFERENCE  the literal printed formulas of roadmap 2.3, hardcoded here and
             deliberately independent of the JSON `formulas` contract;
  LITERAL    the same radicals, but using the JSON inner constants;
  EXPANDED   the power-law form built ONLY from the JSON exponents.

M2 requires all three to agree to machine precision.  Therefore any silent
drift of a prefactor or an exponent is falsifiable, not just cosmetic.

Usage:
  python3 data/check_shale_oil_params.py                 # real files -> verdict
  python3 data/check_shale_oil_params.py --md            # markdown tables A'-C'
  python3 data/check_shale_oil_params.py --json out.json # metric dump
  python3 data/check_shale_oil_params.py --selftest      # positive + mutations
"""

import argparse
import copy
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
PARAM_PATH = os.path.join(HERE, 'shale_oil_params.json')
ROADMAP_PATH = os.path.join(ROOT, 'RATE_FRACABILITY_ROADMAP.md')
B0_DECK_PATH = os.path.join(HERE, 'terzaghi415_inject.yaml')
BUILD_PATH = os.path.join(ROOT, '.github', 'workflows', 'build.yml')

# Published tables A-C carry only 2-3 significant digits.  0.0014 has two,
# so a correctly recomputed 0.0013552 is 3.2% away purely from rounding.
# 4% still rejects every exponent / order-of-magnitude drift by a wide margin.
TOL_LEGACY = 0.04
TOL_MACHINE = 1.0e-10    # literal vs expanded vs reference
TOL_BAG = 1.0e-12        # frozen B0 constants
RATE_WINDOW_DISCLAIMER = (
    'The literature parameter recomputation does NOT produce a field rate '
    'optimization table. B0/B1/B2 are pre-initiation, storage-dominated, '
    'non-converged-damage experiments; beta = f/S is the only identifiable '
    'quantity on p(t) there. K and C are propagation-regime dimensionless '
    'numbers and are NOT identifiable before initiation (Gate 4.15 B1 A4). '
    'Any field rate recommendation requires DFIT / mini-frac calibration and '
    'a genuine fracture-propagation window (Stage B3).'
)

# Exponent/prefactor contract derived from roadmap 2.3.  Hardcoded on purpose:
# this is the expected value, NOT read from the JSON under test.
EXPECTED = {
    'mu_prime_factor': 12.0,
    'K_prefactor': 4.0 / math.sqrt(math.pi),
    'K_inner_constant': 3.0,
    'K_inner_exp': -0.25,
    'K_exp_q': -0.25,
    'K_exp_Eprime': -0.75,
    'K_exp_mu_prime': -0.25,
    'C_prefactor': 2.0,
    'C_inner_constant': 12.0,
    'C_inner_exp': -1.0 / 6.0,
    'C_exp_q': -0.5,
    'C_exp_t': 1.0 / 6.0,
    'C_exp_mu': -1.0 / 6.0,
    'C_exp_Eprime': 1.0 / 6.0,
    'CL_exp_k': 0.5,
    'CL_exp_S': 0.5,
    'CL_exp_mu': -0.5,
    'T_exp_k': 1.0,
    'T_exp_mu': -1.0,
    'T_exp_S': -1.0,
    'T_exp_H': -2.0,
    'T_exp_t': 1.0,
}

PARAM_KEYS = (
    'E_static', 'nu_static', 'KIC', 'k_core', 'k_eff', 'phi', 'c_pore',
    'c_fluid', 'S', 'mu', 'delta_p', 'delta_sigma', 'H_f', 't_f', 'Q',
    'clusters_per_stage', 'cluster_spacing', 'fluid_per_stage',
    'proppant_per_stage', 'horizontal_length', 'well_spacing',
    'pressure_coefficient', 'brittleness_index', 'KIC_anisotropy_ratio',
    'closure_stress_min',
)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1.0e-300)


def close(a, b, tol):
    return rel(a, b) <= tol


# ---------------------------------------------------------------------------
# REFERENCE implementation (independent of the JSON formula contract)
# ---------------------------------------------------------------------------

def ref_Eprime(E_GPa, nu):
    return E_GPa * 1.0e9 / (1.0 - nu * nu)


def ref_Gc(KIC_MPa, Eprime):
    return (KIC_MPa * 1.0e6) ** 2 / Eprime


def ref_sigma_c(Eprime, Gc, ell):
    return math.sqrt(Eprime * Gc / ell)


def ref_cv(k, mu, S):
    return k / (mu * S)


def ref_T(k, t, mu, S, H):
    return k * t / (mu * S * H * H)


def ref_CL(dp, k, S, mu):
    return dp * math.sqrt(k * S / (math.pi * mu))


def ref_K(KIC, q, Eprime, mup):
    return (4.0 * KIC / math.sqrt(math.pi)) * \
           (1.0 / (3.0 * q * Eprime ** 3 * mup)) ** 0.25


def ref_C(CL, Eprime, t, mu, q):
    return 2.0 * CL * (Eprime * t / (12.0 * mu * q ** 3)) ** (1.0 / 6.0)


def ref_knum(T, mu_num, S_num, Hn, tn):
    return T * mu_num * S_num * Hn * Hn / tn


# ---------------------------------------------------------------------------
# LITERAL / EXPANDED implementations driven by the JSON contract
# ---------------------------------------------------------------------------

def lit_K(KIC, q, Eprime, mup, F):
    return F['K_prefactor'] * KIC * \
           (1.0 / (F['K_inner_constant'] * q * Eprime ** 3 * mup)) ** (-F['K_exp_q'])


def exp_K(KIC, q, Eprime, mup, F):
    return (F['K_prefactor'] * F['K_inner_constant'] ** F['K_inner_exp'] * KIC *
            q ** F['K_exp_q'] * Eprime ** F['K_exp_Eprime'] *
            mup ** F['K_exp_mu_prime'])


def lit_C(CL, Eprime, t, mu, q, F):
    return F['C_prefactor'] * CL * \
           (Eprime * t / (F['C_inner_constant'] * mu * q ** 3)) ** F['C_exp_t']


def exp_C(CL, Eprime, t, mu, q, F):
    return (F['C_prefactor'] * F['C_inner_constant'] ** F['C_inner_exp'] * CL *
            Eprime ** F['C_exp_Eprime'] * t ** F['C_exp_t'] *
            mu ** F['C_exp_mu'] * q ** F['C_exp_q'])


def exp_cv(k, mu, S, F):
    return k ** F['T_exp_k'] * mu ** F['T_exp_mu'] * S ** F['T_exp_S']


def exp_T(k, t, mu, S, H, F):
    return (k ** F['T_exp_k'] * mu ** F['T_exp_mu'] * S ** F['T_exp_S'] *
            H ** F['T_exp_H'] * t ** F['T_exp_t'])


def exp_CL(dp, k, S, mu, F):
    return dp * k ** F['CL_exp_k'] * S ** F['CL_exp_S'] * \
           mu ** F['CL_exp_mu'] / math.sqrt(math.pi)


# ---------------------------------------------------------------------------
# Load / compute
# ---------------------------------------------------------------------------

def read_text(path):
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as fh:
            return fh.read()
    except OSError:
        return ''


def load_payload():
    with open(PARAM_PATH, 'r', encoding='utf-8') as fh:
        params = json.load(fh)
    with open(__file__, 'r', encoding='utf-8') as fh:
        script = fh.read()
    return {
        'P': params,
        'texts': {
            'roadmap': read_text(ROADMAP_PATH),
            'deck': read_text(B0_DECK_PATH),
            'script': script,
            'build': read_text(BUILD_PATH),
        },
    }


def legacy_rows(P):
    """Reproduce first-edition tables A-C with the REFERENCE formulas."""
    W = P['window']['legacy']
    F = P['formulas']
    out = []
    for c in P['window']['legacy_classes']:
        Ep = c['E_prime_Pa']
        Gc = c['Gc_N_per_m']
        cv = c['c_v_m2_per_s']
        KIC = math.sqrt(Ep * Gc)
        k = cv * W['mu_Pa_s'] * W['S0_Pa_inv']
        mup = F['mu_prime_factor'] * W['mu_Pa_s']
        T = ref_T(k, W['t_field_s'], W['mu_Pa_s'], W['S0_Pa_inv'], W['H_field_m'])
        CL = ref_CL(W['dp_Pa'], k, W['S0_Pa_inv'], W['mu_Pa_s'])
        sig = ref_sigma_c(Ep, Gc, W['ell_m'])
        Ks, Cs = [], []
        for Qm in W['Q_scan_m3_min']:
            q = (Qm / 60.0) / W['H_field_m']
            Ks.append(ref_K(KIC, q, Ep, mup))
            Cs.append(ref_C(CL, Ep, W['t_field_s'], W['mu_Pa_s'], q))
        knum = ref_knum(T, W['mu_num_Pa_s'], W['S_num_Pa_inv'],
                        W['H_n_m'], W['t_n_s'])
        out.append({
            'name': c['name'], 'Eprime': Ep, 'Gc': Gc, 'cv': cv, 'KIC': KIC,
            'k': k, 'T': T, 'CL': CL, 'sigma_c': sig, 'K': Ks, 'C': Cs,
            'knum': knum, 'k_ratio': knum / W['k_deck_m2'],
            'S_ratio': W['S_num_Pa_inv'] / W['S_deck_Pa_inv'],
            'published': c['published'],
        })
    return out


def literature_point(P, tier):
    """One literature-consistent point at low / rep / high monotonic corner."""
    V = P['params']
    W = P['window']['literature']
    F = P['formulas']
    v = lambda name: V[name][tier]
    E = v('E_static'); nu = v('nu_static')
    Ep = ref_Eprime(E, nu)
    KIC = v('KIC') * 1.0e6
    Gc = ref_Gc(v('KIC'), Ep)
    k = v('k_eff'); S = v('S'); mu = v('mu')
    H = v('H_f'); t = v('t_f'); dp = v('delta_p')
    cv = ref_cv(k, mu, S)
    T = ref_T(k, t, mu, S, H)
    CL = ref_CL(dp, k, S, mu)
    sig = ref_sigma_c(Ep, Gc, W['ell_m'])
    mup = F['mu_prime_factor'] * mu
    Ks, Cs, qs = [], [], []
    for Qm in W['Q_scan_m3_min']:
        q = (Qm / 60.0) / H
        qs.append(q)
        Ks.append(ref_K(KIC, q, Ep, mup))
        Cs.append(ref_C(CL, Ep, t, mu, q))
    knum = ref_knum(T, W['mu_num_Pa_s'], W['S_num_Pa_inv'],
                    W['H_n_m'], W['t_n_s'])
    return {
        'Eprime': Ep, 'Gc': Gc, 'KIC': KIC, 'sigma_c': sig, 'cv': cv,
        'T': T, 'CL': CL, 'K': Ks, 'C': Cs, 'q': qs, 'knum': knum,
        'k': k, 'S': S, 'mu': mu, 'H': H, 't': t, 'dp': dp,
    }


def build_metrics(P):
    F = P['formulas']
    N = P['numerical']
    W = P['window']['literature']
    V = P['params']
    leg = legacy_rows(P)
    pts = {tier: literature_point(P, tier) for tier in ('low', 'rep', 'high')}
    rep = pts['rep']
    m = {}

    def put(k, val):
        m[k] = float(val)

    for tier, p in pts.items():
        for key in ('Eprime', 'Gc', 'KIC', 'sigma_c', 'cv', 'T', 'CL', 'knum'):
            put('so_%s_%s' % (tier, key), p[key])
        put('so_%s_K_q10' % tier, p['K'][0])
        put('so_%s_K_q14' % tier, p['K'][1])
        put('so_%s_C_q10' % tier, p['C'][0])
        put('so_%s_C_q14' % tier, p['C'][1])

    # Monotone-corner envelope of T: each input is chosen at the end that
    # actually minimises / maximises T = k*t/(mu*S*H^2).  The three
    # low/rep/high rows in the tables are CONSISTENT parameter tiers, NOT
    # extrema of any target quantity; keeping the two concepts separate is
    # essential for an honest uncertainty statement.
    pk = lambda name, tier: V[name][tier]
    T_min = ref_T(pk('k_eff', 'low'), pk('t_f', 'low'), pk('mu', 'high'),
                  pk('S', 'high'), pk('H_f', 'high'))
    T_max = ref_T(pk('k_eff', 'high'), pk('t_f', 'high'), pk('mu', 'low'),
                  pk('S', 'low'), pk('H_f', 'low'))
    put('so_T_corner_min', T_min)
    put('so_T_corner_max', T_max)
    put('so_T_span_decades', math.log10(T_max / T_min))
    put('so_knum_rep_over_b0', rep['knum'] / N['k_num_m2'])
    put('so_b0_knum_m2', N['k_num_m2'])
    # Sensitivity to the one unconstrained engineering assumption.
    hc = pts['high']
    dp10 = 1.0e7
    CL10 = ref_CL(dp10, hc['k'], hc['S'], hc['mu'])
    put('so_C_dp10_high_corner',
        ref_C(CL10, hc['Eprime'], hc['t'], hc['mu'], hc['q'][0]))
    put('so_kcore_over_keff_rep', V['k_core']['rep'] / V['k_eff']['rep'])
    put('so_S_derived_rep', V['phi']['rep'] *
        (V['c_pore']['rep'] + V['c_fluid']['rep']))
    put('so_q_field_rep_m2_per_s', (V['Q']['rep']) / rep['H'])
    put('so_q_b0_m2_per_s', N['Q_total_b0_m2_per_s'])
    put('so_q_ratio_field_over_b0',
        (V['Q']['rep'] / rep['H']) / N['Q_total_b0_m2_per_s'])
    put('so_legacy_k_field_I_m2', leg[0]['k'])
    put('so_legacy_k_field_IV_m2', leg[3]['k'])
    put('so_sigma_c_identity_check',
        rep['sigma_c'] / (rep['KIC'] / math.sqrt(W['ell_m'])))
    put('so_b0_f_times_area', N['f_b0_inv_s'] * N['Area_m2'])
    return m, leg, pts


# ---------------------------------------------------------------------------
# M1-M8
# ---------------------------------------------------------------------------

def check_m1(P, texts, m, leg, pts):
    """Parameter-table structural integrity and frozen-constant consistency."""
    b = []
    V = P.get('params', {})
    S = P.get('sources', [])
    ids = [s.get('id') for s in S]
    if len(ids) != 20 or len(set(ids)) != 20:
        b.append('expected 20 unique source ids, got %d / %d' % (len(ids), len(set(ids))))
    for key in PARAM_KEYS:
        if key not in V:
            b.append('missing param %s' % key)
    for key, v in V.items():
        if not all(x in v for x in ('unit', 'low', 'high', 'rep', 'sources', 'note')):
            b.append('param %s incomplete' % key)
            continue
        if not (v['low'] <= v['rep'] <= v['high']):
            b.append('param %s violates low<=rep<=high' % key)
        if not v['unit'] or not v['note']:
            b.append('param %s empty unit/note' % key)
        for sid in v['sources']:
            if sid not in ids:
                b.append('param %s cites unknown source %s' % (key, sid))
    honesty = P.get('honesty', '')
    for token in ('NOT block-specific', 'k_core', 'k_eff', 'DFIT'):
        if token not in honesty:
            b.append('honesty disclaimer missing %r' % token)
    if P.get('schema') != 'shale-oil-literature-params/1':
        b.append('schema version drift')
    if 'window' not in P or 'formulas' not in P or 'numerical' not in P:
        b.append('window / formulas / numerical block missing')
        return b
    W = P['window']
    if len(W['legacy_classes']) != 4:
        b.append('legacy_classes must have exactly 4 rows')
    env, pub = P['legacy_envelope'], W['legacy_published_envelope']
    common = [k for k in env if k in pub and isinstance(env[k], list)]
    for k in common:
        if len(env[k]) != 2 or any(not close(a, c, 1e-9) for a, c in zip(env[k], pub[k])):
            b.append('legacy_envelope vs published_envelope mismatch on %s' % k)
    lit = W['literature']
    for field, pname in lit['map'].items():
        if pname not in V:
            b.append('literature.map %s -> unknown param %s' % (field, pname))
    qlo, qhi = V['Q']['low'] * 60.0, V['Q']['high'] * 60.0
    scan = lit['Q_scan_m3_min']
    if not (close(qlo, scan[0], 1e-3) and close(qhi, scan[1], 1e-3)):
        b.append('Q literature window does not match Q param bounds (%.4f..%.4f vs %s)'
                 % (qlo, qhi, scan))
    N = P['numerical']
    frozen = {
        'S_num_Pa_inv': [N['S_num_Pa_inv'], W['legacy']['S_num_Pa_inv'],
                         lit['S_num_Pa_inv'], env['S_num_Pa_inv']],
        'mu_num_Pa_s': [N['mu_num_Pa_s'], W['legacy']['mu_num_Pa_s'],
                        lit['mu_num_Pa_s']],
        'H_n_m': [N['H_n_m'], W['legacy']['H_n_m'], lit['H_n_m']],
        't_n_s': [N['t_n_s'], W['legacy']['t_n_s'], lit['t_n_s']],
        'k_num_m2': [N['k_num_m2'], None],
    }
    for key, vals in frozen.items():
        vals = [x for x in vals if x is not None]
        if any(not close(vals[0], x, TOL_BAG) for x in vals[1:]):
            b.append('frozen constant %s is inconsistent across blocks' % key)
    if not close(N['f_b0_inv_s'] * N['Area_m2'], N['Q_total_b0_m2_per_s'], 1e-9):
        b.append('B0 f*Area != Q_total')
    dp = V.get('delta_p', {})
    if not (dp.get('low') == dp.get('high') == dp.get('rep') == 3.0e6):
        b.append('delta_p must remain a single-point ENGINEERING ASSUMPTION '
                 '(low==high==rep==3e6 Pa); it has no literature constraint')
    return b


def check_m2(P, texts, m, leg, pts):
    """Three independent formula implementations must agree; contract exact."""
    b = []
    F = P.get('formulas', {})
    for key, want in EXPECTED.items():
        if key not in F:
            b.append('formulas.%s missing' % key)
        elif not close(F[key], want, 1e-12):
            b.append('formulas.%s = %r, expected %r' % (key, F[key], want))
    if b:
        return b
    p = pts['rep']
    Ep, KIC, mu = p['Eprime'], p['KIC'], p['mu']
    mup = F['mu_prime_factor'] * mu
    q, t, CL = p['q'][0], p['t'], p['CL']
    k, S, H, dp = p['k'], p['S'], p['H'], p['dp']
    triads = [
        ('cv', ref_cv(k, mu, S), exp_cv(k, mu, S, F), None),
        ('T', ref_T(k, t, mu, S, H), exp_T(k, t, mu, S, H, F), None),
        ('CL', ref_CL(dp, k, S, mu), exp_CL(dp, k, S, mu, F), None),
        ('K', ref_K(KIC, q, Ep, mup), exp_K(KIC, q, Ep, mup, F),
         lit_K(KIC, q, Ep, mup, F)),
        ('C', ref_C(CL, Ep, t, mu, q), exp_C(CL, Ep, t, mu, q, F),
         lit_C(CL, Ep, t, mu, q, F)),
    ]
    for name, ref, ex, li in triads:
        if not close(ref, ex, TOL_MACHINE):
            b.append('%s reference vs expanded differ (rel %.3e)' % (name, rel(ref, ex)))
        if li is not None and not close(ref, li, TOL_MACHINE):
            b.append('%s reference vs literal differ (rel %.3e)' % (name, rel(ref, li)))
    # sigma_c must be independent of E under Gc = KIC^2/E'
    if not close(m['so_sigma_c_identity_check'], 1.0, 1e-12):
        b.append('sigma_c != KIC/sqrt(ell): E-prime did not cancel')
    return b


def check_m3(P, texts, m, leg, pts):
    """Reproduce every published number of first-edition tables A/B/C."""
    b = []
    for r in leg:
        pub = r['published']
        got = {
            'sigma_c_Pa': r['sigma_c'], 'T_field': r['T'],
            'CL_m_per_sqrts': r['CL'], 'K_q4': r['K'][0], 'K_q20': r['K'][1],
            'C_q4': r['C'][0], 'C_q20': r['C'][1],
            'k_num_m2': r['knum'], 'k_ratio': r['k_ratio'],
            'S_ratio': r['S_ratio'],
        }
        for key, want in pub.items():
            if key not in got:
                b.append('class %s unexpected published key %s' % (r['name'], key))
            elif not close(got[key], want, TOL_LEGACY):
                b.append('class %s %s: recomputed %.6g vs published %.6g (rel %.2f%%)'
                         % (r['name'], key, got[key], want, 100 * rel(got[key], want)))
        if not close(ref_sigma_c(r['Eprime'], r['Gc'], P['window']['legacy']['ell_m']),
                     r['sigma_c'], 1e-12):
            b.append('class %s sigma_c double-evaluation mismatch' % r['name'])
    # min/max of the recomputed columns must sit inside the published envelope
    E = P['window']['legacy_published_envelope']
    cols = {
        'T_field': [r['T'] for r in leg],
        'CL_m_per_sqrts': [r['CL'] for r in leg],
        'K_toughness': [r['K'][1] for r in leg] + [r['K'][0] for r in leg],
        'C_leakoff': [r['C'][0] for r in leg] + [r['C'][1] for r in leg],
        'k_num_m2': [r['knum'] for r in leg],
    }
    for key, vals in cols.items():
        lo, hi = min(vals), max(vals)
        if not (close(lo, E[key][0], 0.05) and close(hi, E[key][1], 0.05)):
            b.append('%s envelope [%.4g,%.4g] != published [%.4g,%.4g]'
                     % (key, lo, hi, E[key][0], E[key][1]))
    return b


def check_m4(P, texts, m, leg, pts):
    """k_core / k_eff separation and S derived from phi*(c_pore+c_fluid)."""
    b = []
    V = P['params']
    lit = P['window']['literature']
    if lit['map'].get('k_m2') != 'k_eff':
        b.append('literature permeability must map to k_eff, not %r'
                 % lit['map'].get('k_m2'))
    if m['so_kcore_over_keff_rep'] < 10.0:
        b.append('k_core / k_eff rep ratio is only %.2f (expected >=10)'
                 % m['so_kcore_over_keff_rep'])
    if not close(V['S']['rep'], m['so_S_derived_rep'], 0.20):
        b.append('S rep %.3e is not phi*(c_pore+c_fluid) %.3e'
                 % (V['S']['rep'], m['so_S_derived_rep']))
    Vk = V['k_eff']['rep']; Vkc = V['k_core']['rep']
    rep = pts['rep']
    T_core = ref_T(Vkc, rep['t'], rep['mu'], rep['S'], rep['H'])
    CL_core = ref_CL(rep['dp'], Vkc, rep['S'], rep['mu'])
    if not (T_core / rep['T'] > 3.0 and CL_core / rep['CL'] > 3.0):
        b.append('substituting k_core would NOT be detectable (T %.2fx, CL %.2fx)'
                 % (T_core / rep['T'], CL_core / rep['CL']))
    if 'k_core' in lit['map'].values():
        b.append('k_core leaked into the literature field map')
    # the recomputation must actually use k_eff
    if not close(rep['k'], Vk, 1e-12):
        b.append('literature point used k=%.3e instead of k_eff=%.3e' % (rep['k'], Vk))
    return b


def check_m5(P, texts, m, leg, pts):
    """Numerical monotonicity of every scaling vs the exponent contract."""
    b = []
    F = P['formulas']
    rep = pts['rep']
    KIC, Ep, mu, CL = rep['KIC'], rep['Eprime'], rep['mu'], rep['CL']
    mup = F['mu_prime_factor'] * mu
    q, t, k, S, H, dp = rep['q'][0], rep['t'], rep['k'], rep['S'], rep['H'], rep['dp']
    cases = [
        ('cv/k', lambda: ref_cv(2 * k, mu, S), F['T_exp_k']),
        ('cv/mu', lambda: ref_cv(k, 2 * mu, S), F['T_exp_mu']),
        ('cv/S', lambda: ref_cv(k, mu, 2 * S), F['T_exp_S']),
        ('T/t', lambda: ref_T(k, 2 * t, mu, S, H), F['T_exp_t']),
        ('T/H', lambda: ref_T(k, t, mu, S, 2 * H), F['T_exp_H']),
        ('CL/k', lambda: ref_CL(dp, 2 * k, S, mu), F['CL_exp_k']),
        ('CL/S', lambda: ref_CL(dp, k, 2 * S, mu), F['CL_exp_S']),
        ('CL/mu', lambda: ref_CL(dp, k, S, 2 * mu), F['CL_exp_mu']),
        ('K/KIC', lambda: ref_K(2 * KIC, q, Ep, mup), 1.0),
        ('K/q', lambda: ref_K(KIC, 2 * q, Ep, mup), F['K_exp_q']),
        ('K/Eprime', lambda: ref_K(KIC, q, 2 * Ep, mup), F['K_exp_Eprime']),
        ('C/CL', lambda: ref_C(2 * CL, Ep, t, mu, q), 1.0),
        ('C/q', lambda: ref_C(CL, Ep, t, mu, 2 * q), F['C_exp_q']),
        ('C/t', lambda: ref_C(CL, Ep, 2 * t, mu, q), F['C_exp_t']),
        ('C/mu', lambda: ref_C(CL, Ep, t, 2 * mu, q), F['C_exp_mu']),
        ('C/Eprime', lambda: ref_C(CL, 2 * Ep, t, mu, q), F['C_exp_Eprime']),
    ]
    baselines = [
        ref_cv(k, mu, S), ref_cv(k, mu, S), ref_cv(k, mu, S),
        ref_T(k, t, mu, S, H), ref_T(k, t, mu, S, H),
        ref_CL(dp, k, S, mu), ref_CL(dp, k, S, mu), ref_CL(dp, k, S, mu),
        ref_K(KIC, q, Ep, mup), ref_K(KIC, q, Ep, mup), ref_K(KIC, q, Ep, mup),
        ref_C(CL, Ep, t, mu, q), ref_C(CL, Ep, t, mu, q), ref_C(CL, Ep, t, mu, q),
        ref_C(CL, Ep, t, mu, q), ref_C(CL, Ep, t, mu, q),
    ]
    for (name, fn, exponent), base in zip(cases, baselines):
        new = fn()
        observed = 1 if new > base else (-1 if new < base else 0)
        expected = 1 if exponent > 0 else (-1 if exponent < 0 else 0)
        if observed != expected:
            b.append('%s: doubling gave sign %d, exponent contract says %d'
                     % (name, observed, expected))
    return b


def check_m6(P, texts, m, leg, pts):
    """Literature regime: pre-initiation storage dominance, no toughness/leakoff.

    K>3 is toughness-dominated, C~O(1) is leakoff-dominated, T~O(1) would
    mean pressure diffusion has filled the domain.  All three consistent
    parameter tiers must stay below those thresholds.  The representative
    tier gets a stricter margin so the headline conclusion is not resting on
    a boundary case.
    """
    b = []
    for tier in ('low', 'rep', 'high'):
        p = pts[tier]
        tag = 'so_%s_' % tier
        if not p['T'] < 1.0e-2:
            b.append('%s T=%.3e is not storage-dominated (need <1e-2)' % (tier, p['T']))
        if not all(x < 3.0 for x in p['K']):
            b.append('%s K=%s enters toughness-dominated regime' % (tier, p['K']))
        if not all(x < 0.1 for x in p['C']):
            b.append('%s C=%s approaches leakoff activation (>=0.1)' % (tier, p['C']))
        if not (p['K'][1] < p['K'][0] and p['C'][1] < p['C'][0]):
            b.append('%s K/C did not decrease with Q (wrong q exponents)' % tier)
        if not close(m[tag + 'K_q14'] / m[tag + 'K_q10'],
                     (10.0 / 14.0) ** 0.25, 1e-9):
            b.append('%s K ratio across Q window violates q^-1/4' % tier)
        if not close(m[tag + 'C_q14'] / m[tag + 'C_q10'],
                     (10.0 / 14.0) ** 0.5, 1e-9):
            b.append('%s C ratio across Q window violates q^-1/2' % tier)
    rep = pts['rep']
    if not (max(rep['K']) < 0.5 and max(rep['C']) < 0.05 and rep['T'] < 1.0e-4):
        b.append('representative tier lost its safety margin '
                 '(K=%s C=%s T=%.3e)' % (rep['K'], rep['C'], rep['T']))
    # The one unconstrained engineering assumption must be reported as a
    # sensitivity, never silently absorbed into the headline conclusion.
    probe = m['so_C_dp10_high_corner']
    if not (0.1 < probe < 1.0):
        b.append('dp=10 MPa sensitivity probe C=%.4f is outside the documented '
                 'boundary band (0.1,1.0); the R14 caveat is stale' % probe)
    return b


def check_m7(P, texts, m, leg, pts):
    """Similitude mapping: literature T -> (k_num,S_num), B0 comparison, q caveat."""
    b = []
    N = P['numerical']
    W = P['window']['literature']
    rep = pts['rep']
    if not close(rep['knum'], ref_knum(rep['T'], W['mu_num_Pa_s'],
                                      W['S_num_Pa_inv'], W['H_n_m'], W['t_n_s']), 1e-12):
        b.append('k_num is not the T-similitude inverse')
    # A wide T envelope is the FINDING, not a defect: it is exactly why no
    # single field rate table can be justified from these published ranges.
    # If someone silently narrows S / k / mu to fake a crisp answer, this fails.
    if m['so_T_span_decades'] < 4.0:
        b.append('literature T corner envelope spans only %.2f decades; '
                 'the published uncertainty was silently narrowed' % m['so_T_span_decades'])
    if not (m['so_T_corner_min'] < rep['T'] < m['so_T_corner_max']):
        b.append('representative T=%.3e fell outside its own corner envelope '
                 '[%.3e, %.3e]' % (rep['T'], m['so_T_corner_min'], m['so_T_corner_max']))
    if not (m['so_knum_rep_over_b0'] > 2.0):
        b.append('literature rep does not move k_num away from B0 (%.2fx)'
                 % m['so_knum_rep_over_b0'])
    # S_num is a SIMILITUDE constant, not a reservoir storativity
    if close(N['S_num_Pa_inv'], P['params']['S']['rep'], 1e-6):
        b.append('S_num was confused with the reservoir S')
    if not (m['so_q_ratio_field_over_b0'] > 1.0):
        b.append('q field/B0 ratio <= 1: injection-rate gap not visible')
    # Gate 4.15 B0 frozen values must still be recoverable from the deck text
    deck = texts['deck']
    for token in ('2.35e-18', '1.50e-10', '1.0e-03', '1.5e-5', '27000.0'):
        if token not in deck:
            b.append('B0 deck lost frozen token %s' % token)
    if not close(P['window']['legacy']['S_num_Pa_inv'] / P['window']['legacy']['S_deck_Pa_inv'],
                 0.326, TOL_LEGACY):
        b.append('S_num/S_deck != published 0.326')
    return b


def check_m8(P, texts, m, leg, pts):
    """State discipline: no field-rate claim, no production perturbation."""
    b = []
    road = texts['roadmap']
    for token in ('起裂前', '存储主导', '不可辨识', 'DFIT'):
        if token not in road:
            b.append('roadmap lost state-discipline token %r' % token)
    if 'RATE_WINDOW_DISCLAIMER' not in texts['script']:
        b.append('checker lost its explicit no-rate-optimization disclaimer')
    if 'NOT produce a field rate' not in RATE_WINDOW_DISCLAIMER:
        b.append('disclaimer text was weakened')
    # Assemble the forbidden signatures at runtime: if they appeared as
    # literal strings here, the checker would flag its own source text.
    sc = texts['script']
    forbidden = ["open(B0_DECK_PATH, " + "'w')",
                 'subprocess' + '.call',
                 'os.remove' + '(ROADMAP']
    for bad in forbidden:
        if bad in sc:
            b.append('checker contains forbidden write/solver call %r' % bad)
    build = texts['build']
    if 'check_shale_oil_params.py' not in build:
        b.append('checker is not wired into CI')
    if 'gate415_b2_verdict.log' not in build:
        b.append('B2 verdict artifact was removed from CI')
    if P['numerical'].get('beta_b0_Pa_per_s') != 1.0e5:
        b.append('B0 beta frozen constant drifted')
    return b


CHECKS = (
    ('M1', 'parameter table integrity + frozen constants', check_m1),
    ('M2', 'three-way formula identity + exponent contract', check_m2),
    ('M3', 'first-edition tables A/B/C reproduction at printed precision', check_m3),
    ('M4', 'k_core / k_eff separation + S derived from phi*c', check_m4),
    ('M5', 'scaling monotonicity vs exponent signs', check_m5),
    ('M6', 'literature regime: storage-dominated, K<3, C<<1', check_m6),
    ('M7', 'T envelope -> (k_num,S_num) similitude + B0 gap + q caveat', check_m7),
    ('M8', 'state discipline + zero production perturbation', check_m8),
)


def run_checks(payload):
    P = payload['P']
    texts = payload['texts']
    try:
        m, leg, pts = build_metrics(P)
    except (KeyError, TypeError, ValueError, ZeroDivisionError) as exc:
        return [(code, name, ['metric computation failed: %s' % exc])
                for code, name, _ in CHECKS], {}, None, None
    results = []
    for code, name, fn in CHECKS:
        try:
            issues = fn(P, texts, m, leg, pts)
        except (KeyError, TypeError, ValueError, ZeroDivisionError) as exc:
            issues = ['exception: %s' % exc]
        results.append((code, name, issues))
    return results, m, leg, pts


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def report(results, m):
    print('=== shale-oil literature parameter recomputation (M1-M8) ===')
    for k in sorted(m):
        print('%-34s = %.6e' % (k, m[k]))
    print('--- verdicts ---')
    failed = []
    for code, name, issues in results:
        if issues:
            failed.append(code)
            print('=== %s %s | FAIL' % (code, name))
            for it in issues:
                print('    FAIL %s | %s' % (code, it))
        else:
            print('=== %s %s | PASS' % (code, name))
    if failed:
        print('SHALE_OIL_PARAMS_RESULT: FAIL (%s)' % ','.join(failed))
    else:
        print('SHALE_OIL_PARAMS_RESULT: PASS')
    return 1 if failed else 0


def markdown_tables(m, leg, pts):
    L = []
    L.append('**表 A′：文献参数反算的相场材料参数（三行是「一致参数档」，'
             '即每个参数同取 low / rep / high，**不是**目标量的极值包络）**')
    L.append('')
    L.append('| 档位 | $E\'$ (GPa) | $G_c=K_{Ic}^2/E\'$ (N/m) | '
             '$\\sigma_c=K_{Ic}/\\sqrt{\\ell}$ (MPa) | $c_v$ (m²/s) | '
             '$T_{field}$ | $C_L$ (m/√s) |')
    L.append('|---|---|---|---|---|---|---|')
    for tier in ('low', 'rep', 'high'):
        L.append('| %s | %.3f | %.3f | %.3f | %.3e | %.3e | %.3e |' % (
            tier, m['so_%s_Eprime' % tier] / 1e9,
            m['so_%s_Gc' % tier], m['so_%s_sigma_c' % tier] / 1e6,
            m['so_%s_cv' % tier], m['so_%s_T' % tier],
            m['so_%s_CL' % tier]))
    L.append('')
    L.append('**表 B′：文献参数下的无量纲分区（Q = 10 → 14 m³/min）**')
    L.append('')
    L.append('| 档位 | $\\mathcal{K}$ @10 | $\\mathcal{K}$ @14 | '
             '$\\mathcal{C}$ @10 | $\\mathcal{C}$ @14 | 结论 |')
    L.append('|---|---|---|---|---|---|')
    for tier in ('low', 'rep', 'high'):
        L.append('| %s | %.4f | %.4f | %.5f | %.5f | '
                 '$\\mathcal{K}<3$、$\\mathcal{C}\\ll1$、$T\\ll1$ |' % (
                     tier, m['so_%s_K_q10' % tier], m['so_%s_K_q14' % tier],
                     m['so_%s_C_q10' % tier], m['so_%s_C_q14' % tier]))
    L.append('')
    L.append('**表 C′：由文献 $T$ 反解数值 deck 参数，并与 B0 冻结值比较**')
    L.append('')
    L.append('| 档位 | $T_{field}$ | $k_{num}$ (m²) | $k_{num}/k_{num}^{B0}$ |')
    L.append('|---|---|---|---|')
    for tier in ('low', 'rep', 'high'):
        L.append('| %s | %.3e | %.3e | %.2f |' % (
            tier, m['so_%s_T' % tier], m['so_%s_knum' % tier],
            m['so_%s_knum' % tier] / m['so_b0_knum_m2']))
    L.append('')
    L.append('**首版假设 → 文献更正版对照**')
    L.append('')
    L.append('| 量 | 首版 4.2 设定 | 文献区间 (rep) | 影响 |')
    L.append('|---|---|---|---|')
    L.append('| $\\mu$ | $1\\times10^{-3}$ Pa·s（清水） | '
             '$1\\times10^{-3}$–$6\\times10^{-2}$ ($5\\times10^{-3}$) | '
             '变黏滑溜水，$c_v$、$C_L$、$\\mathcal{K}$、$\\mathcal{C}$ 全变 |')
    L.append('| $k_{field}$ | $1.0\\times10^{-20}$–$5.0\\times10^{-18}$ m²（反推） | '
             '$k_{eff}=1\\times10^{-18}$–$1\\times10^{-17}$ ($3.2\\times10^{-18}$) m² | '
             '首版 I–III 类低 2–3 个量级 |')
    L.append('| $S$ | $4.6\\times10^{-10}$ Pa⁻¹ | '
             '$1\\times10^{-10}$–$1\\times10^{-8}$ ($1\\times10^{-9}$) Pa⁻¹，'
             '由 $\\phi(c_{pore}+c_{fluid})$ 推导 | 约 2×，且是最宽的不确定源 |')
    L.append('| $H_f$ / $t_f$ | 50 m / 3600 s | 20–100 (42) m / 3600–7200 (7200) s | '
             '$T\\propto t/H^2$，$q=Q/H_f$ |')
    L.append('| Q 扫描 | 4–20 m³/min | 10–14 (12) m³/min | 窗口窄 3.5× |')
    L.append('| $E\'$ / $K_{Ic}$ | 21.7–42 GPa / 0.6–2.0 MPa√m | '
             '8.3–64.6 (24.6) GPa / 0.2–1.2 (0.6) MPa√m | '
             '$\\sigma_c=K_{Ic}/\\sqrt{\\ell}$ 与 $E\'$ 无关 |')
    return '\n'.join(L)


# ---------------------------------------------------------------------------
# Selftest: one positive payload + mutations
# ---------------------------------------------------------------------------

def mutation(name, fn, expect=None):
    return (name, fn, expect)


def selftest():
    base = load_payload()
    pos, _, _, _ = run_checks(base)
    if any(x[2] for x in pos):
        print('=== SELFTEST FAIL: the positive payload does not pass ===')
        report(pos, build_metrics(base['P'])[0])
        return 1

    def setf(path, value):
        def _f(p):
            obj = p['P']
            for key in path[:-1]:
                obj = obj[key]
            obj[path[-1]] = value
        return _f

    M = [
        mutation('drop a source id', lambda p: p['P']['sources'].pop()),
        mutation('param rep > high', setf(['params', 'mu', 'rep'], 1.0)),
        mutation('unresolvable citation', setf(['params', 'mu', 'sources'], ['S99'])),
        mutation('weaken honesty disclaimer', setf(['honesty'], 'fine')),
        mutation('schema drift', setf(['schema'], 'shale-oil-literature-params/0')),
        mutation('wrong schema version', setf(['schema'], 'v2')),
        mutation('K_exp_q -0.25 -> -0.5', setf(['formulas', 'K_exp_q'], -0.5)),
        mutation('T_exp_H -2 -> -1', setf(['formulas', 'T_exp_H'], -1.0)),
        mutation('CL_exp_k 0.5 -> -0.5', setf(['formulas', 'CL_exp_k'], -0.5)),
        mutation('C_exp_q -0.5 -> +0.5', setf(['formulas', 'C_exp_q'], 0.5)),
        mutation('K_prefactor drift', setf(['formulas', 'K_prefactor'], 2.0)),
        mutation('inner constant 3 -> 4', setf(['formulas', 'K_inner_constant'], 4.0)),
        mutation('missing exponent key', lambda p: p['P']['formulas'].pop('C_exp_t')),
        mutation('legacy Gc drift', setf(['window', 'legacy_classes', 1, 'Gc_N_per_m'], 99.0)),
        mutation('legacy c_v drift', setf(['window', 'legacy_classes', 2, 'c_v_m2_per_s'], 1e-3)),
        mutation('published K_q20 drift',
                 setf(['window', 'legacy_classes', 0, 'published', 'K_q20'], 0.5)),
        mutation('legacy envelope mismatch',
                 setf(['window', 'legacy_published_envelope', 'T_field', 1], 9.9)),
        mutation('k_core leaked into field map', setf(['window', 'literature', 'map', 'k_m2'], 'k_core')),
        mutation('k_eff shrunk onto k_core',
                 setf(['params', 'k_eff', 'rep'], 1.2e-16)),
        mutation('S detached from phi*c', setf(['params', 'S', 'rep'], 1e-6)),
        mutation('delta_p widened into a fake range', setf(['params', 'delta_p', 'high'], 1.0e7)),
        mutation('T envelope silently narrowed', lambda p: (
            setf(['params', 'S', 'low'], 1e-9)(p),
            setf(['params', 'S', 'high'], 1e-9)(p),
            setf(['params', 'k_eff', 'low'], 3.2e-18)(p),
            setf(['params', 'k_eff', 'high'], 3.2e-18)(p))),
        mutation('dp=10MPa sensitivity probe stale', setf(['params', 'mu', 'high'], 6.0e-1)),
        mutation('Q window widened', setf(['window', 'literature', 'Q_scan_m3_min'], [4.0, 20.0])),
        mutation('KIC pushed into toughness regime', setf(['params', 'KIC', 'rep'], 100.0)),
        mutation('mu pushed into leakoff regime', setf(['params', 'mu', 'rep'], 1e-12)),
        mutation('t_f pushed out of storage regime', setf(['params', 't_f', 'rep'], 1e9)),
        mutation('frozen S_num drift', setf(['numerical', 'S_num_Pa_inv'], 2.0e-10)),
        mutation('frozen k_num drift', setf(['numerical', 'k_num_m2'], 1.0e-10)),
        mutation('B0 beta drift', setf(['numerical', 'beta_b0_Pa_per_s'], 2.0e5)),
        mutation('deck lost perm_min',
                 lambda p: p['texts'].__setitem__('deck', 'storativity: 1.50e-10\n')),
        mutation('roadmap lost state discipline',
                 lambda p: p['texts'].__setitem__('roadmap', 'no discipline here')),
        mutation('B2 artifact removed from CI',
                 lambda p: p['texts'].__setitem__('build', 'check_shale_oil_params.py\n')),
    ]

    failed = []
    for name, fn, _ in M:
        p = copy.deepcopy(base)
        fn(p)
        res, _, _, _ = run_checks(p)
        codes = [c for c, _, issues in res if issues]
        if not codes:
            failed.append((name, 'mutation was NOT detected'))
        print('  mutation %-42s -> %s' % (name, ','.join(codes) or 'UNDETECTED'))

    # invariance: edits outside the contract must stay green
    invariance = [
        mutation('as_of / citation text only', lambda p: (
            p['P'].__setitem__('as_of', '2099-01-01'),
            p['P']['sources'][0].__setitem__('cite', 'typo-fixed citation'))),
        mutation('extra unrelated source note',
                 lambda p: p['P']['sources'][0].__setitem__('extra_note', 'harmless')),
        mutation('append-only new param',
                 lambda p: p['P']['params'].__setitem__('new_metric', {
                     'unit': '-', 'low': 0.0, 'high': 1.0, 'rep': 0.5,
                     'sources': [], 'note': 'append-only'})),
    ]
    for name, fn, _ in invariance:
        p = copy.deepcopy(base)
        fn(p)
        res, _, _, _ = run_checks(p)
        codes = [c for c, _, issues in res if issues]
        if codes:
            failed.append((name, 'invariance broken: %s' % ','.join(codes)))
        print('  invariance %-40s -> %s' % (name, ','.join(codes) or 'still PASS'))

    if failed:
        print('=== SELFTEST FAIL ===')
        for name, why in failed:
            print('    %s: %s' % (name, why))
        return 1
    print('=== SELFTEST PASS: 1 positive + %d mutations (all rejected) + %d invariances '
          '(all preserved) ===' % (len(M), len(invariance)))
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--md', action='store_true')
    ap.add_argument('--json', dest='json_out', default=None)
    args = ap.parse_args()
    if args.selftest:
        return selftest()
    payload = load_payload()
    results, m, leg, pts = run_checks(payload)
    if args.json_out:
        with open(args.json_out, 'w', encoding='utf-8') as fh:
            json.dump(m, fh, indent=2, sort_keys=True)
    if args.md:
        print(markdown_tables(m, leg, pts))
    rc = report(results, m)
    return rc


if __name__ == '__main__':
    sys.exit(main())
