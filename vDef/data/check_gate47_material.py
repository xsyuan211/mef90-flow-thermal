#!/usr/bin/env python3
"""Gate 4.7 -- TIME-DEPENDENT MATERIAL PARAMETERS (Darcy A=k/mu and S).

Pure standard library (the CI runner has no numpy / scipy / PyYAML).
Run with ``--gate`` for the hard verdict (``GATE47_RESULT: PASS|FAIL``, exit
0/1), ``--selftest`` for the offline positive + negative replay, or with four
paths for a report-only run.

WHAT THE GATE FIXES / PROVES
----------------------------
Before Gate 4.7 ``MEF90FlowSetTransients`` scaled ONLY the volumetric source
and the pressure BC. The two material fields that actually enter the Darcy
operator -- the conductivity ``A = k(v)/mu`` (``FlowCtx%permField``) and the
storativity ``S`` (``FlowCtx%storativity``) -- had no time entry point: A was
rebuilt from damage only, and S was assembled ONCE at init. Compaction,
pressure-dependent permeability and fluid compressibility evolution could not
be represented.

Gate 4.7 adds six option-bag fields (perm_scaling / perm_rate /
stor_scaling / stor_rate / material_perm_offset / material_stor_offset), two
persistent per-cell factor Vecs (permTimeFactor / storTimeFactor) and three
laws for each channel:

    constant : factor = 1                  (strict no-op)
    linear   : factor = 1 + rate*t        (t=0 -> 1, no singularity)
    file     : factor = EXO cell field at material_*_offset / step

``SetTransients`` fills the factor Vecs, and ``UpdatePerm`` /
``UpdateStorativity`` multiply the freshly REBUILT base fields by them.
``vDef.F90`` calls ``SetTransients`` BEFORE both rebuilds and rebuilds S on
EVERY outer step, so there is no one-step lag and no cumulative VecScale
drift. Because mu enters the assembled operator only through A = k/mu, a
mu(t) law is expressed by the same A channel; Gate 4.6's spatial mu(T) field
is untouched.

FALSIFIABLE EXPERIMENT
----------------------
``terzaghi1D_material.yaml`` is ``terzaghi1D.yaml`` plus exactly FOUR lines:

    perm_scaling: linear   perm_rate: 0.8
    stor_scaling: linear   stor_rate: 0.4

The rates are deliberately OPPOSED: A(t) rises (faster drainage) and S(t)
rises too (more storage, slower drainage). With fixed-dt backward Euler the
modal amplitudes are

    a_m(n) = a_m(0) * PROD_{j=1..n} (1 + c_v(j*dt) lam_m^2 dt)^-1
    c_v(t) = c_v0 * (1 + 0.8 t) / (1 + 0.4 t)
    lam_m  = (2m+1) pi / (2H),  H = 10 m, c_v0 = (5e-11/1e-3)/4.6e-10

so the final closed-end field is ~0.664 of the constant-c_v baseline. A bug
that wires only ONE channel gives a different field (relative L2 0.85 / 0.56
against the full solution at step 50), and a completely dead entry point
gives the production field (0.336). X7 separates all four hypotheses.

CHECKS
------
  X1  run integrity: 51 records, identical geometry/time axis, no crash
  X2  static anti-regression: CST/Linear/File branches, factor Vecs read by
      both coefficient assemblers, bag declaration order == registration
      order, constructor arity, vDef call order
  X3  deck isolation: exactly the four expected added lines
  X4  bag visibility: linear + rates in the material log, constant + 0 in the
      production log
  X5  EXO permField (slot 4) == 5e-8*(1+0.8 t) on EVERY record; control flat
  X6  51 ``Gate 4.7:`` rows read back from the factor Vecs match the law
  X7  pressure field vs the time-varying BE product series, plus separation
      from no-scaling / single-channel hypotheses
  X8  production isolation: control field still matches the constant-c_v
      series (Gate 3.3 P6), control permField and factors flat
"""

import difflib
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

# --- Physics / discretisation (terzaghi1D.yaml, CI #66) -------------------
P0 = 1.0e5
K_MIN = 5.0e-11
MU = 1.0e-3
A0 = K_MIN / MU                 # 5e-8: base Darcy conductivity k/mu
STOR0 = 4.6e-10
CV0 = A0 / STOR0                # 108.6957 m^2/s
H = 10.0                        # drainage path, x in [-5,5], z = 5 - x
DT = 0.02                       # fixed backward-Euler step
NSERIES = 300
MIN_RECORDS = 50

PERM_RATE = 0.8                 # 1/s, terzaghi1D_material.yaml
STOR_RATE = 0.4                 # 1/s, terzaghi1D_material.yaml

# --- Hard thresholds (calibrated offline against the CI #109 real field) --
L2_FULL_MAX = 0.08              # FE vs time-varying BE series (predicted ~0.03-0.05)
COS_MIN = 0.99
NORM_LO, NORM_HI = 0.88, 1.12
SEP_CONST_MIN = 0.20           # vs no scaling (analytic 0.336)
SEP_ONLYK_MIN = 0.40           # vs perm-only bug (analytic 0.850)
SEP_ONLYS_MIN = 0.30           # vs stor-only bug (analytic 0.561)
CTRL_L2_MAX = 0.05             # Gate 3.3 P6 = 0.03206
PERM_REL_TOL = 2.0e-6          # on the record-0-normalised A(t) law
PERM_BASE_TOL = 2.0e-3         # loose absolute sanity on k(v)/mu at t=0
FACTOR_TOL = 1.0e-8
GEOM_TOL = 1.0e-9

DECK_CTRL = os.path.join(_HERE, 'terzaghi1D.yaml')
DECK_MAT = os.path.join(_HERE, 'terzaghi1D_material.yaml')
SRC_ASM = os.path.join(_HERE, '..', '..', 'm_Flow',
                       'm_MEF90_FlowAssembly.F90')
SRC_FLOW = os.path.join(_HERE, '..', '..', 'm_Flow', 'm_MEF90_Flow.F90')
SRC_CTX = os.path.join(_HERE, '..', '..', 'm_Flow', 'm_MEF90_FlowCtx.F90')
SRC_DEFAULT = os.path.join(_HERE, '..', 'm_vDefDefault.F90')
SRC_VDEF = os.path.join(_HERE, '..', 'vDef.F90')

FATAL_RE = re.compile(
    r'Caught signal|SEGV|MPI_ABORT|PETSC ERROR|Unknown option|'
    r'attempt to deallocate|Error opening EXO file|TSSolve failed|'
    r'fortran runtime error', re.I)
G47_RE = re.compile(
    r'Gate\s+4\.7:\s*step\s+(\d+)\s+t\s+([-+0-9.eEdD]+)\s+'
    r'permF\s+([-+0-9.eEdD]+)\s+storF\s+([-+0-9.eEdD]+)', re.I)
BAG_HEADER_RE = re.compile(r'^\s*PetscBag Object:\s*(.*?)\s*$')
BAG_ITEM_RE = re.compile(r'^\s{2}(\w+)\s*=\s*(.*?)\s*;')
BAG_KEYS = ('perm_scaling', 'perm_rate', 'stor_scaling', 'stor_rate',
            'material_perm_offset', 'material_stor_offset')
EXPECTED_DECK_LINES = (
    'perm_scaling: linear', 'perm_rate: 0.8',
    'stor_scaling: linear', 'stor_rate: 0.4',
)

_LAM = [(2.0 * m + 1.0) * math.pi / (2.0 * H) for m in range(NSERIES)]
_AMP = [(4.0 / math.pi) * ((-1.0) ** m) / (2.0 * m + 1.0)
        for m in range(NSERIES)]


# --------------------------------------------------------------------------
# Analytic backward-Euler product series with time-dependent c_v
# --------------------------------------------------------------------------
def decay_products(rk, rs):
    """[n][m] modal amplitude multiplier after n BE steps."""
    dec = [1.0] * NSERIES
    out = [list(dec)]
    for n in range(1, MIN_RECORDS + 1):
        t = n * DT
        cv = CV0 * (1.0 + rk * t) / (1.0 + rs * t)
        dec = [d / (1.0 + cv * l * l * DT)
               for d, l in zip(dec, _LAM)]
        out.append(list(dec))
    return out


_DEC = {
    'full': decay_products(PERM_RATE, STOR_RATE),
    'const': decay_products(0.0, 0.0),
    'onlyk': decay_products(PERM_RATE, 0.0),
    'onlys': decay_products(0.0, STOR_RATE),
}
_COS_CACHE = {}


def cos_table(x):
    key = tuple(x)
    if key not in _COS_CACHE:
        z = [5.0 - float(xi) for xi in x]
        tab = [[math.cos(l * zi) for l in _LAM] for zi in z]
        _COS_CACHE[key] = (z, tab)
    return _COS_CACHE[key]


def be_field(tag, n, x):
    """Terzaghi field after n BE steps for hypothesis ``tag``."""
    z, ctab = cos_table(x)
    dec = _DEC[tag][min(n, MIN_RECORDS)]
    out = []
    for zi, ct in zip(z, ctab):
        if zi >= H - 1e-12:
            out.append(0.0)
        else:
            out.append(P0 * sum(a * c * d for a, c, d
                                in zip(_AMP, ct, dec)))
    return out


# --------------------------------------------------------------------------
# Small numeric helpers
# --------------------------------------------------------------------------
def _num(s):
    return float(s.replace('D', 'E').replace('d', 'e'))


def _nrm(a):
    return math.sqrt(sum(v * v for v in a))


def _dot(a, b):
    return sum(u * v for u, v in zip(a, b))


def _cos(a, b):
    d = _nrm(a) * _nrm(b)
    return _dot(a, b) / d if d > 0 else float('nan')


def _rel_l2(a, b):
    nb = _nrm(b)
    if nb <= 0:
        return float('inf')
    return math.sqrt(sum((u - v) ** 2 for u, v in zip(a, b))) / nb


def read_text(path):
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as fh:
            return fh.read()
    except Exception:                                  # noqa: BLE001
        return ''


def strip_yaml(text):
    out = []
    for line in text.splitlines():
        body = re.sub(r'#.*$', '', line).rstrip()
        if body.strip():
            out.append(body)
    return out


def sub_body(src, name):
    pat = re.compile(
        r'^\s*Subroutine\s+' + re.escape(name) + r'\b.*?'
        r'^\s*End\s+Subroutine\s+' + re.escape(name) + r'\s*$',
        re.I | re.S | re.M)
    m = pat.search(src)
    return m.group(0) if m else ''


# --------------------------------------------------------------------------
# IO
# --------------------------------------------------------------------------
def load_exo(path):
    if CDF2 is None:
        raise ValueError('cannot import parse_cdf2: %s' % _CDF2_IMPORT_ERROR)
    if not path or not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise ValueError('missing/empty mesh: %s' % path)
    c = CDF2(path)
    x = [float(v) for v in c.read_static('coordx')]
    nn = len(x)
    if nn <= 0 or not all(math.isfinite(v) for v in x):
        raise ValueError('bad coordx in: %s' % path)
    times, p, perm = [], [], []
    for i in range(c.numrecs):
        times.append(float(c.read_record('time_whole', i)[0]))
        p.append([float(v) for v in c.read_record('vals_nod_var1', i)])
        perm.append([float(v) for v in c.read_record('vals_elem_var4eb1', i)])
    if not all(math.isfinite(t) for t in times):
        raise ValueError('non-finite time axis in: %s' % path)
    for row in p + perm:
        if not row or not all(math.isfinite(v) for v in row):
            raise ValueError('non-finite / empty field in: %s' % path)
    for row in p:
        if len(row) != nn:
            raise ValueError('nodal field length mismatch in: %s' % path)
    return dict(path=path, nn=nn, nrec=c.numrecs, x=x,
                times=times, p=p, perm=perm)


def parse_bag_global(text):
    """Return the Flow Global Options bag as {key: value}."""
    glob, current = {}, None
    for raw in text.splitlines():
        mh = BAG_HEADER_RE.match(raw)
        if mh:
            current = glob if 'flow Flow Global Options' in mh.group(1) else None
            continue
        if current is None:
            continue
        mi = BAG_ITEM_RE.match(raw)
        if mi:
            current[mi.group(1)] = mi.group(2).strip()
    return glob


def parse_factor_rows(text):
    rows = []
    for m in G47_RE.finditer(text):
        try:
            rows.append((int(m.group(1)), _num(m.group(2)),
                         _num(m.group(3)), _num(m.group(4))))
        except ValueError:
            continue
    return rows


# --------------------------------------------------------------------------
# X1 - run integrity
# --------------------------------------------------------------------------
def check_x1(mat, ctrl, mat_log, ctrl_log):
    issues = []
    if not mat_log or not ctrl_log:
        issues.append('a solver log is missing/empty')
    for tag, text in (('material', mat_log), ('control', ctrl_log)):
        if text and FATAL_RE.search(text):
            issues.append('%s log has a fatal signature' % tag)
    for tag, d in (('material', mat), ('control', ctrl)):
        if d['nrec'] < MIN_RECORDS + 1:
            issues.append('%s has only %d records' % (tag, d['nrec']))
    if mat['nn'] != ctrl['nn'] or any(
            abs(a - b) > GEOM_TOL for a, b in zip(mat['x'], ctrl['x'])):
        issues.append('geometry differs between the two meshes')
    if len(mat['times']) != len(ctrl['times']) or any(
            abs(a - b) > 1e-9 for a, b in zip(mat['times'], ctrl['times'])):
        issues.append('time axes differ')
    return (not issues, 'nrec=%d/%d nn=%d fatal=%s; issues=%s'
            % (mat['nrec'], ctrl['nrec'], mat['nn'],
               bool(FATAL_RE.search(mat_log + ctrl_log)),
               ';'.join(issues) or 'none'))


# --------------------------------------------------------------------------
# X2 - static anti-regression
# --------------------------------------------------------------------------
def check_x2(sources):
    asm, flow, ctx, dflt, vdef = sources
    issues = []
    st = sub_body(asm, 'MEF90FlowSetTransients')
    if not st:
        issues.append('SetTransients missing')
    else:
        for token in ('gOpt%permeabilityScaling', 'gOpt%storativityScaling',
                      'MEF90Scaling_Linear', 'MEF90Scaling_File',
                      'VecLoadExodusCell', 'FlowCtx%permTimeFactor',
                      'FlowCtx%storTimeFactor', "Gate 4.7:"):
            if token not in st:
                issues.append('SetTransients lacks ' + token)
    up = sub_body(flow, 'MEF90FlowUpdatePerm')
    us = sub_body(flow, 'MEF90FlowUpdateStorativity')
    if 'FlowCtx%permTimeFactor' not in up or 'facPtr(1)' not in up:
        issues.append('UpdatePerm does not consume permTimeFactor')
    if 'FlowCtx%storTimeFactor' not in us or 'facPtr(1)' not in us:
        issues.append('UpdateStorativity does not consume storTimeFactor')
    # Gate 4.6 V1: every coefficient routine must stay alias-free.
    for name, body in (('UpdatePerm', up), ('UpdateStorativity', us)):
        if 'SectionRealDuplicate(FlowCtx%cellDMScalSec' not in body:
            issues.append(name + ': no private section')
        if re.search(r'SectionReal(?:Restrict|Restore)\s*\(\s*'
                     r'FlowCtx%cellDMScalSec', body, re.I):
            issues.append(name + ': direct shared-section access')
    # PetscBag declaration order must equal registration order (CI #50).
    tb = re.search(r'^\s*Type MEF90FlowGlobalOptions_Type\b(.*?)'
                   r'^\s*End Type MEF90FlowGlobalOptions_Type',
                   ctx, re.S | re.M)
    reg = re.search(r'Subroutine PetscBagRegisterMEF90FlowCtxGlobalOptions'
                    r'\b(.*?)End Subroutine', ctx, re.S)
    if not tb or not reg:
        issues.append('cannot parse the option bag')
    else:
        decl = re.findall(r'^\s*(?:PetscEnum|PetscBool|PetscReal|PetscInt)'
                          r'\s*::\s*(\w+)', tb.group(1), re.M)
        regf = re.findall(r'PetscBagRegister(?:Enum|Bool|Real|Int)\s*\(\s*'
                          r'bag,\s*FlowGlobalOptions%(\w+)', reg.group(1))
        if decl != regf:
            issues.append('bag declaration/registration order mismatch')
        if len(decl) != 36:
            issues.append('expected 36 option fields, got %d' % len(decl))
        for key in BAG_KEYS:
            if "'%s'" % key not in ctx:
                issues.append('option %s not registered' % key)
    # Default constructor must supply all 36 values.
    if not re.search(r'36 components', dflt):
        issues.append('default constructor comment not updated')
    if dflt.count('MEF90Scaling_CST, &') < 4:
        issues.append('default constructor missing CST defaults')
    # vDef: SetTransients BEFORE UpdatePerm, plus a per-step S rebuild.
    m = re.search(r"Gate 3\.1: step',step,' enter transient block", vdef)
    seg = vdef[m.start():m.start() + 6500] if m else ''
    calls = re.findall(r'Call (MEF90Flow\w+)\(', seg)
    for name in ('MEF90FlowSetTransients', 'MEF90FlowUpdatePerm',
                 'MEF90FlowUpdateStorativity'):
        if name not in calls:
            issues.append('transient loop lacks ' + name)
    if ('MEF90FlowSetTransients' in calls and 'MEF90FlowUpdatePerm' in calls
            and calls.index('MEF90FlowSetTransients')
            > calls.index('MEF90FlowUpdatePerm')):
        issues.append('SetTransients still runs after UpdatePerm (one-step lag)')
    return (not issues, 'static anti-regression; issues=%s'
            % (';'.join(issues) or 'none'))


# --------------------------------------------------------------------------
# X3 - deck isolation
# --------------------------------------------------------------------------
def check_x3(mat_deck, ctrl_deck):
    a, b = strip_yaml(ctrl_deck), strip_yaml(mat_deck)
    diff = [l for l in difflib.unified_diff(a, b, lineterm='')
            if l.startswith(('+', '-')) and not l.startswith(('+++', '---'))]
    ok = (len(diff) == 4 and all(l.startswith('+') for l in diff)
          and [l[1:].strip() for l in diff] == list(EXPECTED_DECK_LINES))
    return (ok, '%d stripped lines, %d diffs: %s'
            % (len(b), len(diff), ' | '.join(d.strip() for d in diff[:6])))


# --------------------------------------------------------------------------
# X4 - bag visibility
# --------------------------------------------------------------------------
def _bag_value(bag, key):
    return bag.get(key, '').strip()


def check_x4(mat_log, ctrl_log):
    mb, cb = parse_bag_global(mat_log), parse_bag_global(ctrl_log)
    issues = []
    for key in ('perm_scaling', 'perm_rate', 'stor_scaling', 'stor_rate'):
        if key not in mb:
            issues.append('material bag lacks ' + key)
        if key not in cb:
            issues.append('control bag lacks ' + key)
    want_m = {'perm_scaling': 'linear', 'stor_scaling': 'linear'}
    want_c = {'perm_scaling': 'constant', 'stor_scaling': 'constant'}
    for key, want in want_m.items():
        if key in mb and _bag_value(mb, key).lower() != want:
            issues.append('material %s=%s (want %s)'
                          % (key, mb[key], want))
    for key, want in want_c.items():
        if key in cb and _bag_value(cb, key).lower() != want:
            issues.append('control %s=%s (want %s)'
                          % (key, cb[key], want))
    for key, want in (('perm_rate', PERM_RATE), ('stor_rate', STOR_RATE)):
        try:
            got = float(_bag_value(mb, key))
            if abs(got - want) > 1e-9:
                issues.append('material %s=%g (want %g)' % (key, got, want))
        except (ValueError, KeyError):
            issues.append('material %s not numeric' % key)
        try:
            got = float(_bag_value(cb, key))
            if abs(got) > 1e-12:
                issues.append('control %s=%g (want 0)' % (key, got))
        except (ValueError, KeyError):
            issues.append('control %s not numeric' % key)
    return (not issues, 'material linear %.1f/%.1f, control constant 0/0; issues=%s'
            % (PERM_RATE, STOR_RATE, ';'.join(issues) or 'none'))


# --------------------------------------------------------------------------
# X5 - EXO permField factor
# --------------------------------------------------------------------------
def check_x5(mat, ctrl):
    """A(t) is checked RELATIVE to each run's own record 0.

    Normalising by the run's own t=0 field cancels any small damage-dependent
    spatial pattern in k(v), so the tight 2e-6 tolerance tests ONLY the time
    law and cannot be confused by a slightly different mechanical field.
    A loose absolute check keeps the base value honest.
    """
    issues, worst = [], 0.0
    mbase, cbase = mat['perm'][0], ctrl['perm'][0]
    if len(mbase) != len(cbase):
        issues.append('cell counts differ between the two meshes')
    base_rel = max(abs(v / A0 - 1.0) for v in mbase)
    if base_rel > PERM_BASE_TOL:
        issues.append('material base permField rel=%.3g (want ~5e-8)'
                      % base_rel)
    # Control must stay flat over ALL 51 records (CST law is a no-op).
    ctrl_worst = 0.0
    for row in ctrl['perm']:
        ctrl_worst = max(ctrl_worst,
                         max(abs(v / c - 1.0) for v, c in zip(row, cbase)))
    if ctrl_worst > PERM_REL_TOL:
        issues.append('control permField drifts (rel=%.3g)' % ctrl_worst)
    # Material must follow exactly (1 + 0.8 t) cell by cell.
    for n, t in enumerate(mat['times']):
        fac = 1.0 + PERM_RATE * t
        for c, v in enumerate(mat['perm'][n]):
            want = mbase[c] * fac
            rel = abs(v / want - 1.0) if want else abs(v - want)
            worst = max(worst, rel)
        if worst > PERM_REL_TOL:
            issues.append('material record %d rel=%.3g' % (n, worst))
            break
    ratio = (mat['perm'][-1][0] / mbase[0])
    if abs(ratio - (1.0 + PERM_RATE * mat['times'][-1])) > 1e-6:
        issues.append('final/initial perm ratio=%.6f (want %.6f)'
                      % (ratio, 1.0 + PERM_RATE * mat['times'][-1]))
    return (not issues,
            'max rel|A(t)/A(0)(1+0.8t)-1|=%.3g, base rel=%.3g, '
            'control drift=%.3g, ratio(final/initial)=%.6f; issues=%s'
            % (worst, base_rel, ctrl_worst, ratio,
               ';'.join(issues) or 'none'))


# --------------------------------------------------------------------------
# X6 - factor rows read back from the Vecs
# --------------------------------------------------------------------------
def _check_rows(rows, rk, rs, tag):
    if len(rows) < MIN_RECORDS + 1:
        return ['%s has only %d Gate 4.7 rows' % (tag, len(rows))]
    issues = []
    for step, t, pf, sf in rows:
        if step < 1 or step > MIN_RECORDS + 1:
            issues.append('%s step %d out of range' % (tag, step))
            continue
        want_t = (step - 1) * DT
        if abs(t - want_t) > 1e-9:
            issues.append('%s step %d t=%.6g (want %.6g)'
                          % (tag, step, t, want_t))
        if abs(pf - (1.0 + rk * want_t)) > FACTOR_TOL:
            issues.append('%s step %d permF=%.8g' % (tag, step, pf))
        if abs(sf - (1.0 + rs * want_t)) > FACTOR_TOL:
            issues.append('%s step %d storF=%.8g' % (tag, step, sf))
    return issues[:6]


def check_x6(mat_log, ctrl_log):
    mr, cr = parse_factor_rows(mat_log), parse_factor_rows(ctrl_log)
    issues = (_check_rows(mr, PERM_RATE, STOR_RATE, 'material')
              + _check_rows(cr, 0.0, 0.0, 'control'))
    return (not issues,
            'rows material=%d control=%d, permF(51)=%.4g storF(51)=%.4g; issues=%s'
            % (len(mr), len(cr),
               mr[-1][2] if mr else float('nan'),
               mr[-1][3] if mr else float('nan'),
               ';'.join(issues) or 'none'))


# --------------------------------------------------------------------------
# X7 - analytic time-varying BE series + channel separation
# --------------------------------------------------------------------------
def check_x7(mat):
    x = mat['x']
    issues, met = [], {}
    norms = [_nrm(row) for row in mat['p']]
    if not all(norms[i] >= norms[i + 1] - 1e-6 * norms[0]
               for i in range(len(norms) - 1)):
        issues.append('pressure norm is not monotonically decaying')
    for n in (10, 25, 50):
        if n >= mat['nrec']:
            issues.append('record %d missing' % n)
            continue
        full = be_field('full', n, x)
        got = mat['p'][n]
        l2 = _rel_l2(got, full)
        cs = _cos(got, full)
        ratio = _nrm(got) / _nrm(full)
        met['l2_%d' % n] = l2
        met['cos_%d' % n] = cs
        met['ratio_%d' % n] = ratio
        if l2 > L2_FULL_MAX:
            issues.append('n=%d relL2=%.4f > %.3f' % (n, l2, L2_FULL_MAX))
        if not (cs >= COS_MIN):
            issues.append('n=%d cos=%.5f < %.3f' % (n, cs, COS_MIN))
        if not (NORM_LO <= ratio <= NORM_HI):
            issues.append('n=%d norm ratio=%.4f outside [%g,%g]'
                          % (n, ratio, NORM_LO, NORM_HI))
    if 50 < mat['nrec']:
        got = mat['p'][50]
        for tag, lim in (('const', SEP_CONST_MIN),
                         ('onlyk', SEP_ONLYK_MIN),
                         ('onlys', SEP_ONLYS_MIN)):
            sep = _rel_l2(got, be_field(tag, 50, x))
            met['sep_' + tag] = sep
            if sep < lim:
                issues.append('n=50 cannot reject %s (relL2=%.4f < %.2f)'
                              % (tag, sep, lim))
    return (not issues,
            ('relL2(50)=%.4f cos=%.5f ratio=%.4f | separation '
             'const=%.3f onlyk=%.3f onlys=%.3f; issues=%s')
            % (met.get('l2_50', float('nan')),
               met.get('cos_50', float('nan')),
               met.get('ratio_50', float('nan')),
               met.get('sep_const', float('nan')),
               met.get('sep_onlyk', float('nan')),
               met.get('sep_onlys', float('nan')),
               ';'.join(issues) or 'none')), met


# --------------------------------------------------------------------------
# X8 - production isolation
# --------------------------------------------------------------------------
def check_x8(ctrl):
    x = ctrl['x']
    issues, worst = [], 0.0
    for n in (10, 25, 50):
        if n >= ctrl['nrec']:
            issues.append('control record %d missing' % n)
            continue
        l2 = _rel_l2(ctrl['p'][n], be_field('const', n, x))
        worst = max(worst, l2)
        if l2 > CTRL_L2_MAX:
            issues.append('control n=%d relL2=%.4f > %.3f'
                          % (n, l2, CTRL_L2_MAX))
    cbase = ctrl['perm'][0]
    flat = max(abs(v / c - 1.0)
               for row in ctrl['perm'] for v, c in zip(row, cbase))
    if flat > PERM_REL_TOL:
        issues.append('control permField not flat (rel=%.3g)' % flat)
    return (not issues,
            'control max relL2 vs constant-c_v BE=%.4f, perm flat rel=%.3g; issues=%s'
            % (worst, flat, ';'.join(issues) or 'none'))


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------
def evaluate(sources, decks, logs, meshes):
    chk, met = {}, {}
    try:
        mat, ctrl = meshes
        mat_log, ctrl_log = logs
        mat_deck, ctrl_deck = decks
        chk['X1'] = check_x1(mat, ctrl, mat_log, ctrl_log)
        chk['X2'] = check_x2(sources)
        chk['X3'] = check_x3(mat_deck, ctrl_deck)
        chk['X4'] = check_x4(mat_log, ctrl_log)
        chk['X5'] = check_x5(mat, ctrl)
        chk['X6'] = check_x6(mat_log, ctrl_log)
        x7, met = check_x7(mat)
        chk['X7'] = x7
        chk['X8'] = check_x8(ctrl)
    except Exception as exc:                           # noqa: BLE001
        return ({'X0': (False, 'cannot build Gate 4.7 inputs: %r' % (exc,))},
                {})
    return chk, met


def print_report(chk, met):
    print('=== Gate 4.7 -- time-dependent material parameters verdict ===')
    labels = {
        'X1': 'run integrity (51 records, same geometry, no crash)',
        'X2': 'static anti-regression (branches, factor Vecs, bag order, call order)',
        'X3': 'deck isolation (exactly four added lines)',
        'X4': 'bag visibility (linear + rates, control constant)',
        'X5': 'EXO permField == 5e-8*(1+0.8t) on every record',
        'X6': '51 Gate 4.7 rows read back from factor Vecs',
        'X7': 'time-varying BE product series + channel separation',
        'X8': 'production isolation (constant-c_v series unchanged)',
    }
    fails = []
    for key in sorted(chk):
        ok, msg = chk[key]
        print('  %s %-3s %-70s | %s'
              % (key, 'PASS' if ok else 'FAIL', labels.get(key, ''), msg))
        if not ok:
            fails.append(key)
    if met:
        print('  --- metrics ---')
        for key in ('l2_10', 'l2_25', 'l2_50', 'cos_50', 'ratio_50',
                    'sep_const', 'sep_onlyk', 'sep_onlys'):
            if key in met:
                print('  %-12s = %.6g' % (key, met[key]))
    return fails


def main(argv):
    if argv and argv[0] == '--selftest':
        return selftest()
    gate = bool(argv) and argv[0] == '--gate'
    args = argv[1:] if gate else argv
    usage = '<material.log> <control.log> <material.gen> <control.gen>'
    if len(args) != 4:
        if gate:
            print('GATE47_RESULT: FAIL (usage: --gate %s)' % usage)
            return 1
        print('usage: check_gate47_material.py [--gate|--selftest] %s' % usage)
        return 0
    mat_log_p, ctrl_log_p, mat_gen_p, ctrl_gen_p = args
    sources = tuple(read_text(p) for p in
                    (SRC_ASM, SRC_FLOW, SRC_CTX, SRC_DEFAULT, SRC_VDEF))
    decks = (read_text(DECK_MAT), read_text(DECK_CTRL))
    logs = (read_text(mat_log_p), read_text(ctrl_log_p))
    err = None
    meshes = None
    try:
        meshes = (load_exo(mat_gen_p), load_exo(ctrl_gen_p))
        if not all(logs) or not all(decks) or not all(sources):
            raise ValueError('a log, deck or source file is unreadable')
    except Exception as exc:                           # noqa: BLE001
        err = 'cannot build Gate 4.7 inputs: %r' % (exc,)
    if err:
        chk, met = {'X0': (False, err)}, {}
    else:
        chk, met = evaluate(sources, decks, logs, meshes)
    fails = print_report(chk, met)
    if gate:
        if err or fails or len(chk) < 8:
            print('GATE47_RESULT: FAIL (%s)'
                  % (err if err else ','.join(fails) or 'incomplete criteria'))
            return 1
        print('GATE47_RESULT: PASS (X1-X8: perm_scaling=linear 0.8, '
              'stor_scaling=linear 0.4; EXO A(t), factor rows, '
              'time-varying BE series relL2(50)=%.4f, '
              'onlyk/onlys rejected, production unchanged)'
              % met.get('l2_50', float('nan')))
        return 0
    print('VERDICT: %s (report mode, informational only)'
          % ('FAIL ' + ','.join(fails) if fails else 'PASS'))
    return 0


# --------------------------------------------------------------------------
# Offline selftest: analytic replay + falsification cases
# --------------------------------------------------------------------------
def _synthetic_log(rk, rs, bag_mat, nrec=51):
    """Fortran-list-directed-style Gate 4.7 rows + a Flow Global bag view."""
    lines = ['PetscBag Object:  flow Flow Global Options']
    if bag_mat:
        lines += ['  perm_scaling = linear; (MEF90scaling) law',
                  '  perm_rate = 0.8 ; rate',
                  '  stor_scaling = linear; (MEF90scaling) law',
                  '  stor_rate = 0.4 ; rate',
                  '  material_perm_offset = 0 ; slot',
                  '  material_stor_offset = 0 ; slot']
    else:
        lines += ['  perm_scaling = constant; (MEF90scaling) law',
                  '  perm_rate = 0 ; rate',
                  '  stor_scaling = constant; (MEF90scaling) law',
                  '  stor_rate = 0 ; rate',
                  '  material_perm_offset = 0 ; slot',
                  '  material_stor_offset = 0 ; slot']
    lines.append('PetscBag Object:  cs1 Flow CellSet Options')
    for step in range(1, nrec + 1):
        t = (step - 1) * DT
        lines.append(' Gate 4.7: step %6d  t %20.16f  permF %20.16f'
                     '  storF %20.16f  '
                     % (step, t, 1.0 + rk * t, 1.0 + rs * t))
    return '\n'.join(lines) + '\n'


def _synthetic_mesh(tag, x, times, fe_error=None, scale=1.0,
                    perm_rk=PERM_RATE):
    """Synthetic field = analytic(tag) + scaled real-FE error (or zero).

    ``perm_rk`` is deliberately independent of ``tag``: the EXO permField
    slot is controlled by the permeability channel alone, so a single-channel
    pressure bug can still have a perfectly live A(t).
    """
    p, perm = [], []
    for n, t in enumerate(times):
        row = be_field(tag, n, x)
        if fe_error is not None:
            row = [v + scale * e for v, e in zip(row, fe_error[n])]
        p.append(row)
        perm.append([A0 * (1.0 + perm_rk * t)] * 600)
    return dict(path='synthetic-' + tag, nn=len(x), nrec=len(times),
                x=list(x), times=list(times), p=p, perm=perm)


def selftest():
    print('=== Gate 4.7 offline selftest (analytic replay, no solver) ===')
    sources = tuple(read_text(p) for p in
                    (SRC_ASM, SRC_FLOW, SRC_CTX, SRC_DEFAULT, SRC_VDEF))
    decks = (read_text(DECK_MAT), read_text(DECK_CTRL))
    if not all(sources) or not all(decks):
        print('SELFTEST_RESULT: FAIL (source/deck files unreadable)')
        return 1
    # Real CI #109 geometry / time axis and its measured FE error, so the
    # positive replay uses the ACTUAL discretisation error, not a guess.
    x, times, fe_error = None, None, None
    art = os.path.join(_HERE, '..', '..', 'ci109_diag.zip')
    try:
        import zipfile
        import tempfile
        z = zipfile.ZipFile(art)
        tmp = os.path.join(tempfile.gettempdir(), 'g47_selftest_main.gen')
        if not os.path.exists(tmp):
            with open(tmp, 'wb') as fh:
                fh.write(z.read('out_main.gen'))
        real = load_exo(tmp)
        x, times = real['x'], real['times']
        fe_error = [[u - v for u, v in zip(real['p'][n],
                                          be_field('const', n, x))]
                    for n in range(len(times))]
        print('  replay basis: real CI #109 mesh (%d nodes, %d records)'
              % (len(x), len(times)))
    except Exception as exc:                           # noqa: BLE001
        print('  WARN: no local artifact (%r); using the analytic grid' % exc)
        x = [-5.0 + 10.0 * i / 355.0 for i in range(356)]
        times = [n * DT for n in range(51)]
    mat = _synthetic_mesh('full', x, times, fe_error, 1.0, PERM_RATE)
    ctrl = _synthetic_mesh('const', x, times, fe_error, 1.0, 0.0)
    mat_log = _synthetic_log(PERM_RATE, STOR_RATE, True)
    ctrl_log = _synthetic_log(0.0, 0.0, False)

    cases = [
        ('positive: both channels live',
         (sources, decks, (mat_log, ctrl_log), (mat, ctrl)), set()),
        ('NEG dead entry point (production field)',
         (sources, decks, (mat_log, ctrl_log),
          (_synthetic_mesh('const', x, times, fe_error, 1.0, 0.0), ctrl)),
         {'X5', 'X7'}),
        ('NEG only permeability wired',
         (sources, decks, (mat_log, ctrl_log),
          (_synthetic_mesh('onlyk', x, times, None, 1.0), ctrl)), {'X7'}),
        ('NEG only storativity wired',
         (sources, decks, (mat_log, ctrl_log),
          (_synthetic_mesh('onlys', x, times, None, 1.0), ctrl)), {'X7'}),
        ('NEG permField stays flat',
         (sources, decks, (mat_log, ctrl_log),
          (_flat_perm(mat), ctrl)), {'X5'}),
        ('NEG factor rows stay at 1 (dead VecSet)',
         (sources, decks, (_synthetic_log(0.0, 0.0, True), ctrl_log),
          (mat, ctrl)), {'X6'}),
        ('NEG bag still constant',
         (sources, decks, (_synthetic_log(PERM_RATE, STOR_RATE, False),
                           ctrl_log), (mat, ctrl)), {'X4'}),
        ('NEG deck has an extra change',
         (sources, (decks[0] + '\np0: 2.0e5\n', decks[1]),
          (mat_log, ctrl_log), (mat, ctrl)), {'X3'}),
        ('NEG source branch removed (static regression)',
         ((sources[0].replace('gOpt%permeabilityScaling',
                              'gOpt%sourceScaling', 1),)
          + sources[1:], decks, (mat_log, ctrl_log), (mat, ctrl)), {'X2'}),
        ('NEG control physics disturbed',
         (sources, decks, (mat_log, ctrl_log),
          (mat, _synthetic_mesh('full', x, times, None, 1.0, 0.0))),
         {'X8'}),
    ]
    all_ok = True
    for label, args, want_fail in cases:
        chk, _ = evaluate(*args)
        failing = {k for k, (ok, _) in chk.items() if not ok}
        ok = (failing == want_fail)
        all_ok &= ok
        print('  %-46s -> %s  failing=%s (want %s) %s'
              % (label, 'PASS' if ok else 'FAIL',
                 ','.join(sorted(failing)) or 'none',
                 ','.join(sorted(want_fail)) or 'none',
                 'OK' if ok else '<-- MISMATCH'))
    print('SELFTEST_RESULT: %s' % ('PASS' if all_ok else 'FAIL'))
    return 0 if all_ok else 1


def _flat_perm(mesh):
    d = dict(mesh)
    d['perm'] = [[A0] * len(row) for row in mesh['perm']]
    d['path'] = 'synthetic-flat-perm'
    return d


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
