#!/usr/bin/env python3
"""Gate 4.6 -- mu(T) / permField ALIAS removal verdict.

Pure standard library + the committed ``parse_cdf2.py`` (CI has no numpy /
scipy / PyYAML).  Combines STATIC source checks (the alias can always be
reintroduced by a small edit) with a FALSIFIABLE physical diagnostic run.

Modes
-----
  report : check_gate46_alias.py <visc.log> <ctrl.log> <visc.gen> <ctrl.gen>
      INFORMATIONAL, always exits 0.

  gate   : ... --gate <visc.log> <ctrl.log> <visc.gen> <ctrl.gen>
      HARD Gate 4.6, prints ``GATE46_RESULT: PASS|FAIL`` and exits 0 / 1.
      Missing / empty / unreadable / non-finite data is ALWAYS FAIL.

  selftest : ... --selftest
      Offline replay: real source text + synthetic fields, 1 positive and
      7 negative controls. No PETSc, no Exodus, no CI.

The bug
-------
``FlowCtx%cellDMScalSec`` is ONE SectionReal shared by three routines in
``m_Flow/m_MEF90_Flow.F90``:

  UpdateMuFromTemperature : writes mu(T)      into cellDMScalSec[cell]
  UpdatePerm              : READS muPtr       from cellDMScalSec[cell]
                            then writes k/mu  into the SAME slot
  UpdateStorativity       : writes S         into cellDMScalSec[cell]

So mu and permField ALIASED the same storage and were distinguished only by
call order.  ``cellMuField`` (the persistent Vec that actually holds mu(T))
was only used in an ``Associated()`` guard and never read back.

The fix
-------
All three routines now own a private ``SectionRealDuplicate``; the shared
section is a layout template only.  ``UpdateMuFromTemperature`` scatters its
private section FORWARD to ``cellMuField``; ``UpdatePerm`` scatters that Vec
REVERSE into its OWN private ``muSec``.  A non-positive mu (Vec is
zero-initialised before the first mu(T) assembly) falls back to mu_ref.

Physical anchor (terzaghi1D_viscosity.yaml)
------------------------------------------
The diagnostic deck is terzaghi1D_advection.yaml plus ONE option
``flow.temp_dep_viscosity: true``.  Both have a LIVE HeatXfer field
(temp_Offset 1, pressure_offset 0) with T in [0, 100] K.  Defaults are
T_ref=293 K, Sutherland C=600 K, mu_ref=1e-3 Pa.s, and the solver clamps
T to TMin = 0.5*T_ref = 146.5 K, so EVERY cell sees the same

    mu(T) = mu_ref * (TMin/T_ref)^1.5 * (T_ref+C)/(TMin+C)
          = 4.2293795e-4 Pa.s

Hence permField must be k_min/mu(T) = 1.1822065e-7, i.e. EXACTLY
mu_ref/mu(T) = 2.3644130 times the control deck's k_min/mu_ref = 5e-8.
If mu(T) is not read back from cellMuField (the alias returns, or mu falls
back to mu_ref), the ratio collapses to 1 or blows up -- both hard FAILs.

Criteria
--------
  V1  static: all 3 routines use private duplicate sections, no direct
      writes/reads of the shared section
  V2  static: mu is read back from cellMuField (SCATTER_REVERSE in
      UpdatePerm, SCATTER_FORWARD in UpdateMuFromTemperature); the legacy
      alias expression is absent
  V3  static: no SectionRealComplete/ToVec on the shared section anywhere
  V4  deck isolation (one added line) + solver bag really says TRUE
  V5  live, non-trivial temperature field in the diagnostic mesh
  V6  permField(visc)/permField(control) == 2.3644130 within 5%
  V7  per-cell permField * mu(T) == k_min within 3%
  V8  opt-out path untouched (control permField == 5e-8), 51 steps,
      no crash signature
"""
import difflib
import math
import os
import re
import sys

try:
    from parse_cdf2 import CDF2
except ImportError as exc:                                   # pragma: no cover
    CDF2 = None
    _CDF2_IMPORT_ERROR = exc

# --- Deck-physics constants (vDef/m_vDefDefault.F90 + terzaghi1D.yaml) -----
MU_REF = 1.0e-3                 # dynamicViscosity [Pa.s]
T_REF = 293.0                   # referenceTemperature [K]
SUTH_C = 600.0                  # sutherlandConstant [K]
K_MIN = 5.0e-11                 # perm_min [m^2]
T_MIN = 0.5 * T_REF             # solver clamp TMin = TRef*0.5
MU_T = (MU_REF * (T_MIN / T_REF) ** 1.5
        * (T_REF + SUTH_C) / (T_MIN + SUTH_C))
PERM_VISC = K_MIN / MU_T
PERM_CTRL = K_MIN / MU_REF
PERM_RATIO = MU_REF / MU_T
RATIO_TOL = 0.05
CONS_REL_MAX = 0.03
CTRL_REL_MAX = 0.02
TEMP_SPAN_MIN = 50.0
MIN_RECORDS = 50

BAG_RE = re.compile(r'^\s*temp_dep_viscosity\s*=\s*(TRUE|T)\b', re.I | re.M)
FATAL_RE = re.compile(
    r'Caught signal|SEGV|MPI_ABORT|PETSC ERROR|Unknown option|'
    r'attempt to deallocate|TSSolve failed|fortran runtime error', re.I)
LEGACY_ALIAS_RE = re.compile(
    r'SectionRealRestrict\s*\(\s*FlowCtx%cellDMScalSec\s*,\s*setIdx\(cell\)\s*,'
    r'\s*muPtr', re.I)
SHARED_COMPLETE_RE = re.compile(
    r'SectionReal(?:Complete|ToVec)\s*\(\s*FlowCtx%cellDMScalSec', re.I)
SUBS = ('MEF90FlowUpdatePerm', 'MEF90FlowUpdateMuFromTemperature',
        'MEF90FlowUpdateStorativity')


def read_text(path):
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as fh:
            return fh.read()
    except Exception:                                  # noqa: BLE001 - gate fails
        return ''


def strip_yaml(text):
    """Comments + blank lines removed; used for deck isolation."""
    out = []
    for line in text.splitlines():
        body = re.sub(r'#.*$', '', line).rstrip()
        if body.strip():
            out.append(body)
    return out


def deck_isolation(visc_text, ctrl_text):
    a, b = strip_yaml(ctrl_text), strip_yaml(visc_text)
    diff = [l for l in difflib.unified_diff(a, b, lineterm='')
            if l.startswith(('+', '-')) and not l.startswith(('+++', '---'))]
    ok = (len(diff) == 1 and diff[0].startswith('+')
          and re.match(r'^\+\s*temp_dep_viscosity\s*:\s*true\s*$',
                       diff[0], re.I) is not None)
    return ok, '%d stripped lines, %d differing: %s' % (
        len(b), len(diff), ' | '.join(d.strip() for d in diff[:3]))


def sub_body(src, name):
    pat = re.compile(
        r'^\s*Subroutine\s+' + re.escape(name) + r'\b.*?'
        r'^\s*End\s+Subroutine\s+' + re.escape(name) + r'\s*$',
        re.I | re.S | re.M)
    m = pat.search(src)
    return m.group(0) if m else ''


def static_checks(src):
    """V1, V2, V3 -- the anti-regression guard against re-aliasing."""
    chk = {}
    if not src:
        return {'V1': (False, 'source unreadable'),
                'V2': (False, 'source unreadable'),
                'V3': (False, 'source unreadable')}
    # V1: private sections everywhere
    ok, miss = True, []
    for name in SUBS:
        b = sub_body(src, name)
        if not b:
            ok = False
            miss.append(name + ':missing')
            continue
        if 'SectionRealDuplicate(FlowCtx%cellDMScalSec' not in b:
            ok = False
            miss.append(name + ':no-duplicate')
        if 'SectionRealDestroy(' not in b:
            ok = False
            miss.append(name + ':no-destroy')
        if re.search(r'SectionReal(?:Restrict|Restore)\s*\(\s*FlowCtx%cellDMScalSec',
                     b, re.I):
            ok = False
            miss.append(name + ':direct-shared-access')
    chk['V1'] = (ok, '3 routines private sections; problems=%s'
                 % (','.join(miss) or 'none'))

    # V2: mu really travels through the persistent cellMuField Vec
    bp = sub_body(src, 'MEF90FlowUpdatePerm')
    bm = sub_body(src, 'MEF90FlowUpdateMuFromTemperature')
    ok = (bool(bp) and bool(bm)
          and 'SCATTER_REVERSE' in bp and 'FlowCtx%cellMuField' in bp
          and 'SCATTER_FORWARD' in bm and 'FlowCtx%cellMuField' in bm
          and 'muVal' in bp and 'useMuT' in bp
          and not LEGACY_ALIAS_RE.search(src))
    chk['V2'] = (ok, 'REVERSE(cellMuField)->muSec in UpdatePerm, '
                 'FORWARD->cellMuField in UpdateMu, legacy alias present=%s'
                 % bool(LEGACY_ALIAS_RE.search(src)))

    # V3: shared section is complete/ToVec-free
    hits = SHARED_COMPLETE_RE.findall(src)
    chk['V3'] = (not hits, 'shared SectionRealComplete/ToVec occurrences=%d'
                 % len(hits))
    return chk


def load_exo(path):
    if CDF2 is None:
        raise ValueError('cannot import parse_cdf2: %s' % _CDF2_IMPORT_ERROR)
    if not path or not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise ValueError('missing/empty mesh: %s' % path)
    c = CDF2(path)
    temp, perm = [], []
    for n in range(c.numrecs):
        temp.append([float(v) for v in c.read_record('vals_nod_var1', n)])
        perm.append([float(v) for v in c.read_record('vals_elem_var4eb1', n)])
    for arr in (temp, perm):
        if any(not all(math.isfinite(v) for v in row) for row in arr):
            raise ValueError('non-finite Exodus field in %s' % path)
    return dict(nrec=c.numrecs, temp=temp, perm=perm,
                ncell=len(perm[0]) if perm else 0)


def field_checks(visc, ctrl):
    """V5..V8 from the two Exodus result meshes."""
    chk, met = {}, {}
    n = min(visc['nrec'], ctrl['nrec'])
    met.update(steps=n, cells=visc['ncell'])
    ok = (n >= MIN_RECORDS and ctrl['nrec'] == visc['nrec']
          and visc['ncell'] > 0 and ctrl['ncell'] == visc['ncell'])
    chk['V5'] = (ok and max(abs(v) for v in visc['temp'][n - 1]) >= TEMP_SPAN_MIN,
                 'steps=%d cells=%d, max|T| at last step=%.5g (want >=%.1f)'
                 % (n, visc['ncell'],
                    max(abs(v) for v in visc['temp'][n - 1]) if visc['temp'] else 0,
                    TEMP_SPAN_MIN))
    if not ok:
        for k in ('V6', 'V7', 'V8'):
            chk[k] = (False, 'mesh data incomplete')
        return chk, met

    pv, pc = visc['perm'][n - 1], ctrl['perm'][n - 1]
    ratios = [a / b for a, b in zip(pv, pc) if b > 0]
    ratio = sum(ratios) / len(ratios) if ratios else 0.0
    spread = (max(ratios) - min(ratios)) / ratio if ratios and ratio else float('inf')
    met.update(perm_ratio=ratio, perm_ratio_spread=spread,
               perm_visc=pv[0], perm_ctrl=pc[0])
    chk['V6'] = (len(ratios) == visc['ncell']
                 and abs(ratio - PERM_RATIO) / PERM_RATIO <= RATIO_TOL
                 and spread <= 0.05,
                 'perm(visc)/perm(ctrl)=%.6f (analytic %.6f, tol %.0f%%), '
                 'per-cell spread=%.3g'
                 % (ratio, PERM_RATIO, 100 * RATIO_TOL, spread))

    cons = max(abs(v * MU_T - K_MIN) / K_MIN for v in pv)
    met.update(mu_T=MU_T, consistency_max_rel=cons)
    chk['V7'] = (cons <= CONS_REL_MAX and min(pv) > 0,
                 'max|perm*mu(T)/k_min - 1|=%.5g (<= %.2f), '
                 'mu(T)=%.7g, min perm=%.5g'
                 % (cons, CONS_REL_MAX, MU_T, min(pv)))

    ctr = max(abs(v - PERM_CTRL) / PERM_CTRL for v in pc)
    met['control_perm_rel'] = ctr
    chk['V8'] = (ctr <= CTRL_REL_MAX,
                 'opt-out control perm=%.6g, max rel dev from %.6g = %.5g '
                 '(<= %.2f); %d steps both meshes'
                 % (pc[0], PERM_CTRL, ctr, CTRL_REL_MAX, n))
    return chk, met


def evaluate(src, visc_log, ctrl_log, visc_deck, ctrl_deck, visc, ctrl, err=None):
    chk, met = {}, {}
    if err:
        return chk, met, err
    chk.update(static_checks(src))

    iso, detail = deck_isolation(visc_deck, ctrl_deck)
    bag = bool(BAG_RE.search(visc_log))
    chk['V4'] = (iso and bag and bool(visc_log) and bool(ctrl_log),
                 'deck %s; bag temp_dep_viscosity=TRUE %s'
                 % (detail, bag))

    fatal = bool(FATAL_RE.search(visc_log) or FATAL_RE.search(ctrl_log))
    if fatal:
        chk['V4b'] = (False, 'crash/unknown-option signature present')
        return chk, met, None
    chk['V4b'] = (True, 'no crash/unknown-option signature')

    if visc is None or ctrl is None:
        for k in ('V5', 'V6', 'V7', 'V8'):
            chk[k] = (False, 'Exodus mesh unavailable')
        return chk, met, None
    fchk, fmet = field_checks(visc, ctrl)
    chk.update(fchk)
    met.update(fmet)
    return chk, met, None


def print_report(chk, met, err):
    print('=== Gate 4.6 -- mu(T) alias removal verdict ===')
    if err:
        print('  DATA ERROR: %s' % err)
    for name in sorted(chk):
        ok, detail = chk[name]
        print('  %-4s %-4s %s' % (name, 'PASS' if ok else 'FAIL', detail))
    if met:
        print('  --- metrics ---')
        for key in ('steps', 'cells', 'mu_T', 'perm_ratio',
                    'perm_ratio_spread', 'perm_visc', 'perm_ctrl',
                    'consistency_max_rel', 'control_perm_rel'):
            if key in met:
                print('    %-22s = %s' % (key, met[key]))
    return [k for k, v in chk.items() if not v[0]]


# --------------------------------------------------------------------------
# Offline selftest
# --------------------------------------------------------------------------
def _deck_text(flag):
    base = ['flow:', '    viscosity: 1.0e-03', '    storativity: 4.6e-10',
            '    temp_Offset: 1', 'heatxfer:', '    type: steadystate']
    if flag:
        base.insert(3, '    temp_dep_viscosity: true')
    return '\n'.join(base) + '\n'


def _fields(ratio, temp_span=100.0, ctrl_perm=None):
    ncell, nrec = 12, 51
    perm_v = [PERM_VISC * ratio] * ncell
    perm_c = [PERM_CTRL if ctrl_perm is None else ctrl_perm] * ncell
    temp = [[temp_span * (i + 1) / ncell * (n / (nrec - 1))
             for i in range(ncell)] for n in range(nrec)]
    return (dict(nrec=nrec, temp=temp, perm=[list(perm_v) for _ in range(nrec)],
                 ncell=ncell),
            dict(nrec=nrec, temp=temp, perm=[list(perm_c) for _ in range(nrec)],
                 ncell=ncell))


def selftest():
    src_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            '..', '..', 'm_Flow', 'm_MEF90_Flow.F90')
    src = read_text(src_path)
    good_log = '  temp_dep_viscosity =  TRUE; Enable temperature-dependent viscosity\n'
    bad_log = '  temp_dep_viscosity =  FALSE; Enable temperature-dependent viscosity\n'
    visc_deck, ctrl_deck = _deck_text(True), _deck_text(False)
    cases = [
        ('positive (fixed source, mu active)', src, good_log,
         visc_deck, ctrl_deck, 1.0, None, None, True),
        ('legacy alias reintroduced',
         src.replace('Call SectionRealRestrict(muSec,setIdx(cell),muPtr',
                     'Call SectionRealRestrict(FlowCtx%cellDMScalSec,setIdx(cell),muPtr'),
         good_log, visc_deck, ctrl_deck, 1.0, None, None, False),
        ('mu not read back (perm==control)', src, good_log,
         visc_deck, ctrl_deck, PERM_CTRL / PERM_VISC, None, None, False),
        ('perm blown up by 1e9 (mu<=0 guard)', src, good_log,
         visc_deck, ctrl_deck, 1e9, None, None, False),
        ('temperature field dead', src, good_log,
         visc_deck, ctrl_deck, 1.0, 0.0, None, False),
        ('bag still FALSE', src, bad_log,
         visc_deck, ctrl_deck, 1.0, None, None, False),
        ('deck not one-line isolated', src, good_log,
         ctrl_deck, ctrl_deck, 1.0, None, None, False),
        ('opt-out control perm disturbed', src, good_log,
         visc_deck, ctrl_deck, 1.0, None, PERM_CTRL * 1.5, False),
    ]
    print('=== Gate 4.6 offline selftest (real source + synthetic fields) ===')
    print('  source: %s (%d chars)' % (os.path.normpath(src_path), len(src)))
    print('  analytic mu(T)=%.8g perm=%.8g ratio=%.7f'
          % (MU_T, PERM_VISC, PERM_RATIO))
    all_ok = True
    for (label, s, log, vd, cd, ratio, span, cp, want) in cases:
        v, c = _fields(ratio, temp_span=(100.0 if span is None else span),
                       ctrl_perm=cp)
        chk, met, err = evaluate(s, log, log, vd, cd, v, c)
        fails = [k for k, x in chk.items() if not x[0]]
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
    usage = '<visc.log> <ctrl.log> <visc.gen> <ctrl.gen>'
    if len(args) != 4:
        if gate:
            print('GATE46_RESULT: FAIL (usage: --gate %s)' % usage)
            return 1
        print('usage: check_gate46_alias.py [--gate|--selftest] %s' % usage)
        return 0
    visc_log, ctrl_log, visc_gen, ctrl_gen = args
    here = os.path.dirname(os.path.abspath(__file__))
    src = read_text(os.path.join(here, '..', '..', 'm_Flow',
                                 'm_MEF90_Flow.F90'))
    vd = read_text(os.path.join(here, 'terzaghi1D_viscosity.yaml'))
    cd = read_text(os.path.join(here, 'terzaghi1D_advection.yaml'))
    vl, cl = read_text(visc_log), read_text(ctrl_log)
    err = None
    visc = ctrl = None
    try:
        visc, ctrl = load_exo(visc_gen), load_exo(ctrl_gen)
        if not vl or not cl:
            raise ValueError('diagnostic or control log unreadable/empty')
        if not vd or not cd or not src:
            raise ValueError('deck or source unreadable')
    except Exception as exc:                               # noqa: BLE001
        err = 'cannot build Gate 4.6 inputs: %r' % (exc,)
    chk, met, err = evaluate(src, vl, cl, vd, cd, visc, ctrl, err)
    fails = print_report(chk, met, err)
    if gate:
        if err or fails or len(chk) < 8:
            print('GATE46_RESULT: FAIL (%s)'
                  % (err if err else ','.join(fails) or 'incomplete criteria'))
            return 1
        print('GATE46_RESULT: PASS (V1-V3 private sections, V4 one-line deck '
              '+ bag TRUE, V5 live T, V6 ratio=%.5f vs %.5f, '
              'V7 max|perm*mu/k-1|=%.4g, V8 control perm=%.4g)'
              % (met.get('perm_ratio', 0), PERM_RATIO,
                 met.get('consistency_max_rel', 0), met.get('perm_ctrl', 0)))
        return 0
    print('VERDICT: %s (report mode, informational only)'
          % ('FAIL ' + ','.join(fails) if fails else 'PASS'))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
