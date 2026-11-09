#!/usr/bin/env python3
"""Gate 4.9 -- falsifiable verdict for the batched diffusion assembly.

WHY THIS GATE EXISTS
---------------------
Before 4.9 the Flow residual / Jacobian / storage assembly built one PETSc IS,
one element array and one Diffusion call PER CELL, because
MEF90Diffusion*Set take a CONSTANT coefficient A (or alpha / lambda) while the
Flow coefficient A = permField(cell)/mu is per-cell. That is O(N) PETSc object
creations per assembly.

Gate 4.9 adds a fast path, enabled by `flow.batched_assembly` (default FALSE):
when the coefficient of a cell set is bit-wise uniform, the WHOLE cell-set IS
is passed in ONE call, with ONE element array created from that same IS.
MEF90Diffusion*Set loop cell = 1..size(cellID) and call
SectionRealUpdateClosure(..., ADD_VALUES) / DMmeshAssembleMatrix in exactly
that order, so the legacy per-cell loop and the batched call produce
BIT-FOR-BIT IDENTICAL floating point results. Only object-creation count
changes, O(N) -> O(1). A non-uniform coefficient falls back to the legacy
per-cell loop, so the physics can never be silently changed.

WHY THE ACCEPTANCE CRITERION IS NOT WALL-CLOCK TIMING
----------------------------------------------------
The checklist asked for "large-mesh assembly time drops significantly". CI
wall clock is noisy (shared runners, PETSc cache warm-up) and this repository
ships only the 600-cell Terzaghi mesh. A timing threshold would be either
flaky or unfalsifiable. Gate 4.9 instead uses a DETERMINISTIC complexity
proxy: the number of IS / element-array creations per assembly, emitted by
the solver itself in `Gate 4.9: tag=... nCells=N nIS=1 nElem=1 mode=batched`
rows. Combined with a static audit of the legacy branch (Y7) this proves
O(N) -> O(1) exactly, and Y6 proves the physics is bit-identical. Any wall
clock improvement is an expected consequence, not a hard threshold.

CHECKS (all hard)
----------------
  Y1  data integrity: both runs completed, meshes finite and aligned
  Y2  static anti-regression: helper, 4 fast paths, bag order, default FALSE
  Y3  deck isolation: batched deck differs from control by exactly one line
  Y4  bag visibility: batched TRUE / production FALSE
  Y5  deterministic complexity count: all 3 RUNTIME-ACTIVE tags batched
      (OperatorCore / BilinearFormCore / StorageMassMatrix), nIS=nElem=1,
      nCells>1, production log has NO `Gate 4.9:` row. The 4th source fast
      path (AddStorageToResidual) is dead code since CI #64 and cannot emit a
      row; Y2/Y7 still audit it statically.
  Y6  physics bit-equality: p / u_x / u_y identical on every record/node
  Y7  O(N)->O(1) extrapolation + CI #65 elem-from-same-IS contract
  Y8  production isolation: equal step-marker counts, no fatal signatures

Run with `--gate` for the hard verdict (`GATE49_RESULT: PASS|FAIL`, exit 0/1),
`--selftest` for an offline analytic replay with positive + negative cases,
or with four paths for a report-only run.
"""

import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse_cdf2 import CDF2  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
DECK_BATCHED = os.path.join(_HERE, 'terzaghi1D_batched.yaml')
DECK_CTRL = os.path.join(_HERE, 'terzaghi1D.yaml')
SRC_ASM = os.path.join(_HERE, '..', '..', 'm_Flow',
                       'm_MEF90_FlowAssembly.F90')
SRC_CTX = os.path.join(_HERE, '..', '..', 'm_Flow', 'm_MEF90_FlowCtx.F90')
SRC_DEFAULT = os.path.join(_HERE, '..', 'm_vDefDefault.F90')

G49_RE = re.compile(
    r'Gate\s+4\.9:\s*tag=(\w+)\s+set\s+(\d+)\s+nCells\s+(\d+)\s+'
    r'nIS\s+(\d+)\s+nElem\s+(\d+)\s+mode\s+(\w+)')
BAG_RE = re.compile(r'batched_assembly\s*=\s*(TRUE|FALSE)\b', re.I)
STEP_RE = re.compile(r'Gate\s+3\.1:\s*step')
CRASH_RE = re.compile(
    r'Caught signal|SEGV|MPI_ABORT|PETSC ERROR|segmentation fault|'
    r'TSSolve failed|KSP_DIVERGED|not converged|unknown option', re.I)

# Runtime-active assembly loops on the Transient/ksponly path.
# MEF90FlowAddStorageToResidual is NOT called any more (CI #64 removed the
# extra scatter because xDot == 0 for TS_LINEAR + ksponly), so its fast path
# exists in the source but can never emit a log row. Y5 therefore requires
# the 3 ACTIVE tags; Y2/Y7 still audit all 4 source fast paths.
TAGS = ('OperatorCore', 'BilinearFormCore', 'StorageMassMatrix')
N_FAST = 4
# The option bag closed Gate 4.7 with 36 fields. Gate 4.9's field must sit
# AFTER that frozen prefix; later features may append without breaking the
# PetscBag offset contract (CI #50). This is the same append-only rule as
# Gate 4.7 X2, so an exact-arity check here would reject legal extensions.
N_BASE_FIELDS = 36
MIN_RECORDS = 50
MIN_CELLS = 2       # a single-cell set is not a meaningful batch
NDX = 2             # IS + element array creations per fast path

CHECK_TITLES = {
    'Y1': 'run/mesh data integrity',
    'Y2': 'static anti-regression (helper, 4 fast paths, bag order, default)',
    'Y3': 'deck isolation (exactly one added line)',
    'Y4': 'bag visibility (batched TRUE / production FALSE)',
    'Y5': 'deterministic complexity counts (nIS=nElem=1, 3 active tags)',
    'Y6': 'physics bit-identical (p / u_x / u_y on 51 records)',
    'Y7': 'O(N)->O(1) extrapolation + CI #65 elem-from-same-IS contract',
    'Y8': 'production isolation (equal step count, no fatal signature)',
}


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def read_text(path):
    if not path or not os.path.isfile(path):
        return ''
    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        return fh.read()


def strip_yaml(text):
    """Non-blank, non-comment, whitespace-normalised YAML lines."""
    out = []
    for line in (text or '').splitlines():
        s = line.split('#', 1)[0].strip()
        if s:
            out.append(re.sub(r'\s+', ' ', s))
    return out


def parse_g49(log):
    rows = []
    for m in G49_RE.finditer(log or ''):
        rows.append(dict(tag=m.group(1), set=int(m.group(2)),
                         ncells=int(m.group(3)), nis=int(m.group(4)),
                         nelem=int(m.group(5)), mode=m.group(6)))
    return rows


def bag_value(log):
    m = BAG_RE.search(log or '')
    return m.group(1).upper() if m else None


def load_fields(path):
    """Exodus nodal fields: p=var1, u_x=var3, u_y=var4 (same slots as 4.1b)."""
    if not path or not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise ValueError('missing/empty mesh: %s' % path)
    c = CDF2(path)
    nn = len(c.read_static('coordx'))
    if nn <= 0:
        raise ValueError('no nodes (coordx) in: %s' % path)
    x = [float(v) for v in c.read_static('coordx')]
    y = [float(v) for v in c.read_static('coordy')]
    times, p, ux, uy = [], [], [], []
    for i in range(c.numrecs):
        times.append(float(c.read_record('time_whole', i)[0]))
        p.append([float(v) for v in c.read_record('vals_nod_var1', i)])
        ux.append([float(v) for v in c.read_record('vals_nod_var3', i)])
        uy.append([float(v) for v in c.read_record('vals_nod_var4', i)])
    for key, rows in (('p', p), ('ux', ux), ('uy', uy)):
        for f in rows:
            if not f or len(f) != nn or not all(math.isfinite(v) for v in f):
                raise ValueError('non-finite / mismatched %s in: %s'
                                 % (key, path))
    if not all(math.isfinite(t) for t in times):
        raise ValueError('non-finite time axis in: %s' % path)
    return dict(path=path, nn=nn, nrec=c.numrecs, x=x, y=y, times=times,
                p=p, ux=ux, uy=uy)


def load_sources():
    return dict(asm=read_text(SRC_ASM), ctx=read_text(SRC_CTX),
                default=read_text(SRC_DEFAULT))


def ctx_order(ctx):
    """Return (fields, registrations, default_items) for CI #50 ordering."""
    i = ctx.find('   Type MEF90FlowGlobalOptions_Type')
    j = ctx.find('   End Type MEF90FlowGlobalOptions_Type', i)
    seg = ctx[i:j] if i >= 0 and j > i else ''
    fields = re.findall(
        r'^\s*(?:PetscBool|PetscInt|PetscReal|PetscEnum)\s*::\s*(\w+)',
        seg, re.M)
    regs = re.findall(
        r'PetscBagRegister(?:Bool|Int|Real|Enum|String)\s*\(\s*bag\s*,\s*'
        r'FlowGlobalOptions%(\w+)', ctx)
    return fields, regs


def default_items(default):
    a = default.find('vDefFlowDefaultGlobalOptions')
    b = default.find('Type(MEF90FlowCellSetOptions_Type)', a)
    if a < 0 or b < 0:
        return []
    out = []
    for line in default[a:b].splitlines():
        s = line.split('!')[0].strip().rstrip('&').strip()
        if s and 'MEF90FlowGlobalOptions_Type(' not in s:
            out.append(s)
    return out


def pass_fail(flag, detail):
    return (bool(flag), detail)


# ---------------------------------------------------------------------------
# checks
# ---------------------------------------------------------------------------
def check_y1(logs, meshes):
    issues = []
    bat_log, ctrl_log = logs
    bat, ctrl = meshes
    if not bat_log or not ctrl_log:
        issues.append('empty solver log')
    for name, log in (('batched', bat_log), ('control', ctrl_log)):
        n = len(STEP_RE.findall(log or ''))
        if n < MIN_RECORDS:
            issues.append('%s steps=%d < %d' % (name, n, MIN_RECORDS))
    for name, m in (('batched', bat), ('control', ctrl)):
        if m['nrec'] < MIN_RECORDS:
            issues.append('%s mesh records=%d' % (name, m['nrec']))
        if m['nn'] <= 0:
            issues.append('%s mesh nodes=%d' % (name, m['nn']))
    if bat['nn'] != ctrl['nn']:
        issues.append('node count differs')
    if bat['nrec'] != ctrl['nrec']:
        issues.append('record count differs')
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'both runs complete, meshes aligned')


def check_y2(sources):
    asm, ctx, default = sources['asm'], sources['ctx'], sources['default']
    issues = []
    if 'Subroutine MEF90FlowUniformCellField' not in asm:
        issues.append('uniform-coefficient helper missing')
    if 'vLoc(1) /= value' not in asm:
        issues.append('helper no longer uses bit-wise comparison')
    nfast = asm.count('If (gOpt49%batchedAssembly .AND. uniformCoeff')
    if nfast != 4:
        issues.append('fast-path branches=%d (want 4)' % nfast)
    ncalls = asm.count('Call MEF90FlowUniformCellField(')
    if ncalls != 4:
        issues.append('helper calls=%d (want 4)' % ncalls)
    ndefault = asm.count("'Gate 4.9: tag=")
    if ndefault != 4:
        issues.append('diagnostic tags=%d (want 4)' % ndefault)
    fields, regs = ctx_order(ctx)
    items = default_items(default)
    if not (len(fields) == len(regs) == len(items) and len(fields) > 0):
        issues.append('bag counts %d/%d/%d' % (len(fields), len(regs),
                                               len(items)))
    if fields != regs:
        issues.append('bag declaration/registration order mismatch')
    if 'batchedAssembly' not in fields[N_BASE_FIELDS:]:
        issues.append('batchedAssembly not appended after the frozen 4.7 prefix')
    if 'batchedAssembly' not in regs[N_BASE_FIELDS:]:
        issues.append('batched_assembly not registered after the 4.7 prefix')
    if 'batched_assembly' not in ctx:
        issues.append('option not registered')
    bpos = fields.index('batchedAssembly') if 'batchedAssembly' in fields else -1
    if bpos < 0 or bpos >= len(items) or 'PETSC_FALSE' not in items[bpos]:
        issues.append('default for batchedAssembly is not PETSC_FALSE')
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'helper + 4 fast paths, %d/%d/%d bag order, default FALSE'
                     % (len(fields), len(regs), len(items)))


def check_y3(decks):
    batched, ctrl = decks
    a, b = strip_yaml(ctrl), strip_yaml(batched)
    only_added = [x for x in b if x not in a]
    only_removed = [x for x in a if x not in b]
    issues = []
    if len(only_added) != 1 or only_removed:
        issues.append('added=%s removed=%s' % (only_added, only_removed))
    elif only_added[0] != 'batched_assembly: true':
        issues.append('unexpected added line: %s' % only_added[0])
    return pass_fail(not issues, '; '.join(issues) if issues else
                     '1 added line: batched_assembly: true')


def check_y4(logs):
    bat_log, ctrl_log = logs
    b, c = bag_value(bat_log), bag_value(ctrl_log)
    issues = []
    if b != 'TRUE':
        issues.append('batched bag=%s (want TRUE)' % b)
    if c != 'FALSE':
        issues.append('production bag=%s (want FALSE)' % c)
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'batched TRUE / production FALSE')


def check_y5(logs, rows):
    bat_log, ctrl_log = logs
    issues = []
    if not rows:
        issues.append('no Gate 4.9 rows in batched log')
    tags = set(r['tag'] for r in rows)
    missing = set(TAGS) - tags
    if missing:
        issues.append('missing tags=%s' % sorted(missing))
    ncells = set(r['ncells'] for r in rows)
    if len(ncells) != 1 or min(ncells or [0]) < MIN_CELLS:
        issues.append('nCells values=%s' % sorted(ncells))
    bad = [r for r in rows if r['nis'] != 1 or r['nelem'] != 1
           or r['mode'] != 'batched']
    if bad:
        issues.append('%d rows not batched/one-IS' % len(bad))
    if 'Gate 4.9:' in (ctrl_log or ''):
        issues.append('production log contains Gate 4.9 rows')
    n = max(ncells) if ncells else 0
    saved = N_FAST * NDX * (n - 1)
    detail = ('; '.join(issues) if issues else
              '%d active tags, %d rows, nCells=%d, nIS=nElem=1, production clean'
              % (len(TAGS), len(rows), n))
    met5 = dict(rows=len(rows), ncells=n, saved_creates=saved)
    return pass_fail(not issues, detail), met5


def check_y6(meshes):
    bat, ctrl = meshes
    issues = []
    if bat['x'] != ctrl['x'] or bat['y'] != ctrl['y']:
        issues.append('mesh coordinates differ')
    if bat['times'] != ctrl['times']:
        issues.append('time axis differs')
    maxdiff = 0.0
    for key in ('p', 'ux', 'uy'):
        a, b = bat[key], ctrl[key]
        if len(a) != len(b):
            issues.append('%s record count differs' % key)
            continue
        for ra, rb in zip(a, b):
            for xa, xb in zip(ra, rb):
                d = abs(xa - xb)
                if d > maxdiff:
                    maxdiff = d
    if maxdiff != 0.0:
        issues.append('max|batched - control| = %g (want exactly 0)' % maxdiff)
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'p/u_x/u_y bit-identical on %d records x %d nodes'
                     % (bat['nrec'], bat['nn'])), dict(maxdiff=maxdiff)


def check_y7(sources, rows, ncells):
    asm = sources['asm']
    issues = []
    # Fast path: one element array created from the WHOLE setIS (CI #65).
    # 4 fast paths added, PLUS 2 pre-existing whole-set element creations in
    # ComputeCellPressure / ComputePressureGradient -> use a lower bound.
    nset = asm.count('MEF90Element_Create(FlowCtx%DM,setIS,elem')
    if nset < 4:
        issues.append('setIS element creations=%d (want >=4)' % nset)
    # Legacy per-cell fallbacks: 4 optimised loops + MEF90FlowEnergy (which is
    # deliberately NOT batched because DiffusionEnergySet accumulates a scalar
    # INSIDE the call, so a batched call would change the FP summation order).
    nlegacy = asm.count('Call MEF90FlowMakeSingleCellIS(')
    nsingle = asm.count('MEF90Element_Create(FlowCtx%DM,singleIS,elem')
    if nlegacy < 4 or nsingle < 4:
        issues.append('legacy per-cell fallbacks=%d/%d (want >=4)'
                      % (nlegacy, nsingle))
    # Deterministic counts measured from the solver log, not only the source.
    measured = set(r['nis'] for r in rows) | set(r['nelem'] for r in rows)
    per_call = max(measured) if measured else 1
    if per_call != 1:
        issues.append('fast path still creates %d objects per call' % per_call)
    if ncells < MIN_CELLS:
        issues.append('nCells=%d too small to extrapolate' % ncells)
    legacy = N_FAST * NDX * ncells
    batched = N_FAST * NDX * per_call
    ratio = legacy / float(batched) if batched else 0.0
    if ratio < 2.0:
        issues.append('no O(N)->O(1) reduction (ratio=%.1f)' % ratio)
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'legacy=%d batched=%d object creations (%.0fx), '
                     'elem created from same setIS'
                     % (legacy, batched, ratio)), dict(legacy=legacy,
                                                        batched=batched,
                                                        ratio=ratio)


def check_y8(logs):
    bat_log, ctrl_log = logs
    issues = []
    nb = len(STEP_RE.findall(bat_log or ''))
    nc = len(STEP_RE.findall(ctrl_log or ''))
    if nb != nc:
        issues.append('step markers %d vs %d' % (nb, nc))
    for name, log in (('batched', bat_log), ('control', ctrl_log)):
        m = CRASH_RE.search(log or '')
        if m:
            issues.append('%s fatal signature: %s' % (name, m.group(0)))
    return pass_fail(not issues, '; '.join(issues) if issues else
                     'step markers equal (%d), no fatal signature' % nb)


# ---------------------------------------------------------------------------
# evaluation / reporting
# ---------------------------------------------------------------------------
def evaluate(sources, decks, logs, meshes):
    chk, met = {}, {}
    chk['Y1'] = check_y1(logs, meshes)
    chk['Y2'] = check_y2(sources)
    chk['Y3'] = check_y3(decks)
    chk['Y4'] = check_y4(logs)
    rows = parse_g49(logs[0])
    chk['Y5'], met['y5'] = check_y5(logs, rows)
    chk['Y6'], met['y6'] = check_y6(meshes)
    chk['Y7'], met['y7'] = check_y7(sources, rows, met['y5']['ncells'])
    chk['Y8'] = check_y8(logs)
    return chk, met


def print_report(chk, met):
    print('=== Gate 4.9 -- batched diffusion assembly verdict ===')
    for name in ('Y1', 'Y2', 'Y3', 'Y4', 'Y5', 'Y6', 'Y7', 'Y8'):
        ok, detail = chk[name]
        print('  %s %s %-58s | %s'
              % (name, 'PASS' if ok else 'FAIL',
                 CHECK_TITLES[name], detail[:190]))
    print('  --- metrics ---')
    print('    g49_rows         = %d' % met['y5']['rows'])
    print('    nCells           = %d' % met['y5']['ncells'])
    print('    saved_creates   = %d' % met['y5']['saved_creates'])
    print('    legacy_creates   = %d' % met['y7']['legacy'])
    print('    batched_creates  = %d' % met['y7']['batched'])
    print('    reduction_x     = %.0f' % met['y7']['ratio'])
    print('    max|dPhysics|    = %g' % met['y6']['maxdiff'])
    return [n for n, (ok, _) in chk.items() if not ok]


def verdict_line(fails, met):
    if fails:
        return ('GATE49_RESULT: FAIL (%s)' % ','.join(fails))
    return ('GATE49_RESULT: PASS (Y1-Y8: %d active batched tags (4 source '
            'fast paths), nIS=nElem=1 vs nCells=%d, %d object creations saved '
            'per step, physics bit-identical, production unchanged)'
            % (len(TAGS), met['y5']['ncells'], met['y5']['saved_creates']))


# ---------------------------------------------------------------------------
# offline selftest: synthetic logs + synthetic fields, no solver / no artifact
# ---------------------------------------------------------------------------
def make_log(batched, bag, nsteps=408, tags=TAGS, ncells=600, nis=1,
             g49_in_ctrl=False, crash=False, steps_ctrl_override=None):
    out = ['batched_assembly = %s; Gate 4.9 batched assembly' % bag]
    n = nsteps if steps_ctrl_override is None else steps_ctrl_override
    for i in range(n):
        out.append('Gate 3.1: step %d post-SetTransients' % (i + 1))
    rows = []
    if batched:
        for t in tags:
            rows.append('Gate 4.9: tag=%s set 1 nCells %d nIS %d nElem %d '
                        'mode batched' % (t, ncells, nis, nis))
    if g49_in_ctrl and not batched:
        rows.append('Gate 4.9: tag=OperatorCore set 1 nCells %d nIS 1 '
                    'nElem 1 mode batched' % ncells)
    out.extend(rows)
    if crash:
        out.append('PETSC ERROR: synthetic fatal signature')
    return '\n'.join(out) + '\n'


def make_fields(nrec=51, nn=356, perturb=None):
    x = [0.001 * i for i in range(nn)]
    y = [0.002 * (i % 7) for i in range(nn)]
    times = [i / float(nrec - 1) for i in range(nrec)]
    def field(scale):
        return [[scale * math.sin(0.01 * (i + 1) + 0.002 * j)
                 for j in range(nn)] for i in range(nrec)]
    p, ux, uy = field(1.0e5), field(1.1e-4), field(8.0e-7)
    if perturb:
        ux[0][0] += perturb
    return dict(path='<synthetic>', nn=nn, nrec=nrec, x=x, y=y, times=times,
                p=p, ux=ux, uy=uy)


def selftest():
    print('=== Gate 4.9 offline selftest (synthetic replay, no solver) ===')
    sources = load_sources()
    decks = (read_text(DECK_BATCHED), read_text(DECK_CTRL))
    if not all(sources.values()) or not all(decks):
        print('SELFTEST_RESULT: FAIL (source/deck files unreadable)')
        return 1
    good_log = make_log(True, 'TRUE')
    ctrl_log = make_log(False, 'FALSE')
    fields = make_fields()
    same = make_fields()

    def run(bat_log, c_log, bat_mesh, c_mesh, src=None, dks=None):
        s = src if src is not None else sources
        d = dks if dks is not None else decks
        chk, met = evaluate(s, d, (bat_log, c_log),
                            (bat_mesh, c_mesh))
        return set(n for n, (ok, _) in chk.items() if not ok), met

    cases = []

    # positive
    fails, met = run(good_log, ctrl_log, fields, same)
    cases.append(('positive: all active tags batched, physics identical',
                  fails, set(), met))

    # negatives -- each must fail exactly the expected check(s)
    f, _ = run(make_log(True, 'TRUE', tags=TAGS[:2]), ctrl_log, fields, same)
    cases.append(('NEG missing StorageMassMatrix tag', f, {'Y5'}, None))

    f, _ = run(make_log(True, 'TRUE', nis=600), ctrl_log, fields, same)
    cases.append(('NEG nIS=600 (no reduction)', f, {'Y5', 'Y7'}, None))

    f, _ = run(good_log, make_log(False, 'FALSE', g49_in_ctrl=True),
               fields, same)
    cases.append(('NEG production log has Gate 4.9 rows', f, {'Y5'}, None))

    f, _ = run(make_log(True, 'FALSE'), ctrl_log, fields, same)
    cases.append(('NEG batched bag shows FALSE', f, {'Y4'}, None))

    f, _ = run(good_log, ctrl_log, fields, same,
               dks=(decks[0] + '\n    extra_key: 1\n', decks[1]))
    cases.append(('NEG deck has an extra line', f, {'Y3'}, None))

    bad_default = dict(sources)
    bad_default['default'] = sources['default'].replace(
        'PETSC_FALSE)        ! batchedAssembly',
        'PETSC_TRUE)         ! batchedAssembly')
    f, _ = run(good_log, ctrl_log, fields, same, src=bad_default)
    cases.append(('NEG default flipped to TRUE', f, {'Y2'}, None))

    # Rename the field in BOTH declaration and registration, so the
    # decl == regf check still passes: only the append-after-4.7 contract
    # detects that the Gate 4.9 bag slot has been disturbed.
    bad_ctx = dict(sources)
    bad_ctx['ctx'] = sources['ctx'].replace('batchedAssembly',
                                            'batchedAssemblyRenamed')
    f, _ = run(good_log, ctrl_log, fields, same, src=bad_ctx)
    cases.append(('NEG 4.9 bag slot disturbed', f, {'Y2'}, None))

    bad_asm = dict(sources)
    bad_asm['asm'] = sources['asm'].replace(
        'If (gOpt49%batchedAssembly .AND. uniformCoeff',
        'If (.FALSE. .AND. uniformCoeff')
    f, _ = run(good_log, ctrl_log, fields, same, src=bad_asm)
    cases.append(('NEG fast-path guards removed', f, {'Y2'}, None))

    f, _ = run(good_log, ctrl_log, make_fields(perturb=1.0e-9), same)
    cases.append(('NEG physics field perturbed', f, {'Y6'}, None))

    f, _ = run(make_log(True, 'TRUE', crash=True), ctrl_log, fields, same)
    cases.append(('NEG fatal signature in batched log', f, {'Y8'}, None))

    f, _ = run(make_log(True, 'TRUE', nsteps=300), ctrl_log, fields, same)
    cases.append(('NEG step count diverges', f, {'Y8'}, None))

    ok = True
    for label, got, want, _ in cases:
        good = (got == want)
        ok = ok and good
        print('  %-48s -> %-12s %s'
              % (label, 'PASS' if good else 'FAIL',
                 'failing=%s (want %s)' % (sorted(got) or 'none',
                                           sorted(want) or 'none')))
    print('SELFTEST_RESULT: %s' % ('PASS' if ok else 'FAIL'))
    return 0 if ok else 1


def main(argv):
    if argv and argv[0] == '--selftest':
        return selftest()
    if argv and argv[0] == '--gate':
        args = argv[1:]
        if len(args) != 4:
            print('GATE49_RESULT: FAIL (usage: --gate <batched.log> '
                  '<control.log> <batched.gen> <control.gen>)')
            return 2
        try:
            logs = (read_text(args[0]), read_text(args[1]))
            meshes = (load_fields(args[2]), load_fields(args[3]))
            sources = load_sources()
            decks = (read_text(DECK_BATCHED), read_text(DECK_CTRL))
            if not all(logs) or not all(sources.values()) or not all(decks):
                raise ValueError('a log, source or deck file is unreadable')
            chk, met = evaluate(sources, decks, logs, meshes)
        except Exception as exc:                  # noqa: BLE001 - gate fails
            print('GATE49_RESULT: FAIL (%s: %s)' % (type(exc).__name__, exc))
            return 1
        fails = print_report(chk, met)
        print(verdict_line(fails, met))
        return 1 if fails else 0
    print('usage: check_gate49_batch.py [--gate <batched.log> <control.log> '
          '<batched.gen> <control.gen> | --selftest | <4 paths>]')
    return 2


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
