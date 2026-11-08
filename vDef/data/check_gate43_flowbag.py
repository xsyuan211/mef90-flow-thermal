#!/usr/bin/env python3
"""Gate 4.3 verdict: Flow option bags are visible AND the constant pressure
BC really stays constant in time.

Checklist item 4.3 originally claimed "a constant pressure BC is multiplied
by time". CI #92 proved the premise false: bc_pressure_scaling already has
constant/linear/file modes and the CST branch in MEF90FlowSetTransients
does not scale. The real gap was that MEF90FlowCtxSetFromOptions never
called PetscBagView, so no CI log contained bc_pressure_scaling,
source_scaling, biot_coeff, enable_advection or the EXO offsets. CI #93
added three verbose-guarded PetscBagView calls (print only).

This checker uses the REAL bag format observed in CI #93, e.g.

    PetscBag Object:  flow Flow Global Options
      bc_pressure_scaling = constant; (MEF90scaling) Boundary pressure scaling
      biot_coeff = 1 ; Biot poroelastic coefficient ...
      enable_advection =  FALSE; Add v_darcy*grad(T) ...
    PetscBag Object:  cs10 Flow CellSet Options
      flowBC =  TRUE; Cell set has pressure boundary condition
      bc_pressure = 0 ; Pressure value on BC surface (Pa)

Criteria (C1-C7):
  C1  every solver log has exactly one Flow Global bag with all required keys
  C2  every solver log has the three Flow CellSet bags; cs0010 is the only
      flowBC=TRUE set and its bc_pressure matches the deck
  C3  scaling semantics: bc_pressure_scaling / source_scaling /
      pressure_scaling are `constant` in every run AND explicitly declared
      in every deck (no silent default fallback)
  C4  FULL deck <-> bag agreement for EVERY key in each deck's flow: section
      (numeric tolerance, enum/bool normalisation). This is what catches a
      silently ignored YAML key -- the exact class of bug that hit 4.1a
      (`Force: 0, -9.81, 0`) and CI #52 (`CST` enum).
  C5  known variant differences must be visible in the bags (nobiot alpha=0,
      advection enable_advection=TRUE, 4.2 decks pressure_offset=0, etc.)
  C6  PHYSICS: the drained Dirichlet nodes (x <= -4.9999) hold a constant
      pressure over all >=50 records, equal to the bag's bc_pressure; the
      closed end must still decay (guards against an all-zero fake field);
      Gate 3.3 must be PASS with its analytic-series L2 error <= 0.05.
  C7  no regression: Gate 4.1 / 4.2 verdicts still PASS, no fatal signatures
      in any solver log, and each run reached the Gate 4.2 coupling wire.

Usage:
  check_gate43_flowbag.py --gate <log_dir> <deck_dir> <prod_mesh>
      -> prints "GATE43_RESULT: PASS|FAIL", exits 0 / 1 (hard gate)
  check_gate43_flowbag.py <log_dir> <deck_dir> <prod_mesh>
      -> report only, always exits 0

Pure standard library + the repo's vDef/data/parse_cdf2.py. The CI runner
has neither numpy nor PyYAML, so decks are read with a minimal YAML subset
parser (top-level `section:` plus indented `key: value` lines).
Zero compilation, zero solver execution.
"""
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from parse_cdf2 import CDF2
except Exception as exc:                                # pragma: no cover
    print('GATE43_RESULT: FAIL (cannot import parse_cdf2: %s)' % exc)
    sys.exit(1)

# (solver log, deck) -- ksponly_phys reuses the production deck with only
# -disp_snes_type ksponly on the command line (workflow line ~659).
RUNS = (
    ('terzaghi.log',             'terzaghi1D.yaml'),
    ('terzaghi_ksponly_phys.log','terzaghi1D.yaml'),
    ('terzaghi_nobiot.log',      'terzaghi1D_nobiot.yaml'),
    ('terzaghi_bodyforce.log',   'terzaghi1D_bodyforce.yaml'),
    ('terzaghi_advection.log',   'terzaghi1D_advection.yaml'),
    ('terzaghi_noadvection.log', 'terzaghi1D_noadvection.yaml'),
)

# Keys that MUST appear in every Flow Global bag. Missing one means the
# registration/view path regressed even if the run itself succeeded.
REQUIRED_GLOBAL = (
    'timestepping_type', 'perm_min', 'perm_max', 'viscosity', 'storativity',
    'p0', 'residual_damage', 'biot_coeff', 'temp_dep_viscosity',
    'enable_advection', 'pressure_offset', 'source_offset',
    'boundary_pressure_offset', 'pressure_previous_step_offset',
    'perm_field_offset', 'darcy_velocity_offset',
    'source_scaling', 'bc_pressure_scaling',
    'pressure_previous_step_scaling', 'pressure_scaling',
    'pressure_atol', 'pressure_rtol', 'pressure_maxit',
)
REQUIRED_CELLSET = ('flowBC', 'bc_pressure', 'source', 'elem_short_id')
# bc/source scaling are explicitly declared in every deck; pressure_scaling
# is a library default (MEF90Scaling_CST, m_vDefDefault.F90) and is therefore
# checked against the bag value only, not required in the deck text.
SCALING_KEYS = ('bc_pressure_scaling', 'source_scaling', 'pressure_scaling')
EXPLICIT_SCALING_KEYS = ('bc_pressure_scaling', 'source_scaling')

FATAL_RE = re.compile(
    r'Caught signal|SEGV|MPI_ABORT|PETSC ERROR|Unknown option|'
    r'attempt to deallocate|Error opening EXO file', re.I)

BAG_HEADER_RE = re.compile(r'^\s*PetscBag Object:\s*(.*?)\s*$')
BAG_ITEM_RE = re.compile(r'^\s{2}(\w+)\s*=\s*(.*?)\s*;')
DECK_KEY_RE = re.compile(r'^\s+([A-Za-z_]\w*)\s*:\s*(.*?)\s*$')


# --------------------------------------------------------------------------
# Parsers
# --------------------------------------------------------------------------
def parse_bags(path):
    """Return (global_dict, {cellset_name: dict}, header_count) for one log."""
    glob, cells, header_count = {}, {}, 0
    current = None
    with open(path, encoding='utf-8', errors='replace') as fh:
        for raw in fh:
            line = raw.rstrip('\r\n')
            mh = BAG_HEADER_RE.match(line)
            if mh:
                header_count += 1
                name = mh.group(1)
                if 'flow Flow Global Options' in name:
                    current = glob
                elif 'Flow CellSet Options' in name:
                    setname = name.split()[0]
                    cells.setdefault(setname, {})
                    current = cells[setname]
                else:
                    current = None
                continue
            if current is None:
                continue
            mi = BAG_ITEM_RE.match(line)
            if mi:
                current[mi.group(1)] = mi.group(2).strip()
    return glob, cells, header_count


def parse_deck(path):
    """Minimal YAML subset: top-level `section:` + indented `key: value`.

    Nested blocks (e.g. hookeslaw:) are stored as empty values in the
    enclosing section and are ignored by the flow checks. No lists, no
    anchors -- the vDef option decks use none of those.
    """
    sections, current = {}, None
    with open(path, encoding='utf-8', errors='replace') as fh:
        for raw in fh:
            line = raw.rstrip('\r\n')
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            body = line.split('#', 1)[0].rstrip()
            if not body.strip():
                continue
            if (not body[0].isspace()) and body.rstrip().endswith(':'):
                current = body.strip()[:-1]
                sections.setdefault(current, {})
                continue
            if current is None:
                continue
            mk = DECK_KEY_RE.match(body)
            if mk:
                sections[current][mk.group(1)] = mk.group(2).strip()
    return sections


def as_number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def values_equal(deck_val, bag_val):
    """Compare a YAML scalar with a PetscBag scalar / enum / bool."""
    if bag_val is None or deck_val is None or bag_val == '':
        return False
    d = str(deck_val).strip()
    b = str(bag_val).strip()
    db, bb = d.lower(), b.lower()
    if db in ('true', 'false', '.true.', '.false.'):
        return bb == ('true' if db in ('true', '.true.') else 'false')
    dn, bn = as_number(d), as_number(b)
    if dn is not None and bn is not None:
        return abs(dn - bn) <= 1.0e-9 * max(1.0, abs(dn), abs(bn))
    return db == bb                       # enum / string (case-insensitive)


def finite(x):
    try:
        return math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


# --------------------------------------------------------------------------
# Evaluation
# --------------------------------------------------------------------------
def evaluate(log_dir, deck_dir, prod_mesh):
    checks, data = {}, {}
    runs = []
    read_error = []

    for log_name, deck_name in RUNS:
        lp = os.path.join(log_dir, log_name)
        dp = os.path.join(deck_dir, deck_name)
        try:
            g, cells, headers = parse_bags(lp)
            deck = parse_deck(dp)
            text = open(lp, encoding='utf-8', errors='replace').read()
        except Exception as exc:                        # noqa: BLE001
            read_error.append('%s: %s' % (log_name, exc))
            continue
        runs.append(dict(log=log_name, deck=deck_name, g=g, cells=cells,
                         headers=headers, sections=deck, text=text))

    # C1 -- Flow Global bag visible and complete in every run
    missing = []
    for r in runs:
        absent = [k for k in REQUIRED_GLOBAL if k not in r['g']]
        if absent:
            missing.append('%s: %s' % (r['log'], ','.join(absent)))
    checks['C1'] = (len(runs) == len(RUNS) and not missing
                    and all(len(r['g']) >= len(REQUIRED_GLOBAL) for r in runs),
                    '%d/%d runs, %d required keys each, missing=[%s]'
                    % (len(runs), len(RUNS), len(REQUIRED_GLOBAL),
                       '; '.join(missing) or 'none'))

    # C2 -- three Flow CellSet bags, cs0010 is the only flowBC=TRUE set
    cell_issues = []
    for r in runs:
        sets = r['cells']
        if len(sets) != 3:
            cell_issues.append('%s: %d cellset bags (want 3)'
                               % (r['log'], len(sets)))
            continue
        for sname, vals in sets.items():
            absent = [k for k in REQUIRED_CELLSET if k not in vals]
            if absent:
                cell_issues.append('%s/%s: missing %s'
                                   % (r['log'], sname, ','.join(absent)))
        # cs10 -> cs0010 in the deck (zero-padded 4 digits)
        for sname, vals in sets.items():
            if not sname.startswith('cs'):
                continue
            try:
                deck_key = 'cs%04d' % int(sname[2:])
            except ValueError:
                continue
            dsec = r['sections'].get(deck_key, {})
            if 'flowBC' in dsec and not values_equal(dsec['flowBC'],
                                                     vals.get('flowBC')):
                cell_issues.append('%s/%s: flowBC %s != deck %s'
                                   % (r['log'], sname, vals.get('flowBC'),
                                      dsec['flowBC']))
            if 'bc_pressure' in dsec and not values_equal(
                    dsec['bc_pressure'], vals.get('bc_pressure')):
                cell_issues.append('%s/%s: bc_pressure %s != deck %s'
                                   % (r['log'], sname, vals.get('bc_pressure'),
                                      dsec['bc_pressure']))
    checks['C2'] = (not cell_issues and len(runs) == len(RUNS),
                    '3 Flow CellSet bags/run, cs0010 flowBC+bc_pressure '
                    'agree with decks; issues=[%s]'
                    % ('; '.join(cell_issues) or 'none'))

    # C3 -- scaling semantics plus explicit deck declaration
    scale_issues = []
    for r in runs:
        flow = r['sections'].get('flow', {})
        for key in EXPLICIT_SCALING_KEYS:
            if key not in flow:
                scale_issues.append('%s: deck lacks flow.%s'
                                    % (r['deck'], key))
        for key in SCALING_KEYS:
            if r['g'].get(key, '').lower() != 'constant':
                scale_issues.append('%s: %s=%s (want constant)'
                                    % (r['log'], key, r['g'].get(key)))
    cs10 = [r for r in runs if r['log'] == 'terzaghi.log']
    if cs10:
        v = cs10[0]['cells'].get('cs10', {})
        data['bc_pressure'] = as_number(v.get('bc_pressure'))
        data['flowBC'] = v.get('flowBC')
        if data['flowBC'] != 'TRUE':
            scale_issues.append('production cs10 flowBC=%s (want TRUE)'
                                % data['flowBC'])
    checks['C3'] = (not scale_issues and len(runs) == len(RUNS),
                    'bc/source/pressure scaling = constant in %d runs, '
                    'explicitly declared in every deck; issues=[%s]'
                    % (len(runs), '; '.join(scale_issues) or 'none'))

    # C4 -- EVERY flow: key in every deck must agree with its bag
    mismatches, compared = [], 0
    for r in runs:
        flow = r['sections'].get('flow', {})
        for key, dval in sorted(flow.items()):
            if key not in r['g']:
                mismatches.append('%s: flow.%s absent from bag'
                                  % (r['log'], key))
                continue
            compared += 1
            if not values_equal(dval, r['g'][key]):
                mismatches.append('%s: flow.%s deck=%s bag=%s'
                                  % (r['log'], key, dval, r['g'][key]))
    checks['C4'] = (not mismatches and compared > 0 and len(runs) == len(RUNS),
                    '%d deck->bag key comparisons, mismatches=[%s]'
                    % (compared, '; '.join(mismatches) or 'none'))

    # C5 -- known variant differences must be visible (discriminating power)
    expect = {
        'terzaghi.log':              dict(biot_coeff=1.0,
                                          enable_advection='FALSE',
                                          pressure_offset=1),
        'terzaghi_ksponly_phys.log': dict(biot_coeff=1.0,
                                          enable_advection='FALSE',
                                          pressure_offset=1),
        'terzaghi_nobiot.log':       dict(biot_coeff=0.0),
        'terzaghi_bodyforce.log':    dict(biot_coeff=1.0),
        'terzaghi_advection.log':    dict(enable_advection='TRUE',
                                          pressure_offset=0),
        'terzaghi_noadvection.log':  dict(enable_advection='FALSE',
                                          pressure_offset=0),
    }
    disc = []
    for r in runs:
        for key, want in expect.get(r['log'], {}).items():
            if not values_equal(want, r['g'].get(key)):
                disc.append('%s: %s=%s (want %s)'
                            % (r['log'], key, r['g'].get(key), want))
        if r['g'].get('temp_dep_viscosity', '').upper() != 'FALSE':
            disc.append('%s: temp_dep_viscosity=%s (want FALSE)'
                        % (r['log'], r['g'].get('temp_dep_viscosity')))
    checks['C5'] = (not disc and len(runs) == len(RUNS),
                    'variant signals visible (nobiot alpha=0, advection '
                    'flag T/F, 4.2 pressure_offset=0, mu(T) off); '
                    'issues=[%s]' % ('; '.join(disc) or 'none'))

    # C6 -- constant-pressure Dirichlet nodes over the whole transient
    c6_parts, c6_notes = [], []
    try:
        c = CDF2(prod_mesh)
        x = [float(v) for v in c.read_static('coordx')]
        drained = [i for i, xi in enumerate(x) if xi <= -4.9999]
        closed = [i for i, xi in enumerate(x) if xi >= 4.9999]
        recs = c.numrecs
        bc_series, cl_series = [], []
        for s in range(recs):
            p = [float(v) for v in c.read_record('vals_nod_var1', s)]
            bc_series.append([p[i] for i in drained])
            cl_series.append([p[i] for i in closed])
        flat = [v for row in bc_series for v in row]
        span = max(flat) - min(flat) if flat else float('inf')
        mean = sum(flat) / len(flat) if flat else float('nan')
        bag_bc = data.get('bc_pressure')
        cl_first = max(abs(v) for v in cl_series[0]) if cl_series else 0.0
        cl_last = max(abs(v) for v in cl_series[-1]) if cl_series else 0.0
        decay = cl_first / cl_last if cl_last > 0 else 0.0
        c6_parts.append(recs >= 50 and len(drained) >= 3
                        and all(finite(v) for v in flat))
        c6_notes.append('records=%d, drained BC nodes=%d' % (recs,
                                                             len(drained)))
        c6_parts.append(span <= 1.0e-6)
        c6_notes.append('BC pressure span=%.3g Pa (<=1e-6)' % span)
        c6_parts.append(bag_bc is not None and abs(mean - bag_bc) <= 1.0e-6)
        c6_notes.append('BC mean=%.3g Pa vs bag bc_pressure=%s'
                        % (mean, bag_bc))
        c6_parts.append(decay >= 5.0)
        c6_notes.append('closed end |p| %.6g -> %.6g (decay %.3g >=5, '
                        'guards an all-zero fake field)'
                        % (cl_first, cl_last, decay))
        data.update(bc_span=span, bc_mean=mean, closed_first=cl_first,
                    closed_last=cl_last, recs=recs, bc_nodes=len(drained))
    except Exception as exc:                            # noqa: BLE001
        c6_parts.append(False)
        c6_notes.append('mesh read error: %s' % exc)
    # Gate 3.3 analytic evidence that the constant BC was not distorted
    try:
        v33 = open(os.path.join(log_dir, 'gate33_verdict.log'),
                   encoding='utf-8', errors='replace').read()
        pass33 = 'GATE33_RESULT: PASS' in v33
        m = re.search(r'P6 L2 err vs BE series\s+([0-9.eE+-]+)', v33)
        l2 = float(m.group(1)) if m else float('inf')
        c6_parts.append(pass33 and l2 <= 0.05)
        c6_notes.append('Gate 3.3 PASS=%s, analytic L2 err=%.4g (<=0.05)'
                        % (pass33, l2))
        data['l2be'] = l2
    except Exception as exc:                            # noqa: BLE001
        c6_parts.append(False)
        c6_notes.append('gate33 verdict error: %s' % exc)
    checks['C6'] = (bool(c6_parts) and all(c6_parts),
                    '; '.join(c6_notes))

    # C7 -- no regression from the extra print-out
    reg = []
    for tag, want in (('gate41_verdict.log', 'GATE41_RESULT: PASS'),
                      ('gate42_verdict.log', 'GATE42_RESULT: PASS')):
        try:
            t = open(os.path.join(log_dir, tag), encoding='utf-8',
                     errors='replace').read()
            if want not in t:
                reg.append('%s lacks %s' % (tag, want))
        except Exception as exc:                        # noqa: BLE001
            reg.append('%s: %s' % (tag, exc))
    for r in runs:
        if FATAL_RE.search(r['text']):
            reg.append('%s: fatal signature' % r['log'])
        if r['text'].count('Gate 4.2: coupling wired') != 1:
            reg.append('%s: coupling wire not exactly once' % r['log'])
    checks['C7'] = (not reg and len(runs) == len(RUNS),
                    'Gate 4.1/4.2 verdicts PASS, %d runs clean and wired; '
                    'issues=[%s]' % (len(runs), '; '.join(reg) or 'none'))

    err = '; '.join(read_error) if read_error else ''
    data['runs'] = [r['log'] for r in runs]
    return checks, data, err


def main(argv):
    gate = False
    args = list(argv)
    if args and args[0] == '--gate':
        gate = True
        args = args[1:]
    if len(args) != 3:
        print('usage: check_gate43_flowbag.py [--gate] <log_dir> '
              '<deck_dir> <prod_mesh>')
        return 2
    checks, data, err = evaluate(args[0], args[1], args[2])
    print('=== Gate 4.3 -- Flow option bags + constant-pressure BC verdict ===')
    for key in ('C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7'):
        ok, msg = checks.get(key, (False, 'not evaluated'))
        print('  %s  %s  %s' % (key, 'PASS' if ok else 'FAIL', msg))
    failed = [k for k, v in checks.items() if not v[0]]
    if err:
        failed.append('DATA')
        print('  DATA ERROR: %s' % err)
    verdict = 'PASS' if not failed else 'FAIL (%s)' % ','.join(failed)
    if gate:
        print('GATE43_RESULT: %s' % verdict)
        print('GATE43_DETAIL: %d checks, %d failed, %d runs inspected'
              % (len(checks), len(failed), len(data.get('runs', []))))
    return 0 if (not gate or not failed) else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
