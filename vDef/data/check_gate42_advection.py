#!/usr/bin/env python3
"""Gate 4.2 physics verdict: does cell-wise v_darcy reach the heat field?

Two identical Terzaghi runs (terzaghi1D_advection.yaml /
terzaghi1D_noadvection.yaml) differ in EXACTLY one key,
flow.enable_advection. Both solve the same steady heat problem with the same
Dirichlet temperature BCs. The ONLY possible source of a temperature-field
difference is the v_darcy.grad(T) term assembled by
m_HeatXfer/m_MEF90_HeatXferAssembly.F90.

Verdict deliberately does NOT trust solver messages, fnorm placeholders,
SNESConvergedReason or failure counts (four previous false convergences).
It compares PHYSICAL FIELDS in two Exodus meshes.

Usage:
  check_gate42_advection.py --gate <adv.log> <noadv.log> <adv.gen> <noadv.gen> <main.gen>
      -> prints "GATE42_RESULT: PASS|FAIL", exits 0 / 1 (hard gate)
  check_gate42_advection.py <adv.log> <noadv.log> <adv.gen> <noadv.gen> <main.gen>
      -> report only, always exits 0

Pure standard library + the repo's vDef/data/parse_cdf2.py (the CI runner
has no numpy). Zero compilation, zero solver execution.
"""
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from parse_cdf2 import CDF2
except Exception as exc:                                # pragma: no cover
    print('GATE42_RESULT: FAIL (cannot import parse_cdf2: %s)' % exc)
    sys.exit(1)

NUM = r'[+-]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][+-]?\d+)?'
ROW_RE = re.compile(
    r'Gate\s+4\.2:\s*enabled=\s*([A-Za-z.]+).*?\|v_darcy\|inf=\s*(' + NUM + ')',
    re.I)
WIRE_RE = re.compile(r'Gate\s+4\.2:\s*coupling wired', re.I)
# CI #90: the EXO nodal-variable name table is fixed by MEF90DefMechFormatEXO
# (Pressure, Damage, Displacement_X, Displacement_Y) and CANNOT name the
# temperature field. Lock the slot from the solver's own bag view instead.
TEMP_OFFSET_RE = re.compile(r'(?<![\w])temp_Offset\s*=\s*(\d+)', re.I)
FATAL_RE = re.compile(
    r'segmentation fault|signal 11|Caught signal|MPI_ABORT|PETSC ERROR|'
    r'attempt to deallocate|error while loading shared libraries|'
    r'TSSolve failed|Unknown option|Unrecognized option', re.I)

MIN_RECORDS = 50
MIN_PROBES = 50
REL_DIFF_MIN = 1.0e-4      # |T_on - T_off| / max|T_off| must exceed this
SPATIAL_MIN = 0.05         # de-meaned part of the difference, fraction of max
TEMP_LO, TEMP_HI = 1.0, 1.0e4
MIN_ABS_DIFF = 1.0e-4      # K: reject round-off masquerading as advection
MAX_DISP = 1.0e-3          # m: production |u| ~ 1e-5
MAX_DAMAGE = 0.1           # production damage ~ 4e-5


def finite(x):
    return isinstance(x, float) and math.isfinite(x)


def num(text):
    """Parse a Fortran real, including D exponents."""
    try:
        return float(text.strip().replace('D', 'E').replace('d', 'e'))
    except (ValueError, AttributeError):
        return float('nan')


def truthy(tok):
    t = tok.strip().upper().lstrip('.')
    return t.startswith('T') or t.startswith('1')


def read_log(path):
    """Return (rows, probes, wired, fatal_hit, error)."""
    if not os.path.isfile(path):
        return [], 0, 0, '', 'missing solver log: %s' % path
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as fh:
            lines = fh.read().splitlines()
    except Exception as exc:
        return [], 0, 0, '', 'unreadable solver log: %s' % exc
    rows, probes, wired, fatal = [], 0, 0, ''
    for line in lines:
        if 'Gate 3.2: step' in line:
            probes += 1
        if WIRE_RE.search(line):
            wired += 1
        m = ROW_RE.search(line)
        if m:
            rows.append((truthy(m.group(1)), num(m.group(2))))
        f = FATAL_RE.search(line)
        if f and not fatal and 'echo "' not in line and not line.lstrip().startswith('#'):
            fatal = f.group(0)
    return rows, probes, wired, fatal, ''


def read_mesh(path):
    """Return (CDF2, record_var_names, error)."""
    if not os.path.isfile(path):
        return None, [], 'missing result mesh: %s' % path
    try:
        c = CDF2(path)
        names = [v['name'] for v in c.vars
                 if v.get('is_record') and v['name'] != 'time_whole']
        if not names:
            return None, [], 'no record variables in %s' % path
        return c, names, ''
    except Exception as exc:
        return None, [], 'unreadable mesh %s: %s' % (path, exc)


def series(c, name, step):
    return [float(v) for v in c.read_record(name, step)]


def temp_slot_from_log(path):
    """Lock the temperature EXO slot from the solver's temp_Offset bag value.

    Both 4.2 diagnostic decks set temp_Offset: 1 and flow.pressure_offset: 0,
    so HeatXferViewEXO owns nodal slot 1. Auto-discovering the slot by largest
    |T_on - T_off| is UNSAFE (CI #90: the diverging Displacement_Y field had
    a 3e18 difference and was mistaken for the temperature response).
    """
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as fh:
            for line in fh:
                m = TEMP_OFFSET_RE.search(line)
                if m:
                    return 'vals_nod_var%d' % int(m.group(1))
    except Exception:
        pass
    return None


def evaluate(adv_log, noadv_log, adv_gen, noadv_gen, main_gen):
    checks, metrics, data_error = {}, {}, ''
    finite_ok = False             # defined even when a mesh is missing

    arows, aprobe, awire, afatal, aerr = read_log(adv_log)
    nrows, nprobe, nwire, nfatal, nerr = read_log(noadv_log)
    if aerr or nerr:
        data_error = aerr or nerr
    A, anames, ea = read_mesh(adv_gen)
    N, nnames, en = read_mesh(noadv_gen)
    M, mnames, em = read_mesh(main_gen)
    if ea or en or em:
        data_error = data_error or ea or en or em

    metrics['adv_rows'] = len(arows)
    metrics['noadv_rows'] = len(nrows)
    metrics['adv_probes'] = aprobe
    metrics['noadv_probes'] = nprobe

    # D1 -- both diagnostic runs completed with no fatal signature
    ok = (not data_error and aprobe >= MIN_PROBES and nprobe >= MIN_PROBES
          and len(arows) >= MIN_RECORDS and len(nrows) >= MIN_RECORDS
          and awire >= 1 and nwire >= 1 and not afatal and not nfatal)
    checks['D1'] = (ok, 'adv %d rows/%d probes, noadv %d rows/%d probes, wired %d/%d, fatal=%r'
                    % (len(arows), aprobe, len(nrows), nprobe, awire, nwire,
                       afatal or nfatal or None))
    # CI #89: when either variant crashed early its result mesh is whatever
    # the PREVIOUS variant left behind (the workflow reuses one EXO file).
    # D4/D5/D6 then "PASS" on stale pressure fields while D6 sees a zero
    # difference. Do not let mesh-based checks outrank a failed D1.
    d1_ok = ok

    # D2 -- the enable flag really differs between the two decks
    if arows and nrows:
        aflags = [r[0] for r in arows]
        nflags = [r[0] for r in nrows]
        ok2 = all(aflags) and not any(nflags)
        metrics['flags'] = (sum(aflags), len(aflags), sum(nflags), len(nflags))
        checks['D2'] = (ok2, 'adv enabled T on %d/%d rows, noadv T on %d/%d rows'
                        % (sum(aflags), len(aflags), sum(nflags), len(nflags)))
    else:
        checks['D2'] = (False, 'no Gate 4.2 enabled= rows')

    # D3 -- velocity field is alive, decays, and is IDENTICAL in both runs
    av = [r[1] for r in arows]
    nv = [r[1] for r in nrows]
    if av and nv and len(av) == len(nv):
        allfinite = all(finite(v) for v in av + nv)
        peak = max(av)
        first, last = av[0], av[-1]
        same = all(finite(a) and finite(b) and abs(a - b) <= 1e-9 * max(1.0, abs(a))
                   for a, b in zip(av, nv))
        decaying = finite(first) and finite(last) and first > last > 0.0
        metrics.update(vpeak=peak, vfirst=first, vlast=last, vsame=same)
        checks['D3'] = (allfinite and peak > 0 and decaying and same,
                        '|v|inf %.6g -> %.6g (peak %.6g), same in both logs=%s'
                        % (first, last, peak, same))
    else:
        checks['D3'] = (False, 'velocity rows missing / length mismatch')

    # D4 -- mesh integrity (records, node count, finite values)
    if A is not None and N is not None and M is not None:
        recs_ok = (A.numrecs >= MIN_RECORDS and N.numrecs >= MIN_RECORDS
                   and M.numrecs >= MIN_RECORDS and A.numrecs == N.numrecs)
        shared = [k for k in anames if k in nnames]
        node_ok, finite_ok, nd = True, True, 0
        for k in shared:
            for s in range(min(A.numrecs, N.numrecs)):
                try:
                    a, b = series(A, k, s), series(N, k, s)
                except Exception:
                    finite_ok = False
                    break
                if len(a) != len(b):
                    node_ok = False
                nd = max(nd, len(a))
                if not all(finite(x) for x in a + b):
                    finite_ok = False
        metrics['records'] = (A.numrecs, N.numrecs, M.numrecs)
        metrics['nodes'] = nd
        metrics['record_vars'] = len(shared)
        checks['D4'] = (d1_ok and recs_ok and node_ok and finite_ok and bool(shared),
                        'records %s, nodes=%d, shared slots=%d, all finite=%s, runs_complete=%s'
                        % (str(metrics['records']), nd, len(shared), finite_ok, d1_ok))
    else:
        checks['D4'] = (False, 'mesh unavailable')

    # D5/D6 -- LOCK the temperature slot (CI #90); the per-slot difference
    # table below is kept for diagnostics only, never for the verdict.
    if A is not None and N is not None:
        shared = [k for k in anames if k in nnames]
        diff_by_slot, detail = {}, {}
        for k in shared:
            mx, first_d, last_d, first_span, last_span = 0.0, 0.0, 0.0, 0.0, 0.0
            d = []
            for s in range(min(A.numrecs, N.numrecs)):
                try:
                    a, b = series(A, k, s), series(N, k, s)
                except Exception:
                    continue
                if len(a) != len(b):
                    continue
                d = [u - v for u, v in zip(a, b)]
                m = max(abs(x) for x in d) if d else 0.0
                mx = max(mx, m)
                if s == 0:
                    first_d = m
                    first_span = max(b) - min(b) if b else 0.0
                if s == min(A.numrecs, N.numrecs) - 1:
                    last_d = m
                    last_span = max(b) - min(b) if b else 0.0
                    mean = sum(d) / len(d) if d else 0.0
                    spatial = max(abs(x - mean) for x in d) if d else 0.0
            diff_by_slot[k] = mx
            detail[k] = (first_d, last_d, last_span, spatial, max(abs(x) for x in d) if d else 0.0)
        metrics['slot_diff'] = diff_by_slot
        auto_slot = max(diff_by_slot, key=lambda k: diff_by_slot[k]) if diff_by_slot else None
        metrics['auto_max_slot'] = auto_slot
        # CI #90: auto-discovery selected vals_nod_var4 (Displacement_Y,
        # |d| ~ 3e18) because enabling HeatXfer with the default
        # LinearThermalExpansion=I made the mechanical field diverge. The EXO
        # name table cannot identify temperature, so trust the deck's
        # temp_Offset bag value echoed by the solver instead.
        temp_slot = temp_slot_from_log(adv_log)
        metrics['temp_slot'] = temp_slot
        slot = temp_slot if temp_slot in shared else None
        metrics['response_slot'] = slot
        if slot:
            f_d, l_d, span, spatial, fmax = detail[slot]
            # reference field magnitude in the NOADVECTION mesh (last step)
            try:
                ref = series(N, slot, N.numrecs - 1)
                refmax = max(abs(x) for x in ref)
                ref_finite = bool(ref) and all(finite(x) for x in ref)
                first_ref = max(abs(x) for x in series(N, slot, 0))
            except Exception:
                refmax, ref_finite, first_ref = 0.0, False, 0.0
            grows = ref_finite and finite(first_ref) and refmax > first_ref + 1.0
            metrics.update(diff_first=f_d, diff_last=l_d, temp_span=span,
                           spatial=spatial, refmax=refmax, first_ref=first_ref)
            rel = l_d / refmax if refmax > 0 else 0.0
            spatial_frac = spatial / l_d if l_d > 0 else 0.0
            checks['D5'] = (d1_ok and ref_finite and grows
                            and TEMP_LO <= refmax <= TEMP_HI and span > 1.0,
                            'slot %s (temp_Offset): max|T_off|=%.6g K, first=%.6g K, '
                            'span=%.6g K, grows=%s, finite=%s, runs_complete=%s'
                            % (slot, refmax, first_ref, span, grows, ref_finite, d1_ok))
            checks['D6'] = (d1_ok
                            and all(finite(x) for x in (f_d, l_d, rel, spatial_frac))
                            and l_d >= MIN_ABS_DIFF
                            and rel >= REL_DIFF_MIN
                            and spatial_frac >= SPATIAL_MIN
                            and l_d > f_d,
                            'slot %s: |dT| first=%.6g last=%.6g K (min %.3g), '
                            'rel=%.4g (min %.3g), spatial=%.3g of diff'
                            % (slot, f_d, l_d, MIN_ABS_DIFF, rel, REL_DIFF_MIN,
                               spatial_frac))
        else:
            checks['D5'] = (False, 'temperature slot %r missing from result mesh '
                            '(%d shared slots)' % (temp_slot, len(shared)))
            checks['D6'] = (False, 'temperature slot %r missing' % temp_slot)
    else:
        checks['D5'] = (False, 'mesh unavailable')
        checks['D6'] = (False, 'mesh unavailable')

    # D7 -- production deck is isolated (its slot 1 is still the pressure
    # field, ~1e5 Pa, NOT a ~100 K temperature field)
    if M is not None and 'vals_nod_var1' in mnames:
        try:
            p0 = series(M, 'vals_nod_var1', 0)
            pmax = max(abs(x) for x in p0)
        except Exception:
            pmax = float('nan')
        metrics['prod_slot1_max'] = pmax
        checks['D7'] = (finite(pmax) and pmax > 1.0e4,
                        'production slot1 |max|=%.6g Pa (heatxfer stays null)' % pmax)
    else:
        checks['D7'] = (False, 'production slot1 unavailable')

    # D8 -- mechanical fields must stay at production magnitude.
    # CI #90: with HeatXfer enabled and the default alpha=1/K, thermal strain
    # saturated damage to 1 and |u| diverged to 1e18 m; D4 still saw only
    # FINITE values and auto-discovery mistook that displacement for the
    # temperature response. A valid heat-advection diagnostic must NOT blow
    # up the poroelastic solve.
    if A is not None and N is not None:
        worst_disp, worst_dmg = 0.0, 0.0
        for c in (A, N):
            for k in ('vals_nod_var3', 'vals_nod_var4'):
                try:
                    v = max(abs(x) for x in series(c, k, c.numrecs - 1))
                except Exception:
                    v = float('inf')
                worst_disp = max(worst_disp, v)
            try:
                v = max(abs(x) for x in series(c, 'vals_nod_var2', c.numrecs - 1))
            except Exception:
                v = float('inf')
            worst_dmg = max(worst_dmg, v)
        metrics['defmech_health'] = (worst_disp, worst_dmg)
        checks['D8'] = (d1_ok and finite(worst_disp) and finite(worst_dmg)
                        and worst_disp <= MAX_DISP and worst_dmg <= MAX_DAMAGE,
                        'max|u|=%.3g m (<=%.3g), max damage=%.3g (<=%.3g), '
                        'runs_complete=%s'
                        % (worst_disp, MAX_DISP, worst_dmg, MAX_DAMAGE, d1_ok))
    else:
        checks['D8'] = (False, 'mesh unavailable')

    return checks, metrics, data_error


def main(argv):
    gate = bool(argv) and argv[0] == '--gate'
    args = argv[1:] if gate else argv
    usage = ('[--gate] <adv.log> <noadv.log> <adv.gen> <noadv.gen> <main.gen>')
    if len(args) != 5:
        if gate:
            print('GATE42_RESULT: FAIL (usage: --gate %s)' % usage)
            return 1
        print('usage: check_gate42_advection.py %s' % usage)
        return 0

    checks, metrics, data_error = evaluate(*args)
    print('=== Gate 4.2 -- cell-wise Darcy velocity -> HeatXfer verdict ===')
    if data_error:
        print('  DATA ERROR: %s' % data_error)
    for key in ('D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8'):
        ok, msg = checks.get(key, (False, 'not evaluated'))
        print('  %s  %s  %s' % (key, 'PASS' if ok else 'FAIL', msg))
    slot = metrics.get('response_slot')
    auto_slot = metrics.get('auto_max_slot')
    if metrics.get('slot_diff'):
        print('  --- per-slot max |field_on - field_off| (diagnostic only) ---')
        for k, v in sorted(metrics['slot_diff'].items(), key=lambda kv: -kv[1]):
            mark = ''
            if k == slot:
                mark = ' <== temperature slot (verdict)'
            elif k == auto_slot:
                mark = ' <== largest raw diff (NOT used)'
            print('    %-20s %12.6g%s' % (k, v, mark))
    failed = [k for k, v in checks.items() if not v[0]]
    if data_error:
        failed.append('DATA')
    verdict = 'PASS' if not failed else 'FAIL (%s)' % ','.join(failed)
    if gate:
        print('GATE42_RESULT: %s' % verdict)
        print('GATE42_DETAIL: %d checks, %d failed, temp slot=%s (from temp_Offset), '
              'largest raw diff slot=%s'
              % (len(checks), len(failed), slot, auto_slot))
        return 0 if not failed else 1
    print('  VERDICT: %s' % verdict)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
