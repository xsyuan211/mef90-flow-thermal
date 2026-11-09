#!/usr/bin/env python3
"""Gate 4.8 -- discrete CONSERVATION post-processing for the Darcy solve.

Pure standard library + the committed ``parse_cdf2.py`` (the CI runner has no
numpy / scipy / PyYAML).  No Fortran changes: every quantity is reconstructed
from artefacts the solver already writes.

Modes
-----
  report : check_gate48_conservation.py <main.log> <source.log> <main.gen> <source.gen>
      INFORMATIONAL. Always exits 0, never a gate.

  gate   : ... --gate <main.log> <source.log> <main.gen> <source.gen>
      HARD Gate 4.8. Prints ``GATE48_RESULT: PASS|FAIL`` and exits 0 / 1.
      Missing / empty / short / unreadable / non-finite data is ALWAYS FAIL
      (CI #73 lesson: "no data" must never read as a pass).

  selftest : ... --selftest
      Offline algebraic replay of BOTH conservation laws (exact analytic
      first eigenmode + backward-Euler decay) and 7 negative controls.
      No PETSc, no Exodus, no CI.

Physics
-------
The Flow solve is the transient Darcy / storage equation

        -div(A grad p) + S dp/dt = f ,     A = k(v)/mu ,  S = 4.6e-10 1/Pa

advanced with one backward-Euler step per outer step (dt = 0.02 s).  With the
homogeneous drained Dirichlet face p = 0 at x = -5, two GLOBAL conservation
laws follow from the weak form.

(1) ENERGY / DISSIPATION balance  (test function w = p^{n+1}; the boundary
    term vanishes because p = 0 there):

        S M (p^{n+1}-p^n)/dt . p^{n+1}  +  int A |grad p^{n+1}|^2
                                       -  int f p^{n+1}  = 0
    =>  S1_n/dt + 2 E(n+1) - W(n+1) = 0

    E and W are EXACTLY the Gate 4.5 diagnostics (E = 1/2 int A|grad p|^2,
    W = int f p), so this also validates that 4.5 is wired to the same field.

(2) MASS / VOLUME balance  (test function w = 1; the Dirichlet reaction IS
    the net outflow):

        [ int S p^{n+1} - int S p^n ]/dt + Q_out(n+1) - F = 0
    =>  dM_n/dt + Q(n+1) - F = 0

    M = int S p dOmega is integrated from the Exodus nodal pressure, F =
    int f dOmega is the deck source, and Q = int_{x=-5} q.n dS is recovered
    INDEPENDENTLY from the solver's own per-cell Darcy velocity field
    (Exodus elem var 5 = cellDarcyVelocity_x, Gate 2 output), NOT from the
    equation itself.  U2 proves that field equals -A dp/dx to 1e-5, so this
    is not a circular argument.

Why the residual is not exactly zero
------------------------------------
CI #102 offline replay gives a stable |R|/scale ~ 0.5%-1.0% for BOTH laws,
and both imply the same S_eff/S ~ 0.985-0.989.  That is the FE/BE
discretisation error of the same order as the Gate 3.3 analytic L2 error
(0.032), not a sign or wiring bug.  The hard band is therefore 3%, which is
still violently falsifiable: an unwired E, a flipped source sign, a wrong S
or dt, or a missing outflow all drive |R|/scale to O(1).

Isolation
---------
terzaghi1D_source.yaml differs from terzaghi1D.yaml in EXACTLY ONE option
(cs0001.source 0.0 -> 1.0e-6, Gate 4.4).  Hence F_main = W_main = 0 while
F_src, W_src > 0, and both laws must close on both decks.

Criteria
--------
  U1  data integrity (51 Exodus steps + 51 Gate 4.5 rows per deck, finite,
      monotone time to t ~ 1.0, no crash signature)
  U2  drained face found at x = -5, boundary length ~ 1 m, and the Exodus
      Darcy velocity really equals -A dp/dx (rel <= 1e-5)
  U3  main deck energy balance   max |R_E|/scale <= 3%
  U4  source deck energy balance max |R_E|/scale <= 3% with W_src > 0
  U5  main deck mass balance     max |R_M|/scale <= 3% with F_main = 0
  U6  source deck mass balance   max |R_M|/scale <= 3% with F_src > 0
  U7  transient closure: outflow/dissipation decay, residual does not grow
  U8  cross-check W_src * S * Area == F_src * M_src (both are int f p / S)
      and W_main == 0
"""
import math
import os
import re
import sys

try:
    from parse_cdf2 import CDF2
except ImportError as exc:                                   # pragma: no cover
    CDF2 = None
    _CDF2_IMPORT_ERROR = exc

# --- Deck constants (terzaghi1D.yaml / terzaghi1D_source.yaml) -------------
S = 4.6e-10                    # storativity [1/Pa]
DT = 0.02                      # one BE step per outer step [s]
F_SOURCE = 1.0e-6              # cs0001.source in terzaghi1D_source.yaml
H = 10.0                       # drainage path [m]
WIDTH = 1.0                    # domain width in y [m]
MIN_RECORDS = 50

# --- Hard thresholds -----------------------------------------------------
RES_MAX = 0.03                 # |R|/scale for BOTH laws (CI #102: <= 0.010)
VEL_REL_MAX = 1.0e-5           # Exodus Darcy velocity vs -A dp/dx
BOUNDARY_LEN_LO, BOUNDARY_LEN_HI = 0.90, 1.10
CROSS_REL_MAX = 0.02           # W*S*Area vs F*M
TOL_T_END = 0.05
TAIL_GROWTH_MAX = 2.5          # late residual / early residual
DECAY_FRACTION_MAX = 0.25      # Q(49)/Q(1) upper bound

G45_RE = re.compile(
    r'Gate\s+4\.5:\s*step\s+(\d+)\s+t\s+([-+0-9.eEdD]+)\s+'
    r'E\s+([-+0-9.eEdD]+)\s+W\s+([-+0-9.eEdD]+)', re.I)
FATAL_RE = re.compile(
    r'Caught signal|SEGV|MPI_ABORT|PETSC ERROR|Unknown option|'
    r'attempt to deallocate|TSSolve failed', re.I)


def _num(s):
    """Parse a Fortran real (E or legacy D exponent)."""
    return float(s.replace('D', 'E').replace('d', 'e'))


def read_text(path):
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as fh:
            return fh.read()
    except Exception:                                  # noqa: BLE001 - gate fails
        return ''


def parse_rows(text):
    """Return [(step, t, E, W)] from the Gate 4.5 diagnostic rows."""
    rows = []
    for m in G45_RE.finditer(text):
        try:
            rows.append((int(m.group(1)), _num(m.group(2)),
                         _num(m.group(3)), _num(m.group(4))))
        except ValueError:
            continue
    return rows


def _tri(c):
    """Geometry of all TRI3 cells: (i,j,k,area,b[],cc[]) with 0-based nodes."""
    xs = [float(v) for v in c.read_static('coordx')]
    ys = [float(v) for v in c.read_static('coordy')]
    conn = c.read_static('connect1')
    out = []
    for q in range(len(conn) // 3):
        i, j, k = [int(v) - 1 for v in conn[q * 3:(q + 1) * 3]]
        two_a = ((xs[j] - xs[i]) * (ys[k] - ys[i])
                 - (xs[k] - xs[i]) * (ys[j] - ys[i]))
        if abs(two_a) <= 0.0:
            raise ValueError('degenerate TRI3 cell %d' % q)
        area = abs(two_a) / 2.0
        b = [(ys[j] - ys[k]) / two_a, (ys[k] - ys[i]) / two_a,
             (ys[i] - ys[j]) / two_a]
        cc = [(xs[k] - xs[j]) / two_a, (xs[i] - xs[k]) / two_a,
              (xs[j] - xs[i]) / two_a]
        out.append((i, j, k, area, b, cc))
    return xs, ys, out


def _lin(p, cell):
    """Exact integral of the P1 field p over one triangle."""
    i, j, k, area = cell[:4]
    return area * (p[i] + p[j] + p[k]) / 3.0


def _mquad(u, v, cell):
    """Consistent-mass bilinear form  u^T M v / S  over one triangle.

    M_ij = area/12 * (1 + delta_ij) in the P1 basis.
    """
    i, j, k, area = cell[:4]
    su, sv = u[i] + u[j] + u[k], v[i] + v[j] + v[k]
    return area / 12.0 * (su * sv + u[i] * v[i] + u[j] * v[j] + u[k] * v[k])


def build_terms(gen_path, log_text, f0):
    """Reconstruct every conservation term from an Exodus file + Gate 4.5 log."""
    if CDF2 is None:
        raise ValueError('cannot import parse_cdf2: %s' % _CDF2_IMPORT_ERROR)
    if not gen_path or not os.path.isfile(gen_path) or os.path.getsize(gen_path) == 0:
        raise ValueError('missing/empty mesh: %s' % gen_path)
    c = CDF2(gen_path)
    xs, ys, cells = _tri(c)
    nrec, nn = c.numrecs, len(xs)

    times, press, afield, vxfield = [], [], [], []
    for n in range(nrec):
        times.append(float(c.read_record('time_whole', n)[0]))
        press.append([float(v) for v in c.read_record('vals_nod_var1', n)])
        afield.append([float(v) for v in c.read_record('vals_elem_var4eb1', n)])
        vxfield.append([float(v) for v in c.read_record('vals_elem_var5eb1', n)])
    if len(cells) != len(afield[0]) or len(cells) != len(vxfield[0]):
        raise ValueError('cell-field size mismatch in %s' % gen_path)
    for p in press:
        if len(p) != nn or not all(math.isfinite(v) for v in p):
            raise ValueError('bad pressure field in %s' % gen_path)
    if not all(math.isfinite(t) for t in times):
        raise ValueError('bad time axis in %s' % gen_path)
    for arr in (afield, vxfield):
        if any(not all(math.isfinite(v) for v in row) for row in arr):
            raise ValueError('non-finite cell field in %s' % gen_path)

    # --- drained Dirichlet face: the p=0 nodes of the FIRST step ----------
    pmax = max(abs(v) for v in press[0]) or 1.0
    bc = set(i for i, v in enumerate(press[0]) if abs(v) < 1e-9 * pmax)
    if not bc:
        raise ValueError('no Dirichlet (p=0) node found in %s' % gen_path)
    bx = [xs[i] for i in bc]
    edge_x = sum(bx) / len(bx)
    x_spread = max(bx) - min(bx)
    bedge = []
    for q, cell in enumerate(cells):
        bn = [v for v in cell[:3] if v in bc]
        if len(bn) == 2:
            length = math.hypot(xs[bn[0]] - xs[bn[1]], ys[bn[0]] - ys[bn[1]])
            bedge.append((q, length))
    boundary_len = sum(L for _, L in bedge)

    # --- outflow from the solver's OWN Darcy velocity, cross-checked ------
    # n = (-1, 0) on x = -5, so q.n = -q_x.  dp/dx is reconstructed from the
    # P1 gradient and permField A, which is INDEPENDENT of the velocity Vec.
    vel_rel = 0.0
    qout = []
    for n in range(nrec):
        p, av, vx = press[n], afield[n], vxfield[n]
        tot = 0.0
        for q, L in bedge:
            i, j, k, area, b, cc = cells[q]
            dpdx = b[0] * p[i] + b[1] * p[j] + b[2] * p[k]
            qx_grad = -av[q] * dpdx
            denom = max(abs(vx[q]), abs(qx_grad), 1e-30)
            vel_rel = max(vel_rel, abs(vx[q] - qx_grad) / denom)
            tot += (-vx[q]) * L
        qout.append(tot)

    # --- storage / dissipation terms --------------------------------------
    mass, mass2, s1, area_tot = [], [], [], sum(cell[3] for cell in cells)
    for n in range(nrec):
        p = press[n]
        mass.append(S * sum(_lin(p, cell) for cell in cells))
        mass2.append(S * sum(_mquad(p, p, cell) for cell in cells))
        if n > 0:
            dp = [p[i] - press[n - 1][i] for i in range(nn)]
            s1.append(S * sum(_mquad(dp, p, cell) for cell in cells))

    rows = parse_rows(log_text)
    terms = dict(
        nrec=nrec, times=times, log_rows=len(rows),
        E=[r[2] for r in rows], W=[r[3] for r in rows],
        log_steps=[r[0] for r in rows], log_times=[r[1] for r in rows],
        M=mass, M2=mass2, S1=s1, Q=qout,
        F=f0 * area_tot, area=area_tot,
        boundary_len=boundary_len, edge_x=edge_x, x_spread=x_spread,
        n_bc_nodes=len(bc), n_bedge=len(bedge), vel_rel=vel_rel,
        fatal=bool(FATAL_RE.search(log_text)),
        x_min=min(xs), x_max=max(xs),
    )
    return terms


# --------------------------------------------------------------------------
# Verdict
# --------------------------------------------------------------------------
def _safe_rel(num, den):
    return abs(num) / den if den > 0 else float('inf')


def residuals(t):
    """(rE, rM, scaleE, scaleM) for n = 1 .. N-1."""
    rE, rM, sE, sM = [], [], [], []
    for n in range(1, min(t['nrec'], len(t['E']), len(t['M']))):
        dt = t['times'][n] - t['times'][n - 1]
        if dt <= 0:
            return None
        st = t['S1'][n - 1] / dt
        e2 = 2.0 * t['E'][n]
        w = t['W'][n]
        re_ = st + e2 - w
        dm = (t['M'][n] - t['M'][n - 1]) / dt
        q = t['Q'][n]
        f = t['F']
        rm = dm + q - f
        rE.append(re_)
        rM.append(rm)
        sE.append(abs(st) + abs(e2) + abs(w))
        sM.append(abs(dm) + abs(q) + abs(f))
    return rE, rM, sE, sM


def evaluate(main, src, err=None):
    """Return (checks, metrics, load_error); checks: name -> (ok, detail)."""
    chk, met = {}, {}
    if err:
        return chk, met, err
    for label, t in (('main', main), ('source', src)):
        if t is None:
            return chk, met, '%s deck data unavailable' % label

    # --- U1: data integrity ------------------------------------------------
    ok = True
    for label, t in (('main', main), ('source', src)):
        n = t['nrec']
        ok &= (n >= MIN_RECORDS and t['log_rows'] >= MIN_RECORDS
               and t['log_rows'] == n)
        ok &= (len(t['E']) == t['log_rows'] and len(t['W']) == t['log_rows']
               and len(t['M']) == n and len(t['Q']) == n
               and len(t['S1']) == n - 1)
        ok &= (t['log_steps'] == list(range(1, t['log_rows'] + 1)))
        ok &= all(math.isfinite(v) and v >= 0 for v in t['E'] + t['W'])
        ok &= all(math.isfinite(v) for v in t['M'] + t['Q'] + t['S1'])
        ok &= abs(t['times'][-1] - 1.0) <= TOL_T_END
        ok &= all(b > a for a, b in zip(t['times'], t['times'][1:]))
        ok &= not t['fatal']
    chk['U1'] = (ok, 'main %d steps/%d log rows, source %d/%d; '
                 'finite, monotone t->1.0, no crash signature'
                 % (main['nrec'], main['log_rows'],
                    src['nrec'], src['log_rows']))
    met.update(main_steps=main['nrec'], source_steps=src['nrec'])
    if not ok:
        return chk, met, None

    # --- U2: the outflow is a real, independently-reconstructed field -----
    b = main
    face_ok = (b['n_bc_nodes'] >= 2 and b['x_spread'] < 1e-8
               and abs(b['edge_x'] - b['x_min']) < 1e-8
               and b['n_bedge'] >= 2
               and BOUNDARY_LEN_LO <= b['boundary_len'] <= BOUNDARY_LEN_HI)
    chk['U2'] = (face_ok and b['vel_rel'] <= VEL_REL_MAX
                 and src['vel_rel'] <= VEL_REL_MAX,
                 'drained face x=%.3f (xmin=%.3f), %d nodes, %d edges, '
                 'L=%.4f m; max rel|vel - (-A dp/dx)| = %.3g (<= %.1e)'
                 % (b['edge_x'], b['x_min'], b['n_bc_nodes'], b['n_bedge'],
                    b['boundary_len'], max(b['vel_rel'], src['vel_rel']),
                    VEL_REL_MAX))
    met['vel_rel'] = max(b['vel_rel'], src['vel_rel'])

    # --- U3..U6: both laws on both decks -----------------------------------
    for label, t, fzero, wzero in (
            ('U3', main, True, True), ('U4', src, False, False),
            ('U5', main, True, True), ('U6', src, False, False)):
        r = residuals(t)
        if r is None:
            chk[label] = (False, 'non-positive dt / short series')
            continue
        rE, rM, sE, sM = r
        if label in ('U3', 'U4'):
            rel = [_safe_rel(a, b) for a, b in zip(rE, sE)]
            extra = (max(t['W']) > 0) if label == 'U4' else (max(t['W']) == 0)
            what = 'energy   S1/dt + 2E - W'
        else:
            rel = [_safe_rel(a, b) for a, b in zip(rM, sM)]
            extra = (t['F'] > 0) if label == 'U6' else (t['F'] == 0)
            what = 'mass     dM/dt + Q - F'
        mx = max(rel)
        met['%s_max_rel' % label] = mx
        chk[label] = (mx <= RES_MAX and extra,
                      '%s: max|R|/scale=%.5f (<= %.2f), F=%.3g, '
                      'maxW=%.4g, source/frozen-term ok=%s'
                      % (what, mx, RES_MAX, t['F'], max(t['W']), extra))
        met['%s_res' % label] = (rE if label in ('U3', 'U4') else rM)
        met['%s_scale' % label] = (sE if label in ('U3', 'U4') else sM)

    # --- U7: transient closure (no growing residual, real decay) ----------
    q = main['Q'][2:]
    decay_ok = (all(b <= a * 1.02 for a, b in zip(q, q[1:]))
                and q[-1] <= DECAY_FRACTION_MAX * q[0])
    growth = []
    for lab in ('U3', 'U4', 'U5', 'U6'):
        res = [abs(a) / b for a, b in zip(met['%s_res' % lab],
                                          met['%s_scale' % lab])]
        head = sum(res[:10]) / 10.0
        tail = sum(res[-10:]) / 10.0
        # Machine-precision closure (the analytic selftest) is a PASS, not a
        # 0/0 "growth" artefact.  Only compare a real, finite residual.
        growth.append(1.0 if max(head, tail) < 1e-9
                      else (tail / head if head > 0 else float('inf')))
    grow_ok = max(growth) <= TAIL_GROWTH_MAX
    chk['U7'] = (decay_ok and grow_ok,
                 'Q decay %.4g -> %.4g (ratio %.3f <= %.2f), '
                 'late/early residual growth %s (<= %.1f)'
                 % (q[0], q[-1], q[-1] / q[0] if q[0] else 0,
                    DECAY_FRACTION_MAX,
                    ','.join('%.2f' % g for g in growth), TAIL_GROWTH_MAX))
    met['Q_first'], met['Q_last'] = q[0], q[-1]

    # --- U8: W = f*int p, M = S*int p  =>  W*S*Area == F*M --------------
    cross = []
    for n in range(1, src['nrec']):
        lhs = src['W'][n] * S * src['area']
        rhs = src['F'] * src['M'][n]
        den = max(abs(rhs), 1e-30)
        cross.append(abs(lhs - rhs) / den)
    cross_max = max(cross)
    frozen = max(abs(v) for v in main['W']) == 0 and main['F'] == 0
    chk['U8'] = (cross_max <= CROSS_REL_MAX and frozen and src['F'] > 0,
                 'max|W*S*Area/(F*M) - 1| = %.5f (<= %.2f); '
                 'W_main==0 & F_main==0 = %s, F_src=%.3g'
                 % (cross_max, CROSS_REL_MAX, frozen, src['F']))
    met['cross_rel'] = cross_max
    return chk, met, None


def print_report(chk, met, err):
    print('=== Gate 4.8 -- Darcy discrete conservation verdict ===')
    if err:
        print('  DATA ERROR: %s' % err)
    for name in sorted(chk):
        ok, detail = chk[name]
        print('  %-4s %-4s %s' % (name, 'PASS' if ok else 'FAIL', detail))
    if met:
        print('  --- metrics ---')
        for key in ('main_steps', 'source_steps', 'vel_rel', 'Q_first',
                    'Q_last', 'cross_rel', 'U3_max_rel', 'U4_max_rel',
                    'U5_max_rel', 'U6_max_rel'):
            if key in met:
                print('    %-20s = %s' % (key, met[key]))
    return [k for k, v in chk.items() if not v[0]]


# --------------------------------------------------------------------------
# Offline selftest: exact analytic replay of both laws + negative controls
# --------------------------------------------------------------------------
def analytic_terms(f0, mode='ok'):
    """Algebraic replay of the CONSERVATION LAWS (not of the FE field).

    First eigenmode p(s,n) = a_n cos(lambda s), lambda = pi/(2H), with the
    EXACT backward-Euler decay a_n = a_{n-1}/(1 + A lambda^2 dt / S).
    Because K phi = A lambda^2 phi and M phi = S phi, both discrete laws close
    to machine precision.  Negative controls break exactly one ingredient.
    """
    A = 5.0e-8
    lam = math.pi / (2.0 * H)
    beta = A * lam * lam * DT / S
    rho = 1.0 / (1.0 + beta)
    i1 = 1.0 / lam                     # int_0^H cos(lambda s) ds
    i2 = H / 2.0                       # int_0^H cos^2(lambda s) ds
    area = WIDTH * H
    nrec = 51
    a = [1.0e5 * rho ** n for n in range(nrec)]
    F = f0 * area
    E, W, M, M2, S1, Q = [], [], [], [], [], []
    for n in range(nrec):
        E.append(0.5 * A * WIDTH * a[n] ** 2 * lam ** 2 * i2)
        W.append(f0 * WIDTH * a[n] * i1)
        M.append(S * WIDTH * a[n] * i1)
        M2.append(S * WIDTH * a[n] ** 2 * i2)
        Q.append(A * WIDTH * a[n] * lam)
        if n > 0:
            S1.append(S * WIDTH * (a[n] - a[n - 1]) * a[n] * i2)
    # The source deck adds F and W; keep both laws exactly closed by defining
    # the source-work-consistent dissipation (algebraic replay, not a field).
    if f0 > 0:
        for n in range(1, nrec):
            E[n] = 0.5 * (-S1[n - 1] / DT + W[n])
        for n in range(1, nrec):
            Q[n] = -(M[n] - M[n - 1]) / DT + F

    if mode == 'E_zero':
        E = [0.0] * nrec
    elif mode == 'W_zero':
        W = [0.0] * nrec
    elif mode == 'Q_zero':
        Q = [0.0] * nrec
    elif mode == 'S_wrong':
        M = [10.0 * v for v in M]
        M2 = [10.0 * v for v in M2]
        S1 = [10.0 * v for v in S1]
    elif mode == 'F_sign':
        F = -F
        W = [-v for v in W]
    elif mode == 'short':
        nrec = 12
        E, W, M, Q = E[:nrec], W[:nrec], M[:nrec], Q[:nrec]
        S1 = S1[:nrec - 1]
    return dict(
        nrec=nrec, times=[n * DT for n in range(nrec)], log_rows=nrec,
        E=E, W=W, log_steps=list(range(1, nrec + 1)),
        log_times=[n * DT for n in range(nrec)],
        M=M, M2=M2, S1=S1, Q=Q, F=F, area=area,
        boundary_len=WIDTH, edge_x=-H / 2.0, x_spread=0.0,
        n_bc_nodes=6, n_bedge=5, vel_rel=0.0, fatal=False,
        x_min=-H / 2.0, x_max=H / 2.0)


def selftest():
    cases = [
        ('positive (both laws closed)', 'ok', 'ok', True),
        ('E all zero (Gate 4.5 unwired)', 'E_zero', 'E_zero', False),
        ('W all zero (sourceSec unwired)', 'W_zero', 'W_zero', False),
        ('Q all zero (outflow unwired)', 'Q_zero', 'Q_zero', False),
        ('S off by 10x', 'S_wrong', 'S_wrong', False),
        ('source sign flipped (G+F)', 'ok', 'F_sign', False),
        ('short series (missing rows)', 'short', 'short', False),
    ]
    print('=== Gate 4.8 offline selftest (algebraic replay, no solver) ===')
    all_ok = True
    for label, mm, sm, want in cases:
        # main deck is f=0; but analytic_terms(0, 'E_zero') etc. also apply
        main = analytic_terms(0.0, mm)
        src = analytic_terms(F_SOURCE, sm)
        chk, met, err = evaluate(main, src)
        fails = [k for k, v in chk.items() if not v[0]]
        got = (not fails) and not err
        ok = got == want
        all_ok &= ok
        print('  %-38s -> %-4s (want %-4s) %s  failing=%s'
              % (label, 'PASS' if got else 'FAIL', 'PASS' if want else 'FAIL',
                 'OK' if ok else 'SELFTEST ERROR', ','.join(fails) or 'none'))
        if label.startswith('positive'):
            print('      max|R|/scale: U3=%.3g U4=%.3g U5=%.3g U6=%.3g'
                  % tuple(met.get('%s_max_rel' % u, float('nan'))
                          for u in ('U3', 'U4', 'U5', 'U6')))
    print('SELFTEST_RESULT: %s' % ('PASS' if all_ok else 'FAIL'))
    return 0 if all_ok else 1


def main(argv):
    if argv and argv[0] == '--selftest':
        return selftest()
    gate = bool(argv) and argv[0] == '--gate'
    args = argv[1:] if gate else argv
    usage = '<main.log> <source.log> <main.gen> <source.gen>'
    if len(args) != 4:
        if gate:
            print('GATE48_RESULT: FAIL (usage: --gate %s)' % usage)
            return 1
        print('usage: check_gate48_conservation.py [--gate|--selftest] %s' % usage)
        return 0
    main_log, src_log, main_gen, src_gen = args
    err = None
    main_t = src_t = None
    try:
        mt, st = read_text(main_log), read_text(src_log)
        if not mt or not st:
            raise ValueError('production or source log unreadable/empty')
        main_t = build_terms(main_gen, mt, 0.0)
        src_t = build_terms(src_gen, st, F_SOURCE)
    except Exception as exc:                               # noqa: BLE001
        err = 'cannot build conservation terms: %r' % (exc,)
    chk, met, err = evaluate(main_t, src_t, err)
    fails = print_report(chk, met, err)
    if gate:
        if err or fails or len(chk) < 8:
            print('GATE48_RESULT: FAIL (%s)'
                  % (err if err else ','.join(fails) or 'incomplete criteria'))
            return 1
        print('GATE48_RESULT: PASS (U1 data, U2 outflow field, '
              'U3/U4 energy balance <=%.3g/<=%.3g, '
              'U5/U6 mass balance <=%.3g/<=%.3g, U7 transient closure, '
              'U8 W*S*Area==F*M rel=%.4g)'
              % (met['U3_max_rel'], met['U4_max_rel'],
                 met['U5_max_rel'], met['U6_max_rel'], met['cross_rel']))
        return 0
    print('VERDICT: %s (report mode, informational only)'
          % ('FAIL ' + ','.join(fails) if fails else 'PASS'))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
