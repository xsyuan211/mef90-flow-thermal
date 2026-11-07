"""Post-process the Gate 3.2 Terzaghi result mesh.

Reads vDef/data/Beam1x10-tri3_out.gen (Exodus / NetCDF CDF2), extracts the
nodal pressure history, compares it pointwise with the 1D Terzaghi series
(continuum and backward-Euler-corrected), and writes CSV tables + a PNG chart.

Pure numpy + PIL (no scipy / matplotlib, matching the local environment).
"""
import csv
import math
import os
import struct
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
from parse_cdf2 import CDF2

from PIL import Image, ImageDraw, ImageFont

OUT = _HERE
GEN = os.path.join(OUT, 'Beam1x10-tri3_out.gen')

# ---- Physics / discretisation (terzaghi1D.yaml, CI #66) -------------------
P0 = 1.0e5                 # initial pressure [Pa]
K_OVER_MU = 5.0e-11 / 1.0e-3   # perm_min / viscosity [m^2/(Pa s)]
S = 4.6e-10                # storativity [1/Pa]
CV = K_OVER_MU / S         # consolidation coefficient [m^2/s]
H = 10.0                   # drainage path length [m] (x in [-5,5])
DT = 0.02                  # fixed BE step [s]
NSERIES = 400              # Terzaghi series modes


_MODES = np.arange(NSERIES, dtype=float)
_N = 2.0 * _MODES + 1.0
_LAM = CV * (_N * np.pi / (2.0 * H)) ** 2
_AMP = (4.0 / np.pi) * ((-1.0) ** _MODES) / _N
_COEF = None            # (nn, NSERIES) spatial modes, filled by setup()
_Z = None


def setup(z):
    """Precompute the z-dependent part of the Terzaghi series."""
    global _COEF, _Z
    _Z = np.asarray(z, dtype=float)
    _COEF = _AMP[None, :] * np.cos(np.outer(_Z, _N * np.pi / (2.0 * H)))


def _decay(t, be):
    if be:
        return (1.0 + _LAM * DT) ** (-round(t / DT))
    return np.exp(-_LAM * t)


def analytic_vec(t, be=False):
    """Vectorised p(z_i,t) for all nodes."""
    if t <= 0.0:
        return P0 * np.where(_Z >= H - 1e-12, 0.0, 1.0)
    return P0 * (_COEF @ _decay(t, be))


def analytic_vec_steps(nsteps):
    """BE-series field after exactly nsteps steps (independent of time label).

    nsteps=0 is the initial field WITH the drained Dirichlet BC applied
    (interior p0, drained face 0). The vDef output at t=0 is still the
    raw p0 field without BC, so its fields lag this series by one step.
    """
    if nsteps <= 0:
        return P0 * np.where(_Z >= H - 1e-12, 0.0, 1.0)
    return P0 * (_COEF @ ((1.0 + _LAM * DT) ** (-nsteps)))


def analytic(z, t, be=False):
    """Single-point p(z,t) [Pa]; z=0 closed (x=5), z=H drained (x=-5)."""
    if t <= 0.0:
        return 0.0 if z >= H - 1e-12 else P0
    spatial = _AMP * np.cos(_N * np.pi * z / (2.0 * H))
    return float(P0 * np.dot(spatial, _decay(t, be)))


def l2_rel(a, b):
    num = math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))
    den = math.sqrt(sum(y * y for y in b))
    return num / den if den > 0 else 0.0


def main():
    c = CDF2(GEN)
    x = list(c.read_static('coordx'))
    y = list(c.read_static('coordy'))
    nn = len(x)
    nrec = c.numrecs
    print('nodes %d  records %d  x[%.3f,%.3f]  y[%.3f,%.3f]' %
          (nn, nrec, min(x), max(x), min(y), max(y)))
    names = c.read_static('name_nod_var')
    print('nodal var names (raw chars):', repr(b''.join(names[:33])))
    print('cv = %.6f m^2/s, H = %.1f m, dt = %.3f s' % (CV, H, DT))

    z = [5.0 - xi for xi in x]
    setup(z)
    times, fields = [], []
    for step in range(nrec):
        t = c.read_record('time_whole', step)[0]
        p = list(c.read_record('vals_nod_var1', step))
        times.append(t)
        fields.append(p)

    # ---- CSV 1: compact time series --------------------------------------
    ts_path = os.path.join(OUT, 'terzaghi_pressure_timeseries.csv')
    ts_rows = []
    for step, (t, p) in enumerate(zip(times, fields), start=1):
        ac = analytic_vec(t)
        ab = analytic_vec(t, be=True)
        ab1 = analytic_vec_steps(max(round(t / DT) - 1, 0))
        closed = max(range(nn), key=lambda i: x[i])  # x=5, z=0
        ts_rows.append([
            step, '%.10g' % t,
            '%.10g' % min(p), '%.10g' % max(p),
            '%.10g' % (max(p) / P0),
            '%.10g' % ac[closed], '%.10g' % ab[closed],
            '%.8f' % l2_rel(p, ac), '%.8f' % l2_rel(p, ab),
            '%.10g' % ab1[closed], '%.8f' % l2_rel(p, ab1),
        ])
    with open(ts_path, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['step', 'time_s', 'pmin_Pa', 'pmax_Pa', 'pmax_over_p0',
                    'analytic_continuum_pmax_Pa', 'analytic_BE_pmax_Pa',
                    'L2_rel_err_continuum', 'L2_rel_err_BE',
                    'analytic_BE_shift1_pmax_Pa', 'L2_rel_err_BE_shift1'])
        w.writerows(ts_rows)

    # ---- CSV 2: full pointwise field -------------------------------------
    fld_path = os.path.join(OUT, 'terzaghi_pressure_field.csv')
    with open(fld_path, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['step', 'time_s', 'node', 'x_m', 'y_m',
                    'p_fem_Pa', 'p_analytic_continuum_Pa',
                    'p_analytic_BE_Pa', 'p_analytic_BE_shift1_Pa',
                    'p_fem_minus_BE_Pa', 'p_fem_minus_BE_shift1_Pa'])
        for step, (t, p) in enumerate(zip(times, fields), start=1):
            ac = analytic_vec(t)
            ab = analytic_vec(t, be=True)
            ab1 = analytic_vec_steps(max(round(t / DT) - 1, 0))
            for i in range(nn):
                w.writerow([step, '%.10g' % t, i + 1,
                            '%.10g' % x[i], '%.10g' % y[i],
                            '%.10g' % p[i], '%.10g' % ac[i], '%.10g' % ab[i],
                            '%.10g' % ab1[i],
                            '%.10g' % (p[i] - ab[i]),
                            '%.10g' % (p[i] - ab1[i])])

    # ---- Console verification table --------------------------------------
    print('\n step  time     pmax_FEM   continuum  BE_series  L2err(cont) L2err(BE) L2err(BE-1)')
    for step in (1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51):
        r = ts_rows[step - 1]
        print(' %4s %6.3f  %9.1f %10.1f %10.1f   %7.3f%%   %7.3f%%    %7.3f%%' %
              (r[0], float(r[1]), float(r[3]), float(r[5]), float(r[6]),
               float(r[7]) * 100, float(r[8]) * 100, float(r[10]) * 100))

    # ---- PNG: four panels -------------------------------------------------
    draw_chart(times, fields, x, y, z, ts_rows)

    print('\nwrote:')
    for p in (ts_path, fld_path, os.path.join(OUT, 'terzaghi_validation.png')):
        print('  %s (%d bytes)' % (p, os.path.getsize(p)))


# ---------------------------------------------------------------------------
# Plotting (pure PIL)
# ---------------------------------------------------------------------------
W, H_IMG = 1640, 1040
BG = (250, 250, 252)
FG = (28, 32, 40)
GRID = (220, 224, 232)
C_FEM = (208, 48, 48)       # red
C_CONT = (60, 110, 190)     # blue
C_BE = (30, 150, 110)       # green
C_SHIFT = (145, 75, 170)    # purple: BE series shifted by one dt
PANEL = [(60, 55, 790, 490), (850, 55, 1580, 490),
         (60, 550, 790, 985), (850, 550, 1580, 985)]


def font(sz, bold=False):
    try:
        return ImageFont.truetype('arialbd.ttf' if bold else 'arial.ttf', sz)
    except OSError:
        return ImageFont.load_default()


def frame(dr, box, title, xlabel, ylabel, xmin, xmax, ymin, ymax, logy=False):
    x0, y0, x1, y1 = box
    L, R, T, B = x0 + 78, x1 - 18, y0 + 48, y1 - 52
    dr.rectangle([L, T, R, B], fill=(255, 255, 255), outline=GRID)
    dr.text((L, y0 + 14), title, fill=FG, font=font(17, True))
    dr.text(((L + R) // 2 - 70, y1 - 30), xlabel, fill=FG, font=font(13))
    dr.text((x0 + 6, (T + B) // 2 - 8), ylabel, fill=FG, font=font(13))

    def X(v):
        return L + (v - xmin) / (xmax - xmin) * (R - L)

    def Y(v):
        if logy:
            lo, hi = math.log10(ymin), math.log10(ymax)
            v = max(min(math.log10(max(v, ymin)), hi), lo)
            return B - (v - lo) / (hi - lo) * (B - T)
        return min(max(B - (v - ymin) / (ymax - ymin) * (B - T), T), B)

    tick = (ymax - ymin) / 4
    tickfmt = '%.3f' if tick < 0.02 else ('%.2f' if tick < 0.5 else
                                         ('%.1f' if tick < 2 else '%.0f'))
    for k in range(5):
        gy = T + k * (B - T) / 4
        dr.line([(L, gy), (R, gy)], fill=GRID)
        if logy:
            val = 10 ** (math.log10(ymax) - k * (math.log10(ymax) - math.log10(ymin)) / 4)
            lab = '%.0e' % val
        else:
            val = ymax - k * (ymax - ymin) / 4
            lab = tickfmt % val
        dr.text((L - 66, gy - 8), lab, fill=(110, 116, 128), font=font(11))
    for k in range(6):
        gx = L + k * (R - L) / 5
        dr.line([(gx, T), (gx, B)], fill=GRID)
        val = xmin + k * (xmax - xmin) / 5
        dr.text((gx - 12, B + 6), '%.2f' % val, fill=(110, 116, 128), font=font(11))
    return X, Y


def legend(dr, box, items):
    x0, y0, x1, y1 = box
    px, py = x1 - 220, y0 + 16
    for label, col, dash in items:
        if dash:
            for d in range(0, 26, 9):
                dr.line([(px + d, py), (px + min(d + 5, 26), py)],
                        fill=col, width=3)
        else:
            dr.line([(px, py), (px + 26, py)], fill=col, width=3)
        dr.text((px + 32, py - 8), label, fill=FG, font=font(12))
        py += 20


def poly(dr, pts, col, width=3, dash=False):
    if not dash:
        dr.line(pts, fill=col, width=width, joint='curve')
    else:
        for (ax, ay), (bx, by) in zip(pts[:-1], pts[1:]):
            n = max(1, int(math.hypot(bx - ax, by - ay) / 9))
            for k in range(0, n, 2):
                u0, u1 = k / n, min((k + 1) / n, 1.0)
                dr.line([(ax + (bx - ax) * u0, ay + (by - ay) * u0),
                         (ax + (bx - ax) * u1, ay + (by - ay) * u1)],
                        fill=col, width=width)


def draw_chart(times, fields, x, y, z, rows):
    im = Image.new('RGB', (W, H_IMG), BG)
    dr = ImageDraw.Draw(im)
    dr.text((60, 14), 'Gate 3.2  |  Terzaghi 1D consolidation: FEM vs analytic solution',
            fill=FG, font=font(22, True))
    dr.text((60, 40), 'CI #66  -  51 steps, dt=0.02 s, cv=108.696 m2/s, H=10 m, p0=1e5 Pa',
            fill=(100, 106, 118), font=font(13))

    # (a) pmax/p0 decay
    X, Y = frame(dr, PANEL[0], '(a) Peak pressure decay  p_max / p0',
                 'time t [s]', 'p / p0', 0, 1.0, 0, 1.05)
    fem = [(X(t), Y(max(f) / P0)) for t, f in zip(times, fields)]
    con = [(X(t), Y(float(r[5]) / P0)) for t, r in zip(times, rows)]
    be = [(X(t), Y(float(r[6]) / P0)) for t, r in zip(times, rows)]
    poly(dr, con, C_CONT, 3, dash=True)
    poly(dr, be, C_BE, 3)
    poly(dr, fem, C_FEM, 3)
    legend(dr, PANEL[0], [('FEM (CI #66)', C_FEM, False),
                          ('BE-corrected series', C_BE, False),
                          ('continuum series', C_CONT, True)])

    # (b) final profile along x, averaged over y
    step = len(fields) - 1
    t = times[step]
    p = fields[step]
    bins = {}
    for xi, pi in zip(x, p):
        key = round(xi, 6)
        bins.setdefault(key, []).append(pi)
    xs = sorted(bins)
    prof = [sum(bins[k]) / len(bins[k]) for k in xs]
    con_p = [analytic(5.0 - k, t) for k in xs]
    be_p = [analytic(5.0 - k, t, be=True) for k in xs]
    X, Y = frame(dr, PANEL[1], '(b) Pressure profile at t = 1.0 s',
                 'x [m]  (drained x=-5  ->  closed x=+5)', 'p / p0',
                 -5.0, 5.0, 0, 0.11)
    poly(dr, [(X(k), Y(v / P0)) for k, v in zip(xs, con_p)], C_CONT, 3, dash=True)
    poly(dr, [(X(k), Y(v / P0)) for k, v in zip(xs, be_p)], C_BE, 3)
    poly(dr, [(X(k), Y(v / P0)) for k, v in zip(xs, prof)], C_FEM, 3)
    for k, v in zip(xs, prof):
        dr.ellipse([X(k) - 2, Y(v / P0) - 2, X(k) + 2, Y(v / P0) + 2], fill=C_FEM)
    legend(dr, PANEL[1], [('FEM', C_FEM, False),
                          ('BE-corrected', C_BE, False),
                          ('continuum', C_CONT, True)])

    # (c) L2 error history
    X, Y = frame(dr, PANEL[2], '(c) Field L2 relative error vs time',
                 'time t [s]', 'L2 rel. error [%]', 0, 1.0, 0, 36.0)
    e_con = [(X(t), Y(float(r[7]) * 100)) for t, r in zip(times, rows)]
    e_be = [(X(t), Y(float(r[8]) * 100)) for t, r in zip(times, rows)]
    e_sh = [(X(t), Y(float(r[10]) * 100)) for t, r in zip(times, rows)]
    poly(dr, e_con, C_CONT, 3)
    poly(dr, e_be, C_BE, 3)
    poly(dr, e_sh, C_SHIFT, 3, dash=True)
    legend(dr, PANEL[2], [('vs continuum', C_CONT, False),
                          ('vs BE (n steps)', C_BE, False),
                          ('vs BE (n-1 steps)', C_SHIFT, True)])
    x0, y0, x1, y1 = PANEL[2]
    dr.text((x0 + 250, y0 + 78),
            't<0.06 spike: initial condition / Dirichlet ramp-in',
            fill=(100, 106, 118), font=font(12))
    dr.text((x0 + 250, y0 + 98),
            'early fields match BE(n-1) to 0.05-2%: one-step time-label offset',
            fill=(100, 106, 118), font=font(12))

    # (d) space-time heat map
    x0, y0, x1, y1 = PANEL[3]
    L, R, T, B = x0 + 78, x1 - 18, y0 + 48, y1 - 52
    dr.text((L, y0 + 14), '(d) Space-time pressure field  p(x,t)/p0',
            fill=FG, font=font(17, True))
    cw, ch = 640, 330
    small = Image.new('RGB', (cw, ch))
    px = small.load()
    nx = 160
    grid_x = [-5.0 + 10.0 * i / (nx - 1) for i in range(nx)]
    prof_t = []
    for tt, ff in zip(times, fields):
        b = {}
        for xi, pi in zip(x, ff):
            b.setdefault(round(xi, 6), []).append(pi)
        ks = sorted(b)
        vals = [sum(b[k]) / len(b[k]) for k in ks]
        prof_t.append([vals[min(range(len(ks)), key=lambda j: abs(ks[j] - gx))] / P0
                       for gx in grid_x])
    for iy in range(ch):
        it = min(len(times) - 1, int(iy / ch * len(times)))
        for ix in range(cw):
            ixv = min(nx - 1, int(ix / cw * nx))
            v = max(0.0, min(1.0, prof_t[it][ixv]))
            # white (0) -> deep blue (1)
            px[ix, iy] = (int(250 - 220 * v), int(250 - 175 * v), int(252 - 90 * v))
    im.paste(small, (L, T))
    dr.rectangle([L, T, L + cw, T + ch], outline=GRID)
    for k in range(6):
        gx = L + k * cw / 5
        dr.text((gx - 12, T + ch + 6), '%.0f' % (-5 + k * 2),
                fill=(110, 116, 128), font=font(11))
    for k in range(5):
        gy = T + k * ch / 4
        dr.text((L - 60, gy - 8), '%.2f' % (k * 0.25),
                fill=(110, 116, 128), font=font(11))
    dr.text((L + cw // 2 - 70, y1 - 30), 'x [m]', fill=FG, font=font(13))
    dr.text((x0 + 6, (T + B) // 2 - 8), 'time t [s]', fill=FG, font=font(13))
    # colour bar
    for k in range(cw):
        v = k / cw
        dr.line([(L + k, T + ch + 28), (L + k, T + ch + 38)],
                fill=(int(250 - 220 * v), int(250 - 175 * v), int(252 - 90 * v)))
    dr.text((L, T + ch + 40), '0', fill=(110, 116, 128), font=font(11))
    dr.text((L + cw - 24, T + ch + 40), '1.0', fill=(110, 116, 128), font=font(11))

    out = os.path.join(OUT, 'terzaghi_validation.png')
    im.save(out)


if __name__ == '__main__':
    main()
