#!/usr/bin/env python3
"""Same-source comparison of two KSP choices for the transient Darcy TS.

Run vDef TWICE with the identical binary, mesh and YAML deck, changing only
-flow_ksp_type (gmres vs cg, both with -flow_ksp_rtol 1.0e-10), then compare
the two Exodus result meshes and both solver logs.

CI #69 used this to decide the final KSP choice: CG reproduces GMRES to
2.54e-5 Pa (2.5e-10 relative) with 35% fewer linear iterations, so CI #70
dropped the GMRES command-line override and restored the CreateTS default
KSPCG. The one-off double run was then removed from the workflow; this
script is kept so the comparison can be repeated at any time. The #69 inputs
are archived in the 'terzaghi-ksp-comparison' artifact.

Pure standard library (plus the committed parse_cdf2.py); no numpy / scipy /
matplotlib on the CI runner.

Usage:
    python3 data/compare_ksp_runs.py <gen_gmres> <gen_cg> <log_gmres> <log_cg>

Exit codes:
    0  comparison complete. Either the two fields agree, or CG was found NOT
       VIABLE (it diverged / did not finish) - a valid verdict that says to
       keep the GMRES command-line override.
    1  both runs converged but the pressure fields disagree, or an input is
       broken. This is a real numerical regression and must fail the CI job.
"""

import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse_cdf2 import CDF2  # noqa: E402

P0 = 1.0e5                 # initial pressure [Pa]
TOL_ABS = 1.0e-2           # max allowed |p_gmres - p_cg| [Pa] (1e-7 relative)
TOL_RATIO = 1.0e-5         # max allowed difference of the decay ratios
RATIO_BAND = (0.05, 0.15)  # Gate 3.2 analytic band

ENVELOPE = re.compile(
    r"Gate 3\.2:\s+step\s+(\d+)\s+pmin\s+([-+0-9.eEdD]+)\s+pmax\s+([-+0-9.eEdD]+)")
CONVERGED = re.compile(
    r"Linear solve converged due to ([A-Z_]+) iterations (\d+)")
DIVERGED = re.compile(
    r"Linear solve did not converge|KSP_DIVERGED_|DIVERGED_[A-Z_]+")
# -flow_ts_view prints a "KSP Object: (flow_)" block on every step. Checking
# its type proves the -flow_ksp_type option really took effect (the DefMech
# KSP objects are also CG, so the (flow_) tag must be matched explicitly).
KSPTYPE = re.compile(
    r"KSP Object:\s+\(flow_\)[^\n]*\n\s*type:\s*([A-Za-z0-9_]+)")
TIMESTAMP = re.compile(r"(?m)^\d{4}-\d\d-\d\dT[0-9:.]+Z\s*")


def as_float(s):
    return float(s.replace('D', 'E').replace('d', 'e'))


def read_log(path):
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        # GitHub Actions job logs prefix every line with a UTC timestamp.
        # Strip it so archived job logs can be replayed locally as well.
        text = TIMESTAMP.sub('', f.read())

    # De-duplicate by step number. An archived GitHub job log also contains
    # the shell's own `tail` / `grep` echoes of these lines, out of order;
    # taking the textual last line would then report step 50 instead of 51.
    env = {}
    for m in ENVELOPE.finditer(text):
        env[int(m.group(1))] = (as_float(m.group(2)), as_float(m.group(3)))
    steps = sorted(env)
    conv = [(m.group(1), int(m.group(2))) for m in CONVERGED.finditer(text)]
    iters = [n for _, n in conv]
    reasons = {}
    for why, _ in conv:
        reasons[why] = reasons.get(why, 0) + 1

    ksp_types = {}
    for name in KSPTYPE.findall(text):
        ksp_types[name] = ksp_types.get(name, 0) + 1
    ksp_type = max(ksp_types, key=ksp_types.get) if ksp_types else 'unknown'

    first = env[steps[0]][1] if steps else 0.0
    last = env[steps[-1]][1] if steps else 0.0
    return dict(
        path=path,
        ksp_type=ksp_type,
        ksp_types=ksp_types,
        n_probe=len(steps),
        step_first=steps[0] if steps else 0,
        step_last=steps[-1] if steps else 0,
        pmax_first=first,
        pmax_last=last,
        ratio=(last / first) if first else 0.0,
        solves=len(conv),
        iter_total=sum(iters),
        iter_mean=(sum(iters) / len(iters)) if iters else 0.0,
        iter_max=max(iters) if iters else 0,
        reasons=reasons,
        diverged=len(DIVERGED.findall(text)),
    )


def read_pressure(path):
    c = CDF2(path)
    key = 'vals_nod_var1'
    fields = [c.read_record(key, i) for i in range(c.numrecs)]
    return c.numrecs, len(fields[0]) if fields else 0, fields


def viable(log):
    return (log['diverged'] == 0 and log['solves'] >= 48 and
            log['n_probe'] >= 50 and
            RATIO_BAND[0] <= log['ratio'] <= RATIO_BAND[1])


def log_report(tag, log):
    why = ', '.join('%s x%d' % (k, v) for k, v in sorted(log['reasons'].items())) or 'none'
    print("  %-6s : flow KSP type '%s' %s" %
          (tag, log['ksp_type'], log['ksp_types']))
    print("           %2d envelopes (steps %d..%d), pmax %.6g -> %.6g, ratio %.7f" %
          (log['n_probe'], log['step_first'], log['step_last'],
           log['pmax_first'], log['pmax_last'], log['ratio']))
    print("           %d converged solves (%s), iterations total=%d mean=%.2f max=%d, diverged=%d" %
          (log['solves'], why, log['iter_total'], log['iter_mean'],
           log['iter_max'], log['diverged']))


def main(gen_gmres, gen_cg, log_gmres, log_cg):
    print("=== Same-source KSP comparison: GMRES vs CG, same binary / mesh / deck ===")
    print("  rtol = 1.0e-10 for both; only -flow_ksp_type differs")

    gm = read_log(log_gmres)
    cg = read_log(log_cg)
    print("\n--- solver logs ---")
    log_report('GMRES', gm)
    log_report('CG', cg)

    if gm['ksp_type'] != 'gmres' or cg['ksp_type'] != 'cg':
        print("\n=== FAIL: -flow_ksp_type did not take effect "
              "(GMRES log='%s', CG log='%s'); the comparison is meaningless ===" %
              (gm['ksp_type'], cg['ksp_type']))
        return 1
    if not viable(gm):
        print("\n=== FAIL: the GMRES reference run itself is not healthy ===")
        return 1

    n_g, nn_g, pg = read_pressure(gen_gmres)
    n_c, nn_c, pc = read_pressure(gen_cg)
    print("\n--- pressure field (vals_nod_var1) ---")
    print("  GMRES : %d records x %d nodes" % (n_g, nn_g))
    print("  CG    : %d records x %d nodes" % (n_c, nn_c))

    if (n_g, nn_g) != (n_c, nn_c):
        print("\n=== FAIL: the two result meshes have different shapes ===")
        return 1

    max_abs = 0.0
    max_pmax_delta = 0.0
    worst = (0, 0)
    for i in range(n_g):
        a, b = pg[i], pc[i]
        for j in range(nn_g):
            d = abs(a[j] - b[j])
            if d > max_abs:
                max_abs, worst = d, (i, j)
        max_pmax_delta = max(max_pmax_delta, abs(max(a) - max(b)))

    print("  max |p_GMRES - p_CG|        = %.6g Pa  (record %d, node %d)" %
          (max_abs, worst[0], worst[1]))
    print("  max |d| / p0                = %.6g" % (max_abs / P0))
    print("  max |pmax_GMRES - pmax_CG|  = %.6g Pa" % max_pmax_delta)
    print("  t=1 pmax : GMRES %.6f   CG %.6f" % (max(pg[-1]), max(pc[-1])))

    print("\n--- decay ratio (Gate 3.2 analytic band [0.05, 0.15]) ---")
    print("  GMRES %.7f   CG %.7f   |delta| = %.6g" %
          (gm['ratio'], cg['ratio'], abs(gm['ratio'] - cg['ratio'])))

    print("\n--- linear solver cost ---")
    if cg['iter_total'] and gm['iter_total']:
        delta = 100.0 * (cg['iter_total'] - gm['iter_total']) / gm['iter_total']
        print("  CG total iterations vs GMRES: %+.2f%% (%d vs %d)" %
              (delta, cg['iter_total'], gm['iter_total']))

    print("\n--- verdict ---")
    if not viable(cg):
        print("  CG NOT VIABLE (diverged=%d, converged solves=%d, envelopes=%d, ratio=%.7f)" %
              (cg['diverged'], cg['solves'], cg['n_probe'], cg['ratio']))
        print("  KEEP the -flow_ksp_type gmres command-line override; the "
              "CreateTS default KSPCG must not be used unguarded.")
        print("  (informational field delta above is not used as a pass/fail "
              "criterion when CG did not converge)")
        return 0

    if max_abs > TOL_ABS:
        print("  FAIL: both runs converged but the fields differ by %.6g Pa > %.6g Pa" %
              (max_abs, TOL_ABS))
        return 1
    if abs(gm['ratio'] - cg['ratio']) > TOL_RATIO:
        print("  FAIL: decay ratios differ by more than %.6g" % TOL_RATIO)
        return 1

    print("  PASS: CG reproduces the GMRES pressure field to %.6g Pa "
          "(<= %.6g Pa); the KSP choice does not affect the physics." %
          (max_abs, TOL_ABS))
    print("  -> KSPCG is a valid default for this problem; no GMRES override "
          "is required (CI #70 dropped it).")
    return 0


if __name__ == '__main__':
    if len(sys.argv) != 5:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(*sys.argv[1:5]))
