#!/usr/bin/env python3
"""Summarise DefMech AltMin displacement-SNES convergence from vDef logs.

Usage:
    python3 data/diagnose_snes.py label=/path/to/run.log [label=path ...]

Pure standard library.  CI #71 uses it to compare the production run with
two same-source variants (original L2 line search / SNESKSPONLY) without
recompiling anything.  The fnorm trace comes from ``-disp_snes_monitor``;
only the FIRST outer "Mechanics: step" block is isolated so the numbers are
not mixed across the 51 time steps.

Interpretation (PETSc 3.3, src/snes/impls/ls/ls.c):
  * reason -5 = SNES_DIVERGED_MAX_IT: 50 Newton iterations were spent and
    the convergence test never passed.
  * A fnorm trace that plateaus after one or two iterations points at an
    inconsistent F/J pair or an inaccurate inner KSP solve; a trace that
    oscillates points at the line search.  SNESLINESEARCHL2 fits a secant
    parabola in ||F||^2, while SNESLINESEARCHBASIC takes the full step.
"""

import re
import sys
from collections import Counter

RE_FAIL = re.compile(r"SNESSolve failed with SNESConvergedReason\s+(-?\d+)")
RE_REASON = re.compile(r"SNESConvergedReason\s+(-?\d+)")
RE_STEP = re.compile(r"Mechanics:\s*step\s+(\d+)")
RE_ALTMIN = re.compile(r"Alt\.\s*Min\.\s*step\s+(\d+)")
RE_NORM = re.compile(r"^\s*(\d+)\s+SNES Function norm\s+(\S+)")


def analyse(path):
    with open(path, "r", errors="replace") as handle:
        lines = handle.read().splitlines()
    text = "\n".join(lines)
    fails = Counter(int(x) for x in RE_FAIL.findall(text))
    reasons = Counter(int(x) for x in RE_REASON.findall(text))
    outer = [int(x) for x in RE_STEP.findall(text)]
    altmin = [int(x) for x in RE_ALTMIN.findall(text)]

    starts = [i for i, line in enumerate(lines) if RE_STEP.search(line)]
    if starts:
        stop = starts[1] if len(starts) > 1 else len(lines)
        block = lines[starts[0]:stop]
    else:
        block = lines
    norms = []
    for line in block:
        match = RE_NORM.match(line)
        if match:
            try:
                norms.append((int(match.group(1)), float(match.group(2))))
            except ValueError:
                pass
    return {
        "fails": fails,
        "reasons": reasons,
        "outer": len(outer),
        "altmin": len(altmin),
        "norms": norms,
    }


def fmt(value, width=10):
    return ("%.*e" % (3, value)).rjust(width)


def main(argv):
    if not argv:
        print(__doc__)
        return 2
    rows = []
    for spec in argv:
        if "=" not in spec:
            print("=== SKIP: expected label=path, got %r ===" % spec)
            continue
        label, path = spec.split("=", 1)
        try:
            rows.append((label, path, analyse(path)))
        except IOError as exc:
            print("=== SKIP: %s (%s) ===" % (label, exc))

    print("=== DefMech AltMin SNES diagnosis (CI #71) ===")
    header = "%-14s %-8s %-6s %-6s %-5s %-11s %-11s %-11s %s" % (
        "label", "fails", "outer", "altMin", "iters", "fnorm[0]", "fnorm[min]", "fnorm[-1]", "verdict")
    print(header)
    print("-" * len(header))
    for label, path, info in rows:
        norms = info["norms"]
        values = [v for _, v in norms]
        total_fail = sum(info["fails"].values())
        verdict = "CONVERGED" if total_fail == 0 else "DIVERGED"
        print("%-14s %-8s %-6d %-6d %-5d %s %s %s %s" % (
            label,
            str(total_fail) + ("(" + ",".join("%d:%d" % kv for kv in sorted(info["fails"].items())) + ")" if total_fail else ""),
            info["outer"],
            info["altmin"],
            len(norms),
            fmt(values[0]) if values else fmt(0.0),
            fmt(min(values)) if values else fmt(0.0),
            fmt(values[-1]) if values else fmt(0.0),
            verdict,
        ))
    print("")
    for label, path, info in rows:
        norms = info["norms"]
        print("--- %s: %s" % (label, path))
        if info["reasons"]:
            print("    SNESConvergedReason histogram: %s" %
                  ", ".join("%d x %d" % (count, reason) for reason, count in sorted(info["reasons"].items())))
        if not norms:
            print("    (no -disp_snes_monitor output)")
            continue
        head = norms[:8]
        tail = norms[-4:] if len(norms) > 12 else []
        print("    first norms: %s" % ", ".join("%d:%s" % (i, fmt(v).strip()) for i, v in head))
        if tail:
            print("    last  norms: %s" % ", ".join("%d:%s" % (i, fmt(v).strip()) for i, v in tail))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
