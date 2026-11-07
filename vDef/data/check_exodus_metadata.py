#!/usr/bin/env python3
"""Gate 3.2 metadata check for the Terzaghi Exodus output (CI #68).

Pure standard library (uses the committed parse_cdf2.py). Verifies the whole
nodal variable table and reports each slot's value range at the first and last
record, so a silent overwrite of a DefMech field by a Flow field is visible in
the CI log instead of being discovered during post-processing.

Layout for vDef/data/terzaghi1D.yaml (2D, HeatXfer disabled):
  slot 1  pressureOffset            -> Pressure        (Flow, nodal Vec)
  slot 2  damageOffset              -> Damage          (DefMech)
  slot 3  displacementOffset        -> Displacement_X  (DefMech, vector component 1)
  slot 4  displacementOffset + 1    -> Displacement_Y  (DefMech, vector component 2)

Before CI #68 the YAML also set boundary_pressure_offset: 3. ViewEXO is called
DefMech first and Flow second, so the constant boundaryPressure Vec overwrote
Displacement_X on every step. That offset is now 0 (no output slot).
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse_cdf2 import CDF2  # noqa: E402

EXPECTED = ["Pressure", "Damage", "Displacement_X", "Displacement_Y"]
NAME_WIDTH = 33  # MXSTLN in the Exodus/NetCDF CDF2 header


def main(path):
    c = CDF2(path)
    raw = b"".join(c.read_static("name_nod_var"))
    names = [raw[i:i + NAME_WIDTH].rstrip(b"\x00").decode("ascii", "replace")
             for i in range(0, len(raw), NAME_WIDTH)]

    print("numrecs = %d" % c.numrecs)
    print("name_nod_var = %s" % " | ".join(names))
    if names != EXPECTED:
        print("FAIL: expected %s" % " | ".join(EXPECTED))
        return 1

    last = c.numrecs - 1
    for k, name in enumerate(names, 1):
        key = "vals_nod_var%d" % k
        try:
            a = c.read_record(key, 0)
            b = c.read_record(key, last)
            print("slot %d %-18s t0 [%+.6g, %+.6g]  t1 [%+.6g, %+.6g]" %
                  (k, name, min(a), max(a), min(b), max(b)))
            if name == "Displacement_X" and min(b) == 0.0 and max(b) == 0.0:
                print("[WARN] Displacement_X is identically zero at t=1; "
                      "check that Flow is not overwriting nodal slot 3 "
                      "(boundary_pressure_offset must be 0).")
        except KeyError:
            print("FAIL: %s missing from the Exodus file" % key)
            return 1

    print("PASS: nodal variable table matches the vDef/DefMech/Flow layout")
    return 0


if __name__ == "__main__":
    default = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "Beam1x10-tri3_out.gen")
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else default))
