#!/usr/bin/env python3
"""Check that anelastic advection over terrain leaves a uniform scalar and theta uniform.

Usage: check_uniform_scalar.py <plotfile> <tolerance> [--require-z-split] [--require-level]

On every level of the plotfile, scalar must stay within the given relative tolerance of
its value at the start (prob.A_0 = 1), and theta of prob.T_0 = 300.  So that the check
cannot pass on a flow that never crosses the terrain-following faces, the largest |w|
must also be at least 0.5 m/s (the wind over the hill's flanks).  With --require-z-split
level 0 must also hold a box that starts above the bottom of the domain, and with
--require-level the plotfile must hold a refined level.  A NaN anywhere fails.
"""
import math
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, "..", "..", "Tools"))
import erf_plotfile  # noqa: E402

EXPECTED = {"scalar": 1.0, "theta": 300.0}
MIN_W = 0.5


def extrema(plotfile, level, names, name):
    """Smallest and largest value of a field over the valid cells of a level."""
    level_dir = os.path.join(plotfile, "Level_%d" % level)
    comp = names.index(name)
    ncomp, _, _, fabs = erf_plotfile._read_cell_h(level_dir)
    lo, hi = float("inf"), float("-inf")
    for fname, offset in fabs:
        flo, fhi, data = erf_plotfile._read_fab(os.path.join(level_dir, fname), offset, ncomp)
        npts = 1
        for a, b in zip(flo, fhi):
            npts *= b - a + 1
        values = data[comp * npts:(comp + 1) * npts]
        if any(math.isnan(v) for v in values):
            return float("nan"), float("nan")
        lo, hi = min(lo, min(values)), max(hi, max(values))
    return lo, hi


def main(argv):
    args = [a for a in argv[1:] if a not in ("--require-z-split", "--require-level")]
    if len(args) != 2:
        print(__doc__)
        return 2
    plotfile, tol = args[0], float(args[1])
    hdr = erf_plotfile._read_header(plotfile)
    failed = 0
    largest_w = 0.0
    print("%6s %8s %14s %14s %12s %6s" % ("level", "field", "min", "max", "rel. error", "pass"))
    for lev in range(hdr["finest"] + 1):
        for name, value in EXPECTED.items():
            lo, hi = extrema(plotfile, lev, hdr["names"], name)
            err = max(abs(lo - value), abs(hi - value)) / value
            ok = err <= tol
            failed += 0 if ok else 1
            print("%6d %8s %14.8g %14.8g %12.3e %6s" % (lev, name, lo, hi, err, "yes" if ok else "NO"))
        lo, hi = extrema(plotfile, lev, hdr["names"], "z_velocity")
        largest_w = max(largest_w, abs(lo), abs(hi))
    print("largest |w| = %.3f m/s, tolerance %.1e" % (largest_w, tol))
    if largest_w < MIN_W:
        print("the wind does not cross the hill: |w| below %.1f m/s" % MIN_W)
        failed += 1
    if "--require-level" in argv and hdr["finest"] < 1:
        print("the plotfile has no refined level")
        failed += 1
    if "--require-z-split" in argv:
        _, _, boxes, _ = erf_plotfile._read_cell_h(os.path.join(plotfile, "Level_0"))
        if not any(lo[2] > hdr["lo"][2] for lo, _ in boxes):
            print("the boxes of level 0 are not split in z")
            failed += 1
    print("%d check(s) failed" % failed)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
