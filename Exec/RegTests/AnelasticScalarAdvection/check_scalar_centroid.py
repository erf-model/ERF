#!/usr/bin/env python3
"""Check that an anelastic passive-scalar blob moves at the wind speed on every level.

Usage: check_scalar_centroid.py <plotfile>

The deck advects an x-y cosine blob (uniform in z) in a uniform 10 m/s wind on a
stretched vertical mesh. On every level the blob's x-centroid, taken as a circular
mean over the periodic domain, must sit at x_c + U t within one cell, and its peak
must keep at least 80 % of the initial peak (prob_type 11 gives A_0/2 at the centre).
"""
import math
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, "..", "..", "Tools"))
import erf_plotfile  # noqa: E402

L, U, XC, A0 = 1000.0, 10.0, 500.0, 1.0
PEAK0 = 0.5 * A0        # A_0 (1 + cos(0)) / 4


def main(argv):
    args = argv[1:]
    if len(args) != 1:
        print(__doc__)
        return 2
    hdr, d = erf_plotfile.read_fields(args[0], ["scalar"])
    t = hdr["time"]
    s = d["scalar"]
    nx, ny, nz = len(s), len(s[0]), len(s[0][0])
    dx = L / nx
    expected = (XC + U * t) % L
    failed = 0
    print("t = %.1f s, expected x-centroid %.1f m, tolerance %.2f m (one cell)" % (t, expected, dx))
    print("%6s %14s %12s %10s %6s" % ("level", "centroid [m]", "error [m]", "peak", "pass"))
    for k in range(nz):
        smin = min(s[i][j][k] for i in range(nx) for j in range(ny))
        cs = sn = 0.0
        peak = 0.0
        for i in range(nx):
            th = 2.0 * math.pi * (i + 0.5) * dx / L
            w = sum(s[i][j][k] - smin for j in range(ny))
            cs += w * math.cos(th)
            sn += w * math.sin(th)
            peak = max(peak, max(s[i][j][k] for j in range(ny)))
        xbar = (math.atan2(sn, cs) * L / (2.0 * math.pi)) % L
        err = abs((xbar - expected + 0.5 * L) % L - 0.5 * L)
        ok = err <= dx and peak >= 0.8 * PEAK0
        failed += 0 if ok else 1
        print("%6d %14.1f %12.1f %10.3f %6s" % (k, xbar, err, peak, "yes" if ok else "NO"))
    print("%d level(s) failed" % failed)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
