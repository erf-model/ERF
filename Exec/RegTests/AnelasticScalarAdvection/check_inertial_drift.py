#!/usr/bin/env python3
"""Check a passive-scalar blob against its exact position in an inertial oscillation.

Usage: check_inertial_drift.py <plotfile>

The blob starts at x = y = 500 m and the wind is u = U cos(f t), v = -U sin(f t), so it
reaches x = 500 + U sin(f t) / f, y = 500 + U (cos(f t) - 1) / f.  The centroid, a circular
mean over the periodic domain, must be within 0.5 m of that (a cell is 15.6 m); the wind
must also have turned by the exact angle within 1 % of U, and the blob must have moved.
"""
import math
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, "..", "..", "Tools"))
import erf_plotfile  # noqa: E402

U, PERIOD, L, X0, K = 10.0, 100.0, 1000.0, 500.0, 1
F = 4.0 * math.pi / PERIOD
TOL = 0.5


def centroid(s, axis):
    nx, ny = len(s), len(s[0])
    n = nx if axis == 0 else ny
    smin = min(s[i][j][K] for i in range(nx) for j in range(ny))
    cs = sn = 0.0
    for a in range(n):
        th = 2.0 * math.pi * (a + 0.5) / n
        if axis == 0:
            w = sum(s[a][j][K] - smin for j in range(ny))
        else:
            w = sum(s[i][a][K] - smin for i in range(nx))
        cs += w * math.cos(th)
        sn += w * math.sin(th)
    return (math.atan2(sn, cs) * L / (2.0 * math.pi)) % L


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    hdr, d = erf_plotfile.read_fields(argv[1], ["scalar", "x_velocity", "y_velocity"])
    t = hdr["time"]
    xe = X0 + U * math.sin(F * t) / F
    ye = X0 + U * (math.cos(F * t) - 1.0) / F
    x, y = centroid(d["scalar"], 0), centroid(d["scalar"], 1)
    err = math.hypot(x - xe, y - ye)
    u, v = d["x_velocity"][3][3][K], d["y_velocity"][3][3][K]
    werr = math.hypot(u - U * math.cos(F * t), v + U * math.sin(F * t))
    print("t = %.2f s: centroid (%.3f, %.3f) m, exact (%.3f, %.3f) m, error %.4f m (limit %.1f)"
          % (t, x, y, xe, ye, err, TOL))
    print("wind (%.4f, %.4f) m/s, error %.2e m/s (limit %.2f)" % (u, v, werr, 0.01 * U))
    failed = 0
    if not err <= TOL:
        print("the blob is not at its exact position")
        failed += 1
    if not werr <= 0.01 * U:
        print("the wind is not the inertial oscillation")
        failed += 1
    if t <= 0.0 or math.hypot(xe - X0, ye - X0) < 10.0:
        print("the blob has not moved far enough to test anything")
        failed += 1
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
