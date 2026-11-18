#!/usr/bin/env python3
"""Check that anelastic scalar advection with map factors matches the compressible run.

Usage: check_mapfac_parity.py <anelastic plotfile> <compressible plotfile>

Both decks advect the same x-y cosine blob in a uniform wind with erf.test_mapfactor.
On every level the anelastic blob's x- and y-centroids (circular means over the
periodic domain, in cells) must lie within a quarter cell of the compressible
blob's, and its peak within 5 % of the compressible peak. Both blobs must also
have travelled the analytic distance from the domain centre, U t m / dx and
V t m / dy with the map factor m = 0.5 of erf.test_mapfactor, so
the check does not rest on the compressible run alone and cannot pass on a blob
that has not moved.
"""
import math
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, "..", "..", "Tools"))
import erf_plotfile  # noqa: E402

TOL_CELLS = 0.25
PEAK_TOL = 0.05
U, V, L, MF = 10.0, 5.0, 1000.0, 0.5     # decks: prob.U_0, prob.V_0, domain side, test map factor


def circular_centroid(s, k, axis):
    """Centroid in cells of level k along axis 0 (x) or 1 (y), minimum removed."""
    nx, ny = len(s), len(s[0])
    n = nx if axis == 0 else ny
    smin = min(s[i][j][k] for i in range(nx) for j in range(ny))
    cs = sn = 0.0
    for a in range(n):
        th = 2.0 * math.pi * (a + 0.5) / n
        if axis == 0:
            w = sum(s[a][j][k] - smin for j in range(ny))
        else:
            w = sum(s[i][a][k] - smin for i in range(nx))
        cs += w * math.cos(th)
        sn += w * math.sin(th)
    return (math.atan2(sn, cs) * n / (2.0 * math.pi)) % n, n


def wrapped(d, n):
    return abs((d + 0.5 * n) % n - 0.5 * n)


def main(argv):
    args = argv[1:]
    if len(args) != 2:
        print(__doc__)
        return 2
    ha, da = erf_plotfile.read_fields(args[0], ["scalar"])
    hc, dc = erf_plotfile.read_fields(args[1], ["scalar"])
    sa, sc = da["scalar"], dc["scalar"]
    nx, ny, nz = len(sc), len(sc[0]), len(sc[0][0])
    if (len(sa), len(sa[0]), len(sa[0][0])) != (nx, ny, nz):
        print("the two plotfiles have different grids")
        return 1
    ta, tc = ha["time"], hc["time"]
    if abs(ta - tc) > 1e-6 * max(1.0, tc):
        print("the two plotfiles are at different times: %g and %g s" % (ta, tc))
        return 1
    # blob centre starts at the domain centre; the analytic centroid in cells
    ex = (0.5 * nx + U * tc * MF / (L / nx)) % nx
    ey = (0.5 * ny + V * tc * MF / (L / ny)) % ny
    failed = 0
    print("t = %.1f s; analytic centroid x %.2f, y %.2f cells; tolerance %.2f cells, %.0f %% in the peak"
          % (tc, ex, ey, TOL_CELLS, 100 * PEAK_TOL))
    print("%6s %22s %22s %9s %9s %6s" % ("level", "x anel / comp [cells]", "y anel / comp [cells]",
                                       "peak an", "peak co", "pass"))
    for k in range(nz):
        xa, _ = circular_centroid(sa, k, 0)
        xc, _ = circular_centroid(sc, k, 0)
        ya, _ = circular_centroid(sa, k, 1)
        yc, _ = circular_centroid(sc, k, 1)
        pa = max(sa[i][j][k] for i in range(nx) for j in range(ny))
        pc = max(sc[i][j][k] for i in range(nx) for j in range(ny))
        parity = wrapped(xa - xc, nx) <= TOL_CELLS and wrapped(ya - yc, ny) <= TOL_CELLS and abs(pa - pc) <= PEAK_TOL * pc
        analytic = all(wrapped(c - e, n) <= TOL_CELLS for c, e, n in ((xa, ex, nx), (xc, ex, nx), (ya, ey, ny), (yc, ey, ny)))
        ok = parity and analytic
        failed += 0 if ok else 1
        print("%6d %10.2f / %9.2f %10.2f / %9.2f %9.3f %9.3f %6s%s"
              % (k, xa, xc, ya, yc, pa, pc, "yes" if ok else "NO",
                 "" if ok else "  (" + ", ".join(x for x, b in (("parity", parity), ("analytic", analytic)) if not b) + ")"))
    print("%d level(s) failed" % failed)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
