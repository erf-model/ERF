#!/usr/bin/env python3
"""Numerical checks for the 3D axisymmetric Witch of Agnesi hill k-eqn RANS case.

Usage: check_hill3d.py [--smoke | --physics] <plotfile>

--smoke (default, the CTest entry) checks the wall distance against the
exact distance to the surface of revolution on a sub-sample of columns,
the length-scale bounds built on it, positivity and the wall-cell k
retention over the whole surface; --physics adds the crest speed-up and
the upstream log law after the 4 h run.
"""

import math
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, "..")]
import erf_plotfile  # noqa: E402
import rans_checks as rc  # noqa: E402

Z0 = 0.1
L_G_MAX = 30.0
HMAX = 100.0
L = 500.0
FIELDS = ["x_velocity", "y_velocity", "KE", "Kmv", "Lturb", "walldist", "diss", "density", "z_phys"]


def hill(r):
    return HMAX / (1.0 + (r / L) ** 2)


def exact_distance(r, z):
    """Shortest distance from a point at radius r, height z to the surface
    of revolution z = hill(r): the foot point lies in the meridional plane,
    so this is a 1D search over the foot radius."""
    def dist(s):
        return math.hypot(r - s, z - hill(s))
    best_s, best_d = r, dist(r)
    span = 4.0 * L
    n = 80
    for m in range(-n, n + 1):
        s = max(r + span * m / n, 0.0)
        d = dist(s)
        if d < best_d:
            best_s, best_d = s, d
    a, b = max(best_s - span / n, 0.0), best_s + span / n
    gr = 0.5 * (math.sqrt(5.0) - 1.0)
    c, d_ = b - gr * (b - a), a + gr * (b - a)
    for _ in range(50):
        if dist(c) < dist(d_):
            b, d_ = d_, c
            c = b - gr * (b - a)
        else:
            a, c = c, d_
            d_ = a + gr * (b - a)
    return dist(0.5 * (a + b))


def main(argv):
    mode, plt, _ = rc.parse_args(argv, __doc__)
    hdr, f = erf_plotfile.read_fields(plt, FIELDS)
    nx, ny, nz = [hdr["hi"][d] - hdr["lo"][d] + 1 for d in range(3)]
    dx, dy = hdr["dx"][0], hdr["dx"][1]
    x_lo, y_lo = hdr["prob_lo"][0], hdr["prob_lo"][1]
    xcen = 0.5 * (hdr["prob_lo"][0] + hdr["prob_hi"][0])
    ycen = 0.5 * (hdr["prob_lo"][1] + hdr["prob_hi"][1])
    rep = rc.Report()

    # --- wall distance on every second column ---------------------------
    max_rel, sum_rel, n_rel, max_abs_near = 0.0, 0.0, 0, 0.0
    for i in range(0, nx, 2):
        x = x_lo + (i + 0.5) * dx - xcen
        for j in range(0, ny, 2):
            y = y_lo + (j + 0.5) * dy - ycen
            r = math.hypot(x, y)
            for k in range(nz):
                z = f["z_phys"][i][j][k]
                d_exact = exact_distance(r, z)
                err = f["walldist"][i][j][k] - d_exact
                if d_exact < 100.0:
                    max_abs_near = max(max_abs_near, abs(err))
                rel = abs(err) / d_exact
                max_rel = max(max_rel, rel)
                sum_rel += rel
                n_rel += 1
    # Tucker's Poisson distance is weakest at the convex crest (about 10 %
    # in the first cell there); the absolute error near the surface is
    # judged in units of the vertical cell size.
    rep.check("max rel err walldist vs exact hill distance", max_rel, 0.0, 0.15)
    rep.check("mean rel err walldist vs exact hill distance", sum_rel / n_rel, 0.0, 0.03)
    rep.check("max abs err walldist, d < 100 m [cells]", max_abs_near / hdr["dx"][2], 0.0, 0.2)

    # --- structural checks on the whole field -----------------------------
    finite = all(math.isfinite(f[fl][i][j][k]) for fl in FIELDS for i in range(nx) for j in range(ny) for k in range(nz))
    rep.check("all fields finite", 1.0 if finite else 0.0, 1.0, 0.0)
    rep.check("min KE [m2/s2]", min(f["KE"][i][j][k] for i in range(nx) for j in range(ny) for k in range(nz)), 0.0, 0.0, "min")
    rep.check("min Kmv [kg/m/s]", min(f["Kmv"][i][j][k] for i in range(nx) for j in range(ny) for k in range(nz)), 0.0, 0.0, "min")
    worst = 0.0
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                bound = rc.unstable_bound(rc.geom_length(f["walldist"][i][j][k], Z0, L_G_MAX))
                worst = max(worst, f["Lturb"][i][j][k] / bound)
    rep.check("max Lturb / bound (on walldist)", worst, 1.0, 1e-6, "max")
    worst = 0.0
    for i in range(nx):
        for j in range(ny):
            rho, ke, lt, ds = f["density"][i][j][0], f["KE"][i][j][0], f["Lturb"][i][j][0], f["diss"][i][j][0]
            worst = max(worst, abs((ds / (rho * rc.CMU0 ** 3 * ke ** 1.5 / lt)) ** (2.0 / 3.0) - 1.0))
    rep.check("max |wall cell k_start/k_end - 1|", worst, 0.0, 0.01)

    # --- physics: speed-up at the crest and the upstream log law ---------
    if mode == "physics":
        i_top = min(range(nx), key=lambda i: abs(x_lo + (i + 0.5) * dx - xcen))
        j_top = min(range(ny), key=lambda j: abs(y_lo + (j + 0.5) * dy - ycen))
        i_up = min(range(nx), key=lambda i: abs(x_lo + (i + 0.5) * dx - (xcen - 2000.0)))
        speed = lambda i, j, k: math.hypot(f["x_velocity"][i][j][k], f["y_velocity"][i][j][k])
        # axisymmetric hill: fractional speed-up near the crest about 1.6 h/L
        # (Jackson & Hunt 1975 for the 2D ridge gives 2 h/L; flow around the
        # sides lowers it for a hill of revolution)
        est = 1.6 * HMAX / L
        for k in (0, 1, 2):
            rep.check("crest speed-up (U_top/U_up - 1), k=%d" % k, speed(i_top, j_top, k) / speed(i_up, j_top, k) - 1.0,
                      (0.5 * est, 2.0 * est), 0.0, "range")
        rep.check("speed-up positive in the lowest 8 cells", min(speed(i_top, j_top, k) - speed(i_up, j_top, k) for k in range(8)), 0.0, 0.0, "min")
        ustar = math.sqrt(f["KE"][i_up][j_top][0]) * rc.CMU0
        rep.check("upstream u* from wall k [m/s]", ustar, (0.25, 0.55), 0.0, "range")
        for k in range(3):
            d = f["walldist"][i_up][j_top][k]
            rep.check("upstream |U|(k=%d) vs log law" % k, speed(i_up, j_top, k), ustar / rc.KAPPA * math.log((d + Z0) / Z0), 0.15, "rel")

    rep.dump()
    return 1 if rep.failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
