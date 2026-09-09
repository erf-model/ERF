#!/usr/bin/env python3
"""Numerical checks for the 2D Witch of Agnesi ridge k-eqn RANS case.

Usage: check_hill.py [--smoke | --physics] [--flat] <plotfile>

--smoke (default, the CTest entry) checks the wall distance against the
exact distance to the ridge, the length-scale bounds built on it, positivity
and the wall-cell k retention; --physics adds the hill-top speed-up and
the far-field log law after the 6 h run. --flat is for the flat-fitted
variant (prob.hmax = 1e-6), where the Poisson wall distance must equal
the height above the surface to solver tolerance.
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
UG = 10.0
FIELDS = ["x_velocity", "y_velocity", "KE", "Kmv", "Lturb", "walldist", "diss", "density", "z_phys"]


def ridge(x, hmax, L):
    return hmax / (1.0 + (x / L) ** 2)


def exact_distance(x, z, hmax, L):
    """Shortest distance from (x, z) to the curve z = ridge(x): a scan over
    the foot point with a golden-section refinement around the best sample."""
    best_s, best_d = x, abs(z - ridge(x, hmax, L))
    span = 4.0 * L
    n = 200
    for m in range(-n, n + 1):
        s = x + span * m / n
        d = math.hypot(x - s, z - ridge(s, hmax, L))
        if d < best_d:
            best_s, best_d = s, d
    a, b = best_s - span / n, best_s + span / n
    gr = 0.5 * (math.sqrt(5.0) - 1.0)
    c, d_ = b - gr * (b - a), a + gr * (b - a)
    for _ in range(60):
        fc = math.hypot(x - c, z - ridge(c, hmax, L))
        fd = math.hypot(x - d_, z - ridge(d_, hmax, L))
        if fc < fd:
            b, d_ = d_, c
            c = b - gr * (b - a)
        else:
            a, c = c, d_
            d_ = a + gr * (b - a)
    s = 0.5 * (a + b)
    return math.hypot(x - s, z - ridge(s, hmax, L))


def main(argv):
    flat = "--flat" in argv
    argv = [a for a in argv if a != "--flat"]
    mode, plt, _ = rc.parse_args(argv, __doc__)
    hdr, f = erf_plotfile.read_fields(plt, FIELDS)
    nx, ny, nz = [hdr["hi"][d] - hdr["lo"][d] + 1 for d in range(3)]
    dx = hdr["dx"][0]
    x_lo = hdr["prob_lo"][0]
    xcen = 0.5 * (hdr["prob_lo"][0] + hdr["prob_hi"][0])
    hmax = 1e-6 if flat else 100.0
    L = 500.0
    j = 0
    rep = rc.Report()

    def col(field, i):
        return [f[field][i][j][k] for k in range(nz)]

    # --- wall distance against the exact distance to the ridge -----------
    # The geometric length is capped at max_geom_lscale (30 m here), so the
    # wall distance only shapes the closure within roughly 100 m of the
    # surface; the absolute error is judged there, the relative error
    # everywhere.
    max_rel, sum_rel, n_rel, max_abs_near = 0.0, 0.0, 0, 0.0
    for i in range(nx):
        x = x_lo + (i + 0.5) * dx
        for k in range(nz):
            z = f["z_phys"][i][j][k]
            d_exact = exact_distance(x - xcen, z, hmax, L)
            err = f["walldist"][i][j][k] - d_exact
            if d_exact < 100.0:
                max_abs_near = max(max_abs_near, abs(err))
            rel = abs(err) / d_exact
            max_rel = max(max_rel, rel)
            sum_rel += rel
            n_rel += 1
    if flat:
        # exact to 1e-6 m above the first cell; the first cell carries the
        # odd-reflection Dirichlet ghost of the Poisson solve (dz^2/8 in phi,
        # 1.5 cm in distance here), 0.2 % of its 7.8 m
        rep.check("max rel err walldist vs z - h (flat fitted)", max_rel, 0.0, 3e-3)
        rep.check("max abs err walldist, d < 100 m (flat) [m]", max_abs_near, 0.0, 0.1)
    else:
        # Tucker's Poisson distance is weakest at the convex crest: about
        # 5-10 % in the first cells there. The absolute error near the
        # surface is judged in units of the vertical cell size.
        rep.check("max rel err walldist vs exact ridge distance", max_rel, 0.0, 0.15)
        rep.check("mean rel err walldist vs exact ridge distance", sum_rel / n_rel, 0.0, 0.03)
        rep.check("max abs err walldist, d < 100 m [cells]", max_abs_near / hdr["dx"][2], 0.0, 0.2)

    # --- structural checks on the whole field -----------------------------
    vals = [f[fl][i][j][k] for fl in FIELDS for i in range(nx) for k in range(nz)]
    rep.check("all fields finite", 1.0 if all(math.isfinite(v) for v in vals) else 0.0, 1.0, 0.0)
    rep.check("min KE [m2/s2]", min(f["KE"][i][j][k] for i in range(nx) for k in range(nz)), 0.0, 0.0, "min")
    rep.check("min Kmv [kg/m/s]", min(f["Kmv"][i][j][k] for i in range(nx) for k in range(nz)), 0.0, 0.0, "min")
    worst = 0.0
    for i in range(nx):
        for k in range(nz):
            bound = rc.unstable_bound(rc.geom_length(f["walldist"][i][j][k], Z0, L_G_MAX))
            worst = max(worst, f["Lturb"][i][j][k] / bound)
    rep.check("max Lturb / bound (on walldist)", worst, 1.0, 1e-6, "max")

    # wall-cell k retention along the whole surface
    ratios = []
    for i in range(nx):
        rho, ke, lt, ds = f["density"][i][j][0], f["KE"][i][j][0], f["Lturb"][i][j][0], f["diss"][i][j][0]
        ratios.append((ds / (rho * rc.CMU0 ** 3 * ke ** 1.5 / lt)) ** (2.0 / 3.0))
    rep.check("max |wall cell k_start/k_end - 1|", max(abs(r - 1.0) for r in ratios), 0.0, 0.01)

    # --- physics: speed-up at the crest and the far-field log law ---------
    if mode == "physics" and not flat:
        i_top = min(range(nx), key=lambda i: abs(x_lo + (i + 0.5) * dx - xcen))
        i_up = min(range(nx), key=lambda i: abs(x_lo + (i + 0.5) * dx - (xcen - 2000.0)))
        u_top = [math.hypot(a, b) for a, b in zip(col("x_velocity", i_top), col("y_velocity", i_top))]
        u_up = [math.hypot(a, b) for a, b in zip(col("x_velocity", i_up), col("y_velocity", i_up))]
        # Jackson & Hunt (1975): fractional speed-up near the crest about 2 h/L
        for k in (0, 1, 2):
            rep.check("crest speed-up (U_top/U_up - 1), k=%d" % k, u_top[k] / u_up[k] - 1.0,
                      (0.5 * 2 * hmax / L, 2.0 * 2 * hmax / L), 0.0, "range")
        rep.check("speed-up positive in the lowest 10 cells", min(u_top[k] - u_up[k] for k in range(10)), 0.0, 0.0, "min")
        # far upstream the surface layer follows the log law (u* from the
        # AL01 wall value, k = u*^2 / Cmu0^2)
        ustar = math.sqrt(f["KE"][i_up][j][0]) * rc.CMU0
        rep.check("upstream u* from wall k [m/s]", ustar, (0.25, 0.55), 0.0, "range")
        for k in range(3):
            d = f["walldist"][i_up][j][k]
            rep.check("upstream |U|(k=%d) vs log law" % k, u_up[k], ustar / rc.KAPPA * math.log((d + Z0) / Z0), 0.15, "rel")

    rep.dump()
    return 1 if rep.failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
