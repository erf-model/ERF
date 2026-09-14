#!/usr/bin/env python3
"""Numerical checks for the stable (GABLS1-style) flat-ABL k-eqn RANS case.

Usage: check_stable.py [--smoke | --physics] <plotfile> [surf_hist.dat]

--smoke (default, the CTest entry) runs the structural checks; --physics
adds the checks on the 9 h state: surface cooling, friction velocity,
low-level jet, boundary-layer depth, stable stratification and the AL01
wall value of k.
"""

import math
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, "..")]
import erf_plotfile  # noqa: E402
import rans_checks as rc  # noqa: E402

Z0 = 0.1
L_G_MAX = 30.0          # deck: max_geom_lscale (ceiling of the PBL-height cap)
UG = 8.0
THETA_SURF_0 = 265.0
COOLING = 0.25          # K/h


def main(argv):
    mode, plt, surf = rc.parse_args(argv, __doc__)
    z, p, hdr = erf_plotfile.planar_averages(plt, rc.FIELDS)
    rep = rc.Report()

    # the neutral 265 K layer below 100 m can carry a slightly negative
    # N^2 from noise, so allow the unstable factor on the length bound
    rc.structural_checks(rep, z, p, Z0, L_G_MAX, allow_unstable=True, mode=mode)

    if mode == "physics":
        t = hdr["time"]
        sh = erf_plotfile.read_surf_hist(surf)
        ustar = sh["u_star"]
        rep.check("u_star [m/s]", ustar, (0.20, 0.35), 0.0, "range")

        # surface cooling reaches the first cell
        theta_sfc = THETA_SURF_0 - COOLING * t / 3600.0
        rep.check("theta(k=0) - imposed surface theta [K]", p["theta"][0] - theta_sfc, (0.0, 1.5), 0.0, "range")

        # stable stratification through the boundary layer (dtheta/dz >= 0 below 200 m)
        dth = min((p["theta"][k + 1] - p["theta"][k]) / (z[k + 1] - z[k]) for k in range(len(z) - 1) if z[k] < 200.0)
        rep.check("min dtheta/dz below 200 m [K/m]", dth, 0.0, 1e-3, "min")

        # low-level jet: super-geostrophic maximum between 50 and 300 m
        speed = [math.hypot(u, v) for u, v in zip(p["x_velocity"], p["y_velocity"])]
        kmax = max(range(len(z)), key=lambda k: speed[k])
        rep.check("max |U| / Ug (low-level jet)", speed[kmax] / UG, 1.02, 0.0, "min")
        rep.check("height of |U| max [m]", z[kmax], (50.0, 300.0), 0.0, "range")

        # boundary-layer depth from the TKE profile (GABLS1 LES: 150-200 m)
        rep.check("BL depth from KE [m]", rc.bl_height_from_tke(z, p["KE"]), (80.0, 300.0), 0.0, "range")

        # AL01 wall value (stable: no buoyancy term)
        rep.check("KE(k=0)/u*^2 vs 1/Cmu0^2", p["KE"][0] / ustar ** 2, 1.0 / rc.CMU0 ** 2, 0.05, "rel")

        # length shortened by stability: below the neutral geometric length
        # in the stratified part (above 100 m the profile is stable from the start)
        worst = max(p["Lturb"][k] / rc.geom_length(z[k], Z0, L_G_MAX) for k in range(len(z)) if 120.0 < z[k] < 300.0)
        rep.check("max Lturb / neutral length, 120-300 m", worst, 1.0, 1e-6, "max")

    rep.dump()
    return 1 if rep.failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
