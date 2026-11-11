#!/usr/bin/env python3
"""Check 2D sampled earth-relative winds against the Lambert-conformal rotation.

Usage: check_2d_earth_winds.py <plt2d> <plt3d> --k K --truelat1 T1 --truelat2 T2
                               --stand-lon L0 [--tol-rad R] [--min-alpha DEG]

The 2D plotfile must hold lon_m, u_east_k_<K> and v_north_k_<K>; the 3D plotfile
x_velocity and y_velocity (cell centred). On a Lambert-conformal grid the grid x
axis is at

    alpha = -n (lon - L0),   n = ln(cos T1 / cos T2) / ln(tan(45 - T1/2) / tan(45 - T2/2))

degrees counterclockwise from east (n = sin T1 when T1 = T2), which is the angle
WRF stores as COSALPHA/SINALPHA. ERF measures alpha from its lat/lon arrays
instead, so this is an independent reference. Every cell must satisfy

    |(u_east, v_north) - R(alpha) (x_velocity, y_velocity)| <= R * |wind| + 1e-6

with R = --tol-rad (default 0.003 rad, 0.17 deg: ERF's one-sided difference at
the first and last i columns differs from the closed form by up to 0.12 deg).
The check fails if the grid turns by less than --min-alpha degrees anywhere
(default 5), since then an unrotated field could pass. Standard library only.
Exit status 0 on pass, 1 on fail.
"""

import argparse
import math
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, "..", "Exec", "CanonicalTests", "Canonical_RANS"))
import erf_plotfile  # noqa: E402


def cone_factor(t1, t2):
    r = math.radians
    if abs(t1 - t2) < 1.0e-9:
        return math.sin(r(abs(t1)))
    return (math.log(math.cos(r(t1)) / math.cos(r(t2))) /
            math.log(math.tan(r(45.0 - abs(t1) / 2.0)) / math.tan(r(45.0 - abs(t2) / 2.0))))


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("plt2d")
    ap.add_argument("plt3d")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--truelat1", type=float, required=True)
    ap.add_argument("--truelat2", type=float, required=True)
    ap.add_argument("--stand-lon", type=float, required=True)
    ap.add_argument("--tol-rad", type=float, default=3.0e-3)
    ap.add_argument("--min-alpha", type=float, default=5.0)
    a = ap.parse_args(argv)

    ue_name, vn_name = "u_east_k_%d" % a.k, "v_north_k_%d" % a.k
    _, f2 = erf_plotfile.read_fields(a.plt2d, ["lon_m", ue_name, vn_name])
    _, f3 = erf_plotfile.read_fields(a.plt3d, ["x_velocity", "y_velocity"])
    n = cone_factor(a.truelat1, a.truelat2)

    lon, u_e, v_n = f2["lon_m"], f2[ue_name], f2[vn_name]
    u, v = f3["x_velocity"], f3["y_velocity"]
    worst = 0.0
    worst_unrotated = 0.0
    max_alpha = 0.0
    ncell = 0
    for i in range(len(lon)):
        for j in range(len(lon[i])):
            dlon = lon[i][j][0] - a.stand_lon
            dlon = (dlon + 180.0) % 360.0 - 180.0
            alpha = math.radians(-n * dlon)
            max_alpha = max(max_alpha, abs(math.degrees(alpha)))
            ug, vg = u[i][j][a.k], v[i][j][a.k]
            ue_ref = ug * math.cos(alpha) - vg * math.sin(alpha)
            vn_ref = ug * math.sin(alpha) + vg * math.cos(alpha)
            err = math.hypot(u_e[i][j][0] - ue_ref, v_n[i][j][0] - vn_ref)
            scale = math.hypot(ug, vg)
            worst = max(worst, (err - 1.0e-6) / max(scale, 1.0e-30))
            worst_unrotated = max(worst_unrotated,
                                  math.hypot(u_e[i][j][0] - ug, v_n[i][j][0] - vg))
            ncell += 1

    print("cells %d, cone factor %.6f, max |alpha| %.3f deg, worst error / |wind| %.3e rad, "
          "max |earth - grid wind| %.4f m/s" % (ncell, n, max_alpha, worst, worst_unrotated))
    ok = True
    if max_alpha < a.min_alpha:
        print("FAIL: the grid turns by at most %.3f deg; the test cannot see the rotation" % max_alpha)
        ok = False
    if worst > a.tol_rad:
        print("FAIL: worst error / |wind| %.3e > %.1e" % (worst, a.tol_rad))
        ok = False
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
