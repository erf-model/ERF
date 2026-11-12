#!/usr/bin/env python3
"""Checks for the terrain-following inflow profile cases.

Usage: check_inflow_profile.py [--smoke] <plotfile>

The interior starts from a uniform 10 m/s, so the first column next to the
inflow face (i = 0) carries the log law only through the boundary condition.
The height above the ground, the domain top and each column's level heights
come from z_phys (exact for basic terrain following). The checks:

- the wind speed in the first column, above the wall cell and below the
  Rayleigh layer, follows the log law at the height above the local ground
  (median within 1 %, 90th percentile within 5 %);
- in the lowest 200 m it follows that law at least twice as closely as the log
  law at the height above the domain floor (a profile applied at absolute
  heights), and more closely than the log law at the flat-mesh level height (a
  profile applied by level, as xlo.dirichlet_file is);
- the wind blows from 270 degrees;
- tke in the first column is 0.8 to 1.25 of the profile's;
- every field is finite, and the speed and tke stay bounded over the domain.

xlo.dirichlet_file with the same log law fails the median, by-level and tke
checks on the plateau, incline and cross ridge: see README.md.
"""

import math
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, "..")]
import erf_plotfile  # noqa: E402
import rans_checks as rc  # noqa: E402

SPEED, HREF, Z0, TKE_ZSCALE = 10.0, 10.0, 0.1, 700.0
USTAR = rc.KAPPA * SPEED / math.log((HREF + Z0) / Z0)
TKE0 = USTAR ** 2 / rc.CMU0 ** 2
RAYLEIGH_DEPTH = 200.0
LOW = 200.0
FIELDS = ["x_velocity", "y_velocity", "z_velocity", "KE", "z_phys"]


def log_law(z):
    return USTAR / rc.KAPPA * math.log((max(z, 0.0) + Z0) / Z0)


def tke_profile(z):
    return TKE0 * max((USTAR * TKE_ZSCALE - z) / (max(USTAR, 0.01) * TKE_ZSCALE), 0.01)


def rel(s, z):
    return abs(s - log_law(z)) / log_law(z)


def median(vals):
    s = sorted(vals)
    n = len(s)
    return 0.5 * (s[(n - 1) // 2] + s[n // 2])


def percentile(vals, q):
    s = sorted(vals)
    return s[min(len(s) - 1, max(0, int(round(q * (len(s) - 1)))))]


def main(argv):
    args = [a for a in argv[1:] if not a.startswith("--")]
    if len(args) != 1:
        print(__doc__)
        return 2
    _, d = erf_plotfile.read_fields(args[0], FIELDS)
    u, v, w, ke, zp = (d[f] for f in FIELDS)
    nx, ny, nz = len(u), len(u[0]), len(u[0][0])
    rep = rc.Report()

    finite = all(math.isfinite(f[i][j][k]) for f in (u, v, w, ke)
                 for i in range(nx) for j in range(ny) for k in range(nz))
    rep.check("all fields finite", 1.0 if finite else 0.0, 1.0, 0.0)
    if not finite:
        rep.dump()
        return 1
    top = zp[0][0][nz - 1] + 0.5 * (zp[0][0][nz - 1] - zp[0][0][nz - 2])
    umax = max(math.hypot(u[i][j][k], v[i][j][k]) for i in range(nx) for j in range(ny) for k in range(nz))
    wmax = max(abs(w[i][j][k]) for i in range(nx) for j in range(ny) for k in range(nz))
    kmax = max(ke[i][j][k] for i in range(nx) for j in range(ny) for k in range(nz))
    rep.check("max horizontal speed [m/s]", umax, 2.0 * log_law(top), 0.0, "max")
    rep.check("max |w| [m/s]", wmax, SPEED, 0.0, "max")
    rep.check("max tke [m2/s2]", kmax, 10.0 * TKE0, 0.0, "max")

    err, low_agl, low_floor, low_level, ke_ratio, vs, grounds = [], [], [], [], [], [], []
    i = 0
    for j in range(ny):
        col = zp[i][j]
        ground = col[0] - 0.5 * (col[1] - col[0])
        grounds.append(ground)
        for k in range(1, nz):
            agl = col[k] - ground
            if agl > top - ground - RAYLEIGH_DEPTH:
                break
            s = math.hypot(u[i][j][k], v[i][j][k])
            err.append(rel(s, agl))
            if agl <= LOW:
                low_agl.append(rel(s, agl))
                low_floor.append(rel(s, col[k]))
                low_level.append(rel(s, (k + 0.5) * top / nz))
            ke_ratio.append(ke[i][j][k] / tke_profile(agl))
            vs.append(v[i][j][k])

    print("first column: ground %.1f to %.1f m, %d cells checked, %d in the lowest %.0f m"
          % (min(grounds), max(grounds), len(err), len(low_agl), LOW))
    rep.check("first column: median rel error, log law above ground", median(err), 0.0, 0.01, "max")
    rep.check("first column: 90th pct rel error, log law above ground", percentile(err, 0.9), 0.0, 0.05, "max")
    rep.check("lowest 200 m: error above ground / above floor",
              median(low_agl) / max(median(low_floor), 1e-300), 0.0, 0.5, "max")
    rep.check("lowest 200 m: error above ground / by level",
              median(low_agl) / max(median(low_level), 1e-300), 0.0, 1.0, "max")
    rep.check("first column: mean v [m/s]", sum(vs) / len(vs), 0.0, 0.5, "abs")
    rep.check("first column: median tke / profile tke", median(ke_ratio), (0.8, 1.25), 0.0, "range")
    rep.dump()
    return 1 if rep.failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
