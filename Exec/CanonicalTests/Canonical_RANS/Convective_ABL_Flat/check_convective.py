#!/usr/bin/env python3
"""Numerical checks for the convective flat-ABL k-eqn RANS case.

Usage: check_convective.py [--smoke | --physics] <plotfile> [surf_hist.dat]

--smoke (default, the CTest entry) runs the structural checks; --physics
adds the checks on the 4 h state: the column heat budget against the
imposed surface flux, a well-mixed layer, the inversion height, the
mixed-layer warming and the AL01 wall value of k with its buoyancy term.
"""

import math
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, "..")]
import erf_plotfile  # noqa: E402
import rans_checks as rc  # noqa: E402

Z0 = 0.16
L_G_MAX = 100.0         # deck: max_geom_lscale (ceiling of the PBL-height cap)
SURF_FLUX = 0.24        # K m/s (kinematic)
SOUNDING = [(0.0, 300.0), (937.0, 300.0), (1062.0, 308.0), (2062.0, 311.0)]


def theta_init(zc):
    for (z1, t1), (z2, t2) in zip(SOUNDING[:-1], SOUNDING[1:]):
        if z1 <= zc <= z2:
            return t1 + (t2 - t1) * (zc - z1) / (z2 - z1)
    return SOUNDING[-1][1] + 0.003 * (zc - SOUNDING[-1][0])


def main(argv):
    mode, plt, surf = rc.parse_args(argv, __doc__)
    z, p, hdr = erf_plotfile.planar_averages(plt, rc.FIELDS)
    nz = len(z)
    dz = z[1] - z[0]
    rep = rc.Report()

    rc.structural_checks(rep, z, p, Z0, L_G_MAX, allow_unstable=True, mode=mode)

    if mode == "physics":
        t = hdr["time"]
        sh = erf_plotfile.read_surf_hist(surf)
        ustar = sh["u_star"]
        rep.check("u_star [m/s]", ustar, (0.30, 0.80), 0.0, "range")

        # column heat budget: sum rho (theta - theta_init) dz = rho_sfc F t
        gained = sum(p["density"][k] * (p["theta"][k] - theta_init(z[k])) * dz for k in range(nz))
        rep.check("column heat gain / (rho_sfc F t)", gained / (p["density"][0] * SURF_FLUX * t), 1.0, 0.10)

        # inversion height: strongest dtheta/dz
        grad = [(p["theta"][k + 1] - p["theta"][k]) / dz for k in range(nz - 1)]
        kinv = max(range(nz - 1), key=lambda k: grad[k])
        zi = 0.5 * (z[kinv] + z[kinv + 1])
        rep.check("inversion height [m]", zi, (900.0, 1250.0), 0.0, "range")

        # mixed layer between 0.2 zi and 0.7 zi. A local-K closure carries a
        # superadiabatic lapse -F(z)/K_h through the mixed layer (about
        # -2 K/km here for K_h near 35 m2/s), where LES gives under 0.3 K of
        # spread; 2 K bounds the local-closure value without hiding a
        # broken profile.
        ml = [p["theta"][k] for k in range(nz) if 0.2 * zi < z[k] < 0.7 * zi]
        mean = sum(ml) / len(ml)
        spread = max(ml) - min(ml)
        rep.check("theta spread in 0.2-0.7 zi [K]", spread, 0.0, 2.0)
        grad = [(p["theta"][k + 1] - p["theta"][k]) / dz for k in range(nz - 1) if 0.2 * zi < z[k] < 0.7 * zi]
        rep.check("max dtheta/dz in 0.2-0.7 zi [K/m]", max(grad), 0.0, 1e-3, "max")

        # mixed-layer warming against the encroachment estimate F t / zi
        rep.check("mixed-layer warming / (F t / zi)", (mean - 300.0) / (SURF_FLUX * t / zi), 1.0, 0.30)

        # AL01 wall value with the destabilising term: k >= u*^2 / Cmu0^2
        rep.check("KE(k=0)/u*^2 >= 1/Cmu0^2", p["KE"][0] / ustar ** 2, 1.0 / rc.CMU0 ** 2 * 0.99, 0.0, "min")

        # convective TKE: positive through the mixed layer and small above zi
        ke_ml = min(p["KE"][k] for k in range(nz) if 0.1 * zi < z[k] < 0.8 * zi)
        rep.check("min KE in 0.1-0.8 zi [m2/s2]", ke_ml, 0.05, 0.0, "min")
        ke_above = max(p["KE"][k] for k in range(nz) if z[k] > 1.3 * zi)
        rep.check("max KE above 1.3 zi [m2/s2]", ke_above, 0.05, 0.0, "max")

    rep.dump()
    return 1 if rep.failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
