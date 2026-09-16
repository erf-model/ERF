#!/usr/bin/env python3
"""Compare the turbulence kinetic energy of two runs of the convective case,
one with the implicit vertical diffusion solve and one with explicit vertical
diffusion.

Usage: check_implicit_explicit_ke.py --tol TOL <implicit plotfile> <explicit plotfile>

The buoyancy production of k is g/theta_0 times the vertical heat flux of the
theta diffusion, averaged from the two faces of the cell. Those are the full
face fluxes whether the vertical diffusion is explicit or implicit, so the
source must not depend on that choice. The two runs therefore differ only by the time
discretisation of the diffusion, and the largest planar-mean difference in KE,
relative to the largest KE of the explicit run, must stay below TOL. If the
buoyancy term depended on the diffusion of theta (for example on the face flux
scaled by the explicit fraction, which is zero with the implicit solve), the
implicit run would lose the buoyancy production above the first cell.

Guards against a vacuous pass: both runs finite, turbulence present
(max KE above 0.05 m^2/s^2), and the two runs actually different (the theta
profiles must not be identical, which they would be if both had run the same
vertical diffusion).
"""

import math
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, "..")]
import erf_plotfile  # noqa: E402
import rans_checks as rc  # noqa: E402

FIELDS = ["KE", "theta"]


def main(argv):
    tol = None
    args = []
    it = iter(argv)
    for a in it:
        if a == "--tol":
            try:
                tol = float(next(it))
            except (StopIteration, ValueError):
                print(__doc__)
                return 2
        else:
            args.append(a)
    if tol is None or len(args) != 2:
        print(__doc__)
        return 2
    plt_impl, plt_expl = args

    z, a, _ = erf_plotfile.planar_averages(plt_impl, FIELDS)
    _, b, _ = erf_plotfile.planar_averages(plt_expl, FIELDS)
    nz = len(z)
    rep = rc.Report()

    finite = all(math.isfinite(v) for p in (a, b) for f in FIELDS for v in p[f])
    rep.check("all planar-averaged fields finite", 1.0 if finite else 0.0, 1.0, 0.0)

    ke_max = max(b["KE"])
    rep.check("max KE of the explicit run [m2/s2]", ke_max, 0.05, 0.0, "min")

    dtheta = max(abs(a["theta"][k] - b["theta"][k]) for k in range(nz))
    rep.check("max |theta_impl - theta_expl| > 0 [K]", dtheta, 1e-300, 0.0, "min")

    dke = [abs(a["KE"][k] - b["KE"][k]) for k in range(nz)]
    kmax = max(range(nz), key=lambda k: dke[k])
    rel = dke[kmax] / ke_max if ke_max > 0 else float("inf")
    print("largest KE difference %.4g m2/s2 at z = %.1f m" % (dke[kmax], z[kmax]))
    rep.check("max |KE_impl - KE_expl| / max KE_expl", rel, 0.0, tol, "max")

    rep.dump()
    return 1 if rep.failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
