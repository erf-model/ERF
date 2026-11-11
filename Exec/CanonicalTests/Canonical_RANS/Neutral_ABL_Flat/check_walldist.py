#!/usr/bin/env python3
"""Wall distance on a flat mesh whose bottom is not at z = 0.

Usage: check_walldist.py [--smoke] <plotfile>

The plotfile must hold walldist, and z_phys (cell-centre height) on a
stretched mesh; a constant-dz mesh does not plot z_phys, and its cell-centre
height is prob_lo[2] + (k + 1/2) dz. On a flat mesh the wall distance is
the height above the bottom wall, z_cc - prob_lo[2], on a constant-dz mesh
and on a stretched one alike, whatever the box layout.
The check needs prob_lo[2] > 0: with the bottom at z = 0 the height and the
wall distance are the same number and the check could not tell them apart.
Exits non-zero if a check fails.
"""

import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, ".."), os.path.join(_here, "..", "..", "..", "Tools")]
import erf_plotfile  # noqa: E402
from rans_checks import Report  # noqa: E402


def main(argv):
    args = [a for a in argv if a != "--smoke"]
    if not args:
        print(__doc__)
        return 2
    try:
        _, p, hdr = erf_plotfile.planar_averages(args[0], ["walldist", "z_phys"])
    except KeyError:
        z, p, hdr = erf_plotfile.planar_averages(args[0], ["walldist"])
        p["z_phys"] = z
    zlo = hdr["prob_lo"][2]
    rep = Report()
    rep.check("prob_lo[2] > 0 (else the check is vacuous)", 1.0 if zlo > 0.0 else 0.0, 1.0, 0.0)
    nz = len(p["walldist"])
    err = max(abs(p["walldist"][k] - (p["z_phys"][k] - zlo)) for k in range(nz))
    rep.check("max |walldist - (z_cc - prob_lo[2])| [m]", err, 0.0, 1e-8)
    rep.dump()
    return 1 if rep.failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
