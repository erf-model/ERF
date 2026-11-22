#!/usr/bin/env python3
"""Write input_sounding_log: the neutral case's sounding with a logarithmic wind.

Usage: make_log_sounding.py [output file]

u = (u_* / kappa) ln((z + z0) / z0), v = 0, with u_* = 0.4 m/s, kappa = 0.41 and
z0 = 0.1 m (erf.most.z0 of inputs_neutral). The levels are the cell centres of
the deck (dz = 15.625 m, 64 cells), so the interpolation of the sounding to the
mesh returns the profile itself; theta is that of input_sounding. Used by the
CTest RANS_Neutral_ABL_Flat_LogShear (check_log_shear.py).
"""

import math
import sys

USTAR, KAPPA, Z0 = 0.4, 0.41, 0.1
DZ, NZ = 15.625, 64


def theta(z):
    return 300.0 if z <= 700.0 else 300.0 + 0.03 * (z - 700.0)


def u_log(z):
    return USTAR / KAPPA * math.log((z + Z0) / Z0)


def main(argv):
    out = argv[0] if argv else "input_sounding_log"
    zs = [0.0] + [(k + 0.5) * DZ for k in range(NZ)] + [NZ * DZ]
    with open(out, "w") as fh:
        fh.write("1000.0 300.0 0.0\n")
        for z in zs:
            fh.write("%10.4f %8.3f 0.0 %.12f 0.0\n" % (z, theta(z), u_log(z)))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
