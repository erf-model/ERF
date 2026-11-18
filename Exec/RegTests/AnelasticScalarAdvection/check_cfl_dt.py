#!/usr/bin/env python3
"""Check the time step erf.cfl gives on the map-factor deck.

Usage: check_cfl_dt.py <plotfile>

Steady wind U = 10, V = 5 m/s on dx = dy = 31.25 m with map factor m = 0.5, so the physical
cell is dx / m wide and dt = cfl / (m (U + V) / dx) = 2.0833 s at cfl = 0.5.  The plotfile
comes after 5 steps; its time must be 5 dt within 1 %.
"""
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, "..", "..", "Tools"))
import erf_plotfile  # noqa: E402

NSTEPS, CFL, M, U, V, DX = 5, 0.5, 0.5, 10.0, 5.0, 31.25


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    t = erf_plotfile._read_header(argv[1])["time"]
    expected = NSTEPS * CFL / (M * (U + V) / DX)
    print("t = %.5f s, expected %.5f s" % (t, expected))
    return 0 if abs(t - expected) <= 0.01 * expected else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
