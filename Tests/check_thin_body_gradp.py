#!/usr/bin/env python3
"""Check that the thin-body projection stores a pressure gradient consistent with
its pressure increment.

Usage: check_thin_body_gradp.py <plotfile> [--tol T]

project_velocity_tb adds each corrector pass's pressure increment phi to pp_inc
(plotted as pert_pres when erf.anelastic = 1) and its face gradient to gradp
(plotted, averaged to cell centres, as dpdx, dpdy, dpdz). The average of the two
face gradients around a cell is the centred difference of the cell values, so in
every interior cell

    dpdx(i,j,k) = (pert_pres(i+1,j,k) - pert_pres(i-1,j,k)) / (2 dx)

and the same in y and z. The check reports the largest mismatch over the interior
cells, relative to the largest |dpdx|, |dpdy|, |dpdz|, and fails above --tol.

It also stops when the gradient is zero everywhere: with nothing to compare, the
check could not fail. Exit status 0 on pass, 1 on a mismatch, 2 on a vacuous run.
"""

import argparse
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, "..", "Exec", "Tools"))
import erf_plotfile  # noqa: E402


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("plotfile")
    ap.add_argument("--tol", type=float, default=1e-8)
    a = ap.parse_args(argv)

    names = ["pert_pres", "dpdx", "dpdy", "dpdz"]
    hdr, d = erf_plotfile.read_fields(a.plotfile, names)
    lo, hi, dx = hdr["lo"], hdr["hi"], hdr["dx"]
    n = [hi[i] - lo[i] + 1 for i in range(3)]
    p = d["pert_pres"]

    gmax = 0.0
    err = [0.0, 0.0, 0.0]
    for i in range(1, n[0] - 1):
        for j in range(1, n[1] - 1):
            for k in range(1, n[2] - 1):
                cd = [(p[i + 1][j][k] - p[i - 1][j][k]) / (2.0 * dx[0]),
                      (p[i][j + 1][k] - p[i][j - 1][k]) / (2.0 * dx[1]),
                      (p[i][j][k + 1] - p[i][j][k - 1]) / (2.0 * dx[2])]
                for c, f in enumerate(("dpdx", "dpdy", "dpdz")):
                    g = d[f][i][j][k]
                    gmax = max(gmax, abs(g))
                    err[c] = max(err[c], abs(g - cd[c]))

    print("largest |dp/dx_i| in the interior: %.6e" % gmax)
    if gmax == 0.0:
        print("FAIL: the pressure gradient is zero everywhere; the check is vacuous")
        return 2
    rel = [e / gmax for e in err]
    for c, f in enumerate(("dpdx", "dpdy", "dpdz")):
        print("%s vs centred difference of pert_pres: max error %.3e (relative %.3e)"
              % (f, err[c], rel[c]))
    worst = max(rel)
    if worst > a.tol:
        print("FAIL: relative mismatch %.3e > tol %.1e" % (worst, a.tol))
        return 1
    print("PASS: relative mismatch %.3e <= tol %.1e" % (worst, a.tol))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
