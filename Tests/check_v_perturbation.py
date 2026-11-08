#!/usr/bin/env python3
"""Check whether the initial y-velocity perturbation of the Poiseuille channel is present.

Usage: check_v_perturbation.py <plotfile> --expect {zero,nonzero} [--min M]

The Couette/Poiseuille setup has no mean y-velocity, so y_velocity in the initial
plotfile is the perturbation alone, of amplitude prob.pert_delta_v. With --expect zero
every value must be 0; with --expect nonzero the largest |y_velocity| must be at least
--min (default 0.5 for pert_delta_v = 2). Exit status 0 on pass, 1 on fail.
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
    ap.add_argument("--expect", choices=("zero", "nonzero"), required=True)
    ap.add_argument("--min", type=float, default=0.5)
    a = ap.parse_args(argv)
    _, d = erf_plotfile.read_fields(a.plotfile, ["y_velocity"])
    v = d["y_velocity"]
    vmax = max(abs(x) for plane in v for col in plane for x in col)
    print("max |y_velocity| = %.6e (expected %s)" % (vmax, a.expect))
    ok = (vmax == 0.0) if a.expect == "zero" else (vmax >= a.min)
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
