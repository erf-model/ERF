#!/usr/bin/env python3
"""The flat-fitted variant of the ridge case (prob.hmax = 1e-6): same deck,
same mesh machinery, but the Poisson wall distance must reproduce the
height above the surface to solver tolerance. Delegates to check_hill.py
with --flat."""
import sys
import check_hill

if __name__ == "__main__":
    sys.exit(check_hill.main(sys.argv[1:] + ["--flat"]))
