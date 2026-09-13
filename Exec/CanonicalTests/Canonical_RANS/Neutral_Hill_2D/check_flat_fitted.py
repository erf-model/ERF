#!/usr/bin/env python3
"""The flat-fitted variant of the ridge case (prob.hmax = 1e-6): same deck,
same mesh machinery, but the wall distance must reproduce the height above
the surface to solver tolerance -- with either wall_dist_type, the deck's
terrain_height default or the Poisson solve of the _Poisson CTest variant.
Delegates to check_hill.py with --flat."""
import sys
import check_hill

if __name__ == "__main__":
    sys.exit(check_hill.main(sys.argv[1:] + ["--flat"]))
