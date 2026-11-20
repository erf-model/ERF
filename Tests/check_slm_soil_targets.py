#!/usr/bin/env python3
"""Check that a refined SLM level nudges its soil toward the same targets as level 0.

Usage: check_slm_soil_targets.py <plt_lsm plotfile> [--tol T]

SLM nudges the soil toward the targets soilt_obs and soilw_obs, set from the soil state at the
start of the run; soilt_nudge and soilw_nudge are the nudging tendencies, proportional to the
soil's departure from them.  A refined level made from level 0 partway through the run, or
remade by a regrid, must keep level 0's targets (issue 5178): its tendencies are then close to
level 0's in every column.  Had it taken the soil state at the time it was made as its targets,
its tendencies would be far smaller in the columns made then.

For soilt_nudge and soilw_nudge, the tendency of each column (that of its most nudged layer) on
level 1 is compared with the mean over level 0, and must be within a fraction --tol of it.  The
check stops if level 0 is not nudged at all, since it could then not tell the targets apart.
Exit status 0 on pass, 1 on a mismatch, 2 when the check is vacuous.
"""

import argparse
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, "..", "Exec", "Tools"))
import erf_plotfile  # noqa: E402

FIELDS = ("soilt_nudge", "soilw_nudge")


def column_tendencies(plotfile, level):
    """The tendency of each column of a level, that of its most nudged layer, per field."""
    hdr, data = erf_plotfile.read_fields(plotfile, list(FIELDS), level=level)
    lo, hi = hdr["domains"][level]
    nz = hi[2] - lo[2] + 1
    out = {n: [] for n in FIELDS}
    for blo, bhi in erf_plotfile.level_boxes(plotfile, level):
        for i in range(blo[0], bhi[0] + 1):
            for j in range(blo[1], bhi[1] + 1):
                for n in FIELDS:
                    col = data[n][i - lo[0]][j - lo[1]]
                    out[n].append(max((col[k] for k in range(nz)), key=abs))
    return out


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("plotfile")
    ap.add_argument("--tol", type=float, default=0.1)
    a = ap.parse_args(argv)

    if erf_plotfile.finest_level(a.plotfile) < 1:
        print("FAIL %s has no level 1" % a.plotfile)
        return 1
    coarse = column_tendencies(a.plotfile, 0)
    fine = column_tendencies(a.plotfile, 1)

    status = 0
    for n in FIELDS:
        ref = sum(coarse[n]) / len(coarse[n])
        if ref == 0.0:
            print("VACUOUS %s: level 0 is not nudged" % n)
            return 2
        ratios = [v / ref for v in fine[n]]
        bad = [r for r in ratios if abs(r - 1.0) > a.tol]
        print("%s %s: level 1 / level 0 between %.4f and %.4f over %d columns"
              % ("FAIL" if bad else "ok  ", n, min(ratios), max(ratios), len(ratios)))
        if bad:
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
