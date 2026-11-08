#!/usr/bin/env python3
"""Check the geostrophic wind read from erf.abl_geo_wind_table on every level.

Usage: check_geowind_profile.py <plotfile> [--f F] [--dt DT] [--tol T]

The deck (Tests/test_files/GeoWindTable_TwoLevel) starts at rest with Vg = 0 and
Ug = 2 + z/100 m/s. Over the first step the only force is the geostrophic one, so
v = f dt Ug(z) at every cell centre, on the coarse level and on the fine level,
which is refined in z as well. The neglected terms are of relative size (f dt)^2.
The check needs a level 1 with twice the coarse cells in z (otherwise the level's
own z levels and level 0's would be the same and the check could not fail).
Exit status 0 on pass, 1 on a mismatch, 2 if the plotfile cannot test it.
"""

import argparse
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, "..", "Exec", "Tools"))
import erf_plotfile  # noqa: E402


def level_values(plotfile, lev, name):
    """(k, v) for every cell of the level, from the plotfile's Level_<lev>."""
    hdr = erf_plotfile._read_header(plotfile)
    comp = hdr["names"].index(name)
    level_dir = os.path.join(plotfile, "Level_%d" % lev)
    ncomp, _, boxes, fabs = erf_plotfile._read_cell_h(level_dir)
    out = []
    for (blo, bhi), (fname, off) in zip(boxes, fabs):
        flo, fhi, data = erf_plotfile._read_fab(os.path.join(level_dir, fname), off, ncomp)
        nx, ny, nz = [fhi[d] - flo[d] + 1 for d in range(3)]
        base = comp * nx * ny * nz
        for k in range(blo[2], bhi[2] + 1):
            for j in range(blo[1], bhi[1] + 1):
                for i in range(blo[0], bhi[0] + 1):
                    q = ((k - flo[2]) * ny + (j - flo[1])) * nx + (i - flo[0])
                    out.append((k, data[base + q]))
    return out


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("plotfile")
    ap.add_argument("--f", type=float, default=1.0e-4)
    ap.add_argument("--dt", type=float, default=1.0)
    ap.add_argument("--tol", type=float, default=1.0e-6)
    a = ap.parse_args(argv)

    hdr = erf_plotfile._read_header(a.plotfile)
    if hdr["finest"] < 1:
        print("FAIL: the plotfile has no level 1")
        return 2
    with open(os.path.join(a.plotfile, "Header")) as fh:
        lines = [ln.split() for ln in fh]
    nvars = int(lines[1][0])
    # header: version, nvars, the names, then dim, time, finest level, prob_lo, prob_hi,
    # ref ratios, level domains, step numbers, and one line of cell sizes per level
    dx_line = 2 + nvars + 8
    dz = [float(lines[dx_line + lev][2]) for lev in range(2)]
    zlo = hdr["prob_lo"][2]
    nz0 = round((hdr["prob_hi"][2] - zlo) / dz[0])
    nz1 = round((hdr["prob_hi"][2] - zlo) / dz[1])
    if nz1 != 2 * nz0:
        print("FAIL: level 1 is not refined in z (%d against %d cells)" % (nz1, nz0))
        return 2

    worst = 0.0
    for lev in range(2):
        vals = level_values(a.plotfile, lev, "y_velocity")
        err = 0.0
        for k, v in vals:
            z = zlo + (k + 0.5) * dz[lev]
            expect = a.f * a.dt * (2.0 + z / 100.0)
            err = max(err, abs(v - expect) / abs(expect))
        print("level %d: %d cells, dz = %g m, max relative |v - f dt Ug(z)| = %.3e"
              % (lev, len(vals), dz[lev], err))
        worst = max(worst, err)
    if worst > a.tol:
        print("FAIL: %.3e > tol %.1e" % (worst, a.tol))
        return 1
    print("PASS: %.3e <= tol %.1e" % (worst, a.tol))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
