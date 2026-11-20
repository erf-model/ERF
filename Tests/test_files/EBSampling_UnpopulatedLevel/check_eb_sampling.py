#!/usr/bin/env python3
"""Check that line and plane sampling blank the cells covered by EB terrain.

Usage: check_eb_sampling.py LOG NZ NCOVERED PLANE LINE

LOG is an erf.sample_line log of a full-depth column of NZ cells through the hill. On every
row, after the time, every variable over the column must be erf.sampling_covered_value
(-9999) in its lowest NCOVERED cells, which lie inside the hill, and in none of the cells
above them.

PLANE is a plane plotfile sampled with density, theta and volfrac. Where volfrac is 0 the
plane lies in covered cells, and density and theta must be -9999; where volfrac is 1 both
neighbouring cells are fluid, and neither may be -9999. volfrac itself is never blanked.
The plane must meet both, so the check cannot pass on a plane that misses the hill.

LINE is the name of an erf.do_line_sampling line along the same column, written as text with
volfrac, density and theta. density and theta must be blanked in the same NCOVERED cells,
volfrac must be 0 there and positive above. The heights of level 0 must have been written
(plt_line_lev0.zcc), and none for level 1, which does not exist.
"""
import os
import sys

sys.path.insert(0, os.environ.get("ERF_TOOLS", os.path.join(os.path.dirname(
    os.path.abspath(__file__)), "..", "..", "..", "Exec", "Tools")))
import erf_plotfile  # noqa: E402

COVERED = -9999.0


def is_covered(v):
    return abs(v - COVERED) < 1.0e-6 * abs(COVERED)


def check_log(log, nz, ncovered):
    failed = 0
    with open(log) as f:
        rows = [[float(x) for x in line.split()] for line in f if line.strip()]
    if not rows:
        print("%s: no rows" % log)
        return 1
    for row in rows:
        vals = row[1:]
        if len(vals) == 0 or len(vals) % nz != 0:
            print("%s: a row has %d values, not a multiple of %d" % (log, len(vals), nz))
            return 1
        nvar = len(vals) // nz
        for iv in range(nvar):
            col = vals[iv * nz:(iv + 1) * nz]
            bad_low = [k for k in range(ncovered) if not is_covered(col[k])]
            bad_high = [k for k in range(ncovered, nz) if is_covered(col[k])]
            if bad_low or bad_high:
                print("%s, t = %g, variable %d: not blanked at k = %s, blanked at k = %s"
                      % (log, row[0], iv, bad_low, bad_high))
                failed += 1
    print("%s: %d rows of %d variables over %d cells, lowest %d blanked"
          % (log, len(rows), len(rows[0][1:]) // nz, nz, ncovered))
    return failed


def read_rows(path):
    with open(path) as f:
        return [[float(x) for x in line.split()] for line in f if line.strip()]


def check_line(line, nz, ncovered):
    failed = 0
    for name in ("density", "theta"):
        for row in read_rows(line + "." + name):
            col = row[1:]
            if len(col) != nz:
                print("%s.%s: a row has %d values, expected %d" % (line, name, len(col), nz))
                return failed + 1
            if (not all(is_covered(v) for v in col[:ncovered]) or
                    any(is_covered(v) for v in col[ncovered:])):
                print("%s.%s, t = %g: blanked cells are not the lowest %d"
                      % (line, name, row[0], ncovered))
                failed += 1
    for row in read_rows(line + ".volfrac"):
        col = row[1:]
        if any(v != 0.0 for v in col[:ncovered]) or any(v <= 0.0 for v in col[ncovered:]):
            print("%s.volfrac, t = %g: not 0 in the lowest %d cells and positive above"
                  % (line, row[0], ncovered))
            failed += 1
    if not os.path.exists("plt_line_lev0.zcc"):
        print("plt_line_lev0.zcc was not written")
        failed += 1
    if os.path.exists("plt_line_lev1.zcc"):
        print("plt_line_lev1.zcc was written for a level that does not exist")
        failed += 1
    print("%s: density, theta and volfrac checked" % line)
    return failed


def check_plane(plane):
    hdr, data = erf_plotfile.read_fields(plane, ["density", "theta", "volfrac"])
    failed = 0
    n_cov = n_fluid = 0
    vf = data["volfrac"]
    for i in range(len(vf)):
        for j in range(len(vf[i])):
            for k in range(len(vf[i][j])):
                v = vf[i][j][k]
                if is_covered(v):
                    print("%s: volfrac blanked at (%d,%d,%d)" % (plane, i, j, k))
                    failed += 1
                    continue
                for name in ("density", "theta"):
                    blanked = is_covered(data[name][i][j][k])
                    if v == 0.0 and not blanked:
                        print("%s: %s not blanked at covered (%d,%d,%d)" % (plane, name, i, j, k))
                        failed += 1
                    if v == 1.0 and blanked:
                        print("%s: %s blanked at fluid (%d,%d,%d)" % (plane, name, i, j, k))
                        failed += 1
                n_cov += v == 0.0
                n_fluid += v == 1.0
    if n_cov == 0 or n_fluid == 0:
        print("%s: %d covered and %d fluid samples; the plane must meet both"
              % (plane, n_cov, n_fluid))
        failed += 1
    print("%s: %d covered and %d fluid samples" % (plane, n_cov, n_fluid))
    return failed


def main(argv):
    if len(argv) != 6:
        print(__doc__)
        return 2
    nz, ncovered = int(argv[2]), int(argv[3])
    failed = (check_log(argv[1], nz, ncovered) + check_plane(argv[4]) +
              check_line(argv[5], nz, ncovered))
    print("FAILED" if failed else "PASSED")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
