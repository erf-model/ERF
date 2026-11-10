#!/usr/bin/env python3
"""Check the line coordinates of a level made by a regrid after the start.

Usage: check_line_zcc.py ZCC NZ DZ LINE NROWS_COARSE

ZCC is the plt_line_lev*.zcc file of the level: it must hold one row, the heights of the NZ
cell centres of a full-depth column on a mesh of spacing DZ. LINE is the text file of a
variable of the line along that column: its first NROWS_COARSE rows, sampled before the
level was made, hold the coarse column, and every later row must hold NZ values, sampled on
the new level, so the check cannot pass if the level was never made.
"""
import sys


def read_rows(path):
    with open(path) as f:
        return [[float(x) for x in line.split()] for line in f if line.strip()]


def main(argv):
    if len(argv) != 6:
        print(__doc__)
        return 2
    zcc, nz, dz, line, ncoarse = argv[1], int(argv[2]), float(argv[3]), argv[4], int(argv[5])
    failed = 0
    try:
        rows = read_rows(zcc)
    except OSError:
        print("%s was not written" % zcc)
        return 1
    expected = [(k + 0.5) * dz for k in range(nz)]
    if len(rows) != 1 or len(rows[0]) != nz or any(
            abs(a - b) > 1.0e-5 * abs(b) for a, b in zip(rows[0], expected)):
        print("%s: expected one row of the %d cell-centre heights (k + 1/2) * %g"
              % (zcc, nz, dz))
        failed += 1
    try:
        fine = [len(r) - 1 for r in read_rows(line)[ncoarse:]]
    except OSError:
        print("%s was not written" % line)
        fine = []
    if not fine or any(n != nz for n in fine):
        print("%s: the rows after the first %d hold %s values, expected %d each"
              % (line, ncoarse, fine, nz))
        failed += 1
    print("%s: one row of %d heights; %s: %d rows on the new level"
          % (zcc, nz, line, len(fine)))
    print("FAILED" if failed else "PASSED")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
