#!/usr/bin/env python3
"""Check vertical lines sampled by height on a terrain-fitted mesh against the plotfile.

Usage: check_line_heights.py PLOTFILE AGL PREFIX LINE:i:j:zlo:zhi [LINE:i:j:zlo:zhi ...]

For each line, PREFIX+LINE.z (the height of each sample, last row) must equal z_phys of
column (i,j) over consecutive cells, and those cells must run from the cell holding zlo to the
cell holding zhi, each clamped to the column.  With AGL = 1, zlo and zhi are heights above the
ground.  A LINE:x:j:z spec is a line in x through row j: it must be the cells of every column
of the row at the level k of cells that holds z (or z above the ground) in column (0,j).

The faces are estimated from the cell centres: the deck's levels start at 20 m and are
stretched by 1.2, so the two lowest centres are 10 s and 32 s above the ground (s is the
squeeze of the column), and the faces above are midway between centres.
"""
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
for d in (os.path.join(_here, "..", "..", "..", "Exec", "Tools"), os.environ.get("ERF_TOOLS", "")):
    if d:
        sys.path.insert(0, d)
import erf_plotfile  # noqa: E402

TOL = 1.0e-4


def column_ground(col):
    """The ground under a column, from its two lowest cell centres (see the docstring)."""
    return col[0] - (col[1] - col[0]) * 10.0 / 22.0


def holds_height(col, k, h, ground):
    """Whether cell k of a column holds height h (clamped to the column); faces are estimated."""
    nz = len(col)
    below = ground if k == 0 else 0.5 * (col[k - 1] + col[k])
    above = 2.0 * col[k] - below if k == nz - 1 else 0.5 * (col[k] + col[k + 1])
    return below <= h < above or (k == 0 and h < below) or (k == nz - 1 and h >= above)


def check_row(fields, zp, nz, agl, prefix):
    name, _, j, z = fields[0], fields[1], int(fields[2]), float(fields[3])
    col0 = [zp[0][j][k] for k in range(nz)]
    target = z + column_ground(col0) if agl else z
    with open(prefix + name + ".z") as f:
        rows = [[float(x) for x in line.split()] for line in f if line.strip()]
    zs = rows[-1][1:]
    k0 = min(range(nz), key=lambda k: abs(col0[k] - zs[0]))
    row = [zp[i][j][k0] for i in range(len(zp))]
    match = len(zs) == len(row) and all(abs(a - b) <= TOL * max(1.0, abs(b)) for a, b in zip(zs, row))
    right = match and holds_height(col0, k0, target, column_ground(col0))
    print("%s: row j = %d, wanted the cells holding %.1f m at column (0,%d), sampled k = %d "
          "from %.1f to %.1f m (max %.1f): heights %s, level %s"
          % (prefix + name, j, target, j, k0, zs[0], zs[-1], max(zs),
             "match z_phys" if match else "DIFFER from z_phys", "right" if right else "WRONG"))
    return (not match) + (not right)


def main(argv):
    if len(argv) < 5:
        print(__doc__)
        return 2
    plotfile, agl, prefix = argv[1], argv[2] == "1", argv[3]
    _, data = erf_plotfile.read_fields(plotfile, ["z_phys"])
    zp = data["z_phys"]
    nz = len(zp[0][0])
    failed = 0
    for spec in argv[4:]:
        fields = spec.split(":")
        if fields[1] == "x":
            failed += check_row(fields, zp, nz, agl, prefix)
            continue
        name, i, j, zlo, zhi = fields
        i, j, zlo, zhi = int(i), int(j), float(zlo), float(zhi)
        col = [zp[i][j][k] for k in range(nz)]
        ground = col[0] - (col[1] - col[0]) * 10.0 / 22.0
        lo, hi = (zlo + ground, zhi + ground) if agl else (zlo, zhi)

        with open(prefix + name + ".z") as f:
            rows = [[float(x) for x in line.split()] for line in f if line.strip()]
        z = rows[-1][1:]
        k0 = min(range(nz), key=lambda k: abs(col[k] - z[0]))
        ks = list(range(k0, k0 + len(z)))
        match = ks[-1] < nz and all(abs(z[m] - col[k]) <= TOL * max(1.0, abs(col[k]))
                                    for m, k in enumerate(ks))

        def holds(k, h):
            below = ground if k == 0 else 0.5 * (col[k - 1] + col[k])
            above = 2.0 * col[k] - below if k == nz - 1 else 0.5 * (col[k] + col[k + 1])
            return below <= h < above or (k == 0 and h < below) or (k == nz - 1 and h >= above)

        in_range = match and holds(ks[0], lo) and holds(ks[-1], hi)
        print("%s: column (%d,%d), ground %.1f m, wanted %.1f to %.1f m, sampled k = %d..%d "
              "at %.1f to %.1f m: heights %s, cells %s"
              % (prefix + name, i, j, ground, lo, hi, ks[0], ks[-1], z[0], z[-1],
                 "match z_phys" if match else "DIFFER from z_phys",
                 "right" if in_range else "WRONG"))
        failed += (not match) + (not in_range)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
