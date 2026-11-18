#!/usr/bin/env python3
"""Check the point logs of SamplePoint_LatLon.i against the column logs.

Each argument is point_log:column_log:k. Every row of the point log must be
the time and then each conserved component of cell k of the column, as the
column log writes them at the same time: the column log holds, after the time,
each of its variables (the conserved components first) over the 8 cells of the
column, bottom to top. qv = 0.001 (k+1) makes the cells of a column differ.
"""

import sys

NZ = 8


def rows(name):
    out = [[float(v) for v in r.split()] for r in open(name) if r.strip() and not r.startswith("#")]
    if not out:
        raise SystemExit(f"{name} is empty")
    return out


def main():
    for arg in sys.argv[1:]:
        point, column, k = arg.split(":")
        k = int(k)
        prows, crows = rows(point), rows(column)
        if len(prows) != len(crows):
            raise SystemExit(f"{point} has {len(prows)} rows and {column} {len(crows)}")
        for prow, crow in zip(prows, crows):
            ncons = len(prow) - 1
            if prow[0] != crow[0]:
                raise SystemExit(f"{point}: time {prow[0]} against {crow[0]} in {column}")
            if (len(crow) - 1) % NZ != 0 or (len(crow) - 1) // NZ < ncons:
                raise SystemExit(f"{column}: {len(crow) - 1} values, not a column of {NZ} cells")
            for c in range(ncons):
                want = crow[1 + c * NZ + k]
                if prow[1 + c] != want:
                    raise SystemExit(f"{point}: component {c} at t = {prow[0]} is {prow[1 + c]}, "
                                     f"but cell {k} of {column} has {want}")
        print(f"{point}: {len(prows)} rows of {ncons} components match cell {k} of {column}")


if __name__ == "__main__":
    main()
