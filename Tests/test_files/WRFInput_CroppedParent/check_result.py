#!/usr/bin/env python3
"""Check that a cropped parent patch does not shift its WRF grandchild, in i, j or k."""

import argparse
from pathlib import Path
import re
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("plotfile_reader", type=Path)
    # The lowest level of the level-2 grids: 0 by default, higher for a nest that does not
    # start at the bottom of the domain
    parser.add_argument("--klo", type=int, default=0)
    parser.add_argument("--log", default="run.log")
    args = parser.parse_args()

    sys.path.insert(0, str(args.plotfile_reader))
    from erf_plotfile import _read_cell_h, _read_fab  # noqa: E402

    log = (args.run_dir / args.log).read_text()
    match = re.search(r"WRFInput subdomain 0 at level 2 is\s+\(\((\d+),(\d+),(\d+)\)\s+"
                      r"\((\d+),(\d+),(\d+)\)", log)
    if not match or tuple(map(int, match.groups())) != (78, 78, 0, 95, 95, 7):
        raise AssertionError("level-2 WRF file was not placed at global box 78..95: " +
                             (match.group(0) if match else "no box printed"))

    plotfile = args.run_dir / "plt00000"
    header = (plotfile / "Header").read_text().splitlines()
    names = header[2:2 + int(header[1])]
    if "theta" not in names:
        raise AssertionError("initial plotfile is missing theta")
    theta_component = names.index("theta")
    level_dir = plotfile / "Level_2"
    ncomp, _, boxes, fabs = _read_cell_h(str(level_dir))
    if not boxes:
        raise AssertionError("level 2 has no plotfile boxes")

    cells_checked = 0
    klo_found = min(lo[2] for lo, _ in boxes)
    if klo_found != args.klo:
        raise AssertionError(f"the level-2 grids start at k = {klo_found}, expected {args.klo}")
    for (lo, hi), (filename, offset) in zip(boxes, fabs):
        if lo[0] < 78 or lo[1] < 78 or hi[0] > 95 or hi[1] > 95 or hi[2] > 7:
            raise AssertionError(f"level-2 plotfile box {lo}..{hi} exceeds d03")
        flo, fhi, data = _read_fab(str(level_dir / filename), offset, ncomp)
        nx, ny, nz = [fhi[d] - flo[d] + 1 for d in range(3)]
        base = theta_component * nx * ny * nz
        for k in range(lo[2], hi[2] + 1):
            for j in range(lo[1], hi[1] + 1):
                for i in range(lo[0], hi[0] + 1):
                    index = base + ((k - flo[2]) * ny + j - flo[1]) * nx + i - flo[0]
                    expected = 300.0 + 0.01 * i + 0.001 * j + 0.1 * k
                    actual = data[index]
                    if abs(actual - expected) > 1.0e-4:
                        raise AssertionError(f"theta at level-2 ({i},{j},{k}) = {actual}, "
                                             f"expected {expected}")
                    cells_checked += 1
    expected_cells = 18 * 18 * (8 - args.klo)
    if cells_checked != expected_cells:
        raise AssertionError(f"checked {cells_checked} level-2 cells, expected {expected_cells}")
    print(f"checked level-2 box and theta in {cells_checked} cells")

if __name__ == "__main__":
    main()
