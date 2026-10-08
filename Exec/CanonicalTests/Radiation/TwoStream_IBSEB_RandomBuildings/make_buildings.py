#!/usr/bin/env python3
"""Write a seeded random set of box buildings for erf.buildings_file_name.

    python3 make_buildings.py                 # the case's set (seed 7)
    python3 make_buildings.py --seed 11 --n 20 --out other

Writes three files:
  <out>.txt         the nodal height map in ERF's terrain text format: nx, ny,
                    the x node coordinates, the y node coordinates, then the
                    heights z[ix*ny + iy], one value per line;
  <out>.csv         one row per building, numbered as the balance numbers
                    them (see below): footprint, height and material;
  inputs_buildings  the deck lines that depend on the set (the height map and
                    erf.ibseb.material_by_building), read by the deck through
                    FILE = inputs_buildings.

The buildings are axis-aligned boxes placed by rejection sampling inside the
region given by --region, which must lie inside the refined level with a
margin (see README.md). Footprint corners fall on multiples of --align (the
coarse cell size, so both levels resolve the same footprints) and heights
on multiples of --dz (the cell height, so every roof is a cell face). Two
buildings are never closer than --gap, so every street is resolved and the
4-connected labelling of the balance sees each box as its own building.

The balance numbers buildings in scan order over the columns (i outer, j
inner), so a building's number is its rank in (x_lo, y_lo) among the set;
the CSV and erf.ibseb.material_by_building follow that order.

The heights are sampled at the map's nodes (--node, the fine cell size).
ERF's reader interpolates the map to the nodes of the finest level and the
embedded boundary ramps from the roof to the ground over the one cell
outside each footprint, which the deck's erf.if_snap_partial_cells turns
into a half-height rim on that level (see Exec/RegTests/ImmersedForcingTest/
PartialCells).
"""
import argparse
import numpy as np


def place(a, rng, sides, heights, xlo, xhi, ylo, yhi):
    """One greedy layout: up to --tries draws, each kept if it clears the others by --gap."""
    boxes = []
    for _ in range(a.tries):
        if len(boxes) == a.n:
            break
        wx, wy = rng.choice(sides), rng.choice(sides)
        nxs = int(np.floor((xhi - xlo - wx) / a.align)) + 1
        nys = int(np.floor((yhi - ylo - wy) / a.align)) + 1
        if nxs < 1 or nys < 1:
            continue
        x0 = xlo + a.align * rng.integers(nxs)
        y0 = ylo + a.align * rng.integers(nys)
        x1, y1 = x0 + wx, y0 + wy
        # Keep the streets: the boxes grown by the gap must not overlap.
        if any(x0 < b[1] + a.gap and b[0] < x1 + a.gap and y0 < b[3] + a.gap and b[2] < y1 + a.gap
               for b in boxes):
            continue
        h = rng.choice(heights)
        m = int(rng.integers(1, a.materials + 1))
        boxes.append((x0, x1, y0, y1, h, m))
    return boxes


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, default=7, help="random seed (default 7)")
    ap.add_argument("--n", type=int, default=8, help="number of buildings (default 8)")
    ap.add_argument("--extent", type=float, nargs=2, default=[640.0, 640.0],
                    help="domain extent in x and y [m] (default 640 640)")
    ap.add_argument("--region", type=float, nargs=4, default=[200.0, 440.0, 200.0, 440.0],
                    metavar=("XLO", "XHI", "YLO", "YHI"),
                    help="area the footprints must lie in [m] (default 200 440 200 440)")
    ap.add_argument("--node", type=float, default=10.0, help="node spacing of the map [m] (default 10)")
    ap.add_argument("--align", type=float, default=20.0, help="footprint corners on multiples of this [m] (default 20)")
    ap.add_argument("--dz", type=float, default=10.0, help="heights on multiples of this [m] (default 10)")
    ap.add_argument("--side", type=float, nargs=2, default=[20.0, 60.0], help="footprint side range [m] (default 20 60)")
    ap.add_argument("--height", type=float, nargs=2, default=[10.0, 60.0], help="height range [m] (default 10 60)")
    ap.add_argument("--gap", type=float, default=60.0, help="smallest street width [m] (default 60)")
    ap.add_argument("--materials", type=int, default=3, help="material ids 1..N to draw from (default 3)")
    ap.add_argument("--out", default="buildings", help="output prefix (default buildings)")
    ap.add_argument("--tries", type=int, default=2000, help="draws per layout (default 2000)")
    ap.add_argument("--layouts", type=int, default=500, help="layouts to try before giving up (default 500)")
    a = ap.parse_args()

    xlo, xhi, ylo, yhi = a.region
    if not (0.0 <= xlo < xhi <= a.extent[0] and 0.0 <= ylo < yhi <= a.extent[1]):
        ap.error("--region must lie inside --extent")
    for v, name in ((a.node, "--node"), (a.align, "--align"), (a.dz, "--dz")):
        if v <= 0.0:
            ap.error(f"{name} must be positive")
    if a.align % a.node != 0.0:
        ap.error("--align must be a multiple of --node, so the footprint corners are map nodes")

    rng = np.random.default_rng(a.seed)
    sides = np.arange(np.ceil(a.side[0] / a.align), np.floor(a.side[1] / a.align) + 1) * a.align
    heights = np.arange(np.ceil(a.height[0] / a.dz), np.floor(a.height[1] / a.dz) + 1) * a.dz
    if sides.size == 0 or heights.size == 0:
        ap.error("--side and --height leave no size on the --align / --dz multiples")

    # Greedy placement jams once the first boxes are badly placed, so a layout
    # that falls short is thrown away and a new one drawn from the same
    # generator: still one set per seed.
    for _ in range(a.layouts):
        boxes = place(a, rng, sides, heights, xlo, xhi, ylo, yhi)
        if len(boxes) == a.n:
            break
    if len(boxes) < a.n:
        raise SystemExit(f"no layout of {a.n} buildings in {a.layouts} tries; "
                         "widen --region or lower --n, --side or --gap")

    # The balance's numbering: scan order over the columns, i outer.
    boxes.sort(key=lambda b: (b[0], b[2]))

    nx = int(round(a.extent[0] / a.node)) + 1
    ny = int(round(a.extent[1] / a.node)) + 1
    xs = a.node * np.arange(nx)
    ys = a.node * np.arange(ny)
    z = np.zeros((nx, ny))
    # A building covers the nodes of its footprint, its edges included; the
    # surface then ramps to the ground over the cell outside it.
    for x0, x1, y0, y1, h, _ in boxes:
        ix = (xs >= x0 - 1e-6) & (xs <= x1 + 1e-6)
        iy = (ys >= y0 - 1e-6) & (ys <= y1 + 1e-6)
        z[np.ix_(ix, iy)] = h

    with open(a.out + ".txt", "w") as f:
        f.write(f"{nx}\n{ny}\n")
        for v in xs:
            f.write(f"{v:.3f}\n")
        for v in ys:
            f.write(f"{v:.3f}\n")
        for v in z.ravel():
            f.write(f"{v:.3f}\n")
    with open(a.out + ".csv", "w") as f:
        f.write("building,x_lo_m,x_hi_m,y_lo_m,y_hi_m,height_m,material\n")
        for n, (x0, x1, y0, y1, h, m) in enumerate(boxes, start=1):
            f.write(f"{n},{x0:g},{x1:g},{y0:g},{y1:g},{h:g},{m}\n")
    with open("inputs_buildings", "w") as f:
        f.write(f"# Written by make_buildings.py --seed {a.seed} --n {a.n}; do not edit by hand.\n")
        f.write(f"erf.buildings_file_name = {a.out}.txt\n")
        f.write("erf.ibseb.material_by_building = " + " ".join(str(b[5]) for b in boxes) + "\n")

    area = sum((b[1] - b[0]) * (b[3] - b[2]) for b in boxes)
    print(f"seed {a.seed}: {len(boxes)} buildings, plan area {area:g} m2 "
          f"({100.0 * area / ((xhi - xlo) * (yhi - ylo)):.1f} % of the region), "
          f"heights {min(b[4] for b in boxes):g}-{max(b[4] for b in boxes):g} m")


if __name__ == "__main__":
    main()
