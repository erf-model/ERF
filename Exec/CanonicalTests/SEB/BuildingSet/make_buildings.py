#!/usr/bin/env python3
"""Write buildings_10m_48x48.txt: a small building set on 10 m nodes.

ERF terrain text format: nx, ny, then the x coordinates, the y coordinates
and the heights z[ix*ny + iy], one value per line. Four buildings on a
480 m domain, arranged so they shade one another through the morning
(numbered as the balance numbers them, in scan order over (i, j)):
  1. a 60 m slab, 30 x 60 m, at x = 150-180, y = 200-260 (the tall one);
  2. a 20 m block, 40 x 40 m, at x = 150-190, y = 300-340, north of the slab;
  3. a 40 m cube, 40 x 40 m, at x = 230-270, y = 210-250, east of the slab,
     in its morning shadow;
  4. a 20 m block, 40 x 40 m, at x = 320-360, y = 150-190, off on its own.
The reader ramps every edge over one cell, so the set has sliver cells at
every corner: this case runs with erf.if_snap_partial_cells = true.
"""
import numpy as np
n = 48; c = 10.0 * np.arange(n)
z = np.zeros((n, n))
def box(x0, x1, y0, y1, h):
    for ix, x in enumerate(c):
        for iy, y in enumerate(c):
            if x0 <= x <= x1 and y0 <= y <= y1: z[ix, iy] = max(z[ix, iy], h)
box(150, 180, 200, 260, 60.0)
box(230, 270, 210, 250, 40.0)
box(150, 190, 300, 340, 20.0)
box(320, 360, 150, 190, 20.0)
with open("buildings_10m_48x48.txt", "w") as f:
    f.write(f"{n}\n{n}\n")
    for v in c: f.write(f"{v:.3f}\n")
    for v in c: f.write(f"{v:.3f}\n")
    for ix in range(n):
        for iy in range(n): f.write(f"{z[ix, iy]:.3f}\n")
print("nodes above ground:", int((z > 0).sum()))
