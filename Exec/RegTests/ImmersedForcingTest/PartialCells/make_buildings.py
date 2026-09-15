#!/usr/bin/env python3
"""Write cube_40m_10m_32x32.txt: one 40 m cube at the centre of a 320 m domain.

ERF terrain text format on 10 m nodes: nx, ny, then the x coordinates, the y
coordinates and the heights z[ix*ny + iy], one value per line. The cube
spans x, y = 140-180 m and is 40 m tall, so it is four cells on a side on
the 10 m grid (the embedded-boundary reader steps each edge over one cell,
giving a 4-cell core with a half-height rim; the checks read the geometry
from the face dump).
"""
import numpy as np
n = 32; c = 10.0 * np.arange(n)
z = np.zeros((n, n))
for ix, x in enumerate(c):
    for iy, y in enumerate(c):
        if 140.0 <= x <= 180.0 and 140.0 <= y <= 180.0: z[ix, iy] = 40.0
with open("cube_40m_10m_32x32.txt", "w") as f:
    f.write(f"{n}\n{n}\n")
    for v in c: f.write(f"{v:.3f}\n")
    for v in c: f.write(f"{v:.3f}\n")
    for ix in range(n):
        for iy in range(n): f.write(f"{z[ix, iy]:.3f}\n")
print("nodes above ground:", int((z > 0).sum()))
