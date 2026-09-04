#!/usr/bin/env python3
"""Write two_buildings_5m_128x128.txt: a tall box and a short box east of it.

ERF terrain text format on 5 m nodes: nx, ny, then the x coordinates, the y
coordinates and the heights z[ix*ny + iy], one value per line. The tall box
is 80 x 80 m and 140 m high at x = 280-360 m; the short box is 40 x 40 m and
40 m high at x = 400-440 m, both centred on y = 320 m. With the sun in the
west the tall box shadows the short one's west wall to a height of
140 - 40 tan(elevation) metres.
"""
import numpy as np
n = 128; c = 5.0 * np.arange(n)
z = np.zeros((n, n))
def box(x0, x1, y0, y1, h):
    for ix, x in enumerate(c):
        for iy, y in enumerate(c):
            if x0 <= x <= x1 and y0 <= y <= y1: z[ix, iy] = max(z[ix, iy], h)
box(280.0, 360.0, 280.0, 360.0, 140.0)
box(400.0, 440.0, 300.0, 340.0, 40.0)
with open("two_buildings_5m_128x128.txt", "w") as f:
    f.write(f"{n}\n{n}\n")
    for v in c: f.write(f"{v:.3f}\n")
    for v in c: f.write(f"{v:.3f}\n")
    for ix in range(n):
        for iy in range(n): f.write(f"{z[ix, iy]:.3f}\n")
print("nodes above ground:", int((z > 0).sum()))
