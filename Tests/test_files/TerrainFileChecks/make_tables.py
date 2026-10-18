#!/usr/bin/env python3
"""Write the height tables of the TerrainFileChecks CTests (ERF text format:
nx, ny, the nx x values, the ny y values, then nx*ny heights with y varying
fastest, one value per line).

The domain of TerrainFileChecks.i is 800 m x 400 m, so a table whose x and y
lists are swapped does not cover it.  The ridge runs north-south (along y), so
its slope is all in x.

  ridge_middle.txt     ridge at x = 400 m: gentle next to the inflow and outflow faces
  ridge_outflow.txt    ridge at x = 700 m: |dh/dx| up to about 0.9 next to the outflow face
  swapped_lists.txt    ridge_middle with nx/ny and the x/y lists swapped
  x_not_increasing.txt ridge_middle with two x values exchanged
  too_few_values.txt   ridge_middle without its last height
  extra_values.txt     ridge_middle with one height too many
  two_per_line.txt     ridge_middle with two heights on one line
  short_of_domain.txt  ridge_middle with the x list ending at 700 m, 100 m short of the domain
                       face and more than its 43.75 m spacing
"""
import math

X = [50.0 * i for i in range(17)]   # 0 .. 800
Y = [50.0 * j for j in range(9)]    # 0 .. 400


def ridge(xc):
    return [[120.0 * math.exp(-((x - xc) / 120.0) ** 2) for _ in Y] for x in X]


def write(name, nx, ny, xs, ys, heights):
    with open(name, "w") as f:
        for v in [nx, ny]:
            f.write(f"{v}\n")
        for v in list(xs) + list(ys) + list(heights):
            f.write(f"{v:.6f}\n" if isinstance(v, float) else f"{v}\n")


def flat(h):
    return [h[i][j] for i in range(len(X)) for j in range(len(Y))]   # y fastest


mid = flat(ridge(400.0))
write("ridge_middle.txt", len(X), len(Y), X, Y, mid)
write("ridge_outflow.txt", len(X), len(Y), X, Y, flat(ridge(700.0)))
write("swapped_lists.txt", len(Y), len(X), Y, X, mid)
xs = list(X); xs[5], xs[6] = xs[6], xs[5]
write("x_not_increasing.txt", len(X), len(Y), xs, Y, mid)
write("too_few_values.txt", len(X), len(Y), X, Y, mid[:-1])
write("extra_values.txt", len(X), len(Y), X, Y, mid + [0.0])
write("short_of_domain.txt", len(X), len(Y), [v * 700.0 / 800.0 for v in X], Y, mid)
with open("two_per_line.txt", "w") as f:
    lines = [str(len(X)), str(len(Y))] + [f"{v:.6f}" for v in X + Y + mid]
    lines[40] = lines[40] + " " + lines[41]
    del lines[41]
    f.write("\n".join(lines) + "\n")
