# Write erf_terrain_def for inputs_terrain: a flat surface 108 m high over the whole
# 2048 m x 2048 m domain.  ERF text format (see erf.terrain_file_name in Inputs.rst):
# one value per line, nx, ny, the x values, the y values, then the heights with y
# varying fastest.  The table must cover the domain: the reader stops on a table that
# ends more than one table spacing inside it.
import numpy as np

x = np.arange(0, 2048 + 64, 64)
y = np.arange(0, 2048 + 64, 64)

with open("erf_terrain_def", "w") as target:
    target.write("%d\n%d\n" % (len(x), len(y)))
    for v in x:
        target.write("%g\n" % v)
    for v in y:
        target.write("%g\n" % v)
    for _ in x:
        for _ in y:
            target.write("%g\n" % 108)
