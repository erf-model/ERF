"""Column-by-column check of the example flood engine against the rain.

The example engine keeps all the water that falls on each flood cell, so on
every flood cell its depth [m] times 1000 kg/m3 must equal the surface rain
the microphysics accumulated [kg/m2] in the atmospheric column above it, over
the same time. Flood cell (i, j) lies in column (i // C + i0, j // C + j0) of
the flood grid's level, where C = erf.flood.grid_ratio and (i0, j0) is the
column at the flood grid's lower corner (zero on level 0).

    python3 check_depth.py --flood plt_flood_00120 --rain plt2d00120 \
        [--rain0 plt2d00060] [--level 1] [--ratio 2]

--rain0 is the 2D plotfile at the flood's start, for a flood that started
after t = 0 (a restart from an atmosphere-only checkpoint). The check needs
the rain to vary over the columns (max / min > 1.5), so that a flood cell
reading the wrong column cannot pass. Relative tolerance 1e-9 (the plotfiles
are written in double precision). Uses Exec/Tools/erf_plotfile.py, which the
CTest copies next to this script.
"""

import argparse
import sys

import erf_plotfile as pf


def surface(plotfile, field, level):
    """field at k = 0 on one level, as data[i][j] over the level's domain."""
    hdr, data = pf.read_fields(plotfile, [field], level=level)
    return hdr, [[col[0] for col in row] for row in data[field]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--flood", required=True, help="the engine's flood plotfile")
    ap.add_argument("--rain", required=True, help="the 2D plotfile with precip_rain_accum")
    ap.add_argument("--rain0", help="the 2D plotfile at the flood's start")
    ap.add_argument("--level", type=int, default=0, help="the flood grid's AMR level")
    ap.add_argument("--ratio", type=int, required=True, help="erf.flood.grid_ratio")
    ap.add_argument("--rtol", type=float, default=1.0e-9)
    a = ap.parse_args()

    fhdr, depth = surface(a.flood, "flood_depth", 0)
    rhdr, rain = surface(a.rain, "precip_rain_accum", a.level)
    if a.rain0:
        _, rain0 = surface(a.rain0, "precip_rain_accum", a.level)
        rain = [[r - r0 for r, r0 in zip(row, row0)] for row, row0 in zip(rain, rain0)]

    # The flood grid's lower corner in columns of the level
    lo_l, hi_l = rhdr["domains"][a.level]
    dx_l = (rhdr["prob_hi"][0] - rhdr["prob_lo"][0]) / (hi_l[0] - lo_l[0] + 1)
    dy_l = (rhdr["prob_hi"][1] - rhdr["prob_lo"][1]) / (hi_l[1] - lo_l[1] + 1)
    i0 = int(round((fhdr["prob_lo"][0] - rhdr["prob_lo"][0]) / dx_l))
    j0 = int(round((fhdr["prob_lo"][1] - rhdr["prob_lo"][1]) / dy_l))

    C = a.ratio
    nfx, nfy = len(depth), len(depth[0])
    worst = 0.0
    used = []
    for i in range(nfx):
        for j in range(nfy):
            r = rain[i // C + i0][j // C + j0]
            used.append(r)
            if r <= 0.0:
                print("FAIL: no rain in column (%d, %d)" % (i // C + i0, j // C + j0))
                return 1
            err = abs(depth[i][j] * 1000.0 - r) / r
            if err > worst:
                worst, wi, wj = err, i, j
    spread = max(used) / min(used)
    print("flood cells %d x %d, columns from (%d, %d), rain %.6g .. %.6g kg/m2 (max/min %.3f), "
          "worst relative error %.3g at flood cell (%d, %d)"
          % (nfx, nfy, i0, j0, min(used), max(used), spread, worst,
             wi if worst > 0 else 0, wj if worst > 0 else 0))
    if spread < 1.5:
        print("FAIL: the rain hardly varies over the columns, so the check would not see a wrong column")
        return 1
    if worst > a.rtol:
        print("FAIL: depth * 1000 differs from the rain of its column by more than %g" % a.rtol)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
