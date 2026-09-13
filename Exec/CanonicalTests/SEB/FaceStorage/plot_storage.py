#!/usr/bin/env python3
"""Plot the face storage of the immersed-boundary surface energy balance with yt.

    python3 plot_storage.py plt00004 [--out plots] [--z 70] [--x 320]

Makes four PNGs from one plotfile:
  mask_z.png    horizontal slice of terrain_IB_mask at height z with the
                wall-face count per cell drawn on top
  nfaces_z.png  horizontal slice of ibseb_nfaces (faces touching each cell)
  tskin_x.png   vertical slice of ibseb_tskin (mean face skin temperature)
                through x, with the mask outline
  tskin_y.png   the same through y
The fluid cells next to a wall are the only ones with a nonzero value, so a
building shows as a one-cell ring around its footprint and a cap over its
roof. Needs yt >= 4 and matplotlib.
"""
import argparse, os
import numpy as np
import yt

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("plotfile")
    ap.add_argument("--out", default="plots")
    ap.add_argument("--z", type=float, default=None, help="height [m] of the horizontal slice (default: half the tallest building)")
    ap.add_argument("--x", type=float, default=None, help="x [m] of the y-z slice (default: building centroid)")
    ap.add_argument("--y", type=float, default=None, help="y [m] of the x-z slice (default: building centroid)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    ds = yt.load(a.plotfile)
    ad = ds.all_data()
    mask = ad["boxlib", "terrain_IB_mask"].value
    solid = mask >= 0.5
    if not solid.any():
        raise SystemExit("no solid cells in terrain_IB_mask; nothing to plot")
    # ERF plotfiles carry no unit metadata: code lengths are metres.
    xs, ys, zs = (ad["boxlib", c].in_units("code_length").value for c in ("x", "y", "z"))
    xc, yc = float(xs[solid].mean()), float(ys[solid].mean())
    ztop = float(zs[solid].max())
    z_slice = a.z if a.z is not None else 0.5 * ztop
    x_slice = a.x if a.x is not None else xc
    y_slice = a.y if a.y is not None else yc
    print(f"building centroid x={xc:.1f} y={yc:.1f} m, top {ztop:.1f} m; slices at z={z_slice:.1f} x={x_slice:.1f} y={y_slice:.1f}")

    def slice_plot(axis, field, coord, fname, cmap, title, log=False, contour_mask=True, zlim=None):
        p = yt.SlicePlot(ds, axis, ("boxlib", field), center=center_for(axis, coord))
        p.set_axes_unit("code_length")
        p.set_origin("native")          # absolute coordinates, not domain-centred
        h, v = {"x": ("y", "z"), "y": ("x", "z"), "z": ("x", "y")}[axis]
        if axis == "y":
            # yt draws a y-normal slice as (z, x); swap to (x, z). The labels
            # are applied to the axes before the swap, so they go in reversed.
            p.swap_axes()
            p.set_xlabel(f"{v} [m]")
            p.set_ylabel(f"{h} [m]")
        else:
            p.set_xlabel(f"{h} [m]")
            p.set_ylabel(f"{v} [m]")
        p.set_cmap(("boxlib", field), cmap)
        p.set_log(("boxlib", field), log)
        if zlim is not None:
            p.set_zlim(("boxlib", field), *zlim)
        if contour_mask:
            p.annotate_contour(("boxlib", "terrain_IB_mask"), levels=1, factor=1, take_log=False,
                               clim=(0.5, 0.5), plot_args={"colors": "white", "linewidths": 1.0})
        p.annotate_title(title)
        p.save(os.path.join(a.out, fname))

    def center_for(axis, coord):
        c = ds.domain_center.in_units("code_length").value.copy()
        c["xyz".index(axis)] = coord
        return ds.arr(c, "code_length")

    slice_plot("z", "terrain_IB_mask", z_slice, "mask_z.png", "Greys",
               f"terrain_IB_mask at z = {z_slice:.0f} m", contour_mask=False)
    slice_plot("z", "ibseb_nfaces", z_slice, "nfaces_z.png", "viridis",
               f"ibseb_nfaces (wall faces per cell) at z = {z_slice:.0f} m")
    # Colour scale over the face values only (zero means "no face here").
    ts = ad["boxlib", "ibseb_tskin"].value
    tv = ts[ts > 0]
    tlim = (float(tv.min()) - 1.0, float(tv.max()) + 1.0) if tv.size else None
    slice_plot("x", "ibseb_tskin", x_slice, "tskin_x.png", "inferno",
               f"ibseb_tskin [K] at x = {x_slice:.0f} m", zlim=tlim)
    slice_plot("y", "ibseb_tskin", y_slice, "tskin_y.png", "inferno",
               f"ibseb_tskin [K] at y = {y_slice:.0f} m", zlim=tlim)

    nf = ad["boxlib", "ibseb_nfaces"].value
    print(f"cells with faces: {(nf > 0).sum()}  faces total: {nf.sum():.0f}  "
          f"max faces on one cell: {nf.max():.0f}")
    print(f"wrote {a.out}/mask_z.png nfaces_z.png tskin_x.png tskin_y.png")

if __name__ == "__main__":
    main()
