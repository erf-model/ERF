#!/usr/bin/env python3
"""Plot the phase 3 view fractions and longwave on the building faces with yt.

    python3 plot_longwave.py plt00002 [--out plots] [--z 20] [--y 320]

  f_sky_y.png    x-z centreline slice of ibseb_f_sky: 1 on the tall roof,
                 about one half on open walls, less on walls facing the other box
  lw_net_y.png   x-z centreline slice of ibseb_lw_net (positive into the face)
  f_sky_z.png    horizontal slice of ibseb_f_sky through both boxes
Needs yt >= 4 and matplotlib.
"""
import argparse, os
import yt

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("plotfile"); ap.add_argument("--out", default="plots")
    ap.add_argument("--z", type=float, default=20.0); ap.add_argument("--y", type=float, default=320.0)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    ds = yt.load(a.plotfile)
    def center(axis, coord):
        c = ds.domain_center.in_units("code_length").value.copy(); c["xyz".index(axis)] = coord
        return ds.arr(c, "code_length")
    def slc(axis, field, coord, fname, cmap, title, zlim=None):
        p = yt.SlicePlot(ds, axis, ("boxlib", field), center=center(axis, coord))
        p.set_axes_unit("code_length"); p.set_origin("native")
        h, v = {"x": ("y", "z"), "y": ("x", "z"), "z": ("x", "y")}[axis]
        if axis == "y":
            p.swap_axes(); p.set_xlabel(f"{v} [m]"); p.set_ylabel(f"{h} [m]")
        else:
            p.set_xlabel(f"{h} [m]"); p.set_ylabel(f"{v} [m]")
        p.set_cmap(("boxlib", field), cmap); p.set_log(("boxlib", field), False)
        if zlim: p.set_zlim(("boxlib", field), *zlim)
        p.annotate_contour(("boxlib", "terrain_IB_mask"), levels=1, factor=1, take_log=False,
                           clim=(0.5, 0.5), plot_args={"colors": "white", "linewidths": 1.0})
        p.annotate_title(title); p.save(os.path.join(a.out, fname))
    slc("y", "ibseb_f_sky", a.y, "f_sky_y.png", "viridis", f"sky view fraction at y = {a.y:.0f} m", zlim=(0, 1))
    slc("y", "ibseb_lw_net", a.y, "lw_net_y.png", "coolwarm", f"net longwave [W/m2] at y = {a.y:.0f} m")
    slc("z", "ibseb_f_sky", a.z, "f_sky_z.png", "viridis", f"sky view fraction at z = {a.z:.0f} m", zlim=(0, 1))
    print(f"wrote {a.out}/f_sky_y.png lw_net_y.png f_sky_z.png")

if __name__ == "__main__":
    main()
