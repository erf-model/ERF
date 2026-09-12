#!/usr/bin/env python3
"""Plot the phase 4 wall heat flux and the warm wake with yt.

    python3 plot_sensible.py plt00080 [--out plots] [--z 22.5] [--y 320]

  H_y.png       x-z centreline slice of ibseb_H, the sensible flux out of the
                faces (largest on the windward wall, smallest in the lee)
  theta_z.png   horizontal slice of the potential temperature perturbation
                (theta - 300 K) through the cube: the warm wake downstream
  theta_y.png   the same on the x-z centreline
Needs yt >= 4 and matplotlib.
"""
import argparse, os
import numpy as np
import yt

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("plotfile"); ap.add_argument("--out", default="plots")
    ap.add_argument("--z", type=float, default=22.5); ap.add_argument("--y", type=float, default=320.0)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    ds = yt.load(a.plotfile)
    def _pert(field, data): return data["boxlib", "theta"] - 300.0
    ds.add_field(("boxlib", "theta_pert"), function=_pert, sampling_type="cell", units="")
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
    slc("y", "ibseb_H", a.y, "H_y.png", "inferno", f"sensible heat flux from the faces [W/m2] at y = {a.y:.0f} m")
    slc("z", "theta_pert", a.z, "theta_z.png", "RdBu_r", f"theta - 300 K at z = {a.z:.1f} m", zlim=(-0.5, 0.5))
    slc("y", "theta_pert", a.y, "theta_y.png", "RdBu_r", f"theta - 300 K at y = {a.y:.0f} m", zlim=(-0.5, 0.5))
    print(f"wrote {a.out}/H_y.png theta_z.png theta_y.png")

if __name__ == "__main__":
    main()
