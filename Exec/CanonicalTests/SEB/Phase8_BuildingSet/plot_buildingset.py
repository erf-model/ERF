#!/usr/bin/env python3
"""Plot the building-set morning.

    python3 plot_buildingset.py ibseb_set.csv faces/set [--plotfile plt43200] [--out plots] [--dt 0.5]

  tskin_buildings.png   mean skin temperature of each building against solar time
  shadow_buildings.png  shadow fraction of each building (mutual shading)
  faces_map.png         the faces at the last dump coloured by skin temperature,
                        seen from above (walls drawn at their positions)
  tskin_slice.png       yt slice of ibseb_tskin at 25 m from the plotfile (needs yt)
"""
import argparse, glob, os, re
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

def load(prefix):
    rows = []
    for fn in sorted(glob.glob(prefix + ".rank*.csv")):
        with open(fn) as f:
            hdr = f.readline().strip().split(",")
            for line in f: rows.append([float(v) for v in line.strip().split(",")])
    a = np.array(rows); return {h: a[:, n] for n, h in enumerate(hdr)}

NAMES = {1: "slab 60 m concrete", 2: "north block 20 m timber", 3: "cube 40 m brick", 4: "far block 20 m timber"}

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("csv"); ap.add_argument("prefix")
    ap.add_argument("--plotfile"); ap.add_argument("--out", default="plots"); ap.add_argument("--dt", type=float, default=0.5)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    c = np.genfromtxt(a.csv, delimiter=",", names=True)
    fig, ax = plt.subplots(figsize=(7, 4))
    for b in sorted(set(c["building"].astype(int))):
        m = c["building"] == b; ax.plot(5 + c["time_s"][m] / 3600, c["T_skin_mean_K"][m], label=NAMES.get(b, f"building {b}"))
    ax.set_xlabel("solar time [h]"); ax.set_ylabel("mean skin temperature [K]"); ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "tskin_buildings.png"), dpi=130)
    fig, ax = plt.subplots(figsize=(7, 3.5))
    for b in sorted(set(c["building"].astype(int))):
        m = c["building"] == b; ax.plot(5 + c["time_s"][m] / 3600, c["shadow_frac"][m], label=NAMES.get(b, f"building {b}"))
    ax.set_xlabel("solar time [h]"); ax.set_ylabel("shadow fraction of the faces"); ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "shadow_buildings.png"), dpi=130)
    files = sorted(set(re.sub(r"\.rank\d+\.csv$", "", fn) for fn in glob.glob(a.prefix + ".step*.rank*.csv")))
    d = load(files[-1]); step = int(re.search(r"step(\d+)", files[-1]).group(1))
    fig, ax = plt.subplots(figsize=(6.5, 6))
    sc = ax.scatter(d["x_m"], d["y_m"], c=d["T_skin"], s=14, cmap="inferno", marker="s")
    ax.set_aspect("equal"); ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
    ax.set_title(f"face skin temperature at {5 + (step + 1) * a.dt / 3600:.1f} h solar time (all faces, from above)")
    fig.colorbar(sc, label="T_skin [K]"); fig.tight_layout(); fig.savefig(os.path.join(a.out, "faces_map.png"), dpi=130)
    print(f"wrote {a.out}/tskin_buildings.png shadow_buildings.png faces_map.png")
    if a.plotfile:
        import yt
        ds = yt.load(a.plotfile)
        cen = ds.domain_center.in_units("code_length").value.copy(); cen[2] = 25.0
        p = yt.SlicePlot(ds, "z", ("boxlib", "ibseb_tskin"), center=ds.arr(cen, "code_length"))
        p.set_axes_unit("code_length"); p.set_origin("native"); p.set_log(("boxlib", "ibseb_tskin"), False)
        p.set_zlim(("boxlib", "ibseb_tskin"), 295, 335); p.set_cmap(("boxlib", "ibseb_tskin"), "inferno")
        p.annotate_title("mean skin temperature of the faces touching each cell, z = 25 m")
        p.save(os.path.join(a.out, "tskin_slice.png")); print(f"wrote {a.out}/tskin_slice.png")

if __name__ == "__main__":
    main()
