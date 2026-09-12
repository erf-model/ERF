#!/usr/bin/env python3
"""Plot the isolated-building day.

    python3 plot_isolated.py ibseb_day.csv faces/day [--plotfile plt100800] [--out plots] [--dt 0.5]

  tskin_day.png        mean skin temperature by orientation and the air at the roof
  roof_budget_day.png  roof-mean SW_abs, LW_net, -H, -G through the day
  sun_path.png         solar elevation and azimuth from the CSV
  slab_roof_day.png    temperature through the roof slab (layers x time)
  tskin_slice.png, theta_slice.png   yt slices from the plotfile (needs yt)
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

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("csv"); ap.add_argument("prefix")
    ap.add_argument("--plotfile"); ap.add_argument("--out", default="plots"); ap.add_argument("--dt", type=float, default=0.5)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    files = sorted(set(re.sub(r"\.rank\d+\.csv$", "", fn) for fn in glob.glob(a.prefix + ".step*.rank*.csv")))
    steps = np.array([int(re.search(r"\.step(\d+)$", f).group(1)) for f in files]); th = (steps + 1) * a.dt / 3600.0
    sel = None; T = {}; Ta = []; roof = {k: [] for k in ("SW_abs", "LW_net", "H", "G")}; slab = []
    for f in files:
        d = load(f)
        if sel is None:
            sel = {"roof": d["dir"] == 2, "east": (d["dir"] == 0) & (d["side"] == -1), "west": (d["dir"] == 0) & (d["side"] == 1),
                   "north": (d["dir"] == 1) & (d["side"] == -1), "south": (d["dir"] == 1) & (d["side"] == 1)}
            core = sel["roof"] & (d["z_m"] > d["z_m"][sel["roof"]].max() - 1.0)
            N = 0
            while f"T_slab_{N}" in d: N += 1
            T = {k: [] for k in sel}
        for k, m in sel.items(): T[k].append(d["T_skin"][m].mean())
        Ta.append(d["T_air"][sel["roof"]].mean())
        for k in roof: roof[k].append(d[k][core].mean())
        slab.append([d[f"T_slab_{l}"][core].mean() for l in range(N)])
    fig, ax = plt.subplots(figsize=(7, 4))
    for k in sel: ax.plot(th, T[k], label=k)
    ax.plot(th, Ta, "k--", label="air at the roof")
    ax.set_xlabel("solar time [h]"); ax.set_ylabel("mean skin temperature [K]"); ax.set_xlim(0, 24); ax.set_xticks(range(0, 25, 3))
    ax.legend(fontsize=8); ax.grid(alpha=0.3); fig.tight_layout(); fig.savefig(os.path.join(a.out, "tskin_day.png"), dpi=130)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(th, roof["SW_abs"], label="SW absorbed"); ax.plot(th, roof["LW_net"], label="LW net")
    ax.plot(th, -np.array(roof["H"]), label="-H (to the air)"); ax.plot(th, -np.array(roof["G"]), label="-G (into the roof)")
    ax.set_xlabel("solar time [h]"); ax.set_ylabel("core roof flux [W/m2]"); ax.set_xlim(0, 24); ax.set_xticks(range(0, 25, 3))
    ax.legend(fontsize=8); ax.grid(alpha=0.3); fig.tight_layout(); fig.savefig(os.path.join(a.out, "roof_budget_day.png"), dpi=130)
    c = np.genfromtxt(a.csv, delimiter=",", names=True); tc = c["time_s"] / 3600.0
    fig, ax = plt.subplots(figsize=(7, 3.5)); ax2 = ax.twinx()
    ax.plot(tc, 90 - c["sun_zenith_deg"], "C1", label="elevation"); ax2.plot(tc, c["sun_azimuth_deg"], "C0", label="azimuth")
    ax.set_xlabel("solar time [h]"); ax.set_ylabel("elevation [deg]", color="C1"); ax2.set_ylabel("azimuth from north [deg]", color="C0")
    ax.set_xlim(0, 24); ax.set_xticks(range(0, 25, 3)); ax.grid(alpha=0.3); fig.tight_layout(); fig.savefig(os.path.join(a.out, "sun_path.png"), dpi=130)
    slab = np.array(slab); L = load(files[0])["thickness"][0]
    fig, ax = plt.subplots(figsize=(7, 3.5))
    im = ax.pcolormesh(th, (np.arange(N) + 0.5) * L / N * 100, slab.T, shading="nearest", cmap="inferno")
    ax.invert_yaxis(); ax.set_xlabel("solar time [h]"); ax.set_ylabel("depth below the roof skin [cm]"); ax.set_xlim(0, 24); ax.set_xticks(range(0, 25, 3))
    fig.colorbar(im, label="slab temperature [K]"); fig.tight_layout(); fig.savefig(os.path.join(a.out, "slab_roof_day.png"), dpi=130)
    print(f"wrote {a.out}/tskin_day.png roof_budget_day.png sun_path.png slab_roof_day.png")
    if a.plotfile:
        import yt
        ds = yt.load(a.plotfile)
        def center(axis, coord):
            cen = ds.domain_center.in_units("code_length").value.copy(); cen["xyz".index(axis)] = coord; return ds.arr(cen, "code_length")
        p = yt.SlicePlot(ds, "z", ("boxlib", "ibseb_tskin"), center=center("z", 25.0))
        p.set_axes_unit("code_length"); p.set_origin("native"); p.set_log(("boxlib", "ibseb_tskin"), False)
        p.set_zlim(("boxlib", "ibseb_tskin"), 295, 345); p.set_cmap(("boxlib", "ibseb_tskin"), "inferno")
        p.annotate_title("mean skin temperature of the faces touching each cell, z = 25 m, 14:00")
        p.save(os.path.join(a.out, "tskin_slice.png"))
        def _pert(field, data): return data["boxlib", "theta"] - 300.0
        ds.add_field(("gas", "theta_pert"), function=_pert, sampling_type="cell", units="")
        p = yt.SlicePlot(ds, "y", ("gas", "theta_pert"), center=center("y", 160.0))
        p.set_axes_unit("code_length"); p.set_origin("native"); p.set_log(("gas", "theta_pert"), False); p.swap_axes()
        # yt swaps the labels along with the axes, so the names are given crosswise.
        p.set_xlabel("z [m]"); p.set_ylabel("x [m]"); p.set_cmap(("gas", "theta_pert"), "RdBu_r"); p.set_zlim(("gas", "theta_pert"), -0.5, 0.5)
        p.annotate_title("theta - 300 K through the cube's centre at 14:00, wind from the left")
        p.save(os.path.join(a.out, "theta_slice.png")); print(f"wrote {a.out}/tskin_slice.png theta_slice.png")

if __name__ == "__main__":
    main()
