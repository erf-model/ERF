#!/usr/bin/env python3
"""Plot the phase 6 prognostic balance from the face dumps, the per-building
CSV and a plotfile.

    python3 plot_prognostic.py faces_solar ibseb_solar.csv [--plotfile plt05400] [--out plots] [--dt 1.0]

  tskin_orientation.png  mean skin temperature of the roof and the four walls
                         against time, from the step-tagged dumps
  budget.png             building-mean SW_abs, LW_net, -H, -G against time
                         from the CSV (the four sum to the residual, zero)
  tskin_slice.png        yt slice of ibseb_tskin through the cube (needs yt)
Reads <prefix>.step*.rank*.csv; needs matplotlib.
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
    ap = argparse.ArgumentParser(); ap.add_argument("prefix"); ap.add_argument("csv")
    ap.add_argument("--plotfile"); ap.add_argument("--out", default="plots"); ap.add_argument("--dt", type=float, default=1.0)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    files = sorted(set(re.sub(r"\.rank\d+\.csv$", "", fn) for fn in glob.glob(a.prefix + ".step*.rank*.csv")))
    steps = np.array([int(re.search(r"\.step(\d+)$", f).group(1)) for f in files])
    sel = {"roof": lambda d: d["dir"] == 2, "east": lambda d: (d["dir"] == 0) & (d["side"] == -1),
           "west": lambda d: (d["dir"] == 0) & (d["side"] == 1), "north": lambda d: (d["dir"] == 1) & (d["side"] == -1),
           "south": lambda d: (d["dir"] == 1) & (d["side"] == 1)}
    series = {k: [] for k in sel}
    for f in files:
        d = load(f)
        for k, fn in sel.items(): series[k].append(d["T_skin"][fn(d)].mean())
    t = (steps + 1) * a.dt / 60.0
    fig, ax = plt.subplots(figsize=(6.5, 4))
    for k in sel: ax.plot(t, series[k], label=k)
    ax.set_xlabel("time [min]"); ax.set_ylabel("mean skin temperature [K]"); ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "tskin_orientation.png"), dpi=130)
    print(f"wrote {a.out}/tskin_orientation.png")

    c = np.genfromtxt(a.csv, delimiter=",", names=True)
    tt = c["time_s"] / 60.0
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(tt, c["SW_abs_mean_Wm2"], label="SW absorbed"); ax.plot(tt, c["LW_net_mean_Wm2"], label="LW net")
    ax.plot(tt, -c["H_mean_Wm2"], label="-H (to the air)"); ax.plot(tt, -c["G_mean_Wm2"], label="-G (into the wall)")
    ax.plot(tt, c["SW_abs_mean_Wm2"] + c["LW_net_mean_Wm2"] - c["H_mean_Wm2"] - c["G_mean_Wm2"], "k--", label="sum")
    ax.set_xlabel("time [min]"); ax.set_ylabel("building-mean flux [W/m2]"); ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "budget.png"), dpi=130)
    print(f"wrote {a.out}/budget.png")

    if a.plotfile:
        import yt
        ds = yt.load(a.plotfile)
        cen = ds.domain_center.in_units("code_length").value.copy(); cen[2] = 25.0
        p = yt.SlicePlot(ds, "z", ("boxlib", "ibseb_tskin"), center=ds.arr(cen, "code_length"))
        p.set_axes_unit("code_length"); p.set_origin("native"); p.set_log(("boxlib", "ibseb_tskin"), False)
        p.set_zlim(("boxlib", "ibseb_tskin"), 295, None); p.set_cmap(("boxlib", "ibseb_tskin"), "inferno")
        p.annotate_title("mean skin temperature of the faces touching each cell, z = 25 m")
        p.save(os.path.join(a.out, "tskin_slice.png")); print(f"wrote {a.out}/tskin_slice.png")

if __name__ == "__main__":
    main()
