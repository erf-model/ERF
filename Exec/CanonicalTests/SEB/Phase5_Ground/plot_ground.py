#!/usr/bin/env python3
"""Plot the phase 5 slab: the flux into the walls per face and the slab
temperature profile of one face over time, from the face dumps.

    python3 plot_ground.py faces_thin [--out plots]

  G_faces.png     G of every face against its height, coloured by building
Reads the per-rank dump(s) <prefix>.rank*.csv; needs matplotlib.
"""
import argparse, glob, os
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
    ap = argparse.ArgumentParser(); ap.add_argument("prefix"); ap.add_argument("--out", default="plots")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    d = load(a.prefix)
    fig, ax = plt.subplots(figsize=(6, 4))
    for b in np.unique(d["bid"]):
        m = d["bid"] == b
        ax.scatter(d["G"][m], d["z_m"][m], s=8, label=f"building {int(b)}")
    ax.set_xlabel("G into the wall [W/m2]"); ax.set_ylabel("face height [m]"); ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "G_faces.png"), dpi=130)
    print(f"wrote {a.out}/G_faces.png")

if __name__ == "__main__":
    main()
