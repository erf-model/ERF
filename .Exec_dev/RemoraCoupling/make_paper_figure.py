#!/usr/bin/env python3
"""
Publication figure for the Perlin et al. (2007) uncoupled validation.

Six panels: potential temperature, the two wind components, turbulent kinetic
energy, vertical eddy viscosity and the MYNN master length scale, each at the
72-hour forecast time.

The reference curves are pixel-traced from a 400 dpi render of Figs. 7 and 8 of
the paper at the 25 km offshore location, uncoupled case.  Calibration comes
from the printed tick labels (each figure uses three side-by-side panels with
repeated axes).  The two points of the southward-wind curve that the "25km
offshore" annotation crosses are omitted.  The master length scale is not
plotted in the paper; it is back-computed from Figs. 8a and 8b using ERF's own
stability function, so it carries the uncertainty of both.

Usage:  python make_paper_figure.py [run_dir]   (default run_taper1.35)
"""
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yt

yt.set_log_level("error")

RUN = sys.argv[1] if len(sys.argv) > 1 else "run_taper1.35"
OUT = "perlin_validation.pdf"

# --- traced from Perlin et al. (2007), 25 km offshore, uncoupled ------------
P = {
    "theta": ([286.12, 286.14, 286.18, 286.20, 286.23, 286.36, 286.63,
               286.93, 287.45, 288.26, 289.30],
              [10, 50, 100, 150, 200, 300, 400, 450, 500, 550, 600]),
    "u":     ([2.03, 3.98, 4.17, 4.17, 3.06, 1.20, 0.20],
              [10, 50, 100, 200, 300, 400, 500]),
    "v":     ([8.4, 11.9, 13.4, 14.9, 15.3, 14.7],
              [10, 100, 200, 300, 500, 600]),
    "tke":   ([0.409, 0.369, 0.334, 0.306, 0.257, 0.216, 0.152, 0.118,
               0.077, 0.015],
              [25, 50, 75, 100, 150, 200, 300, 350, 400, 450]),
    "Km":    ([4.37, 4.45, 4.71, 4.58, 4.35, 3.60, 1.80, 0.86, 0.64, 0.41, 0.11],
              [50, 75, 100, 125, 150, 200, 300, 350, 400, 450, 500]),
    "Lm":    ([12.1, 11.8, 10.5, 6.0, 3.2, 3.0], [100, 150, 200, 300, 350, 400]),
}

PANELS = [
    ("theta", r"$\theta$  [K]",                    (285.8, 290.5)),
    ("u",     r"$u$  [m s$^{-1}$]",                (0, 5)),
    ("v",     r"$-v$  [m s$^{-1}$]",               (6, 17)),
    ("tke",   r"$e$  [m$^2$ s$^{-2}$]",            (0, 0.5)),
    ("Km",    r"$K_m$  [m$^2$ s$^{-1}$]",          (0, 6)),
    ("Lm",    r"$L_m$  [m]",                       (0, 25)),
]


def load(run):
    ds = yt.load(f"{run}/plt86400")
    cg = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)

    def f(n):
        return np.asarray(cg[("boxlib", n)].d)[0, 0, :].astype(float)

    return {"z": f("z_phys"), "theta": f("theta"), "u": f("x_velocity"),
            "v": -f("y_velocity"), "tke": f("KE"), "Km": f("nut"),
            "Lm": f("Lturb")}


def main():
    d = load(RUN)
    plt.rcParams.update({
        "font.size": 9, "axes.labelsize": 9, "legend.fontsize": 8,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "axes.linewidth": 0.7, "xtick.direction": "in", "ytick.direction": "in",
    })
    fig, axes = plt.subplots(2, 3, figsize=(7.1, 5.4), sharey=True)
    for ax, (key, xlabel, xlim) in zip(axes.ravel(), PANELS):
        ax.plot(P[key][0], P[key][1], "s", ms=4, mfc="none", mec="0.15",
                mew=0.9, label="Perlin et al. (2007)", zorder=3,
                linestyle="none")
        ax.plot(d[key], d["z"], "-", lw=1.6, color="#0d6b78", label="ERF",
                zorder=2)
        ax.set_xlabel(xlabel)
        ax.set_xlim(*xlim)
        ax.set_ylim(0, 650)
        ax.grid(alpha=0.25, lw=0.5)
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$z$  [m]")
    # legend in the K_m panel, whose upper right is empty
    axes[1, 1].legend(loc="upper right", frameon=False, handlelength=1.6)
    fig.tight_layout(pad=0.4, w_pad=0.8, h_pad=0.9)
    fig.savefig(OUT, bbox_inches="tight")
    fig.savefig(OUT.replace(".pdf", ".png"), dpi=200, bbox_inches="tight")
    print("wrote", OUT, "and", OUT.replace(".pdf", ".png"))

    # agreement summary quoted in the text
    print("\nagreement at the traced levels:")
    for key in ("Km", "tke"):
        v, z = P[key]
        e = np.array([d[key][int(np.argmin(abs(d["z"] - zz)))] for zz in z])
        v = np.array(v)
        m = v > 0.05 * v.max()
        print(f"  {key:5s} max |rel err| = {np.max(abs((e[m]-v[m])/v[m])):.1%}")


if __name__ == "__main__":
    main()
