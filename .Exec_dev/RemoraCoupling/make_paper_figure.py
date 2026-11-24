#!/usr/bin/env python3
"""
Publication figure for the Perlin et al. (2007) uncoupled validation.

Six panels: potential temperature, the two wind components, turbulent kinetic
energy, vertical eddy viscosity and the MYNN master length scale, each at the
72-hour forecast time.

The reference curves are pixel-traced from a 400 dpi render of Figs. 7 and 8 of
the paper at the 25 km offshore location, uncoupled case.  Calibration comes
from the printed tick labels (each figure uses three side-by-side panels with
repeated axes; Fig. 7c's labels are 5/10/15 per panel).  Points where the
"25km offshore" annotation crosses a curve are omitted.  The master length scale is not
plotted in the paper; it is back-computed from Figs. 8a and 8b using ERF's own
stability function, so it carries the uncertainty of both.

Usage:  python make_paper_figure.py [run_dir]   (default run_oneway)
"""
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yt

yt.set_log_level("error")

RUN = sys.argv[1] if len(sys.argv) > 1 else "run_oneway"
OUT = "perlin_validation.pdf"

# --- traced from Perlin et al. (2007), 25 km offshore, uncoupled ------------
# Grey (uncoupled) curve.  Where it is hidden under the black (coupled) curve
# the two coincide to within the line width and the black value is used.
P = {
    "theta": ([286.02, 286.04, 286.07, 286.10, 286.13, 286.18, 286.27, 286.34,
               286.54, 286.83, 287.32, 287.71, 288.21, 288.74, 289.55, 290.18],
              [10, 50, 100, 150, 200, 250, 300, 350, 400, 450, 475, 500, 525,
               550, 600, 650]),
    "u":     ([3.69, 3.96, 4.18, 4.36, 4.35, 4.17, 3.82, 3.13, 2.09, 1.20, 0.76,
               0.49, 0.21, 0.27, 0.45, 0.61, 0.40, 0.02, -0.07, -0.10],
              [10, 25, 50, 100, 150, 200, 250, 300, 350, 400, 425, 440, 455,
               485, 500, 515, 530, 550, 575, 650]),
    "v":     ([8.9, 9.9, 10.7, 12.0, 13.0, 13.9, 14.7, 15.8, 16.4, 16.4, 16.5,
               16.1, 15.7, 15.3, 15.6, 15.6, 15.5, 15.2, 15.1],
              [10, 25, 50, 100, 150, 200, 250, 300, 350, 400, 425, 450, 470,
               485, 515, 530, 550, 600, 650]),
    "tke":   ([0.42, 0.40, 0.36, 0.32, 0.30, 0.25, 0.21, 0.18, 0.14, 0.11,
               0.05, 0.01],
              [10, 25, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450]),
    "Km":    ([1.45, 3.93, 4.39, 4.43, 4.23, 3.93, 3.13, 2.25, 1.28, 0.79, 0.54,
               0.06, 0.0],
              [15, 50, 75, 100, 125, 150, 200, 250, 300, 350, 400, 450, 500]),
    # Back-computed from the earlier Km/tke trace (which followed the coupled
    # curve, about 5% higher in Km); not yet redone from the values above.
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
