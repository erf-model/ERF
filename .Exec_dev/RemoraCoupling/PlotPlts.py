#!/usr/bin/env python3
"""
Plot x_velocity, y_velocity, and Kmv versus z_phys from one or more ERF plotfiles.

Usage examples:
    python PlotPlts.py plt00500
    python PlotPlts.py plt00100 plt00300 plt00500
    python PlotPlts.py plt00500 --mode column --i 32 --j 32
    python PlotPlts.py plt00100 plt00500 --mode scatter --zmax 3000 -o prof.png

When several plotfiles are given, their curves are overlaid on the same axes
with one color per plotfile and a legend identifying each.

Modes:
    mean    : horizontal (i,j) average at each k index (default)
    column  : a single vertical column at index (i, j)
    scatter : every cell, u/v/Kmv vs its own z_phys (shows terrain/spread)

Requires: yt, numpy, matplotlib
"""
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import yt

yt.set_log_level("error")

FIELDS = ["x_velocity", "y_velocity", "Kmv"]
LABELS = {
    "x_velocity": r"$u$ [m/s]",
    "y_velocity": r"$v$ [m/s]",
    "Kmv": r"$K_{mv}$",
}


def load_uniform(pltfile, level):
    """Return dict of 3D numpy arrays (nx, ny, nz) on a uniform grid at `level`."""
    ds = yt.load(pltfile)
    level = min(level, ds.max_level)
    dims = ds.domain_dimensions * ds.refine_by**level
    cg = ds.covering_grid(level=level, left_edge=ds.domain_left_edge, dims=dims)
    data = {}
    for f in FIELDS + ["z_phys"]:
        data[f] = np.asarray(cg[("boxlib", f)].d)
    return ds, data, level


def basename(pltfile):
    return os.path.basename(pltfile.rstrip("/"))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("pltfiles", nargs="+",
                   help="one or more ERF plotfile directories, e.g. plt00100 plt00500")
    p.add_argument("--mode", choices=["mean", "column", "scatter"], default="mean")
    p.add_argument("--i", type=int, default=None, help="i index for column mode (default: center)")
    p.add_argument("--j", type=int, default=None, help="j index for column mode (default: center)")
    p.add_argument("--level", type=int, default=0, help="AMR level to sample (default: 0)")
    p.add_argument("--zmax", type=float, default=None, help="upper z limit for the plot [m]")
    p.add_argument("-o", "--output", default=None, help="output image (default: <pltfile>_profiles.png)")
    p.add_argument("--show", action="store_true", help="also open an interactive window")
    args = p.parse_args()

    fig, axes = plt.subplots(1, 3, figsize=(12, 5), sharey=True)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    zmin = np.inf
    title = None
    levels = []
    last_time = None

    for n, pltfile in enumerate(args.pltfiles):
        ds, d, lev = load_uniform(pltfile, args.level)
        levels.append(lev)
        last_time = float(ds.current_time)
        nx, ny, nz = d["z_phys"].shape
        z = d["z_phys"]
        zmin = min(zmin, float(z.min()))
        color = colors[n % len(colors)]
        label = f"{basename(pltfile)}  (t = {last_time:.2f} s)"

        if args.mode == "mean":
            zprof = z.mean(axis=(0, 1))
            for ax, f in zip(axes, FIELDS):
                ax.plot(d[f].mean(axis=(0, 1)), zprof, "-o", ms=3, color=color, label=label)
            title = "horizontal mean"
        elif args.mode == "column":
            i = nx // 2 if args.i is None else args.i
            j = ny // 2 if args.j is None else args.j
            if not (0 <= i < nx and 0 <= j < ny):
                raise SystemExit(f"(i, j) = ({i}, {j}) out of range for grid {nx} x {ny}")
            for ax, f in zip(axes, FIELDS):
                ax.plot(d[f][i, j, :], z[i, j, :], "-o", ms=3, color=color, label=label)
            title = f"column (i={i}, j={j})"
        else:  # scatter
            zf = z.ravel()
            for ax, f in zip(axes, FIELDS):
                ax.scatter(d[f].ravel(), zf, s=1, alpha=0.3, color=color, label=label)
            title = "all cells"

    for ax, f in zip(axes, FIELDS):
        ax.set_xlabel(LABELS[f])
        ax.grid(True, alpha=0.3)
    axes[0].set_ylabel(r"$z_{phys}$ [m]")
    if args.zmax is not None:
        axes[0].set_ylim(top=args.zmax)
    axes[0].set_ylim(bottom=max(0.0, zmin) if args.mode != "scatter" else zmin)

    if len(args.pltfiles) > 1:
        leg = axes[0].legend(fontsize="small",
                             markerscale=4 if args.mode == "scatter" else 1)
        for h in leg.legend_handles:
            h.set_alpha(1.0)

    lev_str = str(levels[0]) if len(set(levels)) == 1 else ",".join(str(l) for l in levels)
    if len(args.pltfiles) == 1:
        hdr = f"{args.pltfiles[0]}  |  t = {last_time:.2f} s  |  "
    else:
        hdr = f"{len(args.pltfiles)} plotfiles  |  "
    fig.suptitle(f"{hdr}level {lev_str}, {title}")
    fig.tight_layout()

    out = args.output or f"{basename(args.pltfiles[0])}_profiles.png"
    fig.savefig(out, dpi=150)
    print(f"Saved {out}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
