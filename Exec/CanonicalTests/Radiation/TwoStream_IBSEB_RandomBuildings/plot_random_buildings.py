#!/usr/bin/env python3
"""Plot the random-building morning (numpy and matplotlib).

    python3 plot_random_buildings.py [run directory]

Writes, in the run directory:
  random_buildings_map.png     the last 2D plotfile seen from above: the
                               force-restore skin of the open ground, level 1
                               over the refined area and level 0 around it
                               (grey under and around the buildings), the building
                               footprints, the level-1 box, and the roofs
                               coloured by their skin temperature;
  random_buildings_tskin.png   the mean skin temperature of every building on
                               both levels against time, the open ground's
                               force-restore skin over the refined area (from
                               the 2D plotfiles, footprints left out), and the
                               sun's zenith;
  random_buildings_theta.png   potential temperature 15 m above the ground on
                               level 1 at the last 3D plotfile, buildings
                               masked, with the wind.
"""
import glob
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_random_buildings import read_csv, read_dump, dump_steps   # noqa: E402


def read_plotfile(path):
    """Every level of an AMReX plotfile: (names, [[(lo, array[ncomp, ni, nj, nk]) per box] per level], dx)."""
    with open(os.path.join(path, "Header")) as f:
        lines = f.read().split("\n")
    nvar = int(lines[1])
    names = lines[2:2 + nvar]
    p = 2 + nvar + 2
    finest = int(lines[p]); p += 1
    p += 5
    dxs = []
    for _ in range(finest + 1):
        dxs.append([float(v) for v in lines[p].split()]); p += 1
    levels = []
    for lev in range(finest + 1):
        boxes = []
        with open(os.path.join(path, f"Level_{lev}", "Cell_H")) as f:
            ch = f.read()
        for name, off in re.findall(r"FabOnDisk: (\S+) (\d+)", ch):
            with open(os.path.join(path, f"Level_{lev}", name), "rb") as f:
                f.seek(int(off))
                head = b""
                while not head.endswith(b"\n"):
                    head += f.read(1)
                m = re.search(rb"\(\((-?\d+),(-?\d+),(-?\d+)\) \((-?\d+),(-?\d+),(-?\d+)\) \(\d+,\d+,\d+\)\) (\d+)", head)
                lo = [int(m.group(k)) for k in (1, 2, 3)]
                hi = [int(m.group(k)) for k in (4, 5, 6)]
                nc = int(m.group(7))
                n = [hi[d] - lo[d] + 1 for d in range(3)]
                order = "<" if b"(8 7 6 5 4 3 2 1)" in head else ">"
                data = np.frombuffer(f.read(8 * n[0] * n[1] * n[2] * nc), dtype=order + "f8")
                boxes.append((lo, data.reshape(nc, n[2], n[1], n[0]).transpose(0, 3, 2, 1)))
        levels.append(boxes)
    return names, levels, dxs


def plotfile_time(path):
    """Simulation time [s] in an AMReX plotfile's Header."""
    with open(os.path.join(path, "Header")) as f:
        lines = f.read().split("\n")
    return float(lines[2 + int(lines[1]) + 1])


def assemble(boxes, comp, shape, k=None):
    """One component on the level's index space (a k slice, or k = 0 for a 2D file), NaN outside its boxes."""
    f = np.full(shape, np.nan)
    for lo, arr in boxes:
        kk = 0 if k is None else k - lo[2]
        if k is not None and not 0 <= kk < arr.shape[3]:
            continue
        f[lo[0]:lo[0] + arr.shape[1], lo[1]:lo[1] + arr.shape[2]] = arr[comp, :, :, kk]
    return f


def numeric(s):
    return int(re.sub(r"\D", "", s))


def main():
    run = sys.argv[1] if len(sys.argv) > 1 else "."
    os.chdir(run)
    bld = read_csv("buildings.csv")
    nb = len(bld["building"])
    c = read_csv("ibseb_buildings.csv")
    L = 640.0

    def keep_out(shape, dxl, g=40.0):
        """Columns within g of a footprint: the buildings, their ramp cells and the cell around them."""
        x = (np.arange(shape[0]) + 0.5) * dxl[0]; y = (np.arange(shape[1]) + 0.5) * dxl[1]
        ko = np.zeros(shape, dtype=bool)
        for b in range(nb):
            ko |= np.outer((x > bld["x_lo_m"][b] - g) & (x < bld["x_hi_m"][b] + g),
                           (y > bld["y_lo_m"][b] - g) & (y < bld["y_hi_m"][b] + g))
        return ko

    # ---- map --------------------------------------------------------------
    p2 = sorted(glob.glob("plt2d*"), key=numeric)[-1]
    names, levels, dxs = read_plotfile(p2)
    k = names.index("seb_t_sfc")
    n0 = (int(round(L / dxs[0][0])), int(round(L / dxs[0][1])))
    f0 = assemble(levels[0], k, n0)
    f1 = assemble(levels[1], k, (2 * n0[0], 2 * n0[1]))
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    # The two-stream ground balance also covers the footprints, which run much
    # hotter; blank them and their surroundings and scale to the open ground.
    f0[keep_out(f0.shape, dxs[0])] = np.nan
    f1[keep_out(f1.shape, dxs[1])] = np.nan
    vmin, vmax = np.nanmin(np.concatenate([f0.ravel(), f1.ravel()])), np.nanmax(np.concatenate([f0.ravel(), f1.ravel()]))
    ax.set_facecolor("0.85")
    ax.imshow(np.ma.masked_invalid(f0).T, origin="lower", extent=(0, L, 0, L), cmap="inferno", vmin=vmin, vmax=vmax)
    im = ax.imshow(np.ma.masked_invalid(f1).T, origin="lower", extent=(0, L, 0, L), cmap="inferno", vmin=vmin, vmax=vmax)
    fig.colorbar(im, ax=ax, label="force-restore skin, open ground [K]", shrink=0.8)
    lo1 = min(lo[0] for lo, _ in levels[1]) * dxs[1][0]; hi1 = max(lo[0] + a.shape[1] for lo, a in levels[1]) * dxs[1][0]
    lo1y = min(lo[1] for lo, _ in levels[1]) * dxs[1][1]; hi1y = max(lo[1] + a.shape[2] for lo, a in levels[1]) * dxs[1][1]
    ax.add_patch(Rectangle((lo1, lo1y), hi1 - lo1, hi1y - lo1y, fill=False, ec="cyan", lw=1.5, ls="--", label="level 1"))
    s1 = dump_steps("faces/set.lev1")
    d1 = read_dump("faces/set.lev1", s1[-1])
    roof = d1["dir"] == 2
    sc = ax.scatter(d1["x_m"][roof], d1["y_m"][roof], c=d1["T_skin"][roof], s=18, marker="s", cmap="viridis")
    fig.colorbar(sc, ax=ax, label="roof skin, level 1 [K]", shrink=0.8)
    for b in range(nb):
        ax.add_patch(Rectangle((bld["x_lo_m"][b], bld["y_lo_m"][b]), bld["x_hi_m"][b] - bld["x_lo_m"][b],
                               bld["y_hi_m"][b] - bld["y_lo_m"][b], fill=False, ec="white", lw=0.8))
        ax.text(bld["x_hi_m"][b] + 3, bld["y_hi_m"][b] + 3, f"{int(bld['building'][b])}: {bld['height_m'][b]:g} m",
                color="white", fontsize=7)
    t_end = c["time_s"].max()
    ax.set_title(f"{p2}: {t_end / 3600:.1f} h after 15:00 UTC")
    ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]"); ax.legend(loc="lower left", fontsize=8)
    fig.tight_layout(); fig.savefig("random_buildings_map.png", dpi=130); plt.close(fig)

    # ---- skin temperatures -------------------------------------------------
    fig, axs = plt.subplots(2, 1, figsize=(8, 7), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    cmap = plt.get_cmap("tab10")
    lev = c["level"].astype(int)
    for b in range(1, nb + 1):
        for L_, ls in ((0, "--"), (1, "-")):
            m = (lev == L_) & (c["building"] == b)
            axs[0].plot(c["time_s"][m] / 3600.0, c["T_skin_mean_K"][m], ls, color=cmap((b - 1) % 10),
                        label=f"{b} ({bld['height_m'][b - 1]:g} m)" if L_ == 1 else None)
    # The open ground's force-restore skin over level 1's area, away from the
    # footprints and the cell around them, from the 2D plotfiles.
    tg, g0, g1 = [], [], []
    for p2d in sorted(glob.glob("plt2d*"), key=numeric):
        nm, lv, dx2 = read_plotfile(p2d)
        if len(lv) < 2:
            continue
        kk = nm.index("seb_t_sfc")
        n0_ = (int(round(L / dx2[0][0])), int(round(L / dx2[0][1])))
        a0 = assemble(lv[0], kk, n0_)
        a1 = assemble(lv[1], kk, (2 * n0_[0], 2 * n0_[1]))
        on1 = ~np.isnan(a1)
        on0 = np.zeros(n0_, dtype=bool)
        for lo, arr in lv[1]:
            on0[lo[0] // 2:(lo[0] + arr.shape[1]) // 2, lo[1] // 2:(lo[1] + arr.shape[2]) // 2] = True
        means = [f[on & ~keep_out(f.shape, dxl)].mean() for f, on, dxl in ((a0, on0, dx2[0]), (a1, on1, dx2[1]))]
        tg.append(plotfile_time(p2d)); g0.append(means[0]); g1.append(means[1])
    tg = np.array(tg)
    axs[0].plot(tg / 3600.0, g0, "--", color="k", lw=2)
    axs[0].plot(tg / 3600.0, g1, "-", color="k", lw=2, label="open ground")
    axs[0].set_ylabel("mean skin temperature [K]")
    axs[0].legend(fontsize=7, ncol=3, title="building (height); solid level 1, dashed level 0", title_fontsize=7)
    first = (lev == 0) & (c["building"] == 1)
    axs[1].plot(c["time_s"][first] / 3600.0, c["sun_zenith_deg"][first], "C3")
    axs[1].set_ylabel("sun zenith [deg]"); axs[1].set_xlabel("hours after 15:00 UTC (08:20 local solar)")
    fig.tight_layout(); fig.savefig("random_buildings_tskin.png", dpi=130); plt.close(fig)

    # ---- air above the ground ---------------------------------------------
    p3 = sorted([p for p in glob.glob("plt*") if re.fullmatch(r"plt\d+", p)], key=numeric)[-1]
    names, levels, dxs = read_plotfile(p3)
    kz = int(15.0 / dxs[1][2])
    nf = (2 * n0[0], 2 * n0[1])
    th = assemble(levels[1], names.index("theta"), nf, kz)
    mask = assemble(levels[1], names.index("terrain_IB_mask"), nf, kz) if "terrain_IB_mask" in names else None
    u = assemble(levels[1], names.index("x_velocity"), nf, kz)
    v = assemble(levels[1], names.index("y_velocity"), nf, kz)
    if mask is not None:
        solid = mask > 0.5
        th[solid] = np.nan; u[solid] = np.nan; v[solid] = np.nan
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(np.ma.masked_invalid(th).T, origin="lower", extent=(0, L, 0, L), cmap="RdBu_r")
    fig.colorbar(im, ax=ax, label="theta [K]", shrink=0.8)
    x = (np.arange(nf[0]) + 0.5) * dxs[1][0]; y = (np.arange(nf[1]) + 0.5) * dxs[1][1]
    X, Y = np.meshgrid(x, y, indexing="ij")
    sl = (slice(None, None, 2), slice(None, None, 2))
    ok = ~np.isnan(u[sl])
    ax.quiver(X[sl][ok], Y[sl][ok], u[sl][ok], v[sl][ok], scale=80, width=0.002)
    for b in range(nb):
        ax.add_patch(Rectangle((bld["x_lo_m"][b], bld["y_lo_m"][b]), bld["x_hi_m"][b] - bld["x_lo_m"][b],
                               bld["y_hi_m"][b] - bld["y_lo_m"][b], fill=False, ec="k", lw=0.8))
    ax.set_xlim(lo1, hi1); ax.set_ylim(lo1y, hi1y)
    ax.set_title(f"{p3}: level 1, z = {(kz + 0.5) * dxs[1][2]:g} m")
    ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
    fig.tight_layout(); fig.savefig("random_buildings_theta.png", dpi=130); plt.close(fig)
    print("wrote random_buildings_map.png, random_buildings_tskin.png, random_buildings_theta.png")


if __name__ == "__main__":
    main()
