#!/usr/bin/env python3
"""Plot the building-set morning.

    python3 plot_buildingset.py ibseb_set.csv faces/set [--plotfile plt43200] [--out plots] [--dt 0.5]

  tskin_buildings.png   mean skin temperature of each building against solar time
  shadow_buildings.png  shadow fraction of each building (mutual shading)
  faces_map.png         the faces at the last dump coloured by skin temperature,
                        seen from above (walls drawn at their positions)
  tskin_slice.png       yt slice of ibseb_tskin at 25 m from the plotfile (needs yt)
  shading_maps.png      every face coloured by the incident direct beam at three solar
                        times, seen from the south-east, with the sun position in the
                        titles (zero = cast shadow or facing away from the sun)
  shading_top.png       the same from above: roofs as squares, wall columns as the
                        column-mean incident direct beam
  wind_near_surface.png horizontal wind speed and vectors in the first two cells above
                        the ground (5 m and 15 m) at the plotfile time, buildings masked,
                        and a vertical slice of u through the slab and the cube (needs yt)
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

def sun_at(c, t_s):
    """Sun zenith and azimuth from the summary CSV at the row nearest a time [s]."""
    n = np.argmin(np.abs(c["time_s"] - t_s))
    return c["sun_zenith_deg"][n], c["sun_azimuth_deg"][n]

def shading_maps(files, c, dt, out):
    """Absorbed shortwave on every face at three solar times, 3-D view and top-down view."""
    picks = []
    for hour in (6.0, 8.0, 10.0):
        want = int(round((hour - 5.0) * 3600.0 / dt))
        fn = min(files, key=lambda f: abs(int(re.search(r"step(\d+)", f).group(1)) - want))
        picks.append((hour, fn))
    vmax = 0.0
    data = []
    for hour, fn in picks:
        d = load(fn); data.append(d); vmax = max(vmax, d["SW_direct_in"].max())
    fig = plt.figure(figsize=(16, 5.5))
    for n, ((hour, fn), d) in enumerate(zip(picks, data)):
        ax = fig.add_subplot(1, 3, n + 1, projection="3d")
        sc = ax.scatter(d["x_m"], d["y_m"], d["z_m"], c=d["SW_direct_in"], cmap="inferno", vmin=0, vmax=vmax, s=18, marker="s", depthshade=False)
        ax.view_init(elev=35, azim=-50)
        ax.set_xlim(100, 400); ax.set_ylim(100, 400); ax.set_zlim(0, 70)
        ax.set_box_aspect((1, 1, 0.25))
        zen, az = sun_at(c, (hour - 5.0) * 3600.0)
        sh = d["shadow"].mean(); dark = (d["SW_direct_in"] <= 0.0).mean()
        ax.set_title(f"{hour:04.1f} h solar time, sun zenith {zen:.0f}°, azimuth {az:.0f}°\n{100 * dark:.0f} % of faces without direct sun\n({100 * sh:.0f} % by cast shadow, the rest facing away)", fontsize=8.5)
        ax.set_xlabel("x [m]", fontsize=8); ax.set_ylabel("y [m]", fontsize=8); ax.set_zlabel("z [m]", fontsize=8)
    fig.subplots_adjust(left=0.02, right=0.88, wspace=0.05)
    fig.colorbar(sc, ax=fig.axes, shrink=0.6, pad=0.02, label="incident direct beam [W/m2]")
    fig.savefig(os.path.join(out, "shading_maps.png"), dpi=130, bbox_inches="tight"); plt.close(fig)
    fig, axs = plt.subplots(1, 3, figsize=(16, 5.2))
    for ax, (hour, fn), d in zip(axs, picks, data):
        roof = d["dir"] == 2
        sc = ax.scatter(d["x_m"][roof], d["y_m"][roof], c=d["SW_direct_in"][roof], cmap="inferno", vmin=0, vmax=vmax, s=60, marker="s", edgecolors="none")
        # Wall columns: all faces of a column share (x, y); colour by the column-mean direct beam.
        cols = {}
        for x, y, dr, sw in zip(d["x_m"][~roof], d["y_m"][~roof], d["dir"][~roof], d["SW_direct_in"][~roof]):
            cols.setdefault((x, y, int(dr)), []).append(sw)
        for (x, y, dr), v in cols.items():
            col = plt.cm.inferno(min(np.mean(v) / vmax, 1.0))
            if dr == 0: ax.plot([x, x], [y - 5, y + 5], color=col, lw=5, solid_capstyle="butt")
            else:       ax.plot([x - 5, x + 5], [y, y], color=col, lw=5, solid_capstyle="butt")
        zen, az = sun_at(c, (hour - 5.0) * 3600.0)
        ax.set_aspect("equal"); ax.set_xlim(100, 400); ax.set_ylim(100, 400)
        ax.set_title(f"{hour:04.1f} h: roofs (squares) and wall columns (bars) by incident direct beam\nsun zenith {zen:.0f}°, azimuth {az:.0f}°", fontsize=9)
        ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
    fig.subplots_adjust(left=0.05, right=0.9, wspace=0.25)
    fig.colorbar(sc, ax=axs, shrink=0.8, pad=0.02, label="incident direct beam [W/m2]"); fig.savefig(os.path.join(out, "shading_top.png"), dpi=130); plt.close(fig)
    print(f"wrote {out}/shading_maps.png shading_top.png")

def wind_near_surface(plotfile, out):
    """Horizontal wind in the first two cells above the ground and a vertical slice of u."""
    import yt
    ds = yt.load(plotfile)
    g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    u = g["boxlib", "x_velocity"].value; v = g["boxlib", "y_velocity"].value; w = g["boxlib", "z_velocity"].value
    m = g["boxlib", "terrain_IB_mask"].value
    dx = (ds.domain_right_edge - ds.domain_left_edge).value / ds.domain_dimensions
    nx, ny, nz = u.shape
    xc = (np.arange(nx) + 0.5) * dx[0]; yc = (np.arange(ny) + 0.5) * dx[1]; zc = (np.arange(nz) + 0.5) * dx[2]
    t_h = 5.0 + float(ds.current_time) / 3600.0
    solid = m > 0.5
    spd = np.sqrt(u * u + v * v)
    vmax = np.nanmax(np.where(solid, np.nan, spd)[:, :, :2])
    fig, axs = plt.subplots(1, 3, figsize=(16, 5))
    for ax, k in zip(axs[:2], (0, 1)):
        S = np.ma.masked_where(solid[:, :, k], spd[:, :, k]).T
        pc = ax.pcolormesh(xc - dx[0] / 2, yc - dx[1] / 2, S, cmap="viridis", vmin=0, vmax=vmax, shading="auto")
        ax.pcolormesh(xc - dx[0] / 2, yc - dx[1] / 2, np.ma.masked_where(~solid[:, :, k], solid[:, :, k]).T, cmap="Greys", vmin=0, vmax=1.4, shading="auto")
        uu = np.ma.masked_where(solid[:, :, k], u[:, :, k]).T; vv = np.ma.masked_where(solid[:, :, k], v[:, :, k]).T
        ax.quiver(xc, yc, uu, vv, color="w", scale=60, width=0.0035)
        ax.set_aspect("equal"); ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
        ax.set_title(f"horizontal wind at z = {zc[k]:.0f} m (cell {k + 1} above the ground), {t_h:.1f} h solar time", fontsize=9)
        fig.colorbar(pc, ax=ax, shrink=0.8, label="|U_h| [m/s]")
    # Vertical slice of u through the middle of the slab and the cube (y = 230 m).
    j = int(np.argmin(np.abs(yc - 230.0)))
    ax = axs[2]
    U = np.ma.masked_where(solid[:, j, :], u[:, j, :]).T
    pc = ax.pcolormesh(xc - dx[0] / 2, zc - dx[2] / 2, U, cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto")
    ax.pcolormesh(xc - dx[0] / 2, zc - dx[2] / 2, np.ma.masked_where(~solid[:, j, :], solid[:, j, :]).T, cmap="Greys", vmin=0, vmax=1.4, shading="auto")
    ww = np.ma.masked_where(solid[:, j, :], w[:, j, :]).T
    ax.quiver(xc, zc, np.ma.masked_where(solid[:, j, :], u[:, j, :]).T, ww, color="k", scale=60, width=0.003)
    ax.set_xlabel("x [m]"); ax.set_ylabel("z [m]"); ax.set_ylim(0, 100)
    ax.set_title(f"u and (u, w) vectors in the x-z plane at y = {yc[j]:.0f} m (slab and cube)", fontsize=9)
    fig.colorbar(pc, ax=ax, shrink=0.8, label="u [m/s]")
    fig.tight_layout(); fig.savefig(os.path.join(out, "wind_near_surface.png"), dpi=130); plt.close(fig)
    print(f"wrote {out}/wind_near_surface.png")

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
    shading_maps(files, c, a.dt, a.out)
    if a.plotfile:
        wind_near_surface(a.plotfile, a.out)
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
