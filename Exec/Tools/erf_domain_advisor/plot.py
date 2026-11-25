"""--plot: domain.png, the nested boxes over the terrain, and slopes.png, the slope on each level.

Needs matplotlib; the rest of the tool does not.
"""

import math


def _boxes(ax, plan, upto=None, scale=1e-3):
    from matplotlib.patches import Rectangle
    colors = ["black", "tab:blue", "tab:orange", "tab:red", "tab:purple", "tab:brown", "tab:pink"]
    for lv in plan.levels[1:upto]:
        ax.add_patch(Rectangle((lv.x0 * scale, lv.y0 * scale), (lv.x1 - lv.x0) * scale,
                               (lv.y1 - lv.y0) * scale, fill=False, lw=1.6,
                               ec=colors[lv.lev % len(colors)],
                               label="level %d, dx %.4g m" % (lv.lev, lv.dx)))


def _point(ax, plan, scale=1e-3):
    ax.plot(plan.point_xy[0] * scale, plan.point_xy[1] * scale, marker="*", ms=13, mfc="yellow",
            mec="black", ls="none", label="point of interest")


def _wind(ax, plan, lx, ly, scale=1e-3):
    """An arrow blowing toward the domain centre from the side the wind comes from (north up;
    the grid's north differs by the Lambert rotation, a few degrees)."""
    wf = plan.settings.wind_from
    if wf is None:
        return
    th = math.radians(wf)
    ux, uy = math.sin(th), math.cos(th)        # toward where the wind comes from
    r, length = 0.40 * min(lx, ly), 0.16 * min(lx, ly)
    x0, y0 = 0.5 * lx + r * ux, 0.5 * ly + r * uy
    ax.annotate("", xy=((x0 - length * ux) * scale, (y0 - length * uy) * scale),
                xytext=(x0 * scale, y0 * scale),
                arrowprops=dict(arrowstyle="-|>", lw=2.5, color="black"))
    ax.text(x0 * scale, y0 * scale, "wind from %g deg " % wf, fontsize=9,
            ha="right" if ux > 0 else "left", va="bottom" if uy > 0 else "top")


def _credit(fig, dem_name):
    if dem_name and "Terrain Tiles" in dem_name:
        fig.text(0.005, 0.005, "Terrain: AWS Terrain Tiles (Mapzen, Tilezen); 3DEP, GMTED2010, SRTM "
                 "courtesy of the U.S. Geological Survey, and other providers listed in "
                 "erf_domain_advisor/README.md", fontsize=6.5, color="0.3")
    elif dem_name:
        fig.text(0.005, 0.005, "Terrain: %s" % dem_name, fontsize=6.5, color="0.3")


def write(plan, terrain, caution, limit, outdir, dem_name=None):
    """terrain: {lev: (h, water, dx)} from terrain.check (dx is the sample spacing, larger than the
    cell size when a big level was strided), or {} without elevation data."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import Rectangle
    from .terrain import cell_slope

    f = plan.frame
    rw = plan.settings.real_width * plan.base_dx
    land = matplotlib.colors.ListedColormap(plt.cm.terrain(np.linspace(0.25, 1.0, 256)))
    names = []

    # ---- domain.png: level 0 with every box, and a zoom on the refined boxes ----
    fig, axes = plt.subplots(1, 2, figsize=(13, 6), constrained_layout=True)
    for k, ax in enumerate(axes):
        lv0 = plan.levels[0] if k == 0 else (plan.levels[1] if len(plan.levels) > 1 else plan.levels[0])
        show = 0 if k == 0 else lv0.lev
        if show in terrain:
            h, water, _ = terrain[show]
            lv = plan.levels[show]
            ext = [lv.x0 / 1e3, lv.x1 / 1e3, lv.y0 / 1e3, lv.y1 / 1e3]
            hm = np.ma.masked_where(water > 0.5, h)
            ax.imshow(np.zeros_like(h), extent=ext, origin="lower", cmap="Blues", vmin=-1, vmax=1)
            im = ax.imshow(hm, extent=ext, origin="lower", cmap=land,
                           vmin=0, vmax=max(1.0, float(np.nanmax(h))))
            fig.colorbar(im, ax=ax, shrink=0.8, label="terrain height [m] (water in blue)")
        if k == 0:
            ax.add_patch(Rectangle((0, 0), f.lx / 1e3, f.ly / 1e3, fill=False, lw=2, ec="black",
                                   label="level 0, dx %.4g m" % plan.base_dx))
            ax.add_patch(Rectangle((rw / 1e3, rw / 1e3), (f.lx - 2 * rw) / 1e3, (f.ly - 2 * rw) / 1e3,
                                   fill=False, lw=1, ls="--", ec="black",
                                   label="inner edge of the relaxation zone"))
            ax.set_xlim(0, f.lx / 1e3)
            ax.set_ylim(0, f.ly / 1e3)
            ax.set_title("Level 0 (%.0f x %.0f km) and the refined boxes" % (f.lx / 1e3, f.ly / 1e3))
            _wind(ax, plan, f.lx, f.ly)
        else:
            ax.set_xlim(lv0.x0 / 1e3, lv0.x1 / 1e3)
            ax.set_ylim(lv0.y0 / 1e3, lv0.y1 / 1e3)
            ax.set_title("Zoom on level %d" % lv0.lev)
        _boxes(ax, plan)
        _point(ax, plan)
        ax.set_xlabel("x [km]")
        ax.set_ylabel("y [km]")
        ax.set_aspect("equal")
        ax.legend(loc="lower right", fontsize=8, framealpha=0.85)
    fig.suptitle("%s at %.4f, %.4f: %d refined level(s), finest %.4g m%s" % (
        plan.settings.phenomenon, plan.settings.lat, (plan.settings.lon + 180) % 360 - 180,
        len(plan.levels) - 1, plan.target_dx,
        ", fire grid %.3g m" % (plan.target_dx / plan.fire_ratio) if plan.fire else ""))
    _credit(fig, dem_name)
    names.append("domain.png")
    fig.savefig("%s/domain.png" % outdir, dpi=110)
    plt.close(fig)

    # ---- slopes.png: slope on each level, coloured against the caution and limit values ----
    if terrain:
        n = len(plan.levels)
        fig, axes = plt.subplots(1, n, figsize=(4.2 * n + 1, 4.4), constrained_layout=True,
                                 squeeze=False)
        cmap = matplotlib.colors.ListedColormap(["#e8f1e4", "#a6d49a", "#f2c45a", "#d9442b"])
        norm = matplotlib.colors.BoundaryNorm([0, 0.5 * caution, caution, limit, 1e9], cmap.N)
        for lv, ax in zip(plan.levels, axes[0]):
            h, _, dxs = terrain[lv.lev]
            s = cell_slope(h, dxs)
            ext = [lv.x0 / 1e3, lv.x1 / 1e3, lv.y0 / 1e3, lv.y1 / 1e3]
            im = ax.imshow(s, extent=ext, origin="lower", cmap=cmap, norm=norm,
                           interpolation="nearest")
            _boxes(ax, plan)
            ax.set_xlim(ext[0], ext[1])
            ax.set_ylim(ext[2], ext[3])
            ax.set_title("level %d, dx %.4g m: max %.2f" % (lv.lev, lv.dx, float(np.nanmax(s))),
                         fontsize=10)
            ax.set_xlabel("x [km]")
            ax.set_aspect("equal")
        axes[0][0].set_ylabel("y [km]")
        cb = fig.colorbar(im, ax=list(axes[0]), shrink=0.85, ticks=[0.25 * caution, 0.75 * caution,
                                                                  0.5 * (caution + limit), limit + 0.2])
        cb.ax.set_yticklabels(["< %.2g" % (0.5 * caution), "< %.2g" % caution,
                               "%.2g-%.2g (note)" % (caution, limit), "> %.2g (warning)" % limit])
        cb.set_label("terrain slope |grad h|")
        fig.suptitle("Terrain slope on each level, as the terrain-following grid sees it")
        _credit(fig, dem_name)
        names.append("slopes.png")
        fig.savefig("%s/slopes.png" % outdir, dpi=110)
        plt.close(fig)
    return names
