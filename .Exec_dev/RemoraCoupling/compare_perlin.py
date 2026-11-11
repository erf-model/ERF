#!/usr/bin/env python3
"""
Overlay ERF profiles against Perlin et al. (2007, JPO 37, 2081) Figs. 7 and 8,
uncoupled case, 25 km offshore, 72-h forecast.

The reference curves below are hand-digitised from the published figures and
are good to roughly +/-10% in the vertical viscosity and +/-0.3 K in theta.
They are a target, not data.

Usage:
    python compare_perlin.py plt86400 [plt86400_chen ...] [--labels a,b] [-o out.png]

Plots theta, u, v, TKE, K_m and the MYNN master length scale L_m against z.
K_m is taken from "nut" when present (already m^2/s); otherwise from
"Kmv"/density, since ERF writes Kmv as rho*K.
"""
import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yt

yt.set_log_level("error")

# --- Perlin et al. (2007), uncoupled, 25 km offshore, t = 72 h --------------
# Hand-digitised from Figs. 7a-c and 8a-b.
PERLIN = {
    # Pixel-traced from a 400 dpi render of Fig. 7 and Fig. 8, 25 km offshore,
    # uncoupled (grey) curve.  Fig. 8 has three separate 0-6 (8a) and 0-0.6 (8b)
    # panels; calibration is 46.7 px per m^2/s from 12 evenly spaced tick labels.
    # theta/u/v are still read by eye from Fig. 7 and are the least precise here.
    "theta": (
        [286.3, 286.3, 286.3, 286.3, 286.3, 286.4, 287.0, 288.0, 288.6],
        [0, 100, 200, 300, 400, 470, 530, 620, 700]),
    "u": (
        [2.3, 3.6, 4.1, 4.0, 3.9, 3.8, 3.6, 1.2, 0.2],
        [0, 100, 200, 300, 400, 470, 520, 580, 650]),
    "v": (
        [11.0, 13.2, 14.4, 15.0, 15.4, 15.9, 15.3, 15.0, 15.0],
        [0, 100, 200, 300, 400, 470, 540, 620, 700]),
    "tke": (   # Fig. 8b, traced
        [0.409, 0.369, 0.334, 0.306, 0.257, 0.216, 0.152, 0.118, 0.077, 0.015],
        [25, 50, 75, 100, 150, 200, 300, 350, 400, 450]),
    "Km": (    # Fig. 8a, traced
        [4.37, 4.45, 4.71, 4.58, 4.35, 3.60, 1.80, 0.86, 0.64, 0.41, 0.11],
        [50, 75, 100, 125, 150, 200, 300, 350, 400, 450, 500]),
}

# panel key, axis label, x-range (None = autoscale)
PANELS = [
    ("theta", r"$\theta$ [K]", (285.5, 294)),
    ("u", r"$u$ [m/s]", (-3, 5)),
    ("v", r"$-v$ (southward) [m/s]", (8, 17)),
    ("tke", r"TKE [m$^2$/s$^2$]", (0, 0.8)),
    ("Km", r"$K_m$ [m$^2$/s]", (0, 25)),
    ("Lm", r"$L_m$ [m]", (0, 120)),
]


_FLATAU = [6.11239921, 0.443987641, 0.142986287e-1, 0.264847430e-3,
           0.302950461e-5, 0.206739458e-7, 0.640689451e-10,
           -0.952447341e-13, -0.976195544e-15]


def _q_surf_over_sea(theta_surf, p_pa):
    """q_sat over water, as SurfaceLayer::fill_qsurf_with_qsat computes it."""
    T = theta_surf * (p_pa / 1.0e5) ** (287.0 / 1004.5)
    dtt, pm = T - 273.16, p_pa * 0.01
    e = _FLATAU[8]
    for c in reversed(_FLATAU[:8]):
        e = c + dtt * e
    return (287.0 / 461.505) * e / max(e, pm - e)


def load(pltfile, i=0, j=0):
    ds = yt.load(pltfile)
    cg = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    have = {str(n[1]) for n in ds.field_list}

    def f(name):
        return np.asarray(cg[("boxlib", name)].d)[i, j, :].astype(float)

    d = {"z": f("z_phys"), "theta": f("theta"), "u": f("x_velocity"),
         "tke": f("KE"), "qv": f("qv") if "qv" in have else None,
         "rho": f("density"), "p": f("pressure")}
    d["v"] = -f("y_velocity")                      # plot southward as positive
    d["Km"] = f("nut") if "nut" in have else f("Kmv") / f("density")
    d["Lm"] = f("Lturb") if "Lturb" in have else None
    d["t"] = float(ds.current_time)
    return d


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("pltfiles", nargs="+")
    p.add_argument("--labels", default=None,
                   help="comma-separated legend labels, one per plotfile")
    p.add_argument("--i", type=int, default=0)
    p.add_argument("--j", type=int, default=0)
    p.add_argument("--zmax", type=float, default=800.0)
    p.add_argument("-o", "--output", default="perlin_compare.png")
    args = p.parse_args()

    labels = (args.labels.split(",") if args.labels
              else [os.path.basename(f.rstrip("/")) for f in args.pltfiles])
    if len(labels) != len(args.pltfiles):
        raise SystemExit("--labels must give one label per plotfile")

    fig, axes = plt.subplots(1, len(PANELS), figsize=(3.0 * len(PANELS), 5.2),
                             sharey=True)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for n, (pf, lab) in enumerate(zip(args.pltfiles, labels)):
        d = load(pf, args.i, args.j)
        color = colors[n % len(colors)]
        tag = f"{lab}  (t = {d['t']/3600:.1f} h)"
        for ax, (key, _, _xl) in zip(axes, PANELS):
            if d.get(key) is None:
                ax.text(0.5, 0.5, "not in plotfile", ha="center",
                        transform=ax.transAxes, fontsize="small", color="0.5")
                continue
            ax.plot(d[key], d["z"], "-o", ms=3, color=color,
                    label=tag if key == "theta" else None)

    for ax, (key, xlabel, xlim) in zip(axes, PANELS):
        if key in PERLIN:
            x, z = PERLIN[key]
            ax.plot(x, z, "k--s", ms=4, lw=1.5, alpha=0.8,
                    label="Perlin 2007 (digitised)" if key == "theta" else None)
        ax.set_xlabel(xlabel)
        if xlim is not None:
            ax.set_xlim(*xlim)
        ax.grid(True, alpha=0.3)

    axes[0].set_ylabel(r"$z$ [m]")
    axes[0].set_ylim(0, args.zmax)
    axes[0].legend(fontsize="x-small", loc="upper left")
    fig.suptitle("ERF vs Perlin et al. (2007), uncoupled case, 25 km offshore")
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f"Saved {args.output}")

    # Scalar summary against the paper's headline numbers.
    print()
    hdr = (f"{'case':<14}{'peakKm':>8}{'z@pk':>6}{'maxLm':>7}{'MBL':>6}"
           f"{'TKE0':>7}{'th_ml':>8}{'qv0':>8}{'H':>7}{'LE':>7}{'L_ob':>8}")
    print(hdr)
    print("-" * len(hdr))
    print(f"{'Perlin 2007':<14}{4.5:>8.2f}{100:>6.0f}{12:>7.0f}{475:>6.0f}"
          f"{0.42:>7.3f}{286.3:>8.2f}{'':>8}{-10:>7.0f}{-76:>7.0f}{'<0':>8}")
    for pf, lab in zip(args.pltfiles, labels):
        d = load(pf, args.i, args.j)
        k = int(np.argmax(d["Km"]))
        tke = d["tke"]
        above = np.where(tke > 0.01 * max(tke.max(), 1e-12))[0]
        mbl = d["z"][above.max()] if len(above) else 0.0
        lm = d["Lm"].max() if d["Lm"] is not None else float("nan")
        H = LE = Lob = float("nan")
        if d["qv"] is not None:
            # bulk_coeff surface fluxes, with q_surf = q_sat over sea
            wsp = np.hypot(d["u"][0], -d["v"][0])
            qs = _q_surf_over_sea(287.15, d["p"][0])
            wth = 0.0012 * wsp * (287.15 - d["theta"][0])
            wq = 0.0012 * wsp * (qs - d["qv"][0])
            H = d["rho"][0] * 1004.5 * wth
            LE = d["rho"][0] * 2.5e6 * wq
            th0 = d["theta"][0] * (1 + 0.608 * d["qv"][0])
            shf = wth * (1 + 0.608 * d["qv"][0]) + 0.608 * th0 * wq
            ust = np.sqrt(0.0014) * wsp
            Lob = -(th0 * ust ** 3) / (0.41 * 9.81 * shf) if shf != 0 else np.inf
        print(f"{lab:<14}{d['Km'].max():>8.2f}{d['z'][k]:>6.0f}{lm:>7.1f}"
              f"{mbl:>6.0f}{tke[0]:>7.3f}{d['theta'][0]:>8.2f}"
              f"{d['qv'][0]*1000 if d['qv'] is not None else 0:>8.3f}"
              f"{H:>7.1f}{LE:>7.1f}{Lob:>8.0f}")


if __name__ == "__main__":
    main()
