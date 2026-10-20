#!/usr/bin/env python3
"""
Mixed-layer heat budget and entrainment diagnostic for the Perlin comparison.

The question this answers: ERF's marine boundary layer keeps deepening (36-139 m
per 12 h at hour 72) while Perlin et al. (2007) cap theirs at 475 m.  Which term
in the heat budget is responsible?

The problem is 1-D and has no resolved eddies, so the total turbulent flux IS the
subgrid flux, -K_h d(theta_v)/dz, evaluated on z-faces.  That makes the budget
closable from the plotfiles alone:

    d/dt int_0^zt rho*cp*theta dz  =  H_surface  +  F_top  +  R_radiative

with zt chosen below erf.rayleigh_zdamp so the Rayleigh damping term (which
targets the INITIAL profile and would otherwise act as an unaccounted source)
is identically zero.

Usage:
    python heat_budget.py base D_both E_all
"""
import argparse
import pathlib
import sys

import numpy as np
import yt

yt.set_log_level("error")

CP = 1004.5
LV = 2.5e6
R_D = 287.0
R_V = 461.505
EPSV = R_V / R_D - 1.0
KAPPA = 0.41
GRAV = 9.81

CD, CH, CQ = 0.0014, 0.0012, 0.0012
THETA_SFC = 287.15
Q_RAD = -1.1574074074e-5          # K/s, the 1 K/day cooling
ZTOP = 900.0                      # integrate below erf.rayleigh_zdamp = 1000

_FLATAU = [6.11239921, 0.443987641, 0.142986287e-1, 0.264847430e-3,
           0.302950461e-5, 0.206739458e-7, 0.640689451e-10,
           -0.952447341e-13, -0.976195544e-15]


def q_sat_sea(theta_surf, p_pa):
    """q_surf as SurfaceLayer::fill_qsurf_with_qsat sets it over water."""
    T = theta_surf * (p_pa / 1.0e5) ** (R_D / CP)
    dtt, pm = T - 273.16, p_pa * 0.01
    e = _FLATAU[8]
    for c in reversed(_FLATAU[:8]):
        e = c + dtt * e
    return (R_D / R_V) * e / max(e, pm - e)


def q_rad_of(case):
    """Radiative cooling rate [K/s] this run was given, read from its log.

    The run echoes prob.advection_heating_rate; falling back to zero is right
    for a run that was not given the forcing at all.
    """
    import re
    for name in (f"run_{case}/run.log", f"{case}/run.log"):
        try:
            txt = pathlib.Path(name).read_text()
        except OSError:
            continue
        m = re.search(r"advection_heating_rate\s*=\s*(-?[\d.eE+-]+)", txt)
        if m:
            return float(m.group(1))
    return 0.0


def load(case, plt):
    ds = yt.load(f"run_{case}/{plt}")
    cg = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    have = {str(n[1]) for n in ds.field_list}

    def f(n):
        return np.asarray(cg[("boxlib", n)].d)[0, 0, :].astype(float)

    d = {n: f(n) for n in ("z_phys", "theta", "density", "pressure", "KE")}
    d["qv"] = f("qv") if "qv" in have else np.zeros_like(d["theta"])
    d["Khv"] = f("Khv")                       # stored as rho*K
    d["Km"] = f("nut") if "nut" in have else f("Kmv") / d["density"]
    d["thv"] = d["theta"] * (1.0 + EPSV * d["qv"])
    d["t"] = float(ds.current_time)
    # cell faces reconstructed from the cell-center heights
    nz = len(d["z_phys"])
    zf = np.empty(nz + 1)
    zf[0] = 0.0
    for k in range(nz):
        zf[k + 1] = 2.0 * d["z_phys"][k] - zf[k]
    d["zf"] = zf
    d["dz"] = np.diff(zf)
    return d


def surface_fluxes(d):
    """bulk_coeff surface layer, with q_surf = q_sat over sea."""
    # the plotfile has no face velocities; reconstruct |U| from the first cell
    return None


def flux_profile(d):
    """Turbulent flux of theta_v on z-faces, -K_h d(thv)/dz  [K m/s]."""
    z, thv, rho, Khv = d["z_phys"], d["thv"], d["density"], d["Khv"]
    nz = len(z)
    F = np.zeros(nz + 1)
    for k in range(nz - 1):
        Kh_face = 0.5 * (Khv[k] / rho[k] + Khv[k + 1] / rho[k + 1])
        F[k + 1] = -Kh_face * (thv[k + 1] - thv[k]) / (z[k + 1] - z[k])
    return F


def bl_depth(d, frac=0.05):
    km = d["Km"]
    i = np.where(km > frac * max(km.max(), 1e-12))[0]
    return d["z_phys"][i.max()] if len(i) else 0.0


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("cases", nargs="+")
    p.add_argument("--t0", default="plt72000")
    p.add_argument("--t1", default="plt86400")
    args = p.parse_args()

    print(f"Mixed-layer heat budget, integrated 0 to {ZTOP:.0f} m "
          f"(below erf.rayleigh_zdamp = 1000, so the damping term is zero)")
    print(f"Interval: {args.t0} -> {args.t1}.  All terms in W/m^2, positive = warms the column.\n")

    hdr = (f"{'case':<9}{'storage':>9}{'H_sfc':>8}{'F_top':>8}{'R_rad':>8}"
           f"{'resid':>8} | {'h':>6}{'dh/dt':>8}{'wthv_s':>9}{'wthv_h':>9}{'A':>7}{'dthv':>7}")
    print(hdr)
    print("-" * len(hdr))

    for c in args.cases:
        try:
            a, b = load(c, args.t0), load(c, args.t1)
        except Exception as exc:                      # noqa: BLE001
            print(f"{c:<9}  could not load: {exc}")
            continue
        dt = b["t"] - a["t"]
        m = a["z_phys"] < ZTOP

        # storage
        I0 = np.sum((a["density"] * a["theta"] * a["dz"])[m])
        I1 = np.sum((b["density"] * b["theta"] * b["dz"])[m])
        storage = CP * (I1 - I0) / dt

        # radiative sink over the same column (theta tendency x rho x cp)
        rho_dz = np.sum((0.5 * (a["density"] + b["density"]) * a["dz"])[m])
        R = CP * q_rad_of(c) * rho_dz

        # surface sensible flux, bulk_coeff with the mid-interval state
        th0 = 0.5 * (a["theta"][0] + b["theta"][0])
        qv0 = 0.5 * (a["qv"][0] + b["qv"][0])
        rho0 = 0.5 * (a["density"][0] + b["density"][0])
        p0 = 0.5 * (a["pressure"][0] + b["pressure"][0])
        # |U| is not in the plotfile as a face value; the cell-centered pair is exact here
        import yt as _yt                                    # noqa: F401
        wsp = _wsp(c, args.t1)
        H = rho0 * CP * CH * wsp * (THETA_SFC - th0)
        qs = q_sat_sea(THETA_SFC, p0)
        LE = rho0 * LV * CQ * wsp * (qs - qv0)

        # flux through the top of the control volume
        Fa, Fb = flux_profile(a), flux_profile(b)
        kt = int(np.argmin(np.abs(a["zf"] - ZTOP)))
        F_top = -0.5 * (Fa[kt] + Fb[kt]) * rho0 * CP      # into the column

        resid = storage - (H + F_top + R)

        # entrainment diagnostic
        h = bl_depth(b)
        dh = (bl_depth(b) - bl_depth(a)) / dt
        Fm = flux_profile(b)
        zf = b["zf"]
        inbl = (zf > 0.2 * h) & (zf < 1.6 * h)
        wthv_h = Fm[inbl].min() if inbl.any() else np.nan
        # surface virtual flux from the bulk formulae
        wth_s = CH * wsp * (THETA_SFC - b["theta"][0])
        wq_s = CQ * wsp * (q_sat_sea(THETA_SFC, b["pressure"][0]) - b["qv"][0])
        thv0 = b["thv"][0]
        wthv_s = wth_s * (1 + EPSV * b["qv"][0]) + EPSV * thv0 * wq_s
        A = -wthv_h / wthv_s if wthv_s != 0 else np.nan
        # inversion jump across the BL top
        kh = int(np.argmin(np.abs(b["z_phys"] - h)))
        k2 = min(kh + 3, len(b["z_phys"]) - 1)
        dthv = b["thv"][k2] - b["thv"][max(kh - 3, 0)]

        print(f"{c:<9}{storage:>9.1f}{H:>8.1f}{F_top:>8.1f}{R:>8.1f}{resid:>8.1f} | "
              f"{h:>6.0f}{dh*43200:>8.0f}{wthv_s:>9.5f}{wthv_h:>9.5f}{A:>7.2f}{dthv:>7.2f}")

    print()
    print("h      = BL depth [m] (K_m > 5% of peak);  dh/dt per 12 h")
    print("wthv_s = surface virtual heat flux [K m/s];  wthv_h = most negative flux near the BL top")
    print("A      = entrainment ratio -wthv_h/wthv_s   (classic convective BL: ~0.2)")
    print("dthv   = theta_v jump across the inversion [K]")


_WSP_CACHE = {}


def _wsp(case, plt):
    """|U| at the first cell center, which is what bulk_coeff_flux uses."""
    key = (case, plt)
    if key not in _WSP_CACHE:
        ds = yt.load(f"run_{case}/{plt}")
        cg = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
        u = np.asarray(cg[("boxlib", "x_velocity")].d)[0, 0, 0]
        v = np.asarray(cg[("boxlib", "y_velocity")].d)[0, 0, 0]
        _WSP_CACHE[key] = float(np.hypot(u, v))
    return _WSP_CACHE[key]


if __name__ == "__main__":
    sys.exit(main())
