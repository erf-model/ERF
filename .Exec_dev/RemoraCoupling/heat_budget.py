#!/usr/bin/env python3
"""
Mixed-layer heat budget and entrainment diagnostic for the Perlin comparison.

The problem is 1-D and has no resolved eddies, so the total turbulent flux IS the
subgrid flux, -K_h d(theta_v)/dz, evaluated on z-faces.  That makes the budget
closable from the plotfiles alone:

    d/dt int_0^zt rho*cp*theta dz  =  H_sfc + F_top + F_adv + R_rad

    H_sfc  surface sensible heat flux (bulk_coeff, absolute SST converted to theta)
    F_top  turbulent flux through the top of the column
    F_adv  the uniform large-scale heating, prob.advection_heating_rate, read from
           the run's job_info (the last value there is the one the run used)
    R_rad  RRTMGP heating, from the qsrc_sw + qsrc_lw plot variables averaged over
           every plotfile in the interval; reported as n/a (and left in the
           residual) when the run did not write them

with zt chosen below erf.rayleigh_zdamp so the Rayleigh damping term (which
targets the INITIAL profile and would otherwise act as an unaccounted source)
is identically zero.  R_rad is a time average of instantaneous heating rates, so
write plotfiles often (hourly: erf.plot_int_1 = 1200 at dt = 3 s) for the
residual to close; with 12-hourly output the diurnal cycle is badly sampled.

Usage:
    python heat_budget.py <run_dir> [<run_dir> ...] [--t0 plt72000 --t1 plt86400]

Each argument is a run directory, or <case> for a directory named run_<case>.
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
SST = 287.0                       # erf.most.abs_surf_temp, an absolute temperature [K]
ZTOP = 900.0                      # integrate below erf.rayleigh_zdamp = 1000

_FLATAU = [6.11239921, 0.443987641, 0.142986287e-1, 0.264847430e-3,
           0.302950461e-5, 0.206739458e-7, 0.640689451e-10,
           -0.952447341e-13, -0.976195544e-15]


def theta_sfc(p_pa):
    """Surface potential temperature ERF uses for an absolute SST at pressure p_pa."""
    return SST * (1.0e5 / p_pa) ** (R_D / CP)


def q_sat_sea(theta_surf, p_pa):
    """q_surf as SurfaceLayer::fill_qsurf_with_qsat sets it over water."""
    T = theta_surf * (p_pa / 1.0e5) ** (R_D / CP)
    dtt, pm = T - 273.16, p_pa * 0.01
    e = _FLATAU[8]
    for c in reversed(_FLATAU[:8]):
        e = c + dtt * e
    return (R_D / R_V) * e / max(e, pm - e)


def run_dir(case):
    """The run directory: <case> itself, or run_<case> for the old naming."""
    p = pathlib.Path(case)
    return p if p.is_dir() else pathlib.Path(f"run_{case}")


def heating_rate_of(case, plt):
    """prob.advection_heating_rate [K/s] the run used, from the plotfile's job_info.

    job_info lists the inputs file and then the command line, so the last value is
    the one in effect; zero is right for a run not given the forcing at all.
    """
    import re
    try:
        txt = (run_dir(case) / plt / "job_info").read_text()
    except OSError:
        return 0.0
    vals = re.findall(r"advection_heating_rate\s*=\s*(-?[\d.eE+-]+)", txt)
    return float(vals[-1]) if vals else 0.0


def interval_means(case, t0, t1, mask):
    """Time averages over every plotfile from t0 to t1 (trapezoid rule in time).

    Returns (rad, hsfc): rad is the rho-weighted column integral of qsrc_sw + qsrc_lw
    [K/s * kg/m^2], None when the run did not write them; hsfc is the bulk_coeff
    surface sensible heat flux [W/m^2].  The surface flux changes sign between day and
    night, so with frequent plotfiles this is much better than the endpoint value.
    """
    s0, s1 = int(t0[3:]), int(t1[3:])
    plts = sorted((p for p in run_dir(case).glob("plt*") if p.name[3:].isdigit()
                   and s0 <= int(p.name[3:]) <= s1), key=lambda p: int(p.name[3:]))
    times, vals, hs = [], [], []
    have_rad = True
    for p in plts:
        ds = yt.load(str(p))
        have_rad = have_rad and (("boxlib", "qsrc_sw") in ds.field_list)
        cg = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)

        def f(n):
            return np.asarray(cg[("boxlib", n)].d)[0, 0, :].astype(float)

        z = f("z_phys")
        zf = np.empty(len(z) + 1)
        zf[0] = 0.0
        for k in range(len(z)):
            zf[k + 1] = 2.0 * z[k] - zf[k]
        rho, pres = f("density"), f("pressure")
        if have_rad:
            q = (f("qsrc_sw") + f("qsrc_lw")) * rho * np.diff(zf)
            vals.append(np.sum(q[mask]))
        wsp = np.hypot(f("x_velocity")[0], f("y_velocity")[0])
        hs.append(rho[0] * CP * CH * wsp * (theta_sfc(pres[0]) - f("theta")[0]))
        times.append(float(ds.current_time))
    if len(times) < 2:
        return None, None
    trapezoid = getattr(np, "trapezoid", None) or np.trapz      # numpy < 2 has only trapz
    span = times[-1] - times[0]
    rad = trapezoid(vals, times) / span if have_rad else None
    return rad, trapezoid(hs, times) / span


def load(case, plt):
    ds = yt.load(str(run_dir(case) / plt))
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

    hdr = (f"{'case':<14}{'storage':>9}{'H_sfc':>8}{'F_top':>8}{'F_adv':>8}{'R_rad':>8}"
           f"{'resid':>8} | {'h':>6}{'dh/dt':>8}{'wthv_s':>9}{'wthv_h':>9}{'A':>7}{'dthv':>7}")
    print(hdr)
    print("-" * len(hdr))

    for c in args.cases:
        try:
            a, b = load(c, args.t0), load(c, args.t1)
        except Exception as exc:                      # noqa: BLE001
            print(f"{c:<14}  could not load: {exc}")
            continue
        dt = b["t"] - a["t"]
        m = a["z_phys"] < ZTOP

        # storage
        I0 = np.sum((a["density"] * a["theta"] * a["dz"])[m])
        I1 = np.sum((b["density"] * b["theta"] * b["dz"])[m])
        storage = CP * (I1 - I0) / dt

        # large-scale heating over the same column (theta tendency x rho x cp)
        rho_dz = np.sum((0.5 * (a["density"] + b["density"]) * a["dz"])[m])
        F_adv = CP * heating_rate_of(c, args.t1) * rho_dz

        # RRTMGP heating and the surface sensible flux, time-averaged over the
        # plotfiles in the interval
        q_rad, H_mean = interval_means(c, args.t0, args.t1, m)
        R = CP * q_rad if q_rad is not None else np.nan

        # surface sensible flux, bulk_coeff with the mid-interval state
        th0 = 0.5 * (a["theta"][0] + b["theta"][0])
        qv0 = 0.5 * (a["qv"][0] + b["qv"][0])
        rho0 = 0.5 * (a["density"][0] + b["density"][0])
        p0 = 0.5 * (a["pressure"][0] + b["pressure"][0])
        # |U| is not in the plotfile as a face value; the cell-centered pair is exact here
        wsp = _wsp(c, args.t1)
        H = H_mean if H_mean is not None else rho0 * CP * CH * wsp * (theta_sfc(p0) - th0)
        qs = q_sat_sea(theta_sfc(p0), p0)
        LE = rho0 * LV * CQ * wsp * (qs - qv0)

        # flux through the top of the control volume
        Fa, Fb = flux_profile(a), flux_profile(b)
        kt = int(np.argmin(np.abs(a["zf"] - ZTOP)))
        F_top = -0.5 * (Fa[kt] + Fb[kt]) * rho0 * CP      # into the column

        resid = storage - (H + F_top + F_adv + (0.0 if np.isnan(R) else R))

        # entrainment diagnostic
        h = bl_depth(b)
        dh = (bl_depth(b) - bl_depth(a)) / dt
        Fm = flux_profile(b)
        zf = b["zf"]
        inbl = (zf > 0.2 * h) & (zf < 1.6 * h)
        wthv_h = Fm[inbl].min() if inbl.any() else np.nan
        # surface virtual flux from the bulk formulae
        p_b = b["pressure"][0]
        wth_s = CH * wsp * (theta_sfc(p_b) - b["theta"][0])
        wq_s = CQ * wsp * (q_sat_sea(theta_sfc(p_b), p_b) - b["qv"][0])
        thv0 = b["thv"][0]
        wthv_s = wth_s * (1 + EPSV * b["qv"][0]) + EPSV * thv0 * wq_s
        A = -wthv_h / wthv_s if wthv_s != 0 else np.nan
        # inversion jump across the BL top
        kh = int(np.argmin(np.abs(b["z_phys"] - h)))
        k2 = min(kh + 3, len(b["z_phys"]) - 1)
        dthv = b["thv"][k2] - b["thv"][max(kh - 3, 0)]

        R_str = f"{R:>8.1f}" if not np.isnan(R) else f"{'n/a':>8}"
        print(f"{c:<14}{storage:>9.1f}{H:>8.1f}{F_top:>8.1f}{F_adv:>8.1f}{R_str}{resid:>8.1f} | "
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
        ds = yt.load(str(run_dir(case) / plt))
        cg = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
        u = np.asarray(cg[("boxlib", "x_velocity")].d)[0, 0, 0]
        v = np.asarray(cg[("boxlib", "y_velocity")].d)[0, 0, 0]
        _WSP_CACHE[key] = float(np.hypot(u, v))
    return _WSP_CACHE[key]


if __name__ == "__main__":
    sys.exit(main())
