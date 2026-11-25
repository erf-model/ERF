#!/usr/bin/env python3
"""Write a small wrfinput_d01 / wrfbdy_d01 pair on WRF's hybrid vertical coordinate.

The pair is made the way real.exe makes it, so that the wrfbdy arrays at the first time are
exactly the wrfinput state coupled with WRF's column mass and map factors:

  U_B*      = u   * (c1h mu_u + c2h) / msfuy        (couple(..., 'u', msfuy, ...))
  V_B*      = v   * (c1h mu_v + c2h) / msfvx        (couple(..., 'v', msfvx, ...))
  T_B*      = (theta_m - 300) * (c1h mu + c2h)
  QVAPOR_B* = qv  * (c1h mu + c2h)
  PH_B*     = ph' * (c1f mu + c2f)                  (couple(..., 'h', ...): full levels)
  MU_B*     = mu'

with mu = MU + MUB the dry column mass and mu_u, mu_v its average to the u and v points.

The state:
  - hybrid_opt = 2 with eta_c = 0.2 (c1 = B(eta), c2 = (1 - B)(p0 - p_top), WRF's B(eta)),
    or with --hybrid0 the terrain-following sigma coordinate (c1 = 1, c2 = 0);
  - a hill of 1200 m next to the west boundary and a perturbation of the column mass of
    up to 400 Pa, so that c1f mu + c2f differs from mu by tens of percent above eta_c;
  - map factors that vary smoothly over the domain, from 1 at the center to about 1.06 in
    the corners (isotropic: MAPFAC_UY = MAPFAC_U, MAPFAC_VX = MAPFAC_V);
  - physical (grid-relative) winds with shear, moist potential temperature (USE_THETA_M = 1)
    and water vapor varying in x, y and height;
  - a geopotential integrated hydrostatically from ALT, so that the density the boundary
    conversion takes from the geopotential, -(c1h mu + c2h) / (dph * rdnw), is 1/ALT.

  --no-cf             leave C1F and C2F out of wrfinput_d01 (ERF then falls back to c1f = 1,
                      c2f = 0)
  --hybrid0           the sigma coordinate, hybrid_opt = 0
  --bdy-without-msf   couple U_B* and V_B* without the 1/msfuy and 1/msfvx; ERF, which undoes
                      them, then gets boundary winds off by the map factor
"""

import argparse

import netCDF4 as nc
import numpy as np

NX, NY, NZ = 16, 12, 20
DX = 3000.0
BDY_WIDTH = 3
DATES = ["2026-01-01_00:00:00", "2026-01-01_06:00:00"]
NEXT_DATES = ["2026-01-01_06:00:00", "2026-01-01_12:00:00"]

GRAV = 9.81
R_D = 287.0
CP_D = 1004.5
RV_O_RD = 461.6 / 287.0
P0 = 1.0e5
P_TOP = 1.0e4
ETA_C = 0.2


def hybrid_b(eta, hybrid0):
    """WRF's B(eta) for hybrid_opt = 2 (dyn_em/module_initialize_real.F), or 1 for sigma."""
    if hybrid0:
        return np.ones_like(eta)
    d = (1.0 - ETA_C) ** 3
    b1 = 2.0 * ETA_C**2 / d
    b2 = -ETA_C * (4.0 + ETA_C + ETA_C**2) / d
    b3 = 2.0 * (1.0 + ETA_C + ETA_C**2) / d
    b4 = -(1.0 + ETA_C) / d
    b = b1 + b2 * eta + b3 * eta**2 + b4 * eta**3
    return np.where(eta > ETA_C, b, 0.0)


def build_state(hybrid0):
    lx, ly = NX * DX, NY * DX
    # Mass, x-face (u) and y-face (v) points
    xm = (np.arange(NX) + 0.5) * DX
    ym = (np.arange(NY) + 0.5) * DX
    xu = np.arange(NX + 1) * DX
    yv = np.arange(NY + 1) * DX

    def mapfac(x, y):
        return 1.0 + 0.06 * (((x - 0.5 * lx) / (0.5 * lx))**2 + ((y - 0.5 * ly) / (0.5 * ly))**2) / 2.0

    def hill(x, y):
        return 1200.0 * np.exp(-((x - 0.15 * lx)**2 + (y - 0.5 * ly)**2) / (2.0 * (9.0e3)**2))

    Xm, Ym = np.meshgrid(xm, ym)          # (NY, NX)
    Xu, Yu = np.meshgrid(xu, ym)          # (NY, NX+1)
    Xv, Yv = np.meshgrid(xm, yv)          # (NY+1, NX)

    s = {}
    s["MAPFAC_M"] = mapfac(Xm, Ym)
    s["MAPFAC_U"] = mapfac(Xu, Yu)
    s["MAPFAC_V"] = mapfac(Xv, Yv)
    s["HGT"] = hill(Xm, Ym)

    # Vertical coordinate: znw from 1 at the surface to 0 at the top, finer near the surface
    znw = 1.0 - (np.arange(NZ + 1) / NZ) ** 0.8
    znw[0], znw[-1] = 1.0, 0.0
    znu = 0.5 * (znw[:-1] + znw[1:])
    dnw = znw[1:] - znw[:-1]                 # negative
    c1f = hybrid_b(znw, hybrid0)
    c2f = (1.0 - c1f) * (P0 - P_TOP)
    c1h = hybrid_b(znu, hybrid0)
    c2h = (1.0 - c1h) * (P0 - P_TOP)
    s.update(ZNW=znw, ZNU=znu, DNW=dnw, RDNW=1.0 / dnw, C1F=c1f, C2F=c2f, C1H=c1h, C2H=c2h)

    # Column masses: the base state from the terrain, plus a perturbation
    mub = P0 * np.exp(-s["HGT"] / 7500.0) - P_TOP
    mu = 400.0 * np.sin(2.0 * np.pi * Xm / lx) * np.cos(2.0 * np.pi * Ym / ly)
    mut = mub + mu
    s["MUB"], s["MU"] = mub, mu
    s["PSFC"] = mut + P_TOP

    # Half-level dry hydrostatic pressure, full and base
    mu_h = c1h[:, None, None] * mut[None] + c2h[:, None, None]      # (NZ, NY, NX)
    mub_h = c1h[:, None, None] * mub[None] + c2h[:, None, None]
    p_h = mu_h + P_TOP
    pb_h = mub_h + P_TOP

    # Potential temperature and water vapor
    theta = (285.0 + 45.0 * (1.0 - p_h / P0)
             + 2.0 * np.cos(2.0 * np.pi * Xm / lx)[None] * (p_h / P0))
    theta_b = 285.0 + 45.0 * (1.0 - pb_h / P0)
    qv = 0.012 * (p_h / P0) ** 4 * (1.0 + 0.2 * np.sin(2.0 * np.pi * Ym / ly)[None])
    theta_m = theta * (1.0 + RV_O_RD * qv)

    # Inverse dry density from the equation of state at the hydrostatic pressure
    kappa_inv = (CP_D - R_D) / CP_D
    alt = (R_D / P0) * theta_m * (p_h / P0) ** (-kappa_inv)
    alb = (R_D / P0) * theta_b * (pb_h / P0) ** (-kappa_inv)

    # Geopotential: d(phi) = -alpha * (c1h mu + c2h) * dnw, from the terrain up
    phi = np.zeros((NZ + 1, NY, NX))
    phib = np.zeros((NZ + 1, NY, NX))
    phi[0] = GRAV * s["HGT"]
    phib[0] = GRAV * s["HGT"]
    for k in range(NZ):
        phi[k + 1] = phi[k] - alt[k] * mu_h[k] * dnw[k]
        phib[k + 1] = phib[k] - alb[k] * mub_h[k] * dnw[k]
    s["PHB"] = phib
    s["PH"] = phi - phib
    s["ALT"] = alt
    s["PB"] = pb_h
    s["P"] = p_h - pb_h
    s["THM"] = theta_m - 300.0
    s["QVAPOR"] = qv
    s["mu_h"] = mu_h

    # Physical winds with shear
    kk = np.arange(NZ)[:, None, None]
    s["U"] = 6.0 + 0.9 * kk + 2.0 * np.sin(2.0 * np.pi * Xu / lx)[None] + (Yu / ly)[None]
    s["V"] = -3.0 + 0.5 * kk + 1.5 * np.cos(2.0 * np.pi * Yv / ly)[None] + 0.5 * (Xv / lx)[None]
    s["W"] = np.zeros((NZ + 1, NY, NX))
    s["mut"] = mut
    return s


def write_wrfinput(path, s, no_cf, hybrid0):
    with nc.Dataset(path, "w", format="NETCDF4_CLASSIC") as ds:
        for name, size in (
            ("Time", 1), ("DateStrLen", len(DATES[0])),
            ("west_east", NX), ("west_east_stag", NX + 1),
            ("south_north", NY), ("south_north_stag", NY + 1),
            ("bottom_top", NZ), ("bottom_top_stag", NZ + 1),
        ):
            ds.createDimension(name, size)
        for name, value in (
            ("WEST-EAST_GRID_DIMENSION", NX + 1), ("SOUTH-NORTH_GRID_DIMENSION", NY + 1),
            ("BOTTOM-TOP_GRID_DIMENSION", NZ + 1), ("I_PARENT_START", 1), ("J_PARENT_START", 1),
            ("PARENT_GRID_RATIO", 1), ("USE_THETA_M", 1), ("HYBRID_OPT", 0 if hybrid0 else 2),
        ):
            ds.setncattr(name, np.int32(value))
        for name, value in (("DX", DX), ("DY", DX), ("ETAC", ETA_C)):
            ds.setncattr(name, np.float32(value))
        ds.setncattr("START_DATE", DATES[0])
        ds.setncattr("SIMULATION_START_DATE", DATES[0])

        t = ("Time",)
        xy = ("Time", "south_north", "west_east")
        xyz = ("Time", "bottom_top", "south_north", "west_east")
        zf = ("Time", "bottom_top_stag", "south_north", "west_east")

        def put(name, dims, data, dtype="f4"):
            var = ds.createVariable(name, dtype, dims)
            var[:] = np.broadcast_to(data, var.shape)

        put("Times", ("Time", "DateStrLen"), np.array(list(DATES[0]), dtype="S1"), "S1")
        put("U", ("Time", "bottom_top", "south_north", "west_east_stag"), s["U"][None])
        put("V", ("Time", "bottom_top", "south_north_stag", "west_east"), s["V"][None])
        put("W", zf, s["W"][None])
        for name in ("ALT", "THM", "PB", "P", "QVAPOR"):
            put(name, xyz, s[name][None])
        put("QCLOUD", xyz, 0.0)
        put("PH", zf, s["PH"][None])
        put("PHB", zf, s["PHB"][None])
        for name in ("MUB", "MU", "PSFC", "HGT", "MAPFAC_M"):
            put(name, xy, s[name][None])
        put("MAPFAC_U", ("Time", "south_north", "west_east_stag"), s["MAPFAC_U"][None])
        put("MAPFAC_V", ("Time", "south_north_stag", "west_east"), s["MAPFAC_V"][None])
        put("MAPFAC_UY", ("Time", "south_north", "west_east_stag"), s["MAPFAC_U"][None])
        put("MAPFAC_VX", ("Time", "south_north_stag", "west_east"), s["MAPFAC_V"][None])
        for name in ("C1H", "C2H", "RDNW", "ZNU", "DNW"):
            put(name, ("Time", "bottom_top"), s[name][None])
        put("ZNW", ("Time", "bottom_top_stag"), s["ZNW"][None])
        if not no_cf:
            put("C1F", ("Time", "bottom_top_stag"), s["C1F"][None])
            put("C2F", ("Time", "bottom_top_stag"), s["C2F"][None])
        put("P_TOP", t, P_TOP)
        for name, value in (("SST", 285.0), ("TSK", 285.0), ("LANDMASK", 1.0),
                            ("XLAT", 40.0), ("XLONG", -100.0)):
            put(name, xy, value)
        put("IVGTYP", xy, 10, "i4")
        put("ISLTYP", xy, 6, "i4")


def side_slices(s, with_msf=True):
    """The boundary arrays of each side at the first time, as WRF's couple() makes them.
    Each is indexed [w, k, lateral], w counting in from the side."""
    mut, mu_h = s["mut"], s["mu_h"]
    c1h, c2h, c1f, c2f = s["C1H"], s["C2H"], s["C1F"], s["C2F"]

    # Column mass at the u and v points; at the domain edge WRF takes the edge column
    muu = np.empty((NY, NX + 1))
    muu[:, 1:NX] = 0.5 * (mut[:, 1:] + mut[:, :-1])
    muu[:, 0], muu[:, NX] = mut[:, 0], mut[:, NX - 1]
    muv = np.empty((NY + 1, NX))
    muv[1:NY, :] = 0.5 * (mut[1:, :] + mut[:-1, :])
    muv[0, :], muv[NY, :] = mut[0, :], mut[NY - 1, :]

    msfuy = s["MAPFAC_U"] if with_msf else np.ones_like(s["MAPFAC_U"])
    msfvx = s["MAPFAC_V"] if with_msf else np.ones_like(s["MAPFAC_V"])
    Uc = s["U"] * (c1h[:, None, None] * muu[None] + c2h[:, None, None]) / msfuy[None]
    Vc = s["V"] * (c1h[:, None, None] * muv[None] + c2h[:, None, None]) / msfvx[None]
    Tc = s["THM"] * mu_h
    Qc = s["QVAPOR"] * mu_h
    PHc = s["PH"] * (c1f[:, None, None] * mut[None] + c2f[:, None, None])

    out = {}
    for w in range(BDY_WIDTH):
        # (x sides index the lateral direction by j, y sides by i)
        sides = {
            "XS": dict(U=Uc[:, :, w], V=Vc[:, :, w], T=Tc[:, :, w], Q=Qc[:, :, w],
                       PH=PHc[:, :, w], MU=s["MU"][:, w]),
            "XE": dict(U=Uc[:, :, NX - w], V=Vc[:, :, NX - 1 - w], T=Tc[:, :, NX - 1 - w],
                       Q=Qc[:, :, NX - 1 - w], PH=PHc[:, :, NX - 1 - w], MU=s["MU"][:, NX - 1 - w]),
            "YS": dict(U=Uc[:, w, :], V=Vc[:, w, :], T=Tc[:, w, :], Q=Qc[:, w, :],
                       PH=PHc[:, w, :], MU=s["MU"][w, :]),
            "YE": dict(U=Uc[:, NY - 1 - w, :], V=Vc[:, NY - w, :], T=Tc[:, NY - 1 - w, :],
                       Q=Qc[:, NY - 1 - w, :], PH=PHc[:, NY - 1 - w, :], MU=s["MU"][NY - 1 - w, :]),
        }
        for side, fields in sides.items():
            for name, data in fields.items():
                out.setdefault((side, name), []).append(data)
    return {key: np.array(val) for key, val in out.items()}


def write_wrfbdy(path, s, with_msf=True):
    slices = side_slices(s, with_msf)
    nt = len(DATES)
    with nc.Dataset(path, "w", format="NETCDF4_CLASSIC") as ds:
        for name, size in (
            ("Time", None), ("DateStrLen", len(DATES[0])), ("bdy_width", BDY_WIDTH),
            ("west_east", NX), ("west_east_stag", NX + 1),
            ("south_north", NY), ("south_north_stag", NY + 1),
            ("bottom_top", NZ), ("bottom_top_stag", NZ + 1),
        ):
            ds.createDimension(name, size)
        ds.setncattr("USE_THETA_M", np.int32(1))
        ds.setncattr("START_DATE", DATES[0])

        def chars(dates):
            return np.array([list(d) for d in dates], dtype="S1")

        ds.createVariable("Times", "S1", ("Time", "DateStrLen"))[:] = chars(DATES)
        ds.createVariable("md___thisbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_", "S1",
                          ("Time", "DateStrLen"))[:] = chars(DATES)
        ds.createVariable("md___nextbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_", "S1",
                          ("Time", "DateStrLen"))[:] = chars(NEXT_DATES)

        lateral = {"X": ("south_north", "south_north_stag"), "Y": ("west_east", "west_east_stag")}
        for side in ("XS", "XE", "YS", "YE"):
            mass, stag = lateral[side[0]]
            # U is staggered along y sides, V along x sides
            u_lat = stag if side[0] == "Y" else mass
            v_lat = stag if side[0] == "X" else mass
            fields = (
                ("U", ("bottom_top", u_lat), slices[(side, "U")]),
                ("V", ("bottom_top", v_lat), slices[(side, "V")]),
                ("T", ("bottom_top", mass), slices[(side, "T")]),
                ("QVAPOR", ("bottom_top", mass), slices[(side, "Q")]),
                ("PH", ("bottom_top_stag", mass), slices[(side, "PH")]),
                ("MU", (mass,), slices[(side, "MU")]),
                ("PC", (mass,), np.zeros_like(slices[(side, "MU")])),
            )
            for name, dims, data in fields:
                full_dims = ("Time", "bdy_width") + dims
                # The state does not change in time: the same values at both times, and
                # zero tendencies
                ds.createVariable(f"{name}_B{side}", "f4", full_dims)[:] = np.broadcast_to(
                    data[None], (nt,) + data.shape)
                ds.createVariable(f"{name}_BT{side}", "f4", full_dims)[:] = 0.0


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--no-cf", action="store_true", help="leave C1F and C2F out of wrfinput_d01")
    parser.add_argument("--hybrid0", action="store_true", help="the sigma coordinate (hybrid_opt = 0)")
    parser.add_argument("--bdy-without-msf", action="store_true",
                        help="couple U_B*, V_B* without the map factors")
    args = parser.parse_args()

    s = build_state(args.hybrid0)
    write_wrfinput("wrfinput_d01", s, args.no_cf, args.hybrid0)
    write_wrfbdy("wrfbdy_d01", s, with_msf=not args.bdy_without_msf)
    ratio = (s["C1F"][:, None, None] * s["mut"][None] + s["C2F"][:, None, None]) / s["mut"][None]
    print(f"hybrid_opt = {0 if args.hybrid0 else 2}; (c1f mu + c2f)/mu from {ratio.min():.3f} "
          f"to {ratio.max():.3f}; map factors from {s['MAPFAC_M'].min():.4f} to {s['MAPFAC_M'].max():.4f}")


if __name__ == "__main__":
    main()
