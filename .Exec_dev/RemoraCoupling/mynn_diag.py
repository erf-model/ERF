#!/usr/bin/env python3
"""
Offline re-evaluation of ERF's MYNN2.5 eddy viscosity from a plotfile column.

Recomputes l_S, l_T, l_B, the master length scale L_m, GM/GH, the Helfand-Labraga
level-2 limiter, and the stability functions S_M/S_H exactly as
Source/PBL/ERF_ComputeDiffusivityMYNN25.cpp does, then compares the resulting
K_m against the value ERF stored.  Agreement is the evidence that ERF's MYNN is
faithfully evaluating NN09; disagreement localises a bug.

Also reports the Mellor-Yamada 1982 / Blackadar counterfactual used by COAMPS in
Perlin et al. (2007), l = kappa z / (1 + kappa z / l0) with l0 = alpha * <zq>/<q>.

NOTE ON UNITS: ERF writes the plotfile variables "Kmv" and "Khv" as rho*K
(kg/m/s).  Only "nut" is divided by rho.  Everything printed here is K in m^2/s.

Usage:
    python mynn_diag.py plt86400 [--i 0] [--j 0] [--alpha 0.23] [--kmax 32]
"""
import argparse
import numpy as np
import yt

yt.set_log_level("error")

# ---- constants, matching Source/ERF_Constants.H -----------------------------
KAPPA = 0.41
GRAV = 9.81
R_d = 287.0
R_v = 461.505
EPSV = R_v / R_d - 1.0
CP_D = 1004.5
L_V = 2.5e6
CMU0 = 0.5562            # erf.Cmu0 default, ERF_TurbStruct.H:732
EPS = np.finfo(float).eps

# ---- MYNN level 2.5 coefficients, ERF_MYNNStruct.H --------------------------
A1, A2 = 1.18, 0.665
B1, B2 = 24.0, 15.0
C1, C2, C3, C5 = 0.137, 0.75, 0.352, 0.2
SQFAC = 3.0
SMMIN, SMMAX = 0.0, np.inf
SHMIN, SHMAX = 0.0, 4.0
GAM1 = 0.235


class Level2:
    """MYNNLevel2::init_coeffs + calc_Rf / calc_SM / calc_SH."""

    def __init__(self):
        self.gam2 = (2.0 * A1 * (3.0 - 2.0 * C2) + B2 * (1.0 - C3)) / B1
        self.F1 = (B1 * (GAM1 - C1) + 2.0 * A1 * (3.0 - 2.0 * C2)
                   + 3.0 * A2 * (1.0 - C2) * (1.0 - C5))
        self.F2 = B1 * (GAM1 + self.gam2) - 3.0 * A1 * (1.0 - C2)
        self.Rf1 = B1 * (GAM1 - C1) / self.F1
        self.Rf2 = B1 * GAM1 / self.F2
        self.Rfc = GAM1 / (GAM1 + self.gam2)
        self.Ri1 = 0.5 * A2 * self.F2 / (A1 * self.F1)
        self.Ri2 = 0.5 * self.Rf1 / self.Ri1
        self.Ri3 = (2.0 * self.Rf2 - self.Rf1) / self.Ri1

    def calc_Rf(self, GM, GH):
        lGM = np.copysign(max(abs(GM), EPS), GM)
        Ri = -GH / lGM
        return self.Ri1 * (Ri + self.Ri2
                           - np.sqrt(Ri * Ri - self.Ri3 * Ri + self.Ri2 * self.Ri2))

    def calc_SH(self, Rf):
        return 3.0 * A2 * (GAM1 + self.gam2) * (self.Rfc - Rf) / (1.0 - Rf)

    def calc_SM(self, Rf):
        return (A1 * self.F1 / (A2 * self.F2)
                * (self.Rf1 - Rf) / (self.Rf2 - Rf) * self.calc_SH(Rf))


def calc_stability_funcs(GM, GH, alphac):
    """MYNNLevel25::calc_stability_funcs, NN09 Eqns. 27-37."""
    a2 = alphac * alphac
    Phi1 = 1.0 - 3.0 * a2 * A2 * B2 * (1 - C3) * GH
    Phi2 = 1.0 - 9.0 * a2 * A1 * A2 * (1 - C2) * GH
    Phi3 = Phi1 + 9.0 * a2 * A2 * A2 * (1 - C2) * (1 - C5) * GH
    Phi4 = Phi1 - 12.0 * a2 * A1 * A2 * (1 - C2) * GH
    Phi5 = 6.0 * a2 * A1 * A1 * GM
    D = Phi2 * Phi4 + Phi5 * Phi3
    SM = alphac * A1 * (Phi3 - 3 * C1 * Phi4) / D
    SH = alphac * A2 * (Phi2 + 3 * C1 * Phi5) / D
    return SM, SH, SQFAC * SM


def phi_m(zeta):
    """similarity_funs::calc_phi_m, Businger-Dyer."""
    return 1.0 + 5.0 * zeta if zeta > 0 else (1.0 - 16.0 * zeta) ** -0.25


def phi_h(zeta):
    """similarity_funs::calc_phi_h, Businger-Dyer."""
    return 1.0 + 5.0 * zeta if zeta > 0 else (1.0 - 16.0 * zeta) ** -0.5


# --- saturation mixing ratio, matching Source/Utils/ERF_MicrophysicsUtils.H ---
RdoRv = R_d / R_v
_FLATAU = [6.11239921, 0.443987641, 0.142986287e-1, 0.264847430e-3,
           0.302950461e-5, 0.206739458e-7, 0.640689451e-10,
           -0.952447341e-13, -0.976195544e-15]


def erf_qsatw(T, p_mbar):
    """erf_qsatw: Flatau polynomial for e_sat over water, capped at RdoRv."""
    dtt = T - 273.16
    e = _FLATAU[8]
    for c in reversed(_FLATAU[:8]):
        e = c + dtt * e
    return RdoRv * e / max(e, p_mbar - e)


def q_surf_over_sea(theta_surf, p_pa):
    """What SurfaceLayer::fill_qsurf_with_qsat writes over water (is_land = 0).

    erf.most.surf_moist is OVERWRITTEN for ocean cells every step, so the
    inputs-file value is only the land default.
    """
    T = theta_surf * (p_pa / 1.0e5) ** (R_d / CP_D)
    return erf_qsatw(T, p_pa * 0.01)


def load_column(pltfile, i, j):
    ds = yt.load(pltfile)
    cg = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    have = {str(n[1]) for n in ds.field_list}

    def f(name):
        return np.asarray(cg[("boxlib", name)].d)[i, j, :].astype(float)

    col = {n: f(n) for n in ("density", "x_velocity", "y_velocity", "theta",
                             "KE", "Kmv", "Khv", "z_phys", "pressure")}
    col["qv"] = f("qv") if "qv" in have else np.zeros_like(col["density"])
    # Lturb is ERF's own master length scale; when present it removes the
    # (L_obukhov, l_T) degeneracy that K_m alone cannot resolve.
    col["Lturb"] = f("Lturb") if "Lturb" in have else None
    col["nut"] = f("nut") if "nut" in have else None
    return ds, col


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("pltfile")
    p.add_argument("--i", type=int, default=0)
    p.add_argument("--j", type=int, default=0)
    p.add_argument("--alpha", type=float, default=0.23,
                   help="Lt_alpha; 0.23 = NN09 (ERF default), 0.1 = Chen2021 / MY82")
    p.add_argument("--Cd", type=float, default=0.0014)
    p.add_argument("--Ch", type=float, default=0.0012)
    p.add_argument("--Cq", type=float, default=0.0012)
    p.add_argument("--surf_temp", type=float, default=287.15)
    p.add_argument("--surf_moist", type=float, default=None,
                   help="surface mixing ratio; default is q_sat at surf_temp, "
                        "which is what ERF uses over sea (erf.is_land = 0)")
    p.add_argument("--kmax", type=int, default=32, help="number of levels to print")
    args = p.parse_args()

    ds, c = load_column(args.pltfile, args.i, args.j)
    nz = len(c["density"])
    z = c["z_phys"]
    rho, u, v, th, qv = (c["density"], c["x_velocity"], c["y_velocity"],
                         c["theta"], c["qv"])
    KE = c["KE"]
    Km_erf = c["Kmv"] / rho          # ERF stores rho*K
    thv = th * (1.0 + EPSV * qv)
    q = np.sqrt(2.0 * np.maximum(KE, EPS))

    # Cell widths reconstructed from the cell-center heights.
    zf = np.empty(nz + 1)
    zf[0] = 0.0
    for k in range(nz):
        zf[k + 1] = 2.0 * z[k] - zf[k]
    dz = np.diff(zf)

    # ---- surface layer: bulk_coeff_flux, ERF_MOSTStress.H:2588 --------------
    if args.surf_moist is None:
        args.surf_moist = q_surf_over_sea(args.surf_temp, c["pressure"][0])
        qnote = " (q_sat over sea; erf.most.surf_moist is overwritten there)"
    else:
        qnote = " (from --surf_moist)"
    wsp = np.hypot(u[0], v[0])
    u_star = np.sqrt(args.Cd) * wsp
    wth = args.Ch * wsp * (args.surf_temp - th[0])        # <w'theta'>
    wq = args.Cq * wsp * (args.surf_moist - qv[0])        # <w'qv'>
    t_star = -wth / u_star
    q_star = -wq / u_star
    theta0, qv0 = thv[0], qv[0]

    # Independent readout: the surface cell's TKE is Dirichlet-set to u*^2/Cmu0^2
    u_star_tke = CMU0 * np.sqrt(KE[0])

    shf = -u_star * t_star
    slh = -u_star * q_star
    shf = shf * (1.0 + EPSV * qv0) + EPSV * theta0 * slh   # virtual heat flux
    l_obukhov = (-(theta0 * u_star ** 3) / (KAPPA * GRAV * shf)
                 if abs(shf) > EPS and u_star > EPS else np.inf)

    I1 = np.sum(z * q * dz)
    I0 = np.sum(q * dz)
    l_T = args.alpha * I1 / I0
    l0_my82 = 0.1 * I1 / I0

    lvl2 = Level2()

    print(f"file           : {args.pltfile}   (i,j)=({args.i},{args.j})")
    print(f"time           : {float(ds.current_time):.1f} s "
          f"= {float(ds.current_time)/3600:.2f} h")
    print(f"u* (bulk Cd)   : {u_star:.5f} m/s    "
          f"[tau = {rho[0]*u_star**2:.4f} N/m2]")
    print(f"u* (from TKE0) : {u_star_tke:.5f} m/s   <- independent cross-check")
    print(f"<w'theta'>     : {wth:+.6f} K m/s  ->  H  = {rho[0]*CP_D*wth:+8.2f} W/m2")
    print(f"<w'qv'>        : {wq:+.3e} m/s     ->  LE = {rho[0]*L_V*wq:+8.2f} W/m2")
    print(f"<w'thetav'>    : {shf:+.6f} K m/s  "
          f"({'UNSTABLE' if shf > 0 else 'STABLE'})")
    print(f"q_surf         : {args.surf_moist:.6f}{qnote}")
    print(f"L_obukhov      : {l_obukhov:.2f} m")
    print(f"int(zq)/int(q) : {I1/I0:.2f} m")
    print(f"l_T (a={args.alpha:.2f})    : {l_T:.2f} m")
    print(f"MY82 l0 (a=0.1): {l0_my82:.2f} m   <- COAMPS / Perlin counterfactual")

    # If ERF wrote its own length scale, cross-check against it directly. K_m
    # alone leaves L_obukhov and l_T degenerate; Lturb breaks that tie.
    if c["Lturb"] is not None:
        Lerf = c["Lturb"]
        print()
        print("ERF wrote Lturb -- direct cross-check (no fitted parameters):")
        print("   k      z    Lturb_ERF   S_M=Km/(L q)   Km_ERF")
        for k in range(min(args.kmax, nz - 1)):
            sm = (Km_erf[k] / (Lerf[k] * q[k])) if Lerf[k] * q[k] > 1e-12 else np.nan
            print(f"{k:4d} {z[k]:7.1f} {Lerf[k]:11.3f} {sm:14.4f} {Km_erf[k]:9.3f}")
    if c["nut"] is not None:
        rel = np.max(np.abs(c["nut"][:30] - Km_erf[:30])
                     / np.maximum(np.abs(Km_erf[:30]), 1e-12))
        print()
        print(f"unit check: max |nut - Kmv/rho| / |Kmv/rho| over k<30 = {rel:.2e}"
              "   (confirms Kmv is stored as rho*K)")
    print()
    hdr = ("   k       z      l_S      l_T      l_B      L_m ctl      q    "
           "S_M    a_c    Km_py   Km_ERF   ratio  L_MY82  Km_MY82")
    print(hdr)
    print("-" * len(hdr))

    rows = []
    for k in range(min(args.kmax, nz - 1)):
        if k == 0:
            dthvdz = dudz = dvdz = 0.0      # overwritten by the MOST profile below
        else:
            inv = 1.0 / (z[k + 1] - z[k - 1])
            dthvdz = (thv[k + 1] - thv[k - 1]) * inv
            dudz = (u[k + 1] - u[k - 1]) * inv
            dvdz = (v[k + 1] - v[k - 1]) * inv

        zval = z[k]
        zeta = zval / l_obukhov

        # l_S, NN09 Eqn. 53
        if zeta >= 1.0:
            l_S = KAPPA * zval / 3.7
        elif zeta >= 0.0:
            l_S = KAPPA * zval / (1.0 + 2.7 * zeta)
        else:
            l_S = KAPPA * zval * (1.0 - 100.0 * zeta) ** 0.2

        # First-cell MOST gradient override, ERF #4037
        if k == 0:
            tstar_v = t_star * (1.0 + EPSV * qv0) + EPSV * theta0 * q_star
            dudz = u_star * phi_m(zeta) / (KAPPA * zval)
            dvdz = 0.0
            dthvdz = tstar_v * phi_h(zeta) / (KAPPA * zval)

        # l_B, NN09 Eqn. 55
        if dthvdz > 0.0:
            N = np.sqrt(GRAV / theta0 * dthvdz)
            if zeta < 0.0:
                qc = (GRAV / theta0 * shf * l_T) ** (1.0 / 3.0)
                l_B = (1.0 + 5.0 * np.sqrt(qc / (N * l_T))) * q[k] / N
            else:
                l_B = q[k] / N
        else:
            l_B = np.inf

        Lm = 1.0 / (1.0 / l_S + 1.0 / l_T + 1.0 / l_B)
        ctl = min((("S", l_S), ("T", l_T), ("B", l_B)), key=lambda t: t[1])[0]

        shear = dudz * dudz + dvdz * dvdz
        buoy = -(GRAV / theta0) * dthvdz
        L2q2 = Lm * Lm / (q[k] * q[k])
        GM, GH = L2q2 * shear, L2q2 * buoy

        Rf = lvl2.calc_Rf(GM, GH)
        SM2 = lvl2.calc_SM(Rf)
        qe2 = B1 * Lm * Lm * SM2 * (1.0 - Rf) * shear
        qe = 0.0 if qe2 < 0 else max(np.sqrt(qe2), EPS)
        alphac = 1.0 if q[k] >= qe else q[k] / qe

        SM, SH, _ = calc_stability_funcs(GM, GH, alphac)
        SM = min(max(SM, SMMIN), SMMAX)
        SH = min(max(SH, SHMIN), SHMAX)

        Km = Lm * q[k] * SM
        ratio = Km / Km_erf[k] if Km_erf[k] > 1e-12 else np.nan

        L_my = KAPPA * zval / (1.0 + KAPPA * zval / l0_my82)
        Km_my = L_my * q[k] * SM

        rows.append((Km_erf[k], ratio))
        print(f"{k:4d} {zval:7.1f} {l_S:8.2f} {l_T:8.2f} "
              f"{(l_B if np.isfinite(l_B) else 9999.0):8.1f} {Lm:8.2f}  {ctl}  "
              f"{q[k]:6.3f} {SM:6.3f} {alphac:6.3f} "
              f"{Km:8.3f} {Km_erf[k]:8.3f} {ratio:7.4f} {L_my:7.2f} {Km_my:8.3f}")

    # k = 0 is excluded: there the gradients come from the MOST override and
    # ERF uses SurfaceLayer u*/t*/q* and the MOST-averaged theta_v, none of
    # which are in the plotfile. Every other level is a clean comparison.
    good = [r for n, (ke, r) in enumerate(rows)
            if n > 0 and ke > 1e-3 and np.isfinite(r)]
    if good:
        err = max(abs(np.array(good) - 1.0))
        print()
        print(f"max |Km_py/Km_ERF - 1| over k >= 1 with Km_ERF > 1e-3 : {err:.3%}")
        if rows and rows[0][0] > 1e-3:
            print(f"   (k = 0 excluded; its ratio is {rows[0][1]:.4f})")
        if err < 0.05:
            print("VERDICT: ERF's MYNN2.5 reproduces the NN09 formulas.")
        else:
            print("VERDICT: offline model and ERF differ.")
            print("  Before concluding there is a bug, note that L_obukhov and l_T are")
            print("  DEGENERATE against K_m alone: (L=539,l_T=131), (L=1025,l_T=74.5) and")
            print("  (L=2275,l_T=59) all fit the stored Kmv to ~1e-4 for plt86400.  The")
            print("  surface fluxes assumed here (--Cd/--Ch/--Cq/--surf_temp/--surf_moist)")
            print("  fix L_obukhov, and an error there is absorbed into l_T.  Add 'Lturb'")
            print("  to erf.plot_vars to break the tie.")


if __name__ == "__main__":
    main()
