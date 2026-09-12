#!/usr/bin/env python3
"""Numerical checks for the neutral flat-ABL k-eqn RANS case.

Usage: check_neutral.py [--smoke | --physics] <plotfile> [surf_hist.dat]

Every check prints measured value, target and tolerance; the exit code is
non-zero if any enabled check fails. --smoke (default) runs the structural
checks that must hold after a few steps; --physics adds the checks that
need a converged Ekman layer (the 12 h run in the README).
"""

import math
import os
import sys

# erf_plotfile.py lives in Canonical_RANS/ next to the case directories; the
# CTest copies it beside this script, so look in both places.
_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, "..")]
import erf_plotfile  # noqa: E402

KAPPA = 0.41
CMU0 = 0.5562
Z0 = 0.1
L_G_MAX = 30.0
FIELDS = ["x_velocity", "y_velocity", "theta", "KE", "Kmv", "Khv",
          "Lturb", "walldist", "diss", "density", "Rt", "cmu", "cmu_prime"]


class Report:
    def __init__(self):
        self.rows = []
        self.failed = 0

    def check(self, name, value, target, tol, kind="abs"):
        if kind == "abs":
            err = abs(value - target)
        elif kind == "rel":
            err = abs(value - target) / max(abs(target), 1e-300)
        elif kind == "min":     # value must be >= target
            err = max(target - value, 0.0)
        elif kind == "max":     # value must be <= target
            err = max(value - target, 0.0)
        ok = err <= tol
        if not ok:
            self.failed += 1
        self.rows.append((name, value, target, tol, kind, ok))

    def dump(self):
        print("%-40s %14s %14s %10s %5s %s" % ("check", "measured", "target", "tol", "kind", "pass"))
        for name, v, t, tol, kind, ok in self.rows:
            print("%-40s %14.6g %14.6g %10.3g %5s %s" % (name, v, t, tol, kind, "yes" if ok else "NO"))
        print("%d check(s) failed" % self.failed)


def main(argv):
    mode = "smoke"
    args = []
    for a in argv:
        if a in ("--smoke", "--physics"):
            mode = a[2:]
        else:
            args.append(a)
    if not args:
        print(__doc__)
        return 2
    plt = args[0]
    surf = args[1] if len(args) > 1 else os.path.join(os.path.dirname(plt) or ".", "surf_hist.dat")

    z, p, hdr = erf_plotfile.planar_averages(plt, FIELDS)
    nz = len(z)
    rep = Report()

    # --- structural checks (always) -------------------------------------
    finite = all(math.isfinite(v) for f in FIELDS for v in p[f])
    rep.check("all fields finite", 1.0 if finite else 0.0, 1.0, 0.0)
    rep.check("min KE [m2/s2]", min(p["KE"]), 0.0, 0.0, "min")
    rep.check("min Kmv [kg/m/s]", min(p["Kmv"]), 0.0, 0.0, "min")
    rep.check("min diss [kg/m/s3]", min(p["diss"]), 0.0, 0.0, "min")

    # wall distance on a flat constant-dz mesh is the cell-centre height
    wd_err = max(abs(p["walldist"][k] - z[k]) for k in range(nz))
    rep.check("max |walldist - z_cc| [m]", wd_err, 0.0, 1e-8)

    # AL01 geometric length with the harmonic cap, exact in neutral air
    def l_geom(zc):
        lg = KAPPA * (zc + Z0)
        return L_G_MAX * lg / (L_G_MAX + lg)
    for k in (0, 1, 2):
        rep.check("Lturb(k=%d) vs capped kappa(z+z0)" % k, p["Lturb"][k], l_geom(z[k]), 1e-3, "rel")
    rep.check("max Lturb <= max_geom_lscale", max(p["Lturb"]), L_G_MAX, 1e-9, "max")

    # stability functions written by the closure must be AL01 Eqs. 31-32 of
    # the written (smoothed) Rt, and Rt must sit above Rt_min
    def cmu_of(Rt):
        return (CMU0 + 0.108 * Rt) / (1.0 + 0.308 * Rt + 0.00837 * Rt * Rt)
    def cmu_prime_of(Rt):
        return CMU0 / (1.0 + 0.277 * Rt)
    e1 = max(abs(p["cmu"][k] - cmu_of(p["Rt"][k])) for k in range(nz))
    e2 = max(abs(p["cmu_prime"][k] - cmu_prime_of(p["Rt"][k])) for k in range(nz))
    rep.check("max |cmu - Eq.31(Rt)|", e1, 0.0, 1e-6)
    rep.check("max |cmu_prime - Eq.32(Rt)|", e2, 0.0, 1e-6)
    rep.check("min Rt >= Rt_min", min(p["Rt"]), -3.0, 1e-12, "min")
    # Kmv must be rho cmu sqrt(k) Lturb (mean of the planar averages, so a loose check)
    rep.check("Kmv(k=1) vs rho cmu sqrt(KE) Lturb", p["Kmv"][1], p["density"][1] * p["cmu"][1] * math.sqrt(p["KE"][1]) * p["Lturb"][1], 0.05, "rel")

    # Dissipation follows Cmu0^3 rho k^1.5 / L (AL01 Eq. 19). The plotfile
    # holds diss from the start of the last step and KE from its end, so the
    # check is restricted to interior cells with resolved k (above the floor
    # region at the domain top) and allows for one step of k evolution.
    def diss_expected(k):
        return p["density"][k] * CMU0 ** 3 * p["KE"][k] ** 1.5 / p["Lturb"][k]
    dmax = 0.0
    for k in range(1, nz):
        if p["KE"][k] < 1e-3:
            continue
        dmax = max(dmax, abs(p["diss"][k] - diss_expected(k)) / diss_expected(k))
    rep.check("max rel err diss vs AL01 Eq.19 (interior)", dmax, 0.0, 5e-2)

    # In the wall cell the same ratio measures how much of the Dirichlet k
    # imposed at the start of the step survives to its end:
    #   (diss / diss_expected)^(2/3) = k_start / k_end.
    # A consistent wall condition keeps this at one.
    k_ratio = (p["diss"][0] / diss_expected(0)) ** (2.0 / 3.0)
    rep.check("wall cell k_start/k_end", k_ratio, 1.0, 0.01)

    # --- physics checks (converged run) ---------------------------------
    if mode == "physics":
        sh = erf_plotfile.read_surf_hist(surf)
        ustar = sh["u_star"]
        rep.check("u_star [m/s] plausible", ustar, 0.35, 0.15)
        # Wall k from the neutral AL01 condition, k = u*^2 / Cmu0^2
        rep.check("KE(k=0)/u*^2 vs 1/Cmu0^2", p["KE"][0] / ustar ** 2, 1.0 / CMU0 ** 2, 0.05, "rel")
        # log law in the lowest four cells (z < 0.1 zi for zi ~ 700 m)
        for k in range(4):
            speed = math.hypot(p["x_velocity"][k], p["y_velocity"][k])
            u_log = ustar / KAPPA * math.log((z[k] + Z0) / Z0)
            rep.check("|U|(k=%d) vs log law" % k, speed, u_log, 0.10, "rel")
        # eddy viscosity near the wall: rho kappa u* z within a factor 1.3
        rho0 = p["density"][0]
        rep.check("Kmv(k=1)/(rho kappa u* z)", p["Kmv"][1] / (rho0 * KAPPA * ustar * (z[1] + Z0)), 1.0, 0.3)

    rep.dump()
    return 1 if rep.failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
