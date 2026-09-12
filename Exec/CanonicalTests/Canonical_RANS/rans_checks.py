"""Shared pieces of the Canonical_RANS check scripts: the report table and
the structural checks every k-eqn RANS plotfile must satisfy.

Standard library only. Constants are the AL01 defaults; a deck that changes
them must pass the new values to the functions below.
"""

import math
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, "..")]
import erf_plotfile  # noqa: E402

KAPPA = 0.41
CMU0 = 0.5562
CB = 0.35
RT_MIN = -3.0
FIELDS = ["x_velocity", "y_velocity", "theta", "KE", "Kmv", "Khv",
          "Lturb", "walldist", "diss", "density"]


class Report:
    """Collects (name, measured, target, tolerance, kind, pass) rows."""

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
        elif kind == "range":   # target is (lo, hi); tol unused
            lo, hi = target
            err = max(lo - value, value - hi, 0.0)
            target = 0.5 * (lo + hi)
            tol = 0.5 * (hi - lo)
        ok = err <= tol
        if not ok:
            self.failed += 1
        self.rows.append((name, value, target, tol, kind, ok))

    def dump(self):
        print("%-44s %14s %14s %10s %5s %s" % ("check", "measured", "target", "tol", "kind", "pass"))
        for name, v, t, tol, kind, ok in self.rows:
            print("%-44s %14.6g %14.6g %10.3g %5s %s" % (name, v, t, tol, kind, "yes" if ok else "NO"))
        print("%d check(s) failed" % self.failed)


def geom_length(zc, z0, l_cap):
    """AL01 geometric length kappa (z + z0) with the harmonic cap."""
    lg = KAPPA * (zc + z0)
    return l_cap * lg / (l_cap + lg)


def unstable_bound(l_g):
    """Largest unstable length for the default Rt_min (about 1.31 l_g)."""
    return l_g * math.sqrt(1.0 + CMU0 ** 6 / CB ** 2 * abs(RT_MIN))


def parse_args(argv, doc):
    mode = "smoke"
    args = []
    for a in argv:
        if a in ("--smoke", "--physics"):
            mode = a[2:]
        else:
            args.append(a)
    if not args:
        print(doc)
        sys.exit(2)
    plt = args[0]
    surf = args[1] if len(args) > 1 else os.path.join(os.path.dirname(plt) or ".", "surf_hist.dat")
    return mode, plt, surf


def structural_checks(rep, z, p, z0, l_cap, allow_unstable, mode="smoke"):
    """Finite fields, positivity, exact wall distance, length-scale bounds,
    the Eq. 19 dissipation consistency in interior cells and the wall-cell
    k retention (both from the start-of-step diss and end-of-step KE).

    The dissipation check compares diss from the start of the last step with
    KE from its end, so it measures how much k changed in one step: 5 % on a
    converged state, 10 % during the early transient of a smoke run."""
    nz = len(z)
    finite = all(math.isfinite(v) for f in FIELDS for v in p[f])
    rep.check("all fields finite", 1.0 if finite else 0.0, 1.0, 0.0)
    rep.check("min KE [m2/s2]", min(p["KE"]), 0.0, 0.0, "min")
    rep.check("min Kmv [kg/m/s]", min(p["Kmv"]), 0.0, 0.0, "min")
    rep.check("min diss [kg/m/s3]", min(p["diss"]), 0.0, 0.0, "min")
    rep.check("max |walldist - z_cc| [m]", max(abs(p["walldist"][k] - z[k]) for k in range(nz)), 0.0, 1e-8)

    # The length never exceeds the neutral geometric length times the
    # unstable factor (1 when the case cannot be unstable anywhere).
    fac = unstable_bound(1.0) if allow_unstable else 1.0
    worst = max(p["Lturb"][k] / (fac * geom_length(z[k], z0, l_cap)) for k in range(nz))
    rep.check("max Lturb / bound", worst, 1.0, 1e-6, "max")

    def diss_expected(k):
        return p["density"][k] * CMU0 ** 3 * p["KE"][k] ** 1.5 / p["Lturb"][k]
    dmax = 0.0
    for k in range(1, nz):
        if p["KE"][k] < 1e-3:
            continue
        dmax = max(dmax, abs(p["diss"][k] - diss_expected(k)) / diss_expected(k))
    rep.check("max rel err diss vs AL01 Eq.19 (interior)", dmax, 0.0, 5e-2 if mode == "physics" else 1e-1)
    k_ratio = (p["diss"][0] / diss_expected(0)) ** (2.0 / 3.0)
    rep.check("wall cell k_start/k_end", k_ratio, 1.0, 0.01)


def bl_height_from_tke(z, ke):
    """Lowest height where KE drops below max(5 % of the wall value, 0.02)."""
    thr = max(0.05 * ke[0], 0.02)
    for k in range(1, len(z)):
        if ke[k] < thr:
            return z[k]
    return z[-1]
