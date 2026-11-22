#!/usr/bin/env python3
"""Check erf.rans_shear_stencil of the k-eqn closure on a logarithmic wind profile.

Usage: check_log_shear.py --dt DT --cr CR [--stretched] <plotfile central> <plotfile log_law>

Both runs start from input_sounding_log, u = (u_* / kappa) ln((z + z0) / z0),
and take one step of DT seconds with erf.rans_realizability = CR, a value
large enough that the cap stays off. The two plotfiles are those after the
step; the plotfiles of step 0 beside them (same prefix, step 00000) hold the
velocity both the closure and the production of k saw.

1. The shear S the closure used. The closure writes
   cap_ratio = nu_t S / (C_r k) with nu_t = cmu sqrt(k) l, and sqrt(k) follows
   from diss = rho Cmu0^3 k^1.5 / l, so

       S = cap_ratio C_r sqrt(k) / (cmu l).

   Above the wall cell and below the top cell it must be
       central:  (u_k+1 - u_k-1) / (z_k+1 - z_k-1)
       log_law:  (u_k+1 - u_k-1) / [(z_k + z0) ln((z_k+1 + z0) / (z_k-1 + z0))]
   of the step-0 velocity, and for log_law also the closed form
   u_* / (kappa (z + z0)). The ratio central / log_law in the cells 2, 3 and 4
   is pinned to 1.2046, 1.0588 and 1.0285 (erf-model/ERF#5261). In the wall
   cell and the top cell the two runs must agree.

2. The shear production of k. Both runs have the same eddy viscosity in this
   step, so the difference of rho k after it is the difference of the
   production rho nu_t S^2 times DT, up to terms of order DT^2:

       (rho k_log - rho k_central) / (DT Kmv) = S_log^2 - S_central^2.

With --stretched the layers are not uniform: the heights of the cell centres
are read from the plotfile (z_phys), the central stencil is the mean of the
differences over the two faces of the cell, as ComputeStrain_S forms it, and
the checks against the closed form and the pinned ratios are left out (the
sounding is interpolated to the levels, so the initial wind is not the log
profile to round-off).

Guards against a vacuous pass: all fields finite, the cap off (max cap_ratio
below 1, above 0), and the shear of the two runs different by more than 10 %
in the second cell.
"""

import math
import os
import re
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, ".."), os.path.join(_here, "..", "..", "..", "Tools")]
import erf_plotfile  # noqa: E402
import rans_checks as rc  # noqa: E402

USTAR, KAPPA, Z0 = 0.4, 0.41, 0.1  # make_log_sounding.py, erf.most.z0
F_PINNED = (1.2046, 1.0588, 1.0285)
F0 = ["x_velocity", "y_velocity"]
F1 = ["density", "KE", "Kmv", "Lturb", "diss", "cmu", "cap_ratio"]


def step0_of(plt):
    return re.sub(r"\d{5}$", "00000", plt.rstrip("/"))


def main(argv):
    opts = {"--dt": None, "--cr": None}
    args = []
    stretched = "--stretched" in argv
    it = iter([a for a in argv if a != "--stretched"])
    for a in it:
        if a in opts:
            try:
                opts[a] = float(next(it))
            except (StopIteration, ValueError):
                print(__doc__)
                return 2
        else:
            args.append(a)
    dt, cr = opts["--dt"], opts["--cr"]
    if dt is None or cr is None or len(args) != 2 or not dt > 0 or not cr > 0:
        print(__doc__)
        return 2

    rep = rc.Report()
    runs = {}
    for label, plt in zip(("central", "log_law"), args):
        hdr0, f0 = erf_plotfile.read_fields(step0_of(plt), F0)
        hdr, f1 = erf_plotfile.read_fields(plt, F1 + (["z_phys"] if stretched else []))
        nx, ny, nz = [hdr["hi"][d] - hdr["lo"][d] + 1 for d in range(3)]
        if stretched:
            z = [f1["z_phys"][0][0][k] - hdr["prob_lo"][2] for k in range(nz)]
            ratio_z = (z[nz - 1] - z[nz - 2]) / (z[1] - z[0])
            rep.check("[%s] top over bottom layer spacing > 2 (stretched)" % label, ratio_z, 2.0, 0.0, "min")
        else:
            dz = (hdr["prob_hi"][2] - hdr["prob_lo"][2]) / nz
            z = [(k + 0.5) * dz for k in range(nz)]
        cols = [(i, j) for i in range(nx) for j in range(ny)]

        finite = all(math.isfinite(f1[fl][i][j][k]) for fl in F1 for (i, j) in cols for k in range(nz))
        rep.check("[%s] all fields finite" % label, 1.0 if finite else 0.0, 1.0, 0.0)
        ratio = [f1["cap_ratio"][i][j][k] for (i, j) in cols for k in range(nz)]
        rep.check("[%s] max cap_ratio < 1 (cap off)" % label, max(ratio), 1.0, 0.0, "max")
        rep.check("[%s] max cap_ratio > 0 (cap_ratio written)" % label, max(ratio), 1e-6, 0.0, "min")

        # S the closure used, and the two stencils on the step-0 velocity
        S = {}
        err, err_exact = 0.0, 0.0
        for (i, j) in cols:
            u, v = f0["x_velocity"][i][j], f0["y_velocity"][i][j]
            for k in range(nz):
                rho, lt = f1["density"][i][j][k], f1["Lturb"][i][j][k]
                sqrt_k = (f1["diss"][i][j][k] * lt / (rho * rc.CMU0 ** 3)) ** (1.0 / 3.0)
                S[(i, j, k)] = f1["cap_ratio"][i][j][k] * cr * sqrt_k / (f1["cmu"][i][j][k] * lt)
                if 1 <= k <= nz - 2:
                    du = math.hypot(u[k + 1] - u[k - 1], v[k + 1] - v[k - 1])
                    if label == "central" and stretched:
                        want = 0.5 * math.hypot(
                            (u[k + 1] - u[k]) / (z[k + 1] - z[k]) + (u[k] - u[k - 1]) / (z[k] - z[k - 1]),
                            (v[k + 1] - v[k]) / (z[k + 1] - z[k]) + (v[k] - v[k - 1]) / (z[k] - z[k - 1]))
                    elif label == "central":
                        want = du / (z[k + 1] - z[k - 1])
                    else:
                        want = du / ((z[k] + Z0) * math.log((z[k + 1] + Z0) / (z[k - 1] + Z0)))
                        exact = USTAR / (KAPPA * (z[k] + Z0))
                        err_exact = max(err_exact, abs(S[(i, j, k)] - exact) / exact)
                    err = max(err, abs(S[(i, j, k)] - want) / want)
        rep.check("[%s] max rel err of S vs the stencil on the step-0 wind" % label, err, 0.0, 1e-8, "max")
        if label == "log_law" and not stretched:
            rep.check("[log_law] max rel err of S vs u*/(kappa (z + z0))", err_exact, 0.0, 1e-8, "max")
        runs[label] = (cols, nz, f1, S)

    cols, nz, fa, Sa = runs["central"]
    _, _, fb, Sb = runs["log_law"]

    # the factor of the issue in the cells 2, 3, 4, and the cells the option leaves alone
    if stretched:
        fac = [Sa[(i, j, 1)] / Sb[(i, j, 1)] for (i, j) in cols]
        rep.check("S_central / S_log_law in cell 2 > 1.1 (the stencils differ)", min(fac), 1.1, 0.0, "min")
    for n, target in enumerate(() if stretched else F_PINNED):
        k = n + 1
        fac = [Sa[(i, j, k)] / Sb[(i, j, k)] for (i, j) in cols]
        rep.check("S_central / S_log_law in cell %d (min over columns)" % (k + 1), min(fac), target, 6e-5)
        rep.check("S_central / S_log_law in cell %d (max over columns)" % (k + 1), max(fac), target, 6e-5)
    for k, name in ((0, "wall cell"), (nz - 1, "top cell")):
        d = max(abs(Sa[(i, j, k)] - Sb[(i, j, k)]) / max(abs(Sa[(i, j, k)]), 1e-300) for (i, j) in cols)
        rep.check("max rel difference of S in the %s" % name, d, 0.0, 1e-9, "max")

    # production of k: difference of rho k after the step against Kmv (S_log^2 - S_central^2) DT
    err_p, dkmv = 0.0, 0.0
    for (i, j) in cols:
        for k in range(nz):
            dkmv = max(dkmv, abs(fa["Kmv"][i][j][k] - fb["Kmv"][i][j][k]) / fa["Kmv"][i][j][k])
            if 1 <= k <= nz - 2:
                got = (fb["density"][i][j][k] * fb["KE"][i][j][k] - fa["density"][i][j][k] * fa["KE"][i][j][k]) \
                    / (dt * fa["Kmv"][i][j][k])
                want = Sb[(i, j, k)] ** 2 - Sa[(i, j, k)] ** 2
                # relative to the change in the second cell, so cells far up, where the two
                # stencils agree to 1e-4, do not ask for a difference of round-off size
                scale = abs(Sb[(i, j, 1)] ** 2 - Sa[(i, j, 1)] ** 2)
                if abs(want) > 1e-2 * scale:
                    err_p = max(err_p, abs(got - want) / abs(want))
    rep.check("max rel difference of Kmv between the runs (same viscosity)", dkmv, 0.0, 1e-12, "max")
    rep.check("max rel err of d(rho k) / (DT Kmv) vs S_log^2 - S_central^2", err_p, 0.0, 1e-2, "max")

    rep.dump()
    return 1 if rep.failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
