#!/usr/bin/env python3
"""Check the realizability cap of the k-eqn closure on the ridge case.

Usage: check_realizability.py --cr CR <plotfile without cap> <plotfile with cap>

The first run has erf.rans_realizability = 0, the second erf.rans_realizability
= CR. In every cell the closure stores the capped eddy viscosity and
diffusivity, rho nu_t and rho nu_t', together with the length l, the
stability functions cmu and cmu', the dissipation and the diagnostic
cap_ratio = nu_t / (C_r k / S) of the uncapped closure. All of these come
from the same kernel at the start of the last step, while KE in the plotfile
is the k at the end of it, so the k the closure saw is recovered from the
dissipation, diss = rho Cmu0^3 k^1.5 / l (the density does not change in
time in this anelastic case). The checks pin the stored fields to the closed
form of the cap:

    sqrt(k) = (diss l / (rho Cmu0^3))^(1/3)
    Kmv = rho cmu  sqrt(k) l / max(cap_ratio, 1)
    Khv = rho cmu' sqrt(k) l / max(cap_ratio, 1)

cell by cell to round-off, so a cap that is skipped, applied with max in
place of min, applied to nu_t but not nu_t', or reported inverted, fails.
Without the cap the same relations hold with cap_ratio = 0 everywhere.

Guards against a vacuous pass: both runs finite, the cap clearly active in
the capped run (max cap_ratio > 1.2 and at least one cell above 1.05, well
above the 1.0012 of a neutral equilibrium cell at C_r = 0.309), and the two
runs different.
"""

import math
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, ".."), os.path.join(_here, "..", "..", "..", "Tools")]
import erf_plotfile  # noqa: E402
import rans_checks as rc  # noqa: E402

FIELDS = ["density", "Kmv", "Khv", "Lturb", "diss", "cmu", "cmu_prime", "cap_ratio"]


def main(argv):
    cr = None
    args = []
    it = iter(argv)
    for a in it:
        if a == "--cr":
            try:
                cr = float(next(it))
            except (StopIteration, ValueError):
                print(__doc__)
                return 2
        else:
            args.append(a)
    if cr is None or len(args) != 2:
        print(__doc__)
        return 2

    rep = rc.Report()
    runs = {}
    for label, plt in zip(("off", "on"), args):
        hdr, f = erf_plotfile.read_fields(plt, FIELDS)
        nx, ny, nz = [hdr["hi"][d] - hdr["lo"][d] + 1 for d in range(3)]
        cells = [(i, j, k) for i in range(nx) for j in range(ny) for k in range(nz)]
        ratio = [f["cap_ratio"][i][j][k] for (i, j, k) in cells]
        rep.check("[%s] min cap_ratio" % label, min(ratio), 0.0, 0.0, "min")
        if label == "off":
            rep.check("[off] max cap_ratio (cap off writes 0)", max(ratio), 0.0, 0.0)
        else:
            # C_r = 0.309 sits 0.1 % under Cmu0^2 = 0.30936, so a neutral equilibrium
            # cell alone reads cap_ratio = 1.0012; the guards ask for the crest shear
            # layer, well above that
            rep.check("[on] max cap_ratio (cap active in the crest shear layer)", max(ratio), 1.2, 0.0, "min")
            n_capped = sum(1 for r in ratio if r > 1.05)
            rep.check("[on] cells with cap_ratio > 1.05", float(n_capped), 1.0, 0.0, "min")
            print("[on] cap active in %d of %d cells (ratio > 1.05), max cap_ratio %.4g" % (n_capped, len(cells), max(ratio)))
        # closed form of the stored viscosity and diffusivity, and finiteness, in one pass
        err_m, err_h, finite = 0.0, 0.0, True
        for (i, j, k) in cells:
            finite = finite and all(math.isfinite(f[fl][i][j][k]) for fl in FIELDS)
            rho = f["density"][i][j][k]
            lt = f["Lturb"][i][j][k]
            sqrt_k = (f["diss"][i][j][k] * lt / (rho * rc.CMU0 ** 3)) ** (1.0 / 3.0)
            fac = 1.0 / max(f["cap_ratio"][i][j][k], 1.0)
            base = sqrt_k * lt * fac
            kmv = rho * f["cmu"][i][j][k] * base
            khv = rho * f["cmu_prime"][i][j][k] * base
            err_m = max(err_m, abs(f["Kmv"][i][j][k] - kmv) / max(abs(kmv), 1e-300))
            err_h = max(err_h, abs(f["Khv"][i][j][k] - khv) / max(abs(khv), 1e-300))
        rep.check("[%s] all fields finite" % label, 1.0 if finite else 0.0, 1.0, 0.0)
        rep.check("[%s] max rel err Kmv vs rho cmu sqrt(k) l / max(a,1)" % label, err_m, 0.0, 1e-9, "max")
        rep.check("[%s] max rel err Khv vs rho cmu' sqrt(k) l / max(a,1)" % label, err_h, 0.0, 1e-9, "max")
        runs[label] = (cells, f)

    # the two runs must differ (the cap changed the flow)
    cells, fa = runs["off"]
    _, fb = runs["on"]
    dkmv = max(abs(fa["Kmv"][i][j][k] - fb["Kmv"][i][j][k]) for (i, j, k) in cells)
    rep.check("max |Kmv_off - Kmv_on| > 0 [kg/m/s]", dkmv, 1e-300, 0.0, "min")

    rep.dump()
    return 1 if rep.failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
