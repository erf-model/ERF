#!/usr/bin/env python3
"""Checks for erf.input_sounding_theta_above_ground (inputs_theta).

Usage: check_theta_above_ground.py [--smoke] <plotfile>

The flag sets the initial theta, so the comparisons read the step-0 plotfile
(plt00000 next to the given one). ERF resamples an input sounding onto the flat
mesh's cell-centre heights, and the ground and top of the domain, before it
interpolates; the reference here is the same table, looked up at each cell's
height above the local ground (the column that inflow_profile_inversion.txt
imposes on the inflow face). The checks:

- step 0: theta matches that reference in every cell and in the first column
  next to the inflow face (max error 1e-6 K);
- step 0: every column is in discrete hydrostatic balance,
  |dp/dz + g (rho_k + rho_k+1)/2| <= 1e-9 rho_k g;
- the comparison is only evidence if the inversion moves with the ground: the
  reference at physical height must differ from the one above the ground by
  at least 2 K in at least a tenth of the cells;
- the given plotfile (after a few steps): every field is finite and |w| stays
  below 5 m/s (a 10 m/s wind over the 0.15 slope is lifted at 1.5 m/s).

check_theta_flag_off.py runs the step-0 comparison on the flag-off deck and
requires the error to be large instead.
"""

import math
import os
import re
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, "..")]
import erf_plotfile  # noqa: E402
import rans_checks as rc  # noqa: E402

Z_INV, DTH_INV, DZ_INV, LAPSE = 300.0, 8.0, 100.0, 0.003
GRAVITY = 9.81


def sounding_theta(z):
    """The sounding of gen_theta_inputs.py (piecewise linear in height above the ground)."""
    z = max(z, 0.0)
    if z <= Z_INV:
        return 300.0
    if z <= Z_INV + DZ_INV:
        return 300.0 + DTH_INV * (z - Z_INV) / DZ_INV
    return 300.0 + DTH_INV + LAPSE * (z - Z_INV - DZ_INV)


def resampled_table(nz, top):
    """ERF's InputSoundingData table: ground, the flat cell centres, the top."""
    dz = top / nz
    zs = [0.0] + [(k + 0.5) * dz for k in range(nz)] + [top]
    return zs, [sounding_theta(z) for z in zs]


def lookup(zs, fs, z):
    if z <= zs[0]:
        return fs[0]
    if z >= zs[-1]:
        return fs[-1]
    for n in range(1, len(zs)):
        if z <= zs[n]:
            w = (z - zs[n - 1]) / (zs[n] - zs[n - 1])
            return fs[n - 1] + w * (fs[n] - fs[n - 1])
    return fs[-1]


def median(vals):
    s = sorted(vals)
    n = len(s)
    return 0.5 * (s[(n - 1) // 2] + s[n // 2])


def step0_plotfile(plotfile):
    path = plotfile.rstrip("/")
    return re.sub(r"\d{5}$", "00000", path)


def measure(plotfile):
    """Step-0 theta errors and balance, plus finiteness and |w| of the given plotfile."""
    _, d0 = erf_plotfile.read_fields(step0_plotfile(plotfile), ["theta", "pressure", "density", "z_phys"])
    th, p, r, zp = (d0[f] for f in ("theta", "pressure", "density", "z_phys"))
    nx, ny, nz = len(th), len(th[0]), len(th[0][0])
    top = zp[0][0][nz - 1] + 0.5 * (zp[0][0][nz - 1] - zp[0][0][nz - 2])
    zs, fs = resampled_table(nz, top)
    err, err_first, shift, residual = [], [], [], []
    for i in range(nx):
        for j in range(ny):
            col = zp[i][j]
            ground = col[0] - 0.5 * (col[1] - col[0])
            for k in range(nz):
                agl = col[k] - ground
                ref = lookup(zs, fs, agl)
                e = abs(th[i][j][k] - ref)
                err.append(e)
                if i == 0:
                    err_first.append(e)
                shift.append(abs(lookup(zs, fs, col[k]) - ref))
            for k in range(nz - 1):
                dpdz = (p[i][j][k + 1] - p[i][j][k]) / (col[k + 1] - col[k])
                residual.append(abs(dpdz + GRAVITY * 0.5 * (r[i][j][k] + r[i][j][k + 1])) / (GRAVITY * r[i][j][k]))
    _, d = erf_plotfile.read_fields(plotfile, ["theta", "z_velocity"])
    th_final, w = d["theta"], d["z_velocity"]
    finite = all(math.isfinite(f[i][j][k]) for f in (th_final, w)
                 for i in range(len(f)) for j in range(len(f[0])) for k in range(len(f[0][0])))
    wmax = max(abs(w[i][j][k]) for i in range(len(w)) for j in range(len(w[0])) for k in range(len(w[0][0])))
    return dict(err=err, err_first=err_first, shift=shift, residual=residual, finite=finite, wmax=wmax)


def main(argv):
    args = [a for a in argv[1:] if not a.startswith("--")]
    if len(args) != 1:
        print(__doc__)
        return 2
    m = measure(args[0])
    rep = rc.Report()
    print("step 0 from %s: %d cells, %d in the first column" % (step0_plotfile(args[0]), len(m["err"]), len(m["err_first"])))
    rep.check("step 0: cells where theta at physical height differs by >= 2 K [fraction]",
              sum(1 for s in m["shift"] if s >= 2.0) / len(m["shift"]), 0.1, 0.0, "min")
    rep.check("step 0: max |theta - theta(z above ground)| [K]", max(m["err"]), 0.0, 1e-6, "max")
    rep.check("step 0: first column max |theta - theta(z above ground)| [K]", max(m["err_first"]), 0.0, 1e-6, "max")
    rep.check("step 0: max hydrostatic residual / (rho g)", max(m["residual"]), 0.0, 1e-9, "max")
    rep.check("final plotfile: all fields finite", 1.0 if m["finite"] else 0.0, 1.0, 0.0)
    rep.check("final plotfile: max |w| [m/s]", m["wmax"], 0.0, 5.0, "max")
    rep.dump()
    return 1 if rep.failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
