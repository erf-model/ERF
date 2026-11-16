#!/usr/bin/env python3
"""Check the ABL warm bubble and the imposed bulk_coeff heat flux against closed forms.

Usage: check_warm_bubble.py <plotfile> [--self-test]

The deck (inputs) starts a calm, uniform 300 K atmosphere on a vertically
stretched mesh with one warm bubble (prob.bubble_*) and a bulk_coeff surface
layer that imposes the kinematic heat flux F = erf.most.surf_temp_flux. After a
short run the plotfile must show:

1. theta - 300 = bubble_dT cos^2(pi r / 2) at every cell above the bottom one,
   with r from each cell's own height z_phys (WRF's form). A bubble placed on a
   uniform (k + 1/2) dz misses by 0.48 K.
2. Compressible runs: the bubble is made at constant pressure, so the density is
   rho_hse theta_0 / (theta_0 + theta') and the pressure perturbation in the
   bubble is small (0.5 Pa). A bubble made by raising rho*theta at fixed density
   starts with gamma p theta'/theta = 230 Pa at the centre and no density deficit;
   after the 0.5 s run that leaves 14 Pa, a 1.7e-4 density error, and a heat
   content off by 3.2 times the surface gain.
   Runs with a fixed density (anelastic): the density is rho_hse.
3. The heat content sum(rho theta dV) has grown by sum(rho_1 F dA) t, rho_1 the
   bottom-cell density: the imposed flux reached the state. A layer that kept the
   Ch form with Ch = 0 adds nothing.

The bubble (prob.bubble_*) and F are read from the plotfile's job_info, so an edited
deck is checked against its own bubble; the committed deck's values stand in only
for a key job_info lacks, and the line "closed form from:" says which were used.

The exit code is the verdict.
"""
import math
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, "..", "..", "Tools")]
import erf_plotfile  # noqa: E402

THETA0 = 300.0
# The bubble and the flux are read from the run's job_info (read_run_params); these
# are the committed deck's values, used only for a key job_info does not hold and by
# the self-test.
DT_B, XC, YC, ZC, XR, YR, ZR = 0.5, 160.0, 150.0, 30.0, 80.0, 60.0, 20.0
FLUX = 0.24
_KEYS = {"prob.bubble_dT": "DT_B", "prob.bubble_x_c": "XC", "prob.bubble_y_c": "YC",
         "prob.bubble_z_c": "ZC", "prob.bubble_x_r": "XR", "prob.bubble_y_r": "YR",
         "prob.bubble_z_r": "ZR", "erf.most.surf_temp_flux": "FLUX"}

TOL_THETA = 2.0e-3       # K; the bubble moves by millimetres in 0.5 s
TOL_RHO = 2.0e-5         # relative; 5e-6 from the surface heating in 0.5 s
TOL_PRES = 5.0           # Pa; 0.5 Pa measured, 14 Pa for a fixed-density bubble
TOL_BUDGET = 0.02        # relative to the expected heat gain


def header_info(plotfile):
    with open(os.path.join(plotfile, "Header")) as f:
        lines = f.read().split("\n")
    nvars = int(lines[1])
    t = float(lines[2 + nvars + 1])
    prob_lo = [float(v) for v in lines[2 + nvars + 3].split()]
    prob_hi = [float(v) for v in lines[2 + nvars + 4].split()]
    return t, prob_lo, prob_hi


def job_info_values(plotfile):
    """Every "key = value ..." line of the plotfile's job_info, first token of the value."""
    vals = {}
    path = os.path.join(plotfile, "job_info")
    if not os.path.exists(path):
        return vals
    with open(path) as f:
        for line in f:
            s = line.strip()
            if "=" in s and not s.startswith("#"):
                k, v = s.split("=", 1)
                tok = v.split()
                if tok:
                    vals[k.strip()] = tok[0].strip().strip('"')
    return vals


def read_run_params(plotfile):
    """Set the bubble and flux from the run's own inputs (job_info), so the closed form
    checks the bubble the run made; print what was used and where it came from."""
    vals = job_info_values(plotfile)
    used = []
    for key, name in _KEYS.items():
        if key in vals:
            globals()[name] = float(vals[key])
            used.append("%s = %g" % (key, globals()[name]))
        else:
            used.append("%s = %g (not in job_info: the committed deck's value)" % (key, globals()[name]))
    print("closed form from: " + "; ".join(used))


def is_anelastic(plotfile):
    """Read erf.anelastic from the run's job_info (absent means compressible)."""
    path = os.path.join(plotfile, "job_info")
    if not os.path.exists(path):
        return False
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s.startswith("erf.anelastic"):
                val = s.split("=", 1)[1].split()[0].strip().strip('"')
                return val not in ("0", "false", "False")
    return False


def bubble(x, y, z):
    r = math.sqrt(((x - XC) / XR) ** 2 + ((y - YC) / YR) ** 2 + ((z - ZC) / ZR) ** 2)
    if r > 1.0:
        return 0.0
    return DT_B * math.cos(0.5 * math.pi * r) ** 2


def check(t, prob_lo, prob_hi, d, anelastic):
    rho, th, rh, pp, zc = d["density"], d["theta"], d["dens_hse"], d["pert_pres"], d["z_phys"]
    nx, ny, nz = len(rho), len(rho[0]), len(rho[0][0])
    dx = (prob_hi[0] - prob_lo[0]) / nx
    dy = (prob_hi[1] - prob_lo[1]) / ny

    err_th = err_rho = max_pp = 0.0
    max_bubble = 0.0
    n_inside = 0
    heat = heat0 = gain = 0.0
    for i in range(nx):
        x = prob_lo[0] + (i + 0.5) * dx
        for j in range(ny):
            y = prob_lo[1] + (j + 0.5) * dy
            znd = 0.0                      # flat ground at z = 0
            for k in range(nz):
                z = zc[i][j][k]
                dz = 2.0 * (z - znd)
                znd += dz
                dv = dx * dy * dz
                b = bubble(x, y, z)
                heat += rho[i][j][k] * th[i][j][k] * dv
                # the start: rho*theta = rho_hse*theta_0, plus rho_hse*theta' at fixed density
                heat0 += rh[i][j][k] * (THETA0 + (b if anelastic else 0.0)) * dv
                if k == 0:
                    gain += rho[i][j][k] * FLUX * dx * dy * t
                    continue               # the surface heat stays in the bottom cell
                if b > 0.0:
                    n_inside += 1
                    max_bubble = max(max_bubble, b)
                    if not anelastic:
                        max_pp = max(max_pp, abs(pp[i][j][k]))
                err_th = max(err_th, abs(th[i][j][k] - THETA0 - b))
                rho_exp = rh[i][j][k] * (1.0 if anelastic else THETA0 / (THETA0 + b))
                err_rho = max(err_rho, abs(rho[i][j][k] - rho_exp) / rho_exp)
    budget_err = abs((heat - heat0) - gain) / gain

    rows = [
        ("cells inside the bubble", n_inside, n_inside > 50),
        ("largest closed-form theta' [K]", max_bubble, max_bubble > 0.4),
        ("max |theta - 300 - theta'| [K]", err_th, err_th <= TOL_THETA),
        ("max relative density error", err_rho, err_rho <= TOL_RHO),
        ("heat gained / rho F A t - 1", budget_err, budget_err <= TOL_BUDGET),
    ]
    if not anelastic:
        rows.append(("max |p'| in the bubble [Pa]", max_pp, max_pp <= TOL_PRES))
    print("t = %.3f s, %s" % (t, "fixed density (anelastic)" if anelastic else "compressible"))
    print("%-34s %14s %6s" % ("check", "value", "pass"))
    failed = 0
    for name, val, ok in rows:
        failed += 0 if ok else 1
        print("%-34s %14.6g %6s" % (name, val, "yes" if ok else "NO"))
    print("%d check(s) failed" % failed)
    return failed


def self_test():
    """The checker must fail a fixed-density compressible bubble, a uniform-dz bubble and a lost flux."""
    nx, ny, nz = 16, 16, 12
    prob_lo, prob_hi = [0.0, 0.0, 0.0], [320.0, 320.0, 116.0]
    dx = dy = 20.0
    t = 0.5
    zn = [0.0]
    for k in range(nz):
        zn.append(zn[-1] + 4.0 * 1.15 ** k)

    def fields(kind):
        d = {n: [[[0.0] * nz for _ in range(ny)] for _ in range(nx)]
             for n in ("density", "theta", "dens_hse", "pert_pres", "z_phys")}
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    z = 0.5 * (zn[k] + zn[k + 1])
                    zu = (k + 0.5) * (prob_hi[2] / nz)
                    b = bubble((i + 0.5) * dx, (j + 0.5) * dy, zu if kind == "uniform_dz" else z)
                    rh = math.exp(-z / 8000.0)
                    d["z_phys"][i][j][k] = z
                    d["dens_hse"][i][j][k] = rh
                    if kind == "fixed_rho":
                        d["density"][i][j][k] = rh
                        d["pert_pres"][i][j][k] = 1.4 * 1.0e5 * b / THETA0
                    else:
                        d["density"][i][j][k] = rh * THETA0 / (THETA0 + b)
                    d["theta"][i][j][k] = THETA0 + b
                    if k == 0 and kind != "no_flux":
                        # the surface heat, all in the bottom cell
                        dz = zn[1] - zn[0]
                        d["theta"][i][j][k] += FLUX * t / dz
        return d

    ok = check(t, prob_lo, prob_hi, fields("good"), False) == 0
    bad = [check(t, prob_lo, prob_hi, fields(kind), False) > 0
           for kind in ("fixed_rho", "uniform_dz", "no_flux")]
    # the parameters come from job_info: a run whose bubble sits 10 m higher must be
    # checked against that bubble, and the good fields against it must fail
    import tempfile
    good = fields("good")              # made with the committed deck's bubble
    saved = {n: globals()[n] for n in _KEYS.values()}
    with tempfile.TemporaryDirectory() as tmp:
        with open(os.path.join(tmp, "job_info"), "w") as f:
            f.write("prob.bubble_z_c = %g\nerf.most.surf_temp_flux = %g\n" % (ZC + 10.0, FLUX))
        read_run_params(tmp)
    parsed = (ZC == saved["ZC"] + 10.0 and DT_B == saved["DT_B"])
    moved = check(t, prob_lo, prob_hi, good, False) > 0
    globals().update(saved)
    all_ok = ok and all(bad) and parsed and moved
    print("self-test:", "passed" if all_ok else "FAILED")
    return 0 if all_ok else 1


def main(argv):
    if "--self-test" in argv:
        return self_test()
    args = [a for a in argv[1:] if not a.startswith("--")]
    if len(args) != 1:
        print(__doc__)
        return 2
    t, prob_lo, prob_hi = header_info(args[0])
    if t <= 0.0:
        print("the plotfile is at t = 0: the heat-budget check needs a later one")
        return 1
    read_run_params(args[0])
    _, d = erf_plotfile.read_fields(args[0], ["density", "theta", "dens_hse", "pert_pres", "z_phys"])
    return 1 if check(t, prob_lo, prob_hi, d, is_anelastic(args[0])) else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
