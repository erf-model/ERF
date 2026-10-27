#!/usr/bin/env python3
"""Check the constant mass-flux relaxation against its closed form.

Usage: check_massflux_relaxation.py <log> --dt DT [--tau TAU] [--tol TOL] [--steps N]

<log> is the output of an anelastic run with erf.const_massflux_u set (and
erf.const_massflux_v = 0). make_mom_sources prints the integrated mass flux and its
target at each slow RHS evaluation:

    Integrated mass flux : F_x F_y (target: T_x T_y)

The control adds (T_x - F_x) / (A_x * tau) to the x momentum source, where A_x is the
cross-sectional area, so the mean rho*u relaxes towards T_x / A_x with time scale tau.
In uniform flow each RK2 step has two evaluations: the first at F_n, the second after a
forward-Euler predictor, F_n + (dt/tau) (T_x - F_n). The step average then gives

    T_x - F_{n+1} = (T_x - F_n) (1 - dt/tau + (dt/tau)^2 / 2)

so with F_0 = 0, F_n = T_x (1 - r^n) and r = 1 - dt/tau + (dt/tau)^2 / 2. Every printed
F_x must match these values to a relative --tol (relative to T_x), and F_y must stay
at zero to the same tolerance. A control that does not divide by the area, or divides
twice, changes the slope and fails. The check stops if the log has fewer than 2 * --steps
lines, or if T_x is zero (nothing to relax towards).
Exit status 0 on pass, 1 on a mismatch, 2 when the log cannot be checked.
"""

import argparse
import re
import sys

LINE = re.compile(r"Integrated mass flux\s*:\s*(\S+)\s+(\S+)\s*\(target:\s*(\S+)\s+(\S+)\)")


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("log")
    ap.add_argument("--dt", type=float, required=True)
    ap.add_argument("--tau", type=float, default=1.0)
    ap.add_argument("--tol", type=float, default=1e-8)
    ap.add_argument("--steps", type=int, default=3)
    args = ap.parse_args(argv)

    with open(args.log, errors="replace") as f:
        rows = [tuple(float(v) for v in m.groups()) for m in map(LINE.search, f) if m]
    if len(rows) < 2 * args.steps:
        print("found %d mass-flux lines, need %d" % (len(rows), 2 * args.steps))
        return 2
    target = rows[0][2]
    if target == 0.0:
        print("the x target is zero; the check is vacuous")
        return 2

    a = args.dt / args.tau
    r = 1.0 - a + 0.5 * a * a
    status = 0
    worst = 0.0
    for n in range(args.steps):
        f_n = target * (1.0 - r ** n)
        expected = (f_n, f_n + a * (target - f_n))
        for stage in range(2):
            fx, fy, tx, ty = rows[2 * n + stage]
            err = max(abs(fx - expected[stage]), abs(fy), abs(tx - target), abs(ty))
            rel = err / abs(target)
            worst = max(worst, rel)
            ok = rel <= args.tol
            print("step %d stage %d: F_x %.10g expected %.10g, F_y %.3g, relative error %.3e  %s"
                  % (n + 1, stage + 1, fx, expected[stage], fy, rel, "ok" if ok else "FAIL"))
            if not ok:
                status = 1
    print("largest relative error %.3e (tolerance %.1e)" % (worst, args.tol))
    print("PASS" if status == 0 else "FAIL")
    return status


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
