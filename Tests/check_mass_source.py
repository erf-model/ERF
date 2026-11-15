#!/usr/bin/env python3
"""Check the budget of a run with a mass source (erf.mass_source).

Usage: check_mass_source.py <log> (--rate RATE | --expected MASS) [--theta THETA]
                            [--no-tracer] [--tol TOL]

<log> is the output of a run in a closed domain (nothing enters or leaves), with
erf.sum_interval = 1 and erf.v >= 1, so that ERF prints at every step

    TIME= t
     MASS       = M            (one level), or  MASS       SL/ML = M_sl M  (several)
     RHO THETA  = R
     RHO SCALAR = S

The multilevel (last) value is used when both are printed. The mass injected by the last
time printed, E, is RATE * t for a constant rate, or MASS (the total of a rate table that
ends before that time) with --expected. The mass must grow by E, and, with the initial
scalar zero, the scalar integral must equal E, or stay zero with --no-tracer. With --theta
the rho*theta integral must grow by THETA * E (in a closed domain without diffusion, where
the injected gas carries THETA, or where theta is uniform and equal to THETA). Each to a
relative --tol of E. A source added at the wrong stage time, normalized by the wrong volume
or lost by the substeps changes one or more of these. The check also requires the source to
have been on for more than one step.
Exit status 0 on pass, 1 on a mismatch, 2 when the log cannot be checked.
"""

import argparse
import re
import sys


def last_value(line):
    return float(line.split("=")[1].split()[-1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("log")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--rate", type=float)
    g.add_argument("--expected", type=float)
    ap.add_argument("--theta", type=float)
    ap.add_argument("--no-tracer", action="store_true")
    ap.add_argument("--tol", type=float, default=1e-10)
    a = ap.parse_args()

    times, mass, rhth, scal = [], [], [], []
    with open(a.log) as f:
        for line in f:
            if line.startswith("TIME="):
                times.append(float(line.split("=")[1]))
            elif re.match(r"^ *MASS +(SL/ML +)?=", line):
                mass.append(last_value(line))
            elif re.match(r"^ *RHO THETA +(SL/ML +)?=", line):
                rhth.append(last_value(line))
            elif re.match(r"^ *RHO SCALAR +(SL/ML +)?=", line):
                scal.append(last_value(line))

    n = len(times)
    if not (n == len(mass) == len(rhth) == len(scal)) or n < 3:
        print("check_mass_source: found %d TIME, %d MASS, %d RHO THETA and %d RHO SCALAR lines in %s; "
              "need the same number, at least 3 (erf.sum_interval = 1, erf.v >= 1, transport_scalar on?)"
              % (n, len(mass), len(rhth), len(scal), a.log))
        return 2
    if times[0] != 0.0:
        print("check_mass_source: the first integrals are at t = %g, not 0" % times[0])
        return 2

    expected = a.rate * times[-1] if a.rate is not None else a.expected
    errors = {}
    errors["mass gain"] = (mass[-1] - mass[0], expected)
    errors["tracer"] = (scal[-1], 0.0 if a.no_tracer else expected)
    if a.theta is not None:
        errors["rho theta gain"] = (rhth[-1] - rhth[0], a.theta * expected)

    ok = True
    if scal[0] != 0.0:
        print("check_mass_source: FAIL: the tracer integral is %g at t = 0, not 0" % scal[0])
        ok = False
    for name, (got, want) in errors.items():
        scale = expected * (a.theta if name == "rho theta gain" else 1.0)
        err = abs(got - want) / scale
        print("check_mass_source: t = %.6g s, %s %.17g, expected %.17g (rel err %.2e, tol %.1e)"
              % (times[-1], name, got, want, err, a.tol))
        ok = ok and err <= a.tol
    print("check_mass_source: " + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
