#!/usr/bin/env python3
"""Check the mass budget of a run with a constant-rate mass source (erf.mass_source).

Usage: check_mass_source.py <log> --rate RATE [--tol TOL]

<log> is the output of a run in a closed domain (nothing enters or leaves) with
erf.mass_source.rate = RATE, the tracer on, erf.sum_interval = 1 and erf.v >= 1, so that
ERF prints at every step

    TIME= t
     MASS       = M            (one level), or  MASS       SL/ML = M_sl M  (several)
     RHO SCALAR = S            (one level), or  RHO SCALAR SL/ML = S_sl S  (several)

The multilevel (last) value is used when both are printed. With the initial scalar zero,
the mass must grow by RATE * t and the scalar must equal RATE * t at the last time printed,
each to a relative --tol of RATE * t. A source added at the wrong stage, normalized by
the wrong volume (without the terrain metric, say) or lost by the substeps changes one
or both. The check also requires the source to have been on for more than one step.
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
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--tol", type=float, default=1e-10)
    a = ap.parse_args()

    times, mass, scal = [], [], []
    with open(a.log) as f:
        for line in f:
            if line.startswith("TIME="):
                times.append(float(line.split("=")[1]))
            elif re.match(r"^ *MASS +(SL/ML +)?=", line):
                mass.append(last_value(line))
            elif re.match(r"^ *RHO SCALAR +(SL/ML +)?=", line):
                scal.append(last_value(line))

    if not (len(times) == len(mass) == len(scal)) or len(times) < 3:
        print("check_mass_source: found %d TIME, %d MASS and %d RHO SCALAR lines in %s; need the "
              "same number, at least 3 (erf.sum_interval = 1, erf.v >= 1, transport_scalar on?)"
              % (len(times), len(mass), len(scal), a.log))
        return 2
    if times[0] != 0.0:
        print("check_mass_source: the first integrals are at t = %g, not 0" % times[0])
        return 2

    expected = a.rate * (times[-1] - times[0])
    dm = mass[-1] - mass[0]
    err_m = abs(dm - expected) / expected
    err_s = abs(scal[-1] - expected) / expected
    print("check_mass_source: t = %.6g s, expected %.17g kg; mass gain %.17g (rel err %.2e), "
          "tracer %.17g (rel err %.2e), tol %.1e"
          % (times[-1], expected, dm, err_m, scal[-1], err_s, a.tol))
    if scal[0] != 0.0:
        print("check_mass_source: FAIL: the tracer integral is %g at t = 0, not 0" % scal[0])
        return 1
    if err_m > a.tol or err_s > a.tol:
        print("check_mass_source: FAIL")
        return 1
    print("check_mass_source: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
