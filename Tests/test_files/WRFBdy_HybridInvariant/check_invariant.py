#!/usr/bin/env python3
"""Check the wrfbdy invariant that erf.check_wrfbdy_invariant prints.

For each side (xlo, xhi, ylo, yhi) and field (u, v, theta, qv, rho), ERF prints the largest
relative difference between the decoupled, interpolated boundary data at the first boundary
time and the initial state:

  wrfbdy invariant xlo u: max |bdy - init| = ... (rel ...)

  check_invariant.py LOG            every line present, every value within its tolerance
  check_invariant.py LOG --fail F   a control for a known-wrong input: field F must be out
                                    of tolerance on at least one side; the fields not named
                                    are reported but not judged
  check_invariant.py LOG --log TEXT TEXT must also appear in the log

The tolerances: the winds agree to the float precision of the files; theta, qv and the
density agree less closely because the initial state is remapped as rho theta, rho qv and
rho, and the boundary data as theta, qv and rho.
"""

import argparse
import re
import sys

SIDES = ("xlo", "xhi", "ylo", "yhi")
TOLERANCES = {"u": 1.0e-5, "v": 1.0e-5, "theta": 1.0e-4, "qv": 5.0e-4, "rho": 1.0e-4}
LINE = re.compile(r"wrfbdy invariant (\w+) (\w+): max \|bdy - init\| = (\S+) \(rel (\S+)\)")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("log")
    parser.add_argument("--fail", action="append", default=[], choices=sorted(TOLERANCES))
    parser.add_argument("--log", dest="texts", action="append", default=[])
    args = parser.parse_args()

    with open(args.log) as f:
        text = f.read()
    rel = {}
    for side, field, _, value in LINE.findall(text):
        rel[(side, field)] = float(value)

    ok = True
    for field, tol in TOLERANCES.items():
        values = []
        for side in SIDES:
            if (side, field) not in rel:
                print(f"FAIL: no invariant line for {side} {field}")
                ok = False
                continue
            values.append(rel[(side, field)])
        if not values:
            continue
        worst = max(values)
        if field in args.fail:
            status = "ok (out of tolerance, as expected)" if worst > tol else "FAIL (expected out of tolerance)"
            ok = ok and worst > tol
        elif args.fail:
            status = "(not judged)"
        else:
            status = "ok" if worst <= tol else "FAIL"
            ok = ok and worst <= tol
        print(f"{field:6s} max rel {worst:.3e}  tol {tol:.1e}  {status}")

    for t in args.texts:
        if t not in text:
            print(f"FAIL: '{t}' not in {args.log}")
            ok = False
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
