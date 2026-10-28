#!/usr/bin/env python3
"""Check that two initial perturbations add up when both are turned on.

Usage: check_perturbation_sum.py <base> <A> <B> <AB> --vars V [V ...] [--tol T]

<base> is the initial state without either perturbation, <A> and <B> with one of them
each, and <AB> with both. For every listed field the change from both must be the sum
of the two changes alone, cell by cell:

    AB - base = (A - base) + (B - base)

to a relative --tol (relative to the largest |AB - base|). The temperature and the qv
perturbations of the moist setups are applied as such sums, so a setup that drops one
of them when the other is on (it applied the second from the unperturbed state) fails.
The check also stops if either perturbation alone leaves the field unchanged, since the
sum could then not tell anything apart.
Exit status 0 on pass, 1 on a mismatch, 2 when a perturbation does not show.
"""

import argparse
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, "..", "Exec", "Tools"))
import erf_plotfile  # noqa: E402


def flat(field):
    return [v for plane in field for col in plane for v in col]


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("base")
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("ab")
    ap.add_argument("--vars", nargs="+", required=True)
    ap.add_argument("--tol", type=float, default=1e-10)
    args = ap.parse_args(argv)

    runs = [erf_plotfile.read_fields(p, args.vars)[1] for p in (args.base, args.a, args.b, args.ab)]
    status = 0
    for name in args.vars:
        f0, fa, fb, fab = (flat(r[name]) for r in runs)
        da = [a - z for a, z in zip(fa, f0)]
        db = [b - z for b, z in zip(fb, f0)]
        dab = [c - z for c, z in zip(fab, f0)]
        scale = max(abs(v) for v in dab)
        if max(abs(v) for v in da) == 0.0 or max(abs(v) for v in db) == 0.0 or scale == 0.0:
            print("%s: a perturbation leaves it unchanged; the check is vacuous" % name)
            return 2
        err = max(abs(c - (a + b)) for a, b, c in zip(da, db, dab))
        rel = err / scale
        ok = rel <= args.tol
        print("%-10s max |d(both) - d(A) - d(B)| = %.3e, relative %.3e (largest d(both) %.3e)  %s"
              % (name, err, rel, scale, "ok" if ok else "FAIL"))
        if not ok:
            status = 1
    print("PASS" if status == 0 else "FAIL")
    return status


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
