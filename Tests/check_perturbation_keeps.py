#!/usr/bin/env python3
"""Check what an initial perturbation changes and what it keeps.

Usage: check_perturbation_keeps.py <plotfile A> <plotfile B> [--same V [V ...]]
                                   --differ V [V ...] [--tol T]

A is the initial state without the perturbation and B with it. Every --same field
must agree cell by cell to a relative --tol (relative to its largest |value| in A),
and every --differ field must change somewhere by more than that, so the check
cannot pass because the perturbation was never applied. Exit status 0 on pass,
1 when a --same field changes, 2 when a --differ field does not.
"""

import argparse
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, "..", "Exec", "Tools"))
import erf_plotfile  # noqa: E402


def max_diff(fa, fb):
    amax = 0.0
    err = 0.0
    for pa, pb in zip(fa, fb):
        for ca, cb in zip(pa, pb):
            for va, vb in zip(ca, cb):
                amax = max(amax, abs(va))
                err = max(err, abs(vb - va))
    return err, amax


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("plt_a")
    ap.add_argument("plt_b")
    ap.add_argument("--same", nargs="*", default=[])
    ap.add_argument("--differ", nargs="+", required=True)
    ap.add_argument("--tol", type=float, default=1e-12)
    a = ap.parse_args(argv)

    names = a.same + a.differ
    _, A = erf_plotfile.read_fields(a.plt_a, names)
    _, B = erf_plotfile.read_fields(a.plt_b, names)
    status = 0
    for name in a.differ:
        err, amax = max_diff(A[name], B[name])
        rel = err / amax if amax > 0 else err
        ok = rel > a.tol
        print("%-10s changes by %.3e (relative %.3e)  %s" % (name, err, rel, "ok" if ok else "NOT PERTURBED"))
        if not ok:
            return 2
    for name in a.same:
        err, amax = max_diff(A[name], B[name])
        rel = err / amax if amax > 0 else err
        ok = rel <= a.tol
        print("%-10s changes by %.3e (relative %.3e)  %s" % (name, err, rel, "ok" if ok else "FAIL"))
        if not ok:
            status = 1
    print("PASS" if status == 0 else "FAIL")
    return status


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
