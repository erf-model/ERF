#!/usr/bin/env python3
"""Check that run B is run A seen from a frame moving at a uniform velocity.

Usage: check_frame_change.py <plotfile A> <plotfile B> --field NAME[=D] [--field ...] [--tol T]

Run B is run A with every velocity in the problem (the initial wind, a moving surface,
...) changed by the same uniform amount. Where the model depends only on velocity
differences, as a horizontally uniform flow without rotation does, run B is run A seen
from a moving frame: a velocity component of B is that of A plus D, the change along it,
and a field that does not depend on the frame (theta, a stress, ...) is the same in both.
So B(i,j,k) = A(i,j,k) + D for each field named with --field, D being 0 when omitted.

The fields are compared on level 0, over its whole domain. The check reports max |B - A - D|
for each field, relative to the larger of max |A| and |D|, and fails above --tol. A field
also fails when the two runs' domains differ in size, or when either run holds a value that
is not finite (Python's max would otherwise skip a NaN). The check stops when a field could
not tell the two runs apart: a field with D = 0 that is uniform everywhere in run A. A field
with a nonzero D can, uniform or not, since two runs that both ignored the change would agree.
Exit status 0 on pass, 1 on a mismatch, 2 on a vacuous field.
"""

import argparse
import math
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, "..", "Exec", "Tools"))
import erf_plotfile  # noqa: E402


def parse_field(text):
    name, sep, offset = text.partition("=")
    return name, (float(offset) if sep else 0.0)


def shape(f):
    return len(f), len(f[0]), len(f[0][0])


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("plt_a")
    ap.add_argument("plt_b")
    ap.add_argument("--field", type=parse_field, action="append", required=True,
                    help="NAME or NAME=D: compare B with A + D (D defaults to 0)")
    ap.add_argument("--tol", type=float, default=1e-12)
    a = ap.parse_args(argv)

    names = [name for name, _ in a.field]
    _, A = erf_plotfile.read_fields(a.plt_a, names)
    _, B = erf_plotfile.read_fields(a.plt_b, names)
    status = 0
    for name, off in a.field:
        fa, fb = A[name], B[name]
        if shape(fb) != shape(fa):
            print("%-12s domains differ: %s cells in run A, %s in run B  FAIL"
                  % (name, "x".join(map(str, shape(fa))), "x".join(map(str, shape(fb)))))
            status = 1
            continue
        nx, ny, nz = shape(fa)
        cells = [(i, j, k) for i in range(nx) for j in range(ny) for k in range(nz)]
        bad_a = sum(1 for i, j, k in cells if not math.isfinite(fa[i][j][k]))
        bad_b = sum(1 for i, j, k in cells if not math.isfinite(fb[i][j][k]))
        if bad_a or bad_b:
            print("%-12s values not finite: %d cells in run A, %d in run B  FAIL"
                  % (name, bad_a, bad_b))
            status = 1
            continue
        amax = max(abs(fa[i][j][k]) for i, j, k in cells)
        if off == 0.0:
            var = max(abs(fa[i][j][k] - fa[0][0][0]) for i, j, k in cells)
            if amax == 0.0 or var == 0.0:
                print("%s: uniform everywhere in run A with no offset; the check is vacuous" % name)
                return 2
        err = max(abs(fb[i][j][k] - fa[i][j][k] - off) for i, j, k in cells)
        rel = err / max(amax, abs(off))
        ok = rel <= a.tol
        print("%-12s max |B - A - (%g)| = %.3e, relative %.3e  %s"
              % (name, off, err, rel, "ok" if ok else "FAIL"))
        if not ok:
            status = 1
    print("PASS" if status == 0 else "FAIL")
    return status


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
