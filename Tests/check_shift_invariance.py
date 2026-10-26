#!/usr/bin/env python3
"""Check that moving geometry.prob_lo in x moves the initial state with it.

Usage: check_shift_invariance.py <plotfile A> <plotfile B> --shift S --vars V [V ...] [--tol T]

Run B is run A with the domain moved by a distance that maps the problem onto itself
when the cell index is shifted by S in x (a half wavelength of a sine wave whose sign
flips with it, a whole period of a periodic domain, ...). Then B(i,j,k) = A(i+S,j,k)
for every listed field, periodic in i. A field whose position is computed from the
cell index without geometry.prob_lo, or from prob_lo in a way that is not a shift,
breaks this.

The check reports max |B(i) - A(i+S)| over the domain for each field, relative to
max |A|, and fails above --tol. It stops when a field is uniform in x in run A: such
a field could not tell a shifted state from an unshifted one.
Exit status 0 on pass, 1 on a mismatch, 2 on a vacuous field.
"""

import argparse
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, "..", "Exec", "CanonicalTests", "Canonical_RANS"))
import erf_plotfile  # noqa: E402


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("plt_a")
    ap.add_argument("plt_b")
    ap.add_argument("--shift", type=int, required=True)
    ap.add_argument("--vars", nargs="+", required=True)
    ap.add_argument("--tol", type=float, default=1e-12)
    a = ap.parse_args(argv)

    _, A = erf_plotfile.read_fields(a.plt_a, a.vars)
    _, B = erf_plotfile.read_fields(a.plt_b, a.vars)
    status = 0
    for name in a.vars:
        fa, fb = A[name], B[name]
        nx, ny, nz = len(fa), len(fa[0]), len(fa[0][0])
        amax = max(abs(fa[i][j][k]) for i in range(nx) for j in range(ny) for k in range(nz))
        xvar = max(abs(fa[i][j][k] - fa[0][j][k]) for i in range(nx) for j in range(ny) for k in range(nz))
        if amax == 0.0 or xvar == 0.0:
            print("%s: uniform in x in run A; the check is vacuous" % name)
            return 2
        err = max(abs(fb[i][j][k] - fa[(i + a.shift) % nx][j][k])
                  for i in range(nx) for j in range(ny) for k in range(nz))
        rel = err / amax
        ok = rel <= a.tol
        print("%-12s max |B(i) - A(i%+d)| = %.3e, relative %.3e  %s"
              % (name, a.shift, err, rel, "ok" if ok else "FAIL"))
        if not ok:
            status = 1
    print("PASS" if status == 0 else "FAIL")
    return status


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
