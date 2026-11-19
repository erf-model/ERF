#!/usr/bin/env python3
"""Check the heights written with plane samples and subvolumes against the 3D plotfile.

Usage: check_plane_heights.py [--plt PLT] [--plane FILE:DIR:POS ...] [--flat-plane FILE:POS ...]
                              [--subvol FILE:I0:J0:K0 ...] [--names FILE:NAME[,NAME...] ...]
                              [--no-names FILE:NAME[,NAME...] ...] [--tol T]

--plane      a plane plotfile, the plane normal to DIR at POS (in the computational coordinates).
             Its z_phys and theta must be the 3D plotfile's (--plt) interpolated linearly in DIR
             between the two cell centers that bracket POS, the cell itself within half a cell of
             the domain edge, as the plane sampler does.
--flat-plane a horizontal plane plotfile at POS on a mesh with constant dz: its z_phys must be POS.
--subvol     a subvolume plotfile whose box starts at cell (I0,J0,K0) of the 3D plotfile: its
             z_phys and theta must be the 3D plotfile's in that box.
--names      the variables FILE must hold; --no-names those it must not.

A sample fails when it differs from the expected value by more than --tol times the largest
expected value. Exit status 0 on pass, 1 on a mismatch.
"""

import argparse
import math
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, "..", "Exec", "Tools"))
import erf_plotfile  # noqa: E402


def names_of(plotfile):
    return erf_plotfile._read_header(plotfile)["names"]


def domain_lo(hdr):
    return hdr["domains"][0][0]


def compare(label, got, want, tol):
    """got and want are flat lists of the same length; report and return whether they agree."""
    scale = max(abs(v) for v in want) or 1.0
    worst = max(abs(g - w) for g, w in zip(got, want))
    bad = any(not math.isfinite(g) for g in got) or worst > tol * scale
    print("%s %s: %d samples, max |diff| %.3e (scale %.3e)" %
          ("FAIL" if bad else "ok  ", label, len(got), worst, scale))
    return not bad


def check_plane(spec, plt, tol):
    fname, d, pos = spec.split(":")
    d = int(d)
    pos = float(pos)
    hdr3, f3 = plt
    lo3 = hdr3["prob_lo"][d]
    dx = hdr3["dx"][d]
    n = hdr3["domains"][0][1][d] - hdr3["domains"][0][0][d] + 1

    # The two cells whose centers bracket pos, and the weight of the upper one; within half a
    # cell of the domain edge the cell pos lies in, alone
    s = (pos - lo3) / dx - 0.5
    il = int(math.floor(s))
    w = s - il
    if il < 0:
        il, w = 0, 0.0
    if il >= n - 1:
        il, w = n - 1, 0.0
    ih = min(il + 1, n - 1)

    hdr, fp = erf_plotfile.read_fields(fname, ["z_phys", "theta"])
    plo = domain_lo(hdr)
    ok = True
    for var in ("z_phys", "theta"):
        a = fp[var]
        got, want = [], []
        for i in range(len(a)):
            for j in range(len(a[0])):
                for k in range(len(a[0][0])):
                    idx = [i + plo[0], j + plo[1], k + plo[2]]
                    idx[d] = il
                    lo_v = f3[var][idx[0]][idx[1]][idx[2]]
                    idx[d] = ih
                    hi_v = f3[var][idx[0]][idx[1]][idx[2]]
                    got.append(a[i][j][k])
                    want.append((1.0 - w) * lo_v + w * hi_v)
        ok &= compare("%s %s" % (fname, var), got, want, tol)
    return ok


def check_flat_plane(spec, tol):
    fname, pos = spec.split(":")
    _, fp = erf_plotfile.read_fields(fname, ["z_phys"])
    got = [v for plane in fp["z_phys"] for row in plane for v in row]
    return compare("%s z_phys" % fname, got, [float(pos)] * len(got), tol)


def check_subvol(spec, plt, tol):
    fname, i0, j0, k0 = spec.split(":")
    i0, j0, k0 = int(i0), int(j0), int(k0)
    _, f3 = plt
    ok = True
    for var in ("z_phys", "theta"):
        _, fs = erf_plotfile.read_fields(fname, [var])
        boxes = erf_plotfile.level_boxes(fname, 0)
        got, want = [], []
        for blo, bhi in boxes:
            for i in range(blo[0], bhi[0] + 1):
                for j in range(blo[1], bhi[1] + 1):
                    for k in range(blo[2], bhi[2] + 1):
                        got.append(fs[var][i][j][k])
                        want.append(f3[var][i][j][k])
        lo = [min(b[0][d] for b in boxes) for d in range(3)]
        if lo != [i0, j0, k0]:
            print("FAIL %s: the box starts at %s, not (%d,%d,%d)" % (fname, lo, i0, j0, k0))
            return False
        ok &= compare("%s %s" % (fname, var), got, want, tol)
    return ok


def check_names(spec, present):
    fname, wanted = spec.split(":")
    have = names_of(fname)
    ok = True
    for name in wanted.split(","):
        if (name in have) != present:
            print("FAIL %s %s %s (has: %s)" %
                  (fname, "lacks" if present else "holds", name, " ".join(have)))
            ok = False
    if ok:
        print("ok   %s %s %s" % (fname, "holds" if present else "lacks", wanted))
    return ok


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--plt")
    ap.add_argument("--plane", action="append", default=[])
    ap.add_argument("--flat-plane", action="append", default=[])
    ap.add_argument("--subvol", action="append", default=[])
    ap.add_argument("--names", action="append", default=[])
    ap.add_argument("--no-names", action="append", default=[])
    ap.add_argument("--tol", type=float, default=1e-12)
    a = ap.parse_args(argv)

    plt = None
    if a.plane or a.subvol:
        if not a.plt:
            ap.error("--plane and --subvol need --plt")
        plt = erf_plotfile.read_fields(a.plt, ["z_phys", "theta"])

    ok = True
    for spec in a.names:
        ok &= check_names(spec, True)
    for spec in a.no_names:
        ok &= check_names(spec, False)
    for spec in a.plane:
        ok &= check_plane(spec, plt, a.tol)
    for spec in a.flat_plane:
        ok &= check_flat_plane(spec, a.tol)
    for spec in a.subvol:
        ok &= check_subvol(spec, plt, a.tol)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
