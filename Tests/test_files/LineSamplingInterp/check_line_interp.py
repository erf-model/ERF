#!/usr/bin/env python3
"""Check lines interpolated to requested heights (erf.sample_line_interp) against the plotfile.

Usage: check_line_interp.py PLOTFILE AGL PREFIX COVERED LINE:i:j:zlo:zhi:n [LINE:x:j:z ...]

A LINE:i:j:zlo:zhi:n spec is a vertical line in column (i,j) at n heights evenly spaced from zlo
to zhi; a LINE:x:j:z spec is a line in x through row j at height z, one sample per column.  With
AGL = 1 the heights are above the ground under each column.  LINE@L names a line sampled on
level L, whose column indices are those of level L.

For each line, PREFIX+LINE.z (last row) must be the target heights above sea level, and
PREFIX+LINE.cons_3 (rho times the scalar, rhoadv_0 in the plotfile) must be rhoadv_0 of the
column interpolated linearly in height between the two cell centres (z_phys) that bracket the
target: the lowest cell centre's value between the ground and it, the top cell centre's between
it and the domain top, and COVERED below the ground or above the top.  The interpolation here is
independent of ERF's.

The ground is estimated from the cell centres: the deck's levels start at 20 m and are
stretched by 1.2, so the two lowest centres are 10 s and 32 s above the ground (s is the
squeeze of the column), on every level (the refinement is in x and y only).
"""
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
for d in (os.path.join(_here, "..", "..", "..", "Exec", "Tools"), os.environ.get("ERF_TOOLS", "")):
    if d:
        sys.path.insert(0, d)
import erf_plotfile  # noqa: E402

# The text output has 6 significant digits
TOL = 2.0e-5


def column_ground(zc):
    """The ground under a column, from its two lowest cell centres (see the docstring)."""
    return zc[0] - (zc[1] - zc[0]) * 10.0 / 22.0


def interp(zc, f, z, top, covered):
    """f of a column at height z: linear between the cell centres around it, clamped."""
    if z < column_ground(zc) or z > top:
        return covered
    if z <= zc[0]:
        return f[0]
    if z >= zc[-1]:
        return f[-1]
    for k in range(len(zc) - 1):
        if zc[k] <= z < zc[k + 1]:
            w = (z - zc[k]) / (zc[k + 1] - zc[k])
            return (1.0 - w) * f[k] + w * f[k + 1]
    raise RuntimeError("no bracket for z = %g" % z)


def last_row(fname):
    with open(fname) as f:
        rows = [[float(x) for x in line.split()] for line in f if line.strip()]
    return rows[-1][1:]


def close(a, b, scale):
    return abs(a - b) <= TOL * max(scale, abs(b))


def main(argv):
    if len(argv) < 6:
        print(__doc__)
        return 2
    plotfile, agl, prefix, covered = argv[1], argv[2] == "1", argv[3], float(argv[4])
    levels = {}

    def level(lev):
        if lev not in levels:
            hdr, data = erf_plotfile.read_fields(plotfile, ["z_phys", "rhoadv_0"], level=lev)
            levels[lev] = (data["z_phys"], data["rhoadv_0"], hdr["prob_hi"][2])
        return levels[lev]

    failed = 0
    for spec in argv[5:]:
        fields = spec.split(":")
        name, _, lev = fields[0].partition("@")
        zp, rs, top = level(int(lev) if lev else 0)
        nx = len(zp)
        if fields[1] == "x":
            j, z = int(fields[2]), float(fields[3])
            samples = [(i, j, z) for i in range(nx)]
        else:
            i, j, zlo, zhi, n = int(fields[1]), int(fields[2]), float(fields[3]), float(fields[4]), int(fields[5])
            samples = [(i, j, zlo if n == 1 else zlo + (zhi - zlo) * m / (n - 1)) for m in range(n)]
        got_z = last_row(prefix + name + ".z")
        got_s = last_row(prefix + name + ".cons_3")
        if len(got_z) != len(samples) or len(got_s) != len(samples):
            print("%s: %d samples, %d heights and %d values written"
                  % (prefix + name, len(samples), len(got_z), len(got_s)))
            failed += 1
            continue
        bad, nbelow, worst = 0, 0, 0.0
        for (i, j, h), gz, gs in zip(samples, got_z, got_s):
            zc = zp[i][j]
            f = rs[i][j]
            zt = h + column_ground(zc) if agl else h
            if zc[0] == 0.0 and zc[1] == 0.0:
                print("  sample (%d,%d): the column is not on this level" % (i, j))
                bad += 1
                continue
            want = interp(zc, f, zt, top, covered)
            nbelow += want == covered
            # the heights to 6 digits and a ground estimated to TOL; values of order 0.1 to 1
            if not close(gz, zt, 1.0) or not close(gs, want, 1.0e-2):
                bad += 1
                print("  sample (%d,%d) at %.3f m: z %.6g (want %.6g), cons_3 %.6g (want %.6g)"
                      % (i, j, h, gz, zt, gs, want))
            worst = max(worst, abs(gs - want))
        print("%s (level %s): %d samples, %d missing, worst difference %.2e: %s"
              % (prefix + name, lev or "0", len(samples), nbelow, worst, "match" if bad == 0 else "DIFFER"))
        failed += bad > 0
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
