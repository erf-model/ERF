#!/usr/bin/env python3
"""Check that the level 1 grids cover exactly the cells a refinement indicator tags.

Usage: check_tagging_conditions.py <fields plotfile> <grids plotfile>
           --cond FIELD OP VALUE [--cond ...] [--combine and|or]
           [--box XLO XHI YLO YHI [ZLO ZHI]] [--ref-ratio R]

The expected tags are recomputed here from the level 0 fields of the first plotfile: a cell
passes a condition when FIELD >= VALUE (OP ge, value_greater) or FIELD <= VALUE (OP le,
value_less), and the criterion when it passes all the conditions (combine and) or any one of
them (or). With --box it must also have its center inside the box, bounds included, with z
the height of the center above the flat bottom of the domain. The tagged run is the second
plotfile, written with amr.n_error_buf = 0, a coarse blocking factor of 1 and
amr.grid_eff = 1, so its level 1 grids, coarsened by R, are the tagged cells.

The fields plotfile is from a run of the same deck with amr.max_level = 0: the level 0 data
of the tagged run has been averaged down from level 1, so it is no longer the data that was
tagged. A condition that is always true (e.g. density ge 0) leaves only the box, and then
both plotfiles may be the tagged one.

The check stops as vacuous when nothing is tagged, when the box (if given) excludes no cell
that passes the conditions, or, with several conditions, when the combination tags the same
cells inside the box as one of its conditions alone.
Exit status 0 on pass, 1 on a mismatch, 2 on a vacuous check.
"""

import argparse
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, "..", "Exec", "Tools"))
import erf_plotfile  # noqa: E402


def covered_cells(plotfile, ratio):
    """Level 0 cells covered by the level 1 grids of plotfile."""
    level_dir = os.path.join(plotfile, "Level_1")
    if not os.path.isdir(level_dir):
        sys.exit("FAIL: %s has no level 1" % plotfile)
    _, _, boxes, _ = erf_plotfile._read_cell_h(level_dir)
    cells = set()
    for lo, hi in boxes:
        for d in range(3):
            if lo[d] % ratio != 0 or (hi[d] + 1) % ratio != 0:
                sys.exit("FAIL: level 1 box %s %s does not cover whole level 0 cells" % (lo, hi))
        for i in range(lo[0] // ratio, (hi[0] + 1) // ratio):
            for j in range(lo[1] // ratio, (hi[1] + 1) // ratio):
                for k in range(lo[2] // ratio, (hi[2] + 1) // ratio):
                    cells.add((i, j, k))
    return cells


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("plt_fields")
    ap.add_argument("plt_grids")
    ap.add_argument("--cond", nargs=3, action="append", required=True,
                    metavar=("FIELD", "OP", "VALUE"))
    ap.add_argument("--combine", choices=["and", "or"], default="and")
    ap.add_argument("--box", nargs="+", type=float, default=None)
    ap.add_argument("--ref-ratio", type=int, default=2)
    a = ap.parse_args(argv)

    if a.box is not None and len(a.box) not in (4, 6):
        sys.exit("--box takes XLO XHI YLO YHI [ZLO ZHI]")
    for _, op, _ in a.cond:
        if op not in ("ge", "le"):
            sys.exit("condition operator must be ge or le, not %s" % op)

    fields = sorted({c[0] for c in a.cond})
    hdr, data = erf_plotfile.read_fields(a.plt_fields, fields)
    lo, hi, dx, prob_lo = hdr["lo"], hdr["hi"], hdr["dx"], hdr["prob_lo"]
    n = [hi[d] - lo[d] + 1 for d in range(3)]

    def passes(cond, i, j, k):
        field, op, value = cond
        v = data[field][i][j][k]
        return v >= float(value) if op == "ge" else v <= float(value)

    def in_box(i, j, k):
        if a.box is None:
            return True
        c = [prob_lo[d] + (idx + 0.5) * dx[d] for d, idx in enumerate((i, j, k))]
        c[2] -= prob_lo[2]
        b = a.box
        nd = len(b) // 2
        return all(b[2 * d] <= c[d] <= b[2 * d + 1] for d in range(nd))

    combine = all if a.combine == "and" else any
    expected, outside, single = set(), 0, [set() for _ in a.cond]
    for i in range(n[0]):
        for j in range(n[1]):
            for k in range(n[2]):
                inside = in_box(i, j, k)
                for s, cond in zip(single, a.cond):
                    if inside and passes(cond, i, j, k):
                        s.add((i, j, k))
                if combine(passes(c, i, j, k) for c in a.cond):
                    if inside:
                        expected.add((i, j, k))
                    else:
                        outside += 1

    print("expected tags: %d; passing the conditions outside the box: %d; per condition in box: %s"
          % (len(expected), outside, [len(s) for s in single]))
    if not expected:
        print("VACUOUS: no cell is tagged")
        return 2
    if a.box is not None and outside == 0:
        print("VACUOUS: the box excludes no cell that passes the conditions")
        return 2
    if len(a.cond) > 1 and any(s == expected for s in single):
        print("VACUOUS: the combination tags the same cells as one condition alone")
        return 2

    covered = covered_cells(a.plt_grids, a.ref_ratio)
    missing = sorted(expected - covered)
    extra = sorted(covered - expected)
    print("level 1 covers %d level 0 cells" % len(covered))
    if missing or extra:
        print("FAIL: %d expected cells not refined (first %s), %d refined cells not expected (first %s)"
              % (len(missing), missing[:5], len(extra), extra[:5]))
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
