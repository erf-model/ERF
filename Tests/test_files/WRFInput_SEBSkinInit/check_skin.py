#!/usr/bin/env python3
"""Compare the balance's skin temperature in a 2D plotfile with make_wrfinput.py.

  check_skin.py file PLT [--all-land] [--land VALUE]
                                           land columns hold the file's TSK, water columns
                                           erf.rad_t_sfc; --land as in make_wrfinput.py
  check_skin.py default PLT                every column holds erf.rad_t_sfc
  check_skin.py tagged PLT                 an all-land file with a second level, refined by
                                           2, without a file: level 0 holds TSK and level 1
                                           the same linear field at its own column centers
  check_skin.py nested PLT                 the deck's file with a wrfinput_d02 for level 1:
                                           level 0 as in "file"; a level 1 land column
                                           holds the TSK of wrfinput_d02, a level 1 water
                                           column what level 0 gave it (here the TSK of
                                           the land column of level 0 under it)
  check_skin.py same PLT_A PLT_B           the two plotfiles hold the same skin, bit for bit
  check_skin.py log LOG TEXT [--absent]    LOG holds TEXT (or does not, with --absent)

The plotfile modes also take --log LOG --count N: the log must report N land columns set.
"""

import argparse
import os
import re
import sys

_here = os.path.dirname(os.path.abspath(__file__))
for d in (_here, os.path.join(_here, "..", "..", "..", "Exec", "Tools"), os.environ.get("ERF_TOOLS", "")):
    if d:
        sys.path.insert(0, d)
import erf_plotfile  # noqa: E402

from make_wrfinput import (LAND_COLUMN, NEST_LO, NEST_N, NEST_RATIO, NX, NY, RAD_T_SFC,  # noqa: E402
                           is_water, nest_is_water, nest_tsk, tsk)

FIELD = "seb_t_sfc"
REPORT = re.compile(r"ERF::init_seb_skin_from_file: level (\d+): the skin temperature of (\d+) "
                    r"land column\(s\) starts from the initial file")


def skin(plotfile, level=0):
    hdr, data = erf_plotfile.read_fields(plotfile, [FIELD], level=level)
    lo, hi = hdr["domains"][level]
    return lo, hi, data[FIELD]


def compare(plotfile, expected, level=0, boxes=None):
    """expected(i, j) -> K for the columns of the level (all, or those in boxes)."""
    lo, hi, got = skin(plotfile, level)
    errors = []
    n = 0
    for i in range(lo[0], hi[0] + 1):
        for j in range(lo[1], hi[1] + 1):
            if boxes is not None and not any(b[0][0] <= i <= b[1][0] and b[0][1] <= j <= b[1][1]
                                             for b in boxes):
                continue
            n += 1
            want = expected(i, j)
            have = got[i - lo[0]][j - lo[1]][0]
            if have != want:
                errors.append(f"level {level} column ({i}, {j}): {have!r}, expected {want!r}")
    return n, errors


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("file", "default", "tagged", "nested", "same", "log"))
    parser.add_argument("args", nargs="+")
    parser.add_argument("--all-land", action="store_true")
    parser.add_argument("--absent", action="store_true")
    parser.add_argument("--land", type=float, default=None)
    parser.add_argument("--log")
    parser.add_argument("--count", type=int, action="append", default=[])
    a = parser.parse_args()

    errors = []
    if a.mode == "log":
        log_path, text = a.args
        found = text in open(log_path).read()
        if found == a.absent:
            errors.append(f"{log_path}: '{text}' " + ("is present" if found else "is not present"))
    elif a.mode == "same":
        la, ha, da = skin(a.args[0])
        lb, hb, db = skin(a.args[1])
        if (la, ha) != (lb, hb):
            errors.append("the two plotfiles are not on the same domain")
        elif da != db:
            errors.append(f"{a.args[0]} and {a.args[1]} hold different skin temperatures")
        elif all(v[0] == RAD_T_SFC for row in da for v in row):
            errors.append("both plotfiles hold erf.rad_t_sfc everywhere: the comparison shows nothing")
    else:
        plotfile = a.args[0]
        if a.mode == "default":
            n, errors = compare(plotfile, lambda i, j: RAD_T_SFC)
        else:
            all_land = a.all_land or a.mode == "tagged"

            def expected(i, j):
                if is_water(i, j, all_land):
                    return RAD_T_SFC
                if a.land is not None and (i, j) == LAND_COLUMN:
                    return a.land
                return tsk(i, j)

            n, errors = compare(plotfile, expected)
        if n != NX * NY:
            errors.append(f"level 0 has {n} columns, expected {NX * NY}")
        if a.mode == "nested":
            if erf_plotfile.finest_level(plotfile) != 1:
                errors.append("the plotfile has no level 1")
            else:
                boxes = erf_plotfile.level_boxes(plotfile, 1)
                # The two water columns of the nest lie over the level 0 land column (5, 2)
                # and keep the value the fill from level 0 gave them, without a slope.
                n1, e1 = compare(plotfile, lambda i, j: tsk(i // NEST_RATIO, j // NEST_RATIO)
                                 if nest_is_water(i, j) else nest_tsk(i, j), level=1, boxes=boxes)
                errors += e1
                lo = NEST_RATIO * NEST_LO
                covered = all(any(b[0][0] <= i <= b[1][0] and b[0][1] <= j <= b[1][1] for b in boxes)
                              for i in range(lo, lo + NEST_N) for j in range(lo, lo + NEST_N))
                if n1 != NEST_N * NEST_N or not covered:
                    errors.append(f"level 1 has {n1} columns, expected the {NEST_N * NEST_N} of wrfinput_d02")
        if a.mode == "tagged":
            if erf_plotfile.finest_level(plotfile) != 1:
                errors.append("the plotfile has no level 1")
            else:
                # Refined by 2: the fine column I has its center at the level 0 coordinate
                # (I + 0.5) / 2 - 0.5, and the interpolation with slopes is exact for a
                # linear field away from the domain faces. Bit for bit only because the
                # slopes 0.5 and 0.25 are powers of two: every term is exact.
                boxes = erf_plotfile.level_boxes(plotfile, 1)
                n1, e1 = compare(plotfile, lambda i, j: tsk((i + 0.5) / 2 - 0.5, (j + 0.5) / 2 - 0.5),
                                 level=1, boxes=boxes)
                errors += e1
                if n1 == 0:
                    errors.append("level 1 has no column")
                print(f"level 1: {n1} columns checked")

    if a.log:
        counts = [int(m.group(2)) for m in REPORT.finditer(open(a.log).read())]
        if counts != a.count:
            errors.append(f"{a.log}: land columns set per level {counts}, expected {a.count}")

    for e in errors[:20]:
        print("FAIL:", e)
    if errors:
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
