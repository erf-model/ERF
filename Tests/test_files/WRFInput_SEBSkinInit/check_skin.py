#!/usr/bin/env python3
"""Compare the balance's skin or deep temperature in a 2D plotfile with make_wrfinput.py.

Without --field, or with --field skin, the field is the skin temperature seb_t_sfc: the
file's TSK on land and erf.rad_t_sfc on water. With --field deep it is the deep temperature
seb_t_deep: the file's TMN on land and erf.radiation.seb_t_deep_default on water. The modes
below are written for the skin and read the same way for the deep temperature.

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
                                           With --averaged the level 0 columns under the
                                           nest hold the mean of the nest's land columns
                                           over them, as after a step with two-way coupling
  check_skin.py same PLT_A PLT_B           the two plotfiles hold the same skin, bit for bit,
                                           on every level
  check_skin.py near PLT --tol K           every land column of level 0 is within K of the
                                           file's value, and at least one differs from it
  check_skin.py log LOG TEXT [--absent]    LOG holds TEXT (or does not, with --absent)

The plotfile modes also take --log LOG --count N: the log must report N land columns set
(by ERF::init_seb_skin_from_file, or with --field deep by ERF::init_seb_t_deep_from_file).
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

import make_wrfinput as mw  # noqa: E402
from make_wrfinput import (LAND_COLUMN, NEST_LO, NEST_N, NEST_RATIO, NX, NY,  # noqa: E402
                           is_water, nest_is_water)

# Per field: the plotfile variable, the value a column holds when the file is not used, the
# formulas of the two files, and the line the log reports the land columns in.
FIELDS = {
    "skin": ("seb_t_sfc", mw.RAD_T_SFC, mw.tsk, mw.nest_tsk,
             re.compile(r"ERF::init_seb_skin_from_file: level (\d+): the skin temperature of (\d+) "
                        r"land column\(s\) starts from the initial file")),
    "deep": ("seb_t_deep", mw.T_DEEP_DEFAULT, mw.tmn, mw.nest_tmn,
             re.compile(r"ERF::init_seb_t_deep_from_file: level (\d+): the deep temperature of (\d+) "
                        r"land column\(s\) is the TMN of the initial file")),
}
FIELD, DEFAULT, value, nest_value, REPORT = FIELDS["skin"]


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
    parser.add_argument("mode", choices=("file", "default", "tagged", "nested", "same", "near", "log"))
    parser.add_argument("args", nargs="+")
    parser.add_argument("--all-land", action="store_true")
    parser.add_argument("--absent", action="store_true")
    parser.add_argument("--land", type=float, default=None)
    parser.add_argument("--averaged", action="store_true")
    parser.add_argument("--tol", type=float, default=None)
    parser.add_argument("--log")
    parser.add_argument("--count", type=int, action="append", default=[])
    parser.add_argument("--field", choices=sorted(FIELDS), default="skin")
    a = parser.parse_args()
    global FIELD, DEFAULT, value, nest_value, REPORT
    FIELD, DEFAULT, value, nest_value, REPORT = FIELDS[a.field]

    errors = []
    if a.mode == "log":
        log_path, text = a.args
        found = text in open(log_path).read()
        if found == a.absent:
            errors.append(f"{log_path}: '{text}' " + ("is present" if found else "is not present"))
    elif a.mode == "same":
        finest = erf_plotfile.finest_level(a.args[0])
        if finest != erf_plotfile.finest_level(a.args[1]):
            errors.append("the two plotfiles do not have the same levels")
            finest = -1
        for level in range(finest + 1):
            la, ha, da = skin(a.args[0], level)
            lb, hb, db = skin(a.args[1], level)
            if (la, ha) != (lb, hb):
                errors.append(f"level {level}: the two plotfiles are not on the same domain")
            elif da != db:
                errors.append(f"level {level}: {a.args[0]} and {a.args[1]} hold different {FIELD}")
            elif level == 0 and all(v[0] == DEFAULT for row in da for v in row):
                errors.append(f"both plotfiles hold the default {DEFAULT} K everywhere: the comparison shows nothing")
        print(f"{finest + 1} level(s) compared")
    elif a.mode == "near":
        lo, hi, got = skin(a.args[0])
        differs = 0
        for i in range(lo[0], hi[0] + 1):
            for j in range(lo[1], hi[1] + 1):
                if is_water(i, j):
                    continue
                have = got[i - lo[0]][j - lo[1]][0]
                if not abs(have - value(i, j)) <= a.tol:
                    errors.append(f"column ({i}, {j}): {have!r} is not within {a.tol} K of {value(i, j)!r}")
                differs += have != value(i, j)
        if differs == 0:
            errors.append(f"every land column still holds the file's value: {FIELD} did not evolve")
    else:
        plotfile = a.args[0]
        if a.mode == "default":
            n, errors = compare(plotfile, lambda i, j: DEFAULT)
        else:
            all_land = a.all_land or a.mode == "tagged"

            def expected(i, j):
                if is_water(i, j, all_land):
                    return DEFAULT
                if a.averaged and all(NEST_LO <= c < NEST_LO + NEST_N // NEST_RATIO for c in (i, j)):
                    # Two-way coupling: the mean of the nest's land columns over this one
                    fine = [nest_value(I, J)
                            for I in range(NEST_RATIO * i, NEST_RATIO * (i + 1))
                            for J in range(NEST_RATIO * j, NEST_RATIO * (j + 1))
                            if not nest_is_water(I, J)]
                    return sum(fine) / len(fine)
                if a.land is not None and (i, j) == LAND_COLUMN:
                    return a.land
                return value(i, j)

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
                n1, e1 = compare(plotfile, lambda i, j: value(i // NEST_RATIO, j // NEST_RATIO)
                                 if nest_is_water(i, j) else nest_value(i, j), level=1, boxes=boxes)
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
                n1, e1 = compare(plotfile, lambda i, j: value((i + 0.5) / 2 - 0.5, (j + 0.5) / 2 - 0.5),
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
