#!/usr/bin/env python3
"""Check the force-restore surface fields of a level 1 made over a shore.

The runs of TwoStream_SEBLandMaskNewLevel.i read the land mask of make_wrfinput.py: land in
the west, water in the east, one land column in the water. Level 1 has no wrfinput file, so
it takes its parent's mask, and it covers coarse columns i = 6..11, j = 4..9. Asserts:

1. --regrid-dir, level 1 made by the regrid at t = 5 s (ERF::MakeNewLevelFromCoarse), at every
   2D plotfile from --first-step to --last-step:
   a. level 1 covers exactly the fine columns over the coarse patch;
   b. a level 1 column over a coarse water column holds --skin-temp and --soil-water, and
      the same bits as a level 0 water column that level 1 does not cover;
   c. a level 1 column over a coarse land column has moved at least --min-change from
      --skin-temp and is within --land-rtol of a level 0 land column that level 1 does not
      cover (the deck is uniform over land);
   d. a level 0 water column holds --skin-temp and --soil-water, under level 1 or not;
   e. the log has no "land column(s) lie over a coarse water column" line: no level 1 land
      column lies over coarse water.
2. --scratch-dir, level 1 from the start (ERF::MakeNewLevelFromScratch): 1a, 1b, 1d and 1e at
   steps 0 and --last-step, and 1c at --last-step. At the start land and water hold the same
   values, so 1e is the check that sees the order of the fill and the mask there. The text
   of the line is pinned by RunTwoStreamSEBLandMask.cmake, which expects it.
3. --nestfile-dir, level 1 from the start with its own wrfinput file (wrfinput_d02 of
   make_wrfinput.py, the parent's mask under each cell): the checks of 2.
4. --restart-dir, the regrid run restarted from a checkpoint written after the regrid: its
   last 2D plotfile holds the same values as the uninterrupted run's, on both levels.
"""

import argparse
import os
import sys

NX = NY = 16
PATCH_I = (6, 11)
PATCH_J = (4, 9)
FAR_WATER = (15, 0)   # a level 0 water column level 1 does not cover
FAR_LAND = (0, 15)    # a level 0 land column level 1 does not cover
FIELDS = ('seb_t_sfc', 'seb_q_sfc')
FILL_LINE = 'land column(s) lie over a coarse water column'


def land(i, j):
    """LANDMASK of coarse column (i, j), as make_wrfinput.py writes it."""
    return i <= 7 or (i == 8 and j == 6)


def level_columns(reader, plotfile, level):
    """{(i, j): {field: value}} of the columns the level's grids cover."""
    hdr, data = reader.read_fields(plotfile, FIELDS, level=level)
    lo, _ = hdr["domains"][level]
    boxes = reader.level_boxes(plotfile, level)
    out = {}
    for blo, bhi in boxes:
        for j in range(blo[1], bhi[1] + 1):
            for i in range(blo[0], bhi[0] + 1):
                out[(i, j)] = {f: data[f][i - lo[0]][j - lo[1]][0] for f in FIELDS}
    return out


# What a missing or short plotfile or log raises in the reader: reported as a failure of
# the check, so that the failures collected so far are still printed.
READ_ERRORS = (OSError, IndexError, KeyError, ValueError, StopIteration, AttributeError)


def check_levels(reader, directory, step, args, deck, check_land, failures):
    tag = f"{os.path.basename(os.path.normpath(directory))} step {step}"
    plotfile = os.path.join(directory, f"{args.prefix}{step:0{args.digits}d}")
    try:
        _check_levels(reader, plotfile, tag, args, deck, check_land, failures)
    except READ_ERRORS as err:
        failures.append(f"{tag}: could not read {plotfile}: {type(err).__name__}: {err}")


def _check_levels(reader, plotfile, tag, args, deck, check_land, failures):
    if reader.finest_level(plotfile) != 1:
        failures.append(f"{tag}: the 2D plotfile has no level 1")
        return
    crse = level_columns(reader, plotfile, 0)
    fine = level_columns(reader, plotfile, 1)

    expected = {(i, j) for i in range(2 * PATCH_I[0], 2 * PATCH_I[1] + 2)
                for j in range(2 * PATCH_J[0], 2 * PATCH_J[1] + 2)}
    if set(fine) != expected:
        failures.append(f"{tag}: level 1 covers {len(fine)} columns, not the "
                        f"{len(expected)} over the patch")
        return

    n_water = n_land = 0
    for (i, j), column in sorted(fine.items()):
        if land(i // 2, j // 2):
            n_land += 1
            if not check_land:
                continue
            moved = abs(column['seb_t_sfc'] - deck['seb_t_sfc'])
            if moved < args.min_change:
                failures.append(f"{tag}: level 1 land column ({i},{j}) moved only "
                                f"{moved:.3e} K from {deck['seb_t_sfc']}")
            ref = crse[FAR_LAND]['seb_t_sfc']
            if abs(column['seb_t_sfc'] - ref) > args.land_rtol * abs(ref):
                failures.append(f"{tag}: level 1 land column ({i},{j}) holds seb_t_sfc = "
                                f"{column['seb_t_sfc']!r}, level 0 land holds {ref!r}")
            continue
        n_water += 1
        for name in FIELDS:
            value = column[name]
            if abs(value - deck[name]) > args.start_rtol * abs(deck[name]):
                failures.append(f"{tag}: level 1 water column ({i},{j}) holds {name} = "
                                f"{value!r}, not the held {deck[name]}")
            elif value != crse[FAR_WATER][name]:
                failures.append(f"{tag}: level 1 water column ({i},{j}) {name} = {value!r} "
                                f"is not the level 0 water value {crse[FAR_WATER][name]!r}")
    # 6 x 6 coarse columns, 13 of them land: columns i = 6, 7 and the column (8, 6).
    if (n_land, n_water) != (4 * 13, 4 * 23):
        failures.append(f"{tag}: {n_land} land and {n_water} water columns on level 1, "
                        f"expected 52 and 92")

    for (i, j), column in sorted(crse.items()):
        if land(i, j):
            continue
        for name in FIELDS:
            if abs(column[name] - deck[name]) > args.start_rtol * abs(deck[name]):
                failures.append(f"{tag}: level 0 water column ({i},{j}) holds {name} = "
                                f"{column[name]!r}, not the held {deck[name]}")
    if check_land:
        print(f"{tag}: level 1 land seb_t_sfc {fine[(12, 8)]['seb_t_sfc']!r}, "
              f"level 0 land {crse[FAR_LAND]['seb_t_sfc']!r}, "
              f"water {fine[(23, 8)]['seb_t_sfc']!r}")


def check_log(directory, failures):
    tag = os.path.basename(os.path.normpath(directory))
    try:
        with open(os.path.join(directory, 'simulation.log')) as handle:
            lines = [line.strip() for line in handle if FILL_LINE in line]
    except OSError as err:
        failures.append(f"{tag}: could not read the log: {err}")
        return
    for line in lines:
        failures.append(f"{tag}: the fill counted land columns over water: {line}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--reader-dir', required=True, help='directory of erf_plotfile.py')
    parser.add_argument('--regrid-dir', required=True)
    parser.add_argument('--scratch-dir', required=True)
    parser.add_argument('--nestfile-dir', required=True)
    parser.add_argument('--restart-dir', required=True)
    parser.add_argument('--first-step', type=int, required=True,
                        help='first 2D plotfile of the regrid run with level 1')
    parser.add_argument('--last-step', type=int, required=True)
    parser.add_argument('--prefix', default='plt2d')
    parser.add_argument('--digits', type=int, default=5)
    parser.add_argument('--skin-temp', type=float, default=299.0, help='erf.rad_t_sfc [K]')
    parser.add_argument('--soil-water', type=float, default=0.0095,
                        help='erf.radiation.seb_q_sfc_default')
    parser.add_argument('--min-change', type=float, default=1.0e-3, help='[K]')
    parser.add_argument('--land-rtol', type=float, default=1.0e-5,
                        help='level 1 land against level 0 land; one step moves land by 6e-5')
    parser.add_argument('--start-rtol', type=float, default=1.0e-6,
                        help='a held value against the deck value')
    args = parser.parse_args()

    sys.path.insert(0, args.reader_dir)
    import erf_plotfile as reader  # noqa: E402

    deck = {'seb_t_sfc': args.skin_temp, 'seb_q_sfc': args.soil_water}
    failures = []

    # The regrid run has one level until the regrid.
    before = os.path.join(args.regrid_dir, f"{args.prefix}{args.first_step - 1:0{args.digits}d}")
    try:
        if reader.finest_level(before) != 0:
            failures.append(f"regrid: level 1 exists before step {args.first_step}; the run "
                            f"does not make it by a regrid")
    except READ_ERRORS as err:
        failures.append(f"regrid: could not read {before}: {type(err).__name__}: {err}")
    for step in range(args.first_step, args.last_step + 1):
        check_levels(reader, args.regrid_dir, step, args, deck, True, failures)
    check_log(args.regrid_dir, failures)

    for directory in (args.scratch_dir, args.nestfile_dir):
        for step in (0, args.last_step):
            check_levels(reader, directory, step, args, deck, step > 0, failures)
        check_log(directory, failures)

    name = f"{args.prefix}{args.last_step:0{args.digits}d}"
    for level in (0, 1):
        try:
            ran = level_columns(reader, os.path.join(args.regrid_dir, name), level)
            again = level_columns(reader, os.path.join(args.restart_dir, name), level)
        except READ_ERRORS as err:
            failures.append(f"restart: could not read level {level} of {name}: "
                            f"{type(err).__name__}: {err}")
            continue
        if set(ran) != set(again):
            failures.append(f"restart: level {level} covers other columns than the "
                            f"uninterrupted run")
            continue
        differ = [key for key in sorted(ran) if ran[key] != again[key]]
        if differ:
            i, j = differ[0]
            failures.append(f"restart: {len(differ)} level {level} column(s) differ from the "
                            f"uninterrupted run, first ({i},{j}): {again[(i, j)]} "
                            f"against {ran[(i, j)]}")

    if failures:
        print("TwoStream_SEBLandMaskNewLevel: FAILED")
        for failure in failures:
            print("  " + failure)
        return 1
    print("TwoStream_SEBLandMaskNewLevel: passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
