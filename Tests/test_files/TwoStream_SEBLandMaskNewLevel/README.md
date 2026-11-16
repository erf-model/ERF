# TwoStream_SEBLandMaskNewLevel

## Purpose

The two-stream force-restore surface energy balance holds its state on the water columns of
the land mask. When a new fine level is made, its surface fields are filled from the parent,
and the fill reads the fine level's land mask to choose between a land value and the held
value. With `wrfinput`, a level with no file of its own takes its parent's mask. This case
checks that the fill runs after the mask has reached the level: on a level 1 made over a
shore, every column over coarse water keeps the parent's held value, and the fill counts no
land column over water.

Before the change the fill ran first, against a mask that was land everywhere
(`erf.is_land`, default 1). The level 1 columns over coarse water within one coarse column
of the shore then took a mean of land columns, and the log counted the columns farther out
as land columns with no land source.

`TwoStream_SEBLandMask` covers one mask value per level. This case has both values on
level 0 and on level 1.

## Files

| File | Role |
|---|---|
| `TwoStream_SEBLandMaskNewLevel.i` | The deck: 16 x 16 x 8 cells of 1 km, `wrfinput` start, two-stream radiation with the prognostic balance and fixed turbulent fluxes, two-way coupling, level 1 over coarse columns i = 6..11, j = 4..9 tagged from t = 4.5 s. |
| `make_wrfinput.py` | Writes `wrfinput_d01` and `wrflowinp_d01`: flat, at rest, land for i <= 7, water for i >= 8, and one land column (8, 6) in the water. |
| `check_new_level.py` | The checker (below). |
| `../../RunTwoStreamSEBLandMaskNewLevel.cmake` | Writes the files, runs the three legs and the checker. |

## Data

None. The case uses no external data: `make_wrfinput.py` writes its two input files.

## How to run

With CTest, from a build directory configured with `-DERF_ENABLE_TESTS=ON`,
`-DERF_ENABLE_NETCDF=ON` and a Python 3 with `netCDF4` and `numpy`, after building
`erf_exec`:

```bash
ctest -R TwoStream_SEBLandMaskNewLevel --output-on-failure
```

By hand: run `make_wrfinput.py`, make three directories `regrid/`, `scratch/` and
`restart/`, copy `wrfinput_d01` and `wrflowinp_d01` into each, then:

```bash
cd regrid;  mpiexec -n 1 ../erf_exec ../TwoStream_SEBLandMaskNewLevel.i erf.check_int=6 > simulation.log; cd ..
cd scratch; mpiexec -n 1 ../erf_exec ../TwoStream_SEBLandMaskNewLevel.i erf.patch.start_time=0.0 > simulation.log; cd ..
cd restart; mpiexec -n 1 ../erf_exec ../TwoStream_SEBLandMaskNewLevel.i erf.restart=../regrid/chk00006 > simulation.log; cd ..
python3 check_new_level.py --reader-dir <ERF>/Exec/Tools --regrid-dir regrid \
    --scratch-dir scratch --restart-dir restart --first-step 6 --last-step 8
```

The checker reads `simulation.log` in each directory.

## Expected Results

Level 1 holds 144 columns: 52 over coarse land and 92 over coarse water.

| # | Leg | Check | Value |
|---|---|---|---|
| 1 | regrid | 2D plotfile 5 has one level, plotfiles 6 to 8 have two | level 1 is made by the regrid after step 5 |
| 2 | regrid, steps 6 to 8 | level 1 `seb_t_sfc`, `seb_q_sfc` over coarse water | 299 K and 0.0095, the bits of a level 0 water column |
| 3 | regrid, steps 6 to 8 | level 1 `seb_t_sfc` over coarse land | moved by more than 1e-3 K from 299 K; within a relative 1e-5 of a level 0 land column (one step moves land by 6e-5) |
| 4 | regrid, scratch | level 0 water columns | 299 K and 0.0095 |
| 5 | regrid, scratch | log line "land column(s) lie over a coarse water column" | absent |
| 6 | scratch, steps 0 and 8 | checks 2 and 4 with level 1 from the start; check 3 at step 8 | as above |
| 7 | restart | 2D plotfile 8 against the regrid leg's | the same bits on both levels |

With the fill run before the mask reaches the level, checks 2 and 5 fail on the regrid leg
and check 5 fails on the scratch leg: the log reports 60 level 1 land columns over water,
and the 32 level 1 columns over coarse water next to land hold the land value (299.086 K at
step 6). At the start land and water hold the same values, so on the scratch leg the log
line is the only check that sees the order. `TwoStream_SEBLandMask` expects the same line
with its count, so a change of its text does not go unseen.

The coupling is two-way, so the average-down of the surface fields runs at every step with
both mask values on both levels; check 4 covers the level 0 water columns under level 1.
The deep reservoir, the hydrology and the snow layer are off in this deck.

## CTest

`TwoStream_SEBLandMaskNewLevel`, labels `regression;radiation;wrfinput`, registered in
`Tests/CTestList.cmake` next to `WRFInput_CroppedParent`. It needs a NetCDF build and a
Python 3 with `netCDF4` and `numpy`.
