# WRFInput_SEBSkinInit

## Purpose

The two-stream force-restore surface energy balance evolves a skin temperature `T_s` in
every land column. It started `T_s` at the one value `erf.rad_t_sfc` in every column. A
wrfinput file carries a skin temperature per column, `TSK`, and ERF already reads it for
the surface layer. With `erf.radiation.seb_init_skin_from_file = true` a land column of the
balance now starts at the file's `TSK`. A water column keeps `erf.rad_t_sfc`: the balance
does not run there. The default is `false`, and then every column starts at `erf.rad_t_sfc`
as before.

A finer level that reads a wrfinput file of its own takes that file's `TSK` on land in the
same way. A finer level without a file is filled from its parent, as before.

This case checks the first value of `T_s`, column by column, against the file.

The balance restores `T_s` toward a deep temperature `T_deep`, which was the one value
`erf.radiation.seb_t_deep_default` in every column. A wrfinput file carries `TMN`, the soil
temperature at the lower boundary, per column. With
`erf.radiation.seb_init_t_deep_from_file = true` the `T_deep` of a land column is the file's
`TMN`, and a water column keeps `seb_t_deep_default`. `T_deep` is then kept per column: ERF
no longer resets it before each step, a level without a file takes it from its parent, and
the checkpoint holds it. The default is `false`, and then every column holds
`seb_t_deep_default` as before. The `WRFInput_SEBDeepInit` tests check `T_deep` column by
column against the file, at the start, after two steps and after a restart.

## Files

| File | Role |
| --- | --- |
| `WRFInput_SEBSkinInit.i` | One level, 12 x 8 x 8 cells in two boxes, `max_step = 0`, two-stream radiation with `seb_prognostic_enable`, `erf.rad_t_sfc = 300`, `seb_t_deep_default = 290`, the skin input on and the deep input off (the `WRFInput_SEBDeepInit` tests set it on the command line); writes `seb_t_sfc` and `seb_t_deep` to the 2D plotfile `plt2d00000` |
| `make_wrfinput.py` | Writes `wrfinput_d01` and `wrflowinp_d01`; `--all-land` leaves out the water; `--land VALUE` puts the value in the land column (8, 5) and `--land-west VALUE` in the land column (2, 5) of the other box; `--water VALUE` puts it in the water column (8, 1); `--nest` also writes `wrfinput_d02` and `wrflowinp_d02` for a second level; for `TMN`, `--no-tmn` leaves it out, `--tmn-land VALUE`, `--tmn-land-west VALUE` and `--tmn-water VALUE` put the value in the same three columns |
| `check_skin.py` | Compares `seb_t_sfc` (or, with `--field deep`, `seb_t_deep`) of every column of the plotfile with the formulas in `make_wrfinput.py`; also compares two plotfiles and looks for a message in a log |

## Data and acknowledgements

None is downloaded and none is redistributed. `make_wrfinput.py` writes the files:

- `TSK = 280 + 0.5 i + 0.25 j` [K], with `i`, `j` the column indices from zero, so 280 K to
  287.25 K. Every value is exact in single precision.
- `LANDMASK = 0` (water) where `i >= 6` and `j <= 2`: 18 of the 96 columns, in both boxes
  of the deck (the boxes meet between `i = 7` and `i = 8`). The other 78 are land.
- `SST = 285` [K].
- `TMN = 275 + 0.25 i + 0.5 j` [K], so 275 K to 281.25 K, exact in single precision.
  `erf.radiation.seb_t_deep_default = 290` K is a value `TMN` never takes.
- With `--nest`, `wrfinput_d02` covers the level 0 columns 2..5 x 2..5 with 8 x 8 columns
  refined by 2. Its `TSK = 290 + 0.125 I + 0.0625 J` [K], with `I`, `J` the fine column
  indices of the whole domain (4..11), so it differs from the level 0 field everywhere. Its
  `LANDMASK = 0` in the fine columns (10, 5) and (11, 5): 62 of the 64 are land. Its
  `TMN = 284 + 0.0625 I + 0.125 J` [K].

`erf.rad_t_sfc = 300` K is a value `TSK` never takes, so a column shows which of the two it
holds.

## How to run

Needs ERF built with NetCDF, and Python with `netCDF4` and `numpy`. `check_skin.py` reads
the plotfile with `Exec/Tools/erf_plotfile.py`.

```bash
python3 make_wrfinput.py
./erf_exec WRFInput_SEBSkinInit.i > run.log
python3 check_skin.py file plt2d00000 --log run.log --count 78
```

For the deep temperature:

```bash
python3 make_wrfinput.py
./erf_exec WRFInput_SEBSkinInit.i erf.radiation.seb_init_t_deep_from_file=true max_step=2 > run.log
python3 check_skin.py file plt2d00002 --field deep --log run.log --count 78
```

## Expected Results

The log holds one line per level that took its file's values:

```
ERF::init_seb_skin_from_file: level 0: the skin temperature of 78 land column(s) starts from the initial file (280 K to 287.25 K); water columns keep erf.rad_t_sfc = 300 K. The restoring temperature stays erf.radiation.seb_t_deep_default = 290 K in every column, up to 10 K from the values set; the skin relaxes toward it.
```

| Column | `seb_t_sfc` at step 0 [K] |
| --- | --- |
| land, for example (0, 0), (5, 2), (8, 5), (11, 7) | `280 + 0.5 i + 0.25 j`: 280, 283, 285.25, 287.25 |
| water, `i >= 6` and `j <= 2` | 300 |
| any column with the input `false` | 300 |

The check is exact: no tolerance.

The restoring temperature of the balance stays the one value `erf.radiation.seb_t_deep_default`
(290 K in the deck), so the skin relaxes from the file's value toward it during a run. The log line
prints the largest distance between the two, here 290 - 280 = 10 K. This case stops at step 0
and does not test the relaxation.

With a second level, refined by 2 and without a file of its own, the fine level is filled
from level 0 by the interpolation ERF already used. For this linear field and away from
the domain faces that gives the same formula at the fine column centers,
`280 + 0.5 x + 0.25 y` with `x = (I + 0.5) / 2 - 0.5` for the fine column `I`.

With a second level that reads `wrfinput_d02` and `wrflowinp_d02`, a fine land column holds
`290 + 0.125 I + 0.0625 J`. The two fine water columns are not written: they keep what the
fill from level 0 gave them, here 283 K, the `TSK` of the level 0 land column (5, 2) under
them. The log then holds a second line, for level 1, with 62 land columns.

With `wrfinput_d02` and no `erf.nc_low_file_1`, ERF replaces the skin temperature the level
read by its parent's, interpolated, so that the level follows level 0's `wrflowinp` file.
That is not the file's field, so level 1 keeps the fill from level 0 and the log says so.

Not covered here:

- A metgrid file (`SKINTEMP`). It reaches the balance through the same field as `TSK`, but
  ERF has no metgrid test with a generated file.
- A wrfinput file without `TSK`. The wrfinput reader stops on it before this input is used,
  so the test of "no skin temperature" is a start without a file.
- The stop with a land-surface, urban or surface model: those need soil fields this file
  does not have.

With `erf.radiation.seb_init_t_deep_from_file = true` the log holds one more line per level
that took its file's `TMN`, and the skin's line ends after the water columns:

```
ERF::init_seb_t_deep_from_file: level 0: the deep temperature of 78 land column(s) is the TMN of the initial file (275 K to 281.25 K); water columns keep erf.radiation.seb_t_deep_default = 290 K. The balance does not evolve it.
```

| Column | `seb_t_deep` at step 0 and at step 2 [K] |
| --- | --- |
| land, for example (0, 0), (5, 2), (8, 5), (11, 7) | `275 + 0.25 i + 0.5 j`: 275, 277.25, 279.5, 281.25 |
| water, `i >= 6` and `j <= 2` | 290 |
| any column with the input `false` | 290 |

The second level without a file, and the second level with `wrfinput_d02`, hold `TMN` as
described for the skin: the linear field at the fine column centers, and
`284 + 0.0625 I + 0.125 J` on fine land with 277.25 K, the `TMN` of the level 0 land column
(5, 2), in the two fine water columns. Without `erf.nc_low_file_1` the level still takes the
`TMN` of `wrfinput_d02`: the lower boundary file does not hold a soil temperature.

A restart reads `T_deep` from `TwoStream_TDeep` of the checkpoint and not from the file. A
restart with the input from a checkpoint of a run without it holds 290 K in every column,
which is what that run used.

Not covered for the deep temperature: the stop with a land-surface, urban or surface model,
as for the skin; a regrid during the run, which uses the same carry-over as the skin; and a
nest with a file that takes its radiation from its parent (a nest of partial height), which
is filled from its parent.

The unit test `SEBDeepInit` (`Tests/Unit/Radiation`) checks that the balance keeps a
per-column deep temperature and restores the skin toward it (the two columns' skins differ
by `dt (2 pi / tau) dT_deep` after one step), that it is one of the fields a new level
takes from its parent, and the four checkpoint cases.

The unit test `SEBSkinInit` (`Tests/Unit/Radiation`) checks the copy itself, with water, a
building's column, ghost columns, two boxes and values that are not usable.

## CTest

| Test | What it runs | Passes when |
| --- | --- | --- |
| `WRFInput_SEBSkinInit` | the deck, with `TSK = -999` in the water column (8, 1) | 78 land columns hold `TSK`, 18 water columns hold 300 K (a water value is not read), and the log reports 78 and the largest distance from the restoring temperature, 10 K |
| `WRFInput_SEBSkinInit_Off` | `erf.radiation.seb_init_skin_from_file=false` | every column holds 300 K and the log has no report line |
| `WRFInput_SEBSkinInit_TaggedLevel` | a file without water (`--all-land`); a second level made by tagging at the start on level 0 columns 2..5 x 2..5 | all 96 columns of level 0 hold `TSK`; the 64 fine columns hold the linear field at their centers; the log has no line for level 1 |
| `WRFInput_SEBSkinInit_Restart` | two steps with a checkpoint after the first, then a restart from it | the restart's skin at step 2 is that of the first run, bit for bit; its log has the NOTE that the input is not used on a restart and no report line |
| `WRFInput_SEBSkinInit_NoFile` | `erf.init_type=Uniform` | the run stops; the message names the input and says the file of level 0 gave no skin temperature |
| `WRFInput_SEBSkinInit_NotTwoStream` | `erf.radiation_model=None` | the run stops; the message names the input and `erf.radiation_model = TwoStream` |
| `WRFInput_SEBSkinInit_Nested` | `--nest`; a second level from `wrfinput_d02` and `wrflowinp_d02` | level 0 as the first test; 62 fine land columns hold the `TSK` of `wrfinput_d02`, the 2 fine water columns 283 K; the log reports 78 and 62 |
| `WRFInput_SEBSkinInit_NestedParentLow` | `--nest`; a second level from `wrfinput_d02`, no `erf.nc_low_file_1` | level 0 as the first test; only level 0 is reported; the log says level 1 takes the skin of level 0, interpolated |
| `WRFInput_SEBSkinInit_AboveRange` | `TSK = 350` in the land column (8, 5), above `seb_prognostic_t_max_k = 340` | the run stops; the message names the input, counts 1 land column and prints the range 200 to 340 K |
| `WRFInput_SEBSkinInit_RangeRaised` | `TSK = 350` there and `erf.radiation.seb_prognostic_t_max_k=360` | the run finishes and the column holds 350 K |
| `WRFInput_SEBSkinInit_NotPrognostic` | `erf.radiation.seb_prognostic_enable=false` | the run stops; the message names the input and `seb_prognostic_enable` |
| `WRFInput_SEBSkinInit_MPI` | the deck on 4 ranks (2 boxes, so 2 ranks hold none); MPI builds only | as the first test |
| `WRFInput_SEBSkinInit_FlagOnLand_MPI` | `TSK = -999` in one land column of each box, on 2 ranks; MPI builds only | the run stops and the message counts 2 land columns, the sum over the ranks |
| `WRFInput_SEBSkinInit_Nested_MPI` | the nested test on 3 ranks with `amr.max_grid_size=4`; MPI builds only | as the nested test |
| `WRFInput_SEBSkinInit_Restart_MPI` | the restart test on 3 ranks; MPI builds only | as the restart test |
| `WRFInput_SEBDeepInit` | the deep input on, two steps, with `TMN = 350` in the water column (8, 1) | 78 land columns hold `TMN` at step 0 and at step 2, 18 water columns hold 290 K; the skin is as in the first test; the log reports 78 and no longer says the restoring temperature stays one value |
| `WRFInput_SEBDeepInit_Off` | the deck (the deep input off), two steps | every column holds 290 K at step 0 and step 2; the log has no report line |
| `WRFInput_SEBDeepInit_SkinOff` | the deep input on, `seb_init_skin_from_file=false` | `seb_t_deep` holds `TMN` on land; every skin holds 300 K |
| `WRFInput_SEBDeepInit_TaggedLevel` | as the skin's tagged test, with the deep input | all 96 columns of level 0 hold `TMN`; the 64 fine columns hold the linear field at their centers; no line for level 1 |
| `WRFInput_SEBDeepInit_Nested` | `--nest`, `wrfinput_d02` and `wrflowinp_d02` | 62 fine land columns hold the `TMN` of `wrfinput_d02`, the 2 fine water columns 277.25 K; the log reports 78 and 62 |
| `WRFInput_SEBDeepInit_NestedParentLow` | `--nest`, no `erf.nc_low_file_1` | as the nested test for `TMN`, while the log says the skin of level 1 is that of level 0 |
| `WRFInput_SEBDeepInit_Restart` | two steps with a checkpoint after the first, then a restart | the restart's `seb_t_deep` and `seb_t_sfc` at step 2 are those of the first run, bit for bit; `seb_t_deep` holds `TMN`; the restart's log has no report line |
| `WRFInput_SEBDeepInit_OlderCheckpoint` | the first run without the deep input, the restart with it | every column of the restart holds 290 K; its skin is that of the first run, bit for bit; the log has the NOTE that the checkpoint has no `TwoStream_TDeep` |
| `WRFInput_SEBDeepInit_DeepPrognostic` | the deep input with `seb_deep_prognostic_enable` | `seb_t_deep` holds `TMN` on land at step 0; at step 2 every land column is within 0.01 K of `TMN` and at least one differs from it; the log says it is the first value of the deep reservoir |
| `WRFInput_SEBDeepInit_DeepPrognostic_Restart` | the same through a restart | `seb_t_deep` and `seb_t_sfc` at step 2 are those of the first run, bit for bit, and `seb_t_deep` is within 0.01 K of `TMN` |
| `WRFInput_SEBDeepInit_RestartInputOff` | the first run with the deep input, the restart without it | the first run holds `TMN` at step 2; every column of the restart holds 290 K; the log has the NOTE that the checkpoint's `TwoStream_TDeep` is not read |
| `WRFInput_SEBDeepInit_Nested_Restart` | the nested test through a restart (two-way coupling) | `seb_t_deep` and `seb_t_sfc` at step 2 are those of the first run on both levels, bit for bit; the nest holds its file's `TMN`; a level 0 land column under the nest holds the mean of the nest's land `TMN` over it |
| `WRFInput_SEBDeepInit_NoTMN` | a file without `TMN` (`--no-tmn`) | the run stops; the message names the input and says the file of level 0 gave no usable `TMN` |
| `WRFInput_SEBDeepInit_NoFile` | `erf.init_type=Uniform`, the skin input off | the same stop |
| `WRFInput_SEBDeepInit_FlagOnWater` | `TMN = -999` in the water column (8, 1) | the same stop: the reader does not keep a field with a value below 0 K, on water too |
| `WRFInput_SEBDeepInit_NotTwoStream` | `erf.radiation_model=None`, the skin input off | the run stops; the message names the input and `erf.radiation_model = TwoStream` |
| `WRFInput_SEBDeepInit_AboveRange` | `TMN = 350` in the land column (8, 5) | the run stops; the message names the input, counts 1 land column and prints the range 200 to 340 K |
| `WRFInput_SEBDeepInit_BelowRange` | `TMN = 150` there, below `seb_prognostic_t_min_k = 200` | the same stop |
| `WRFInput_SEBDeepInit_RangeRaised` | `TMN = 350` there and `erf.radiation.seb_prognostic_t_max_k=360` | the run finishes and the column holds 350 K |
| `WRFInput_SEBDeepInit_NotPrognostic` | `erf.radiation.seb_prognostic_enable=false`, the skin input off | the run stops; the message names the input and `seb_prognostic_enable` |
| `WRFInput_SEBDeepInit_MPI` | the first deep test on 4 ranks; MPI builds only | `TMN` on land at step 0 and step 2 |
| `WRFInput_SEBDeepInit_AboveRange_MPI` | `TMN = 350` in one land column of each box, on 2 ranks; MPI builds only | the run stops and the message counts 2 land columns |
| `WRFInput_SEBDeepInit_Nested_MPI` | the nested test on 3 ranks with `amr.max_grid_size=4`; MPI builds only | as the nested test |
| `WRFInput_SEBDeepInit_Restart_MPI` | the restart test on 3 ranks; MPI builds only | as the restart test |

All have the labels `regression;wrfinput` and are registered only with
`ERF_ENABLE_NETCDF` and a Python that has `netCDF4` and `numpy`.
