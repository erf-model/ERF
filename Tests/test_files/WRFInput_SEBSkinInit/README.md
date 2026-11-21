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

## Files

| File | Role |
| --- | --- |
| `WRFInput_SEBSkinInit.i` | One level, 12 x 8 x 8 cells in two boxes, `max_step = 0`, two-stream radiation with `seb_prognostic_enable`, `erf.rad_t_sfc = 300`, the input on; writes `seb_t_sfc` to the 2D plotfile `plt2d00000` |
| `make_wrfinput.py` | Writes `wrfinput_d01` and `wrflowinp_d01`; `--all-land` leaves out the water; `--land VALUE` puts the value in the land column (8, 5) and `--land-west VALUE` in the land column (2, 5) of the other box; `--water VALUE` puts it in the water column (8, 1); `--nest` also writes `wrfinput_d02` and `wrflowinp_d02` for a second level |
| `check_skin.py` | Compares `seb_t_sfc` of every column of the plotfile with the formulas in `make_wrfinput.py`; also compares two plotfiles and looks for a message in a log |

## Data and acknowledgements

None is downloaded and none is redistributed. `make_wrfinput.py` writes the files:

- `TSK = 280 + 0.5 i + 0.25 j` [K], with `i`, `j` the column indices from zero, so 280 K to
  287.25 K. Every value is exact in single precision.
- `LANDMASK = 0` (water) where `i >= 6` and `j <= 2`: 18 of the 96 columns, in both boxes
  of the deck (the boxes meet between `i = 7` and `i = 8`). The other 78 are land.
- `SST = 285` [K].
- With `--nest`, `wrfinput_d02` covers the level 0 columns 2..5 x 2..5 with 8 x 8 columns
  refined by 2. Its `TSK = 290 + 0.125 I + 0.0625 J` [K], with `I`, `J` the fine column
  indices of the whole domain (4..11), so it differs from the level 0 field everywhere. Its
  `LANDMASK = 0` in the fine columns (10, 5) and (11, 5): 62 of the 64 are land.

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

## Expected Results

The log holds one line per level that took its file's values:

```
ERF::init_seb_skin_from_file: level 0: the skin temperature of 78 land column(s) starts from the initial file (280 K to 287.25 K); water columns keep erf.rad_t_sfc = 300 K.
```

| Column | `seb_t_sfc` at step 0 [K] |
| --- | --- |
| land, for example (0, 0), (5, 2), (8, 5), (11, 7) | `280 + 0.5 i + 0.25 j`: 280, 283, 285.25, 287.25 |
| water, `i >= 6` and `j <= 2` | 300 |
| any column with the input `false` | 300 |

The check is exact: no tolerance.

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

The unit test `SEBSkinInit` (`Tests/Unit/Radiation`) checks the copy itself, with water, a
building's column, ghost columns, two boxes and values that are not usable.

## CTest

| Test | What it runs | Passes when |
| --- | --- | --- |
| `WRFInput_SEBSkinInit` | the deck | 78 land columns hold `TSK`, 18 water columns hold 300 K, and the log reports 78 |
| `WRFInput_SEBSkinInit_Off` | `erf.radiation.seb_init_skin_from_file=false` | every column holds 300 K and the log has no report line |
| `WRFInput_SEBSkinInit_AllLand` | a file without water | all 96 columns hold `TSK` |
| `WRFInput_SEBSkinInit_TaggedLevel` | a file without water; a second level made by tagging at the start on level 0 columns 2..5 x 2..5 | level 0 holds `TSK`; the 64 fine columns hold the linear field at their centers; the log has no line for level 1 |
| `WRFInput_SEBSkinInit_Restart` | two steps with a checkpoint after the first, then a restart from it | the restart's skin at step 2 is that of the first run, bit for bit; its log has the NOTE that the input is not used on a restart and no report line |
| `WRFInput_SEBSkinInit_NoFile` | `erf.init_type=Uniform` | the run stops; the message names the input and says the file of level 0 gave no skin temperature |
| `WRFInput_SEBSkinInit_NotTwoStream` | `erf.radiation_model=None` | the run stops; the message names the input and `erf.radiation_model = TwoStream` |
| `WRFInput_SEBSkinInit_Nested` | `--nest`; a second level from `wrfinput_d02` and `wrflowinp_d02` | level 0 as the first test; 62 fine land columns hold the `TSK` of `wrfinput_d02`, the 2 fine water columns 283 K; the log reports 78 and 62 |
| `WRFInput_SEBSkinInit_NestedParentLow` | `--nest`; a second level from `wrfinput_d02`, no `erf.nc_low_file_1` | level 0 as the first test; only level 0 is reported; the log says level 1 takes the skin of level 0, interpolated |
| `WRFInput_SEBSkinInit_FlagOnLand` | `TSK = -999` in the land column (8, 5) | the run stops; the message names the input and counts 1 land column |
| `WRFInput_SEBSkinInit_AboveRange` | `TSK = 350` in that column, above `seb_prognostic_t_max_k = 340` | the same, and the message prints the range 200 to 340 K |
| `WRFInput_SEBSkinInit_RangeRaised` | `TSK = 350` there and `erf.radiation.seb_prognostic_t_max_k=360` | the run finishes and the column holds 350 K |
| `WRFInput_SEBSkinInit_NotPrognostic` | `erf.radiation.seb_prognostic_enable=false` | the run stops; the message names the input and `seb_prognostic_enable` |
| `WRFInput_SEBSkinInit_FlagOnWater` | `TSK = -999` in the water column (8, 1) | the run finishes with the values of the first test: a water value is not read |
| `WRFInput_SEBSkinInit_MPI` | the first test on 4 ranks (2 boxes, so 2 ranks hold none); MPI builds only | as the first test |
| `WRFInput_SEBSkinInit_SmallBoxes_MPI` | the first test on 3 ranks with `amr.max_grid_size=4`; MPI builds only | as the first test |
| `WRFInput_SEBSkinInit_FlagOnLand_MPI` | `TSK = -999` in one land column of each box, on 2 ranks; MPI builds only | the run stops and the message counts 2 land columns, the sum over the ranks |
| `WRFInput_SEBSkinInit_Nested_MPI` | the nested test on 3 ranks with `amr.max_grid_size=4`; MPI builds only | as the nested test |
| `WRFInput_SEBSkinInit_Restart_MPI` | the restart test on 3 ranks; MPI builds only | as the restart test |

All have the labels `regression;wrfinput` and are registered only with
`ERF_ENABLE_NETCDF` and a Python that has `netCDF4` and `numpy`.
