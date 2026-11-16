# WRFInput_SurfaceFields

## Purpose

A wrfinput file carries land and soil fields. ERF asked for five of them (`LAI`, `VEGFRA`,
`TMN`, `TSLB`, `SMOIS`) only when a land model was on, and never read `ZNT`, `ALBEDO`,
`ALBBCK`, `SEAICE`, `SNOW` or `SNOWH`. With the two-stream force-restore surface energy
balance and no land model, the reader now takes all eleven, plus the soil layer depths `ZS`
and thicknesses `DZS`, and keeps them per level in `WRFSurfaceFields`. Each one is optional.
A field with a value outside its range is reported and left out; the run goes on. The
fields are written to the checkpoint and read back on a restart.

This case checks the reader only. Nothing uses the fields yet.

## Files

| File | Role |
| --- | --- |
| `WRFInput_SurfaceFields.i` | One level, 12 x 8 x 8 cells, `max_step = 0`, two-stream radiation with `seb_prognostic_enable` |
| `make_wrfinput.py` | Writes `wrfinput_d01` and `wrflowinp_d01`; `--omit NAME ...` leaves variables out; `--bad NAME VALUE` puts the value in one column of a field; `--fill NAME VALUE` puts it in every column; `--surface-zs` starts `ZS` at zero; `--nz N` sets the atmospheric levels; `--zero-layers` leaves `ZS` and `DZS` at zero |
| `check_log.py` | Compares the reader's report in the log with the formulas in `make_wrfinput.py` |

## Data

None is downloaded. `make_wrfinput.py` writes the files. Every field varies in both
horizontal directions and with the soil layer, for example
`TSLB = 290 - 2 k + 0.5 i + 0.0625 j` and `SMOIS = 0.1 + 0.05 k + 0.003 j + 0.0002 i`
(`i`, `j` column indices, `k` soil layer, all from zero). `SEAICE` is 1 in one column.

## How to run

Needs ERF built with NetCDF, and Python with `netCDF4` and `numpy`.

```bash
python3 make_wrfinput.py
./erf_exec WRFInput_SurfaceFields.i > run.log
python3 check_log.py full run.log
```

## Expected Results

The log holds one line per field, for example

```
   wrfinput surface field TSLB at level 0: 4 components, min 284, max 295.938 [K]
   wrfinput soil layer depths ZS at level 0: 0.05 0.25 0.7 1.5 [m]
```

| Field | Components | Min | Max |
| --- | --- | --- | --- |
| `LAI` | 1 | 0.5 | 3.46875 |
| `VEGFRA` | 1 | 10 | 47.75 |
| `TMN` | 1 | 280 | 285.9375 |
| `ZNT` | 1 | 0.01 | 0.091 |
| `ALBEDO` | 1 | 0.1 | 0.327 |
| `ALBBCK` | 1 | 0.15 | 0.274 |
| `SEAICE` | 1 | 0 | 1 |
| `SNOW` | 1 | 0 | 15.375 |
| `SNOWH` | 1 | 0 | 0.0755 |
| `TSLB` | 4 | 284 | 295.9375 |
| `SMOIS` | 4 | 0.1 | 0.2732 |

The file holds single-precision values, so the log shows for example 0.3269999921 for 0.327. The check allows a relative difference of 2e-5.

With a second level and `erf.v=1` the log also holds one line per field each time the
level is put on new grids, for example

```
   wrfinput surface field TSLB at level 1 on the new grids: 4 components, 64 of 96 columns from level 0, min 286.125, max 294.8125 [K]
```

The level takes the level 0 column each of its columns lies in, so the range is that of
the level 0 columns under it. The regrid legs pin the counts and the ranges; the unit
tests check that a column the old grids held keeps the level's own value.

The log check pins the component count and the range of each field, not where each value
lands. A reversed order of the `TSLB` or `SMOIS` layers would pass it; a reversed `ZS` or
`DZS` would not. The unit test `WRFSurfaceFields` (`Tests/Unit/Initialization`) checks
every column, the ghost columns and the layer order of the copy into the holder.

## CTest

| Test | What it runs | Passes when |
| --- | --- | --- |
| `WRFInput_SurfaceFields` | every field in the file | all eleven fields, `ZS` and `DZS` are reported with the values above |
| `WRFInput_SurfaceFields_Absent` | `LAI`, `SEAICE`, `SNOWH`, `SMOIS`, `DZS` left out | the run finishes, these five are skipped and the others are reported |
| `WRFInput_SurfaceFields_BalanceOff` | `erf.radiation.seb_prognostic_enable=false` | none of the fields is asked for |
| `WRFInput_SurfaceFields_Drop_Albedo` | one `ALBEDO` column is 1.5 | the run finishes; `ALBEDO` is reported as not held with 1 value outside 0 to 1; the other fields are reported |
| `WRFInput_SurfaceFields_Drop_FlagValues` | `SNOWH` is -999 in every column and one `SEAICE` column is -0.001 | the run finishes; both are reported as not held (96 values and 1 value); the other fields are reported |
| `WRFInput_SurfaceFields_ThinColumn` | a file and a grid with 3 atmospheric levels and 4 soil layers | all four soil layers are reported, as in the first test |
| `WRFInput_SurfaceFields_ZeroLayers` | `ZS` and `DZS` declared and left at zero | the run finishes, both are reported as not given and the eleven fields are read |
| `WRFInput_SurfaceFields_SurfaceZS` | `ZS` is 0, 0.1, 0.4, 1.0 | `ZS` is held with these values |
| `WRFInput_SurfaceFields_Restart` | the first test with a checkpoint, then a restart from it | the restart reports every field as restored with the values above, and does not read the file again |
| `WRFInput_SurfaceFields_Restart_OtherGrids` | the same, the restart with `amr.max_grid_size=4` | the same report on the smaller boxes |
| `WRFInput_SurfaceFields_TaggedLevel` | a second level, refined by 2, made by tagging at the start on level 0 columns 2..5 x 2..5; no file of its own | level 1 reports every field with 64 of 64 columns from level 0 and the range of those level 0 columns |
| `WRFInput_SurfaceFields_Regrid` | three steps; the regrid after step 1 makes level 1 on the same columns, the regrid after step 2 moves and widens it to columns 4..9 x 2..5 | level 1 reports every field twice: 64 of 64 columns from level 0, then 64 of 96, each with the range of the level 0 columns under it; no field is dropped |
| `WRFInput_SurfaceFields_Regrid_MPI` | the regrid test on 3 ranks with `amr.max_grid_size=4`; MPI builds only | the same report as on one rank |
| `WRFInput_SurfaceFields_MPI` | the first test on 4 ranks (2 boxes, so 2 ranks hold none); MPI builds only | the same report as the first test |
| `WRFInput_SurfaceFields_Drop_Albedo_MPI` | the albedo of 1.5 on 2 ranks; the second rank holds the column; MPI builds only | the same report as on one rank, printed by the first rank |
| `WRFInput_SurfaceFields_Restart_MPI` | the restart test on 3 ranks; MPI builds only | the same report as on one rank |

All have the labels `regression;wrfinput` and are registered only with
`ERF_ENABLE_NETCDF` and a Python that has `netCDF4` and `numpy`.
