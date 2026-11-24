# TwoStream_SEBLandMask

## Purpose

The two-stream force-restore surface energy balance is a model of soil. This case checks
that it reads the land mask: a water column (mask 0) keeps its skin temperature and soil
water, is left out of the means, and gives the radiation sweep the surface layer's sea
surface temperature as its longwave boundary. Before the change the balance ran on every
column, so water was updated as soil and the sweep used that skin over water.

An idealized run has one mask value per level (`erf.is_land`), so the case runs one deck
all land, all water, and with two levels: water on level 0 and a level 1 of land over its
western half. The third run checks the average-down of the balance's fields between the
levels, which is over land columns: the level 0 water columns under level 1 keep their
values. A level that holds both is covered by the unit test
`Tests/Unit/Radiation/ERF_GTestTwoStreamLandMask.cpp`, and so is the fill of a new level
along a shore. In the third run every level 1 column lies over water with no land next to
it, so the fill has no land value to give: the columns take the parent's held value and the
log reports their number.

## Files

| File | Role |
|---|---|
| `TwoStream_SEBLandMask.i` | The deck: 4 x 4 x 16 cells, Kessler, MRF, surface layer at the bottom, two-stream radiation with the prognostic balance. It is the deck of `TwoStream_SEBSurfaceLayerFluxes` with a surface emissivity of one. |
| `input_sounding` | Ideal sounding, surface pressure 950 hPa. |
| `check_land_mask.py` | The checker (below). |
| `../../RunTwoStreamSEBLandMask.cmake` | Runs the three legs and the checker. |

## Data

None. The case uses no external data.

## How to run

With CTest, from a build directory configured with `-DERF_ENABLE_TESTS=ON` and
`-DERF_ENABLE_FCOMPARE=ON`, after building `erf_exec` and `fextract`:

```bash
ctest -R TwoStream_SEBLandMask --output-on-failure
```

By hand: make two empty directories `land/` and `water/`, copy the deck and
`input_sounding` into each, run each leg in its directory, then run the checker from the
directory above them:

```bash
mpiexec -n 1 ./erf_exec TwoStream_SEBLandMask.i erf.is_land=1 max_step=10      # in land/
mpiexec -n 1 ./erf_exec TwoStream_SEBLandMask.i erf.is_land=0 max_step=10      # in water/
python3 check_land_mask.py --fextract amrex_fextract --land-dir land --water-dir water --steps 10
```

The two-level run takes the extra inputs listed in `RunTwoStreamSEBLandMask.cmake`
(`run_leg("${NEST_DIR}" ...)`); pass its directory with `--nest-dir`.

## Expected Results

| # | Check | Value |
|---|---|---|
| 1 | Water run, every step: `seb_t_sfc` and `seb_q_sfc` | what plotfile 0 holds, bit for bit; plotfile 0 holds 299 K and 0.0095 |
| 2 | Land run: `seb_t_sfc` moves | by more than 1e-3 K in 10 steps |
| 3 | Step 1: `seb_lw_net`(land) - `seb_lw_net`(water) | sigma (T_sea^4 - T_skin^4), relative 1e-6 |
| 4 | `T_s_mean` in `radiation_diag.dat` | finite on land, not reported on water |
| 5 | Two-level run, step 10 | level 0 `seb_t_sfc`, `seb_q_sfc`: what plotfile 0 holds, bit for bit; level 1 `seb_t_sfc` moves by more than 1e-3 K |

In check 3, T_skin = `erf.rad_t_sfc` = 299 K and T_sea = `t_surf` (`surf_pres` / p0)^(R_d/c_p),
297.1122 K in the run for `erf.most.pot_surf_temp` = 301.5 K. The checker also compares
T_sea with 301.5 (0.95)^(R_d/c_p) = 297.1137 K, a number from `input_sounding` alone, to a
relative 1e-4 (0.03 K). The run is 5e-6 away: its surface pressure is 94998.36 Pa, 1.64 Pa
below 950 hPa, because the code takes the surface pressure as the lowest cell's pressure
plus half a cell of dry air. Both runs see the same sky at the first sweep and the
emissivity is one, so the difference is the two surface emissions alone:

| Surface pressure | T_sea [K] | sigma (T_sea^4 - T_skin^4) [W/m^2] |
|---|---|---|
| 950 hPa (sounding) | 297.1137 | -11.3290 |
| 94998.36 Pa (run) | 297.1122 | -11.3377 |

With an emissivity below one the closed form would carry that factor; the checker has none.
With the balance run on every column, checks 1, 3 and 4 fail. With the plain average-down
between levels, check 5 fails.

## CTest

`TwoStream_SEBLandMask`, labels `regression;radiation`, registered in `Tests/CTestList.cmake`
(`add_test_two_stream_seb_land_mask`). It needs Python 3 and `amrex_fextract`.
