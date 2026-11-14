# TwoStream_SEBSnow

## Purpose

Checks the force-restore snow layer of the two-stream surface energy balance
(`erf.radiation.seb_snow_enable`, `Source/Radiation/TwoStream/ERF_ForceRestoreSnow.H`)
against the closed forms of the scheme it implements: the snow routine of ISBA as described
in Sec. 2.1 of Douville, Royer and Mahfouf (1995), for bare ground.

The column is the one of `TwoStream_SEBSurfaceLayerFluxes`: 4 x 4 x 16 cells, Kessler
microphysics, a fixed sun, surface-layer fluxes, and a skin that starts at 299 K. A snowpack
given at the start therefore melts.

## Files

| File | Role |
|---|---|
| `TwoStream_SEBSnow.i` | The deck: a 50 kg/m^2 snowpack with the paper's albedo (0.70) and critical snow water equivalent (10 kg/m^2). |
| `input_sounding` | The sounding of `TwoStream_SEBSurfaceLayerFluxes`. |
| `check_two_stream_seb_snow.py` | Reads the 2D plotfiles of the legs with `amrex_fextract` and checks them per cell. |
| `../../RunTwoStreamSEBSnow.cmake` | Runs the legs, the two `fcompare` comparisons and the checker. |

## Data

None. The test uses no external data set; every expected value is a closed form of the
paper evaluated from the run's own inputs and plotfiles.

## How to run

With CTest, from a build directory configured with `-DERF_ENABLE_TESTS=ON` and
`-DERF_ENABLE_FCOMPARE=ON`, after building `erf_exec`, `fcompare` and `fextract`:

```bash
ctest -R TwoStream_SEBSnow --output-on-failure
```

By hand, the deck runs as it is (the melt leg):

```bash
mpiexec -n 1 ./erf_exec TwoStream_SEBSnow.i
```

## Expected Results

| Leg | Inputs changed from the deck | Check |
|---|---|---|
| `bare` | the `erf.radiation.seb_snow_*` lines removed | reference for the next two legs |
| `no_snow` | `seb_snow_initial_swe_kg_m2 = 0` | last 2D plotfile equal to the `bare` leg's (`fcompare`): the snow layer without snow changes nothing |
| `melt` | none | every step: `seb_t_sfc` = 273.16 K (Eq. 12); `seb_snow_swe - 50 = snowfall - sublimation - melt` (Eq. 6); the step's sublimation = `dt max(LE, 0) / L_v * p_ng` (Eq. 7); `seb_sw_net` = the bare leg's times `(1 - alpha_t) / (1 - alpha_g)` (Eq. 8) |
| `full_melt` | `seb_snow_initial_swe_kg_m2 = 1`, `seb_snow_albedo = 0.3` (the ground's) | after step 1: no snow left; melt + sublimation = 1 kg/m^2; `seb_t_sfc` below the bare leg's by `L_f melt / C_s` (about 16.7 K) and above 273.16 K; `seb_q_sfc` above the bare leg's by `1 / (rho_w d_s)` = 0.01 |
| `restart` | checkpoint at step 5, restart to step 10 | last 2D plotfile equal to the `melt` leg's (`fcompare`) |
| `two_level` | a second level over the middle half of x after step 5 | the melting point and the budget on every cell of the finest data; the snow and the melt total there equal the `melt` leg's; the last plotfile has level 1 |

Equation numbers are those of Douville et al. (1995). `p_ng = W_n / (W_n + W_crn)` uses the
snow at the start of the step, `alpha_t = alpha_g + (alpha_n - alpha_g) p_ng`.

The start-up checks of the snow inputs have one abort test each
(`TwoStream_SEBSnowAbort*` in `Tests/CTestList.cmake`).

Frozen precipitation reaching the snow layer is covered by the unit tests
(`ForceRestoreSnow*` in `Tests/Unit/Radiation/ERF_GTestForceRestoreSnow.cpp`); this deck's
Kessler scheme makes none.

## CTest

`TwoStream_SEBSnow` (labels `regression`, `radiation`), and the abort tests
`TwoStream_SEBSnowAbort*`.

## Reference

Douville, H., J.-F. Royer, and J.-F. Mahfouf, 1995: A new snow parameterization for the
Meteo-France climate model. Part I: validation in stand-alone experiments. *Climate
Dynamics*, 12, 21-35.
