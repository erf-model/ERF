# Anelastic advection of the slow scalars

## Purpose

With the anelastic integrator, `erf_slow_rhs_post` advects the slow scalars (TKE, moisture,
passive scalars) with the momenta that `AdvectionSrcForRho` sets in `erf_slow_rhs_pre` from
the start of each stage, in the form `ax rho u / mf_uy`, `ay rho v / mf_vx` and
`az Omega / (mf_mx mf_my)`; `AdvectionSrcForScalars` divides the flux differences by `detJ`.
These momenta were projected by the previous stage or step. These cases check the scalars
against exact answers. Theta is advected in `erf_slow_rhs_pre` with the same faces, so it
serves as a reference.

## Files

| File | What it is |
|---|---|
| `inputs_mapfac_anelastic` | Blob in a uniform wind, uniform mesh, `erf.test_mapfactor = true` (m = 0.5), MLMG projection |
| `inputs_mapfac_compressible` | The same blob with the compressible integrator, no substepping: the reference for the case above |
| `inputs_fitted` | Blob on a terrain-fitted mesh (`VariableDz`) over a 1e-6 m hill, in effect flat ground, 12 levels from 20 m stretched by 1.2 |
| `inputs_stretched` | The same blob and stretching without terrain (`StretchedDz`); needs the FFT projection |
| `inputs_hill` | Uniform scalar and uniform theta in the wind over a 100 m Witch of Agnesi hill, the same stretching |
| `inputs_inertial` | Blob in an inertial oscillation: a uniform wind that turns in time, with an exact position |
| `check_inertial_drift.py` | Blob centroid and wind against the exact inertial oscillation |
| `check_cfl_dt.py` | The time step from `erf.cfl` with map factors, on `inputs_mapfac_anelastic` |
| `check_mapfac_parity.py` | Blob centroids on every level against the analytic drift `U t m / dx`, `V t m / dy` and the compressible run |
| `check_scalar_centroid.py` | Blob x-centroid on every level against `x_c + U t` |
| `check_uniform_scalar.py` | Scalar and theta on every AMR level stay at their starting values |

The plotfile reader is `Exec/Tools/erf_plotfile.py`.

## Data

None: every case starts from the input deck.

## How to run

```bash
ctest -R AnelasticScalarAdvection
```

or by hand, from a run directory:

```bash
mpiexec -n 2 erf_exec inputs_hill max_step=30 erf.plot_int_1=30
python3 check_uniform_scalar.py plt00030 1.0e-4
```

## Expected Results

Measured on 2 ranks. "Before" is development at 8bdde5f8a, "after" is the current code.

| CTest `AnelasticScalarAdvection_` | Check | Before | After |
|---|---|---|---|
| `MapFactor` | centroid within 0.25 cell | 2.4 cells behind in x, 1.2 in y | on the analytic drift |
| `Fitted`, `Stretched` | centroid within one cell (31.25 m) | 40 to 457 m off | within one cell |
| `Hill`, `Hill_SplitZ` | relative error below 1e-4 | 0.28 | 1.6e-7, the same as theta |
| `Hill_MidPoint` | relative error below 1e-4 | 0.28 | 8.4e-8 |
| `Hill_GmresFFT` | relative error below 1e-4 | 0.28 | 1.1e-8 |
| `Hill_MapFactor` | relative error below 1e-4 | 0.11 | 7.3e-9 |
| `Hill_TwoLevels` | relative error below 1e-4 | 0.31 | 3.0e-5, the same as theta |
| `Hill_PartialDepth` | relative error below 1e-4 | 0.28 | 1.6e-7 |
| `Hill_TwoWay` | relative error below 1e-3 | 0.31 | 4.2e-4; theta 1.5e-4 |
| `Inertial_RK2`, `Inertial_MidPoint` | blob within 0.5 m of the exact position | 3.4 m | 0.072 m |
| `CFL_MapFactor` | time after 5 steps within 1 % of 10.417 s | 5.208 s | 10.417 s |

## CTest

The cases that use multigrid (`erf.terrain_poisson_solver = mlmg` or the uniform-mesh MLMG
projection) are registered in every build with Python and double precision, so they run in
CI. `Stretched` and `Hill_GmresFFT` need `ERF_ENABLE_FFT`.
