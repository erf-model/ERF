# Terrain-following inflow profile

Regression cases for `<face>.inflow_profile` (`log_law` or `file`) and for
`erf.input_sounding_theta_above_ground`. The profile option imposes a vertical
profile of the wind, and of theta and tke where given, at each boundary cell's
own height above the ground beneath it. The inflow follows the terrain along the
face, and the terrain stays on the boundary, as WindNinja's inlets do.
`xlo.dirichlet_file` instead holds one value per grid level along the whole face,
from heights above the domain floor.

## Setup

- **Domain:** 1600 m x 800 m x 800 m, 32 x 16 x 32 cells (50 m horizontally, 25 m
  vertically over flat ground), basic terrain following.
- **Boundaries:** periodic in y, `Inflow` at xlo and `Outflow` at xhi. Surface
  layer with z0 = 0.1 m, slip-wall top with Rayleigh damping of w in the top 200 m.
- **Solver:** anelastic, fixed dt = 1 s, one-equation k RANS with
  `erf.dirichlet_k` and `erf.init_tke_from_ustar` plus
  `erf.init_tke_at_wall_value`. Neutral, theta = 300 K. The terrain-fitted
  projection is GMRES with the FFT preconditioner, so these cases need a build
  with FFT; CTest registers them only when `ERF_ENABLE_FFT` is on.
- **Start:** the interior starts from a uniform 10 m/s westerly, so only the
  inflow face can turn the first column into a log law.
- **Profile:** 10 m/s at 10 m above the local ground over z0 = 0.1 m, from 270
  degrees. The tke is u*^2/Cmu0^2 at the ground, tapering to zero at 700 u*.

`gen_inputs.py` writes the terrain, sounding and profile files. Every sloping
terrain puts the ground on the inflow face at 150 m or more. There, a profile
applied by level or at absolute height differs clearly from one applied at the
height above the ground.

| CTest | Terrain | Ground on the inflow face |
|---|---|---|
| `RANS_InflowProfile_Plateau` | flat, raised to 150 m | 150 m |
| `RANS_InflowProfile_Incline` | rising downstream at 0.15 | 150 m (390 m at the outflow) |
| `RANS_InflowProfile_Decline` | falling downstream at 0.15 | 290 m (50 m at the outflow) |
| `RANS_InflowProfile_CrossRidge` | cosine ridge across the face | 150 to 230 m along y |
| `RANS_InflowProfile_CrossRidge_File` | the same, profile read from `inflow_profile_crossridge.txt` (`# z u v T tke`) | 150 to 230 m |
| `RANS_InflowProfile_Parity_Flat` | flat on the floor | 0 m |

## Checks

`check_inflow_profile.py` reads the first column next to the inflow face (i = 0)
after 60 steps. It takes the height above the ground from `z_phys`, which is exact
for basic terrain following. It requires:

- the wind speed follows the log law at the height above the local ground,
  above the wall cell and below the Rayleigh layer (median error at most 1 %,
  90th percentile at most 5 %);
- in the lowest 200 m, that law fits at least twice as well as the log law at the
  height above the floor (a profile applied at absolute heights);
- in the lowest 200 m, it also fits better than the log law at the flat-mesh
  level heights (a profile applied by level, as `dirichlet_file` is);
- the wind blows from 270 degrees;
- the tke is 0.8 to 1.25 times the profile's;
- every field is finite, and the maximum speed, |w| and tke stay bounded.

`RANS_InflowProfile_Parity_Flat` runs `inputs_parity` twice: once with
`xlo.dirichlet_file = inflow_levels_flat.txt`, once with the same rows through
`xlo.inflow_profile = file`. On flat ground on the floor the two lookups
coincide, so the plotfiles must be identical (fcompare with zero tolerance). The
comparison is sensitive: the same rows with u scaled by 1.001 make fcompare fail
(x_velocity relative difference 1.1e-3).

## Results

Measured on 2026-09-12 (Release with FFT, 2 ranks, 60 steps), after the
anelastic scalar-advection fix (hgopalan/ERF#415). The rows marked
*dirichlet_file* are negative controls, not CTest cases: the same log law applied
by level (`inputs_parity` with `xlo.dirichlet_file=inflow_levels_flat.txt` on the
named terrain).

| Case | Median error | 90th pct | Error ratio vs floor | Error ratio vs by-level | tke / profile | Max tke [m2/s2] | Checks failed |
|---|---|---|---|---|---|---|---|
| Plateau | 2.2e-5 | 4.5e-5 | 0.0004 | 0.0015 | 0.98 | 4.5 | 0 |
| Incline | 0.55 % | 1.3 % | 0.090 | 0.44 | 0.98 | 14.3 | 0 |
| Decline | 0.69 % | 2.2 % | 0.058 | 0.14 | 0.98 | 2.5 | 0 |
| CrossRidge | 2.3e-5 | 5.2e-5 | 0.0003 | 0.0009 | 0.98 | 5.9 | 0 |
| CrossRidge_File | 4.9e-4 | 6.7e-4 | 0.0043 | 0.014 | 0.98 | 5.9 | 0 |
| Plateau, *dirichlet_file* | 2.6 % | | 0.36 | 87 | 0.51 | | 3 |
| Incline, *dirichlet_file* | 3.2 % | | 0.52 | 4.7 | 0.54 | | 4 |
| CrossRidge, *dirichlet_file* | 3.6 % | | 0.43 | 115 | 0.57 | | 3 |

The parity pair gave identical plotfiles.

**Stability over 30 minutes (1800 steps).** Every case is steady by the 5-minute
plotfile and stays there:

| Case | Max speed [m/s] | Max tke [m2/s2] | Max \|w\| [m/s] |
|---|---|---|---|
| Plateau | 19.08 | 2.57 | 0.04 |
| Incline | 28.13 | 12.84 | 3.93 |
| Decline | 18.39 | 2.83 | 1.18 |
| CrossRidge | 19.08 | 3.60 | 0.04 |

The highest tke in each run is at step 0, before the flow adjusts to the
projection.

**Why the incline carries more tke.** The slip-wall top is fixed at 800 m, so the
column thins as the ground rises and the anelastic flow speeds up to carry the
same mass. At the outflow end of the earlier 50 m to 290 m incline, the
column-mean wind was 1.44 times its inflow value, against a depth ratio of 1.43.
The wall-cell tke followed the local equilibrium u*^2/Cmu0^2 to within 0.97 to
1.12 along the whole slope, with an eddy viscosity of only 6.6 m2/s at the
maximum. This is the equilibrium for the faster wind, not a runaway.

## Theta above the ground: `inputs_theta`

`erf.input_sounding_theta_above_ground = true` interpolates the input sounding's
theta at the height above the local ground, as the inflow profile's `T` column
is. The pressure keeps the sounding's value at the physical height, the density
follows from it, and each column is rebalanced hydrostatically.

The case uses the incline under a 1200 m domain (48 levels) with a capping
inversion: theta is 300 K up to 300 m above the ground, rises by 8 K to 400 m and
by 3 K/km above (`gen_theta_inputs.py` writes `input_sounding_inversion` and
`inflow_profile_inversion.txt`). It is compressible, without turbulence closure
and with slip walls, so it runs in every build.

| CTest | Check |
|---|---|
| `InflowProfile_ThetaAboveGround` | `check_theta_above_ground.py`: at step 0, theta in every cell and in the first column equals the sounding at the height above the ground (max error 1e-6 K), each column is in discrete hydrostatic balance (residual at most 1e-9 of rho g), and the inversion does move with the ground (theta at physical height differs by at least 2 K in at least a tenth of the cells); after 4 steps every field is finite and \|w\| stays below 5 m/s |
| `InflowProfile_ThetaAboveGround_FlagOff` | `check_theta_flag_off.py`: the same deck with the flag off must start theta at physical heights (median error at least 0.5 K, max at least 4 K), so the first test cannot pass without the flag |

The reference in both checkers is the sounding as ERF stores it: resampled onto the
ground, the flat mesh's cell-centre heights and the top, then interpolated
linearly. Compared with the exact piecewise-linear sounding, that resampling rounds
the inversion base by up to 0.5 K in the cell next to it, with the flag on or off.

Measured on 2026-09-12: with the flag, the step-0 theta error is 2.3e-13 K
(first column 1.7e-13 K), the hydrostatic residual 3.9e-13, and max |w| is 2.1 m/s
after 4 steps (a 10 m/s wind lifted over the 0.15 slope). The inversion moves
by at least 2 K in 35 % of the cells. With the flag off, the error is 1.0 K median
and 8.8 K max.

On a terrain-fitted mesh started from an input sounding, ERF prints a warning
when an inflow face has a profile and `erf.input_sounding_wind_above_ground` is
off, or the profile has a `T` column and `erf.input_sounding_theta_above_ground`
is off.

## Running by hand

```
python3 gen_inputs.py        # only to regenerate the tables
python3 gen_theta_inputs.py
mpirun -np 2 erf_exec inputs_inflow max_step=60 erf.plot_int_1=60 erf.terrain_file_name=terrain_incline.txt
python3 check_inflow_profile.py plt00060
mpirun -np 2 erf_exec inputs_theta max_step=4 erf.plot_int_1=4
python3 check_theta_above_ground.py plt00004
```

The check scripts need `erf_plotfile.py` and `rans_checks.py` from the parent
`Canonical_RANS` folder, which CTest copies next to them.
