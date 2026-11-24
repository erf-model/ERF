# Canonical RANS: neutral ABL on flat ground

A neutral, barotropic Ekman layer driven by a 10 m/s geostrophic wind at
f = 1e-4 1/s over a z0 = 0.1 m surface, closed with the one-equation k RANS
model of Axell & Liungman (2001, Environ. Fluid Mech. 1, 71-106) under a
MOST surface layer. The initial potential temperature is 300 K up to 700 m
with a 3 K per 100 m inversion above, so the boundary layer stays neutral
and bounded. There are no initial perturbations: the problem is a column
and the coarse 8 x 8 horizontal grid only exists so the deck decomposes
across ranks.

| item | value |
| --- | --- |
| domain | 2560 x 2560 x 1000 m |
| grid | 8 x 8 x 64, dz = 15.625 m, first cell center 7.8 m |
| time step | 5 s fixed, anelastic with FFT |
| closure | `erf.rans_type = kEqn`, AL01 defaults, `sigma_k = 1` |
| wall k | `erf.dirichlet_k = true` (AL01 Eq. 16) |
| initial k | `erf.init_tke_from_ustar = true` |
| physics run | 12 h (8640 steps), about a minute on 2 ranks |
| smoke run | 40 steps (`ctest -R RANS_Neutral_ABL_Flat`) |

## Running

```bash
mpirun -np 2 erf_exec inputs_neutral
python3 check_neutral.py --physics plt08640 surf_hist.dat
```

`check_neutral.py` needs only the Python standard library. It reads the
plotfile with `Exec/Tools/erf_plotfile.py`, averages every field over x and y, and
prints one row per check with the measured value, the target and the
tolerance. The exit code is non-zero if any enabled check fails.

## Checks

Smoke (`--smoke`, the CTest entry, 40 steps):

| check | target | tolerance |
| --- | --- | --- |
| all planar-averaged fields finite | yes | exact |
| min KE, Kmv, diss | >= 0 | exact |
| wall distance vs cell-center height above `geometry.prob_lo[2]` | equal | 1e-8 m |
| Lturb in the three lowest cells vs the capped kappa (z + z0) | equal | 0.1 % |
| max Lturb | <= max_geom_lscale (30 m) | 1e-9 m |
| diss vs Cmu0^3 rho k^1.5 / Lturb (AL01 Eq. 19), interior cells | equal | 5 % |
| wall cell k at step start over step end, from diss and KE | 1 | 1 % |

Physics (`--physics`, adds to the smoke checks, 12 h run):

| check | target | tolerance |
| --- | --- | --- |
| u* | 0.35 m/s | +/- 0.15 m/s |
| KE at the first cell over u*^2 | 1 / Cmu0^2 = 3.23 | 5 % |
| wind speed in the four lowest cells vs u*/kappa ln((z + z0)/z0) | log law | 10 % |
| Kmv at the second cell over rho kappa u* (z + z0) | 1 | +/- 0.3 |


## Wall distance over a raised bottom

`check_walldist.py` checks one thing: on a flat mesh the wall distance is
the height above the bottom wall, z_cc - `geometry.prob_lo[2]`, whatever
the mesh and the box layout. Two CTest entries run the deck for one step
with the bottom at 300 m and `input_sounding_raised` (the same profile
moved up by 300 m):

| test | mesh | boxes |
| --- | --- | --- |
| `RANS_WallDist_RaisedBase` | constant dz | one box |
| `RANS_WallDist_RaisedBase_StretchedZSplit` | stretched, `erf.initial_dz = 4` | split in z (`amr.no_box_split_dir = -1`) |

The script stops if `geometry.prob_lo[2]` is 0, where height and wall
distance are the same number and the check could not fail. Before
erf-model/ERF#4609 and #4610 were fixed, both entries failed: the
constant-dz distance included `prob_lo[2]` (300 m too large), and the upper
box of the split stretched mesh read the surface node at k = 0 outside its
own data, an undefined value (300 m off in the build where this was measured).

## Shear stencil of the production of k

The production of k uses the vertical shear at the cell centre. The default
(`erf.rans_shear_stencil = central`) is the average of the shears on the two
faces of the cell, on uniform layers the central difference
(u_k+1 - u_k-1) / (2 dz). For a logarithmic profile that difference is larger
than du/dz at the centre: by 20 %, 6 % and 3 % in the cells 2, 3 and 4 of
this deck (erf-model/ERF#5261). With `erf.rans_shear_stencil = log_law` the
shear above the wall cell and below the top cell is

    (u_k+1 - u_k-1) / [(d_k + z0) ln((d_k+1 + z0) / (d_k-1 + z0))]

with d the wall distance and z0 the roughness length of the surface layer,
which is exact for the log profile over flat ground and less accurate than
`central` for a profile that is linear in height (Docs: theory/RANS, "Shear in
the production of k"). The stencil follows the `log_law` option of GUST
(erf-model/gust-abl#132).

`RANS_Neutral_ABL_Flat_LogShear` runs the deck for one step of 0.05 s from
`input_sounding_log`, u = (u_* / kappa) ln((z + z0) / z0) with u_* = 0.4 m/s
at the cell centres (written by `make_log_sounding.py`), once with each
stencil. `check_log_shear.py` reads the plotfiles of step 0 and step 1 of
both runs:

| check | target | tolerance |
| --- | --- | --- |
| all fields finite; cap off (0 < max cap_ratio < 1) | yes | exact |
| S of the closure, central run, vs the central stencil on the step-0 wind | equal | 1e-8 |
| S of the closure, log_law run, vs the log_law stencil on the step-0 wind | equal | 1e-8 |
| S of the log_law run vs u_* / (kappa (z + z0)) | equal | 1e-8 |
| S_central / S_log_law in the cells 2, 3, 4 | 1.2046, 1.0588, 1.0285 | 6e-5 |
| S in the wall cell and in the top cell, between the runs | equal | 1e-9 |
| Kmv between the runs | equal | 1e-12 |
| (rho k_log_law - rho k_central) / (dt Kmv) vs S_log_law^2 - S_central^2 | equal | 1 % |

`RANS_Neutral_ABL_Flat_LogShear_Stretched` repeats the pair with
`erf.initial_dz = 6` (top layer 5 times the bottom one) and boxes cut in z;
`check_log_shear.py --stretched` reads the cell heights from the plotfile and
leaves out the closed form and the pinned ratios, because the sounding is
interpolated to the levels. `RANS_Neutral_ABL_Flat_LogShear_Fitted` repeats it
on a terrain-fitted mesh over flat ground with all checks of the table (FFT
builds only).

S is recovered from the plotfile as cap_ratio C_r sqrt(k) / (cmu l); the runs
set `erf.rans_realizability = 10` so the closure writes cap_ratio without the
cap acting.

```bash
python3 make_log_sounding.py input_sounding_log
mpirun -np 1 erf_exec inputs_neutral erf.use_fft=false erf.input_sounding_file=input_sounding_log \
    max_step=1 erf.fixed_dt=0.05 erf.rans_realizability=10 erf.plot_int_1=1 erf.check_int=-1 \
    "erf.plot_vars_1=density x_velocity y_velocity KE Kmv Lturb diss cmu cap_ratio" \
    erf.rans_shear_stencil=central erf.plot_file_1=a_plt
# repeat with erf.rans_shear_stencil=log_law erf.plot_file_1=b_plt
python3 check_log_shear.py --dt 0.05 --cr 10 a_plt00001 b_plt00001
```

The 12 h run of the deck with each stencil (2 ranks, about a minute each):

| cell | z [m] | k / u_*^2, `central` | k / u_*^2, `log_law` | central / log_law shear |
| --- | --- | --- | --- | --- |
| 1 | 7.81 | 3.232 | 3.232 | wall cell, unchanged |
| 2 | 23.44 | 3.334 | 2.887 | 1.2046 |
| 3 | 39.06 | 3.050 | 2.829 | 1.0588 |
| 4 | 54.69 | 2.864 | 2.742 | 1.0285 |
| 5 | 70.31 | 2.720 | 2.645 | 1.0169 |
| 6 | 85.94 | 2.596 | 2.547 | 1.0112 |
| 7 | 101.56 | 2.484 | 2.449 | 1.0080 |
| 8 | 117.19 | 2.380 | 2.353 | 1.0060 |

u_* is 0.3931 m/s with `central` and 0.3897 m/s with `log_law`. With `central`
k rises from the wall cell to the second cell; with `log_law` it falls
monotonically. Both runs pass `check_neutral.py --physics`.

The converged wind is not the log profile at the cell centres, because the
momentum equations keep their central stress divergence. The shear the
production uses, times l / u_*, is

| cell | `central` | `log_law` |
| --- | --- | --- |
| 2 | 1.058 | 0.924 |
| 3 | 0.962 | 0.947 |
| 4 | 0.934 | 0.929 |
| 5 | 0.913 | 0.910 |

so `log_law` removes the excess in the second cell and leaves a deficit of
similar size. Two more numbers move with the option: the wind speed in cell 4
is 8.3 % above the log law of the run's u_* (6.5 % with `central`), and
Kmv(cell 1) / (kappa u_* (z_1 + z0)) is 0.714 (0.767 with `central`).
