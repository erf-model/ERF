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
| grid | 8 x 8 x 64, dz = 15.625 m, first cell centre 7.8 m |
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
| wall distance vs cell-centre height above `geometry.prob_lo[2]` | equal | 1e-8 m |
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
