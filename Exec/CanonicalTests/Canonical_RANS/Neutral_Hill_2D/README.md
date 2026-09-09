# Canonical RANS: neutral flow over a 2D ridge (Witch of Agnesi)

Neutral, periodic flow over a two-dimensional Witch of Agnesi ridge,
hmax = 100 m and L = 500 m (half-width at half height), on a
terrain-fitted mesh under a MOST surface layer, driven by a 10 m/s
geostrophic wind at f = 1e-4 1/s, z0 = 0.1 m. Closed with the one-equation
k RANS model of Axell & Liungman (2001). The wall distance the closure
needs on a fitted mesh comes from the Poisson solve of Tucker (2003,
J. Comput. Phys. 190, 229-248); this case checks it against the exact
distance to the ridge, and the flat-fitted variant (`prob.hmax = 1e-6`)
checks it against the analytic height above the surface on the same mesh
machinery.

| item | value |
| --- | --- |
| domain | 5120 x 40 x 1000 m, one cell wide in y |
| grid | 128 x 1 x 64, dx = 40 m, dz = 15.6 m, basic terrain following |
| time step | 1.5 s fixed (Courant 0.375), anelastic with the MLMG projection |
| closure | `erf.rans_type = kEqn`, AL01 defaults, `dirichlet_k = true` |
| physics run | 6 h (14400 steps), about 10 min on 2 ranks |
| smoke runs | 40 steps (`ctest -R RANS_Neutral_Hill_2D`, `ctest -R RANS_Flat_Fitted_2D`) |

`amr.blocking_factor = 1` because the domain is one cell wide in y. The
wall distance is `erf.wall_dist_type = terrain_height`, the height above
the local surface projected on its normal (exact to 1e-10 on a flat
mesh, 0.02 % mean error on this ridge, no linear solve); the Poisson
distance of Tucker (2003) is exercised by the `_Poisson` CTest variants
of this deck and its flat-fitted version.

## Running

```bash
mpirun -np 2 erf_exec inputs_hill
python3 check_hill.py --physics plt14400
mpirun -np 2 erf_exec inputs_hill prob.hmax=1e-6 max_step=40 erf.plot_int_1=40
python3 check_flat_fitted.py plt00040
```

The scripts read the full 3D fields (`../erf_plotfile.py`,
`read_fields`), so the planar-average logs, which ERF does not write on a
fitted mesh, are not needed.

## Checks

Wall distance (both modes; the field is set at initialisation):

| check | target | tolerance |
| --- | --- | --- |
| ridge: max relative error vs the exact distance to the curve | 0 | 15 % (crest cells, Tucker's convex-corner error) |
| ridge: mean relative error | 0 | 3 % |
| ridge: max absolute error where the distance is under 100 m | 0 | 0.2 dz |
| flat fitted: max relative error vs z - h | 0 | 3e-3 (first cell, see below) |
| flat fitted: max absolute error where the distance is under 100 m | 0 | 0.1 m |

The absolute error is judged below 100 m because the geometric length is
capped at 30 m, so the distance shapes the closure only near the surface.
On the flat fitted mesh the distance is exact to 1e-6 m above the first
cell; the first cell is 1.5 cm long because the Poisson solve represents
the wall Dirichlet value by an odd-reflection ghost, which shifts phi by
dz^2/8.

Structural (both modes): finite fields, KE and Kmv non-negative, Lturb at
most the unstable bound of the geometric length built on the wall
distance, wall-cell k retention within 1 % along the whole surface.

Physics (`--physics`, 6 h):

| check | target | tolerance |
| --- | --- | --- |
| crest speed-up U_top/U_up - 1 in the three lowest cells | 2 h/L = 0.4 (Jackson & Hunt 1975) | 0.2 to 0.8 |
| speed-up positive in the lowest ten cells | > 0 | exact |
| upstream u* from the wall k | 0.25 to 0.55 m/s | range |
| upstream wind in the three lowest cells vs the log law | equal | 15 % |

The upstream reference column sits 2 km, four half-widths, ahead of the
crest; with periodic hills every 5.1 km it still carries a little of the
previous ridge's wake, which the 15 % log-law tolerance allows for.
