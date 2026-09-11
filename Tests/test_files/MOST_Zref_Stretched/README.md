# MOST_Zref_Stretched

MOST reference height on a flat periodic column with a stretched vertical
grid: 4x4x40 cells, `erf.initial_dz = 10`, `erf.grid_stretching_ratio = 1.1`
(top 4425.9 m), `zlo.type = surface_layer`, `erf.most.z0 = 0.1`, MRF,
anelastic, `erf.most.zref` unset. The initial wind is a uniform 15 m/s, so at
t = 0 the plane-averaged u* must equal the log law
`0.41 * 15 / ln(zref / 0.1)` at the height the surface layer uses.

`run_most_zref.py` runs the deck through the three MOST lookups a flat
stretched mesh can take, plus a uniform 10 m column (443 levels):

| case | mesh / terrain | lookup | zref |
|---|---|---|---|
| `fitted_plane` | StretchedDz, fitted (set by the stretching) | `set_k_indices_T` | 15.5 m |
| `fitted_interp` | VariableDz, `erf.terrain_type = StaticFittedMesh` | `set_norm_positions_T` + interpolation | 10 m |
| `immersed_stretched` | StretchedDz, `erf.terrain_type = ImmersedForcing` | `set_k_indices_N` | 5 m |
| `uniform_reference` | ConstantDz, dz = 10 m | `set_k_indices_N` | 5 m |

It passes when every case reports the expected height, u*(0) matches the log
law to 2e-5 (hist.dat holds 6 digits), and the stretched no-terrain column
stays within 1e-3 of the uniform column over the 10 steps.

## Old and new

Measured on 2026-09-11 against development 48f0fe0cd (`--old-exe`):

| case | old | new zref | new u*(0) | log law |
|---|---|---|---|---|
| `fitted_plane` | abort, `ERF_MOSTAverage.cpp:768` "zref not found with terrain!" | 15.5 | 1.21941 | 1.21941 |
| `fitted_interp` | abort, `ERF_MOSTAverage.H:371` "Height above terrain not found" | 10 | 1.33546 | 1.33546 |
| `immersed_stretched` | zref 55.32, u*(0) 0.973749 | 5 | 1.57208 | 1.57208 |
| `uniform_reference` | zref 5, u*(0) 1.57208 | 5 | 1.57208 | 1.57208 |

- **Fitted cases:** the default 10 m query lies exactly on the top face of
  the first cell. The height search accepted only `z_lo < z < z_hi`, so it
  found no cell. The search is now `z_lo <= z < z_hi`: a face belongs to the
  cell above it, and the interpolation is continuous across the face.
  `erf.most.zref = 10.01` was the workaround and gives the same answer as the
  new default.
- **No-terrain stretched case:** the default height was half of
  `(prob_hi - prob_lo)/nz` (110.6 m), paired with the first-cell wind. The
  cell centers and indices now come from the staggered z levels. The given-zref
  and `k_arr_in` branches use them too.

Stretched against uniform, new binary:

| t [s] | stretched u* | uniform u* | rel diff |
|---|---|---|---|
| 0 | 1.57208 | 1.57208 | 0 |
| 5 | 1.47586 | 1.47594 | 5.4e-05 |
| 10 | 1.37895 | 1.37967 | 5.2e-04 |

## Running by hand

```
python3 run_most_zref.py --exe /path/to/erf_exec --mpi-cmd "mpiexec -n 1" \
    [--old-exe /path/to/erf_exec_before_fix]
```
