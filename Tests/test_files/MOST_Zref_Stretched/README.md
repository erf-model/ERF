# MOST_Zref_Stretched

MOST reference height on a flat periodic column with a stretched vertical
grid: 4x4x40 cells, `erf.initial_dz = 10`, `erf.grid_stretching_ratio = 1.1`
(top 4425.9 m), `zlo.type = surface_layer`, `erf.most.z0 = 0.1`, MRF,
anelastic, `erf.most.zref` unset.

The initial wind is sheared on purpose. `sounding_most_zref` ramps it linearly
from 10 m/s at the ground to 15 m/s at 50 m, and holds it constant above. ERF
interpolates the sounding linearly to the cell centers, so every cell below
50 m holds the exact ramp value at its height. So does any interpolation that
is linear in z. At t = 0 the plane-averaged u* must therefore equal
`0.41 * U(zref) / ln(zref / 0.1)`, with the wind taken at the reference height
itself. If a lookup reports one height but takes the wind from another cell,
u* comes out different. From u*(0) and the reported zref, the script inverts
the ramp and prints the height the wind actually came from. With a uniform
wind this check could not tell the cells apart.

`run_most_zref.py` runs the deck through the three MOST lookups a flat
stretched mesh can take, plus a uniform 10 m column (443 levels):

| case | mesh / terrain | lookup | zref |
|---|---|---|---|
| `fitted_plane` | StretchedDz, fitted (set by the stretching) | `set_k_indices_T` | 15.5 m |
| `fitted_interp` | VariableDz, `erf.terrain_type = StaticFittedMesh` | `set_norm_positions_T` + `trilinear_interp_T` | 10 m |
| `immersed_stretched` | StretchedDz, `erf.terrain_type = ImmersedForcing` | `set_k_indices_N` | 5 m |
| `uniform_reference` | ConstantDz, dz = 10 m | `set_k_indices_N` | 5 m |

The test passes when:
- every case reports the expected height;
- u*(0) matches the log law with the wind at that height to 2e-5 (hist.dat
  holds 6 digits);
- the stretched no-terrain column stays within 1e-3 of the uniform column
  over the 10 steps.

The script also aborts if the sounding is ever made too uniform to tell a
cell from its upper neighbour.

## Old and new

Measured on 2026-09-11 against development 48f0fe0cd (`--old-exe`):

| case | old | new zref | new u*(0) | log law | wind taken at |
|---|---|---|---|---|---|
| `fitted_plane` | abort, `ERF_MOSTAverage.cpp:768` "zref not found with terrain!" | 15.5 | 0.938945 | 0.938945 | 15.50 m |
| `fitted_interp` | abort, `ERF_MOSTAverage.H:371` "Height above terrain not found" | 10 | 0.979334 | 0.979334 | 10.00 m |
| `immersed_stretched` | zref 55.32, u*(0) 0.681625, wind taken at 5.00 m | 5 | 1.10045 | 1.10045 | 5.00 m |
| `uniform_reference` | zref 5, u*(0) 1.10045 | 5 | 1.10045 | 1.10045 | 5.00 m |

- **Fitted cases:** the default 10 m query lies exactly on the top face of the
  first cell. The height search accepted only `z_lo < z < z_hi`, so it found no
  cell. The search is now `z_lo <= z < z_hi`: a face belongs to the cell above
  it.
- **Interpolation weight:** `trilinear_interp_T` used to weight the two
  bracketing cells by the fraction of the containing cell. That is exact only
  for equal cell heights. On this column, at the 10 m face it returned the plain
  mean of the 5 m and 15.5 m cells, which is the wind at 10.25 m (u*(0)
  0.98156, 2.3e-3 high). It now weights by the physical heights of the two cell
  centers. On equal cell heights the new weight is the same as the old one, up
  to round-off.
- **No-terrain stretched case:** the default height was half of
  `(prob_hi - prob_lo)/nz` (110.6 m), but the wind still came from the first
  cell at 5 m. The cell centers and indices now come from the staggered z
  levels, and the given-zref and `k_arr_in` branches use them too.

**Checking the checker.** A build whose k lookups deliberately return the cell
above the right one fails every index case:
- `fitted_plane`: wind taken at 27.05 m;
- `immersed_stretched`: wind taken at 15.50 m;
- `uniform_reference`: wind taken at 15.00 m.

Each is about 10% off in u*(0).

Stretched against uniform, new binary:

| t [s] | stretched u* | uniform u* | rel diff |
|---|---|---|---|
| 0 | 1.10045 | 1.10045 | 0 |
| 5 | 1.0636 | 1.06323 | 3.5e-04 |
| 10 | 1.03113 | 1.03083 | 2.9e-04 |

## Running by hand

```
python3 run_most_zref.py --exe /path/to/erf_exec --mpi-cmd "mpiexec -n 1" \
    [--old-exe /path/to/erf_exec_before_fix]
```
