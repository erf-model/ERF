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
plotfile with `../erf_plotfile.py`, averages every field over x and y, and
prints one row per check with the measured value, the target and the
tolerance. The exit code is non-zero if any enabled check fails.

## Checks

Smoke (`--smoke`, the CTest entry, 40 steps):

| check | target | tolerance |
| --- | --- | --- |
| all planar-averaged fields finite | yes | exact |
| min KE, Kmv, diss | >= 0 | exact |
| wall distance vs cell-centre height | equal | 1e-8 m |
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

Results for the current code are recorded in `RESULTS.md` at the top of
`Canonical_RANS/` as the phases of the RANS plan land.
