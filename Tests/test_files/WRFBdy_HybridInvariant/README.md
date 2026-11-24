# WRFBdy_HybridInvariant

## Purpose

real.exe writes the wrfbdy arrays as the WRF state coupled with the dry column mass and the
map factors (`couple()` in `dyn_em/module_big_step_utilities_em.F`):

| Array | Coupling |
| --- | --- |
| `U_B*` | `u (c1h mu_u + c2h) / msfuy` |
| `V_B*` | `v (c1h mu_v + c2h) / msfvx` |
| `T_B*`, `QVAPOR_B*` | `(theta_m - 300) (c1h mu + c2h)`, `qv (c1h mu + c2h)` |
| `PH_B*` | `ph' (c1f mu + c2f)` (full levels) |

At the first boundary time these are the wrfinput state, so once ERF has decoupled them and
interpolated them onto its own heights they must give back the initial state that
`init_from_wrfinput` builds from the same columns. This case checks that, on each side, for
u, v, theta, qv and the density.

It covers issue #5149 (PH was divided by `mu`, which is right only for `hybrid_opt = 0`) and
the map factors of the boundary winds (ERF's velocities are the physical ones, so the
boundary has to undo WRF's `1/msfuy` and `1/msfvx`).

`erf.check_wrfbdy_invariant = true` makes ERF print, at the first boundary time, the
largest difference on each side, over the points whose ERF height lies inside the WRF
column (outside it the two remaps fill differently):

    wrfbdy invariant xlo u: max |bdy - init| = 2.76e-06 (rel 1.10e-07)

## Files

| File | Role |
| --- | --- |
| `make_wrf_hybrid.py` | Writes `wrfinput_d01` and `wrfbdy_d01`: 16 x 12 x 20 cells, `hybrid_opt = 2` with `eta_c = 0.2`, a 1200 m hill next to the west boundary, a column-mass perturbation of up to 400 Pa, map factors from 1 to about 1.05, sheared physical winds, `USE_THETA_M = 1`, water vapor, and a geopotential integrated hydrostatically from `ALT`; three-cell boundary zone, two times. `--no-cf` leaves out C1F/C2F, `--hybrid0` uses the sigma coordinate, `--bdy-without-msf` couples the winds without the map factors |
| `WRFBdy_HybridInvariant.i` | One level in four boxes, `max_step = 0`, `erf.real_width = 3`, `erf.avg_grid_faces_to_nodes = true` and no rebalance (so the initial density is `1/ALT`), Kessler moisture, the check on |
| `check_invariant.py` | Reads the invariant lines of the log and compares them with the tolerances; `--fail FIELD` requires that field out of tolerance (a control), and the fields not named are then reported but not judged |

## Tests

| Test | Input | Expected |
| --- | --- | --- |
| `WRFBdy_HybridInvariant` | the pair as made | all fields within tolerance (u, v to about 1e-7) |
| `WRFBdy_HybridInvariant_NoCF` | no C1F/C2F in wrfinput | u, v and the density out of tolerance (about 3e-3 to 6e-3) under the c1f = 1, c2f = 0 fallback |
| `WRFBdy_HybridInvariant_Sigma_NoCF` | sigma coordinate, no C1F/C2F | within tolerance: the fallback is exact for `hybrid_opt = 0` |
| `WRFBdy_HybridInvariant_NoMapFac` | winds coupled without the map factors | u and v out of tolerance (about 5%) |
| `WRFBdy_HybridInvariant_MPI` | the pair as made, four ranks | as the first |

## Tolerances

u and v: 1e-5 (they agree to the float precision of the files). theta 1e-4, qv 5e-4,
density 1e-4: the initial state is remapped as rho theta, rho qv and rho, the boundary data
as theta, qv and rho, which differ at second order in the layer thickness.

## Data

None is downloaded or redistributed; `make_wrf_hybrid.py` writes both files.
