# Canonical RANS results by phase

Numbers from the check scripts, recorded as each phase of `PLAN.md` lands.
"Before" is the code as branched from `upstream/development` (c001e148d).

## Neutral_ABL_Flat

12 h physics run, 2 ranks, `check_neutral.py --physics plt08640`.

| check | target | tol | before (phase 1) | after phase 2 |
| --- | --- | --- | --- | --- |
| all fields finite; KE, Kmv, diss >= 0 | yes | exact | pass | pass |
| max abs(walldist - z_cc) [m] | 0 | 1e-8 | 0 | 0 |
| Lturb(k=0,1,2) vs capped kappa (z + z0) | equal | 0.1 % | pass (2.93, 7.30, 10.46 m) | pass |
| max Lturb [m] | <= 30 | 1e-9 | 26.9 | 25.4 |
| diss vs AL01 Eq. 19, interior cells | equal | 5 % | 2.2 % | 2.1 % |
| u* [m/s] | 0.35 | +/- 0.15 | 0.393 | 0.393 |
| KE(k=0) / u*^2 | 3.23 | 5 % | **2.60 (fail, 20 % low)** | 3.232 |
| wall cell k_start / k_end | 1 | 1 % | **1.244 (fail)** | 1.000 |
| abs(U)(k=0..3) vs log law | equal | 10 % | -0.3 %, +1.7 %, +4.2 %, +6.7 % | -0.3 %, +1.5 %, +4.1 %, +6.5 % |
| Kmv(k=1) / (rho kappa u* (z + z0)) | 1 | +/- 0.3 | 0.76 | 0.77 |

Smoke run (40 steps): every structural check passes on 1 and 2 ranks, and
the two rank counts agree to 2e-15 in every planar-averaged field.

Reading of the two phase-1 failures: the Dirichlet wall value of k was
written into the first cell once per step, then the one-sided diffusive
flux against a zero ghost cell drained a fifth of it before the step ended,
so the converged wall k sat 20 % below the AL01 value.

Phase 2 keeps the logical BC for RhoKE at foextrap (ghost cell equals the
first cell, zero diffusive flux through the wall face via the surface-layer
branch), re-imposes the first-cell value after every RK stage, and pins the
bottom row of the implicit vertical diffusion solve. Both failures clear;
the interior profile is unchanged to the digits shown. A 40-step run with a
checkpoint at step 20 and a restart from it matches the straight run to
1e-14 in every field. The anelastic integrator disables the implicit
vertical solve, so the pinned row was exercised with a compressible variant
of the deck (dt 0.5 s, 20 acoustic substeps, 2 h): see the note below.

Compressible check of the implicit path (2 h, dt 0.5 s, 20 acoustic
substeps, 2 ranks): `vert_implicit_fac 1 1 0` with `tke = 1` in the banner,
wall cell k_start/k_end = 1.000 and KE(0)/u*^2 = 3.232 with the implicit
KE solve on and off; the two runs differ by 1e-9 in KE, so the pinned row
is exercised and holds the wall value.

## Phase 3: robustness

Neutral_ABL_Flat after phase 3 (12 h, 2 ranks): every physics check passes
with the same numbers as phase 2 to the digits in the table above. The
closure refactor onto `ERF_RANSClosure.H` was bit-identical to phase 2
before the unstable-length bound went in; with the bound, KE differs by
2e-6 and Lturb by 4e-4 m at most, in cells where dtheta/dz noise makes N^2
slightly negative.

| run | dt [s] | integrator | KE diffusion | dissipation | outcome |
| --- | --- | --- | --- | --- | --- |
| deck | 5 | anelastic | explicit (forced by anelastic) | explicit | all checks pass |
| deck | 5 | anelastic | explicit | implicit | all checks pass, same numbers |
| 4x dt | 20 | anelastic | explicit | explicit | all checks pass |
| 4x dt | 20 | anelastic | explicit | implicit | **abort at 10.6 h, negative theta at 290 m** |
| large dt | 60 | compressible, 2400 substeps | implicit | explicit | all checks pass |
| large dt | 60 | compressible, 2400 substeps | implicit | implicit | all checks pass |

Reading: the anelastic integrator switches every vertical diffusion to
explicit, and dz^2 / (2 K/rho) is about 19 s mid-layer for this deck, so
dt = 20 s sits on the explicit diffusion limit. The explicit-dissipation
run survived it by a small margin (dissipation damps k and so K); the
implicit-dissipation run kept slightly more k and crossed it. Neither
result is about the dissipation itself: the compressible pair at dt = 60 s,
where vertical diffusion is implicit, passes both ways. In this deck the
dissipation time scale is never the binding one because the wall cell is
held by the Dirichlet condition and the second cell's scale is about a
minute. The option stays opt-in and verified equivalent; the practical
limit for RANS under the anelastic integrator is the explicit vertical
diffusion, which is a dycore matter outside this plan.

Unit tests (`erf_unit_tests --gtest_filter=RANSClosure*`, 7 tests): the
first version caught the Burchard & Petersen smoothing returning -2 for
Rt = -1e16 and +1.4e14 for Rt = -1e30 through cancellation; the
rearranged form `Rt_crit + a x / (x + a)` matches the original to
round-off up to |Rt| of 1e8 and holds Rt_min beyond.

Input validation: `erf.Rt_min = -4` now aborts with the pole message;
`erf.tke_floor = 1e-4` runs.

