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

## Stable_ABL_Flat (phase 4)

9 h GABLS1 run, 2 ranks, `check_stable.py --physics plt16200`.

| check | target | measured |
| --- | --- | --- |
| u* [m/s] | 0.20 to 0.35 | 0.244 |
| theta(k=0) minus imposed surface theta [K] | 0 to 1.5 | 0.25 |
| min dtheta/dz below 200 m [K/m] | >= 0 | 0.011 |
| max wind over Ug (low-level jet) | >= 1.02 | 1.23 |
| height of the wind maximum [m] | 50 to 300 | 154 |
| BL depth from KE [m] | 80 to 300 | 134 |
| KE(k=0)/u*^2 | 3.23 within 5 % | 3.232 |
| Lturb over neutral length, 120 to 300 m | <= 1 | 0.12 |

GABLS1 LES ensemble for reference: u* 0.26 to 0.30, jet near 150 to 200 m,
depth 150 to 200 m.

## Convective_ABL_Flat (phase 4)

4 h run at dt = 2 s, 2 ranks, `check_convective.py --physics plt07200`.

| check | target | measured |
| --- | --- | --- |
| u* [m/s] | 0.30 to 0.80 | 0.485 |
| column heat gain over rho_sfc F t | 1 within 10 % | 0.9998 |
| inversion height [m] | 900 to 1250 | 1020 |
| theta spread in 0.2 to 0.7 zi [K] | <= 2 (local closure) | 1.12 |
| max dtheta/dz in 0.2 to 0.7 zi [K/m] | <= 0 | -0.0012 |
| mixed-layer warming over F t / zi | 1 within 30 % | 1.075 |
| KE(k=0)/u*^2 | >= 3.23 (buoyancy adds) | 3.82 |
| min KE in 0.1 to 0.8 zi [m2/s2] | >= 0.05 | 0.55 |
| max KE above 1.3 zi [m2/s2] | <= 0.05 | 2e-16 |
| max Lturb over bound | <= 1 | 0.79 |

At dt = 5 s the run aborts after 1.7 h with a negative theta at 290 m:
Kmv reaches 42 kg/m/s, i.e. K/rho at the explicit limit dz^2/(2 dt) = 40
m2/s of the anelastic integrator. The profile at 4 h keeps a
superadiabatic lapse of -1 to -3 K/km through the mixed layer (K_h 40 to
67 kg/m/s, Lturb 24 to 45 m under the PBL-height cap), the expected
behaviour of a local-K closure without countergradient transport.

## Neutral_Hill_2D (phase 5)

6 h run at dt = 1.5 s, 2 ranks, `check_hill.py --physics plt14400`.

| check | target | measured |
| --- | --- | --- |
| walldist vs exact ridge distance, max relative | <= 10 % | 5.2 % |
| walldist vs exact ridge distance, mean relative | <= 3 % | 1.0 % |
| walldist abs error within 100 m of the surface [m] | <= 3 | 2.36 |
| crest speed-up, k = 0, 1, 2 | 0.2 to 0.8 (2 h/L = 0.4) | 0.55, 0.45, 0.39 |
| speed-up positive in the lowest ten cells | > 0 | min 1.43 m/s |
| upstream u* from the wall k [m/s] | 0.25 to 0.55 | 0.346 |
| upstream wind vs log law, k = 0, 1, 2 | within 15 % | 1.5 %, 1.7 %, 3.7 % |
| max Lturb over bound | <= 1 | 0.76 |
| wall-cell k retention along the surface | within 1 % | 3e-16 |

Flat-fitted variant (prob.hmax = 1e-6, 40 steps): before the gradient fix
the Poisson distance was z (1 - dz/2H), 0.78 % short at every height;
after it, exact to 1e-6 m above the first cell and 1.5 cm long in the
first cell. The ridge mean error went from 1.5 % to 1.0 %; the speed-ups
did not change to three digits because the length is capped at 30 m.

## Neutral_Hill_3D (phase 6)

4 h run at dt = 1.5 s, 64 x 64 x 20 at 40 m, 2 ranks, `check_hill3d.py --physics plt09600`.

| check | target | measured |
| --- | --- | --- |
| walldist vs exact hill distance, max relative | <= 15 % | 7.4 % (terrain_height), 9.5 % (poisson) |
| walldist vs exact hill distance, mean relative | <= 3 % | 0.013 % (terrain_height), 0.31 % (poisson) |
| walldist abs error within 100 m of the surface [cells] | <= 0.2 | 0.03 (terrain_height), 0.12 (poisson) |
| crest speed-up, k = 0, 1, 2 | 0.16 to 0.64 (1.6 h/L = 0.32) | 0.41, 0.28, 0.22 |
| speed-up positive in the lowest eight cells | > 0 | min 0.71 m/s |
| upstream u* from the wall k [m/s] | 0.25 to 0.55 | 0.318 |
| upstream wind vs log law, k = 0, 1, 2 | within 15 % | -5.7 %, +4.0 %, +10.4 % |
| wall-cell k retention over the surface | within 1 % | 3e-16 |

Wall-distance methods on the 2D ridge (40-step run): terrain_height mean
0.02 %, max 2.3 %, 0.01 cells near the surface; poisson mean 1.0 %, max
5.2 %, 0.15 cells. On the flat fitted mesh terrain_height is exact to
2e-10, poisson 0.2 % in the first cell and 1e-6 above it.

Restart (2D and 3D terrain decks, checkpoint at 20, compare at 40):
every field identical to 1e-14.

Askervein (20 steps, 4 ranks, 26 s): walldist 7.9 to 707 m, KE up to
6.9 m2/s2, Kmv up to 9.4 kg/m/s, all finite.

Mesh finding (not a RANS matter, see PLAN phase 6): any 3D fitted mesh
here with dz != dx yields a deterministic pre-projection divergence of
1.788e139 and aborts, independent of the closure and not a memory read
(no trap under `amrex.init_snan` at initialisation); a separate
uninitialised read in the w boundary fill trips the trap in the first
advance.

