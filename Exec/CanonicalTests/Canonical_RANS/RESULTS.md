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
