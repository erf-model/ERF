# Canonical RANS results by phase

Numbers from the check scripts, recorded as each phase of `PLAN.md` lands.
"Before" is the code as branched from `upstream/development` (c001e148d).

## Neutral_ABL_Flat

12 h physics run, 2 ranks, `check_neutral.py --physics plt08640`.

| check | target | tol | before (phase 1) |
| --- | --- | --- | --- |
| all fields finite; KE, Kmv, diss >= 0 | yes | exact | pass |
| max abs(walldist - z_cc) [m] | 0 | 1e-8 | 0 |
| Lturb(k=0,1,2) vs capped kappa (z + z0) | equal | 0.1 % | pass (2.93, 7.30, 10.46 m) |
| max Lturb [m] | <= 30 | 1e-9 | 26.9 |
| diss vs AL01 Eq. 19, interior cells | equal | 5 % | 2.2 % |
| u* [m/s] | 0.35 | +/- 0.15 | 0.393 |
| KE(k=0) / u*^2 | 3.23 | 5 % | **2.60 (fail, 20 % low)** |
| wall cell k_start / k_end | 1 | 5 % | **1.244 (fail)** |
| abs(U)(k=0..3) vs log law | equal | 10 % | -0.3 %, +1.7 %, +4.2 %, +6.7 % |
| Kmv(k=1) / (rho kappa u* (z + z0)) | 1 | +/- 0.3 | 0.76 |

Smoke run (40 steps): every structural check passes on 1 and 2 ranks, and
the two rank counts agree to 2e-15 in every planar-averaged field.

Reading of the two failures: the Dirichlet wall value of k is written into
the first cell once per step, then the one-sided diffusive flux against a
zero ghost cell drains a fifth of it before the step ends. The converged
wall k therefore sits 20 % below the AL01 value. Phase 2 addresses this.
