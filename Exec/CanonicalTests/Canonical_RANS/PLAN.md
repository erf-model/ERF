# Canonical RANS: a minimal working one-equation TKE model

Branch `claude-RANS` off `upstream/development`, worktree `.worktrees/claude-RANS`.
Goal: the Axell & Liungman (2001) k-equation closure (`erf.rans_type = kEqn`) runs
correctly on flat and terrain-fitted meshes under a MOST surface layer, with
regression decks in `Exec/CanonicalTests/Canonical_RANS/` laid out like
`Canonical_LES/` (one folder per case: README, input_sounding, inputs decks,
planar-average data logs, a check script) and CTest entries. Fire coupling
comes afterwards in ERF-Hazard.

Out of scope (agreed 2026-09-09): moist buoyancy, EB and thin-body wall
distance, wall distance after a regrid, hybrid RANS-LES, two-equation closures.

Conventions: every behaviour change is opt-in with the old path as default
unless stated in the phase; each phase ends with its deck(s) passing and a
commit; decks are sized for 40-step CTest smoke runs plus a longer
"physics" run documented in the README.

Regression-test rule (Harish, 2026-09-09): every deck ships with a Python
check script (`check_<case>.py`) that reads the run output (plotfile via yt
or the `data_log` profiles) and compares numbers against stated targets with
tolerances, exiting non-zero on failure. The CTest entry runs the deck and
then the script, so a run that finishes but drifts fails the test. Each
script prints the measured value, the target and the tolerance for every
check, and the README lists the same table. A clean exit plus a plotfile is
never the pass criterion on its own.

## Phase 1: baseline harness (flat, neutral)

- Configure `cmake-build/` (Release, MPI, all warnings, no fire/dust).
- `Canonical_RANS/Neutral_ABL_Flat/`: ABL problem, 2.5 km x 2.5 km x 1 km,
  coarse RANS grid (e.g. 16x16x64 with stretched dz), geostrophic 10 m/s,
  MOST with z0, `rans_type = kEqn`, `init_tke_from_ustar`, `data_log`
  profiles, plot vars `KE Kmv Khv Lturb walldist diss`.
- `check_neutral.py`: reads `mean_profiles.dat` and `sfs_profiles.dat`;
  reports log-law error in the surface layer, k/u*^2 at the first cell,
  min k, max Lturb.
- `add_test_rans` in `Tests/CTestList.cmake` (40 steps, labels
  `rans;regression`); register the deck.
- Record the unfixed behaviour (this is the "before" column for the PR).

Exit: deck runs 40 steps cleanly on 1 and 2 ranks; check script runs.

Status (2026-09-09): done. `RANS_Neutral_ABL_Flat` passes the smoke checks on
1 and 2 ranks (rank counts agree to 2e-15); the 12 h run gives the "before"
column in `RESULTS.md`: log law within 7 %, but the wall cell keeps only
80 % of its Dirichlet k through a step and KE(0)/u*^2 = 2.60 against 3.23.
Two harness lessons: `surf_hist.dat` needs `erf.v > 0`, and the plotfile
holds `diss` from the start of the step and `KE` from its end, so the
Eq. 19 check is restricted to interior cells with a 5 % allowance. The
wall-cell `k_start/k_end` check moves from `--physics` to `--smoke` in
phase 2 once it passes.

## Phase 2: wall boundary for k

- Make `dirichlet_k` self-consistent: the ghost cell below the wall gets the
  same AL01 Eq. 16 value as the first cell, the diffusive KE flux through the
  wall face is zero, and the first-cell value is re-imposed every RK stage
  (zero its RHS) rather than once per step on `S_old`.
- Keep implicit KE diffusion available with Dirichlet k (extend the BC
  sanity check instead of silently switching it off).
- Neutral (no buoyancy) branch of Eq. 16 when `t_star` is positive stays.
- Decide default: `dirichlet_k = true` when `rans_type = kEqn` and
  `zlo.type = surface_layer` (proposed; flagged for Harish).

Exit: neutral deck reaches a log-law profile within a few percent between
2 z0 and 0.1 zi after the physics run; k at the first cell equals
u*^2 / Cmu0^2 to round-off; 40-step smoke deck bit-identical on restart.

Status (2026-09-09): done. Implementation: `init_bcs` no longer switches
RhoKE to ext_dir at a surface-layer wall (foextrap keeps the ghost cell
equal to the first cell and the surface-layer branch of the diffusion gives
a zero wall flux); `erf_slow_rhs_post` restores the first-cell value from
`S_old` after every stage; `ImplicitDiffForStateLU_{N,S,T}` pin the bottom
row for RhoKE. KE(0)/u*^2 = 3.232 (target 3.2325), k_start/k_end = 1.000,
restart bit-exact, interior profile unchanged. Harish's decision
(2026-09-09): keep the default false, warn at startup when kEqn runs under a
surface layer without it (near-wall k settles at half the AL01 value, the
mean wind still follows the log law), and set the flag in every deck that
enables kEqn (Canonical_RANS, Askervein; no Tests/test_files deck uses kEqn).

## Phase 3: robustness and hygiene

- Validate `Rt_crit`, `Rt_min`, `Cmu0`, `Cb`, `max_geom_lscale` in
  `TurbChoice::init_params`; guard the cmu' and smoothing denominators.
- `erf.tke_floor` (default 0, meaning the historical machine-epsilon floor)
  used in the viscosity and the slow RHS. Not `tke_min`: that existing key
  is the initial TKE of the prognostic closures.
- Opt-in `erf.implicit_tke_dissipation`: dissipation linearised as
  (Cmu0^3 sqrt(k) / l) * k and folded into the update.
- Remove the dead grown-box loop in `ComputeTurbulentViscosityRANS`, drop
  unused parameters, fix the `hfx_z(i,j,-1)` comment, use `geom[lev]` for
  periodicity in the wall-distance solve.
- gtest for the closure kernel: cmu(0) = Cmu0, cmu' monotone in Rt, smoothing
  continuous at Rt_crit, length scale limits (neutral, strongly stable).

Exit: unit tests pass; neutral deck at 4x the CFL-limited dt stays bounded
with the implicit dissipation on.

Status (2026-09-09): done, with the exit criterion corrected. Landed:
closure relations factored into `Source/Diffusion/ERF_RANSClosure.H`
(namespace `AL01`, bit-identical refactor); input validation for Cmu0, Cb,
max_geom_lscale, Rt_crit <= 0, Rt_min < Rt_crit, Rt_min > -3.6 (poles of
Eqs. 31-32); `erf.tke_floor`; opt-in `erf.implicit_tke_dissipation`
(source skips the sink, the update divides by 1 + dt c with
c = diss_old / (rho k)_old, half-weighted on the anelastic stage 1);
component fill on the tilebox only; unused parameters dropped; comment
and per-level periodicity fixes; 7 gtests. The unit test caught the
smoothing formula cancelling for |Rt| > 1e15; it is now
`Rt_crit + a x / (x + a)`. Two limiters from Harish's Kynema KLAxell were
reviewed against the paper: the unstable length now evaluates Eq. 28 once
with the smoothed Rt from the geometric length (bounded by about
1.31 l_g), replacing the two-pass corrector, which iterated the fixed-point
map of Eq. 26 that has no fixed point in strong convection (AL01 p. 78);
the extra stable cap sqrt(Cmu k / N^2) is never active because Cb = 0.35 is
below sqrt(Cmu0) = 0.75, so it was not added. Kynema's sigma_k = 0.5
differs from Table I (1.0); ERF keeps 1.0. The 4x-dt criterion turned out
to test the anelastic integrator's explicit vertical diffusion, not the
dissipation (see RESULTS.md); the dissipation option was verified with a
compressible dt = 60 s pair instead. Harish (2026-09-09): anelastic with
FFT is the main use, compressible stays optional.

## Phase 4: stratified flat cases

- Consistent diffusivities: Theta_h follows Theta_v (cmu'/cmu), scalars and
  moisture use the same ratio, behind `erf.rans_consistent_diffusivities`.
- Opt-in `erf.rans_lscale_from_pblh`: cap l_g at kappa * 0.1 * zi using the
  surface-layer PBL height diagnostic instead of the fixed 30 m.
- `Stable_ABL_Flat/` (GABLS1-style surface cooling, 9 h physics run) and
  `Convective_ABL_Flat/` (specified surface heat flux, 4 h) decks with
  check scripts: BL depth, inversion strength, Rt within [Rt_min, ...],
  k > 0 everywhere, well-mixed theta in the convective case.

Exit: both decks pass their checks; 40-step smoke entries registered.

## Phase 5: terrain wall distance

- Verify the Poisson wall distance on a terrain-fitted mesh: with a flat
  `StaticFittedMesh` it must reproduce the analytic distance; on a gentle
  hill it must match the vertical height above ground to a few percent.
- Confirm ghost-cell fill of `walldist` on terrain, the z0 offset, and the
  h_zeta scaling of dtheta/dz in the RANS kernel.
- `Neutral_Hill_2D/`: Witch of Agnesi (from `Idealized_Terrain/`) under
  MOST with `rans_type = kEqn`, 2D (ny = 1) with the nodal Poisson solve
  (check that the 2D case is handled; the solver cannot hide a direction).
- check script: walldist error map, hill-top speed-up, no negative k.

Exit: walldist error under the tolerance on both meshes; hill deck runs
40 steps on 1 and 2 ranks and the physics run gives a smooth speed-up.

## Phase 6: 3D terrain canonical case and restart

- `Neutral_Hill_3D/`: Gaussian hill (h/L about 0.2), 3D, neutral, MOST,
  geostrophic forcing; compare hill-top fractional speed-up with the
  Jackson-Hunt estimate (about 2 h/L) and check lee-side k.
- Restart test on terrain: checkpoint mid-run, restart, compare plotfiles
  (wall distance is recomputed in `InitData_post`, Dirichlet k restored).
- Askervein deck (`Real_Terrain/Askervein/inputs_anel`) updated to the new
  input names and run once as the "real" sanity case (not a CTest).

Exit: speed-up within the expected band; restart bit-identical.

## Phase 7: documentation and diagnostics

- Theory: new `Docs/sphinx_doc/theory/RANS.rst` (equations 15-32 of AL01,
  the length-scale limits, the Poisson wall distance, the wall condition).
- `Inputs.rst`: every new key; `RegressionTests.rst`: the Canonical_RANS
  entries; a README per case in the Canonical_LES style.
- Optional plot variables `Rt`, `cmu`, `cmu_prime` (stored in the eddy
  diffusivity container's spare components or a small diagnostic MultiFab).

Exit: docs build; every input key documented.

## Phase 8: sweep and PR

- Run all Canonical_RANS decks and `ctest -L rans`; single-precision build;
  warning-clean build; codespell.
- PR to `erf-model/ERF` `development` from `hgopalan:claude-RANS` with the
  before/after table from Phase 1, limitations (moist, EB, regrid, hybrid)
  and effort estimates for each.
- Note for ERF-Hazard: fire couples through the surface heat flux at k = 0
  (already the buoyancy source) and, if wanted, a volumetric Q_fire term in
  the buoyancy production; both are follow-ups after the merge.

Exit: PR open; memory note updated.
