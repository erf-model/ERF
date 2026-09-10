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

Reference-implementation rule (Harish, 2026-09-09): every change to the
closure or its TKE sources is cross-checked against the Kynema KLAxell and
KransAxell code (github.com/kynema/kynema-sgf, `src/turbulence/RANS/` and
`src/equation_systems/tke/source_terms/`), and the phase status records what
matches and what differs.

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

Status (2026-09-09): done. `erf.rans_consistent_diffusivities` (heat,
scalar and moisture diffusivities all rho cmu' sqrt(k) L) and
`erf.rans_lscale_from_pblh` with `erf.rans_lscale_min` (cap
kappa 0.1 zi from `erf.most.pblh_calc = MYNN25`, clamped to
[rans_lscale_min, max_geom_lscale]) landed, both opt-in; the neutral deck
is bit-identical with them off. `Stable_ABL_Flat` (GABLS1, 9 h) passes
every physics check first time: u* 0.244, jet 1.23 Ug at 154 m, BL depth
134 m, KE(0)/u*^2 3.23. `Convective_ABL_Flat` (MS94-B sounding, 0.24 K m/s,
4 h) closes the column heat budget to 0.02 %, holds the inversion at
1020 m, warms the mixed layer 1.08 times the encroachment estimate, and
carries the AL01 buoyancy term at the wall (3.82 against 3.23 neutral).
Two findings: the convective deck aborted at dt = 5 s because K/rho
reaches 40 m2/s, exactly the explicit vertical-diffusion limit of the
anelastic integrator (phases 9-10), so the deck runs at 2 s; and a
local-K closure keeps a superadiabatic lapse of about -2 K/km through
the mixed layer (1.1 K spread where LES gives under 0.3 K), which the
check bounds at 2 K with the gradient sign checked separately. Kynema
comparison: no countergradient term there either, and its Prandtl
function gives the same 2.2 heat-to-momentum ratio at Rt = -3 as cmu'/cmu.
The smoke-mode dissipation-lag tolerance is 10 % (early transient),
5 % in physics mode. Shared check code moved to `rans_checks.py`.

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

Status (2026-09-09): done, with a wall-distance bug found and fixed. The
flat-fitted check showed the Poisson distance short by exactly z dz/(2H)
(0.78 % on 64 cells): `poisson_wall_dist` took the cell's low-face fluxes
as its gradient, half a cell off centre in every direction, which
overstates |grad phi| by dz/2. It now forms a cell-centred gradient from
centred differences with the cell-centre metrics (chain rule for a mesh
deformed in z), needing one ghost cell of phi and the cell's own nodes.
After the fix the flat-fitted distance is exact to 1e-6 m above the first
cell (1.5 cm there, the odd-reflection Dirichlet ghost); on the ridge the
mean error is 1.0 %, the max 5.2 % (at 900 m, above the 30 m cap) and
2.4 m within 100 m of the surface on a 40 m by 15.6 m grid. Kynema uses
the vertical height above the terrain, clamped to dz/2, so there is no
counterpart to compare. `Neutral_Hill_2D` (periodic Witch of Agnesi,
h 100 m, L 500 m, ny = 1, anelastic with MLMG since FFT needs a flat
mesh, blocking factor 1, dt 1.5 s for Courant 0.375 at dx 40 m) passes:
crest speed-up 0.55, 0.45, 0.39 in the three lowest cells against the
Jackson-Hunt 2 h/L = 0.4, upstream log law within 4 %. The Poisson solve
with ny = 1 works without the hidden-direction hint. 1 and 2 ranks agree
to 1e-6 in walldist and 1e-11 in the state, the iterative-solver
tolerances. `RANS_Flat_Fitted_2D` runs the same deck with prob.hmax = 1e-6.
The first trial at dt 5 s aborted from a Courant number of 1.25 in the
no-RANS control as well: a deck error, not the closure.

## Phase 6: 3D terrain canonical case and restart

- `Neutral_Hill_3D/`: Gaussian hill (h/L about 0.2), 3D, neutral, MOST,
  geostrophic forcing; compare hill-top fractional speed-up with the
  Jackson-Hunt estimate (about 2 h/L) and check lee-side k.
- Restart test on terrain: checkpoint mid-run, restart, compare plotfiles
  (wall distance is recomputed in `InitData_post`, Dirichlet k restored).
- Askervein deck (`Real_Terrain/Askervein/inputs_anel`) updated to the new
  input names and run once as the "real" sanity case (not a CTest).

Exit: speed-up within the expected band; restart bit-identical.

Status (2026-09-09): done. `Neutral_Hill_3D` (periodic radial Witch of
Agnesi, h 100 m, L 500 m, 64x64x20 at dx = dy = dz = 40 m, 4 h) passes:
crest speed-up 0.41, 0.28, 0.22 in the lowest cells against the
axisymmetric estimate 1.6 h/L = 0.32, upstream log law within 10 %.
Restart is bit-exact on both the 2D and the 3D terrain decks (checkpoint
at step 20, compare at 40). Askervein runs clean for 20 steps (all fields
finite, first-cell wall distance 7.9 to 10 m on a 20 m cell); its deck
gains nothing beyond the phase-2 flag.

New option `erf.wall_dist_type = terrain_height` (Harish's suggestion,
after the amr-wind immersed terrain and Kynema): height above the local
surface projected on its normal, no linear solve. It is exact to 1e-10 on
a flat fitted mesh and closer to the true distance than the Poisson solve
on both hills (mean 0.02 % and 0.01 % against 1.0 % and 0.3 %), so the
hill decks use it; the Poisson path keeps three `_Poisson` CTest
variants. The Poisson solve is posed positive definite now (same
iterates as before).

Finding, outside this plan: on a 3D terrain-fitted mesh with dz different
from dx (20 or 80 m at dx = 40 m, flat or hill, periodic or inflow,
custom or file terrain, stretched or not, any box layout) the divergence
of the initial field is of order 1e139 before the first projection, the
wall-distance MLMG then diverges (residual 18x after one cycle, 1e10 by
iteration 100, unchanged by the sign convention or by semi-coarsening),
and the run aborts within a step. Askervein is unaffected at dx/dz of
0.5, 1 and 2, so the aspect ratio is not the cause; the RANS code is not
involved (same with Smagorinsky). Two separate defects: (1) the
initialisation one is deterministic, since under `amrex.init_snan = 1`
with the invalid-operation trap armed initialisation completes with no
trap and the bit-identical 1.788e139 divergence, so it is an arithmetic
error that depends on dz relative to dx, not a memory read; (2) the trap
then fires in `ERFPhysBCFunct_w` during the first advance, an
uninitialised read in the w boundary fill that the unit-aspect mesh does
not trigger. Reproducer for (1): `inputs_hill3d amr.n_cell="64 64 40"
prob.hmax=1e-6 max_step=0 erf.mg_v=2 erf.v=1` and read the divergence
before the solve. Harish's rule for the decks: mass inflow and
pressure outflow instead of periodic if periodic turns out to be the
issue; it did not, so the hill decks stay periodic at unit aspect ratio.

## Phase 7: documentation and diagnostics

- Theory: new `Docs/sphinx_doc/theory/RANS.rst` (equations 15-32 of AL01,
  the length-scale limits, the Poisson wall distance, the wall condition).
- `Inputs.rst`: every new key; `RegressionTests.rst`: the Canonical_RANS
  entries; a README per case in the Canonical_LES style.
- Optional plot variables `Rt`, `cmu`, `cmu_prime` (stored in the eddy
  diffusivity container's spare components or a small diagnostic MultiFab).

Exit: docs build; every input key documented.

Status (2026-09-09): done. `Docs/sphinx_doc/theory/RANS.rst` (in the
THEORY toctree after DNSvsLES) covers Eqs. 11-32 of AL01 as implemented,
the bounded unstable length, the cancellation-free smoothing, the wall
condition, both wall distances, the limitations and an input table;
`RegressionTests.rst` gains a Canonical RANS section; every new key is in
`Inputs.rst`; a top-level `Canonical_RANS/README.md` states the rules.
Plot variables `Rt`, `cmu`, `cmu_prime` (three new EddyDiff components
written by the closure) let the stability functions be inspected; the
neutral smoke check verifies them against Eqs. 31-32 to 1e-16 and
Kmv = rho cmu sqrt(k) L. `sphinx-build` (7.4, without the math-dollar
extension) reports no error on the touched pages after two table rows
were realigned. Nine CTest entries and seven gtests pass.

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

Status (2026-09-09): done. All five physics decks rerun on the final code
and pass; `ctest -R RANS_` (9 entries) passes; the closure gtests pass in
double and single precision (one tolerance made precision-aware);
Release build with all warnings on is clean in the changed files;
codespell clean with the repo config; the docs build; upstream
development (two commits, MOST fixes and RemakeLevel) merged without
conflict. PR opened against erf-model/ERF development from
hgopalan:claude-RANS. CI runs `ctest -L regression` in a Debug build, so
the nine entries run there. Phases 9 and 10 (implicit vertical diffusion
under anelastic) and the dycore defects of phase 6 stay open.

## Phase 9: implicit vertical diffusion of scalars under anelastic

Added 2026-09-09 after phase 3 showed that the anelastic integrator, the
main use for RANS, runs every vertical diffusion explicitly (dz^2 / 2K,
about 19 s on the neutral deck) while the compressible path has a column
tridiagonal (Thomas) solve that works with the RANS diffusivities
(4 h compressible, explicit vs implicit at the same dt: 3e-6 m/s in wind,
2e-7 in KE). PR erf-model/ERF#3329 zeroed `vert_implicit_fac` under
anelastic without a stated reason; the momentum solve lives in the
substepping path, which anelastic replaces with a projection, so the
switch reads as a guard against a half-wired configuration.

- Assert at grid creation that no level is decomposed in z (the z entry of
  `amr.max_grid_size` below the domain height) whenever any
  `vert_implicit_fac` is nonzero: the column solves take the box bounds as
  the column and would apply domain BCs at interior faces. This closes a
  silent error in the compressible path as it stands.
- Let the post-stage scalar solve (theta, KE, moisture) run under
  anelastic: stop zeroing the factor for scalars; on the second stage of
  the trapezoidal update, which recovers the first-stage tendency from the
  state difference, apply the implicit operator with half the step (the
  same pattern as the implicit dissipation in phase 3).
- Momentum stays explicit in this phase; its factor stays zero under
  anelastic.

Exit: neutral deck under anelastic at dt = 20 s with implicit scalars
passes the physics checks (it aborts today at dt = 20 s); dt = 5 s matches
the explicit answer to the tolerances above; compressible results
unchanged; CTest entry with `erf.vert_implicit_fac = 1 1` on the neutral
deck.

## Phase 10: implicit vertical diffusion of momentum under anelastic

- Call the momentum tridiagonal on the explicitly updated momenta in the
  no-substep path, before the projection (diffuse, then project, so the
  divergence constraint holds), for the constant-dz, stretched and terrain
  variants; the MOST wall stress enters as it does in the compressible
  path.
- Remove the anelastic zeroing entirely; `erf.vert_implicit` keeps its
  meaning.

Exit: neutral deck under anelastic at dt = 60 s passes the physics checks
and matches the compressible implicit run at dt = 60 s (phase 3, RESULTS)
to the same tolerances; 2D hill deck (phase 5) under anelastic at 4x its
explicit limit runs and matches its explicit answer; restart bit-exact.

