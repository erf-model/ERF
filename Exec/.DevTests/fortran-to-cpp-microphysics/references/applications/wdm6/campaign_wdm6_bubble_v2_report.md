# WDM6 Bubble Plotfile-Parity Campaign Report — `campaign_wdm6_bubble_v2`

Serial bridge-vs-native analogue of the nightly
`https://ccse.lbl.gov/pub/RegressionTesting1/ERF/2026-08-12/Bubble_WSM6.html`.

## Verdict

**Milestone A FAIL. Campaign stopped at A; B and C not run** per stop-diff-retreat.
One blocking runtime defect found and fixed en route.

## Provenance

| Field | Value |
| --- | --- |
| ERF SHA | `4ddfb8c42756ad17c91e976d8827d91d78db97c6` |
| Embedded runtime hash | `26.06-287-g4ddfb8c42756` (repaired; matches HEAD) |
| Worktree | dirty — `Source/Microphysics/WDM6/ERF_AdvanceWDM6.cpp` only |
| Build | `build_wdm6`, `CMAKE_BUILD_TYPE=Release`, `ERF_ENABLE_MPI=OFF` |
| Ranks | 1, no `mpirun` (Rule 35 serial scope gate) |
| Input | `inputs_BF02_moist_bubble`, `erf.moisture_model=WDM6`, `fixed_dt=0.5` |
| fcompare | `tools/fcompare_serial`; cross-checked identical against `Submodules/AMReX/Tools/Plotfile/fcompare.gnu.ex` |

Provenance was repaired before any run: build-info tier 1 (object + dep removal)
was **insufficient** because the generated source itself carried the stale hash;
tier 2 (removing `build_wdm6/erf_srclib/AMReX_buildInfo.cpp`) was required and
triggered a CMake regeneration, as `process-operators.md:13` anticipates.

All rows below are **triage-only**: the fix is uncommitted, so no formal
manifest promotion is claimed.

## Blocking defect found and fixed: bridge heap corruption

The Fortran bridge had never completed a single step at any recent SHA — an
abort recorded open in the source-corpus ledger since G5a as
*"free(): invalid pointer downstream of G5a"*. It was plausibly a heavy-diagnostic
buffering artifact (`wsm6_validation_operator.md:337`). **It was not.** Under full
campaign posture — `microphysics_debug=0`, `GFORTRAN_UNBUFFERED_ALL=y`, no
`mpirun` — it still aborted, as
`Fatal glibc error: malloc.c:4376 (_int_malloc): assertion failed`.

`valgrind` memcheck localized it exactly: **28 invalid reads + 15 invalid writes**,
all inside `mp_wdm6_run_c`, all past 6,400-byte blocks allocated by
`FArrayBox::FArrayBox` in `WDM6::Advance`. 6,400 bytes = 800 doubles = 200×4 =
the **active tile** extent.

Root cause: the seven 2D precip buffers (`rainacc`, `rainncv`, `sr`, `snowacc`,
`snowncv`, `graupacc`, `graupelncv`) were allocated on the **active** tile box,
while `mp_wdm6_run_c` declares them `dimension(ims:ime, jms:jme)` — i.e. on
**storage** bounds, which `WDM6::Advance` passes as `imlo,imhi,jmlo,jmhi`. Any
ghost extent makes the Fortran side index past the allocation.

This violates two documented core rules:
- `core-rules.md:37` — "Respect slab/box allocation patterns for 2D accumulators
  and other mixed-dimensional storage."
- `core-rules.md:16` — "Pass storage bounds and active bounds separately when the
  Fortran scheme distinguishes memory extents from compute extents."

Fix mirrors the working WSM6 pattern at `ERF_AdvanceWSM6.cpp:963-994` — allocate
on `fab_box2d` (storage slab) and initialize over `fab_box2d`, keeping `box2d` for
valid-region updates. WSM6 even carries the explanatory comment, which is now
carried across.

Result: `ERROR SUMMARY: 0 errors from 0 contexts`, and the bridge leg completes
step 1 and writes `plt00001` — the first time that has happened at a recent SHA.

## Milestone A — step 1, t = 0.5 s

`plt00000` pair: **AGREE** (identical initial state). So the divergence arises
entirely within the first WDM6 call.

`plt00001`, zero tolerance, 31 variables:

| variable | absolute error | relative error |
| --- | --- | --- |
| temp | 6.631607696 | 2.291e-02 |
| eq_pot_temp | 6.669060511 | 1.986e-02 |
| theta | 6.625510053 | 1.976e-02 |
| rhotheta | 3.01690336 | 8.859e-03 |
| pressure | 867.1722896 | 8.722e-03 |
| pert_pres | 867.1722896 | 1.000524732 |
| qi | 5.347781098e-07 | 2.647e-03 |
| rhoQ3 | 2.198031793e-07 | 2.608e-03 |
| qsnow | 1.225136291e-05 | 6.247e-04 |
| nr | 1.225136291e-05 | 6.247e-04 |
| nc | 1.238390409e-05 | 6.192e-04 |
| rhoQ5 | 5.641599512e-06 | 6.374e-04 |
| qv | 3.151323603e-13 | 2.667e-11 |
| rhoQ1 | 1.434944557e-13 | 1.034e-11 |
| density, qc, qrain, qgraup, rhoQ2/4/6, x/y/z_velocity, scalar, pres_hse, dens_hse, pert_dens, rhoadv_0 | 0 | 0 |

Max-error zone (`-z temp`, `-z rhotheta`, `-z pressure` all agree):
**`(i,j,k) = (199,3,90)`** — the far corner, ~9 km altitude.

Diff artifact archived: `Exec/RegTests/Bubble/campaign_wdm6_bubble_v2/diffs_plt00001_temp/`
(untracked: `campaign_*/` is gitignored, so run artifacts do not survive a clean
checkout -- regenerate from the ledger row's run controls if needed).

This is percent-scale on the thermodynamic state — not epsilon, and not eligible
for variable-aware growth adjudication.

## Rule 36 — restart reproducibility: `MATCH`

Both legs restarted from the shared Fortran-leg `chk00000` and advanced one step
reproduce the divergence **bit-for-bit** (`temp 6.631607696`,
`rhotheta 3.01690336`, `pressure 867.1722896`, `qi 5.347781098e-07`). Not a
checkpoint, ghost-state, or restart-state artifact. Rule 37 trigger satisfied.

## Rule 37 — tag retreat: inconclusive *by construction*

Tag pair re-run at the **failing** column `(199,3)` rather than the `(ilo,jlo)`
default, since `(100,0)` already passes everywhere. Result: **all 33 groups
`P0..G17` report `max_abs_diff = 0.000e+00`** — while the plotfile diverges at
`(199,3,90)`.

The retreat cannot resolve this, because WDM6's group tags index every field at
`kts` with no `do k` loop:

```fortran
'WDM6-FORT_PRE_G7', kts, qrs(i_dbg_local,kts,1), ... t(i_dbg_local,kts), ...
```

The WSM6 contract on branch `wsm6-mpi-cpu` instead validates a full column:

```fortran
if (mpdbg_level >= 1 .and. loop == 1) then
  do k = kts, kte
    write(*,'(A,I3,6E24.16)') 'WSM6-FORT_DENFAC ', k, denfac(its,k), ...
  enddo
endif
```

So **each WDM6 group tag validates one cell; the WSM6 equivalent validated ~100.**
Every WDM6 group closure rests on roughly 1/100th of the evidence its WSM6
counterpart did, and no group evidence exists at the k where Milestone A fails.
Separately, `G10A/G10B/G10C` emit ~20,000 Fortran marker lines (full i-k plane)
against one C++ line per group, so those groups have no matched surface either.

Group statuses are not reopened on this evidence, but they are qualified as
**single-cell closures**.

## Candidate divergence site (not confirmed)

A structural asymmetry sits outside every group tag. The Fortran branch converts
updated absolute temperature back to potential temperature after the bridge call:

```cpp
Real exner = std::pow(p_arr(i,j,k) / p0, rdOcp);
theta_arr(i,j,k) = t_arr(i,j,k) / exner;      // Fortran branch only
```

The native branch never writes `mic_fab_vars[theta]` — it ends at precipitation
accumulation. `Copy_Micro_to_State` then forms `RhoTheta = rho * theta` from that
field, so the native path carries a **pre-microphysics theta** and drops
microphysical latent heating from the thermodynamic state. That matches the
signature precisely: moisture agrees, thermodynamics diverges at percent scale.
The site is after `POST_G17`, outside all group boundaries — which is why no tag
observes it.

Held as **candidate, not root cause**: retreat discipline forbids patching logic
without paired `PRE_<site>`/`POST_<site>` first-line evidence, the design choice
(mirror the bridge conversion vs. have the native kernel update `theta` directly)
is a deliberate one, and the residual `qi`/`qsnow`/`nc`/`nr` differences are not
explained by writeback alone.

## Coverage corrections

- **Bubble is not warm-phase-only.** `qi` and `qsnow` carry nonzero differences
  and the failing zone is at ~9 km, below freezing. Cold-phase groups are active
  at upper levels despite an ice-free initial state. SquallLine remains the
  stronger cold-phase gate, but the earlier "warm-rain only" framing was wrong.
- **`ni` is not being compared.** Despite the `erf.plot_vars_1` override adding
  `nc ni nr`, only 31 variables appear — `nc` and `nr` are present, `ni` is not.
  `RhoQ8`/`nn` (the CCN field) is therefore uncompared. Must be resolved before
  any milestone is treated as full-coverage.

## What the nightly-page equivalent still needs

1. Commit the heap-corruption fix so Milestone A can be re-run at clean SHA.
2. Restore the WSM6 all-k print contract in WDM6 group tags, then re-validate the
   frontier — currently 33 single-cell closures.
3. Resolve the `theta` writeback design question and re-run A.
4. Fix `ni` plotfile emission so `nn` is compared.
5. Then B and C. Milestone C at step 100 is the direct analogue of the nightly page.
6. Deferred lanes: 4-rank MPI (needs an `ERF_ENABLE_MPI=ON` build), GPU, and
   blessing a `Bubble_WDM6` benchmark.
