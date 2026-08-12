# WDM6 Bubble Plotfile-Parity Campaign Report — `plotfile_parity_wdm6_bubble_v1a`

Clean-SHA re-run of Milestone A. Serial bridge-vs-native analogue of the nightly
`https://ccse.lbl.gov/pub/RegressionTesting1/ERF/2026-08-12/Bubble_WSM6.html`.

**Supersedes the triage-only caveat in `campaign_wdm6_bubble_v2_report.md`.** That
report remains an accurate record of the v2 campaign and of the heap-corruption
investigation; it is not amended. What changes here is evidentiary standing, not
physics: v2's rows were `dirty_flag=dirty` / `clean_sha_basis=dirty_triage_only`
because the fix was uncommitted. These rows are `clean`.

## Verdict

**Milestone A FAIL, on clean-SHA evidence. Campaign stopped at A; B and C not run**
per stop-diff-retreat. The result reproduces the v2 triage finding **byte-for-byte**.

## Provenance

| Field | Value |
| --- | --- |
| ERF SHA | `9a74610149d5b5c4609134725c69fe5a736d6ce3` |
| Embedded runtime hash | `26.06-292-g9a74610149d5` — matches worktree describe, **no `-dirty`** |
| Worktree | tracked-clean at build time |
| Build | `build_wdm6`, `CMAKE_BUILD_TYPE=Release`, `ERF_ENABLE_MPI=OFF` |
| Ranks | 1, no `mpirun` (serial scope gate) |
| Input | `inputs_BF02_moist_bubble`, `erf.moisture_model=WDM6`, `fixed_dt=0.5` |
| fcompare | `Submodules/AMReX/Tools/Plotfile/fcompare.gnu.ex` (tracked submodule source) |

The embedded hash is recorded in all four run logs as `ERF git hash: 26.06-292-g9a74610149d5`.

### Provenance repair — tier 2 again required

`AMReXBuildInfo.cmake:181` captures the hash with `git describe --abbrev=12 --dirty
--always --tags` inside `execute_process` at **CMake configure time**. Two consequences
drove the sequencing:

- Every tracked modification had to be committed *before* the refresh, or the run would
  be stamped `-dirty`. Two commits landed first: the pending `source-corpus` TSV rows,
  and the naming convention plus its operator promotion.
- Tier 1 (removing the build-info object and its `.o.d`) was **insufficient a second
  time** — the generated source still carried `26.06-287-g4ddfb8c42756-dirty` at
  `AMReX_buildInfo.cpp:137`. Tier 2 (removing `build_wdm6/erf_srclib/AMReX_buildInfo.cpp`)
  was required and triggered a CMake regeneration, exactly as `process-operators.md:11-13`
  anticipates and as v2 recorded.

Only the build-info translation unit recompiled and `erf_exec` relinked; `erf_srclib`
and `amrex_3d` were already up to date. The binary is therefore **source-identical** to
the v2 triage build, which is why byte-identical results were the expected outcome
rather than a surprise.

## Milestone A — step 1, t = 0.5 s

`plt00000` pair: **AGREE** — all 31 variables exactly `0`. The divergence arises
entirely within the first WDM6 call.

`plt00001`, zero tolerance, 31 variables:

| variable | absolute error | relative error |
| --- | --- | --- |
| temp | 6.631607696 | 2.290772895e-02 |
| eq_pot_temp | 6.669060511 | 1.985579005e-02 |
| theta | 6.625510053 | 1.976311604e-02 |
| rhotheta | 3.01690336 | 8.859488505e-03 |
| pressure | 867.1722896 | 8.721967687e-03 |
| pert_pres | 867.1722896 | 1.000524732 |
| qi | 5.347781098e-07 | 2.647491793e-03 |
| rhoQ3 | 2.198031793e-07 | 2.608285552e-03 |
| rhoQ5 | 5.641599512e-06 | 6.373683213e-04 |
| qsnow | 1.225136291e-05 | 6.246766316e-04 |
| nr | 1.225136291e-05 | 6.246766316e-04 |
| nc | 1.238390409e-05 | 6.191952043e-04 |
| qt | 10 | 5.458145434e-08 |
| qv | 3.151323603e-13 | 2.667406784e-11 |
| rhoQ1 | 1.434944557e-13 | 1.034266105e-11 |
| qp | 0.01 | 1.797693135e+308 |
| density, qc, qrain, qgraup, rhoQ2/4/6, x/y/z_velocity, scalar, pres_hse, dens_hse, pert_dens, rhoadv_0 | 0 | 0 |

`qp`'s relative error is `HUGE` because its reference norm is zero, not because the
absolute error is large. `qt` carries the largest absolute error in the table but the
third-smallest relative error. Both rows were absent from the v2 report's table and are
recorded here for completeness; they come from the archived `fcompare` output, not from
prose.

Max-error zone: **`(i,j,k) = (199,3,90)`** — the far corner, ~9 km altitude. Confirmed
independently by `-z temp`, `-z rhotheta` and `-z pressure`.

Archived under `Exec/RegTests/Bubble/plotfile_parity_wdm6_bubble_v1a/` (gitignored, so
these do not survive a clean checkout — regenerate from the ledger rows' run controls):

- `fcompare_A_plt00000_zero_tol.txt` — the AGREE table
- `fcompare_A_plt00001_zero_tol.txt` — the FAIL table above
- `fcompare_A_zone_temp.txt`, `fcompare_A_diff_temp.txt`
- `diffs_plt00001_temp/` — the `-d temp` diff plotfile

This is percent-scale on the thermodynamic state — not epsilon, and not eligible for
variable-aware growth adjudication (`validation-operators.md:32`).

## Rule 36 — restart reproducibility: `MATCH`

Both legs restarted from the shared Fortran-leg `chk00000` and advanced one step. The
restarted 31-variable table diffs **byte-identical** against the direct clean-SHA table,
and `-z temp` returns the same zone `(199,3,90)`. Not a checkpoint, ghost-state, or
restart-state artifact. Rule 37 trigger satisfied.

## Rule 37 — deliberately not re-run

The standing `20260812T182700Z_wdm6_tag_surface_scope_regression` decision establishes
that the tag retreat is inconclusive *by construction*: WDM6 group tags index every
field at `kts` with no `do k` loop, so each tag validates one cell and the surface
cannot observe k=90. Re-running it at clean SHA would produce a ~53 MB marker log
reproducing a known-inconclusive result. Restoring the WSM6 all-k print contract is the
prerequisite, and it is out of scope here.

## Reproduction against v2

Three independent byte-identical diffs, all against
`campaign_wdm6_bubble_v2/fcompare_A_zero_tol.txt`:

1. clean-SHA direct `plt00001` table
2. clean-SHA Rule 36 restarted `plt00001` table
3. the two clean-SHA tables against each other

The v2 physics finding stands unchanged. The clean SHA converts it from triage into
formal evidence.

## Coverage caveat carried forward

**`ni` is still not compared.** Despite the `erf.plot_vars_1` override adding `nc ni nr`,
only 31 variables appear — `nc` and `nr` are present, `ni` is not. `RhoQ8`/`nn` (the CCN
field) is therefore uncompared, so this milestone is not full-coverage even as a FAIL.
Unchanged from v2; an ERF plotfile emission gap, not a tooling artifact.

## Naming restart

This is the first campaign under `plotfile_parity_*`, adopting the WSM6 convention
(`plotfile_parity_v1c`, `plotfile_parity_v1d`) promoted into `process-operators.md` as
`IMPL_CAMPAIGN_NAMING`. The earlier `campaign_wdm6_bubble_v*` line is retired: it was
invented without a corpus check and its numbering does not index the ledgers — v3, v4
and v5 were run but never recorded, and have been deleted. `v1` and `v2` are retained;
they are cited by existing rows.

## What the nightly-page equivalent still needs

1. ~~Commit the heap-corruption fix so Milestone A can be re-run at clean SHA.~~ Done —
   `38a3f6ead`, confirmed by this campaign.
2. Restore the WSM6 all-k print contract in WDM6 group tags, then re-validate the
   frontier — currently 33 single-cell closures. The `tools/build_ci` harness needs the
   same fix (hardwired `micro_diag_target_column`, first-marker-only comparison).
3. Resolve the `theta` writeback design question and re-run A. Still a CANDIDATE site;
   no logic was patched in this session.
4. Fix `ni` plotfile emission so `nn` is compared.
5. Then B and C. Milestone C at step 100 is the direct analogue of the nightly page.
6. Deferred lanes: 4-rank MPI (needs an `ERF_ENABLE_MPI=ON` build), GPU, and blessing a
   `Bubble_WDM6` benchmark.
