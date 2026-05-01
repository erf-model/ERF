# WSM6 Validation Operator (Agentic, Rule-Gated)

Purpose: move rows in `validation_manifest.tsv` from `PENDING` to `PASS` / `EPSILON_OK` / `FAIL` / `RETREAT` using the rule framework in `fortran_to_cpp_microphysics_skill.md` (Rules 30/30A/30B/31/32/32.6/33/34).

This workflow is intentionally agentic and small-process. Do not introduce a large automation script.

## Execution and Decision Ledgers

`validation_runs.tsv` is the legacy combined ledger
for historical rows recorded before the ledger split.

For all new formal validation runs:
- `executions.tsv` records one executable invocation
  per row (execution facts only).
- `decisions.tsv` records manifest or frontier state
  changes, each citing one or more `execution_id` values.
- `validation_manifest.tsv` remains the current-state
  frontier table.

No manifest change is valid unless a `decisions.tsv`
row cites clean-SHA execution evidence.
`dirty_flag=dirty` rows are acceptable for triage
and debug only — not for promotion or re-open decisions.

---

## 1. Inputs and Outputs

Inputs:
- `Docs/wsm6_deep_dive/validation_manifest.tsv`
- `Docs/wsm6_deep_dive/executions.tsv`
- `Docs/wsm6_deep_dive/decisions.tsv`
- `Docs/wsm6_deep_dive/validation_runs.tsv` (legacy only)
- `Docs/wsm6_deep_dive/fortran_to_cpp_microphysics_skill.md`
- `Exec/CanonicalTests/SquallLine_2D/inputs_moisture_WSM6`

Outputs:
- Updated `validation_manifest.tsv` statuses and divergence variables.
- Appended row(s) in `executions.tsv` for each run leg.
- Appended row(s) in `decisions.tsv` for each manifest/frontier status change.

## 2. Non-Negotiable Rules

1. Runtime toggles only:
- Fortran path: `erf.use_wsm6_cpp_answer=0`
- C++ path: `erf.use_wsm6_cpp_answer=1`
- Tier1 canonical: `erf.microphysics_debug=1`, `erf.micro_diag_mode=canonical|both`
- Tier2 forensic: `erf.microphysics_debug>=2`, `erf.micro_diag_mode=forensic|both`

2. For any host-side print after GPU kernels, call `Gpu::synchronize()` first.

3. Stop, Diff, Retreat:
- On first divergence beyond epsilon, stop forward validation.
- Record first divergent group/tag/step.
- Mark current group `FAIL` and downstream groups `RETREAT` as needed.
- Retreat method is two-stage and fidelity-first:
  1) Block-by-block: move process-group by process-group to isolate the first failing block.
  2) Line-by-line inside that block: prove the first failing statement with paired diagnostics around one line/block:
     `PRE_<site>` (inputs/local vars) and `POST_<site>` (outputs/local vars).
- Per-line fidelity is required before changing logic: inputs must still match across Fortran/C++ at `PRE_<site>` and first diverge at `POST_<site>`.
- Bisection inside the block is allowed, but every proposed fix must be backed by this per-line evidence.
- Record the confirmed site in `validation_runs.tsv` notes (for example `first_line_divergence=<file>:<line>`).
- Retreat print granularity:
  - `erf.microphysics_debug=1` is canonical-tag mode for forward validation.
  - `erf.microphysics_debug>=2` is one-value-per-line forensic mode.

4. Use tracked assets only for formal validation.

5. Tier2 schema/format contract (when forensic mode is enabled):
- Record shape:
  `WSM6-DIAG-T2 diag_schema=... tag=... phase=... source_layer=... path_id=... expr_id=... store_id=... loop=... i_dbg=... j_dbg=... k_dbg=... k_raw=... debug_level=... var=... value=...`
- Field order is contractual.
- `i_dbg,j_dbg` are target-column coordinates in outer C++ space.
- `k_dbg` is normalized 1..nk; `k_raw` is raw storage index.
- Numeric token must preserve ULP-level distinguishability for double precision (minimum 17 significant digits; current contract uses 20 fractional digits with explicit sign and fixed exponent width).

6. When Tier2 is required (skill-level trigger):
- Use Tier1 first for normal frontier progression.
- Enter Tier2 only after Tier1 identifies the first failing group/tag, or when an index/column-label ambiguity prevents trusted interpretation.
- Do not run broad Tier2 traces during forward pass validation; scope to a single failing tag/group and a single target column.

## 3. Setup IDs and Pass Gates

Supported setup IDs in manifest:
- `1-step_dry`
- `1-step_debug|2-step_debug|10-step_debug`

Interpretation:
- `1-step_*`: initial-state and immediate copyback checks.
- `2-step_*`: first sequential accumulation check.
- `10-step_*`: short drift/compounding check.

Acceptance:
- `PASS`: bitwise equal at required tags.
- `EPSILON_OK`: within approved thresholds (Density/RhoTheta) when bitwise exact is not required by rule context.
- `FAIL`: exceeds thresholds or tag mismatch.
- `RETREAT`: blocked by earlier failing dependency.

## 4. Per-Group Operator Loop

For each `group_id` still `PENDING`:

1. Select setup from `setup_id`.
2. Run Path A (Fortran bridge) and capture `log.fortran`.
3. Run Path B (native C++) and capture `log.cpp`.
4. Compare matching `fortran_tag` vs `cpp_tag` streams for the target column (Tier1).
5. In retreat mode, compare Tier2 schema lines first on structure, then value deltas.
6. Determine earliest divergence location (group/tag/step/variable).
7. Update manifest row:
- `status`
- `divergence_variable` (first failing variable; keep concise)
8. Append run metadata to `executions.tsv`; append status/frontier actions to `decisions.tsv`.

If current group passes, proceed to next group. If fails, stop forward progress and apply retreat strategy.

## 5. Command Templates (Minimal)

From `Exec/CanonicalTests/SquallLine_2D` using the actively built executable (default shown is DEBUG+MPI):

```bash
# Path A (Fortran bridge; unbuffered recommended for heavy print modes)
GFORTRAN_UNBUFFERED_ALL=y mpirun -n 1 ../../ERF3d.gnu.DEBUG.MPI.ex inputs_moisture_WSM6 \
  erf.use_wsm6_cpp_answer=0 erf.microphysics_debug=1 erf.micro_diag_mode=canonical max_step=<N> \
  > log.fortran 2>&1

# Path B
mpirun -n 1 ../../ERF3d.gnu.DEBUG.MPI.ex inputs_moisture_WSM6 \
  erf.use_wsm6_cpp_answer=1 erf.microphysics_debug=1 erf.micro_diag_mode=canonical max_step=<N> \
  > log.cpp 2>&1
```

Tag-focused compare (example):

```bash
rg "WSM6-(FORT|CPP)_(DENFAC|QSAT|SLOPE1|NISLFV_R|NISLFV_SG|FALL|SLOPE2|MELT|VICE|PRECIP|PHASE|SLOPE3|PRAUT|PRACW|PREVP|PRACI|PSACI|PRACS|PSEML|PIDEP|PSAUT|PSEVP|UPDATE|QSAT2|PCOND)" log.fortran log.cpp
```

Tier2 retreat pair (target column controlled at runtime):

```bash
# Path A (Fortran bridge; use unbuffered mode for Tier2 verbosity)
GFORTRAN_UNBUFFERED_ALL=y mpirun -n 1 ../../ERF3d.gnu.DEBUG.MPI.ex inputs_moisture_WSM6 \
  erf.use_wsm6_cpp_answer=0 \
  erf.microphysics_debug=2 \
  erf.micro_diag_mode=forensic \
  erf.micro_diag_tags=DENFAC \
  erf.micro_diag_expr=standing \
  erf.micro_diag_store=standing \
  erf.micro_diag_target_column="3 3" \
  max_step=1 > log.fortran.t2 2>&1

# Path B
mpirun -n 1 ../../ERF3d.gnu.DEBUG.MPI.ex inputs_moisture_WSM6 \
  erf.use_wsm6_cpp_answer=1 \
  erf.microphysics_debug=2 \
  erf.micro_diag_mode=forensic \
  erf.micro_diag_tags=DENFAC \
  erf.micro_diag_expr=standing \
  erf.micro_diag_store=standing \
  erf.micro_diag_target_column="3 3" \
  max_step=1 > log.cpp.t2 2>&1
```

Tier2 compare policy:
- First compare schema/order/index tokens and path metadata.
- Then compare numeric values.
- If a provenance token is intentionally path-specific (for example source-layer naming), normalize that token explicitly and document the normalization in `executions.tsv` notes.

## 6. Manifest Update Convention

`validation_manifest.tsv`:
- `divergence_variable` is the first failing variable only.
- Keep `compare_variables` unchanged unless scope itself changes.
- Keep `print_contract` unchanged unless formatting contract changes.

`executions.tsv` append fields (minimum):
- `run_id`, `timestamp_utc`, `git_sha`, `setup_id`, `use_wsm6_cpp_answer`, `microphysics_debug`, `max_step`, `first_divergence_group`, `first_divergence_tag`, `first_divergence_step`, `status`, `notes`.
- `git_sha` policy:
  - `*-dirty` rows are allowed for debug/triage evidence, but are not sufficient for formal status promotion.
  - Any `validation_manifest.tsv` status/frontier change (`PASS`, `EPSILON_OK`, new `FAIL` frontier, downstream `RETREAT` shift) should be backed by a rerun from a clean commit SHA.
  - Do not use `git commit --amend` to update failing-run provenance; prefer a new small checkpoint commit, rerun, then append new rows with the new clean SHA.

`decisions.tsv` append fields (minimum):
- `decision_id`, `timestamp_utc`, `manifest_row_id|group_id`, `old_status`, `new_status`, `reason`, `evidence_execution_ids`, `operator_notes`.
- Every decision row that changes manifest/frontier state must cite one or more clean-SHA execution IDs.

## 7. Dependency Ordering

Validate in process order; do not mark downstream groups PASS when upstream failed:
- G1b -> G1c -> G4 -> G5b -> G5c -> G5d -> G6 -> G7 -> G8 -> G9 -> G10 -> G11 -> G13a..G13j -> G14 -> G15 -> G16

## 8. Practical Guardrails

- Prefer single-rank (`-n 1`) for deterministic diagnostics.
- Keep diagnostics at level 1 unless actively isolating a local issue.
- If heavy Fortran diagnostic output aborts with `corrupted double-linked list`, rerun with `GFORTRAN_UNBUFFERED_ALL=y` and record that env in `executions.tsv` notes.
- Avoid broad code edits during validation runs; isolate one suspected cause at a time.
- Commit manifest/run-ledger updates in small, reviewable chunks.
- Enforce strict Rule 30 group-boundary snapshots: emit each canonical tag at the immediate process-group boundary (`NISLFV_R` right after `nislfv_rain_plm`, `NISLFV_SG` right after `nislfv_rain_plm6`), not from a later aggregate `POST-G*` state.
- If a group boundary value is later copied/reduced in subsequent sub-steps, print from boundary snapshot buffers/derived temporaries captured at the boundary, and note `strict_group_boundary_snapshot` in `validation_runs.tsv` notes.

## 9. Plotfile Parity Lane (Rules 35-38)

### What a campaign is

A campaign is one run pair compared at physics-grounded milestone
checkpoints. Each passed milestone is both a promotion gate and a
knowledge artifact. Rule details are in Rules 35-38 of the skill doc.
Scheme-specific milestone definitions and command templates belong in
implementation notes appendices.

### Trigger condition

Activate this lane after the active tag-frontier closes at clean SHA.
For WSM6, this means the frontier is closed through `G16` with no active
`PENDING`, `FAIL`, or `RETREAT` rows.

### Scope gate: serial parity first

Section 9 is serial-only:
- single-rank
- CPU path
- non-MPI executable mode

MPI and GPU validation are deferred to later scope after serial C++ parity
is verified. Those future lanes require distinct `validation_lane` values,
separate acceptance criteria, and separate operator sections.

### microphysics_debug policy

- Plotfile parity runs: `erf.microphysics_debug=0`
- Tag retreat re-entry: `erf.microphysics_debug=1`
- Line-level retreat isolation: `erf.microphysics_debug>=2`

### Sequential milestone comparison

After both run legs complete:
1. Compare milestones from earliest to latest.
2. If `EPSILON_OK`, record active regime metadata and continue.
3. On first failure, stop immediately.
4. Apply Rule 36 (restart reproducibility), then Rule 37 (tag retreat).
5. Do not compare later milestones while an earlier one fails.

If the failing milestone was detected on coarse plot cadence (`plot_int_1 > 1`),
run a bounded restart refinement from the nearest earlier checkpoint with
`plot_int_1=1`, compare per-step explicit pairs, and identify earliest failing
substep `N*`. Use `N*` (not the coarse milestone step) for Rule 36 and Rule 37.
Use a shared restart source for both legs (same checkpoint path, typically the
Fortran-leg checkpoint) during refinement and restart reproducibility checks.
If single-step checkpoint breadcrumbs are needed for tag retreat handoff, set
`check_int=1` during the bounded refinement run.

### Dependency ordering extension (extends Section 7)

After terminal tag-frontier closure at clean SHA, open plotfile campaign.
On milestone failure:
- Rule 36 -> Rule 37 -> re-enter Section 4 loop at mapped tag group.
- Mark that group `FAIL`; mark downstream groups `RETREAT` as required.
- Resume plotfile campaign only after upstream re-validation at clean SHA.

### Campaign directory discipline

Use a unique `campaign_id` directory root per campaign execution.
Do not reuse campaign output directories.
Archive only after corresponding `validation_runs.tsv` rows are closed.
