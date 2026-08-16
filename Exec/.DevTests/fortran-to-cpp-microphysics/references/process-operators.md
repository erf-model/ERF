# Process Operators

These are campaign habits and operator-level behaviors that matter but are not core physics rules.

## Evidence Capture
- Treat branch history as evidence, not as the primary instruction surface.
- Keep source pointers short and factual in application notes.
- Preserve first-fail artifacts immediately; never rely on overwritable default outputs.
- Record why a decision changed, not just that it changed.
- If clean worktree SHA and runtime embedded git hash disagree, treat the run as triage-only until provenance is repaired.
- Repair provenance with the lightest path first: treat the generated build-info object and dependency file as the first refresh tier, and remove both before the first rebuild check (for example `build_wdm6/Exec/CMakeFiles/buildInfoerf_srclib.dir/__/erf_srclib/AMReX_buildInfo.cpp.o` and `build_wdm6/Exec/CMakeFiles/buildInfoerf_srclib.dir/__/erf_srclib/AMReX_buildInfo.cpp.o.d`). Only remove the generated build-info source as the second-tier fallback if that object+dep refresh still fails to repair the embedded hash (for example `build_wdm6/erf_srclib/AMReX_buildInfo.cpp`).
- After each build-info refresh tier, rebuild and re-check the embedded runtime hash before escalating to the next file-removal tier.
- If removing the generated `AMReX_buildInfo.cpp` forces CMake to reconfigure or regenerate build files, allow that regeneration to complete before judging the provenance repair result.
- Do not run broad reconfigure or near-full rebuild just to refresh embedded hash provenance unless the lightweight build-info refresh fails.
- Batch commits before repairing provenance, and repair immediately before the
  formal run rather than after each commit. On a CMake build the embedded hash
  is captured at configure time, so **every** commit re-stales it, including
  ledger-only ones. Repairing per commit multiplies the cost by the number of
  commits for no benefit. The correct order is: land all code commits, repair
  once, run the formal gate, then commit the ledger rows the gate produced.
  Repairing provenance needs no confirmation; it is a routine step of the gate.

## Debug Discipline
- Keep compile-time bridge toggles separate from runtime diagnostic levels.
- Prefer host-gated diagnostics around kernels over device-side print noise.
- Start with summary prints and escalate to exhaustive dumps only for intentionally small forensic cases.
- When a diagnostic surface becomes contractual, preserve schema and field ordering.

## Campaign Management
- Serial-first acceptance is the default unless a rule explicitly opens MPI/GPU parity scope.
- Re-opening an early group invalidates downstream closure until the earlier scope is re-established.
- Use bounded refinement when a coarse milestone fails; do not jump directly to sprawling manual scans.
- Keep operator workflows small-process and agentic. Avoid replacing judgment with a monolithic automation script.
- Re-check the source corpus when extracted rules are too vague to choose a low-cost operator action safely.
- Re-check the source corpus before inventing new campaign habits around provenance, rerun hygiene, artifact layout, or ledger practice.
- If the same underspecified operator question recurs across turns and the source corpus has a stable answer, promote that answer into the extracted skill package.

## Campaign Artifact Naming
- Use a unique `<campaign_id>` root per campaign, never reused.
- Name the root for the validation lane, not the word "campaign": `plotfile_parity_v1c`, `plotfile_parity_v1d`. Add scheme and case when one ledger serves more than one of them.
- Version the root as `v<number><letter>`. The letter advances within a campaign generation as scope or SHA changes; the number advances for a new generation.
- Derive the ledger ids from the root so artifacts and rows stay greppable together:
  - `run_pair_id` = `<campaign_id>_<pair><scope token>`
  - `execution_id` = `<UTC timestamp>Z_<activity>_<version>_<pair>_<leg>_dbg<level>`
- Encode build type in the id as a `_debug` or `_release` token and record it in notes as `build_type=debug` or `build_type=release`. Do not add build_type as a TSV column — the notes token is enough.
- Bake the clean short SHA into the directory name for retreat and forensic sub-runs, whose provenance varies run to run. Omit it when a whole campaign sits at one SHA that `git_sha` already records.
- Keep the ignore rule aligned with the naming scheme. A root that is not matched by `.gitignore` puts run artifacts into the worktree and breaks the clean-SHA gate.

## Extraction Discipline
- The source corpus includes markdown and TSV files; both must survive the split with clear homes.
- No “good summary” is enough without provenance. Use `source-ledger.md` and `coverage-matrix.md`.
- Default to near-verbatim carry-forward for hard-won rules and operators.
- If normalized wording weakens a requirement, call it out explicitly in the coverage matrix.
