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

## Extraction Discipline
- The source corpus includes markdown and TSV files; both must survive the split with clear homes.
- No “good summary” is enough without provenance. Use `source-ledger.md` and `coverage-matrix.md`.
- Default to near-verbatim carry-forward for hard-won rules and operators.
- If normalized wording weakens a requirement, call it out explicitly in the coverage matrix.
