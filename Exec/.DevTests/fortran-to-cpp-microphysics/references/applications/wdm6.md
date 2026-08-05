# WDM6 Application Note

## Current Branch and PR Lineage
- current local branch: `wdm6-mpi-cpu`
- current head: `f23e571c8`
- fetched PR mirror: `pr-3517-wdm6`
- upstream PR: `erf-model/ERF#3517`
- PR title: `[WDM6] Initial CPU fortran implementation`
- PR head owner/branch: `adamwise95:WDM6`
- PR base: `development`

## Status as of August 3, 2026
- Path A bridge: partial
  - build-system wiring, interfaces, and wrappers exist
  - Fortran module and isohelper are present
  - the branch README still describes the Fortran module as a stub needing WRF-source completion
- Path B native C++: partial
  - infrastructure, state mapping, and helper scaffolding exist
  - `ERF_AdvanceWDM6.cpp` contains dual-mode structure and a partial kernel path
  - the branch README still describes full physics as incomplete

This means WDM6 is not yet a validated reusable pattern the way mature WSM6/Morrison rules are.

## State Mapping and Scheme-Specific Facts
- WDM6 is double-moment for liquid only.
- Mixing-ratio state follows the WSM6-style six-class layout:
  - `RhoQ1..RhoQ6` map to `qv, qc, qi, qr, qs, qg`
- Number concentrations add:
  - `RhoQ7 -> nc`
  - `RhoQ8 -> nn`
  - `RhoQ9 -> nr`
- `nn` is aerosol/CCN number, not Morrison-style ice number.
- WDM6 introduces `wdm6.ccn0` and `xland` behavior that WSM6 does not have.

## Current Working Artifacts
- canonical file-order group map:
  `references/applications/wdm6/group_map.md`
- use the tracked Fortran source as the frozen kernel reference:
  `Source/Microphysics/WDM6/ERF_module_mp_wdm6.F90`
- use the untracked top-level `module_mp_wdm6.F` only as comment/naming guidance,
  not as the authoritative line-number source

## Entry Point Comparison
Framework-level wiring matches WSM6 and Morrison:
- `MoistureType::WDM6` is registered in ERF moisture-model selection.
- `modelType()` classifies it as Eulerian.
- `EulerianMicrophysics` dispatches it through `SetModel<WDM6>()`.
- It implements the expected `Init`, `Copy_State_to_Micro`, `Advance`, and `Update` pattern.

Current divergences from mature patterns:
- `WDM6::Init` seeds `nn` from `wdm6.ccn0` and preserves it through the first state-copy path.
- `WDM6::Advance` hardwires some debug behavior and target-column defaults in the bridge call, unlike the more campaign-configurable WSM6 path.
- The bridge uses WSM6-like argument packing but extends it with `nn/nc/nr`, `ccn0`, and `xland`.

## Bubble vs SquallLine Posture
- Bubble is the preferred newer testing framework in team practice, but this branch does not yet include an in-repo WDM6 Bubble input.
- The current in-repo moist Bubble case still defaults to Kessler and only comments Morrison as an alternative.
- The existing mature WSM6 parity corpus is centered on SquallLine, not Bubble.
- The WDM6 README points to external Bubble-oriented verification assets outside this repo.

Practical consequence:
- Bubble is the likely future WDM6 validation surface.
- SquallLine remains the better-documented in-repo parity precedent because the WSM6 operator corpus is built around it.

## Bubble Lane Conventions
- For formal Bubble validation, use distinct run roots per path and SHA:
  `erf.plot_file_1=<campaign>/<path>/plt`
  `erf.check_file=<campaign>/<path>/chk`
  `> <campaign>/<path>/run.log 2>&1`
- Keep Bubble run artifacts in ignored directories or stash them before clean-SHA reruns so the worktree stays suitable for formal evidence capture.
- When bridge and native runs share the same Bubble input, vary only the intended comparison knobs such as `erf.use_wdm6_cpp_answer`, debug level, and target column.

## Hardwired or Non-General Assumptions
- hardwired diagnostic defaults inside the WDM6 bridge path should not become general rules
- README status text describing “stub/incomplete” state is branch-specific, not reusable skill guidance
- external WDM6 verification paths in the README are evidence pointers only

## What Prevents Promotion to Generic Rules Today
- the branch README still marks both the C++ and Fortran implementations as incomplete
- no mature WDM6 validation operator corpus exists yet
- no in-repo WDM6 Bubble campaign assets are wired the way WSM6 SquallLine assets are
- debug hardwiring has not yet been normalized into a reusable operator contract

## Extraction Policy
- promote WDM6 lessons only when they generalize beyond liquid double-moment specifics
- treat `nn/nc/nr`, `ccn0`, and `xland` handling as scheme-local until a second double-moment case confirms the pattern
- use this note as the first-stop context for future WDM6 work before touching generic rules
