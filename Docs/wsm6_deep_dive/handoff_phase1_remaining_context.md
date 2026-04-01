# WSM6 Phase-1 Handoff (Remaining Context Budget)

## Progress Log
- Section 1 completed on branch `wsm6-moisture-stubs`.
  - Commit: `672df6d0`
  - Summary: added WSM6 scaffold files under `Source/Microphysics/WSM6/` and wired Make/CMake sources/includes.
- Section 2 completed (selection plumbing).
  - Added `WSM6` to `MoistureType` enum in `Source/DataStructs/ERF_DataStruct.H`.
  - Added WSM6 moisture index mapping as a 6-species model (`qv,qc,qi,qr,qs,qg`).
  - Added WSM6 to Eulerian model detection in `Source/Microphysics/ERF_Microphysics.H`.
  - Added `#include "ERF_WSM6.H"` and `SetModel<WSM6>()` in `Source/Microphysics/ERF_EulerianMicrophysics.H`.
  - Added WSM6 handling to accumulation-copy branch in `Source/IO/ERF_Plotfile.cpp` so `rain_accum/snow_accum/graup_accum` are not skipped.

## Current State
- Workspace progress was stashed before deep-dive work.
- Deep-dive artifacts are complete:
  - `Docs/wsm6_deep_dive/entrypoints.md`
  - `Docs/wsm6_deep_dive/args_map.csv`
  - `Docs/wsm6_deep_dive/state_map.csv`
  - `Docs/wsm6_deep_dive/risk_log.md`
  - `Docs/wsm6_deep_dive/implementation_checklist.md`
- Decision baseline is locked:
  - Morrison setup template = initial integration + critical fixes.
  - Canonical external source = `../WRF_model`; `../WRF_NCAR` delta-only.
  - Phase 1 goal = Fortran-in-place runnable WSM6 in ERF.

## Goal for Remaining Context
Implement the highest-leverage subset of Phase 1 that gets us as close as possible to a compile-time integrated WSM6 path, with minimal risk and minimal surface area.

## Recommended Scope to Implement Now (in order)

### 1. Create WSM6 skeleton files and wire into build
Why now: highest leverage, low ambiguity, mostly mechanical.

Files to create:
- `Source/Microphysics/WSM6/ERF_WSM6.H`
- `Source/Microphysics/WSM6/ERF_InitWSM6.cpp`
- `Source/Microphysics/WSM6/ERF_AdvanceWSM6.cpp`
- `Source/Microphysics/WSM6/ERF_UpdateWSM6.cpp`
- `Source/Microphysics/WSM6/ERF_WSM6_Fortran_Interface.H`
- `Source/Microphysics/WSM6/Make.package`

Files to edit:
- `Exec/Make.ERF` (add `ERF_MOISTURE_WSM6_DIR` include block, mirroring existing moisture package pattern)
- CMake files currently used for microphysics source inclusion (same pattern as Morrison integration)

Acceptance criteria:
- New files exist with compilable class/function declarations and stub bodies.
- Build system includes WSM6 source directory/package (no runtime behavior yet required).

### 2. Add moisture-model selection plumbing
Why now: required to make WSM6 addressable from inputs and dispatch path.

Files likely to edit:
- `Source/Microphysics/ERF_Microphysics.H`
- `Source/Microphysics/ERF_EulerianMicrophysics.H`
- Any `MoistureType` enum definition location in ERF data/solver structs

Acceptance criteria:
- `MoistureType::WSM6` exists.
- Eulerian microphysics can instantiate `WSM6` model class.
- No regressions for existing moisture models.

### 3. Add minimal WSM6 state mapping and lifecycle stubs
Why now: enables incremental integration without full Fortran call complexity.

Implement in WSM6 class:
- `Define`
- `Init`
- `Copy_State_to_Micro`
- `Advance` (initially stubbed/no-op or guarded call)
- `Copy_Micro_to_State`
- `Qmoist_Size`, `Qstate_Moist_Size`, restart var hooks

Acceptance criteria:
- ERF compiles with WSM6 model selected at compile time.
- WSM6 model class follows the same interface contracts as other Eulerian models.

Status:
- Completed in scaffold form.
- Remaining in this area is only Fortran bridge insertion into `WSM6::Advance` and any compile-fix fallout.

## Defer to Next Pass (do not spend remaining context here)

### A. Full Fortran argument bridge implementation
- `mp_wsm6_init`/`mp_wsm6_run` call signatures and full argument marshalling.
- Optional arguments: radiation effective radii, chemistry hooks.

Reason to defer: higher complexity and higher bug risk; better done after skeleton + dispatch compile cleanly.

### B. Runtime validation case setup and physics debugging
- New regtest input case for WSM6.
- Numerical debugging/field validation.

Reason to defer: depends on successful compile and basic bridge.

### C. C++ per-process decomposition
- Any split equivalent to Morrison cloud/ice/precip process files.

Reason to defer: out of scope for phase-1 Fortran-in-place objective.

## Implementation Notes
- Keep phase-1 WSM6 simple: single-moment moisture fields only (`qv,qc,qi,qr,qs,qg`) plus accumulators.
- Reuse Morrison lifecycle shape, not Morrison number-concentration state design.
- Avoid backup/temp files (`*~`, `#*#`) when copying references from Morrison directory.

## Stop Conditions for This Context Budget
Stop and hand off once all are true:
1. Build wiring and model dispatch for WSM6 are in place.
2. WSM6 class skeleton compiles (or is one small compile-fix away with explicit note).
3. Remaining work is primarily Fortran call wiring + runtime validation.

Current evaluation:
- (1) Satisfied.
- (2) Not yet verified by build in this pass.
- (3) Expected once targeted compile is run and clean.

## Next-Agent Starting Checklist
1. Confirm stashes before modifying anything:
   - `git stash list -n 3`
2. Implement sections 1-3 above.
3. Run a targeted compile command.
4. Update this document with:
   - what compiled,
   - what remains blocked,
   - exact first failing file/line if build fails.
