# ERF Two-Stream Radiation Module — Development Log

## Overview

The ERF Two-Stream Radiation module implements a column-based radiative
transfer solver (shortwave direct-beam Beer-Lambert + longwave gray-gas
two-stream) as a standalone, runtime-selectable physics module in
`Source/Radiation/`, decoupled from both the SLUCM urban radiation stub
(`Source/UrbanCanopy/ERF_UCMRadiationForcing.H`) and the RRTMGP coupling
path (`Source/PhysicsInterfaces/Radiation/`). It is being developed
incrementally over a 12-phase roadmap, starting from a clear-sky,
non-scattering solver and building toward full moist/cloud/aerosol/coupled
radiative heating.

**Reference scenarios:**
- Beer-Lambert direct-beam shortwave attenuation (analytical clear-sky test)
- Gray-gas two-stream longwave with isothermal analytical verification
  (Kirchhoff's law: F_up = F_down = σT⁴, zero net flux)
- Future phases: real per-column vertical integration, moist/cloud LW
  emissivity, scattering, coupling to `RhoTheta` heating tendency

---

## Phase Roadmap (12 phases, Phase 1 in progress/complete)

| Phase | Title | Key Deliverables | Status |
|-------|-------|-------------------|--------|
| 1 | Core infrastructure + clear-sky/no-scattering solver skeleton | `RadType`/`RadChoice` (ParmParse), `Source/Radiation/` skeleton, SW Beer-Lambert + LW gray-gas isothermal test, 2 RegTests | ✅ COMPLETE (with hotfixes 1b, 1c — see below) |
| 2 | Real per-column vertical integration | Replace Phase 1 analytic stub with real GPU-safe `ParallelFor` sweep over actual grid state (temperature, density) | 🟡 PLANNED |
| 3+ | Clouds, scattering, aerosols, diurnal cycle, moisture coupling, `RhoTheta` heating injection (Phase 8), PBL coupling, Immersed Forcing shadow/SVF/LW | — | ⚪ FUTURE |

**Status legend:** ✅ COMPLETE | 🟡 PLANNED (active) | ⚪ FUTURE (backlog)

---

## Phase 1 — Core Infrastructure + Clear-Sky Solver Skeleton

**Status:** ✅ COMPLETE (required three follow-up hotfixes: 1b wiring, 1c SW
bug fix, plus manual sounding/input-file corrections — see below)

### Overview

Phase 1 delivered:
1. `Source/DataStructs/ERF_RadStruct.H` — `RadType` enum (`None`,
   `TwoStream`, `RRTMGP`) and `RadChoice` struct with ParmParse-driven
   `init_params(...)`, defaulting to `RadType::None` (fully backward
   compatible / opt-in).
2. `Source/Radiation/` module skeleton: `ERF_RadiationParams.H`,
   `ERF_TwoStreamSW.H` (Beer-Lambert direct-beam math),
   `ERF_TwoStreamLW.H` (gray-gas two-stream upward/downward sweep math),
   `ERF_RadiationCoupling.H` (Phase 1 no-op placeholder;
   `RhoTheta` injection deferred to Phase 8), `ERF_RadiationDiagnostics.H`/
   `.cpp` (CSV diagnostics writer, mirrors `UCMDiagnostics` pattern),
   `Make.package`.
3. Two RegTests under `Exec/CanonicalTests/Radiation/`:
   `SW_ClearSky_Analytical/` and `LW_Isothermal/`.

### Contract #R1 (new): Radiation type is opt-in, backward compatible

`erf.radiation_type` defaults to `RadType::None`. All sub-option ParmParse
queries (`erf.radiation.sw_enabled`, `tau_per_layer`, etc.) are gated behind
`if (rad_type == RadType::TwoStream)`, so existing non-radiation tests and
runs are completely unaffected (mirrors UCM Contract #28/#37 "opt-in
preservation" pattern).

### Files touched (Phase 1 primary)

1. `Source/DataStructs/ERF_RadStruct.H` — `RadType` enum + `RadChoice`
   struct, reusing `query_one_or_per_level_enum_case_insensitive` helper
   from `ERF_TurbStruct.H`.
2. `Source/Radiation/ERF_RadiationParams.H`, `ERF_TwoStreamSW.H`,
   `ERF_TwoStreamLW.H`, `ERF_RadiationCoupling.H`,
   `ERF_RadiationDiagnostics.{H,cpp}`, `Make.package`.
3. `Exec/CanonicalTests/Radiation/{README.md, SW_ClearSky_Analytical/,
   LW_Isothermal/}`.

### Known gap discovered post-merge: build/SolverChoice wiring incomplete

**Bug:** After the primary Phase 1 PR merged, verification found that:
1. `Source/Radiation/*.cpp` was NOT registered in `CMake/BuildERFExe.cmake`
   or the GNUmake `Make.package` inclusion chain — the module compiled into
   nothing, silently excluded from the build.
2. `SolverChoice` in `Source/DataStructs/ERF_DataStruct.H` had no
   `radChoice` member, no `#include "ERF_RadStruct.H"`, and no
   `init_params(...)` call site — `ERF_RadStruct.H` was a fully correct,
   completely orphaned file.
3. `RadiationDiagnostics::append(...)` was never called anywhere in the
   codebase, so no diagnostic CSV was ever produced — the RegTests could
   never pass even with correct math.

**Root cause class:** identical to UCM Contract #13-style bugs — a
component was implemented correctly in isolation but never actually wired
into the consumer/build/call chain. "Files exist in the repo" was
mistaken for "the feature works."

### Phase 1b — Wire Radiation Module into Build System and SolverChoice

**Status:** ✅ COMPLETE.

Fixed the three gaps above:
1. Registered `Source/Radiation/ERF_RadiationDiagnostics.cpp` in the build
   (CMake `target_sources(...)` + `Make.package` inclusion), matching the
   pattern used by other physics modules (PBL, UrbanCanopy).
2. Added `#include "ERF_RadStruct.H"` and `RadChoice radChoice;` member to
   `SolverChoice`; wired `radChoice.init_params(lev, max_level, "erf")`
   alongside other physics-choice initialization.
3. Added a real call site — new file
   `Source/Radiation/ERF_AdvanceTwoStreamRadiation.cpp` implementing
   `ERF::compute_twostream_radiation_diagnostics(lev, nstep, time_step)` —
   invoked from the timestep loop, gated on
   `rad_choice.rad_type == RadType::TwoStream`, calling
   `RadiationDiagnostics::append(...)` so the CSV diagnostic file is
   actually produced.

### Contract #R2 (new): Build/wiring verification must be demonstrated, not assumed

Before considering ANY future phase complete, the PR must include actual
confirmation snippets (grep output / diff lines) showing:
(a) the new `.cpp` file appears in `target_sources(...)` and/or
`Make.package`,
(b) the new `*Choice` struct member is declared in `SolverChoice` with its
`init_params(...)` call site,
(c) any new diagnostics writer has a real call site producing output
during a run — not just a class skeleton.
"Files exist" ≠ "feature works." (Mirrors UCM's broader lesson that
compensating/silent gaps are the worst kind of bug — see
`RAD_MPI_SKILLS.md` Lesson R1.)

### Phase 1c — SW Cumulative Optical Depth Bug Fix

**Status:** ✅ COMPLETE.

**Bug:** `ERF_AdvanceTwoStreamRadiation.cpp`'s SW block computed:
```cpp
amrex::Real tau_cum = rad_choice.tau_per_layer;   // single layer, hardcoded!