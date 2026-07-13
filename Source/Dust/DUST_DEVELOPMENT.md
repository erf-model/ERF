# ERF-Dust Development Progress

This file records the completion status of each development phase.
Source code comments are restricted to technical descriptions.
No phase references or progress notes appear in source files.

## Phase Status

| Phase | Title | Status | Date | Notes |
|-------|-------|--------|------|-------|
| 1 | Directory scaffold and CMake integration | Complete | 2026-07-13 | PR #128, #131. |
| 2 | 2D dust grid definition | Complete | 2026-07-13 | PR #129, #130, #131. |
| 3 | Surface property map reader | Complete | 2026-07-13 | PR #132. |
| 4 | PHREEQC output file reader and u*t mapper | Complete | 2026-07-13 | PR #134. Linker fix: PR #136. |
| 5 | Threshold friction velocity computation | Complete | 2026-07-13 | PR #136. |
| 6 | Saltation and vertical dust emission flux | Complete | 2026-07-13 | PR #137. |
| 7 | Blasting schedule emission source | Complete | 2026-07-13 | PR #138. |
| 8 | Suppression agent degradation model | Complete | 2026-07-13 | PR #139. ERF_DustSuppression.H. Retreat flag. amrex.fpe_trap_invalid=0 applied to all tests. AMReX_ParallelFor.H bug fixed (commit 6054023). |
| 9 | Wind extraction from ERF-Atm to ERF-Dust | Not started | | |
| 10 | ERF-Dust to ERF-Atm aerosol injection | Not started | | |
| 11 | Multi-size-bin scalar advection and Stokes settling | Not started | | |
| 12 | Dry deposition lower boundary condition | Not started | | |
| 13 | Two-way coupling: ERF-Atm returns fields to ERF-Dust | Not started | | |
| 14 | MRF nonlocal countergradient extension to dust scalars | Not started | | |
| 15 | MRF high-resolution diffusivity bounds at fine AMR levels | Not started | | |
| 16 | Dust diagnostics and plotfile output | Not started | | |
| 17 | EPA NAAQS PM2.5/PM10 compliance output module | Not started | | |
| 18 | MSHA worker exposure tracking module | Not started | | |
| 19 | Lagrangian super-particle receptor tracking | Not started | | |
| 20 | Multi-mine domain with per-site PHREEQC tables | Not started | | |
| 21 | PHREEQC deposition feedback file writer | Not started | | |
| 22 | Haul road vehicle schedule emission module | Not started | | |
| 23 | DOE Critical Materials Assessment output module | Not started | | |
| 24 | Regression test suite and DUST_DEVELOPMENT.md finalisation | Not started | | |

## Build Rules (apply to every phase)

- New `.cpp` file: add to BOTH `CMake/BuildERFExe.cmake` dust block AND `Source/Dust/Make.package`.
- Header-only file: add to `Make.package` only. Never add to `Source/CMakeLists.txt`.
- Never `#include <AMReX_ParallelFor.H>` in any dust file (not even commented out). Use `<AMReX_MFIter.H>`.

## Test Input Requirements (apply to every test)

All files under `Exec/CanonicalTests/Dust/*/inputs` must have:
- `amrex.fpe_trap_invalid = 0`
- `erf.dust.dust_debug = true`
- `amr.max_level = 0`
- `erf.use_gravity = false`
- `erf.most.zref = 50.0`
- `erf.most.z0 = 0.1`

## Known Limitations

None at this stage.

## Parameter Recommendations

To be populated during Phase 24 based on results from arid Western US test cases.
