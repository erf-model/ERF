# ERF-Dust Development Progress

This file records the completion status of each development phase.
Source code comments are restricted to technical descriptions.
No phase references or progress notes appear in source files.

## Phase Status

| Phase | Title | Status | Date | Notes |
|-------|-------|--------|------|-------|
| 1 | Directory scaffold and CMake integration | Complete | 2026-07-13 | All `#ifdef ERF_USE_DUST` guards in place. PR #131. |
| 2 | 2D dust grid definition | Complete | 2026-07-13 | `ERF.H` includes `ERF_Dust.H` under `#ifdef ERF_USE_DUST`. `ERF.cpp` calls `m_DustLayer->initialize()` from `InitData_post()` and `m_DustLayer->advance()` from `Advance()` under `#ifdef ERF_USE_DUST`. PR #131. |
| 3 | Surface property map reader | Complete | 2026-07-13 | Implements `ERF_DustSurfaceReader.H/cpp` with ESRI ASCII reader, MPI broadcast, and bilinear GPU interpolation. Row reversal convention matches `ERF_FuelMap.H`. Regression test added with test rasters and inputs file. |
| 4 | PHREEQC output file reader and u*t mapper | Complete | 2026-07-13 | Implements `ERF_PhreeqcReader.H/cpp` with CSV reader, rank-0 read + Bcast + GPU fill pattern. Marticorena & Bergametti (1995) u*_t reduction formula applied. Time tracking and phreeqc_update_interval_s logic added to DustLayer::advance. Regression test with PhreeqcReader test case. |
| 5 | Threshold friction velocity computation | Not started | | |
| 6 | Saltation and vertical dust emission flux | Not started | | |
| 7 | Blasting schedule emission source | Not started | | |
| 8 | Suppression agent degradation model | Not started | | |
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

## Known Limitations

None at this stage.

## Parameter Recommendations

To be populated during Phase 24 based on results from arid Western US test cases.
