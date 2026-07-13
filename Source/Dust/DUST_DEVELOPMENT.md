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
| 5 | Threshold friction velocity computation | Complete | 2026-07-13 | Implements header-only `ERF_DustThreshold.H` with GPU device functions and MultiFab kernel. Combines Bagnold base, chemistry (crust/efflor), moisture, and suppression modifiers. Called in DustLayer::initialize and DustLayer::advance. New regression test DustThreshold with comprehensive expected values. |
| 6 | Saltation and vertical dust emission flux | Complete | 2026-07-13 | Implements header-only `ERF_DustEmission.H` and `ERF_DustEmission.cpp` with GPU device functions for horizontal and vertical emission flux. Owen (1964) saltation form with Marticorena & Bergametti (1995) sandblasting efficiency. Phase 6 uses placeholder uniform u* from `test_ustar` parameter; Phase 9 replaces with actual MRF u*. Regression test DustEmission with expected flux calculations. |
| 7 | Blasting schedule emission source | Complete | 2026-07-13 | Implements header-only `ERF_DustBlastSchedule.H` adapting `ERF_IgnitionSchedule.H` from `Source/Fire/`. Timed CSV-based mass injection into `dust_emission_flux`. Rank-0 CSV read + MPI_Bcast + GPU ParallelFor event dispatch. Added `blast_schedule_file` and `blast_reactivity` parameters to DustParams. Integration into DustLayer::initialize and DustLayer::advance. Regression test DustBlast with CSV and inputs. Documentation added to DustModule.rst. |
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

## Recent Fixes and Enhancements

### Phase 4 Linker Error Fix
- Fixed undefined symbol `update_dust_from_phreeqc` on arm64 and x86 architectures
- Root cause: `ERF_PhreeqcReader.cpp` was added to `Source/Dust/Make.package` but NOT to `CMake/BuildERFExe.cmake`
- The authoritative build source for dust module is `CMake/BuildERFExe.cmake` (not `Source/CMakeLists.txt`)
- **Solution**: Added `${SRC_DIR}/Dust/ERF_PhreeqcReader.cpp` to `target_sources()` in `CMake/BuildERFExe.cmake` dust block
- **Rule for all future phases**: Every new `.cpp` dust file must be added to BOTH `CMake/BuildERFExe.cmake` AND `Source/Dust/Make.package`
- Header-only files (like Phase 5's `ERF_DustThreshold.H`) are added only to `Make.package`, not `CMake/BuildERFExe.cmake`

### Debug Printing Enhancement
- Added exhaustive debug printing to `DustLayer::initialize()` and `DustLayer::advance()` when `erf.dust.debug = true`
- Debug output includes grid creation details, MultiFab allocations, initial values, and PHREEQC update triggers
- Added prerequisite check debug logging in `verify_dust_prerequisites()` for detailed diagnostic output
- Output format mirrors Fire module implementation with `[DUST DEBUG]` prefix

### Test Case Fixes
- **DustGrid/inputs**: Fixed `erf.most.zref = 10.0` → `50.0` (must be above first cell center); added `amr.max_level = 0`; disabled gravity (`erf.use_gravity = false`); increased `erf.most.z0 = 0.01` → `0.1` for realistic surface layer
- **DustScaffold/inputs**: Fixed `erf.most.zref = 16.0` → `100.0`; added `amr.max_level = 0`; disabled gravity; increased `erf.most.z0 = 0.1`
- **DustSurfaceReader/inputs**: Fixed `erf.most.zref = 10.0` → `50.0`; added `amr.max_level = 0`; disabled gravity; increased `erf.most.z0 = 0.01` → `0.1`; added `erf.dust.dust_debug = false`
- **PhreeqcReader/inputs**: Fixed `erf.most.zref = 10.0` → `50.0`; added `amr.max_level = 0`; disabled gravity; increased `erf.most.z0 = 0.01` → `0.1`; added `erf.dust.dust_debug = false`
- All test cases now use canonical ABL-like structure with proper MOST parameters

## Known Limitations

None at this stage.

## Parameter Recommendations

To be populated during Phase 24 based on results from arid Western US test cases.
