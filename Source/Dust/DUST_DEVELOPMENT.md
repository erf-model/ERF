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
| 9 | Wind extraction from ERF-Atm to ERF-Dust | Complete | 2026-07-13 | PR #140. ERF_DustWindExtract.H/cpp. Adapts fill_fire_wind_from_interpolation from Source/Fire/. SurfaceLayer accessors: get_u_star(0), get_t_surf(0), get_pblh(0). Replaces test_ustar, test_surf_temp_K, test_wind_speed placeholders. |
| 10 | ERF-Dust to ERF-Atm aerosol injection | Complete | 2026-07-13 | ERF_DustAtmCoupling.H/cpp. Adapts coarsen_fire_flux_to_atm and apply_fire_tendency_to_cc_source from Source/Fire/. Single scalar slot (RhoAdv_comp). transport_bins_separately=false default. Injection at k=0 only. |
| 11 | Multi-size-bin scalar advection and Stokes settling | Complete | 2026-07-14 | ERF_DustSettling.H. Stokes settling with Cunningham correction. First-order upwind flux divergence. Per-bin diameters via bin_diameters parameter. transport_bins_separately flag controls 1 vs n_size_bins 3D scalar slots. |
| 12 | Dry deposition lower boundary condition | Complete | 2026-07-14 | ERF_DustDeposition.H. Zhang et al. (2001) resistance scheme. dep_flux_atm buffer and dust_deposition_rate accumulator. All Phase 9-12 tests updated to neutral ABL (sounding_neutral_abl, 3000x3000x1024 m, 8x8x64, u*~0.56 m/s). erf.transport_scalar=true added to all coupling tests (PR #145 fix). |
| 13 | Two-way coupling: ERF-Atm returns fields to ERF-Dust | Complete | 2026-07-14 | ERF_DustAtmReturn.H. fill_dust_conc_from_atm (Shao 2001 loading feedback). fill_dust_moist_from_atm (Fecan 1999 moisture inhibition, null-safe: Q1fx3 null when moisture_type==None). loading_feedback_coeff and use_dynamic_moisture parameters. dep_total diagnostic. Pattern follows fill_dust_scalar_from_atm (Phase 9). |
| 14 | MRF nonlocal countergradient extension to dust scalars | Complete | 2026-07-14 | Verified EddyDiff::Scalar_v set in ERF_ComputeDiffusivityMRF.cpp. dust_mrf_Sc_t parameter in TurbStruct. gamma_dust=0 confirmed. Debug print of Scalar_v_max and Theta_v_max in ERF_AdvanceDycore.cpp. DustMRFDiffusion regtest. |
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
- `erf.use_gravity = true`
- `erf.most.zref = 24.0`
- `erf.most.z0 = 0.1`
- Neutral ABL configuration from `Exec/CanonicalTests/ABL/MRF_Enhancements/canonical/neutral_abl`
- Sounding file: `sounding_neutral_abl` (copied to each test directory)
- Geostrophic wind: `15.0 0.0 0.0` (u*~0.56 m/s from MRF)
- `erf.transport_scalar = true` for all tests with `erf.dust.atm_feedback >= 1.0`

## Known Limitations

None at this stage.

## Parameter Recommendations

To be populated during Phase 24 based on results from arid Western US test cases.
