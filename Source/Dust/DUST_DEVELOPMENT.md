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
| 16 | Dust diagnostics and plotfile output | Complete | 2026-07-14 | ERF_DustPlotfile.H/cpp (VisMF::Write + WriteGenericPlotfileHeader on native dust grid, pattern: ERF_FirePlotfile.cpp). ERF_DustPlotfileCatalog.H (8 fields). ERF_DustStatsOutput.H (CSV append, pattern: ERF_FireStatsOutput.H). Final-step guard: nstep > m_last_dust_plot_step. DustOutput regtest (max_step=7, dust_plot_int=3, plots at 0,3,6,7). |
| 17 | EPA NAAQS PM2.5/PM10 compliance output module | Complete | 2026-07-14 | ERF_DustPM.H (compute_pm_concentrations, update_running_average, compute_exceedance_flag). ERF_DustNAAQSOutput.H (append_naaqs_stats CSV). ERF_DustPlotfileCatalog.H updated to 14 fields. dust_naaqs.csv. DustNAAQS regtest. |
| 18 | MSHA worker exposure tracking module | Complete | 2026-07-15 | ERF_DustMSHA.H (update_msha_dose, compute_msha_exceed). ERF_DustMSHAOutput.H (append_msha_stats, write_msha_shift_summary, append_receptor_sample). Catalog updated 14->18 fields. msha_exposure.csv, msha_shift_summary.csv, msha_receptor_<name>.csv. Shift reset via floor(t/T_shift). DustMSHA regtest. |
+| 19 | Lagrangian super-particle receptor tracking | Complete | 2026-07-15 | ERF_DustPC.H/cpp (ParticleContainer<0,0,5,0>). ReleaseParticles (LoopOnCpu). AdvanceParticles (nearest-cell Euler + Stokes). dust_source_map in catalog (ncomp=19). Guarded ERF_USE_PARTICLES. DustParticles regtest. |
| 20 | Multi-mine domain with per-site PHREEQC tables | Complete | 2026-07-15 | ERF_DustSiteRegistry.H (populate_dust_site_id, count_site_cells). dust_site_id MultiFab on dust grid. Per-site PHREEQC u*_t factors. last-wins bounding box assignment. Catalog updated 19->20 fields. DustMultiSite regtest (2 sites, 32 cells each).
| 21 | PHREEQC deposition feedback file writer | Complete | 2026-07-15 | ERF_DustPHREEQCWriter.H. write_phreeqc_deposition_file (64 rows, overwrite). append_phreeqc_site_summary (per-site CSV, append). phreeqc_feedback_interval_s, phreeqc_feedback_file, phreeqc_site_summary_file params. Final-step guard m_last_phreeqc_write_step. Called from write_output(). DustPHREEQCFeedback regtest (7 steps, interval=1.0s, writes at t=1,2,3,3.5s).
| 22 | Haul road vehicle schedule emission module | Complete | 2026-07-15 | ERF_DustRoadSchedule.H. load_road_schedule (rank-0 CSV + Bcast, pattern: ERF_DustBlastSchedule.H). apply_road_schedule (GPU ParallelFor, EPA AP-42 Ch.13.2.2). RoadEvent struct, DustRoadConst namespace. road_schedule_file and road_diag_file params. Active_roads debug counter. DustRoadSchedule regtest (2 roads, time-windowed).
| 23 | DOE Critical Materials Assessment output module | Complete | 2026-07-15 | ERF_DustCriticalMaterials.H. compute_cm_flux (per-bin cm_fractions, GPU ParallelFor). append_cm_budget (domain + per-site CSV, LoopOnCpu+ReduceRealSum, pattern: ERF_FireStatsOutput.H). dust_cm_flux MultiFab. Catalog updated 20->21 fields. DustCriticalMaterials regtest (2 sites, cm_fractions=0.001).
| 24 | Regression test suite and DUST_DEVELOPMENT.md finalisation | Complete | 2026-07-15 | DustIntegration end-to-end test. All 23 individual phase tests audited and repaired. Parameter recommendations added. |

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

- Phase 15 (MRF scale-aware blending) is not implemented in ERF-Mining.
  It is implemented separately in the blending_les branch.
- Phase 19 particle advection uses nearest-cell velocity interpolation.
  Full trilinear interpolation (as in ERFPCEvolve.cpp) is a future improvement.
- The PHREEQC coupling is loose (offline file exchange) not tight (runtime).
  See ERF_DustPHREEQCWriter.H for the workflow description.
- Road emission uses the EPA AP-42 Chapter 13.2.2 formula for unpaved roads.
  Paved road emission (AP-42 Chapter 13.2.1) is not yet implemented.
- CM fractions are currently constant per bin. Time-varying fractions from
  PHREEQC mineralogy updates are a planned extension (see ERF_DustCriticalMaterials.H).
- All tests use max_level=0 (no AMR). Multi-level dust transport with AMR
  is not validated.

## Parameter Recommendations

Based on implementation design and standard atmospheric boundary layer values:

### Emission parameters
- `threshold_A_coeff = 0.0123`: Bagnold (1941) coefficient for quartz-dominated mine tailings.
- `silt_fraction = 0.05–0.20`: Typical range for mine tailings (0.10 is a reasonable default).
- `rho_air = 1.225`: Standard sea-level air density [kg/m^3].

### Settling and deposition
- `bin_diameters = 7.0e-6 2.5e-6 50.0e-6`: Covers PM2.5, PM10, and coarse fractions.
- `deposition_E0 = 3.0e-3`: Impaction efficiency for bare mine surface (Slinn 1982).
  Increase to 5e-3 for rough tailings; decrease to 1e-3 for paved surfaces.

### Suppression
- `supp_tau_base_s = 3600.0`: Water half-life ~1 h at 20°C in low wind.
  Use 7200–14400 for MgCl2 (longer-lived). Use 1800 for high temperatures (>35°C).

### PHREEQC coupling
- `phreeqc_update_interval_s = 86400.0`: Daily geochemical updates are appropriate
  for crust formation timescales (Parkhurst & Appelo 2013).

### MSHA / NAAQS output
- `msha_pel_mg_m3 = 5.0`: 30 CFR 56.5001 respirable dust PEL.
- `msha_shift_duration_s = 28800.0`: Standard 8-hour mine shift.

### Critical materials
- `cm_fractions`: Typical REE mass fractions in mine tailings are 0.001–0.01 (0.1–1%).
  Use 0.0 (empty) to disable CM tracking.
  Reference: USGS Mineral Resources Program (https://www.usgs.gov/programs/mineral-resources-program).

### Performance
- `transport_bins_separately = false`: Use for production runs unless per-bin
  concentration is required. Reduces conserved state array size.
- `atm_feedback = 0.0`: Disable atmospheric injection for standalone dust diagnostics
  without full atmosphere coupling.
