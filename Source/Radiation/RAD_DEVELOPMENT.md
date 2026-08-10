# Radiation Development Roadmap & Phase History

This document tracks the development of the two-stream radiation model through phases, including contracts, architectural decisions, known issues, and fixes.

## Roadmap

### Roadmap Status Policy (as of 2026-08-10)

- **Complete**: implemented and validated against designated regression/benchmark cases.
- **Planned (Active)**: in current execution scope and expected to be implemented in sequence.
- **Shelved**: retained for future reference but out of current execution scope.

---

| Phase | Name | Status | PR | Key Feature | Ease | RegTest |
|---|---|---|---|---|---|---|
| **1** | Two-Stream Skeleton with Analytic Stub | ✅ Complete | N/A | Clear-sky SW/LW, diagnostic output, single-layer optical depth | Easy | `SW_ClearSky_Analytical`, `LW_Isothermal` |
| **2** | Real Per-Column Two-Stream Radiation | ✅ Complete | #283 | Per-column vertical integration, actual grid bounds, GPU-safe kernel | Moderate | `SW_ClearSky_Analytical`, `LW_Isothermal` |
| **3** | Cloud Optical Properties | ✅ Complete | N/A (manual) | Height-varying cloud-layer optical depth, cloud fraction masking | Easy | `SW_Cloud_Layer` (+ Phase 1–2 regressions) |
| **4** | Scattering Effects | ✅ Complete | Merged | Diffuse SW scattering via Meador-Weaver two-stream approximation | Moderate | `SW_Scattering_Cloud` (+ Phase 1–3 regressions) |
| **5** | RhoTheta Coupling | ✅ Complete | N/A (manual) | Per-level SW/LW heating written to `qheating_rates` and injected into `RhoTheta` | Moderate | `Phase5_RhoTheta_Coupling` (+ Phase 1–4 regr[...]
| **6** | Time-Stepping Integration | ✅ Complete | TBD | TwoStream call-cadence and temporal consistency with slow-step/source application + call_site diagnostics | Moderate | `Phase6_TimeIntegratio[...]
| **7** | TwoStream Runtime Diagnostics Controls | ✅ Complete | TBD | Runtime diagnostics controls (mode/streams/schema toggles), no physics change | Easy | `TwoStream_DiagControls` |
| **8** | Validation & Benchmarking | ✅ Complete | TBD | Canonical benchmark suite + automated metric checks | Moderate | `Radiation_Benchmark_Suite` |
| **9** | TwoStream Integration Polish I | ✅ Complete | TBD | Cadence/de-dup hardening + nonuniform-`dz` heating framework + finite guards | Easy | `TwoStream_Cadence_NonuniformDZ` |
| **10** | True Nonuniform `dz(k)` Wiring | ✅ Complete | TBD | Wire per-level `dz(k)` from physical vertical geometry (`z_phys_cc`) with uniform fallback retained | Moderate | `TwoStream_NonuniformD[...]
| **11** | Surface Heterogeneity + Fallback (Albedo/Emissivity/`t_sfc`) | ✅ Complete | TBD | TwoStream consumes per-column LSM/Radiation surface fields with robust fallback path | Moderate | `TwoStr[...]
| **12** | Moisture/Cloud-Aware Dynamic Optical Depth | ✅ Complete | TBD | Diagnose SW/LW `tau(k)` from `qv`, `qc`, `rho`, `dz` with safe fallback | Moderate | `TwoStream_DynamicTau_MoistCloud` |
| **13** | PBL Coupling Focus (YSUNew-only) | ✅ Complete | TBD | YSUNew radiative tendency smoothing/limiter + diagnostic hooks; MRF deferred | Moderate | `TwoStream_PBL_MRF_YSU_Coupling` |
| **14** | Prognostic Cloud Fraction for Radiation | ✅ Complete | TBD | RH/`qc`-based diagnosed cloud fraction with bounds and temporal smoothing | Easy–Moderate | `TwoStream_ProgCloudFraction` |
| **14A** | LSM Surface Properties Wiring + TwoStream Bugfixes | ✅ Complete | TBD | Wire LSM surface fields (albedo/emissivity/t_sfc) into TwoStream with standalone fallback MultiFabs; fix cloudy-column dead call; fix pre_dycore time arg | Easy–Moderate | `TwoStream_ProgCloudFraction` (extended) |
| **15** | Bulk Aerosol/Turbidity Option | ✅ Complete | N/A | Prescribed aerosol optical-depth profile (constant/exponential/table), optional LW hook | Easy–Moderate | `TwoStream_Aerosol_Turbidity` |
| **16** | Time-Varying Solar Geometry | ✅ Complete | TBD | Solar zenith evolution with time/lat/day; fixed-angle fallback retained | Easy | `TwoStream_DiurnalSolarGeometry` |
| **17** | Simplified Surface Energy Balance (SEB) — Diagnostic Mode | ⏳ Planned (Active) | TBD | Compute/report SEB residual terms from TwoStream + surface inputs (no prognostic `T_s` update) | Moderate | `TwoStream_SEB_Diagnostic` |
| **18** | Simplified SEB — Prognostic `T_s` Mode | ⏳ Planned (Active) | TBD | Optional explicit `T_s` tendency update with limiter/clamps and fallback-safe behavior | Moderate | `TwoStream_SEB_PrognosticTs` |
| **19** | SEB Coupling Safeguards (Noah-MP/SurfaceLayer Interop) | ⏳ Planned (Active) | TBD | Anti-double-count rules and precedence guards for `T_s`, `H`, `LE`, and radiative terms | Moderate | `TwoStream_SEB_InteropGuards` |
| **20** | SEB Validation & Benchmark Suite | ⏳ Planned (Active) | TBD | Canonical SEB closure/stability tests, tolerances, and CI-ready reports | Moderate | `TwoStream_SEB_BenchmarkSuite` |
| **14A** | LSM Surface Properties Wiring + TwoStream Bugfixes | ✅ Complete | TBD | Wire LSM surface fields (albedo/emissivity/t_sfc) into TwoStream with standalone fallback MultiFabs; fix cloudy-co[...]
| **15** | Bulk Aerosol/Turbidity Option | ✅ Complete | N/A | Prescribed aerosol optical-depth profile (constant/exponential/table), optional LW hook | Easy–Moderate | `TwoStream_Aerosol_Turbidity[...]
| **16** | Time-Varying Solar Geometry | ⏳ Planned (Active) | TBD | Solar zenith evolution with time/lat/day; fixed-angle fallback retained | Easy | `TwoStream_DiurnalSolarGeometry` |
| **17** | Simplified SEB — MultiFab Infrastructure + Noah-MP Passthrough | ⏳ Planned (Active) | TBD | Create all SEB-related MultiFabs (albedo, emissivity, SW, LW, sensible heat, latent heat, ground heat flux, ground/surface temperature, ground/surface moisture, deep soil temperature, deep soil moisture); wire logic to copy values from Noah-MP when active; simplified/standalone SEB computation deferred to Phase 18 | Moderate | `TwoStream_SEB_MultiFabInfra` |
| **18** | Simplified SEB — Diagnostic Mode | ⏳ Planned (Active) | TBD | Compute/report SEB residual terms from TwoStream + surface inputs using Phase 17 MultiFabs (no prognostic `T_s` update) | Moderate | `TwoStream_SEB_Diagnostic` |
| **19** | Simplified SEB — Prognostic Mode | ⏳ Planned (Active) | TBD | Optional explicit `T_s` (and surface moisture) tendency update with limiter/clamps and fallback-safe behavior | Moderate | `TwoStream_SEB_Prognostic` |
| **20** | SEB Coupling Safeguards (Noah-MP/SurfaceLayer Interop) | ⏳ Planned (Active) | TBD | Anti-double-count rules and precedence guards for `T_s`, `H`, `LE`, and radiative terms between simplified SEB and Noah-MP/SurfaceLayer | Moderate | `TwoStream_SEB_CouplingSafeguards` |

**Note (2026-08-10)**: SEB Validation & Benchmark Suite is deferred/deprioritized for now (not currently scheduled as a numbered phase). It may be reinstated as a future phase once Phases 17–20 are complete and stable.

---

## Phase 17 Implementation (Simplified SEB — MultiFab Infrastructure + Noah-MP Passthrough)

**Status**: ⏳ Planned (Active) — not yet implemented  
**Scope**: TwoStream radiation + new SEB surface-field infrastructure; Noah-MP passthrough wiring; no new standalone SEB physics computation yet  
**Key Feature**: Allocate all MultiFabs required for simplified Surface Energy Balance, and wire precedence logic so that when Noah-MP (or another active LSM) provides a field, its value is copied/used; otherwise the field is populated by a placeholder/fallback pending the real simplified SEB solver in Phase 18.

### Planned Fields (per-level 2D MultiFabs, following the `twostream_alb_sw`/`ba2d[lev]` pattern from Phase 14A)

1. `sfc_alb_sw` — surface shortwave albedo
2. `sfc_emis_lw` — surface longwave emissivity
3. `sw_flux_sfc` — net/incident shortwave flux at surface
4. `lw_flux_sfc` — net/incident longwave flux at surface
5. `hfx_sfc` — sensible heat flux (H)
6. `lh_sfc` — latent heat flux (LE)
7. `grdflx_sfc` — ground heat flux (G)
8. `t_sfc` — surface (skin) temperature (already exists from Phase 11/14A; reused, not duplicated)
9. `q_sfc` — surface moisture (specific humidity / moisture availability, `beta`)
10. `t_deep` — deep soil temperature (force-restore reservoir)
11. `q_deep` — deep soil moisture (force-restore reservoir; optional/deferred activation, see below)

### Technical Design (planned)

1. **Allocation**: All fields allocated in `ERF_MakeNewArrays.cpp`, same step/pattern as `qheating_rates`/`twostream_alb_sw`, gated on `solverChoice.radChoice.rad_type == RadType::TwoStream` and a new master switch (see below). Vector members added to `ERF.H` and resized in `ERF_Constructors.cpp::ERF_shared()` (per the Phase 14B lesson — resize immediately alongside existing `twostream_*` vectors to avoid the same undefined-behavior bug).

2. **Precedence / passthrough logic** (mirrors Phase 14A's LSM wiring pattern exactly):
   - If Noah-MP (or other active LSM) exposes a corresponding named field (e.g., `hfx3`, `q1fx3`/latent heat equivalent, `grdflx`, `tg`/`t_sfc`, soil moisture/temperature layers), copy/reference that value into the corresponding SEB MultiFab each step, using `lsm.Get_DataIdx()` / `lsm.Get_Data_Ptr()` as in Phase 14A.
   - If no LSM is active, or a specific field is unavailable, fall back to a constant-filled placeholder (from new `RadChoice` scalar defaults), matching the Phase 14A standalone-MultiFab fallback pattern. No new physics is computed yet in Phase 17 — this phase is infrastructure + passthrough only.

3. **New RadChoice parameters (planned)**:
   - `seb_enable` [bool]: master switch for allocating/using SEB infrastructure; default `false` (fully backward compatible, no new MultiFabs allocated when off).
   - Scalar fallback defaults for each new field (e.g., `seb_hfx_default`, `seb_lh_default`, `seb_grdflx_default`, `seb_q_sfc_default`, `seb_t_deep_default`, `seb_q_deep_default`), all validated/clamped in `init_params()`.

4. **`q_deep` activation**: Per prior design discussion, `q_deep` and a prognostic moisture-reservoir update are optional/deferred — Phase 17 allocates the MultiFab for forward compatibility, but a constant/prescribed `q_sfc`/`beta` is sufficient initially; full force-restore moisture coupling may be added in a future phase if needed.

### Backward Compatibility (required)

- `seb_enable = false` (default): no new MultiFabs allocated, zero behavior change, bitwise-identical to Phase 15/16 baseline.
- When enabled with no active LSM: all fields populated from constant scalar fallbacks; no physics change to existing TwoStream radiation output.
- When enabled with active Noah-MP: fields mirror Noah-MP values; no double-counting yet since no independent SEB solve occurs in this phase (deferred to Phase 18/19, with anti-double-count safeguards formalized in Phase 20).

### GPU Safety (required)
- All allocation/copy logic follows existing `AMREX_GPU_DEVICE`/host-side MFIter patterns established in Phase 14A; no new host-side I/O in device kernels.

### Deferred to later phases
- Actual simplified SEB residual computation (`R_net - H - LE - G`): **Phase 18**.
- Prognostic `T_s` (and surface moisture) update: **Phase 19**.
- Anti-double-count precedence/coupling safeguards vs. Noah-MP/SurfaceLayer: **Phase 20**.
- Full validation/benchmark suite: deferred/deprioritized (see roadmap note above).

---

## Phase 11 Fallback Contract (Surface Heterogeneity)

Per-column surface properties for TwoStream must resolve in this priority order:

1. **Use heterogeneous model field** (LSM/Radiation interface) when available and finite:
   - `t_sfc(i,j)`
   - `sfc_emis(i,j)`
   - SW albedo inputs (`sfc_alb_dir_vis/nir`, `sfc_alb_dif_vis/nir`) or mapped broadband equivalents
2. **Else use configured scalar fallback** from runtime parameters (`RadChoice`).
3. **Else use hard-safe default** (documented constants), with warn-once behavior.

Required guards:
- Clamp albedo/emissivity to `[0,1]`.
- Enforce finite values for all fallback-resolved fields.
- Preserve backward-compatible behavior when heterogeneous fields are absent.
- Fallback logic must not abort simulation solely due to missing LSM fields.

---

## Phase 16 Implementation (Time-Varying Solar Geometry)

**Status**: ✅ Complete (as of 2026-08-10)  
**Replaces**: Phase 15 fixed solar zenith angle  
**Key Feature**: Dynamic solar position computation from astronomical formulas

### Implementation Summary

Phase 16 extends TwoStream to compute solar zenith angle dynamically from simulation time, site latitude, longitude, and day-of-year using standard astronomical formulas. The fixed-angle path (static `solar_zenith_deg` scalar) is **fully retained** as the default fallback behavior — this phase is purely additive via the `solar_geometry_dynamic_enable` master switch (default `false` for backward compatibility).

#### Technical Design

1. **Solar Geometry Module** (`ERF_SolarGeometry.H`):
   - All functions are GPU-safe: `AMREX_GPU_DEVICE AMREX_FORCE_INLINE`
   - Helper functions compute:
     - **Solar declination** from day-of-year (Spencer's Fourier formula, ~0.006 rad accuracy)
     - **Equation of time** correction for solar time (~16 minute range)
     - **Solar hour angle** from UTC time, longitude, and time zone offset
     - **Solar zenith angle** from declination, hour angle, latitude (standard spherical trig formula)
     - **Solar azimuth angle** (for future building-shadow/orientation work)
   - Input guards against invalid/out-of-range values (NaN/Inf/latitude outside [-90,90], longitude outside [-180,180], doy outside [1,366]) with safe fallbacks and clamping consistent with existing `clamp_finite()` conventions

2. **New RadChoice Parameters** (in `ERF_RadStruct.H`):
   - `solar_geometry_dynamic_enable` [bool]: master switch, default `false` (preserves Phase 15 behavior exactly)
   - `latitude_deg` [real]: site latitude [-90,90], default 0.0, clamped in `init_params()`
   - `longitude_deg` [real]: site longitude [-180,180], default 0.0, clamped/wrapped in `init_params()`
   - `day_of_year` [real/int]: reference day-of-year at simulation start [1,366], default 172.0 (summer solstice), clamped in `init_params()`
   - `time_zone_offset_hours` [real]: time zone offset from UTC [hours], default 0.0 (UTC), ensures finite in `init_params()`
   - All parameters validated/clamped in `init_params()` matching Phase 12-15 pattern

3. **Function Changes**:
   - **`vertical_two_stream_sweep()`**: Added `time_utc_seconds` parameter (default 0.0) to accept simulation time
     - When `solar_geometry_dynamic_enable == false` (default): `cos_zenith` computation is **bitwise-identical** to Phase 15 (uses static `solar_zenith_deg`)
     - When enabled: computes `cos_zenith` from solar geometry helpers using current simulation time
   - **`compute_twostream_radiation_diagnostics()`**: Converts absolute simulation time `t_old[lev]` to UTC seconds within day via `std::fmod(t_old[lev], 86400.0)`, then computes dynamic `cos_zenith` for `SW_TOA` and passes time to `vertical_two_stream_sweep()` calls
   - Both `sw_surface_flux` and `SW_TOA` diagnostics reflect dynamic zenith angle when enabled

#### New RadChoice Parameters Table

| Parameter | Type | Default | Range | Phase | Description |
|-----------|------|---------|-------|-------|-------------|
| `solar_geometry_dynamic_enable` | bool | false | bool | 16 | Enable time-varying solar position computation; false → Phase 15 behavior (bitwise-identical) |
| `latitude_deg` | real | 0.0 | [-90, 90] | 16 | Site latitude (degrees, positive north) |
| `longitude_deg` | real | 0.0 | [-180, 180] | 16 | Site longitude (degrees, positive east, wrapped/clamped) |
| `day_of_year` | real | 172.0 | [1, 366] | 16 | Reference day-of-year at sim start (summer solstice default) |
| `time_zone_offset_hours` | real | 0.0 | unconstrained | 16 | Time zone offset from UTC (hours); ensures finite |

#### Backward Compatibility

- **Disabled by default**: `solar_geometry_dynamic_enable=false` preserves Phase 15 behavior exactly
- **Bitwise-identical path**: When disabled, `cos_zenith` computation is **exactly** `std::cos(solar_zenith_deg * M_PI / 180.0)` — no floating-point differences
- **No changes to other parameters**: All existing RadChoice defaults remain unchanged (static `solar_zenith_deg`, etc.)
- **Existing RegTests**: `TwoStream_Aerosol_Turbidity`, `TwoStream_ProgCloudFraction`, etc. continue to pass unmodified (they don't enable the feature)

#### GPU Safety

- All new solar-geometry helpers: `AMREX_GPU_DEVICE AMREX_FORCE_INLINE`
- No host-side I/O, dynamic allocation, or thread-local state inside device code
- All trig/transcendental math uses `std::` functions (e.g., `std::cos`, `std::sin`, `std::acos`, `std::atan2`), matching existing zenith-angle computation
- Input validation uses `std::isfinite()` guards on all inputs, matching `clamp_finite()` conventions already established
- Time modulo arithmetic (`std::fmod()`) is GPU-safe

#### Integration Points & Example Configuration

```ini
# Enable Phase 16 dynamic solar geometry
erf.radiation.solar_geometry_dynamic_enable = true
erf.radiation.latitude_deg = 45.0           # Site at 45°N (e.g., Minneapolis/Seattle)
erf.radiation.longitude_deg = -93.0         # Western hemisphere (central US: -93°)
erf.radiation.day_of_year = 172.0           # June 21 (summer solstice)
erf.radiation.time_zone_offset_hours = -6.0 # Central Daylight Time (UTC-6)

# Static solar zenith still required (used when dynamic disabled):
erf.radiation.solar_zenith = 60.0            # Fallback/reference zenith angle [degrees]
```

**Behavior**:
- Simulation time `t_old[lev]` is converted to UTC seconds within a day via modulo 86400 s
- Assuming simulation starts at 00:00 UTC on the specified `day_of_year`
- For latitude 45°N, longitude -93°, day 172 (June 21), time zone UTC-6:
  - Local solar noon is ~12:45 UTC (due to longitude + equation-of-time correction)
  - Sunrise ~05:45 UTC, sunset ~21:00 UTC (partial range; full day shown in test)
  - `cos_zenith` ranges from ~-0.5 (night) to ~+0.5 (peak near noon at 45°N in summer)

#### New RegTest: `TwoStream_DiurnalSolarGeometry`

Located at `Exec/CanonicalTests/Radiation/TwoStream_DiurnalSolarGeometry/`, with:

- **`inputs_baseline`**: `solar_geometry_dynamic_enable = false`
  - Validates **backward compatibility**: `cos_zenith` and fluxes match static `solar_zenith_deg=60.0` baseline
  - Confirms no regressions in existing TwoStream functionality
  
- **`inputs_dynamic`**: `solar_geometry_dynamic_enable = true`
  - Latitude 45°N, longitude -93°, day 172 (June 21), time zone UTC-6
  - Full 24-hour simulation (86400 s) to observe complete diurnal cycle
  - Validates:
    - `cos_zenith` varies plausibly across the day (lowest at sunrise/sunset, highest at noon)
    - `SW_surface = 0` when sun is below horizon (`cos_zenith <= 0`)
    - `SW_surface` increases as sun rises higher (toward local solar noon)
    - No NaN/Inf in any diagnostic output across full day
    - `SW_TOA` peaks near local solar noon (around 12:45 UTC for this site)

- **`check_solar.py`**: Python validation script
  - Reads diagnostics CSV files from both baseline and dynamic cases
  - Checks finite values, peak magnitudes, and diurnal pattern plausibility
  - Reports pass/fail with detailed error messages

#### Future Extensions

- **Building/urban-canopy shadowing** (deferred to Phase 18+): azimuth output from this phase is preparatory infrastructure only
- **Multi-band solar geometry** (deferred): current implementation assumes broadband solar constant `S0`; future work may extend to per-band declination corrections
- **Simplified Surface Energy Balance** (Phases 17–19, unaffected by this work): SEB fluxes will use dynamic `cos_zenith` when enabled, naturally improving surface temperature feedback

#### Verification & Validation Checklist

- [x] Solar geometry module created with all required helper functions (declination, equation-of-time, hour-angle, zenith, azimuth)
- [x] All helpers are GPU-safe (`AMREX_GPU_DEVICE AMREX_FORCE_INLINE`)
- [x] Input validation guards against NaN/Inf/out-of-range with safe fallbacks
- [x] RadChoice parameters added with defaults and init_params() validation
- [x] vertical_two_stream_sweep() accepts time parameter and computes dynamic cos_zenith when enabled
- [x] compute_twostream_radiation_diagnostics() converts t_old[lev] to UTC seconds and passes to sweep function
- [x] Backward-compatible path (disabled by default) is bitwise-identical to Phase 15
- [x] TwoStream_DiurnalSolarGeometry RegTest directory created with baseline and dynamic inputs
- [x] Validation script (check_solar.py) created to verify finite values, backward compat, and diurnal pattern
- [x] RAD_DEVELOPMENT.md documentation updated with full Phase 16 section
- [x] Code builds cleanly with existing CMake (no new dependencies)
- [x] Existing RegTests continue to pass when feature is disabled (default)
- [x] New RegTest validates both backward-compat and feature-on behavior

---

## References

- Toon et al., 1989: "Rapid calculation of radiative heating rates...", *J. Geophys. Res.*, 94, 16387–16405.
- Beer, A., 1852: "Bestimmung der Absorption des rothen Lichts in farbigen Flüssigkeiten", *Ann. Phys. Chem.*, 86, 78–88.
- Meador, W. E., and W. R. Weaver, 1980: "Two-stream approximations to radiative transfer in planetary atmospheres: A unified description of existing methods and a new improvement", *J. Atmos. Sci.*
- Spencer, J. W., 1971: "Fourier series representation of the position of the sun", *Search*, 2(5), 172–172.
- Duffie, J. A., and W. A. Beckman, 1991: *Solar Engineering of Thermal Processes*, John Wiley & Sons.
- AMReX GPU Guide: https://amrex-codes.github.io/amrex/docs_html/GPU.html
