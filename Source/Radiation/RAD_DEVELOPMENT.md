# Radiation Development Roadmap & Phase History

This document tracks the development of the two-stream radiation model through phases, including contracts, architectural decisions, known issues, and fixes.

## Roadmap

### Roadmap Status Policy (as of 2026-08-08)

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
| **5** | RhoTheta Coupling | ✅ Complete | N/A (manual) | Per-level SW/LW heating written to `qheating_rates` and injected into `RhoTheta` | Moderate | `Phase5_RhoTheta_Coupling` (+ Phase 1–4 regressions) |
| **6** | Time-Stepping Integration | ✅ Complete | TBD | TwoStream call-cadence and temporal consistency with slow-step/source application + call_site diagnostics | Moderate | `Phase6_TimeIntegration` |
| **7** | TwoStream Runtime Diagnostics Controls | ✅ Complete | TBD | Runtime diagnostics controls (mode/streams/schema toggles), no physics change | Easy | `TwoStream_DiagControls` |
| **8** | Validation & Benchmarking | ✅ Complete | TBD | Canonical benchmark suite + automated metric checks | Moderate | `Radiation_Benchmark_Suite` |
| **9** | TwoStream Integration Polish I | ✅ Complete | TBD | Cadence/de-dup hardening + nonuniform-`dz` heating framework + finite guards | Easy | `TwoStream_Cadence_NonuniformDZ` |
| **10** | True Nonuniform `dz(k)` Wiring | ✅ Complete | TBD | Wire per-level `dz(k)` from physical vertical geometry (`z_phys_cc`) with uniform fallback retained | Moderate | `TwoStream_NonuniformDZ` |
| **11** | Surface Heterogeneity + Fallback (Albedo/Emissivity/`t_sfc`) | ✅ Complete | TBD | TwoStream consumes per-column LSM/Radiation surface fields with robust fallback path | Moderate | `TwoStream_SurfaceHeterogeneity` |
| **12** | Moisture/Cloud-Aware Dynamic Optical Depth | ✅ Complete | TBD | Diagnose SW/LW `tau(k)` from `qv`, `qc`, `rho`, `dz` with safe fallback | Moderate | `TwoStream_DynamicTau_MoistCloud` |
| **13** | PBL Coupling Focus (YSUNew-only) | ✅ Complete | TBD | YSUNew radiative tendency smoothing/limiter + diagnostic hooks; MRF deferred | Moderate | `TwoStream_PBL_MRF_YSU_Coupling` |
| **14** | Prognostic Cloud Fraction for Radiation | ⏳ Planned (Active) | TBD | RH/`qc`-based diagnosed cloud fraction with bounds and temporal smoothing | Easy–Moderate | `TwoStream_ProgCloudFraction` |
| **15** | Bulk Aerosol/Turbidity Option | ⏳ Planned (Active) | TBD | Prescribed aerosol optical-depth profile (constant/exponential/table), optional LW hook | Easy–Moderate | `TwoStream_Aerosol_Turbidity` |
| **16** | Time-Varying Solar Geometry | ⏳ Planned (Active) | TBD | Solar zenith evolution with time/lat/day; fixed-angle fallback retained | Easy | `TwoStream_DiurnalSolarGeometry` |
| **17** | Simplified Surface Energy Balance (SEB) — Diagnostic Mode | ⏳ Planned (Active) | TBD | Compute/report SEB residual terms from TwoStream + surface inputs (no prognostic `T_s` update) | Moderate | `TwoStream_SEB_Diagnostic` |
| **18** | Simplified SEB — Prognostic `T_s` Mode | ⏳ Planned (Active) | TBD | Optional explicit `T_s` tendency update with limiter/clamps and fallback-safe behavior | Moderate | `TwoStream_SEB_PrognosticTs` |
| **19** | SEB Coupling Safeguards (Noah-MP/SurfaceLayer Interop) | ⏳ Planned (Active) | TBD | Anti-double-count rules and precedence guards for `T_s`, `H`, `LE`, and radiative terms | Moderate | `TwoStream_SEB_InteropGuards` |
| **20** | SEB Validation & Benchmark Suite | ⏳ Planned (Active) | TBD | Canonical SEB closure/stability tests, tolerances, and CI-ready reports | Moderate | `TwoStream_SEB_BenchmarkSuite` |

---

## Phase 12 Implementation (Moisture/Cloud-Aware Dynamic Optical Depth)

**Status**: ✅ Complete (as of 2026-08-08)  
**Replaces**: Phase 11 static optical depth  
**Key Feature**: Per-level dynamic optical depth diagnosis from atmospheric moisture and cloud content

### Implementation Summary

Phase 12 extends TwoStream to compute per-level optical depth dynamically from available atmospheric state fields (water vapor qv, cloud liquid qc). Instead of using a uniform `tau_per_layer` / `tau_lw_per_layer` for all levels, Phase 12 diagnoses optical depth at each level as:

```
tau_sw(k) = tau_per_layer + tau_sw_coeff_qv * qv(k) + tau_sw_coeff_qc * qc(k)
tau_lw(k) = tau_lw_per_layer + tau_lw_coeff_qv * qv(k) + tau_lw_coeff_qc * qc(k)
```

This allows optical depth to vary with height based on the moisture/cloud profile, enabling more realistic absorption by water vapor and cloud droplets.

#### Technical Design

1. **Dynamic Tau Diagnosis Functions** (GPU-safe device kernels):
   - `diagnose_tau_sw_dynamic()`: Computes SW optical depth from qv/qc at level (i,j,k)
   - `diagnose_tau_lw_dynamic()`: Computes LW optical depth from qv/qc at level (i,j,k)
   - Both retrieve mixing ratios (qv, qc) from state array using safe division
   - Both guard against invalid values (NaN, Inf, negative) by using 0 fallback
   - Both clamp output to physically reasonable range [0, 100]

2. **New RadChoice Parameters** (Phase 12):
   - `tau_sw_dynamic_enable` [bool]: Master switch for dynamic SW tau (default false)
   - `tau_lw_dynamic_enable` [bool]: Master switch for dynamic LW tau (default false)
   - `tau_sw_coeff_qv` [real]: SW absorption coefficient for qv (default 0.0)
   - `tau_sw_coeff_qc` [real]: SW absorption coefficient for qc (default 0.0)
   - `tau_lw_coeff_qv` [real]: LW absorption coefficient for qv (default 0.0)
   - `tau_lw_coeff_qc` [real]: LW absorption coefficient for qc (default 0.0)
   - All coefficients clamped to nonnegative by init_params() validation

3. **Function Changes**:
   - **ERF_RadStruct.H**: Added Phase 12 parameters + init_params queries
   - **ERF_AdvanceTwoStreamRadiation.cpp**:
     - New device-inline functions: `diagnose_tau_sw_dynamic()`, `diagnose_tau_lw_dynamic()`
     - SW downward sweep: After computing static `tau_sw`, call dynamic diagnosis if enabled
     - LW upward sweep: After computing static `tau_lw`, call dynamic diagnosis if enabled
     - LW downward sweep: After computing static `tau_lw`, call dynamic diagnosis if enabled

#### Backward Compatibility

- **Disabled by default**: `tau_sw_dynamic_enable=false` and `tau_lw_dynamic_enable=false` preserve Phase 11 behavior
- **Zero coefficients**: When all coefficients are 0.0 (default), dynamic path returns static tau unchanged
- **Numerical path**: When disabled or with zero coefficients, kernel path is **bitwise-identical** to Phase 11
- **Missing fields**: If qv or qc unavailable/invalid (outside domain, NaN, Inf), code silently uses 0 and falls back to static tau
- **Fallback safety**: No crash or exception; simulation continues with static tau

#### GPU Safety

- All new helpers: `AMREX_GPU_DEVICE AMREX_FORCE_INLINE`
- No host-side I/O in device code
- No dynamic allocations or thread-local state
- Safe mixing ratio extraction: `qv = rho_qv / rho` with guard `if (rho > 0.0)`
- Array bounds checked via `Array4::contains()` before access
- Respects existing AMReX patterns (lambda captures, inline functions)

#### Integration Points

1. **Inputs File** (Phase 12 RegTest):
   ```
   # Disabled by default (Phase 11 compat)
   erf.radiation.tau_sw_dynamic_enable = false
   erf.radiation.tau_lw_dynamic_enable = false
   erf.radiation.tau_sw_coeff_qv = 0.0
   erf.radiation.tau_sw_coeff_qc = 0.0
   erf.radiation.tau_lw_coeff_qv = 0.0
   erf.radiation.tau_lw_coeff_qc = 0.0
   
   # To enable dynamic tau, set:
   # erf.radiation.tau_sw_dynamic_enable = true
   # erf.radiation.tau_lw_dynamic_enable = true
   # erf.radiation.tau_sw_coeff_qv = <value>   # e.g., 10.0
   # erf.radiation.tau_sw_coeff_qc = <value>   # e.g., 100.0
   # erf.radiation.tau_lw_coeff_qv = <value>   # e.g., 20.0
   # erf.radiation.tau_lw_coeff_qc = <value>   # e.g., 200.0
   ```

2. **Diagnostics Output**: Unchanged; heating rates and surface fluxes remain domain-averaged scalars

3. **Future Extensions** (Phase 13+):
   - Spectral band differentiation (SW vs. LW use different coefficients)
   - Cloud optical properties (liquid, ice, different sizes)
   - Aerosol optical depth coupling
   - Prognostic/diagnostic cloud fraction modulation

#### Verification & Validation

1. **Compile Check**: ✅ No new dependencies, minimal code changes (3 device functions, 6 ParmParse queries)
2. **Fallback Safety**: ✅ Invalid/missing qv/qc silently use 0; coefficients guard against negative values
3. **Finite Guards**: ✅ NaN/Inf in qv/qc caught and replaced; output clamped to [0, 100]
4. **Phase 11 Regression**: ✅ Existing tests with dynamic tau disabled show bitwise-identical output
5. **Phase 12 RegTest**: ✅ New test exercises both static (backward compat) and dynamic (qv/qc-dependent) paths

---

## Phase 13 Implementation (PBL Coupling Focus: YSUNew-only)

**Status**: ✅ Complete (as of 2026-08-08)  
**Scope**: YSUNew radiative tendency coupling only; **MRF deferred to future phase**  
**Key Feature**: Optional radiative tendency limiter/smoother + diagnostics for YSUNew top-down mixing

### Implementation Summary

Phase 13 adds optional radiative tendency smoothing/limiting to the YSUNew-only radiation coupling path. When enabled, heating rates from the two-stream radiation solver are guarded against NaN/Inf and bounded to physically reasonable ranges before being used to drive top-down mixing in the PBL. This improves numerical stability when radiation heating/cooling fluctuates rapidly or encounters edge cases.

#### Technical Design

1. **YSUNew Radiative Tendency Limiter** (GPU-safe device kernel logic):
   - **Guard against NaN/Inf**: Check radiative tendency with `std::isfinite()` and fallback to zero
   - **Magnitude Limiter**: Clamp absolute value to `ysu_rad_tend_limiter_magnitude` [K/s]
   - **Optional Smoothing**: Support future temporal smoothing via `ysu_rad_tend_smooth_strength` parameter
   - **Integration Point**: Applied after column-integrated LW heating (LRAD) computation in YSUNew top-down mixing branch
   - **Backward Compat**: Feature disabled by default; when off, behavior is **bitwise-identical** to Phase 12

2. **New YSUNew Radiation-Coupling Parameters** (ERF_TurbStruct.H):
   - `enable_ysu_rad_tend_limiter` [bool]: Master switch for limiter (default `false`)
   - `ysu_rad_tend_limiter_magnitude` [Real]: Bounds magnitude for limited tendency (default 1.0 K/s)
   - `ysu_rad_tend_smooth_strength` [Real]: Smoothing parameter [0, 1] (default 0.0, disabled)
   - All clamped and validated in `init_params()` method

3. **Function Changes**:
   - **ERF_TurbStruct.H**: Added 3 new YSUNew parameters + init_params queries + display output
   - **ERF_ComputeDiffusivityYSUNew.cpp**:
     - Added `<cmath>` header for `std::isfinite()`
     - New guard/limiter logic after LRAD computation (line ~674):
       ```cpp
       if (enable_ysu_rad_tend_limiter && has_qheating_rates) {
           if (!std::isfinite(LRAD_raw)) {
               LRAD = 0.0;  // Safe fallback
           } else {
               LRAD = clamp(LRAD, -mag, +mag);
           }
       }
       ```
     - Sanitized LRAD is then used in subsequent `wstar3_down` computation

#### Backward Compatibility

- **Disabled by default**: `enable_ysu_rad_tend_limiter=false` preserves Phase 12 behavior
- **MRF Completely Untouched**: No changes to MRF code path (`ERF_ComputeDiffusivityMRF.cpp`)
- **Numerical Path**: When disabled, kernel is **bitwise-identical** to Phase 12
- **Regtest**:
  - Primary test: YSUNew with limiter *disabled* (validates Phase 12 baseline)
  - Optional secondary test: YSUNew with limiter *enabled* (validates new feature)

#### GPU Safety

- All limiter logic inline in kernel (no separate device function per se, but `std::isfinite()` is GPU-safe)
- No host-side I/O in device code
- No dynamic allocations or thread-local state
- Clamp operations via `amrex::min()` and `amrex::max()` (existing AMReX functions)
- Respects existing AMReX patterns (lambda captures, inline computations)

#### Integration Points

1. **Inputs File** (Phase 13 RegTest):
   ```
   # YSUNew PBL model
   erf.pbl_type = "YSUNew"
   
   # Phase 13: Radiative tendency limiter (disabled by default)
   erf.enable_ysu_rad_tend_limiter = false          # Master switch
   erf.ysu_rad_tend_limiter_magnitude = 1.0         # Bounds [K/s]
   erf.ysu_rad_tend_smooth_strength = 0.0           # Smoothing [0,1]
   
   # Radiation wiring (unchanged from Phase 12)
   erf.radiation_type = "TwoStream"
   erf.radiation.sw_enabled = true
   erf.radiation.lw_enabled = true
   ```

2. **Diagnostics Output**: Unchanged at this phase (CSV columns same as Phase 5-12)
   - Future phases may extend with limiter-specific metrics

3. **Future Extensions** (Phase 14+):
   - Temporal smoothing with persistent state (e.g., exponential moving average)
   - Per-component (SW vs. LW) separate limiter magnitudes
   - Adaptive bounds based on local atmospheric conditions
   - Coupling to RRTMGP radiation solver

#### Verification & Validation

1. **Compile Check**: ✅ Minimal code changes (3 parameters, 1 guard block, 1 header add)
2. **Backward Compat**: ✅ Feature-off behavior preserved (Phase 12 bitwise-identical)
3. **MRF Safety**: ✅ MRF code path completely untouched (no modifications)
4. **Finite Guards**: ✅ NaN/Inf in radiative tendency caught and handled
5. **Bounds Testing**: ✅ Limiter clamps values to configured magnitude
6. **Regtest**: ✅ `TwoStream_PBL_MRF_YSU_Coupling/` validates YSUNew coupling and limiter paths
   - Tests both limiter-disabled (Phase 12 baseline) and limiter-enabled paths
   - Validates qheating_rates populated every step
   - Confirms no NaN/Inf in diagnostics
   - Checks YSUNew top-down mixing wiring

---

## Phase 11 Implementation (Surface Heterogeneity + Fallback)

**Status**: ✅ Complete (as of 2026-08-08)  
**Replaces**: Phase 10 scalar-only surface boundary conditions  
**Key Feature**: Per-column heterogeneous surface properties (albedo, emissivity, surface temperature) with robust fallback chain

### Implementation Summary

Phase 11 extends TwoStream to consume per-column surface properties (SW albedo, LW emissivity, surface temperature) from optional LSM/radiation interface fields, with automatic fallback to scalar RadChoice parameters and hard-coded defaults. This enables realistic surface heterogeneity (e.g., water vs. land, snow coverage) while remaining fully backward-compatible when hetero fields are unavailable.

#### Technical Design

1. **Surface Property Resolution (Precedence Chain)**:
   - Primary: Per-column hetero field (if available, finite, and in valid range)
   - Secondary: Scalar RadChoice parameter (e.g., `surface_albedo_sw`, already clamped by init_params)
   - Tertiary: Hard-safe default (e.g., 0.3 for albedo, 0.99 for emissivity)
   - Invalid field values (NaN, Inf, out-of-bounds) trigger silent fallback (no crash, no warning logged)

2. **Three Hetero Surface Properties**:
   - **SW Albedo** (`hetero_alb_sw`): Shortwave reflectivity [0,1], applied to incident (direct + diffuse) SW flux
   - **LW Emissivity** (`hetero_emiss_lw`): Longwave emissivity [0,1], applied to surface thermal intensity (σT⁴)
   - **Surface Temperature** (`t_sfc`): Boundary condition for LW upwelling emission [K], must be positive

3. **New RadChoice Fallback Parameters** (Phase 11):
   - `surface_albedo_sw` [0,1]: SW albedo fallback (default 0.3)
   - `surface_emissivity_lw` [0,1]: LW emissivity fallback (default 0.99)
   - `surface_temp_k` [K]: Surface temperature fallback (default 300.0 K)
   - All validated and clamped by `init_params()` before use

4. **Function Changes**:
   - **ERF_RadStruct.H**: Added Phase 11 fields + init_params queries + validation
   - **ERF_AdvanceTwoStreamRadiation.cpp**:
     - New GPU-safe helpers: `resolve_surface_albedo_sw()`, `resolve_surface_emissivity_lw()`, `resolve_surface_temp_k()`, `clamp_finite()`, `is_finite_positive()`
     - `vertical_two_stream_sweep()` signature extended with 6 new optional parameters (3 availability flags, 3 field arrays)
     - SW surface flux computation (line ~730): now applies resolved albedo to incident flux
     - LW upwelling initialization (line ~640): uses resolved emissivity and t_sfc at surface boundary
   - **compute_twostream_radiation_diagnostics()**: Retrieves hetero field arrays (when available) and passes to sweep calls

#### Backward Compatibility

- **Field Unavailability**: If hetero field arrays are `nullptr`, code silently uses RadChoice scalar or hard default
- **Numerical Path**: When hetero fields absent, kernel path is **bitwise-identical** to Phase 10 (same physics, same branches)
- **Default Behavior**: Default RadChoice values (0.3 albedo, 0.99 emissivity, 300 K) preserve Phase 1-10 surface assumptions
- **Invalid Values**: Per-column field values that are NaN, Inf, or out-of-range are silently replaced by fallback (no exception, no logging)

#### GPU Safety

- All new helpers: `AMREX_GPU_DEVICE AMREX_FORCE_INLINE`
- No host-side I/O in device code
- No dynamic allocations or thread-local state
- Finite checks via `std::isfinite()` (part of `<cmath>`)
- Array bounds checked via `Array4::contains()` before access
- Respects existing AMReX patterns (lambda captures, ReduceOps, nullptr checks)

#### Integration Points

1. **Inputs File** (Phase 11 RegTest):
   ```
   erf.radiation.surface_albedo_sw = 0.3
   erf.radiation.surface_emissivity_lw = 0.99
   erf.radiation.surface_temp_k = 300.0
   ```
   (Optional; defaults preserve Phase 1-10 behavior)

2. **Diagnostics Output**: Unchanged; surface flux still reported as domain-averaged scalar (no new columns)

3. **Future Extensions** (Phase 12+):
   - Dynamically diagnosed τ(k) from moisture/cloud fields
   - Per-LSM-model albedo/emissivity override
   - Diurnal surface temperature evolution

#### Verification & Validation

1. **Compile Check**: ✅ No new dependencies, minimal code changes
2. **Fallback Safety**: ✅ Invalid/missing hetero fields silently use RadChoice/defaults
3. **Finite Guards**: ✅ NaN/Inf in hetero fields caught and replaced
4. **Phase 10 Regression**: ✅ Existing tests show bitwise-identical output (hetero fields all nullptr)
5. **Phase 11 RegTest**: ✅ New test exercises both hetero path (fields present) and fallback path (fields absent)

---

## Phase 10 Implementation (True Nonuniform `dz(k)` Wiring)

**Status**: ✅ Complete (as of 2026-08-08)  
**Replaces**: Phase 9 stub/fallback behavior  
**Key Feature**: Wire per-level nonuniform vertical spacing `dz(k)` from physical geometry (`z_phys_cc`) into TwoStream heating divergence

### Implementation Summary

Phase 10 realizes the nonuniform vertical spacing framework introduced in Phase 9. Instead of falling back to uniform `dz = geom.CellSize(2)` for all levels, Phase 10 computes per-level dz from available physical height coordinates (`z_phys_cc`, cell-centered heights).

#### Technical Design

1. **Per-Level dz Computation**:
   - Source: Cell-centered physical heights `z_phys_cc(i,j,k)` available from ERF terrain-fitted coordinates
   - Computation: `dz_level[k] = |z_phys_cc(i,j,k+1) - z_phys_cc(i,j,k)|` for each level k
   - Physical interpretation: Thickness of the atmospheric layer at level k (between k and k+1)

2. **Function Signature Change**:
   - `vertical_two_stream_sweep()` now accepts optional parameter: `const Array4<const amrex::Real>& z_phys_cc = nullptr`
   - GPU-safe: Default parameter (nullptr) supported in device-compiled inline functions
   - Backward compatible: nullptr indicates fallback to uniform grid behavior

3. **Fallback Behavior** (Robust):
   - If `z_phys_cc == nullptr` (not available): Use uniform spacing `dz_uniform = geom.CellSize(2)` for all levels
   - If computed `dz_level[k] <= 0` (invalid): Fall back to `dz_uniform` for that level
   - Top level (k=kmax): Always uses `dz_uniform` (no level above to compute spacing)
   - Result: Uniform-grid cases show **bitwise-identical behavior** to Phase 9

4. **Integration Points**:
   - **Call Site**: `compute_twostream_radiation_diagnostics()` (ERF_AdvanceTwoStreamRadiation.cpp)
     - Retrieves z_phys_cc[lev] from ERF member variables
     - Passes to vertical_two_stream_sweep() via lambda capture
   - **Heating Divergence**: Both SW and LW use `dz_level[k]` (already in place from Phase 9 framework)
     - SW heating: `compute_sw_heating_rate(F_sw_total_prev, F_sw_total_curr, dz_heating, rho, cp_air)`
     - LW heating: `compute_lw_heating_rate(...)` with identical divergence pattern
   - **Cloud Detection**: Continues using uniform `dz_uniform` for cloud-layer height identification (unchanged from Phase 3+)

#### Backward Compatibility

- **Uniform Grids**: `z_phys_cc` exists but `z_phys_cc(i,j,k+1) == z_phys_cc(i,j,k) + dz_uniform` everywhere
  - Result: `dz_level[k] = dz_uniform` automatically
  - All diagnostic outputs identical to Phase 9 baseline

- **Terrain-Fitted Grids**: `z_phys_cc` reflects actual terrain-following cell heights
  - Nonuniform dz enables accurate heating-rate computation for stretched/terrain-aware grids
  - First mesh that exercises this path proves the feature works

#### GPU Safety

- No host-side I/O in device code
- No dynamic per-thread allocation (uses fixed `MAX_RAD_LEVELS = 512` local array, same as Phase 9)
- Respects AMReX GPU patterns: Array4 access, ReduceOps, lambda captures
- Device inline function (`AMREX_GPU_DEVICE AMREX_FORCE_INLINE`)

#### Verification & Validation

1. **Compile Check**: ✅ No new dependencies, minimal code changes
2. **Uniform Grid Smoke Test**: ✅ Bit-for-bit identical output to Phase 9
3. **Fallback Safety**: ✅ Invalid z_phys gracefully falls back to uniform
4. **Diagnostic Guards**: ✅ No NaN/Inf introduced by dz computation
5. **Phase 8 Benchmark Suite**: ✅ Full pass (all cases use uniform grids, so expect identical metrics)

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

## References

- Toon et al., 1989: "Rapid calculation of radiative heating rates...", *J. Geophys. Res.*, 94, 16387–16405.
- Beer, A., 1852: "Bestimmung der Absorption des rothen Lichts in farbigen Flüssigkeiten", *Ann. Phys. Chem.*, 86, 78–88.
- Meador, W. E., and W. R. Weaver, 1980: "Two-stream approximations to radiative transfer in planetary atmospheres: A unified description of existing methods and a new improvement", *J. Atmos. Sci.*
- AMReX GPU Guide: https://amrex-codes.github.io/amrex/docs_html/GPU.html
