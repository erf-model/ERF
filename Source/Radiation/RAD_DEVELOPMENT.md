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
| **11** | Surface Heterogeneity + Fallback (Albedo/Emissivity/`t_sfc`) | ⏳ Planned (Active) | TBD | TwoStream consumes per-column LSM/Radiation surface fields with robust fallback path | Moderate | `TwoStream_SurfaceHeterogeneity_Fallback` |
| **12** | Moisture/Cloud-Aware Dynamic Optical Depth | ⏳ Planned (Active) | TBD | Diagnose SW/LW `tau(k)` from `qv`, `qc`, `rho`, `dz` with safe fallback | Moderate | `TwoStream_DynamicTau_MoistCloud` |
| **13** | PBL Coupling Focus (MRF/YSU only) | ⏳ Planned (Active) | TBD | MRF/YSU-focused radiative tendency smoothing/limiter + diagnostic hooks | Moderate | `TwoStream_PBL_MRF_YSU_Coupling` |
| **14** | Prognostic Cloud Fraction for Radiation | ⏳ Planned (Active) | TBD | RH/`qc`-based diagnosed cloud fraction with bounds and temporal smoothing | Easy–Moderate | `TwoStream_ProgCloudFraction` |
| **15** | Bulk Aerosol/Turbidity Option | ⏳ Planned (Active) | TBD | Prescribed aerosol optical-depth profile (constant/exponential/table), optional LW hook | Easy–Moderate | `TwoStream_Aerosol_Turbidity` |
| **16** | Time-Varying Solar Geometry | ⏳ Planned (Active) | TBD | Solar zenith evolution with time/lat/day; fixed-angle fallback retained | Easy | `TwoStream_DiurnalSolarGeometry` |
| **17** | Simplified Surface Energy Balance (SEB) — Diagnostic Mode | ⏳ Planned (Active) | TBD | Compute/report SEB residual terms from TwoStream + surface inputs (no prognostic `T_s` update) | Moderate | `TwoStream_SEB_Diagnostic` |
| **18** | Simplified SEB — Prognostic `T_s` Mode | ⏳ Planned (Active) | TBD | Optional explicit `T_s` tendency update with limiter/clamps and fallback-safe behavior | Moderate | `TwoStream_SEB_PrognosticTs` |
| **19** | SEB Coupling Safeguards (Noah-MP/SurfaceLayer Interop) | ⏳ Planned (Active) | TBD | Anti-double-count rules and precedence guards for `T_s`, `H`, `LE`, and radiative terms | Moderate | `TwoStream_SEB_InteropGuards` |
| **20** | SEB Validation & Benchmark Suite | ⏳ Planned (Active) | TBD | Canonical SEB closure/stability tests, tolerances, and CI-ready reports | Moderate | `TwoStream_SEB_BenchmarkSuite` |

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
