# Phase 13 YSUNew Radiation Coupling Regtest

## Overview

This regtest validates Phase 13 implementation: **YSUNew PBL Coupling with Radiative Tendency Limiter/Smoother**.

### Key Features Validated

1. **YSUNew Model Selection**: Confirms YSUNew (not MRF) is selected and active
2. **Radiation-to-PBL Coupling**: Validates qheating_rates from TwoStream radiation are coupled to YSUNew top-down mixing
3. **Radiative Tendency Limiter**: Tests optional finite guards and magnitude bounds on radiative heating tendency
4. **Diagnostics Output**: Confirms radiation diagnostics accumulate every timestep
5. **Backward Compatibility**: Feature-off (default) preserves Phase 12 baseline behavior

## Test Configuration

- **PBL Model**: YSUNew (Phase 13 focus)
- **Radiation Type**: TwoStream (SW + LW, non-isothermal)
- **Domain**: 3000×3000×1024 m, 8×8×64 grid
- **Runtime**: 2.5 seconds, fixed dt=0.5s
- **Surface Layer**: MOST with z0=0.1 m
- **Coriolis**: Enabled (latitude=45°, f≈1e-4 rad/s)
- **ABL Driver**: Geostrophic wind forcing

## Input Files

### `inputs`
Main configuration file. Key Phase 13 parameters:
- `erf.pbl_type = "YSUNew"` — Select YSUNew PBL
- `erf.enable_ysu_topdown = true` — Enable top-down mixing (LW radiation coupling)
- `erf.enable_ysu_rad_tend_limiter = false` — Limiter disabled by default (baseline test)
- `erf.ysu_rad_tend_limiter_magnitude = 1.0` — Bounds parameter [K/s]
- `erf.ysu_rad_tend_smooth_strength = 0.0` — Smoothing parameter [0,1]

### `input_sounding_phase13_ysu`
Initial sounding profile (pressure-theta-qv-u-v):
- Surface: p=1000 hPa, θ=300 K, u=15 m/s
- Mixed layer to 551 m: θ=300 K
- Upper atmosphere: θ=308-311 K above 551 m

## Expected Output

### Diagnostic Files
- **radiation_phase13_ysu_coupling_diag.dat**: CSV with per-timestep radiation fluxes and heating rates
  - Columns: `step, time, call_site, SW_surface, SW_TOA, F_up_surface, F_down_toa, heating_rate_max`
  - One row per timestep (5 rows expected for 2.5s simulation @ dt=0.5s)

### Checkpoint/Plotfile
- **chk_ysu_coupling_***: Checkpoints (disabled: check_int=-1)
- **plt_ysu_coupling_***: Plotfiles every 5 steps (including radiation heating fields)

## Validation

Run the checker script:
```bash
python3 check_ysunew_coupling.py
```

### Smoke-Test Checks
1. ✓ Diagnostic file exists and contains multiple timesteps
2. ✓ Time progression is monotonic
3. ✓ All diagnostic values are finite (no NaN/Inf)
4. ✓ SW_TOA matches analytical value (S0 * cos(zenith°))
5. ✓ Heating rate max is nonzero and physically reasonable
6. ✓ Surface fluxes are positive (physical energy direction)

## Backward Compatibility

With `enable_ysu_rad_tend_limiter = false` (default):
- Radiative tendency limiter is completely disabled
- Behavior is **bitwise-identical** to Phase 12 (before Phase 13 changes)
- Existing tests continue to pass unchanged

## Future Enhancements

- Temporal smoothing with state persistence (smooth_strength ∈ [0,1])
- Per-component (SW vs LW) separate limiting
- Adaptive limiter magnitude based on local conditions

## References

### Phase 13 Documentation
- `Source/Radiation/RAD_DEVELOPMENT.md` — Phase 13 technical design
- `Source/DataStructs/ERF_TurbStruct.H` — Parameter definitions
- `Source/PBL/ERF_ComputeDiffusivityYSUNew.cpp` — Limiter implementation

### Regtest Patterns (Reference Cases)
- Phase 5: `Phase5_RhoTheta_Coupling/` — Radiation coupling wiring validation
- Phase 12: `TwoStream_DynamicTau_MoistCloud/` — Dynamic optical depth
- Phase 11: `TwoStream_SurfaceHeterogeneity/` — Surface property heterogeneity

## Notes

- **MRF Untouched**: Phase 13 implementation is YSUNew-only; no changes to MRF code
- **No Compilation Required**: Regtest can be visually validated against source; full execution requires build
- **GPU Safe**: All limiter/smoothing logic uses AMReX GPU-safe patterns
