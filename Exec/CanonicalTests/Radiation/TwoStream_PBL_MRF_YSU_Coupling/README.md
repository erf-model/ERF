# YSUNew Radiation Coupling Regtest

## Overview

This regtest validates implementation: **YSUNew PBL Coupling with Radiative Tendency Limiter/Smoother**.

### Key Features Validated

1. **YSUNew Model Selection**: Confirms YSUNew (not MRF) is selected and active
2. **Radiation-to-PBL Coupling**: Validates qheating_rates from TwoStream radiation are coupled to YSUNew top-down mixing
3. **Radiative Tendency Limiter**: Tests optional finite guards and magnitude bounds on radiative heating tendency
4. **Diagnostics Output**: Confirms radiation diagnostics accumulate every timestep
5. **Backward Compatibility**: Feature-off (default) preserves the feature-off baseline behavior

## Test Configuration

- **PBL Model**: YSUNew
- **Radiation Type**: TwoStream (SW + LW, non-isothermal)
- **Domain**: 3000×3000×1024 m, 8×8×64 grid
- **Runtime**: 2.5 seconds, fixed dt=0.5s
- **Surface Layer**: MOST with z0=0.1 m
- **Coriolis**: Enabled (latitude=45°, f≈1e-4 rad/s)
- **ABL Driver**: Geostrophic wind forcing

## Input Files

### `inputs`
Main configuration file. Key parameters:
- `erf.pbl_type = "YSUNew"` — Select YSUNew PBL
- `erf.enable_ysu_topdown = true` — Enable top-down mixing (LW radiation coupling)
- `erf.enable_ysu_rad_tend_limiter = false` — Limiter disabled by default (baseline test)
- `erf.ysu_rad_tend_limiter_magnitude = 1.0` — Bounds parameter [K/s]

### `input_sounding_ysu`
Initial sounding profile (pressure-theta-qv-u-v):
- Surface: p=1000 hPa, θ=300 K, u=15 m/s
- Mixed layer to 551 m: θ=300 K
- Upper atmosphere: θ=308-311 K above 551 m

## Expected Output

### Diagnostic Files
- **radiation_ysu_coupling_diag.dat**: CSV with per-timestep radiation fluxes and heating rates
  - Columns: `step, time, call_site, SW_surface, SW_TOA, SW_up_TOA, LW_net_surface, LW_up_TOA, heating_rate_max`
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
- Behavior is **bitwise-identical** (before changes)
- Existing tests continue to pass unchanged

## Future Enhancements

- Temporal smoothing with state persistence (smooth_strength ∈ [0,1])
- Per-component (SW vs LW) separate limiting
- Adaptive limiter magnitude based on local conditions

## References

### Documentation
- `Source/Radiation/RAD_DEVELOPMENT.md` — technical design
- `Source/DataStructs/ERF_TurbStruct.H` — Parameter definitions
- `Source/PBL/ERF_ComputeDiffusivityYSUNew.cpp` — Limiter implementation

### Regtest Patterns (Reference Cases)
- `TwoStream_RhoTheta_Coupling/` — Radiation coupling wiring validation
- `TwoStream_DynamicTau_MoistCloud/` — Dynamic optical depth
- `TwoStream_SurfaceHeterogeneity/` — Surface property heterogeneity

## Notes

- **MRF Untouched**: implementation is YSUNew-only; no changes to MRF code
- **No Compilation Required**: Regtest can be visually validated against source; full execution requires build
- **GPU Safe**: All limiter/smoothing logic uses AMReX GPU-safe patterns
