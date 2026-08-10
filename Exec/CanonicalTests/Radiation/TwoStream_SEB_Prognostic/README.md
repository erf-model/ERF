# Phase 19b: Simplified SEB — Prognostic Surface Temperature and Moisture

## Objective

Validate the Phase 19b SEB prognostic evolution feature:
- **Prognostic T_s and q_s evolution** from SEB residual using force-restore formulation
- **Time integration** via explicit Euler with configurable timescales and bounds
- **Noah-MP gating** — update skipped when Noah-MP actively drives LSM fields
- **Backward compatibility** — when disabled (default), output is bitwise-identical to Phase 19a
- **GPU-safe implementation** — time integration via device-side update kernels

## Test Design

### Baseline Scenario (Disabled)
- **SEB prognostic disabled** (`seb_prognostic_enable=false`, default)
- **Phase 17/18 infrastructure active** (`seb_enable=true`, `seb_diagnostic_enable=true`)
- **Expected behavior**: Identical to Phase 19a baseline output
- **Validates**: Full backward compatibility; no new computation when feature is off

### Feature-On Scenario
- **SEB prognostic enabled** (`seb_prognostic_enable=true`)
- **Phase 17/18 infrastructure active** (auto-enabled if needed)
- **No LSM active** — all fields from scalar fallback defaults (deterministic)
- **Expected behavior**: T_s and q_s evolve over time steps according to force-restore equations
- **Validates**: Time integration correct, clamping bounds respected, diagnostics present

## Input Files

### `input_sounding` (Atmospheric Profile)
Physically reasonable mid-latitude sounding with moisture profile (same as Phase 18):
- Surface (0 m): T=300K, qv=0.008 kg/kg (~8 g/kg, typical mid-latitude)
- Upper levels: qv decays to 0.004 kg/kg at 1550 m (realistic moisture gradient)

See Phase 17 README for full sounding documentation.

### Input Configuration Files

#### `inputs_seb_prognostic_disabled`
- **seb_prognostic_enable = false** (default)
- Tests baseline case for backward compatibility
- Should produce bitwise-identical output to Phase 19a

#### `inputs_seb_prognostic_enabled`
- **seb_prognostic_enable = true** (feature on)
- Uses Phase 17 scalar fallback defaults:
  - `seb_sw_flux_default = 50.0` W/m^2
  - `seb_lw_flux_default = -25.0` W/m^2
  - `seb_hfx_default = 10.0` W/m^2
  - `seb_lh_default = 20.0` W/m^2
  - `seb_grdflux_default = 5.0` W/m^2
  - `seb_q_sfc_default = 0.01` kg/kg
  - `seb_t_deep_default = 295.0` K
  - `seb_q_deep_default = 0.20` kg/kg
- Uses Phase 19b prognostic parameters:
  - `seb_surface_heat_capacity = 2.0e4` J/(m^2*K)
  - `seb_restore_timescale_s = 86400.0` s (1 day, weak damping)
  - `seb_moisture_layer_depth_m = 0.1` m
  - `seb_moisture_restore_timescale_s = 86400.0` s
  - `seb_prognostic_t_min_k = 200.0` K (min clamp)
  - `seb_prognostic_t_max_k = 340.0` K (max clamp)
  - `seb_prognostic_q_min = 0.0` kg/kg (min clamp)
  - `seb_prognostic_q_max = 1.0` kg/kg (max clamp)

## Validation Criteria

### Baseline Test
1. **No new CSV columns** — CSV output must have exactly 12 columns (Phase 19a format)
2. **No prognostic diagnostics** — CSV files contain only NaN for `T_s_mean`, `T_s_max`, `q_s_mean`, `q_s_max`
3. **Bitwise-identical output** — All diagnostics match Phase 19a baseline
4. **Finite values** — All reported diagnostics are finite

### Feature-On Test
1. **New CSV columns present** — CSV output has 16 columns including T_s and q_s prognostic values
   - `T_s_mean`: mean surface temperature [K]
   - `T_s_max`: maximum surface temperature [K]
   - `q_s_mean`: mean surface moisture [kg/kg]
   - `q_s_max`: maximum surface moisture [kg/kg]
2. **Temperature evolution**:
   - Expected SEB residual: `R_net - H - LE - G = (50 - 25) - 10 - 20 - 5 = -10 W/m^2`
   - Tendency (first order): `dT_s/dt ≈ -10 / 2.0e4 = -0.0005 K/s = -1.8 K/hr`
   - Temperature should decrease from initial T_s in direction of negative residual
   - Must remain within bounds `[200.0, 340.0]` K throughout run
3. **Moisture evolution**:
   - Tendencies computed from latent heat flux and restoring term
   - Must remain within bounds `[0.0, 1.0]` kg/kg throughout run
4. **Restoring term validation**:
   - With very small timescale (e.g., 100 s), T_s should approach T_deep asymptotically
   - Restoring term sign/magnitude correct (approaching deep value from above/below)
5. **Finite diagnostics** — All T_s and q_s values finite (no NaN/Inf)
6. **No impact on radiation** — SW/LW/heating diagnostics identical to baseline
   (prognostic update occurs after radiation calculation; no feedback)

## Surface Energy Balance Prognostic Equations

### Temperature Evolution (Force-Restore)

```
C_s * dT_s/dt = R_net - H - LE - G - C_s * (2*pi/tau) * (T_s - T_deep)

dT_s/dt = (R_net - H - LE - G) / C_s - (2*pi/tau) * (T_s - T_deep)
        = SEB_residual / C_s - (2*pi/tau) * (T_s - T_deep)
```

Where:
- **SEB_residual** = R_net - H - LE - G
- **C_s** = effective surface heat capacity [J/(m^2*K)]
- **tau** = force-restore timescale [s]
- **T_deep** = deep soil temperature [K]

Euler update:
```
T_s^(n+1) = T_s^n + dt * dT_s/dt^n
           clamped to [T_min, T_max]
```

### Moisture Evolution (Force-Restore, Bucket-Style)

```
dq_s/dt = -(LE / (L_v * rho_w * d_s)) - (1/tau_q) * (q_s - q_deep)
```

Where:
- **LE** = latent heat flux [W/m^2]
- **L_v** = 2.5e6 J/kg (latent heat of vaporization, hardcoded)
- **rho_w** = 1000.0 kg/m^3 (water density, hardcoded)
- **d_s** = effective surface moisture layer depth [m]
- **tau_q** = moisture force-restore timescale [s]
- **q_deep** = deep soil moisture [kg/kg]

Euler update:
```
q_s^(n+1) = q_s^n + dt * dq_s/dt^n
           clamped to [q_min, q_max]
```

## Running the Tests

### Baseline Mode
```bash
erf inputs_seb_prognostic_disabled
python check_seb_prognostic.py baseline
# Verify: 12 columns (Phase 19a format), T_s/q_s columns are NaN, bitwise-identical to Phase 19a
```

### Feature-On Mode
```bash
erf inputs_seb_prognostic_enabled
python check_seb_prognostic.py feature_on
# Verify: 16 columns, T_s/q_s columns finite, temperatures evolve in correct direction,
#         all values within configured bounds, no NaN/Inf
```

## Phase 19b Implementation Summary

### New Files
- `Source/Radiation/ERF_SimplifiedSEB.H` — GPU-safe prognostic tendency kernels
- `Exec/CanonicalTests/Radiation/TwoStream_SEB_Prognostic/` — RegTest directory

### Modified Files
- `Source/DataStructs/ERF_RadStruct.H` — Added 9 new prognostic parameters
- `Source/Radiation/ERF_RadiationDiagnostics.H/.cpp` — Extended CSV output with T_s/q_s columns
- `Source/Radiation/ERF_AdvanceTwoStreamRadiation.cpp` — Integrated prognostic update with Noah-MP gating
- `Source/Radiation/RAD_DEVELOPMENT.md` — Phase 19b section and roadmap update

### Key Design Decisions
1. **Prognostic-only update**: T_s and q_s evolved in place; no feedback to radiation or atmosphere
2. **GPU-safe**: All time integration via `AMREX_GPU_DEVICE AMREX_FORCE_INLINE` kernels
3. **Noah-MP gating**: When Noah-MP drives LSM at a level, Phase 19b update skipped (LSM takes precedence)
4. **Auto-enable prerequisites**: If `seb_prognostic_enable=true`, auto-enable `seb_enable` and `seb_diagnostic_enable`
5. **Safe no-op on non-finite**: If any input NaN/Inf or parameters invalid, return 0.0 tendency (no update)
6. **Clamping**: All updated values clamped to configured bounds to prevent instability

## References

- `Source/Radiation/RAD_DEVELOPMENT.md` — Phase 19b Implementation section
- `Source/DataStructs/ERF_RadStruct.H` — RadChoice parameters documentation
- Oke, T. R., 1987: Boundary Layer Climates (2nd ed.), Routledge. [SEB theory reference]
