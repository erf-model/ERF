# Phase 18: Simplified SEB — Diagnostic Mode

## Objective

Validate the Phase 18 SEB diagnostic residual computation feature:
- **SEB residual diagnosed** from net radiation and turbulent/ground heat fluxes
- **Diagnostic-only** — no prognostic surface temperature or flux update (Phase 19 work)
- **Backward compatibility** — when disabled (default), output is bitwise-identical to Phase 17
- **GPU-safe implementation** — residual computed via device-side reduction kernels

## Test Design

### Baseline Scenario (Disabled)
- **SEB diagnostic disabled** (`seb_diagnostic_enable=false`, default)
- **Phase 17 infrastructure still active** (`seb_enable=true`)
- **Expected behavior**: Identical to Phase 17 baseline output
- **Validates**: Full backward compatibility; no new computation when feature is off

### Feature-On Scenario
- **SEB diagnostic enabled** (`seb_diagnostic_enable=true`)
- **SEB infrastructure active** (`seb_enable=true`)
- **No LSM active** — all fluxes from scalar fallback defaults (deterministic)
- **Expected behavior**: SEB residual computed and reported in diagnostics
- **Validates**: Residual computation correct and diagnostic output present

## Input Files

### `input_sounding` (Atmospheric Profile)
Physically reasonable mid-latitude sounding with moisture profile:
- Surface (0 m): T=300K, qv=0.008 kg/kg (~8 g/kg, typical mid-latitude)
- Upper levels: qv decays to 0.004 kg/kg at 1550 m (realistic moisture gradient)

See Phase 17 README for full sounding documentation.

### Input Configuration Files

#### `inputs_seb_diagnostic_disabled`
- **seb_diagnostic_enable = false** (default)
- Tests baseline case for backward compatibility
- Should produce bitwise-identical output to Phase 17

#### `inputs_seb_diagnostic_enabled`
- **seb_diagnostic_enable = true** (feature on)
- Uses Phase 17 scalar fallback defaults:
  - `seb_sw_flux_default = 50.0` W/m^2
  - `seb_lw_flux_default = -25.0` W/m^2
  - `seb_hfx_default = 10.0` W/m^2
  - `seb_lh_default = 20.0` W/m^2
  - `seb_grdflux_default = 5.0` W/m^2

## Validation Criteria

### Baseline Test
1. **No new CSV columns** — CSV output must have exactly 8 columns (Phase 17 format)
2. **No SEB diagnostics** — CSV files contain no `SEB_residual_*` values
3. **Bitwise-identical output** — Radiation fluxes/heating rates match Phase 17 exactly
4. **Finite values** — All flux diagnostics finite and physically reasonable

### Feature-On Test
1. **New CSV columns present** — CSV output has 10 columns including SEB residuals
   - `SEB_residual_mean`: mean residual across all surface columns [W/m^2]
   - `SEB_residual_max`: maximum |residual| across surface columns [W/m^2]
2. **Residual computation correct** — Expected value calculation:
   ```
   R_net = SW_net + LW_net = 50.0 + (-25.0) = 25.0 W/m^2
   SEB_residual = R_net - H - LE - G = 25.0 - 10.0 - 20.0 - 5.0 = -10.0 W/m^2
   ```
   Residual must be approximately **-10.0 W/m^2** (within ±0.1 tolerance)
3. **Finite diagnostics** — All residual values finite (no NaN/Inf)
4. **No impact on radiation** — SW/LW/heating diagnostics identical to baseline
   (residual is diagnostic-only; does not affect physics)

## Surface Energy Balance Equation

At each surface column (i,j):

```
R_net(i,j) = SW_net(i,j) + LW_net(i,j)

SEB_residual(i,j) = R_net(i,j) - H(i,j) - LE(i,j) - G(i,j)
```

Where:
- **SW_net**: Net shortwave flux at surface (positive down)
- **LW_net**: Net longwave flux at surface (positive up)
- **H**: Sensible heat flux
- **LE**: Latent heat flux
- **G**: Ground heat flux

A perfectly closed budget gives `SEB_residual ≈ 0`. In this test with scalar fallback defaults, the residual reflects the parameterized constant values only — not physically balanced, but deterministic and finite (safe diagnostic).

## Running the Tests

### Baseline Mode
```bash
erf inputs_seb_diagnostic_disabled
python check_seb_diagnostic.py
# Verify: 8 columns, no SEB residual output, bitwise-identical to Phase 17
```

### Feature-On Mode
```bash
erf inputs_seb_diagnostic_enabled
python check_seb_diagnostic.py
# Verify: 10 columns, SEB residual ~-10.0 W/m^2, finite values
```

## Phase 18 Implementation Summary

### New Files
- `Source/Radiation/ERF_SimplifiedSEB.H` — GPU-safe residual diagnostic kernel
- `Exec/CanonicalTests/Radiation/TwoStream_SEB_Diagnostic/` — RegTest directory

### Modified Files
- `Source/DataStructs/ERF_RadStruct.H` — Added `seb_diagnostic_enable` parameter
- `Source/Radiation/ERF_RadiationDiagnostics.H/.cpp` — Extended CSV output with SEB columns
- `Source/Radiation/ERF_AdvanceTwoStreamRadiation.cpp` — Integrated residual computation
- `Source/Radiation/RAD_DEVELOPMENT.md` — Phase 18 section and roadmap update

### Key Design Decisions
1. **Diagnostic-only**: No feedback to T_s, heating rates, or any prognostic fields
2. **GPU-safe**: All computation via `AMREX_GPU_DEVICE AMREX_FORCE_INLINE` kernels
3. **Backward compatible**: When disabled (default), zero overhead and bitwise-identical output
4. **Auto-enable SEB**: If `seb_diagnostic_enable=true` but `seb_enable=false`, auto-enable SEB internally
5. **Safe no-op**: If any input flux is NaN/Inf, return 0.0 residual (safe fallback)

## References

- `Source/Radiation/RAD_DEVELOPMENT.md` — Phase 18 Implementation section
- `Source/DataStructs/ERF_RadStruct.H` — RadChoice parameters documentation
- Oke, T. R., 1987: Boundary Layer Climates (2nd ed.), Routledge. [SEB theory reference]
