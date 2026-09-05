# Phase 12 Dynamic Tau Test: Moisture/Cloud-Aware Optical Depth

## Objective

Validate the Phase 12 dynamic optical depth diagnosis feature for TwoStream radiation:
- **SW/LW optical depth** is computed per-level from atmospheric moisture (qv) and cloud liquid (qc)
- **Backward compatibility** maintained: dynamic tau disabled by default, matching Phase 11 output
- **Robust fallback** when dynamic feature is disabled or moisture fields unavailable

## Test Design

### Baseline Scenario (Default)
- **Dynamic tau disabled** (`tau_sw_dynamic_enable=false`, `tau_lw_dynamic_enable=false`)
- **All coefficients zero** (`tau_sw_coeff_qv=0`, `tau_sw_coeff_qc=0`, etc.)
- **Expected behavior**: Identical to Phase 11 (bitwise-compatible)
- **Validates**: No regression when Phase 12 feature is disabled

### Extended Scenario (Optional)
- **Enable dynamic tau** by uncommenting Phase 12 parameters in inputs
- **Set nonzero coefficients** (e.g., `tau_sw_coeff_qv=10.0`, `tau_sw_coeff_qc=100.0`)
- **Use moist sounding** (qv ~0.005-0.008 kg/kg in BL, decreasing with height)
- **Expected behavior**: Optical depth varies per-level based on moisture content
- **Validates**: Dynamic path exercised and produces reasonable heating rates

## Inputs

- **inputs**: Control file with Phase 12 parameters (disabled by default)
- **input_sounding_phase12_moist**: Atmospheric sounding with moisture profile
  - Surface: qv ~0.008 kg/kg, qc ~0.0001 kg/kg (dry cloud edge)
  - Upper troposphere: qv decreases to ~0.005 kg/kg

## Validation Criteria

1. **Diagnostics file** (`radiation_diag_phase12.dat`) must exist and parse correctly
2. **Flux values** (SW_TOA, SW_surface, LW_net) must be finite and physically reasonable
3. **Heating rates** must be nontrivial and not contain NaN/Inf
4. **Static tau mode** (default): output matches Phase 11 regression baseline
5. **Dynamic tau mode** (when enabled): optical depth varies with height due to moisture gradient

## Running the Test

### Baseline (Static Tau Mode)
```bash
erf inputs
# Check: radiation_diag_phase12.dat has finite values
# Check: qsrc_sw, qsrc_lw fields in plot files show reasonable heating
```

### Extended (Dynamic Tau Mode)
Edit `inputs` and uncomment:
```
erf.radiation.tau_sw_dynamic_enable = true
erf.radiation.tau_lw_dynamic_enable = true
erf.radiation.tau_sw_coeff_qv = 10.0
erf.radiation.tau_sw_coeff_qc = 100.0
erf.radiation.tau_lw_coeff_qv = 20.0
erf.radiation.tau_lw_coeff_qc = 200.0
```

Then run:
```bash
erf inputs
# Check: radiation_diag_phase12.dat shows tau varying with moisture
# Check: heating rates change vs. static tau baseline
```

## Key Features Validated

- **GPU-safe device code**: Dynamic tau functions use AMReX-compatible inline device kernels
- **Numerical guards**: NaN/Inf handling, clamping to [0, 100] range
- **Fallback safety**: When fields unavailable, silently uses static tau
- **Parameter validation**: Coefficients clamped to nonnegative
- **Backward compatibility**: Zero coefficients = no-op, identical to Phase 11

## References

- RAD_DEVELOPMENT.md Phase 12 section
- Source/Radiation/ERF_AdvanceTwoStreamRadiation.cpp (dynamic tau functions)
- Source/DataStructs/ERF_RadStruct.H (RadChoice parameters)
