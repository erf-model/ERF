# Phase 11 Two-Stream Radiation RegTest: Surface Heterogeneity + Fallback

## Overview

This regression test validates the **Phase 11** implementation of per-column heterogeneous surface properties in the TwoStream radiation solver, with a robust fallback chain for missing or invalid data.

## Test Purpose

Phase 11 extends TwoStream to:

1. **Consume Per-Column Surface Properties**: Accept albedo, emissivity, and surface temperature from optional LSM/radiation interface fields
2. **Implement Robust Fallback Chain**: 
   - Primary: Use hetero field value if finite and in valid range
   - Secondary: Fall back to scalar RadChoice parameter (from inputs file)
   - Tertiary: Fall back to hard-coded default value
3. **Maintain Backward Compatibility**: When hetero fields unavailable, produce bitwise-identical results to Phase 10

## Test Configuration

### Inputs File: `inputs`

- **Domain**: 3000m × 3000m × 1024m (small for fast CI testing)
- **Grid**: 8 × 8 × 64 cells
- **Radiation**: TwoStream (SW and LW enabled)
- **Surface Properties**:
  - `surface_albedo_sw = 0.3` (typical land/water)
  - `surface_emissivity_lw = 0.99` (blackbody-like)
  - `surface_temp_k = 300.0` (typical surface)
- **Duration**: 2.5 seconds (5 timesteps at 0.5s each)

### Sounding File: `input_sounding_hetero`

Simple neutral sounding with constant potential temperature (300 K).

### Checker Script: `check_hetero_accuracy.py`

Validates:
- Diagnostics file created and parseable
- Finite values (no NaN/Inf)
- Nonzero heating rates
- Sensible surface flux values
- Fallback path exercised

## What Is Being Tested?

### Scenario A: Fallback Mode (Default)

In this test scenario:
- **Hetero Fields**: All unavailable (nullptr)
- **Surface Properties Used**: RadChoice scalar parameters
  - `surface_albedo_sw = 0.3`
  - `surface_emissivity_lw = 0.99`
  - `surface_temp_k = 300.0`
- **Expected Behavior**: 
  - SW flux = (incident flux) × (1 - 0.3) [i.e., 70% absorbed]
  - LW flux = 0.99 × σ × (300)^4 ≈ 450 W/m²
  - Heating rates computed from flux divergence

### Scenario B: Heterogeneous Mode (Future/Extension)

To test the hetero field path:
1. Modify inputs to provide spatially varying surface properties
2. Create hetero_alb_sw, hetero_emiss_lw, t_sfc MultiFabs
3. Run with hetero fields populated
4. Verify that hetero values override fallback scalars

This scenario is not included in this basic test but can be added by:
- Extending the inputs file with hetero field generation
- Creating a companion test case with spatial heterogeneity

## Expected Outputs

### Radiation Diagnostics File: `radiation_diag_phase11.dat`

CSV format with columns:
```
step, time, call_site, SW_surface, SW_TOA, F_up_surface, F_down_toa, heating_rate_max
```

### Typical Values (Fallback Mode)

- **SW_TOA**: ~962 W/m² (= 1361 × cos(45°), constant in time)
- **SW_surface**: ~300-400 W/m² (absorbed = (1-0.3) × incident)
- **F_up_surface**: ~450 W/m² (= 0.99 × σ × 300^4 from surface temp)
- **heating_rate_max**: ~1-10 K/s (typical clear-sky radiative heating)

### Plotfiles

Generated at every 5 timesteps:
- `plt_hetero/Plt00000/` etc.
- Contains velocity, temperature, and radiative heating fields

### Data Logs

- `hetero_hist.dat` - Time history of domain-averaged quantities
- `hetero_profiles.dat` - Vertical profiles at each output interval

## Running This Test

### Quick Run (CPU)

```bash
cd /path/to/ERF/build
./erf /path/to/ERF/Exec/CanonicalTests/Radiation/TwoStream_SurfaceHeterogeneity/inputs
```

Expected runtime: 30 seconds to 2 minutes

### With Custom Parameters

```bash
./erf \
  /path/to/ERF/Exec/CanonicalTests/Radiation/TwoStream_SurfaceHeterogeneity/inputs \
  max_level=0 \
  amr.n_cell="16 16 64" \
  stop_time=5.0 \
  erf.radiation.surface_albedo_sw=0.5
```

### Validate Results

```bash
python3 check_hetero_accuracy.py
```

This runs the validation script in the current directory (where the test was run).

## Key Phase 11 Features Exercised

1. **Helper Functions**:
   - `resolve_surface_albedo_sw()`: Resolves per-column albedo
   - `resolve_surface_emissivity_lw()`: Resolves per-column emissivity
   - `resolve_surface_temp_k()`: Resolves per-column surface temperature
   - `clamp_finite()`: Safely clamps invalid values
   - `is_finite_positive()`: Validates temperature values

2. **Function Signature Updates**:
   - `vertical_two_stream_sweep()` accepts 6 new optional parameters
   - Backwards compatible (default nullptr for all fields)

3. **Physics Integration**:
   - SW flux *= (1 - albedo) [surface absorption]
   - LW upwelling = emissivity × σ × T^4 [surface emission]
   - Both fallback to RadChoice or hard defaults if hetero field unavailable

4. **Diagnostics**:
   - CSV file accumulates data over multiple timesteps
   - All values finite and physically sensible
   - No crashes or assertion failures

## Backward Compatibility Validation

This test implicitly validates Phase 10 compatibility:
- When hetero fields are nullptr (as they are here), the code path is **identical** to Phase 10
- Domain-averaged surface fluxes should match Phase 10 output exactly
- Heating rates should show identical spatial and temporal patterns

To explicitly verify:
1. Run this test → get diagnostics
2. Run equivalent Phase 10 test (e.g., TwoStream_NonuniformDZ) → get diagnostics
3. Compare CSV files → differences should be < 1e-12 (rounding only)

## Future Extensions

### Heterogeneous Albedo Test

Add LSM-style heterogeneous surface fields:
```
inputs_hetero_spatial:
  - ocean region: albedo = 0.06
  - desert region: albedo = 0.35
  - snow region: albedo = 0.80
```

Expected result: SW surface flux varies by ~5-10x across domain

### Heterogeneous Temperature Test

Add time-varying surface temperature (e.g., diurnal cycle):
```
inputs_hetero_temporal:
  - Morning: t_sfc = 280 K → F_up ≈ 390 W/m²
  - Afternoon: t_sfc = 320 K → F_up ≈ 520 W/m²
```

Expected result: LW upwelling flux varies with local time

## Troubleshooting

### Issue: Test fails to compile

**Solution**: Check that all Phase 11 code changes were applied:
- ERF_RadStruct.H: Three new fields + init_params() queries
- ERF_AdvanceTwoStreamRadiation.cpp: Five helper functions + modified vertical_two_stream_sweep()

### Issue: Checker script reports failures

**Solution**: Inspect radiation_diag_phase11.dat:
- Check that file exists and contains data
- Look for non-finite values (NaN, Inf)
- Verify flux ranges are sensible (0-2000 W/m²)
- Check heating_rate_max > 1e-10 (not zero)

### Issue: Test produces NaN in output

**Solution**:
1. Verify surface temperature is positive (default 300 K)
2. Check that albedo/emissivity are in [0,1] (defaults: 0.3, 0.99)
3. Look for invalid values in sounding file

## References

- **RAD_DEVELOPMENT.md**: Detailed Phase 11 implementation notes
- **ERF_AdvanceTwoStreamRadiation.cpp**: Source code (lines with Phase 11 comments)
- **MANUAL_VERIFICATION.md**: Step-by-step verification procedures

---

**RegTest Version**: 1.0  
**Phase**: 11 (Surface Heterogeneity + Fallback)  
**Created**: 2026-08-08  
**Status**: Production
