# Clear-Sky Shortwave Analytical Test

## Objective

Validate shortwave (solar) radiation transport using the Beer-Lambert direct-beam attenuation model against an exact analytical solution.

## Test Purpose

This test confirms that the two-stream shortwave solver correctly computes:
- Direct-beam solar flux at any height using exponential attenuation
- Zero-diffuse contributions (clear-sky approximation)
- Proper solar geometry handling (zenith angle)
- Consistency with analytical Beer-Lambert prediction

## Test Design

### Configuration

- **Domain**: 1000 m × 1000 m horizontal, 10 km vertical (20 layers)
- **Time**: Single timestep, 0.1 s duration
- **Solar Constant**: S₀ = 1361 W/m²
- **Solar Zenith Angle**: 60°
- **Optical Depth per Layer**: τ = 0.05 (uniform, clear-sky)

### Key Physics

The direct-beam flux at height z follows Beer-Lambert attenuation:

```
F_dir(z) = S₀ · cos(zenith) · exp(-τ_cumulative(z) / cos(zenith))
```

where τ_cumulative = τ_per_layer × (number of layers above z).

### Analytical Solution

- **TOA Flux**: S₀ · cos(60°) ≈ 680.5 W/m²
- **Incident Surface Flux** (after 64 layers, τ = 0.003125 per layer): 680.5 × exp(-0.2 / cos(60°)) ≈ 456.2 W/m²
- **Absorbed Surface Flux** (the `SW_surface` diagnostic, albedo 0.3): 456.2 × (1 − 0.3) ≈ 319.3 W/m²
- **Heating Rate**: 0 (transparent medium, no absorption to temperature)

## Files

- `inputs` — Main control file with radiation parameters
- `sounding_us_standard_atm` — Reference atmospheric sounding (U.S. Standard Atmosphere)
- `check_flux_accuracy.py` — Python validation script

## Running the Test

```bash
cd Exec/CanonicalTests/Radiation/SW_ClearSky_Analytical
mpirun -np 1 erf.ex inputs
python3 check_flux_accuracy.py
```

## Validation Criteria

The checker script verifies:

1. **TOA flux** matches S₀ · cos(zenith) to within 1%
2. **Surface flux** matches analytical attenuation through 20 layers to within 5%
3. **All fluxes are non-negative** (physical constraint)
4. **Diagnostics file created** (`radiation_sw_diag.dat` or similar)
5. **No NaN or Inf values** in output

## Expected Output

- Radiation diagnostics file with columns: step, time, SW_surface, SW_TOA, heating rates
- CHECK PASS message if all validation criteria satisfied
- Heating rates should be negligible (clear-sky, non-absorbing)

## Related Documentation

- `RAD_DEVELOPMENT.md` — Base Two-Stream Solver section
- Beer-Lambert Law references in the main README
