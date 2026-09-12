# Shortwave with Cloud Layer Test

## Objective

Validate shortwave radiation in the presence of a prescribed cloud optical-depth layer.

## Test Purpose

This test confirms that the two-stream solver correctly:
- Applies height-dependent cloud optical-depth contributions
- Computes flux attenuation through cloudy regions
- Maintains backward compatibility when cloud parameters are disabled
- Produces realistic heating rates from cloud absorption/scattering

## Test Design

### Configuration

- **Domain**: 1000 m × 1000 m horizontal, 10 km vertical (20 layers)
- **Time**: Single timestep, 0.1 s duration
- **Solar Constant**: S₀ = 1361 W/m²
- **Solar Zenith Angle**: 60°
- **Background Optical Depth**: τ_bg = 0.05 per layer (clear-sky baseline)
- **Cloud Layer**: Defined between z = 4 km and z = 6 km
- **Cloud Optical Depth**: τ_cloud = 1.0 in cloudy region, 0.0 elsewhere

### Key Physics

Total optical depth in cloudy layers:
```
τ(k) = τ_background(k) + τ_cloud(k)
```

This results in:
- Strong attenuation in the cloud layer (higher optical depth)
- Flux divergence inside the cloud (heating source)
- Reduced surface flux compared to clear-sky case

## Files

- `inputs` — Main control file with cloud optical-depth parameters
- `sounding_us_standard_atm` — Reference atmospheric sounding
- `check_flux_accuracy.py` — Python validation script

## Running the Test

```bash
cd Exec/CanonicalTests/Radiation/SW_Cloud_Layer
mpirun -np 1 erf.ex inputs
python3 check_flux_accuracy.py
```

## Validation Criteria

The checker script verifies:

1. **Cloudy layer flux** is significantly reduced from clear-sky baseline
2. **Maximum heating rate** occurs near cloud top/bottom (strong gradient)
3. **All fluxes remain non-negative** and physically reasonable
4. **Diagnostics file created** with expected structure
5. **No NaN or Inf values** in output
6. **Flux profile** shows signature cloud-layer structure

## Expected Output

- Radiation diagnostics file showing reduced surface flux (cloud opacity)
- CHECK PASS message confirming cloud layer processing
- Heating-rate maximum at cloud boundaries
- Clear distinction between clear and cloudy regions

## Related Documentation

- `RAD_DEVELOPMENT.md` — Cloud Optical Depth section
- Main README table linking to related cloud tests
