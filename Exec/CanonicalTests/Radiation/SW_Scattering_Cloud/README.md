# Shortwave with Scattering Test

## Objective

Validate shortwave radiation scattering using the two-stream approximation for cloudy atmospheres.

## Test Purpose

This test confirms that the two-stream solver correctly:
- Computes diffuse scattering in the presence of cloud particles
- Distributes radiation between direct and diffuse components
- Maintains physically realistic scattering properties
- Produces appropriate cloud heating from scattering/absorption

## Test Design

### Configuration

- **Domain**: 1000 m × 1000 m horizontal, 10 km vertical (20 layers)
- **Time**: Single timestep, 0.1 s duration
- **Solar Constant**: S₀ = 1361 W/m²
- **Solar Zenith Angle**: 60°
- **Cloud Layer**: z = 4–6 km
- **Cloud Optical Depth**: τ_cloud = 1.0 (scattering cloud)
- **Single-Scattering Albedo**: ω₀ ≈ 0.8–0.9 (typical for water clouds)

### Key Physics

In the two-stream approximation, total radiation splits into:
```
F_total(z) = F_direct(z) + F_diffuse(z)
```

with scattering as:
```
F_diffuse ∝ (1 - ω₀) · τ · I
```

where ω₀ is the single-scattering albedo (fraction scattered vs. absorbed).

## Files

- `inputs` — Main control file with scattering optical-depth and albedo parameters
- `sounding_us_standard_atm` — Reference atmospheric sounding
- `check_flux_accuracy.py` — Python validation script

## Running the Test

```bash
cd Exec/CanonicalTests/Radiation/SW_Scattering_Cloud
mpirun -np 1 erf.ex inputs
python3 check_flux_accuracy.py
```

## Validation Criteria

The checker script verifies:

1. **Diffuse flux component** is non-zero in and near cloud layer
2. **Total flux divergence** matches absorption + scattering
3. **All fluxes remain non-negative** (fundamental constraint)
4. **Surface flux** reduced from clear-sky (cloud reflectance)
5. **Diagnostics file created** with expected structure
6. **No NaN or Inf values** in output
7. **Heating rate signature** consistent with cloud scattering/absorption

## Expected Output

- Radiation diagnostics file with significant diffuse component in cloud layer
- CHECK PASS message confirming scattering processing
- Reduced direct flux but elevated diffuse flux in cloudy region
- Reasonable heating-rate profile

## Related Documentation

- `RAD_DEVELOPMENT.md` — Scattering section
- Meador-Weaver (1980) for two-stream scattering formulation
- Main README for related cloud tests
