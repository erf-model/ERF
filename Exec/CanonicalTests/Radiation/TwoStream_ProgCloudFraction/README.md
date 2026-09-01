# Prognostic Cloud Fraction Test

## Objective

Validate that cloud fraction is diagnosed from relative humidity and cloud liquid water, and correctly applied to optical-depth masking.

## Test Purpose

This test confirms that:
- Cloud fraction is computed from RH/qc without explicit input
- Cloud-layer opacity masks only cloudy portions of each column
- Heating rates are reduced in non-cloudy regions
- Temporal smoothing (if enabled) prevents spurious oscillations
- Backward compatibility maintained (disabled by default)

## Test Design

### Configuration

- **Domain**: 1000 m × 1000 m horizontal, 10 km vertical (20 layers)
- **Time**: Multiple timesteps (5+ steps)
- **Sounding**: Includes moisture profile with RH variations
- **Cloud Optical Depth**: τ_cloud = 1.0 (where clouds form)
- **Cloud Fraction**: Diagnosed from RH relative to saturation and qc
- **Prognostic Cloud Fraction**: Enabled with smoothing timescale ~10 min

### Key Physics

Cloud fraction is diagnosed as:
```
cf = min(1.0, max(0.0, (RH - RH_threshold) / (1.0 - RH_threshold)))
```

or from qc directly:
```
cf = min(1.0, qc / qc_threshold)
```

Total optical depth becomes:
```
τ(k) = τ_bg(k) + cf(k) · τ_cloud
```

This ensures that dry layers contribute no cloud opacity.

## Files

- `inputs` — Main control file with prognostic cloud fraction enabled
- `input_sounding_moist` — Atmospheric sounding with moisture profile
- `check_cloud_fraction_accuracy.py` — Python validation script

## Running the Test

```bash
cd Exec/CanonicalTests/Radiation/TwoStream_ProgCloudFraction
mpirun -np 1 erf.ex inputs
python3 check_cloud_fraction_accuracy.py
```

## Validation Criteria

The checker script verifies:

1. **Cloud fraction** is bounded 0 ≤ cf ≤ 1 everywhere
2. **Cloud fraction** is highest where RH is highest (physical consistency)
3. **Heating** reduced in low-cf regions, enhanced in high-cf regions
4. **Temporal smoothing** prevents spurious oscillations between steps
5. **Diagnostics file** includes cloud fraction or heating diagnostics
6. **No NaN or Inf values** in output
7. **Backward compatibility** (disabled by default, no change to clear-sky case)

## Expected Output

- Radiation diagnostics with cloud-layer signature heating
- CHECK PASS message confirming prognostic cloud fraction
- Heating profile showing maximum in cloudy regions
- Smooth evolution of cloud fraction over multiple steps

## Related Documentation

- `RAD_DEVELOPMENT.md` — Prognostic Cloud Fraction section
- Cloud optical-depth formulation in `Source/Radiation/`
