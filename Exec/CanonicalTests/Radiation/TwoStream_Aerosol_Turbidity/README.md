# Aerosol/Turbidity Test

## Objective

Validate prescribed aerosol or atmospheric turbidity optical-depth contributions to shortwave and longwave radiation.

## Test Purpose

This test confirms that:
- Aerosol optical depth is correctly added to background atmosphere
- Vertical aerosol profiles (constant/exponential/table-based) are applied
- Aerosol effects reduce surface flux and enhance heating
- Fallback behavior is robust when aerosol fields are unavailable
- Backward compatibility maintained (aerosol disabled by default)

## Test Design

### Configuration

- **Domain**: 1000 m × 1000 m horizontal, 10 km vertical (20 layers)
- **Time**: Single or multi-step run
- **Solar Constant**: S₀ = 1361 W/m²
- **Solar Zenith Angle**: 60°
- **Background Optical Depth**: τ_bg = 0.05 per layer (clear-sky)
- **Aerosol Profile**: Exponential decay from surface (e-folding scale ~2 km)
- **Aerosol Optical Depth at Surface**: τ_aero(0) = 0.2 (moderate turbidity)

### Key Physics

Total optical depth includes aerosol contribution:
```
τ(k) = τ_bg(k) + τ_aero(k) · exp(-z / H)
```

where H is the scale height (e.g., 2 km).

Enhanced optical depth reduces surface flux and produces heating aloft (aerosol absorption).

## Files

- `inputs` — Main control file with aerosol parameters enabled
- Sounding file — Reference atmospheric profile
- `check_aerosol_accuracy.py` — Python validation script

## Running the Test

```bash
cd Exec/CanonicalTests/Radiation/TwoStream_Aerosol_Turbidity
mpirun -np 1 erf.ex inputs
python3 check_aerosol_accuracy.py
```

## Validation Criteria

The checker script verifies:

1. **Surface flux** is reduced compared to clear-sky baseline
2. **Aerosol heating** is maximum near surface (peak optical depth)
3. **Total optical depth** reflects aerosol + background
4. **Vertical profile** of aerosol optical depth decreases with height
5. **Fallback path** exercises correctly when aerosol disabled
6. **Diagnostics file** includes aerosol contribution or heating
7. **No NaN or Inf values** in output
8. **Backward compatibility** (aerosol disabled by default)

## Expected Output

- Radiation diagnostics with reduced surface flux (aerosol opacity)
- CHECK PASS message confirming aerosol processing
- Heating rate maximum near surface (aerosol absorption zone)
- Clear contrast with clean-air baseline case

## Related Documentation

- `RAD_DEVELOPMENT.md` — Aerosol/Turbidity section
- Aerosol optical-depth formulation in `Source/Radiation/`
