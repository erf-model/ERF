# Longwave Isothermal Consistency Test

## Objective

Validate longwave (thermal) radiation in a special isothermal configuration where net heating is zero. This serves as a stringent self-consistency test for the two-stream solver.

## Test Purpose

This test confirms that the longwave solver correctly:
- Computes upwelling and downwelling fluxes from surface temperature
- Maintains energy conservation in isothermal conditions
- Produces zero heating rates when atmosphere and surface are at same temperature
- Preserves isothermal state (no spurious temperature changes)

## Test Design

### Configuration

- **Domain**: 1000 m × 1000 m horizontal, 10 km vertical (20 layers)
- **Time**: Single timestep, 0.1 s duration
- **Isothermal Temperature**: T_iso = 288.15 K (everywhere)
- **Longwave Optical Depth per Layer**: τ_lw = 1.0 (uniform)
- **Isothermal Test Mode**: ENABLED (forces all T = T_iso)

### Key Physics

In isothermal equilibrium, all layers at the same temperature must have identical upwelling and downwelling fluxes:

```
F_up(z) = σ · T_iso⁴  (constant everywhere)
F_down(z) = σ · T_iso⁴  (constant everywhere)
F_net = F_up - F_down = 0  (perfect balance)
dT/dt = 0  (no heating)
```

where σ = 5.670374419 × 10⁻⁸ W/(m²·K⁴) is the Stefan-Boltzmann constant.

For T_iso = 288.15 K:
```
F = 5.670e-8 × (288.15)⁴ ≈ 384.4 W/m²
```

### Physical Insight

This test is **extremely stringent**: if the two-stream solver has any bug in its vertical integration, coupling, or symmetry handling, isothermal conditions will be violated. Zero heating rates are impossible to achieve by accident.

## Files

- `inputs` — Main control file with isothermal test mode enabled
- `sounding_us_standard_atm` — Reference atmospheric sounding
- `check_flux_accuracy.py` — Python validation script

## Running the Test

```bash
cd Exec/CanonicalTests/Radiation/LW_Isothermal
mpirun -np 1 erf.ex inputs
python3 check_flux_accuracy.py
```

## Validation Criteria

The checker script verifies:

1. **F_up matches σ · T_iso⁴** to machine precision (< 1e-10 relative error)
2. **F_down matches F_up** to machine precision
3. **F_net is exactly zero** (within rounding)
4. **Heating rates are zero** everywhere (< 1e-12 K/s)
5. **Energy conservation holds** globally
6. **Diagnostics file created** with finite values

## Expected Output

- Radiation diagnostics file with identical up and down fluxes at all levels
- CHECK PASS message confirming isothermal consistency
- No spurious heating anywhere in the domain

## Related Documentation

- `RAD_DEVELOPMENT.md` — Base Two-Stream Solver section
- Stefan-Boltzmann Law references in the main README
- Toon et al. (1989) for two-stream formulation
