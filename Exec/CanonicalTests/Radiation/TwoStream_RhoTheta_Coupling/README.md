# RhoTheta Coupling Test

## Objective

Validate that radiative heating computed by the two-stream solver is correctly deposited into the model's thermodynamic equation via the `qheating_rates` pathway.

## Test Purpose

This test confirms that:
- Heating rates are computed by TwoStream
- Heating rates are stored in the `qheating_rates` array
- Heating rates are correctly injected into `RhoTheta` during the slow-step thermodynamic update
- No heating is lost or duplicated during the coupling
- Backward compatibility is maintained (heating only applied when TwoStream is enabled)

## Test Design

### Configuration

- **Domain**: 1000 m × 1000 m horizontal, 10 km vertical (20 layers)
- **Time**: Multiple timesteps (5 × 0.5 s each)
- **Solar Constant**: S₀ = 1361 W/m²
- **Solar Zenith Angle**: 60°
- **Optical Depth**: τ = 0.05 per layer (clear-sky)
- **Thermodynamic Integration**: Enabled (heating applied to RhoTheta)

### Key Physics

The radiative heating rate dT/dt is converted to RhoTheta tendency:
```
(dRhoTheta/dt)_rad = Rho · Cp · dT/dt
```

where Cp is the specific heat at constant pressure.

Over the test duration, cumulative heating should:
1. Increase potential temperature throughout the domain
2. Show maximum heating near the surface (strongest solar flux)
3. Produce a consistent vertical profile each timestep

## Files

- `inputs` — Main control file with TwoStream and thermodynamic coupling enabled
- `input_sounding_phase5_coupling` — Reference atmospheric sounding
- `check_flux_accuracy.py` — Python validation script

## Running the Test

```bash
cd Exec/CanonicalTests/Radiation/TwoStream_RhoTheta_Coupling
mpirun -np 1 erf.ex inputs
python3 check_flux_accuracy.py
```

## Validation Criteria

The checker script verifies:

1. **qheating_rates** field is non-zero in radiation-enabled case
2. **qheating_rates** is zero when radiation disabled (backward compatibility)
3. **Temperature tendency** from heating matches flux divergence
4. **Cumulative temperature change** is monotonically increasing over timesteps
5. **Vertical profile** of heating shows maximum near surface
6. **No spurious oscillations** in heating rates
7. **Energy conservation** check (total heat input ≈ temperature change)

## Expected Output

- Radiation diagnostics with heating rate values
- Plotfile outputs showing temperature increase due to radiation
- CHECK PASS message confirming proper coupling
- Heating-rate profiles sensible (peak near surface, decrease with altitude)

## Related Documentation

- `RAD_DEVELOPMENT.md` — RhoTheta Coupling section
- `Source/SourceTerms/ERF_MakeSources.cpp` — where `qheating_rates` is added to the RhoTheta source
