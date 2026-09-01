# Nonuniform Vertical Spacing Test

## Objective

Validate that radiative heating is computed correctly when physical vertical spacing (Δz) varies with height.

## Test Purpose

This test confirms that:
- Optical-depth computation accounts for varying layer thickness
- Heating rates scale appropriately with layer depth
- Vertical integration remains accurate with nonuniform grid
- Surface flux calculations are correct despite geometric variation
- Backward compatibility preserved (uniform grid produces original results)

## Test Design

### Configuration

- **Domain**: 1000 m × 1000 m horizontal, 10 km vertical (nonuniform layers)
- **Vertical Grid**: Stretched toward surface (finer near ground, coarser aloft)
  - Example: dz = 10–50 m near surface, dz = 500–1000 m at top
- **Time**: Single or multi-step run
- **Solar Constant**: S₀ = 1361 W/m²
- **Solar Zenith Angle**: 60°
- **Optical Depth**: τ_per_meter (interpreted per physical thickness)

### Key Physics

Optical depth scales with physical thickness:
```
τ(k) = τ_coefficient · dz(k)
```

In nonuniform grids, thin layers accumulate less optical depth, and thick layers accumulate more. The two-stream solver must correctly weight fluxes by actual layer thickness.

## Files

- `inputs` — Main control file with nonuniform geometry
- Sounding file — Reference atmospheric profile
- `check_nonuniform_accuracy.py` — Python validation script

## Running the Test

```bash
cd Exec/CanonicalTests/Radiation/TwoStream_NonuniformDZ
mpirun -np 1 erf.ex inputs
python3 check_nonuniform_accuracy.py
```

## Validation Criteria

The checker script verifies:

1. **Fluxes computed correctly** despite nonuniform grid
2. **Heating rate magnitude** scales with layer thickness (thin layers → small dT/dt)
3. **Cumulative optical depth** reflects actual physical path
4. **Surface flux** correct for nonuniform attenuation
5. **No spurious artifacts** at grid refinement boundaries
6. **Smooth heating profile** without jumps/discontinuities
7. **Diagnostics file** created and finite

## Expected Output

- Radiation diagnostics showing heating variations consistent with grid spacing
- CHECK PASS message confirming nonuniform-grid correctness
- Heating rates lower in thin layers, higher in thick layers (when normalized by layer depth)

## Related Documentation

- `RAD_DEVELOPMENT.md` — Dynamic Optical Depth section
- Grid geometry handling in `Source/Radiation/`
