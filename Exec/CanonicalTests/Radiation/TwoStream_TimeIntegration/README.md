# Time Integration Consistency Test

## Objective

Validate that repeated calls to the TwoStream radiation solver maintain temporal consistency and diagnostic fidelity.

## Test Purpose

This test confirms that:
- Radiation solver produces consistent results when called repeatedly
- Diagnostic output tracks call-site identity correctly
- Time-integration cadence does not introduce spurious oscillations or drifts
- Multiple slow-step calls per fast-step produce expected heating accumulation
- Backward compatibility maintained (output independent of call pattern when disabled)

## Test Design

### Configuration

- **Domain**: 1000 m × 1000 m horizontal, 10 km vertical (20 layers)
- **Time**: Multiple timesteps (10+ steps with varying call patterns)
- **Solar Constant**: S₀ = 1361 W/m²
- **Solar Zenith Angle**: 60° (fixed for consistency)
- **Optical Depth**: τ = 0.05 per layer
- **Call Patterns**: Pre-step only, both pre/post, diagnostic tracking enabled

### Key Physics

Heating accumulated per call should follow:
```
Q_total = Q_per_call × N_calls_per_timestep
```

If the solver is called twice per timestep (once pre-step, once post-step), total heating should be approximately 2× that of a single-call case.

## Files

- `inputs` — Main control file with TwoStream enabled
- `input_sounding_phase6_timing` — Reference atmospheric sounding
- `check_timing_consistency.py` — Python validation script for temporal consistency
- `radiation_phase6_timing_diag.dat` — Expected diagnostic output (reference)

## Running the Test

```bash
cd Exec/CanonicalTests/Radiation/TwoStream_TimeIntegration
mpirun -np 1 erf.ex inputs
python3 check_timing_consistency.py
```

## Validation Criteria

The checker script verifies:

1. **Call-site diagnostics** correctly identify pre-step vs. post-step calls
2. **Heating accumulation** scales linearly with number of calls
3. **Flux values** are consistent across all calls at the same time
4. **No temporal drift** in fluxes (same solar geometry, same result)
5. **Diagnostics file** parsing and expected structure
6. **No NaN or Inf values** in output
7. **Multi-call pattern** produces expected total heating (sum of parts)

## Expected Output

- Multiple diagnostic entries per timestep (one per call)
- CHECK PASS message confirming temporal consistency
- Heating and flux values repeatable when called at same time
- Clear diagnostic differentiation between pre and post calls

## Related Documentation

- `RAD_DEVELOPMENT.md` — Time Integration and Diagnostics Cadence section
- Main README for related diagnostic/timing tests
