# FireDustMassConservation

## Purpose
This hazard case verifies that fire-driven dust emission and transport remain mass-consistent in a coupled simulation, especially when dry deposition is disabled so airborne mass should not decrease spuriously.

## Physics / Model Features Exercised
- Fire spread / ignition configuration
- Atmospheric forcing and boundary-condition setup
- Coupled hazard-module configuration
- Cross-module diagnostics or interaction controls

## Expected Results
See the input-file header comments in this directory for the specific validation target. In general, these cases should reproduce the documented analytical trend, qualitative regime change, or engineering diagnostic associated with the scenario.

## Key Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| `erf.prob_name` | `"ABL"` | Primary configuration value taken from `inputs`. |
| `max_step` | `100` | Primary configuration value taken from `inputs`. |
| `amr.max_level` | `0` | Primary configuration value taken from `inputs`. |
| `amrex.fpe_trap_invalid` | `0` | Primary configuration value taken from `inputs`. |
| `geometry.prob_extent` | `3000 3000 1024` | Primary configuration value taken from `inputs`. |
| `amr.n_cell` | `8 8 64` | Primary configuration value taken from `inputs`. |
| `geometry.is_periodic` | `1 1 0` | Primary configuration value taken from `inputs`. |
| `amr.max_grid_size` | `32` | Primary configuration value taken from `inputs`. |

## References
- Rothermel 1972, A Mathematical Model for Predicting Fire Spread in Wildland Fuels.
- Andrews 2018, The Rothermel surface fire spread model and associated developments.
