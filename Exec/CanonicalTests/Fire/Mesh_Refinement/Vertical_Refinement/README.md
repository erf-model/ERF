# Vertical_Refinement

## Purpose
This case documents a fire configuration with enhanced vertical resolution, used to confirm that the model setup and boundary conditions remain valid under refined z discretization.

## Physics / Model Features Exercised
- Fire spread / ignition configuration
- Atmospheric forcing and boundary-condition setup

## Expected Results
See the input-file header comments in this directory for the specific validation target. In general, these cases should reproduce the documented analytical trend, qualitative regime change, or engineering diagnostic associated with the scenario.

## Key Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| `geometry.prob_lo` | `0.0 0.0 0.0` | Primary configuration value taken from `inputs_fire_vertical_refinement`. |
| `geometry.prob_hi` | `2000.0 2000.0 1000.0` | Primary configuration value taken from `inputs_fire_vertical_refinement`. |
| `amr.n_cell` | `40 40 64` | Primary configuration value taken from `inputs_fire_vertical_refinement`. |
| `geometry.is_periodic` | `0 0 0` | Primary configuration value taken from `inputs_fire_vertical_refinement`. |
| `amr.max_level` | `0` | Primary configuration value taken from `inputs_fire_vertical_refinement`. |
| `amr.max_grid_size` | `100` | Primary configuration value taken from `inputs_fire_vertical_refinement`. |
| `amr.max_grid_size_z` | `64` | Primary configuration value taken from `inputs_fire_vertical_refinement`. |
| `erf.initial_dz` | `4.0` | Primary configuration value taken from `inputs_fire_vertical_refinement`. |

## References
- Rothermel 1972, A Mathematical Model for Predicting Fire Spread in Wildland Fuels.
- Andrews 2018, The Rothermel surface fire spread model and associated developments.
