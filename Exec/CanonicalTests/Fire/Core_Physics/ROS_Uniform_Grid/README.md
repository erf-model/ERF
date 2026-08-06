# ROS_Uniform_Grid

## Purpose
This case verifies that fire spread on a flat, uniform grid behaves consistently with the analytical/unit-test benchmark. It is useful for separating physics verification from mesh or terrain complications.

## Physics / Model Features Exercised
- Fire spread / ignition configuration
- Atmospheric forcing and boundary-condition setup

## Expected Results
See the input-file header comments in this directory for the specific validation target. In general, these cases should reproduce the documented analytical trend, qualitative regime change, or engineering diagnostic associated with the scenario.

## Key Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| `geometry.prob_lo` | `0.0 0.0 0.0` | Primary configuration value taken from `inputs_fire_flat_uniform`. |
| `geometry.prob_hi` | `2000.0 2000.0 500.0` | Primary configuration value taken from `inputs_fire_flat_uniform`. |
| `amr.n_cell` | `40 40 50` | Primary configuration value taken from `inputs_fire_flat_uniform`. |
| `geometry.is_periodic` | `0 0 0` | Primary configuration value taken from `inputs_fire_flat_uniform`. |
| `amr.max_level` | `0` | Primary configuration value taken from `inputs_fire_flat_uniform`. |
| `amr.max_grid_size` | `100` | Primary configuration value taken from `inputs_fire_flat_uniform`. |
| `amr.max_grid_size_z` | `50` | Primary configuration value taken from `inputs_fire_flat_uniform`. |
| `zlo.type` | `"surface_layer"` | Primary configuration value taken from `inputs_fire_flat_uniform`. |

## References
- Rothermel 1972, A Mathematical Model for Predicting Fire Spread in Wildland Fuels.
- Andrews 2018, The Rothermel surface fire spread model and associated developments.
