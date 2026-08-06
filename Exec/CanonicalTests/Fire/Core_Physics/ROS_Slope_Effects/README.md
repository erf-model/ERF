# ROS_Slope_Effects

## Purpose
This case checks that the fire model responds correctly to inclined terrain or slope forcing by increasing upslope spread relative to the flat-terrain baseline. It exercises the slope term in the ROS formulation while keeping the rest of the setup simple.

## Physics / Model Features Exercised
- Fire spread / ignition configuration
- Atmospheric forcing and boundary-condition setup
- Idealized terrain effects

## Expected Results
See the input-file header comments in this directory for the specific validation target. In general, these cases should reproduce the documented analytical trend, qualitative regime change, or engineering diagnostic associated with the scenario.

## Key Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| `geometry.prob_lo` | `0.0 0.0 0.0` | Primary configuration value taken from `inputs_fire_phase2_slope`. |
| `geometry.prob_hi` | `1000.0 500.0 1000.0` | Primary configuration value taken from `inputs_fire_phase2_slope`. |
| `amr.n_cell` | `20 10 100` | Primary configuration value taken from `inputs_fire_phase2_slope`. |
| `geometry.is_periodic` | `1 1 0` | Primary configuration value taken from `inputs_fire_phase2_slope`. |
| `amr.max_level` | `0` | Primary configuration value taken from `inputs_fire_phase2_slope`. |
| `amr.max_grid_size` | `100` | Primary configuration value taken from `inputs_fire_phase2_slope`. |
| `amr.max_grid_size_z` | `100` | Primary configuration value taken from `inputs_fire_phase2_slope`. |
| `zlo.type` | `"surface_layer"` | Primary configuration value taken from `inputs_fire_phase2_slope`. |

## References
- Rothermel 1972, A Mathematical Model for Predicting Fire Spread in Wildland Fuels.
- Andrews 2018, The Rothermel surface fire spread model and associated developments.
