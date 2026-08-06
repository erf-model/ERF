# ABL_with_MRF

## Purpose
These cases combine the fire module with multiple atmospheric boundary-layer soundings and MRF settings to test sensitivity across neutral, unstable, diurnal, and shear-driven conditions.

## Physics / Model Features Exercised
- Fire spread / ignition configuration
- Atmospheric forcing and boundary-condition setup

## Expected Results
See the input-file header comments in this directory for the specific validation target. In general, these cases should reproduce the documented analytical trend, qualitative regime change, or engineering diagnostic associated with the scenario.

## Key Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| `erf.prob_name` | `"ABL"` | Primary configuration value taken from `inputs_fire_abl_diurnal`. |
| `max_step` | `20` | Primary configuration value taken from `inputs_fire_abl_diurnal`. |
| `amrex.fpe_trap_invalid` | `0` | Primary configuration value taken from `inputs_fire_abl_diurnal`. |
| `fabarray.mfiter_tile_size` | `1024 1024 1024` | Primary configuration value taken from `inputs_fire_abl_diurnal`. |
| `geometry.prob_extent` | `8000  8000  2048` | Primary configuration value taken from `inputs_fire_abl_diurnal`. |
| `amr.n_cell` | `8   8   64` | Primary configuration value taken from `inputs_fire_abl_diurnal`. |
| `geometry.is_periodic` | `1 1 0` | Primary configuration value taken from `inputs_fire_abl_diurnal`. |
| `zlo.type` | `"surface_layer"` | Primary configuration value taken from `inputs_fire_abl_diurnal`. |

## References
- Rothermel 1972, A Mathematical Model for Predicting Fire Spread in Wildland Fuels.
- Andrews 2018, The Rothermel surface fire spread model and associated developments.
