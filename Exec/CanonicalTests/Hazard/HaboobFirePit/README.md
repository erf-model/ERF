# HaboobFirePit

## Purpose
This case uses bowl-shaped terrain to test confined or sheltered haboob-fire interactions and resulting hazard evolution.

## Physics / Model Features Exercised
- Fire spread / ignition configuration
- Atmospheric forcing and boundary-condition setup
- Coupled hazard-module configuration
- Cross-module diagnostics or interaction controls
- Idealized terrain effects

## Expected Results
See the input-file header comments in this directory for the specific validation target. In general, these cases should reproduce the documented analytical trend, qualitative regime change, or engineering diagnostic associated with the scenario.

## Key Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| `erf.prob_name` | `"ABL"` | Primary configuration value taken from `inputs`. |
| `max_step` | `20` | Primary configuration value taken from `inputs`. |
| `amr.max_level` | `0` | Primary configuration value taken from `inputs`. |
| `amrex.fpe_trap_invalid` | `0` | Primary configuration value taken from `inputs`. |
| `fabarray.mfiter_tile_size` | `1024 1024 1024` | Primary configuration value taken from `inputs`. |
| `geometry.prob_extent` | `4000.0 4000.0 1500.0` | Primary configuration value taken from `inputs`. |
| `amr.n_cell` | `32 32 64` | Primary configuration value taken from `inputs`. |
| `geometry.is_periodic` | `1 1 0` | Primary configuration value taken from `inputs`. |

## References
- Rothermel 1972, A Mathematical Model for Predicting Fire Spread in Wildland Fuels.
- Andrews 2018, The Rothermel surface fire spread model and associated developments.
- Density-current / haboob literature for idealized cold-pool initialization.
- Analytical Gaussian terrain idealization used for topographic verification.
