# DustTerrainFlat

## Purpose
This case provides a terrain-free hazard reference for comparing the effects of hills, pits, and slopes on dust transport and exposure.

## Physics / Model Features Exercised
- Dust emission / transport / deposition controls
- Exposure or air-quality diagnostics as configured
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
| `geometry.prob_extent` | `3000 3000 1024` | Primary configuration value taken from `inputs`. |
| `amr.n_cell` | `8 8 64` | Primary configuration value taken from `inputs`. |
| `geometry.is_periodic` | `1 1 0` | Primary configuration value taken from `inputs`. |
| `amr.max_grid_size` | `32` | Primary configuration value taken from `inputs`. |

## References
- Bagnold 1941, The Physics of Blown Sand and Desert Dunes.
- Marticorena and Bergametti 1995, Modeling the atmospheric dust cycle.
- Analytical Gaussian terrain idealization used for topographic verification.
