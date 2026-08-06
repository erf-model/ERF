# DustAtmReturn

## Purpose
This case documents the return-coupling pathway in which dust can feed back or be reintroduced according to the configured return logic.

## Physics / Model Features Exercised
- Dust emission / transport / deposition controls
- Exposure or air-quality diagnostics as configured

## Expected Results
See the input-file header comments in this directory for the specific validation target. In general, these cases should reproduce the documented analytical trend, qualitative regime change, or engineering diagnostic associated with the scenario.

## Key Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| `erf.prob_name` | `"ABL"` | Primary configuration value taken from `inputs`. |
| `max_step` | `5` | Primary configuration value taken from `inputs`. |
| `amr.max_level` | `0` | Primary configuration value taken from `inputs`. |
| `amrex.fpe_trap_invalid` | `0` | Primary configuration value taken from `inputs`. |
| `fabarray.mfiter_tile_size` | `1024 1024 1024` | Primary configuration value taken from `inputs`. |
| `geometry.prob_extent` | `3000 3000 1024` | Primary configuration value taken from `inputs`. |
| `amr.n_cell` | `8    8   64` | Primary configuration value taken from `inputs`. |
| `geometry.is_periodic` | `1 1 0` | Primary configuration value taken from `inputs`. |

## References
- Bagnold 1941, The Physics of Blown Sand and Desert Dunes.
- Marticorena and Bergametti 1995, Modeling the atmospheric dust cycle.
