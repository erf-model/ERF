# ROS_Slope_Effects

## Purpose
This case puts the fire on a resolved hill and checks that the slope term of the rate-of-spread model accelerates it relative to the flat-ground baseline. It is the only fire canonical case that configures terrain, so it also covers the terrain-following paths: slope computation from `z_phys_nd`, wind extraction above true ground, and the projection of surface spread into map view.

## Terrain
A Gaussian ridge, `z(x) = 100 exp(-((x - 500)/150)^2)` metres, uniform in y and decaying to 1.5 mm at both x boundaries so it is compatible with the periodic domain. The crest is 100 m at x = 500 m and the peak slope is tan(theta) = 0.57, about 30 degrees, at x = 394 m and x = 606 m. Ignition sits on the steepest part of the windward flank, 100 m below the crest.

`terrain_gaussian_ridge.txt` is read twice: by `erf.terrain_file_name` for the atmospheric terrain-fitted mesh at 50 m, and by `erf.fire.terrain_file_name` so the fire grid resolves the same ridge at its own 10 m spacing. The wind extraction datum stays on the atmospheric column in either case, since the wind profile being interpolated belongs to that column.

## Physics / Model Features Exercised
- Rothermel slope factor, `phi_s = 5.275 beta^(-0.3) tan^2(theta)` (Eq. 51)
- Fire-grid slope computation from nodal terrain
- Terrain-following wind extraction and map-view projection of surface spread
- Fire spread / ignition configuration over a terrain-fitted mesh

## Expected Results
With no wind, the flat-ground rate of spread is the no-wind, no-slope value `R0`, and the ridge multiplies it by `1 + phi_s`:

| Configuration | Max fire-grid slope | Max ROS | `phi_s` / tan^2(theta) |
|---|---|---|---|
| Flat (terrain disabled) | 0.000 | 0.0202 m/s | - |
| Half-amplitude ridge (h = 50 m) | 0.285 | 0.0888 m/s | 41.9 |
| This case (h = 100 m) | 0.569 | 0.2942 m/s | 41.8 |

The slope factor tracks `tan^2(theta)` to within 0.2% across the two amplitudes, which is Eq. 51 with the FM1 packing ratio. The full 900 s run burns 1204 fire cells and crosses the crest without instability.

Compare directly against `ROS_Basic_Calculation` for the flat baseline.

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
| `erf.terrain_type` | `StaticFittedMesh` | Terrain-fitted atmospheric mesh over the ridge. |
| `erf.terrain_file_name` | `"terrain_gaussian_ridge.txt"` | Ridge for the atmospheric grid. |
| `erf.fire.terrain_file_name` | `"terrain_gaussian_ridge.txt"` | Same ridge at fire-grid resolution. |
| `erf.fire.coupling_type` | `"passive"` | ROS test; no fire-atmosphere feedback. |
| `erf.fire.use_terrain_wind` | `false` | Terrain flow is resolved, so the empirical corrections would double count. |
| `erf.fire.ignition_x` | `400.0` | On the steepest part of the windward flank. |

## References
- Rothermel 1972, A Mathematical Model for Predicting Fire Spread in Wildland Fuels.
- Andrews 2018, The Rothermel surface fire spread model and associated developments.
