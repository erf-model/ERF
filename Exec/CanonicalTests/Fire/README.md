# Fire Model Tests

This directory contains test cases for the 2D fire model implementation in ERF.

## Fire Model Overview

The fire model implements a 2D wildfire propagation solver using:
- **Rothermel Fire Spread Model**: Computes rate of fire spread based on fuel characteristics, weather, and topography
- **FARSITE Elliptical Expansion**: Models fire expansion as an ellipse with orientation based on wind and slope

## Rothermel C Coefficient Correction (Critical Fix)

### Problem
The original implementation had an incorrect wind factor coefficient C in Rothermel Equation 47:
```cpp
// WRONG (produced C ≈ 1.74e-5 for FM1):
Real C = 7.47 * std::exp(-0.133 * std::pow(sigma, 0.55));
```

This caused:
- C coefficient to be near-zero (1.74e-5) instead of ~7.4
- MEWS wind speed cap to be unrealistically low (~4 ft/min instead of ~300+ ft/min)
- Maximum Rate of Spread (ROS) to be clamped to unrealistic values
- Fire spread rates far too slow compared to physical models

### Solution
The correct Rothermel (1972) Equation 47 uses **negative exponent** on sigma:
```cpp
// CORRECT (produces C ≈ 7.4 for FM1):
Real C = 7.47 * std::exp(-0.8711 * std::pow(sigma, -0.55));
```

### Verification
For standard Anderson FBFM13 fuel models:
- **FM1 (Short Grass)**: σ=3500 ft⁻¹ → C = 7.40
- **FM4 (Chaparral)**: σ=1739 ft⁻¹ → C = 7.37
- **FM3 (Tall Grass)**: σ=1500 ft⁻¹ → C = 7.35

These match published BEHAVE/BehavePlus reference values.

## Quick Validation

### Run Python Unit Tests (No ERF Required)
```bash
cd Exec/CanonicalTests/Fire/
python3 test_rothermel_unit.py              # Pure Python physics check (~1 sec)
python3 test_fire_ros_regression.py         # Reference value regression check (~1 sec)
```

### Run Integration Test (Requires Compiled ERF)
```bash
cd Exec/CanonicalTests/Fire/
python3 test_rothermel_ros.py --erf_exe ./erf --input_file inputs_fire_flat_uniform
```

## Test Cases

### Existing Test Cases

#### inputs_fire_dummy
Dummy test case that verifies fire model initialization and basic function calls.

#### inputs_fire_phase2
Phase 2 implementation test on a flat domain using:
- GR1 fuel model (short grass) from Anderson FBFM13
- 5 m/s uniform wind
- 8% fuel moisture content
- Rothermel fire spread rate computation

#### inputs_fire_phase2_slope
Phase 2 test case with terrain effects:
- 30-degree slope
- 3 m/s upslope wind
- Terrain wind corrections enabled (ridges, valleys, sheltering)
- Demonstrates fire spread on sloped terrain

#### inputs_fire_phase3
Phase 3 test case with FARSITE level-set propagation:
- Elliptical fire expansion
- Level-set method for fire front advancement
- Wind and slope-based propagation direction

#### inputs_fire_abl_mrf_unstable
Atmospheric Boundary Layer test case with Mesoscale Radiation Framework (MRF):
- Couples MRF boundary layer physics with fire model
- Uses realistic atmospheric initialization from input sounding
- Includes Coriolis effects and geostrophic wind forcing
- Flat terrain for clean fire-atmosphere coupling testing
- Uses MRF turbulence closure for stable boundary layer conditions

### New Regression Test Cases

#### inputs_fire_flat_uniform
Minimal flat uniform mesh regression test:
- FM1 (Short Grass) fuel
- Flat terrain, uniform atmospheric grid
- 5 m/s nominal wind
- 8% fuel moisture
- **Expected ROS**: 0.1-1.0 m/s over 1800s simulation
- **Purpose**: Baseline regression test for flat domain fire spread

#### inputs_fire_vertical_refinement
Vertical grid stretching regression test:
- FM1 (Short Grass) fuel
- Stretched vertical grid (grid_stretching_ratio=1.08)
- Wind interpolation at z_ref=6.1 m with non-uniform dz
- **Expected ROS**: 0.1-1.0 m/s
- **Purpose**: Verify wind extraction and ROS computation with stretched grids

#### inputs_fire_phase5
Phase 5 test case with heat flux and diagnostics:
- FM1 (Short Grass) fuel
- 5 m/s nominal wind
- 8% fuel moisture
- Heat flux computation via Albini exponential burnout model
- Fuel load depletion via forward Euler
- Byram fireline intensity and Thomas flame length diagnostics
- CSV output with max heat flux tracking
- **Expected outputs**: fire_heat_flux > 0 in burned cells, fuel_load decreasing, flame length consistent with Byram
- **Purpose**: Regression test for heat flux, fuel depletion, fireline intensity, and flame length

## Test Summary Table

| Test Name | Type | Fuel | Wind | Purpose | Command |
|-----------|------|------|------|---------|---------|
| test_rothermel_unit.py | Unit | FM1, FM4 | Various | Physics equations validation (no ERF) | `python3 test_rothermel_unit.py` |
| test_fire_ros_regression.py | Regression | FM1, FM4 | Various | Reference value checks (no ERF) | `python3 test_fire_ros_regression.py` |
| test_rothermel_ros.py | Integration | FM1 | Uniform | ERF ROS magnitude check | `python3 test_rothermel_ros.py --erf_exe ./erf --input_file inputs_fire_flat_uniform` |
| inputs_fire_flat_uniform | Regression | FM1 | 5 m/s | Flat uniform mesh fire spread | `./erf inputs_fire_flat_uniform` |
| inputs_fire_vertical_refinement | Regression | FM1 | 5 m/s | Stretched vertical grid fire spread | `./erf inputs_fire_vertical_refinement` |
| inputs_fire_phase5 | Regression | FM1 | 5 m/s | Heat flux and diagnostics | `./erf inputs_fire_phase5` |
| inputs_fire_dummy | Smoke | FM1 | Dummy | Initialization and basic calls | `./erf inputs_fire_dummy` |

## Output Variables

All fire test cases output 2D fire-specific variables at each time step:
- `fire_phi`: Level-set function (<0 burned, >0 unburned)
- `fire_ros`: Rate of spread [m/s]
- `fire_wind_eff`: Effective wind after WAF and terrain corrections [m/s]
- `fire_wind_ref`: Reference height wind (6.1 m) [m/s]
- `fire_slopes`: Terrain slope components (dz/dx, dz/dy)
- `fire_fuel_mc`: Fuel moisture content (1hr, 10hr, 100hr)

Phase 5 additionally outputs:
- `fire_heat_flux`: Sensible heat flux [W/m²]
- `fire_fuel_load`: Remaining fuel load [kg/m²]
- `fire_fireline_intensity`: Byram fireline intensity [kW/m]
- `fire_flame_length`: Thomas flame length [m]

Output is written to `plt2d_fire*` directories at intervals specified by `erf.plot2d_int_1`.

## Running Tests

To run a fire model test:

```bash
./erf inputs_fire_phase2
```

Plotfiles will be written to `plt2d_fire*` directories in the run directory.

## Input Parameters

Key fire model parameters:

- `erf.fire.enable`: Enable fire module (true/false)
- `erf.fire.fuel_model_id`: Fuel model number (1=GR1, etc.)
- `erf.fire.moisture_1hr`, `erf.fire.moisture_10hr`, `erf.fire.moisture_100hr`: Fuel moisture content
- `erf.fire.ignition_x`, `erf.fire.ignition_y`, `erf.fire.ignition_r`: Ignition location and radius
- `erf.fire.wind_ref_ht`: Reference height for wind extraction
- `erf.fire.use_waf`: Apply Wind Adjustment Factor for vegetation
- `erf.fire.use_terrain_wind`: Enable terrain wind corrections
- `erf.fire.grid_ratio`: Fire grid refinement ratio relative to atmospheric grid

## Debug Output

When `erf.fire.fire_debug = true`, the fire model outputs:
- **Rothermel parameters**: R0, I_R, C, B, E, β, U_max for each fuel/moisture combination
- **Fire grid resolution**: dx, dy, cell counts for the refined fire grid
- **Wind extraction**: Maximum reference height wind, effective wind after corrections
- **Rate of spread**: Maximum ROS at each time step
- **Ignition state**: Fire front tracking and level-set function statistics

Example debug output:
```
[FIRE] FireLayer initialized: C=5, fuel_model=1, grid=4 boxes
[FIRE DEBUG] Fire grid mesh resolution: dx=8 m, dy=8 m, extent: 200 x 200 cells, grid_ratio=5
Rothermel Params: R0=2.5 ft/min (0.01270 m/s) I_R=342 BTU/ft²/min C=7.40 B=0.0877 E=0.8734 ...
[FIRE DEBUG] Wind extraction completed. Max reference wind: 5.02 m/s
[FIRE DEBUG] Rate-of-spread computed. Max: 0.08523 m/s
```

