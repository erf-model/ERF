# Fire Model Tests

This directory contains test cases for the 2D fire model implementation in ERF.

## Fire Model Overview

The fire model implements a 2D wildfire propagation solver using:
- **Rothermel Fire Spread Model**: Computes rate of fire spread based on fuel characteristics, weather, and topography
- **FARSITE Elliptical Expansion**: Models fire expansion as an ellipse with orientation based on wind and slope

## Test Cases

### inputs_fire_dummy
Dummy test case that verifies fire model initialization and basic function calls.

### inputs_fire_phase2
Phase 2 implementation test on a flat domain using:
- GR1 fuel model (short grass) from Anderson FBFM13
- 5 m/s uniform wind
- 8% fuel moisture content
- Rothermel fire spread rate computation

### inputs_fire_phase2_slope
Phase 2 test case with terrain effects:
- 30-degree slope
- 3 m/s upslope wind
- Terrain wind corrections enabled (ridges, valleys, sheltering)
- Demonstrates fire spread on sloped terrain

### inputs_fire_phase3
Phase 3 test case with FARSITE level-set propagation:
- Elliptical fire expansion
- Level-set method for fire front advancement
- Wind and slope-based propagation direction

### inputs_fire_abl_mrf_unstable
Atmospheric Boundary Layer test case with Mesoscale Radiation Framework (MRF):
- Couples MRF boundary layer physics with fire model
- Uses realistic atmospheric initialization from input sounding
- Includes Coriolis effects and geostrophic wind forcing
- Flat terrain for clean fire-atmosphere coupling testing
- Uses MRF turbulence closure for stable boundary layer conditions

## Output Variables

All fire test cases output 2D fire-specific variables at each time step:
- `fire_phi`: Level-set function (<0 burned, >0 unburned)
- `fire_ros`: Rate of spread [m/s]
- `fire_wind_eff`: Effective wind after WAF and terrain corrections [m/s]
- `fire_wind_ref`: Reference height wind (6.1 m) [m/s]
- `fire_slopes`: Terrain slope components (dz/dx, dz/dy)
- `fire_fuel_mc`: Fuel moisture content (1hr, 10hr, 100hr)

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
