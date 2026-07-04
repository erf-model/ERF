# 2D Fire Model Implementation - Phase 1

## Overview

The 2D fire model in ERF implements a wildfire propagation solver using the Rothermel fire spread model combined with the FARSITE elliptical fire expansion algorithm. This documentation covers Phase 1, which provides the basic framework and dummy implementations.

## Architecture

### Core Components

1. **Fire Model Class** (`ERF_Fire.H`, `ERF_Fire.cpp`)
   - Main fire layer class
   - Stores fire model parameters and state variables
   - Implements dummy functions for core calculations

2. **Integration with ERF**
   - Fire model is instantiated in the main ERF class
   - Initialization occurs during ERF construction
   - Fire state can be advanced during main time loop

### Physical Models

#### Rothermel Fire Spread Model

The Rothermel model computes fire spread rate based on:
- Fuel characteristics (moisture, bed depth, particle density, energy content, load)
- Environmental conditions (wind speed, slope angle)
- Rate of spread calculations for head, flank, and back fires

Reference: Rothermel, R. C. (1972). A mathematical model for predicting fire spread 
in wildland fuels. Res. Paper INT-115, USDA Forest Service, Intermountain Forest 
and Range Experiment Station.

#### FARSITE Elliptical Expansion

The FARSITE algorithm models fire expansion as an ellipse with:
- Length-to-width ratio
- Eccentricity based on wind and slope
- Head, flank, and back fire rates defining ellipse axes

Reference: Finney, M. A. (2004). FARSITE: Fire Area Simulator model development 
and evaluation. Res. Paper RMRS-RP-4 Revised, USDA Forest Service, Rocky Mountain 
Research Station.

## Implementation Status

### Phase 1 Components

- [x] Fire layer class structure
- [x] Dummy Rothermel model function
- [x] Dummy FARSITE elliptical expansion function
- [x] Dummy fire intensity calculation
- [x] Integration with main ERF class
- [x] Input file support
- [x] Dummy regression test
- [x] Basic documentation

### Dummy Variables

#### Rothermel Model Parameters
- `fuel_moisture_content`: Fuel moisture content (%)
- `fuel_bed_depth`: Fuel bed depth (m)
- `fuel_particle_density`: Fuel particle density (kg/m³)
- `fuel_energy_content`: Fuel energy content (J/kg)
- `fuel_load`: Fuel load (kg/m²)
- `wind_speed`: Wind speed (m/s)
- `slope`: Slope angle (degrees)

#### FARSITE Elliptical Expansion Parameters
- `ellipse_length_width_ratio`: Length-to-width ratio of fire ellipse
- `ellipse_eccentricity`: Eccentricity of fire ellipse
- `ellipse_major_axis`: Major axis length (m)
- `ellipse_minor_axis`: Minor axis length (m)

#### Fire State Variables
- `head_fire_rate_of_spread`: Head fire rate of spread (m/s)
- `flank_fire_rate_of_spread`: Flank fire rate of spread (m/s)
- `back_fire_rate_of_spread`: Back fire rate of spread (m/s)
- `fire_line_intensity`: Fire line intensity (W/m)
- `flame_length`: Flame length (m)

## Dummy Functions

### Define(int, SolverChoice&)
Initializes fire model parameters with default values. Reads user-provided parameters 
from input file.

### Init(int, MultiFab, Geometry, Real)
Performs initialization for a given AMR level. Stores grid dimensions.

### Advance(int, Real, Real, MultiFab, Geometry)
Advances fire simulation by one time step by calling Rothermel, elliptical expansion, 
and intensity calculations.

### ComputeRothermellSpreadRate(int, Geometry)
Dummy implementation of Rothermel fire spread rate calculation.

### ComputeEllipticalExpansion(int, Geometry)
Dummy implementation of FARSITE elliptical fire expansion calculation.

### ComputeFireIntensity(int)
Dummy implementation of fire intensity calculation using empirical relationships.

## Input Parameters

Fire model parameters can be set in the input file using the `fire.*` prefix:

```
# Fire model parameters
fire.fuel_moisture_content = 0.08
fire.fuel_bed_depth = 0.3
fire.wind_speed = 5.0
fire.slope = 0.0
fire.ellipse_length_width_ratio = 1.5
```

## Testing

### Running Dummy Test

```bash
python3 test_fire_dummy.py --erf_exe ./erf --input_file inputs_fire_dummy
```

### Test Verification

The dummy test verifies:
1. Model initialization
2. Dummy function calls during time stepping
3. Successful completion without errors

## Future Phases

Phase 2-16 will implement:
1. Full Rothermel model equations
2. Wind and slope factor calculations
3. Fuel model parameterizations
4. Atmospheric coupling
5. Heat release modeling
6. Smoke and particle tracking
7. Multi-scale AMR refinement
8. Performance optimization
9. Visualization and analysis tools
10. Integration with atmospheric chemistry
11. Fire behavior prediction validation
12. Advanced fuel characterization
13. Real-world scenario implementation
14. Model coupling with radiation
15. Final validation and benchmarking

## References

- Rothermel, R. C. (1972). A mathematical model for predicting fire spread in wildland 
  fuels. Res. Paper INT-115, USDA Forest Service.
- Finney, M. A. (2004). FARSITE: Fire Area Simulator model development and evaluation. 
  Res. Paper RMRS-RP-4 Revised, USDA Forest Service.
- Andrews, P. L. (2018). Current status and future needs of the BehavePlus Fire 
  Modeling System. International Journal of Wildland Fire, 27(9), 558-566.
