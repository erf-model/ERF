# Fire Model Tests

This directory contains test cases for the 2D fire model implementation in ERF.

## Fire Model Overview

The fire model implements a 2D wildfire propagation solver using:
- **Rothermel Fire Spread Model**: Computes rate of fire spread based on fuel characteristics, weather, and topography
- **FARSITE Elliptical Expansion**: Models fire expansion as an ellipse with orientation based on wind and slope

## Phase 1 Implementation

Phase 1 provides the basic framework with:
- Fire layer class structure
- Dummy implementations of core functions
- Basic input parameters
- Dummy regression tests

## Test Cases

### inputs_fire_dummy
Dummy test case that verifies:
- Fire model initialization
- Dummy function calls (Rothermel spread rate, elliptical expansion, intensity)
- Basic variable updates

## Running Tests

To run the fire model test:

```bash
# Compile ERF with fire model support
# Run the test
./erf inputs_fire_dummy
```

## Input Parameters

Key fire model parameters that can be set in input files:

- `fire.fuel_moisture_content`: Fuel moisture content (%)
- `fire.fuel_bed_depth`: Fuel bed depth (m)
- `fire.wind_speed`: Wind speed (m/s)
- `fire.slope`: Slope angle (degrees)
- `fire.ellipse_length_width_ratio`: Length-to-width ratio of fire ellipse

## Future Phases

- Phase 2-16: Full implementation of Rothermel model equations
- Integration with atmospheric variables (temperature, humidity, wind)
- Heat release and energy transfer modeling
- Lagrangian particle tracking for smoke and heat
- Multi-scale modeling with AMR
