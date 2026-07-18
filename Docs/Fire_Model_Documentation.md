# 2D Fire Model Implementation - Phase 1, 2, and 3

## Overview

The 2D fire model in ERF implements a wildfire propagation solver using the Rothermel fire spread model combined with the FARSITE elliptical fire expansion algorithm. This documentation covers Phase 1 (framework), Phase 2 (Rothermel implementation), and Phase 3 (level-set propagation).

## Architecture

### Core Components

1. **Fire Layer Class** (`ERF_FireLayer.H`, `ERF_FireLayer.cpp`)
   - Main fire simulation container
   - Stores fire state on refined grid
   - Implements fire computation pipeline

2. **Integration with ERF**
   - Fire layer is instantiated in the main ERF class
   - Initialization occurs in `ERF::InitData_post()`
   - Fire state is advanced each atmospheric timestep via `FireLayer::advance()`

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
- Length-to-width ratio derived from wind speed (Anderson 1983)
- Richards (1990) coefficients for head, flank, and backing fire rates
- Level-set signed-distance representation of the fire front
- Two-pass propagation: GPU computation of spread vectors, then host collection and MPI gather

Reference: 
- Finney, M. A. (2004). FARSITE: Fire Area Simulator model development 
and evaluation. Res. Paper RMRS-RP-4 Revised, USDA Forest Service, Rocky Mountain 
Research Station.
- Richards, G.D. (1990). An elliptical growth model of forest fire fronts
and its numerical solution. Int. J. Numer. Meth. Eng. 30(6):1163-1179.

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

### Phase 2 Components

- [x] Rothermel model implementation
- [x] Rate-of-spread field computation
- [x] Wind extraction from MOST layer
- [x] Wind Adjustment Factor (WAF) application
- [x] Terrain slope and curvature computation
- [x] Terrain wind corrections (FARSITE)
- [x] Anderson FBFM13 fuel model database
- [x] Phase 2 regression test

### Phase 3 Components

- [x] Level-set field initialization (signed-distance circle)
- [x] FARSITE elliptical propagation kernel
- [x] Anderson L/W ratio computation (wind-dependent)
- [x] Richards coefficient conversion
- [x] Two-pass propagation: GPU spread vectors, host MPI gather
- [x] Single-cell and Gaussian stamping modes
- [x] CFL-limited fire subcycling
- [x] Phase 3 regression test

## Model State Variables

### Fire Grid MultiFabs (all C*nx × C*ny × 1)

| Name           | ncomp | Ghosts | Description                                  |
|----------------|-------|--------|----------------------------------------------|
| fire_phi       | 1     | 1      | Level-set: <0 burned, >0 unburned, ≈0 front  |
| fire_wind_ref  | 2     | 0      | u,v at z_ref (e.g., 6.1 m) from MOST        |
| fire_wind_eff  | 2     | 0      | After WAF + terrain corrections              |
| fire_slopes    | 2     | 0      | dz/dx, dz/dy (static, computed at init)      |
| fire_curvature | 1     | 0      | Terrain curvature (static, at init)          |
| fire_ros       | 1     | 0      | Rate of spread [m/s]                         |
| fire_fuel_load | 1     | 0      | Fuel load [kg/m²]                            |
| fire_fuel_mc   | 3     | 0      | M_1hr, M_10hr, M_100hr [fraction]           |
| fire_heat_flux | 1     | 0      | Heat flux [W/m²] (Phase 6, stub)             |
| fire_spread_vec| 2     | 0      | Scratch for FARSITE spread vectors (Phase 3) |

### FARSITE Parameters

Configuration parameters read from `erf.fire.farsite.*`:

- `phi_threshold`: Front detection level in phi (default 0.0)
- `use_anderson_lw`: Derive coefficients from Anderson L/W ratio (default 1)
- `coeff_a`: Richards head coefficient (default 0.5)
- `coeff_b`: Richards flank coefficient (default 0.25)
- `coeff_c`: Richards backing coefficient (default 0.1)
- `gaussian_sigma`: Phi stamping radius [m] (<0 = single-cell, =0 = auto, >0 = user)
- `cfl_fire`: Fire CFL number for subcycle timestep control (default 0.5)

### Level-Set Convention

The level-set field `phi` represents:
- phi < 0: burned region (inside fire)
- phi > 0: unburned fuel (outside fire)
- phi ≈ 0: fire front interface

## Dummy Functions

### Phase 1 Functions

### Define(int, SolverChoice&)
Initializes fire model parameters with default values. Reads user-provided parameters 
from input file.

### Init(int, MultiFab, Geometry, Real)
Performs initialization for a given AMR level. Stores grid dimensions.

### ComputeFireIntensity(int)
Dummy implementation of fire intensity calculation using empirical relationships.

## Input Parameters

Fire model parameters can be set in the input file using the `erf.fire.*` prefix:

### Basic Configuration

```
erf.fire.enable = true
erf.fire.grid_ratio = 5              # Fire grid refinement factor
erf.fire.fuel_model_id = 1           # Anderson FBFM13 model (1-13)
erf.fire.ignition_x = 1000.0         # Ignition center [m]
erf.fire.ignition_y = 1000.0         # Ignition center [m]
erf.fire.ignition_r = 20.0           # Ignition radius [m]
erf.fire.fire_debug = false          # Enable detailed debug output for fire calculations
```

### Moisture and Wind

```
erf.fire.moisture_1hr = 0.08
erf.fire.moisture_10hr = 0.08
erf.fire.moisture_100hr = 0.10
erf.fire.wind_ref_ht = 6.1           # Reference height for wind [m]
erf.fire.use_waf = true              # Apply Wind Adjustment Factor
erf.fire.waf_formula = "andrews"     # "andrews" or "behaviorplus"
```

### Terrain and Wind Corrections

```
erf.fire.use_terrain_wind = true     # Apply FARSITE terrain corrections
erf.fire.k_ridge = 1.5               # Ridge speed-up factor
erf.fire.k_shelter = 0.6             # Sheltered wind factor
erf.fire.k_valley = 0.8              # Valley channeling factor
erf.fire.k_deflect = 0.3             # Wind deflection factor
```

### Phase 3: FARSITE Propagation

```
erf.fire.farsite.phi_threshold = 0.0
erf.fire.farsite.use_anderson_lw = 1
erf.fire.farsite.gaussian_sigma = -1.0    # Single-cell stamping
erf.fire.farsite.cfl_fire = 0.5
```

## Testing

### Phase 2 Regression Test

Input: `inputs_fire_phase2`
- Flat domain, GR1 fuel, 5 m/s wind, 8% moisture
- Verifies Rothermel ROS computation and basic fire initialization
- Expected: max_ROS and mean_ROS values printed to stdout

### Phase 3 Regression Test

Input: `inputs_fire_phase3`
- Flat domain, GR1 fuel, 5 m/s wind, 8% moisture
- Verifies level-set propagation and FARSITE elliptical expansion
- Expected output checks:
  1. `phi_min < 0` at t > 0 — fire front has propagated
  2. `n_substeps` printed each step ≥ 1
  3. At t = 1800 s (30 min), burned ellipse with correct L/W ratio
  4. No NaN values in phi, ros, wind_eff

#### Expected Results at t = 1800 s

With GR1 fuel, 5 m/s wind, flat terrain:
- Head fire ROS: ~0.06–0.10 m/s (WAF-dependent)
- Wind speed in mph: ~5.6 (adjusted by WAF ≈ 0.4)
- Anderson L/W ratio: ~1.5–2.0
- Burned ellipse major axis (downwind): ~200–350 m
- Minor axis (cross-wind): ~100–200 m

## Debugging

The fire module provides debug output to trace calculation steps when enabled. This can be useful for verifying model behavior and investigating issues.

### Enabling Debug Output

Set `erf.fire.fire_debug = true` in the input file to enable step-by-step debug messages. When enabled, the fire module prints information about:
- Wind extraction from MOST boundary layer
- Wind Adjustment Factor application
- FARSITE terrain wind corrections
- Rate-of-spread field computation
- Level-set propagation and number of fire subcycles

Example input configuration:

```
erf.fire.enable = true
erf.fire.fire_debug = true
```

Debug messages are prefixed with `[FIRE DEBUG]` and printed to standard output alongside regular fire diagnostics (prefixed with `[FIRE]`).

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
