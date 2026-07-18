# Fire Model Tests

This directory contains organized test cases for the 2D fire model implementation in ERF, grouped by physical features and behavior being tested.

## Fire Model Overview

The fire model implements a 2D wildfire propagation solver using:
- **Rothermel Fire Spread Model**: Computes rate of fire spread based on fuel characteristics, weather, and topography
- **FARSITE Elliptical Propagation**: Models fire expansion as an ellipse with orientation based on wind and slope

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
cd Exec/CanonicalTests/Fire/Unit_Tests
python3 test_rothermel_unit.py              # Pure Python physics check (~1 sec)
python3 test_fire_ros_regression.py         # Reference value regression check (~1 sec)
```

### Run Integration Test (Requires Compiled ERF)
```bash
cd Exec/CanonicalTests/Fire/Core_Physics/ROS_Uniform_Grid
python3 ../../test_rothermel_ros.py --erf_exe ./erf --input_file inputs_fire_flat_uniform
```

## Test Directory Organization

### Unit_Tests/
Pure Python tests requiring no compiled ERF binary:
- `test_rothermel_unit.py`: Physics equation validation
- `test_fire_ros_regression.py`: Reference value checks

### Core_Physics/
Basic fire rate-of-spread calculations and fundamental fire physics:
- **ROS_Uniform_Grid/**: Rate of spread on flat uniform mesh
- **ROS_Basic_Calculation/**: Basic ROS computation on flat domain
- **ROS_Slope_Effects/**: ROS computation with terrain slope effects
- **Fuel_Moisture_Sensitivity/**: Fire spread at different fuel moisture levels (dry/wet conditions)
- **Wind_Speed_Variation/**: Fire spread at low and high wind speeds
- **Multiple_Fuel_Models/**: Different Anderson FBFM13 fuel types

### FARSITE_Propagation/
FARSITE elliptical fire expansion method:
- **Elliptical_Propagation/**: FARSITE elliptical fire front advancement

### Heat_Flux_Diagnostics/
Heat output and fuel depletion calculations:
- **Heat_Flux_and_Intensity/**: Heat flux, fuel depletion, fireline intensity, and flame length

### Fire_Atmosphere_Coupling/
Fire-to-atmosphere interactions and feedback:
- **Lagged_Coupling/**: Lagged fire-atmosphere feedback mode
- **Synchronous_Coupling/**: Synchronous fire-atmosphere feedback mode
- **Passive_Baseline/**: Passive mode (no fire-atmosphere feedback)

### Atmospheric_Boundary_Layer/
Atmospheric boundary layer physics with fire:
- **ABL_with_MRF/**: Atmospheric coupling with MRF turbulence model
- **Atmospheric_Stability/**: Tests in stable and unstable boundary layers

### Mesh_Refinement/
Grid refinement and mesh structure tests:
- **Vertical_Refinement/**: Stretched vertical grid fire spread

### Fire_Behavior/
Advanced fire behavior patterns:
- **Ignition_Patterns/**: Multiple ignition sources and fire interaction
- **Spotting/**: Albini (1983) stochastic ember spotting (Phase 8)

### Supporting_Files/
Documentation and utility scripts:
- `README.md`: This documentation
- `ANIMATION_GUIDE.md`: Guide for creating fire animations
- `plot_fire_animation.py`: Plotting utility

## Test Cases by Feature

### ROS Computation Tests
- **ROS_Uniform_Grid**: Minimal flat uniform mesh regression test with FM1 (Short Grass) fuel
- **ROS_Basic_Calculation**: Basic ROS computation on flat domain with GR1 fuel
- **ROS_Slope_Effects**: ROS with 30-degree slope and terrain wind corrections

### Fuel and Weather Variation Tests
- **Fuel_Moisture_Sensitivity/inputs_fire_moisture_dry**: Dry condition (3% moisture)
- **Fuel_Moisture_Sensitivity/inputs_fire_moisture_wet**: Wet condition (15% moisture)
- **Wind_Speed_Variation/inputs_fire_wind_low**: Low wind (1 m/s)
- **Wind_Speed_Variation/inputs_fire_wind_high**: High wind (10 m/s)
- **Multiple_Fuel_Models/inputs_fire_fuel_fm4_chaparral**: FM4 (Chaparral) fuel model

### Fire Propagation Tests
- **Elliptical_Propagation/**: FARSITE elliptical fire expansion method

### Heat Output Tests
- **Heat_Flux_and_Intensity/**: Heat flux computation, fuel depletion, Byram fireline intensity, Thomas flame length

### Atmospheric Coupling Tests
- **Lagged_Coupling/**: Fire heat flux injected one step lagged into atmospheric equations
- **Synchronous_Coupling/**: Fire advances using post-dycore wind with synchronized coupling
- **Passive_Baseline/**: Fire spreads under prescribed atmospheric conditions (no feedback)

### Boundary Layer Tests
- **ABL_with_MRF/**: Couples MRF boundary layer physics with fire model
- **Atmospheric_Stability/inputs_fire_stable_atmosphere**: Stable stratification conditions
- **Atmospheric_Stability/inputs_fire_unstable_atmosphere**: Unstable/convective conditions

### Mesh Tests
- **Vertical_Refinement/**: Stretched vertical grid with wind interpolation verification

### Fire Behavior Tests
- **Ignition_Patterns/inputs_fire_multiple_ignitions/**: Multiple ignition source interactions
- **Spotting/inputs_fire_phase8**: Albini (1983) stochastic ember spotting with lofting + trajectory + phi stamp

#### inputs_fire_vertical_refinement (Vertical_Refinement/)
Vertical grid stretching regression test:
- FM1 (Short Grass) fuel
- Stretched vertical grid (grid_stretching_ratio=1.08)
- Wind interpolation at z_ref=6.1 m with non-uniform dz
- **Expected ROS**: 0.1-1.0 m/s
- **Purpose**: Verify wind extraction and ROS computation with stretched grids

#### inputs_fire_heat_flux (Heat_Flux_and_Intensity/)
Heat flux computation test with diagnostics:
- FM1 (Short Grass) fuel
- 5 m/s nominal wind
- 8% fuel moisture
- Heat flux computation via Albini exponential burnout model
- Fuel load depletion via forward Euler
- Byram fireline intensity and Thomas flame length diagnostics
- CSV output with max heat flux tracking
- **Expected outputs**: fire_heat_flux > 0 in burned cells, fuel_load decreasing, flame length consistent with Byram
- **Purpose**: Regression test for heat flux, fuel depletion, fireline intensity, and flame length

#### inputs_fire_lagged_coupling (Lagged_Coupling/)
Lagged fire-to-atmosphere coupling test:
- FM1 (Short Grass) fuel
- 5 m/s nominal wind
- 8% fuel moisture
- Heat flux computation and fuel depletion
- Fire heat flux distributed through atmospheric column via exponential decay profile (WRF-Fire approach)
- Latent heat flux injection alongside sensible heat
- Lagged coupling: fire flux from step n enters dycore at step n+1 via RK source term
- Fire uses pre-dycore wind (vars_old)
- 30-minute simulation to allow coupling effects to develop
- **Expected outputs**: theta anomaly above burned area; qv increase if moisture is active; heat flux decays with altitude per exp(-z/alfg)
- **Purpose**: Regression test for WRF-Fire-parity fire-atmosphere coupling (lagged mode)

#### inputs_fire_synchronous_coupling (Synchronous_Coupling/)
Synchronous fire-to-atmosphere coupling test:
- FM1 (Short Grass) fuel
- 5 m/s nominal wind
- 8% fuel moisture
- Heat flux computation and fuel depletion (Phase 5)
- Fire heat flux distributed through atmospheric column via exponential decay profile
- Latent heat flux injection alongside sensible heat
- **Synchronous coupling**: fire advances using post-dycore wind (vars_new) that reflects atmospheric response to previous fire heating
- Fire heat flux from step n still enters dycore at step n+1 (same injection mechanism as lagged)
- 30-minute simulation
- **Expected outputs**: similar theta anomaly as phase6 but with different fire spread pattern reflecting post-dycore wind feedback
- **Purpose**: Regression test for synchronous fire-atmosphere coupling

#### inputs_fire_passive_coupling (Passive_Baseline/)
Passive fire-to-atmosphere coupling baseline test:
- Identical configuration to synchronous coupling except with passive coupling mode
- Fire spreads and heat flux is computed but NOT injected into atmosphere
- Fire uses pre-dycore wind (unchanged from initial state in absence of coupling)
- 30-minute simulation
- **Expected outputs**: fire spread with zero theta anomaly above burned area; no atmospheric response to fire
- **Purpose**: Baseline reference test for passive mode; demonstrates fire spread under prescribed atmospheric conditions

#### inputs_fire_phase8 (Fire_Behavior/Spotting/)
Phase 8 test case with Albini (1983) ember spotting:
- FM4 (Chaparral) fuel — produces high fireline intensity > 100 kW/m
- 8 m/s nominal wind to drive firebrand trajectories downwind
- 8% fuel moisture
- Passive coupling (no fire-atmosphere thermal feedback) to isolate spotting physics
- Albini lofting height H_z = 12.2 * I_B^(1/3) [m], Albini (1983) INT-309
- Forward-Euler trajectory: 20 sub-steps, terminal velocity 0.5 m/s
- Deterministic spotting with fixed random_seed = 42
- 15-minute simulation (900 s)
- **Expected outputs**: spot_fires_this_step > 0 after t ~ 60 s; fire_spot_active shows
  isolated patches downwind; max_spot_dist_m increases with fire intensity
- **Purpose**: Regression test for Albini spotting: lofting, trajectory, phi stamping, CSV output

#### test_albini_spotting.py (Unit_Tests/)
Unit tests for the Albini (1983) firebrand lofting and spotting distance model.
Verifies H_z = 12.2 × I_B^(1/3) against Albini (1983) INT-309 Table 1 reference values,
spotting distance scaling with wind speed and terminal velocity, and threshold behaviour.

Run: `python3 Exec/CanonicalTests/Fire/Unit_Tests/test_albini_spotting.py`

## Test Summary Table

| Test Category | Test Name | Fuel | Wind | Purpose |
|---------------|-----------|------|------|---------|
| Unit Tests | test_rothermel_unit.py | FM1, FM4 | Various | Physics equations validation (no ERF) |
| Unit Tests | test_fire_ros_regression.py | FM1, FM4 | Various | Reference value checks (no ERF) |
| ROS Uniform Grid | ROS_Flat_Uniform | FM1 | 5 m/s | Flat uniform mesh fire spread |
| ROS Basic | ROS_Basic_Calculation | GR1 | 5 m/s | Basic ROS computation |
| ROS with Terrain | ROS_Slope_Effects | GR1 | 3 m/s | Fire spread on sloped terrain |
| Fuel Moisture | Dry conditions | FM1 | 5 m/s | Fire spread at low moisture |
| Fuel Moisture | Wet conditions | FM1 | 5 m/s | Fire spread at high moisture |
| Wind Variation | Low wind | FM1 | 1 m/s | Fire spread at low wind speed |
| Wind Variation | High wind | FM1 | 10 m/s | Fire spread at high wind speed |
| Multiple Fuels | FM4 Chaparral | FM4 | 5 m/s | Dense, fast-burning fuel model |
| FARSITE Propagation | Elliptical propagation | GR1 | 5 m/s | FARSITE elliptical fire expansion |
| Heat Flux | Heat and diagnostics | FM1 | 5 m/s | Heat flux, fuel depletion, intensity |
| Lagged Coupling | Fire-atmosphere (lagged) | FM1 | 5 m/s | Fire-atmosphere coupling (lagged mode) |
| Synchronous Coupling | Fire-atmosphere (sync) | FM1 | 5 m/s | Fire-atmosphere coupling (synchronous mode) |
| Passive Baseline | Fire-atmosphere (passive) | FM1 | 5 m/s | Fire-atmosphere coupling (passive mode) |
| ABL with MRF | ABL with MRF physics | FM1 | Realistic | MRF boundary layer with fire |
| Stable Atmosphere | Stable conditions | FM1 | 5 m/s | Fire spread in stable stratification |
| Unstable Atmosphere | Unstable conditions | FM1 | 5 m/s | Fire spread in convective conditions |
| Mesh Refinement | Vertical refinement | FM1 | 5 m/s | Stretched vertical grid |
| Fire Behavior | Multiple ignitions | FM1 | 5 m/s | Multiple ignition source interaction |
| Fire Behavior | inputs_fire_phase8 | FM4 | 8 m/s | Albini spotting: lofting + trajectory + phi stamp |
| Unit Tests | test_albini_spotting.py | FM4 | Various | Albini lofting height and spotting distance (no ERF) |

## Output Files

ERF writes separate plot files for atmospheric and fire variables due to their different grid structures:

**Atmospheric 3D Variables** (atmospheric grid - 50m spacing):
- Directory: `plt3d_atm*` 
- Variables: theta, qv, u, v, w
- Configuration: `erf.plot_file`, `erf.plot_int`, `erf.plot_vars`

**Fire 2D Variables** (fire grid - 10m spacing, 5× refined):
- Directory: `plt_fire_XXXXX` (timestep numbered)
- Variables: fire physics on 10m fire cells
- Configuration: `erf.fire_plot_file`, `erf.fire_plot_int`

**Important**: Fire variables use separate dedicated fire plot output parameters because the fire grid (10m cells) is 5× finer than the atmospheric grid (50m cells).

## Fire Output Variables

All fire test cases output 2D fire-specific variables on the fire grid:

**Fire Output Configuration**:
- Parameter: `erf.fire_plot_file = plt_fire_` (output file prefix)
- Parameter: `erf.fire_plot_int = 100` (output interval in timesteps)
- Output directory: `plt_fire_XXXXX` (timestep-numbered directories)

**Fire variables on 10m fire grid**:

Heat flux tests additionally output:
- `fire_heat_flux`: Sensible heat flux [W/m²]
- `fire_fuel_load`: Remaining fuel load [kg/m²]
- `fire_fireline_intensity`: Byram fireline intensity [kW/m]
- `fire_flame_length`: Thomas flame length [m]

Fire-atmosphere coupling tests additionally output:
- `fire_heating_rate`: Heating rate induced by fire [K/s]
- `fire_moisture_flux`: Moisture flux from fire [kg/(m²·s)]

Output is written to separate `plt2d_fire*` directories (intervals) and `plt2d_fire_final*` (final timestep) on the fire grid (10m cells).

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


## Fire-to-Atmosphere Coupling Modes (Phases 6-7)

The fire model supports three fire-to-atmosphere coupling modes via `erf.fire.coupling_type`:

### Coupling Modes Overview

**Passive** (`coupling_type = "passive"`)
- Fire spread and heat flux are computed but the heat flux is NOT injected into the atmospheric equations.
- The atmosphere drives fire spread through wind; fire does not modify the atmosphere.
- Fire uses pre-dycore wind (vars_old) at each timestep.
- Use this mode to study fire spread under prescribed atmospheric conditions or as a baseline to assess coupling effects.

**Lagged** (`coupling_type = "lagged"`, default)
- Fire heat flux from step n is injected into the atmospheric RK source term at step n+1 using the WRF-Fire approach.
- The fire at step n+1 uses wind from the pre-dycore state (vars_old), which is the wind from the end of step n's dycore integration.
- This is the default mode and equivalent to WRF-Fire's standard explicit coupling (`fire_atm_feedback = 1`).
- The heat flux is re-injected inside the RK integrator (ERF_TI_slow_rhs_pre.H) at every RK stage so it survives the make_sources() reset at each stage (PR #84).

**Synchronous** (`coupling_type = "synchronous"`)
- Fire advances after `advance_dycore()`, using the post-dycore wind (vars_new) that reflects the atmospheric momentum response to fire heating from the previous step.
- Fire heat flux is still injected one step later via the same RK source term mechanism as lagged coupling (not immediately applied).
- The difference from lagged coupling is which wind the fire sees for spread rate calculation, not when the heat flux enters the dycore.
- Use this mode to capture the feedback between atmospheric heating and wind changes on fire spread.

### Coupling Mechanism

For both lagged and synchronous modes, the injection mechanism is identical:

1. **Heat Flux Computation**: Fire sensible heat flux Q [W/m²] is computed from fuel consumption rates and heat content (Phase 5).

2. **Vertical Distribution**: The surface heat flux is distributed through the atmospheric column using an exponential decay profile:

   ```
   hfx(z) = (Q_sfc / Cp_d) * exp(-z / alfg)
   ```

   where `z` is height above terrain and `alfg` is the e-folding scale height [m]. At `z = alfg`, the flux decays to 37% of surface value.

3. **Flux Divergence**: The vertical flux divergence is converted to a potential temperature tendency:

   ```
   d(RhoTheta)/dt = -rho * d(hfx)/dz
   ```

   This is added to the atmospheric source term inside the RK integrator at every RK stage.

4. **Latent Heat**: An equivalent latent heat flux is computed from fuel moisture and distributed the same way:

   ```
   qfx(z) = (Q_lat / Lv) * exp(-z / alfg)
   d(RhoQ1)/dt = -rho * d(qfx)/dz
   ```

5. **One-Step Lag**: Fire flux from timestep n is stored in `m_Q_atm_prev` and injected at timestep n+1. This lag is necessary because the fire is advanced AFTER the dycore has already run for the current step.

### Key Parameters (Phases 6-7)

- `erf.fire.coupling_type` (string, default `"lagged"`): Coupling mode — `"passive"`, `"lagged"`, or `"synchronous"`
- `erf.fire.heat_flux_alfg` (real, default 45.0): E-folding height scale for heat distribution [m]
- `erf.fire.fire_atm_feedback` (real, default 1.0): Multiplier on fire flux (0 = no feedback, 1 = full coupling)
- `erf.fire.inject_latent` (bool, default `true`): Inject latent heat alongside sensible heat

### Expected Physical Behavior

With coupling enabled (lagged or synchronous):

1. **Sensible Heat**: Produces a warming anomaly in the atmosphere above the fire, with maximum warming at the surface and exponential decay with height.

2. **Latent Heat**: Increases moisture in the atmospheric boundary layer if a moisture model is active.

3. **Dynamics**: The buoyancy from heating can induce local circulation and enhance updraft near the fire.

4. **Time Evolution**: The atmospheric response grows over the course of the fire spread, with maximum coupling when burning rate is highest.

5. **Wind Feedback** (synchronous only): The wind response to fire heating affects the fire spread rate in subsequent timesteps, creating a feedback loop.

### Example Usage

Run the lagged coupling test case (Phase 6):

```bash
./erf Exec/CanonicalTests/Fire/inputs_fire_phase6
```

Run the synchronous coupling test case (Phase 7):

```bash
./erf Exec/CanonicalTests/Fire/inputs_fire_phase7
```

Run the passive (no-feedback) baseline:

```bash
./erf Exec/CanonicalTests/Fire/inputs_fire_phase7_passive
```

### Customization Examples

To disable fire-to-atmosphere coupling but keep fire spread:

```
erf.fire.coupling_type = "passive"
```

To disable latent heat injection while maintaining sensible heat:

```
erf.fire.inject_latent = false
```

To reduce coupling strength:

```
erf.fire.fire_atm_feedback = 0.5
```

To increase the height scale over which heat is distributed:

```
erf.fire.heat_flux_alfg = 100.0
```

### References

- Mandel, J., et al. (2011). A wildland fire model with data assimilation. Mathematics and Computers in Simulation, 79(3), 584-606.
- WRF-SFIRE: https://wiki.ucar.edu/display/wrf/WRF%27s+Fire+and+Smoke+Module
