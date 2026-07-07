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

#### inputs_fire_phase6
Phase 6 test case with lagged fire-to-atmosphere coupling:
- FM1 (Short Grass) fuel
- 5 m/s nominal wind
- 8% fuel moisture
- Heat flux computation and fuel depletion (Phase 5)
- Fire heat flux distributed through atmospheric column via exponential decay profile (WRF-Fire approach)
- Latent heat flux injection alongside sensible heat
- Lagged coupling: fire flux from step n enters dycore at step n+1 via RK source term
- Fire uses pre-dycore wind (vars_old)
- 30-minute simulation to allow coupling effects to develop
- **Expected outputs**: theta anomaly above burned area; qv increase if moisture is active; heat flux decays with altitude per exp(-z/alfg)
- **Purpose**: Regression test for WRF-Fire-parity fire-atmosphere coupling (lagged mode)

#### inputs_fire_phase7
Phase 7 test case with synchronous fire-to-atmosphere coupling:
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

#### inputs_fire_phase7_passive
Phase 7 baseline test with passive fire-to-atmosphere coupling:
- Identical configuration to inputs_fire_phase7 except with passive coupling mode
- Fire spreads and heat flux is computed but NOT injected into atmosphere
- Fire uses pre-dycore wind (unchanged from initial state in absence of coupling)
- 30-minute simulation
- **Expected outputs**: fire spread with zero theta anomaly above burned area; no atmospheric response to fire
- **Purpose**: Baseline reference test for passive mode; demonstrates fire spread under prescribed atmospheric conditions

## Test Summary Table

| Test Name | Type | Fuel | Wind | Purpose | Command |
|-----------|------|------|------|---------|---------|
| test_rothermel_unit.py | Unit | FM1, FM4 | Various | Physics equations validation (no ERF) | `python3 test_rothermel_unit.py` |
| test_fire_ros_regression.py | Regression | FM1, FM4 | Various | Reference value checks (no ERF) | `python3 test_fire_ros_regression.py` |
| test_rothermel_ros.py | Integration | FM1 | Uniform | ERF ROS magnitude check | `python3 test_rothermel_ros.py --erf_exe ./erf --input_file inputs_fire_flat_uniform` |
| inputs_fire_flat_uniform | Regression | FM1 | 5 m/s | Flat uniform mesh fire spread | `./erf inputs_fire_flat_uniform` |
| inputs_fire_vertical_refinement | Regression | FM1 | 5 m/s | Stretched vertical grid fire spread | `./erf inputs_fire_vertical_refinement` |
| inputs_fire_phase5 | Regression | FM1 | 5 m/s | Heat flux and diagnostics | `./erf inputs_fire_phase5` |
| inputs_fire_phase6 | Regression | FM1 | 5 m/s | Fire-atmosphere coupling (lagged mode) | `./erf inputs_fire_phase6` |
| inputs_fire_phase7 | Regression | FM1 | 5 m/s | Fire-atmosphere coupling (synchronous mode) | `./erf inputs_fire_phase7` |
| inputs_fire_phase7_passive | Regression | FM1 | 5 m/s | Fire-atmosphere coupling (passive mode) | `./erf inputs_fire_phase7_passive` |
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
