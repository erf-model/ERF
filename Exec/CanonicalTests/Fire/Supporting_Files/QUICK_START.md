# Fire Test Suite - Quick Reference Guide

## Running Tests Quickly

### Unit Tests (No ERF Executable Required - ~2 seconds)
```bash
cd Unit_Tests
python3 test_rothermel_unit.py              # Rothermel physics validation
python3 test_fire_ros_regression.py         # Reference value checks
python3 test_farsite_ellipse.py             # FARSITE ellipse propagation
```

### Integration Tests (Requires Compiled ERF)

#### Core ROS Tests
```bash
cd Core_Physics/ROS_Basic_Calculation
../../Core_Physics/ROS_Uniform_Grid/test_rothermel_ros.py --erf_exe ./erf --input_file inputs_fire_phase2

cd Core_Physics/ROS_Uniform_Grid
./erf inputs_fire_flat_uniform

cd Core_Physics/ROS_Slope_Effects
./erf inputs_fire_phase2_slope
```

#### Fuel/Weather Sensitivity Tests
```bash
cd Core_Physics/Fuel_Moisture_Sensitivity
./erf inputs_fire_moisture_dry           # Dry condition (3% moisture)
./erf inputs_fire_moisture_wet           # Wet condition (15% moisture)

cd Core_Physics/Wind_Speed_Variation
./erf inputs_fire_wind_low               # Low wind (1 m/s)
./erf inputs_fire_wind_high              # High wind (10 m/s)

cd Core_Physics/Multiple_Fuel_Models
./erf inputs_fire_fuel_fm4_chaparral     # FM4 Chaparral fuel model
```

#### Fire Propagation Tests
```bash
cd FARSITE_Propagation/Elliptical_Propagation
./erf inputs_fire_phase3                 # FARSITE elliptical expansion
```

#### Heat Flux Tests
```bash
cd Heat_Flux_Diagnostics/Heat_Flux_and_Intensity
./erf inputs_fire_phase5                 # Heat flux, fuel depletion, intensity
```

#### Fire-Atmosphere Coupling Tests
```bash
cd Fire_Atmosphere_Coupling/Lagged_Coupling
./erf inputs_fire_phase6                 # Lagged coupling

cd Fire_Atmosphere_Coupling/Synchronous_Coupling
./erf inputs_fire_phase7                 # Synchronous coupling

cd Fire_Atmosphere_Coupling/Passive_Baseline
./erf inputs_fire_phase7_passive         # Passive (no feedback) baseline
```

#### Atmospheric Boundary Layer Tests
```bash
cd Atmospheric_Boundary_Layer/ABL_with_MRF
./erf inputs_fire_abl_mrf_unstable       # MRF boundary layer with fire

cd Atmospheric_Boundary_Layer/Atmospheric_Stability
./erf inputs_fire_stable_atmosphere      # Stable stratification
./erf inputs_fire_unstable_atmosphere    # Unstable/convective
```

#### Mesh Refinement Tests
```bash
cd Mesh_Refinement/Vertical_Refinement
./erf inputs_fire_vertical_refinement    # Stretched vertical grid
```

#### Fire Behavior Tests
```bash
cd Fire_Behavior/Ignition_Patterns
./erf inputs_fire_multiple_ignitions     # Multiple ignition sources
```

## Test Organization by Feature

| Feature | Category | Tests |
|---------|----------|-------|
| **Rate of Spread** | Core_Physics | ROS_Basic_Calculation, ROS_Uniform_Grid, ROS_Slope_Effects |
| **Fuel Properties** | Core_Physics | Fuel_Moisture_Sensitivity (dry/wet) |
| **Wind Effects** | Core_Physics | Wind_Speed_Variation (low/high), ROS_Slope_Effects |
| **Fuel Models** | Core_Physics | Multiple_Fuel_Models (FM1, FM4) |
| **Propagation** | FARSITE_Propagation | Elliptical_Propagation |
| **Heat Output** | Heat_Flux_Diagnostics | Heat_Flux_and_Intensity |
| **Coupling Modes** | Fire_Atmosphere_Coupling | Lagged_Coupling, Synchronous_Coupling, Passive_Baseline |
| **ABL Physics** | Atmospheric_Boundary_Layer | ABL_with_MRF, Atmospheric_Stability (stable/unstable) |
| **Mesh Handling** | Mesh_Refinement | Vertical_Refinement |
| **Fire Behavior** | Fire_Behavior | Ignition_Patterns |

## Output Variables by Test Type

### Separate Grid Structures
Fire variables are output to **separate 2D plot files** (plt2d_fire*) from atmospheric variables (plt3d_atm*) because they live on different grids:
- **Atmospheric grid**: 50m spacing (40×40 cells in horizontal)
- **Fire grid**: 10m spacing (200×200 cells in horizontal, 5× refined)

### All Standard Tests Output

**At intervals** (every 100 timesteps):
- Fire directory: `plt2d_fire_XXXXX` (on 10m fire grid)
- Atmospheric directory: `plt3d_atm_XXXXX` (on 50m atmospheric grid)

**At final timestep**:
- Fire directory: `plt2d_fire_final_XXXXX` (on 10m fire grid)

**Fire variables** (10m grid):
- `fire_phi` - Fire front status (burned/unburned)
- `fire_ros` - Rate of spread [m/s]
- `fire_wind_eff` - Effective wind [m/s]
- `fire_wind_ref` - Reference height wind [m/s]
- `fire_slopes` - Terrain slope [dz/dx, dz/dy]
- `fire_fuel_mc` - Fuel moisture content [fraction]
- `fire_tsurf` - Surface temperature [K]
- `fire_fuel_burned` - Cumulative fuel burned [kg/m²]
- `fire_ignition_mask` - Ignition state (0=unburned, 1=burning, 2=burned)

**Atmospheric variables** (50m grid):
- theta, qv, u, v, w

### Heat Flux Tests Additionally Output
- `fire_heat_flux` - Sensible heat flux [W/m²]
- `fire_fuel_load` - Remaining fuel load [kg/m²]
- `fire_fireline_intensity` - Byram intensity [kW/m]
- `fire_flame_length` - Thomas flame length [m]

### Coupling Tests Additionally Output
- `fire_heating_rate` - Atmospheric heating rate [K/s]
- `fire_moisture_flux` - Moisture flux [kg/(m²·s)]

## Expected Output Patterns

### Dry Conditions Test
- Higher ROS than baseline (3% vs 8% moisture)
- Faster fire spread expected
- Lower ignition delay

### Wet Conditions Test
- Lower ROS than baseline (15% vs 8% moisture)
- Slower fire spread expected
- Higher ignition delay

### Low Wind Test
- Lower ROS at all points
- More circular fire shape
- Lower wind factor contribution

### High Wind Test
- Higher ROS, especially head fire
- More elongated elliptical shape
- Strong directional effect

### FM4 Chaparral Test
- Significantly higher ROS than GR1 grass
- Denser fuel with higher heat content
- More intense fire behavior

### Stable Atmosphere Test
- Reduced vertical mixing
- Wind speeds closer to ground values
- Different fire-atmosphere interaction pattern

### Unstable Atmosphere Test
- Enhanced vertical mixing
- Higher effective wind speeds aloft
- Stronger updraft response to fire heating

### FARSITE Propagation Test
- Elliptical fire front
- Head fire advances faster than flanks
- Back fire spreads slowly into wind

### Heat Flux Tests
- Non-zero heat flux in burned areas
- Fuel load decreases over time
- Flame length and intensity consistent with ROS

### Coupling Tests
- Lagged: theta anomaly appears one step after burning
- Synchronous: theta anomaly affects next fire advance
- Passive: no atmospheric response to fire heating

## Common Parameters Across Tests

| Parameter | Value | Notes |
|-----------|-------|-------|
| Domain size | 2000×2000 m | Standard fire domain |
| Atm grid | 40×40×50 cells | ~50 m horizontal spacing |
| Fire grid ratio | 5 | 10 m fire cells |
| Time step | 1.0 s | Fixed timestep |
| Simulation time | 1800 s (30 min) | Standard 30-minute run |
| Reference height | 6.1 m | MOST wind extraction |
| Fuel moisture baseline | 8% (1hr, 10hr) | Standard conditions |
| Ignition | Center at (1000, 1000) | Domain center |
| Ignition radius | 20 m | Initial fire size |

## Creating Custom Tests

To create a new fire test:

1. **Copy template input file** from similar test
2. **Modify parameters**:
   - `erf.geostrophic_wind_x/y` - Wind speed/direction
   - `erf.fire.moisture_*hr` - Fuel moisture
   - `erf.fire.fuel_model_id` - Fuel type (1=GR1, 4=FM4, etc.)
   - `erf.dtheta_ref` - Atmospheric stability
   - Domain geometry if needed
3. **Update output variables** in `erf.plot2d_vars_1`
4. **Run**: `./erf inputs_your_test`
5. **Analyze**: Plot files saved to `plt2d_fire*`

## File Organization Reference

```
Fire/
├── Unit_Tests/
│   ├── test_rothermel_unit.py
│   ├── test_fire_ros_regression.py
│   └── test_farsite_ellipse.py
├── Core_Physics/
│   ├── ROS_Basic_Calculation/
│   ├── ROS_Uniform_Grid/
│   ├── ROS_Slope_Effects/
│   ├── Fuel_Moisture_Sensitivity/
│   ├── Wind_Speed_Variation/
│   └── Multiple_Fuel_Models/
├── FARSITE_Propagation/
│   └── Elliptical_Propagation/
├── Heat_Flux_Diagnostics/
│   └── Heat_Flux_and_Intensity/
├── Fire_Atmosphere_Coupling/
│   ├── Lagged_Coupling/
│   ├── Synchronous_Coupling/
│   └── Passive_Baseline/
├── Atmospheric_Boundary_Layer/
│   ├── ABL_with_MRF/
│   └── Atmospheric_Stability/
├── Mesh_Refinement/
│   └── Vertical_Refinement/
├── Fire_Behavior/
│   └── Ignition_Patterns/
└── Supporting_Files/
    ├── README.md
    ├── QUICK_START.md (this file)
    ├── ANIMATION_GUIDE.md
    └── plot_fire_animation.py
```

## Troubleshooting

**Fire doesn't spread**
- Check fuel moisture (too high)
- Check wind speed (too low)
- Verify ignition radius and location
- Check fuel model ID (1-13 valid)

**Very slow spreading**
- Increase wind speed
- Decrease moisture content
- Check mesh resolution (may need finer grid)

**Fire spreads unrealistically fast**
- Decrease wind speed
- Increase moisture content
- Check fuel model (FM4 naturally faster than FM1)

**Errors in output**
- Verify fire grid ratio (typically 5)
- Check domain periodicity (must be 0 0 0)
- Ensure max_grid_size and max_grid_size_z are set
- Verify boundary conditions for fire (outflow on sides)
