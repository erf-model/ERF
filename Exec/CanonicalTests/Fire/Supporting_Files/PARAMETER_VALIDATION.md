# Fire Test Parameter Validation Report

## Overview
All parameters used in fire input files have been validated against the ERF codebase and confirmed to be available during parameter parsing.

## Fire Parameters (erf.fire.*)

All fire-specific parameters are defined in `Source/Fire/ERF_FireParams.H` and properly parsed:

### Core Parameters
✓ `enable` - Enable/disable fire module (bool)
✓ `fire_debug` - Debug output flag (bool)
✓ `grid_ratio` - Fire grid refinement factor (int)
✓ `fuel_model_id` - Anderson FBFM13 fuel model 1-13 (int)

### Fuel Moisture Parameters
✓ `moisture_1hr` - 1-hour fuel moisture content [fraction]
✓ `moisture_10hr` - 10-hour fuel moisture content [fraction]
✓ `moisture_100hr` - 100-hour fuel moisture content [fraction]
✓ `moisture_dynamic` - Enable dynamic moisture evolution (bool)
✓ `moisture_live` - Live fuel moisture [fraction]
✓ `precip_rate_mm_hr` - Precipitation rate [mm/hr]
✓ `use_dynamic_mext` - Compute moisture of extinction (bool)

### Ignition Parameters
✓ `ignition_x` - Ignition center X coordinate [m]
✓ `ignition_y` - Ignition center Y coordinate [m]
✓ `ignition_r` - Ignition radius [m]

### Wind and Terrain Parameters
✓ `wind_ref_ht` - MOST wind reference height [m]
✓ `use_waf` - Apply Wind Adjustment Factor (bool)
✓ `waf_formula` - WAF formula: "andrews" or "behaviorplus" (string)
✓ `use_wind_limit` - Apply MEWS wind speed cap (bool)
✓ `use_terrain_wind` - Apply terrain wind corrections (bool)
✓ `k_ridge` - Ridge speed-up factor
✓ `k_shelter` - Sheltered wind factor
✓ `k_valley` - Valley channeling factor
✓ `k_deflect` - Wind deflection factor

### Heat Flux Parameters (Phase 5+)
✓ `tau_residence_s` - Residence time [s]; 0 = derive from fuel
✓ `f_c_min` - Minimum combustion fraction
✓ `f_c_max` - Maximum combustion fraction
✓ `heat_flux_alfg` - E-folding height for heat distribution [m]

### FARSITE Propagation Parameters
✓ `farsite_phi_threshold` - Front detection level in phi
✓ `farsite_use_anderson_lw` - Derive coefficients from Anderson L/W (bool)
✓ `farsite_coeff_a` - Richards head coefficient
✓ `farsite_coeff_b` - Richards flank coefficient
✓ `farsite_coeff_c` - Richards backing coefficient
✓ `farsite_gaussian_sigma` - Phi stamping radius [m]
✓ `farsite_cfl_fire` - Fire CFL number

### Fire-Atmosphere Coupling Parameters (Phase 6-7)
✓ `coupling_type` - Coupling mode: "passive", "lagged", "synchronous"
✓ `fire_atm_feedback` - Fire flux multiplier (0-1)
✓ `inject_latent` - Inject latent heat (bool)

### Output Parameters
✓ `write_fire_stats_csv` - Write fire statistics CSV (bool)
✓ `fire_stats_csv_file` - CSV output filename (string)
✓ `terrain_file_name` - Fine-grid terrain file path (string)

## Fire Plot Parameters (erf.*)

✓ `erf.fire_plot_file` - Fire plot file name prefix
✓ `erf.fire_plot_int` - Fire plot output interval (timesteps)

**Note**: Fire variables are output to separate files from atmospheric variables because they exist on different grids (fire grid is 5× finer).

## General ERF Parameters (Validated)

### Time Control
✓ `erf.fixed_dt` - Fixed timestep [s]
✓ `stop_time` - Simulation end time [s]
✓ `max_step` - Maximum number of timesteps

### I/O and Checkpoints
✓ `erf.check_file` - Checkpoint file prefix
✓ `erf.check_int` - Checkpoint interval
✓ `erf.plot_file` - Atmospheric plot file prefix
✓ `erf.plot_int` - Atmospheric plot interval
✓ `erf.plot_vars` - Atmospheric variables to output
✓ `erf.sum_interval` - Statistics summary interval

### Physics Configuration
✓ `erf.init_type` - Initialization type ("uniform", etc.)
✓ `erf.theta_ref` - Reference potential temperature [K]
✓ `erf.molec_diff_type` - Molecular diffusion type
✓ `erf.les_type` - LES turbulence model
✓ `erf.rans_type` - RANS turbulence model
✓ `erf.use_gravity` - Include gravity (bool)
✓ `erf.use_coriolis` - Include Coriolis (bool)

### Numerics
✓ `erf.cfl` - CFL number
✓ `erf.substepping_cfl` - Substepping CFL

### Surface Layer (MOST)
✓ `erf.most.z0` - Roughness length [m]
✓ `erf.most.use_monin_obukhov` - Stability corrections (bool)
✓ `erf.most.zref` - Reference height

### ABL/MRF Parameters
✓ `erf.pbl_type` - PBL model type
✓ `erf.abl_driver_type` - ABL driver type
✓ `erf.abl_geo_wind` - Geostrophic wind configuration
✓ `erf.geostrophic_wind_x` - Geostrophic wind X [m/s]
✓ `erf.geostrophic_wind_y` - Geostrophic wind Y [m/s]
✓ `erf.latitude` - Latitude [degrees]
✓ `erf.rotational_time_period` - Coriolis period

### Atmospheric Stratification
✓ `erf.dtheta_ref` - Potential temperature gradient [K/m]

### Grid and Mesh
✓ `erf.grid_stretching_ratio` - Vertical grid stretching factor
✓ `erf.initial_dz` - Initial vertical cell size [m]

## Geometry and AMR (Standard)

✓ `geometry.prob_lo` - Domain lower bound
✓ `geometry.prob_hi` - Domain upper bound
✓ `geometry.is_periodic` - Periodicity (0/1 for each direction)
✓ `amr.n_cell` - Initial grid cells
✓ `amr.max_level` - Maximum AMR level
✓ `amr.max_grid_size` - Maximum grid size
✓ `amr.max_grid_size_z` - Maximum Z grid size

## Boundary Conditions (Standard)

✓ `xlo.type` - X-low boundary condition
✓ `xhi.type` - X-high boundary condition
✓ `ylo.type` - Y-low boundary condition
✓ `yhi.type` - Y-high boundary condition
✓ `zlo.type` - Z-low boundary condition (surface_layer for fire)
✓ `zhi.type` - Z-high boundary condition

## Test Files Summary

All 18 input files have been verified:
- ✓ 15 fire input files (inputs_fire_*)
- ✓ All parameters present and valid
- ✓ Fire plot output configured correctly
- ✓ Atmospheric plot output configured
- ✓ All boundary conditions set appropriately

### Files Checked
```
Core_Physics/
  ├── ROS_Basic_Calculation/inputs_fire_phase2
  ├── ROS_Uniform_Grid/inputs_fire_flat_uniform
  ├── ROS_Slope_Effects/inputs_fire_phase2_slope
  ├── Fuel_Moisture_Sensitivity/inputs_fire_moisture_dry
  ├── Fuel_Moisture_Sensitivity/inputs_fire_moisture_wet
  ├── Wind_Speed_Variation/inputs_fire_wind_low
  ├── Wind_Speed_Variation/inputs_fire_wind_high
  └── Multiple_Fuel_Models/inputs_fire_fuel_fm4_chaparral

FARSITE_Propagation/
  └── Elliptical_Propagation/inputs_fire_phase3

Heat_Flux_Diagnostics/
  └── Heat_Flux_and_Intensity/inputs_fire_phase5

Fire_Atmosphere_Coupling/
  ├── Lagged_Coupling/inputs_fire_phase6
  ├── Synchronous_Coupling/inputs_fire_phase7
  └── Passive_Baseline/inputs_fire_phase7_passive

Atmospheric_Boundary_Layer/
  ├── ABL_with_MRF/inputs_fire_abl_mrf_unstable
  └── Atmospheric_Stability/
      ├── inputs_fire_stable_atmosphere
      └── inputs_fire_unstable_atmosphere

Mesh_Refinement/
  └── Vertical_Refinement/inputs_fire_vertical_refinement

Fire_Behavior/
  └── Ignition_Patterns/inputs_fire_multiple_ignitions
```

## Validation Status: ✓ PASSED

All parameters in fire test input files are valid and available during parameter parsing.
No invalid, undefined, or deprecated parameters detected.

**Last Updated**: 2026-07-07
**ERF Version**: Current main branch
