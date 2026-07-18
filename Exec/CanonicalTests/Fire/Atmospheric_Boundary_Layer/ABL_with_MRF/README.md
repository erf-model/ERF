# ABL-Enhanced Fire Configurations

## Overview

This directory contains fire simulation configurations combined with physically realistic Atmospheric Boundary Layer (ABL) scenarios. These configurations are designed to test fire behavior under different atmospheric conditions that can realistically occur before and during wildfire events.

## Key Features

- **Periodic Boundary Conditions**: All configurations use periodic boundaries in x and y directions (`geometry.is_periodic = 1 1 0`), with surface layer at bottom and slip wall at top
- **Physically Realistic Scenarios**: Only ABL conditions that are conducive to pre-fire atmospheric states are included
- **Matched Wind Parameters**: `erf.abl_geo_wind` values are extracted directly from sounding files to ensure atmospheric consistency
- **Integrated Fire Physics**: Fire module is enabled with appropriate fuel models and moisture parameters

## Available Scenarios

### 1. **Neutral** (10.0 m/s)
- **File**: `inputs_fire_abl_neutral`
- **Sounding**: `input_sounding_MRF_neutral`
- **Description**: Shear-driven turbulence with uniform potential temperature profile
- **Characteristics**:
  - Zero surface heat flux (adiabatic)
  - Turbulence driven by wind shear, not buoyancy
  - PBL height governed by wind shear relation: h ~ u*/f
- **Fire Relevance**: Typical pre-fire state with moderate wind speeds and stable buoyancy
- **Physics**:
  - `erf.most.surf_heating_rate = 0.0` (adiabatic surface)
  - `zhi.theta_grad = 0.0` (adiabatic upper boundary)

### 2. **Unstable** (5.0 m/s)
- **File**: `inputs_fire_abl_unstable`
- **Sounding**: `input_sounding_MRF_unstable`
- **Description**: Daytime convective conditions with buoyancy-driven turbulence
- **Characteristics**:
  - Positive sensible heat flux creating unstable conditions
  - Convective boundary layer with enhanced vertical mixing
  - Lower shear production, buoyancy-driven turbulence
- **Fire Relevance**: Common afternoon conditions when fires are most active; unstable profiles enhance fire-induced convection
- **Physics**:
  - Potential temperature: 303.0 K (warmer than neutral)
  - Enhanced vertical mixing promotes updrafts

### 3. **Diurnal** (8.0 m/s)
- **File**: `inputs_fire_abl_diurnal`
- **Sounding**: `input_sounding_MRF_diurnal`
- **Description**: Diurnal heating cycle typical of fire-prone days
- **Characteristics**:
  - Intermediate wind speeds with transitional stability
  - Heating-driven PBL evolution
  - Moderate wind for fire propagation
- **Fire Relevance**: Represents typical afternoon fire conditions with moderate winds and heating
- **Physics**:
  - Cooler surface at start (295.0 K) evolving through diurnal cycle
  - Transitional stability conditions

### 4. **Extreme Heating** (4.0 m/s)
- **File**: `inputs_fire_abl_extreme_heating`
- **Sounding**: `input_sounding_MRF_extreme_heating`
- **Description**: High sensible heat flux preceding intense fire activity
- **Characteristics**:
  - Extreme atmospheric instability
  - Very strong buoyancy-driven convection
  - Favorable conditions for strong fire-induced updrafts
- **Fire Relevance**: Pre-fire atmospheric conditions with maximum instability; generates strong updrafts that enhance fire intensity
- **Physics**:
  - Aggressive heating profile
  - Maximum buoyancy production
  - Strong convective overturning

### 5. **High Shear** (35.0 m/s)
- **File**: `inputs_fire_abl_high_shear`
- **Sounding**: `input_sounding_MRF_high_shear`
- **Description**: Strong wind shear conditions promoting rapid fire spread
- **Characteristics**:
  - Very high geostrophic wind speed (35 m/s)
  - Shear-dominated turbulence
  - Minimal buoyancy effects
- **Fire Relevance**: Extreme wind scenarios that dramatically enhance fire spread; represents Santa Ana, Föhn, or other high-wind fire events
- **Physics**:
  - Neutral stability maintained at high wind speeds
  - Shear production dominates buoyancy
  - Enhanced horizontal momentum transport

### 6. **Weak Convection Transition** (8.0 m/s)
- **File**: `inputs_fire_abl_weak_convection_transition`
- **Sounding**: `input_sounding_MRF_weak_convection_transition`
- **Description**: Transitional conditions between stable and unstable regimes
- **Characteristics**:
  - Weak convection with modest instability
  - Combined shear and buoyancy effects
  - Transitional PBL characteristics
- **Fire Relevance**: Morning-to-afternoon transition when fires transition from wind-driven to convection-driven
- **Physics**:
  - Mixed-mode turbulence (shear + buoyancy)
  - Transitional Richardson number
  - Evolving vertical structure

## Configuration Details

### Geometry & Domain
- **Domain extent**: 8000 m × 8000 m × 2048 m (x, y, z)
- **Resolution**: 8 × 8 × 64 cells (coarse for quick testing)
- **Boundary conditions**:
  - Periodic in x and y: `geometry.is_periodic = 1 1 0`
  - Surface layer at bottom: `zlo.type = "surface_layer"`
  - Slip wall at top: `zhi.type = "SlipWall"`
  - **No xlo, xhi, ylo, yhi type specifications** (implicit periodic)

### Surface Layer (MOST)
- **Roughness length**: z₀ = 0.1 m
- **Reference height**: 48.0 m
- **Reference temperature**: 300.0 K

### Fire Configuration
- **Fire model**: Enabled (`erf.fire.enable = true`)
- **Fuel model**: FM4 (Chaparral) - suitable for diverse scenarios
- **Grid refinement**: 5x finer in horizontal
- **Fuel moisture** (NFMD standard):
  - 1-hour: 12%
  - 10-hour: 12%
  - 100-hour: 15%
  - Live: 90%
- **Wind reference height**: 10.0 m
- **Coupling**: Lagged coupling between atmosphere and fire

### Time Control
- **Time step**: 0.1 s (fixed)
- **Duration**: 100 steps (~10 seconds simulation)
- **Output interval**: Every 15 time steps

## Wind Validation Matrix

| Scenario | Sounding File | erf.abl_geo_wind | Justification |
|----------|---------------|------------------|---------------|
| Neutral | input_sounding_MRF_neutral | 10.0 0.0 0.0 | Typical pre-fire state |
| Unstable | input_sounding_MRF_unstable | 5.0 0.0 0.0 | Daytime heating |
| Diurnal | input_sounding_MRF_diurnal | 8.0 0.0 0.0 | Afternoon fire conditions |
| Extreme Heating | input_sounding_MRF_extreme_heating | 4.0 0.0 0.0 | Pre-intense fire conditions |
| High Shear | input_sounding_MRF_high_shear | 35.0 0.0 0.0 | Santa Ana / Föhn events |
| Weak Conv. Trans. | input_sounding_MRF_weak_convection_transition | 8.0 0.0 0.0 | Morning-afternoon transition |

## How to Use

### Run a Single Scenario
```bash
cd Exec/CanonicalTests/Fire/Atmospheric_Boundary_Layer/ABL_with_MRF
mpirun -np 4 ../../erf_exec inputs_fire_abl_neutral
```

### Run All Scenarios
```bash
for scenario in neutral unstable diurnal extreme_heating high_shear weak_convection_transition; do
  echo "Running $scenario..."
  mpirun -np 4 ../../erf_exec inputs_fire_abl_$scenario > log_$scenario.log
done
```

### Output Files
For each scenario `<scenario>`, the following files are generated:
- `plt_fire_abl_<scenario>_*` - Fire visualization files
- `fire_stats_abl_<scenario>.csv` - Fire statistics and diagnostics
- Standard ERF output files (checkpoint, plot files, etc.)

## Physical Justification for Scenario Selection

### Why These Scenarios?
All selected scenarios represent **pre-fire or fire-active atmospheric conditions**:

1. **Neutral** - Baseline shear-driven state before heating develops
2. **Unstable** - Daytime convection when most fires occur
3. **Diurnal** - Evolution through typical fire-prone day
4. **Extreme Heating** - Maximum atmospheric instability (pre-intense fire)
5. **High Shear** - Wind-driven extreme fire events
6. **Weak Conv. Trans.** - Transitional morning-to-afternoon regime

### Excluded Scenarios
The following ABL scenarios were **not included** as they are unphysical or non-fire-conducive:
- **Stable**: Suppresses buoyancy-driven fire dynamics
- **Arctic**: Unphysical for wildfire context
- **Calm Conditions**: Insufficient wind for fire propagation
- **Cloud Topped**: Clouds suppress solar heating needed for fire initiation
- **Fine dt Stable**: Stable conditions counter to fire dynamics
- **Marine**: Stable marine layer suppresses fire conditions
- **Rapid Cooling**: Cooling opposes fire conditions
- **Saturated Layer**: High moisture content suppresses fire

## Fire Physics Considerations

### Fire-Atmosphere Coupling
- **Coupling Type**: Lagged coupling for computational efficiency
- **Heat Flux Scaling**: `heat_flux_alfg = 45.0` (moderate to strong)
- **Feedback**: `fire_atm_feedback = 1.0` (full coupling)
- **Latent Heat**: Injected into atmosphere via `inject_latent = true`

### Wind Representation
All scenarios represent **geostrophic wind** at the domain scale. Fire-resolving grid sees wind modulated by:
- Surface layer effects (MOST parameterization)
- Turbulence-induced perturbations
- Fire-induced modifications (when coupled)

### Model Validation
These configurations are suitable for:
- Testing fire behavior across stability regimes
- Validating fire-atmosphere coupling algorithms
- Studying wind-driven vs. buoyancy-driven fire regimes
- Benchmarking spread rate predictions
- Comparing fire models under realistic conditions

## References

### Atmospheric Dynamics
- **Neutral ABL**: Sorbjan, Z. (1989), Structure of the Atmospheric Boundary Layer
- **Convective ABL**: Stull, R. B. (1988), An Introduction to Boundary Layer Meteorology
- **Wind Shear**: Panofsky, H. A., and Dutton, J. A. (1984), Atmospheric Turbulence

### Fire-Atmosphere Coupling
- **Heat Flux Models**: Butler et al. (2004), Influence of increasing wildfire and fuel break size on fire behavior in shrubland
- **Fire Propagation**: Andrews, P. L. (2018), The Rothermel surface fire spread model and associated developments

## Configuration Files

Each scenario includes:
1. **Main input file**: `inputs_fire_abl_<scenario>` - Complete ERF configuration
2. **Sounding file**: `input_sounding_MRF_<scenario>` - Atmospheric initial conditions

### File Structure
- Geometry specification
- Boundary conditions (periodic + surface layer)
- ABL physics (MRF not active in fire runs)
- Fire module parameters
- Geostrophic wind from sounding
- I/O and diagnostics

## Future Enhancements

Potential extensions to this configuration set:
1. **Stable conditions** with weak fires for low-activity regimes
2. **Variable fuel moisture** scenarios
3. **Terrain-influenced wind** configurations
4. **Crown fire** scenarios (requires extended atmosphere domain)
5. **Multiple fire** scenarios (spotting, fire mergers)

## Contact & Support

For questions or issues with these configurations:
- Check individual input file comments for scenario-specific details
- Verify wind values match sounding file (column 4, row 2)
- Ensure `erf.abl_geo_wind` matches extracted sounding wind speed

---
**Last Updated**: July 2026  
**Status**: Production ready  
**Validation**: Wind parameters verified against sounding files
