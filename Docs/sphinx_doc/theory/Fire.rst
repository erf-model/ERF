
 .. role:: cpp(code)
    :language: c++

.. _sec:Fire:

Fire Model
==========

Overview
--------

The fire model in ERF simulates wildfire propagation using a coupled approach combining the Rothermel fire spread model with the FARSITE elliptical fire expansion algorithm. The implementation uses a level-set method to track the fire front interface, allowing for accurate representation of fire behavior in complex terrain.

The fire model operates on a refined grid with adaptive mesh refinement to capture the dynamics of fire spread at appropriate spatial scales. Fire state is advanced each atmospheric timestep, with communications between atmospheric and fire solvers handled through interpolation and mapping functions.

Physical Models
---------------

Rothermel Fire Spread Model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Rothermel model computes the rate of fire spread based on fuel characteristics and environmental conditions. The fire spread rate is determined by:

- **Fuel properties**: Fuel moisture content, bed depth, particle density, energy content, and fuel load
- **Environmental factors**: Wind speed and slope angle
- **Rate calculations**: Separate computation of head fire (downwind), flank fire (cross-wind), and backing fire (upwind) spread rates

The Rothermel model is based on:

Rothermel, R. C. (1972). A mathematical model for predicting fire spread in wildland fuels. Res. Paper INT-115, USDA Forest Service, Intermountain Forest and Range Experiment Station.

FARSITE Elliptical Fire Expansion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The FARSITE algorithm models fire expansion as an elliptical shape with characteristics derived from wind and topography:

- **Ellipticity**: The length-to-width ratio of the fire ellipse is derived from wind speed using the Anderson (1983) formulas
- **Directional spread rates**: The ellipse is oriented according to wind direction and expanding at rates corresponding to head, flank, and backing spread rates computed by Rothermel
- **Coefficients**: The Richards (1990) coefficients relate head fire rate of spread to flank and backing rates through the ellipse geometry
- **Level-set representation**: The fire front is represented as a signed-distance function where negative values indicate burned area, positive values indicate unburned fuel, and values near zero indicate the active fire front

References:

- Finney, M. A. (2004). FARSITE: Fire Area Simulator model development and evaluation. Res. Paper RMRS-RP-4 Revised, USDA Forest Service, Rocky Mountain Research Station.
- Richards, G.D. (1990). An elliptical growth model of forest fire fronts and its numerical solution. Int. J. Numer. Meth. Eng. 30(6):1163-1179.

Level-Set Propagation
~~~~~~~~~~~~~~~~~~~~~~

The fire front is tracked using a level-set method where the signed-distance function :math:`\phi` represents the fire front location:

- :math:`\phi < 0`: Burned region (inside fire)
- :math:`\phi > 0`: Unburned fuel (outside fire)  
- :math:`\phi \approx 0`: Fire front interface

The level-set field is advanced by computing spread vectors at each grid point based on the local rate of spread, elliptical orientation, and elevation. A two-pass algorithm handles GPU computation of spread vectors followed by host-side MPI gather operations for distributed computing.

Wind and Terrain Effects
~~~~~~~~~~~~~~~~~~~~~~~~~

The model includes adjustments for wind speed variation with height and terrain effects:

- **Wind Adjustment Factor (WAF)**: Adjusts wind speed from reference height to the flame zone where fire spread occurs, based on fuel moisture and other properties
- **Terrain corrections**: FARSITE terrain wind corrections account for ridge speed-up, sheltered wind reduction, valley channeling, and wind deflection effects

Model State Variables
---------------------

The fire model maintains several grid-based state variables on the fire-refined mesh:

.. list-table:: Fire Grid MultiFabs (dimensions: Cf × nx × Cf × ny × 1)
   :widths: 20 15 20 45
   :header-rows: 1

   * - Field
     - Components
     - Ghost Cells
     - Description
   * - fire_phi
     - 1
     - 1
     - Level-set field: <0 burned, >0 unburned, ≈0 front
   * - fire_wind_ref
     - 2
     - 0
     - Reference wind components (u, v) at specified height
   * - fire_wind_eff
     - 2
     - 0
     - Effective wind after WAF and terrain corrections
   * - fire_slopes
     - 2
     - 0
     - Terrain slopes (dz/dx, dz/dy)
   * - fire_curvature
     - 1
     - 0
     - Terrain curvature
   * - fire_ros
     - 1
     - 0
     - Rate of spread [m/s]
   * - fire_fuel_load
     - 1
     - 0
     - Fuel load [kg/m²]
   * - fire_fuel_mc
     - 3
     - 0
     - Fuel moisture content (1-hr, 10-hr, 100-hr) [fraction]
   * - fire_heat_flux
     - 1
     - 0
     - Heat flux [W/m²]
   * - fire_spread_vec
     - 2
     - 0
     - Spread vectors for FARSITE expansion

Implementation Components
--------------------------

Core Classes
~~~~~~~~~~~~

**Fire Layer Class** (ERF_FireLayer.H, ERF_FireLayer.cpp)
   - Main fire simulation container
   - Manages fire state on refined grid
   - Implements fire computation pipeline
   - Handles communication between atmospheric and fire solvers

**Fire Parameters** (ERF_FireParams.H)
   - Stores fire model configuration and parameters
   - Reads user settings from input file
   - Manages fuel model specifications

**Wind Extraction** (ERF_FireWindExtract.H, ERF_FireWindExtract.cpp)
   - Extracts wind field from atmospheric MOST layer
   - Interpolates wind to fire grid
   - Applies Wind Adjustment Factor

**Terrain Handling** (ERF_FireTerrainReader.H, ERF_FireTerrainReader.cpp)
   - Reads terrain elevation data
   - Computes slopes and curvature on fire grid
   - Applies FARSITE terrain corrections

**Ignition** (ERF_FireIgnition.H)
   - Initializes fire front at specified location and time
   - Sets up initial level-set field

**Plotfile Output** (ERF_FirePlotfile.H, ERF_FirePlotfile.cpp)
   - Writes fire state variables to plotfiles
   - Enables visualization and analysis of results

Integration with ERF
~~~~~~~~~~~~~~~~~~~~

The fire layer is instantiated and managed by the main ERF class:

- **Initialization**: Fire layer is created in ``ERF::InitData_post()`` after atmospheric grid setup
- **Advance**: Fire state is advanced each atmospheric timestep via ``FireLayer::advance()``
- **Wind coupling**: Wind field is extracted from the MOST boundary layer model
- **Output**: Fire variables are included in plotfile output for visualization

Configuration and Input Parameters
-----------------------------------

Fire model behavior is controlled through input file parameters using the ``erf.fire.*`` prefix.

Basic Configuration
~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   erf.fire.enable = true                # Enable/disable fire model
   erf.fire.grid_ratio = 5               # Fire grid refinement factor relative to atmospheric grid
   erf.fire.fuel_model_id = 1            # Anderson FBFM13 fuel model (1-13)
   erf.fire.ignition_x = 1000.0          # Ignition center x-coordinate [m]
   erf.fire.ignition_y = 1000.0          # Ignition center y-coordinate [m]
   erf.fire.ignition_r = 20.0            # Initial fire radius [m]
   erf.fire.ignition_time = 0.0          # Ignition time [s]

Fuel and Moisture
~~~~~~~~~~~~~~~~~

.. code-block:: text

   erf.fire.fuel_model_id = 1            # Anderson fuel model (1-13)
   erf.fire.moisture_1hr = 0.08          # 1-hour fuel moisture [fraction]
   erf.fire.moisture_10hr = 0.08         # 10-hour fuel moisture [fraction]
   erf.fire.moisture_100hr = 0.10        # 100-hour fuel moisture [fraction]

Wind Parameters
~~~~~~~~~~~~~~~

.. code-block:: text

   erf.fire.wind_ref_ht = 6.1            # Reference height for wind input [m]
   erf.fire.use_waf = true               # Apply Wind Adjustment Factor
   erf.fire.waf_formula = "andrews"      # WAF formula: "andrews" or "behaviorplus"

Terrain Corrections
~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   erf.fire.use_terrain_wind = true      # Apply FARSITE terrain corrections
   erf.fire.k_ridge = 1.5                # Ridge speed-up factor
   erf.fire.k_shelter = 0.6              # Sheltered wind reduction factor
   erf.fire.k_valley = 0.8               # Valley channeling factor
   erf.fire.k_deflect = 0.3              # Wind deflection factor

FARSITE Algorithm Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   erf.fire.farsite.phi_threshold = 0.0     # Level-set detection threshold
   erf.fire.farsite.use_anderson_lw = 1     # Use Anderson L/W ratio (1) or fixed coefficients (0)
   erf.fire.farsite.coeff_a = 0.5           # Richards head fire coefficient
   erf.fire.farsite.coeff_b = 0.25          # Richards flank fire coefficient
   erf.fire.farsite.coeff_c = 0.1           # Richards backing fire coefficient
   erf.fire.farsite.gaussian_sigma = -1.0   # Phi stamping radius [m]
                                            # <0: single-cell mode
                                            # =0: auto mode
                                            # >0: fixed Gaussian radius
   erf.fire.farsite.cfl_fire = 0.5          # Fire CFL number for subcycle timestep

Debugging
~~~~~~~~~

.. code-block:: text

   erf.fire.fire_debug = false           # Enable debug output for fire calculations

When enabled, debug messages provide information about wind extraction, wind adjustment factor application, terrain corrections, rate-of-spread calculations, and level-set propagation.

Development Phases
------------------

The fire model has been developed in phases to progressively implement functionality:

**Phase 1: Framework**
   - Fire layer class structure
   - Dummy Rothermel and FARSITE functions
   - Integration with main ERF class
   - Input file support
   - Basic regression tests

**Phase 2: Rothermel Model**
   - Full Rothermel fire spread equations
   - Wind and slope factor calculations
   - Anderson FBFM13 fuel model database
   - Wind Adjustment Factor application
   - Terrain slope computation

**Phase 3: Level-Set Propagation**
   - Level-set field initialization
   - FARSITE elliptical propagation kernel
   - Anderson L/W ratio computation
   - Two-pass GPU/MPI propagation
   - CFL-limited subcycling
   - Single-cell and Gaussian stamping modes

Testing
-------

Regression tests verify fire model functionality at different development stages:

**Phase 2 Regression Test** (inputs_fire_phase2)
   - Flat domain with GR1 fuel model
   - 5 m/s wind, 8% fuel moisture
   - Verifies Rothermel rate-of-spread computation
   - Output: max_ROS and mean_ROS values

**Phase 3 Regression Test** (inputs_fire_phase3)
   - Flat domain with GR1 fuel model
   - 5 m/s wind, 8% fuel moisture
   - Verifies level-set propagation and FARSITE elliptical expansion
   - Test checks:
      - Fire front propagation (phi_min < 0 at t > 0)
      - Subcycle count verification
      - Burned area geometry at specified time
      - No NaN values in computed fields

Expected Results
~~~~~~~~~~~~~~~~

At t = 1800 seconds with GR1 fuel, 5 m/s wind, and flat terrain:

   - Head fire rate of spread: 0.06-0.10 m/s (WAF-dependent)
   - Adjusted wind speed: 5-6 m/s
   - Anderson L/W ratio: 1.5-2.0
   - Burned ellipse major axis: 200-350 m
   - Burned ellipse minor axis: 100-200 m

Output and Visualization
------------------------

Fire simulation results are written to plotfiles for visualization and analysis:

- **Level-set field** (phi): Shows burned area (negative), unburned fuel (positive), and fire front location
- **Rate of spread** (fire_ros): Spatial distribution of fire spread rates
- **Wind fields** (fire_wind_ref, fire_wind_eff): Reference and effective wind components
- **Fuel properties** (fire_fuel_load, fire_fuel_mc): Fuel load and moisture content
- **Terrain data** (fire_slopes, fire_curvature): Elevation derivatives

These variables enable quantitative analysis of fire behavior and validation against observational data.

References
----------

- Rothermel, R. C. (1972). A mathematical model for predicting fire spread in wildland fuels. Res. Paper INT-115, USDA Forest Service, Intermountain Forest and Range Experiment Station.

- Finney, M. A. (2004). FARSITE: Fire Area Simulator model development and evaluation. Res. Paper RMRS-RP-4 Revised, USDA Forest Service, Rocky Mountain Research Station.

- Richards, G.D. (1990). An elliptical growth model of forest fire fronts and its numerical solution. Int. J. Numer. Meth. Eng. 30(6):1163-1179.

- Andrews, P. L. (2018). Current status and future needs of the BehavePlus Fire Modeling System. International Journal of Wildland Fire, 27(9), 558-566.

- Anderson, H.E. (1983). Predicting wind-driven wild fire size and shape. USDA Forest Service Research Paper INT-305.
