
.. role:: cpp(code)
   :language: c++

.. _sec:FireFuelMoisture:

Fire Model - Fuel Moisture Integration
=======================================

Overview
--------

The fuel moisture model provides dynamic updating of fuel moisture content based on atmospheric conditions. This section documents the implementation of fuel moisture integration into the fire simulation pipeline.

Current Implementation Status
-----------------------------

**Implementation Location:**

- **Header:** ``Source/Fire/ERF_FireLayer.H`` (lines 213-225)
- **Implementation:** ``Source/Fire/ERF_FireLayer.cpp`` (lines 238-312)
- **Status:** Functional implementation with proper physics

Function Signature
~~~~~~~~~~~~~~~~~~

.. code-block:: cpp

   void advance_fuel_moisture(amrex::Real dt_s,
                              const amrex::MultiFab& T_atm_k0,
                              const amrex::MultiFab& RH_atm_k0);

Required Inputs
~~~~~~~~~~~~~~~

**dt_s** (Real): Atmospheric timestep [seconds]
   - Available at call site in ``ERF::Advance()``
   - Source: ``dt_lev`` parameter passed to ``FireLayer::advance()``

**T_atm_k0** (MultiFab): Atmospheric potential temperature at k=0 [K]
   - Source: Derived from ``Theta_prim[lev]`` at k=0
   - Type: Potential temperature (not absolute temperature)
   - Grid: Atmospheric level 0 grid

**RH_atm_k0** (MultiFab): Atmospheric relative humidity at k=0 [fraction 0-1]
   - Source: Computed from conserved variables
   - Type: Fraction (0.0 to 1.0), not percentage
   - Grid: Atmospheric level 0 grid

Implementation Details
~~~~~~~~~~~~~~~~~~~~~~

The ``advance_fuel_moisture()`` function:

- Converts timestep from seconds to hours
- Iterates over all fire grid cells using ``MFIter``
- Maps each fire cell to corresponding atmospheric grid cell (accounting for refinement factor C)
- Reads atmospheric temperature and RH at k=0
- Updates fuel moisture for three fuel classes (1hr, 10hr, 100hr) using time-lag ODEs
- Recomputes moisture of extinction
- Handles hysteresis in moisture equilibrium
- Applies temperature and precipitation corrections

Validation
~~~~~~~~~~

- Properly handles grid refinement (fire grid refined by factor C relative to atmospheric grid)
- Converts K to °C for temperature calculations
- Converts RH fraction to percent for internal calculations
- Updates multiple moisture components
- Recomputes moisture of extinction after update
- Supports GPU-compatible kernels with AMREX_DEVICE annotations

Integration Requirements
------------------------

Atmospheric State Extraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``advance()`` method receives ``dt`` and ``surface_layer`` but requires direct access to atmospheric state variables. Available data in ``ERF::Advance()`` context includes:

- ``S_old`` and ``S_new``: Conserved variables (include RhoTheta, RhoQ1)
- ``Theta_prim[lev]``: Primitive potential temperature
- ``Qv_prim[lev]``: Primitive water vapor mixing ratio
- ``Geom(lev)``: Geometry on level 0

Grid Alignment
~~~~~~~~~~~~~~

The fire layer operates on a refined grid (fire cells refined by factor C relative to atmospheric cells). The ``advance_fuel_moisture()`` method handles this internally through mapping:

.. code-block:: cpp

   // Map fire grid cell to atmospheric grid cell
   int i_a = i_f / C;  // i_f = fire grid index, C = refinement factor
   int j_a = j_f / C;

Data Types and Units
~~~~~~~~~~~~~~~~~~~~

**Temperature:**
   - Input: Potential temperature (K)
   - Grid: Atmospheric grid at k=0 (first vertical level)

**Relative Humidity:**
   - Input: Fraction (0.0 to 1.0)
   - Grid: Atmospheric grid at k=0
   - Avoids floating-point precision issues with percentage representation

Implementation Approach
-----------------------

Modified FireLayer::advance() Signature
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**File:** ``Source/Fire/ERF_FireLayer.H`` (lines 98-145)

Updated method signature to accept atmospheric state at k=0:

.. code-block:: cpp

   void advance(amrex::Real dt, 
                SurfaceLayer& surface_layer,
                const amrex::MultiFab& T_atm_k0,
                const amrex::MultiFab& RH_atm_k0);

Documentation includes:

- Parameter descriptions for T_atm_k0 and RH_atm_k0
- Note that atmospheric fields are on coarser grid than fire grid
- Fire computation pipeline step: "Update fuel moisture from atmospheric state"

Fuel Moisture Update in advance()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**File:** ``Source/Fire/ERF_FireLayer.cpp`` (lines 106-200)

Integrated call to ``advance_fuel_moisture()`` in the fire computation pipeline:

.. code-block:: cpp

   // Update fuel moisture from atmospheric state
   if (m_params.fire_debug) {
       amrex::Print() << "[FIRE DEBUG] Updating fuel moisture from atmospheric state" << std::endl;
   }
   advance_fuel_moisture(dt, T_atm_k0, RH_atm_k0);

**Placement:** After wind processing and before rate-of-spread computation

**Rationale:**
   - Fuel moisture update occurs before ROS computation so updated moisture values affect fire spread rates
   - Wind fields are already available and stable before moisture update
   - Updated moisture values flow directly to Rothermel model calculations

Relative Humidity Computation Utility
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**File:** ``Source/Fire/ERF_FireUtils.H``

Function ``compute_rh_from_conservative()`` converts conserved atmospheric variables to relative humidity.

**Input Variables:**
   - Density (Rho_comp)
   - Potential temperature (RhoTheta_comp)
   - Water vapor mixing ratio from RhoQ1_comp

**Computation Steps:**
   1. Extract conserved variables at k=0
   2. Compute mixing ratio: qv = RhoQ1 / Rho
   3. Compute potential temperature: theta = RhoTheta / Rho
   4. Compute pressure: p = getPgivenRTh(RhoTheta, qv)
   5. Convert to saturation functions units: p_mbar = p * 0.01
   6. Compute temperature: T = theta * (p/p_0)^(R_d/c_p) using Exner function
   7. Compute saturation vapor pressure: e_sat = erf_esatw(T)
   8. Compute saturation mixing ratio: q_sat = erf_qsat_from_esat(e_sat, p)
   9. Compute RH: RH = qv / q_sat (clamped to [0, 1])

**Physical Validity:**
   - Uses standard thermodynamic relationships from atmospheric dynamics
   - Leverages ERF's validated EOS functions
   - Uses Flatau polynomial saturation vapor pressure
   - Handles edge cases with division guards and clamping

Call Site Update in ERF_Advance()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**File:** ``Source/TimeIntegration/ERF_Advance.cpp`` (lines 1-420)

Added include:

.. code-block:: cpp

   #ifdef ERF_ENABLE_FIRE
   #include <ERF_FireUtils.H>
   #endif

Modified fire layer advance call (lines 415-428):

.. code-block:: cpp

   #ifdef ERF_ENABLE_FIRE
       // Advance fire simulation at level 0
       if (lev == 0 && m_fire_layer) {
           // Extract atmospheric state at k=0 for fuel moisture update
           const MultiFab& T_atm_k0 = *Theta_prim[lev];  
           
           // Compute relative humidity at k=0 from conserved state
           MultiFab RH_atm_k0(S_old.boxArray(), S_old.DistributionMap(), 1, 0);
           compute_rh_from_conservative(RH_atm_k0, S_old, Geom(lev));
           
           m_fire_layer->advance(dt_lev, *m_SurfaceLayer, T_atm_k0, RH_atm_k0);
       }
   #endif

Data Flow
---------

The fuel moisture update integrates into the fire computation pipeline as follows:

.. code-block:: text

   ERF::Advance() [lev=0, dt=dt_lev]
       |
       ├─ Available: S_old (conserved), Theta_prim[0]
       |
       ├─ Extract: T_atm_k0 = Theta_prim[0]
       |
       ├─ Compute: RH_atm_k0 from (S_old, Geom)
       |   └─ compute_rh_from_conservative()
       |       ├─ Extract conserved at k=0
       |       ├─ Compute p from RhoTheta, qv
       |       ├─ Compute T from theta, p
       |       ├─ Compute RH = qv / q_sat
       |       └─ Clamp to [0, 1]
       |
       └─ Call: FireLayer::advance(dt, surface_layer, T_atm_k0, RH_atm_k0)
           |
           ├─ Extract wind from MOST
           ├─ Apply WAF (optional)
           ├─ Apply terrain corrections (optional)
           |
           ├─→ advance_fuel_moisture(dt, T_atm_k0, RH_atm_k0)  [NEW]
           |   └─ Update MC_1hr, MC_10hr, MC_100hr
           |
           ├─ Compute ROS (now with updated moisture)
           ├─ Advance level-set
           └─ Output diagnostics

Technical Verification
----------------------

Constants and Functions Used
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``Rho_comp``, ``RhoTheta_comp``, ``RhoQ1_comp`` (from ERF_IndexDefines.H)
- ``R_d``, ``Cp_d``, ``p_0`` (from ERF_Constants.H)
- ``getPgivenRTh()``, ``getExnergivenP()`` (from ERF_EOS.H)
- ``erf_esatw()``, ``erf_qsat_from_esat()`` (from ERF_MicrophysicsUtils.H)
- ``ParallelFor``, ``Array4``, ``MultiFab`` (from AMReX)

GPU Compatibility
~~~~~~~~~~~~~~~~~

- All kernels marked with ``AMREX_GPU_DEVICE``
- Uses AMReX parallel loop constructs
- Compatible with CPU and GPU execution

Error Handling
~~~~~~~~~~~~~~

- Division guard: ``if (q_sat > 1.0e-10_rt)`` prevents NaN
- RH clamping: ``amrex::min/max`` ensures valid range [0, 1]
- Geometry parameter available for future domain-specific logic

Testing Checklist
-----------------

Unit Tests
~~~~~~~~~~

- RH computation: verify output in [0, 1] range
- Temperature conversion: verify against analytical formulas
- Grid refinement mapping: verify fire cells map to correct atmospheric cells
- Hysteresis: verify adsorption/desorption path selection
- Bounds: verify moisture stays in [M_MIN, M_MAX]

Integration Tests
~~~~~~~~~~~~~~~~~

- Compile with ERF_ENABLE_FIRE=true
- Run dummy test: ``Exec/CanonicalTests/Fire/test_fire_dummy.py``
- Check diagnostics output includes fuel moisture updates
- Verify no NaNs or unphysical values in output

Regression Tests
~~~~~~~~~~~~~~~~

- Compare fuel moisture evolution against benchmark
- Verify ROS changes with moisture updates
- Check level-set evolution reflects moisture-dependent fire spread
- Monitor conservation: mass/energy bounds

Future Enhancements
-------------------

1. **Temporal Averaging:** Consider time-averaging RH over sub-steps
2. **Vertical Profile:** Extend to multiple vertical levels
3. **Spatial Filtering:** Apply smoothing to reduce grid-scale noise
4. **Diagnostic Output:** Add RH field to plot file for visualization
5. **Validation:** Compare computed RH against sounding measurements

References
----------

- Fuel moisture model: ``Source/Fire/ERF_FuelMoisture.H``
- EOS functions: ``Source/Utils/ERF_EOS.H``
- Microphysics utilities: ``Source/Utils/ERF_MicrophysicsUtils.H``
- Fire layer: ``Source/Fire/ERF_FireLayer.H``
