ERF-Dust 2D Surface Emission Module
====================================

Overview
--------

The ERF-Dust module provides a 2D surface dust emission layer that operates
on the horizontal domain of the ERF 3D atmospheric solver. The emission layer
exchanges fields with the 3D solver through a two-way coupling interface
modelled on the wildfire coupling in the ERF-Fire branch.

Geochemical inputs from PHREEQC are read as pre-computed files at prescribed
intervals. No runtime calls to PHREEQC are made within the ERF timestepping
loop. This one-way file-based coupling is appropriate because geochemical
processes evolve on timescales (days to months) that are much longer than
the atmospheric dynamical timestep (seconds).

The module is designed for dust emission from domestic critical mineral
extraction sites, including open-pit mines, tailings impoundments, haul
roads, and evaporation ponds.

Coupling Architecture
---------------------

.. code-block:: text

   PHREEQC (offline)  --[file, periodic interval]-->  ERF-Dust 2D
                                                       |        ^
                                               upward  |        | downward
                                               aerosol |        | u*, T_sfc,
                                               flux    |        | H, h_PBL,
                                                       v        | wind
                                                   ERF-Atm 3D

The ERF-Dust 2D layer injects aerosol mass flux into the
``AdvectionSrcForScalars`` source term at the lowest model level. The 3D
solver returns friction velocity :math:`u^*`, surface temperature
:math:`T_\mathrm{sfc}`, sensible heat flux :math:`H`, PBL height :math:`h`
(diagnosed by the MRF scheme), and near-surface wind to the emission layer
at every timestep.

Boundary Layer Parameterisation
---------------------------------

The MRF (Medium Range Forecast) PBL scheme [Hong1996]_ is used. The scheme
provides:

- Bulk Richardson number :math:`Ri_b` diagnosis of PBL height :math:`h`.
- Stability-dependent eddy diffusivity:

  .. math::

     K_m = \rho \, w_* \, \kappa \, z \left(1 - \frac{z}{h}\right)^2

- Nonlocal countergradient flux corrections for heat (HGAMT) and
  moisture (HGAMQ) when ``enable_mrf_countergradient = true``.
- High-resolution grid-dependent diffusivity bounds at fine AMR levels
  when ``pbl_mrf_highres_bounds = true``.

In Phase 14, the countergradient framework is extended to dust scalar
tracers, enabling nonlocal vertical transport through convective PBL
structures.

Dust Emission Physics
---------------------

The threshold friction velocity :math:`u^*_t` below which no emission
occurs is computed using the Bagnold [Bagnold1941]_ formula:

.. math::

   u^*_t = A \sqrt{\frac{\rho_p \, g \, d}{\rho_a}}

where :math:`A` is the threshold coefficient (default 0.0123), :math:`\rho_p`
is the particle density (default 2650 kg/m³ for quartz), :math:`g` is
gravitational acceleration, :math:`d` is the particle diameter, and
:math:`\rho_a` is the air density.

The threshold is modified by surface crust strength, salt efflorescence
fraction, soil moisture, and suppression agent coverage. These modifiers
are provided by PHREEQC output files (Phase 4) and by ERF surface layer
fields (Phases 5 and 8).

Size Bins
---------

Three particle size bins are used by default:

+-------+----------------------------------+---------------------+
| Bin   | Diameter range                   | Representative d    |
+=======+==================================+=====================+
| 0     | 2.5 µm < d ≤ 10 µm (PM10)       | 7.0 µm              |
+-------+----------------------------------+---------------------+
| 1     | 1.0 µm < d ≤ 2.5 µm             | 3.5 µm              |
+-------+----------------------------------+---------------------+
| 2     | d ≤ 1.0 µm (PM2.5)              | 0.7 µm              |
+-------+----------------------------------+---------------------+

The number of bins is controlled by ``erf.dust.n_size_bins``. Each bin
maps to one ``RhoQ*_comp`` scalar component in the ERF conserved state
vector.

Input Parameters
----------------

All parameters are read from the ``erf.dust`` ParmParse prefix.

.. list-table::
   :header-rows: 1
   :widths: 35 15 15 35

   * - Parameter
     - Type
     - Default
     - Description
   * - ``erf.dust.enable``
     - bool
     - false
     - Enable the dust module
   * - ``erf.dust.dust_debug``
     - bool
     - false
     - Write per-timestep debug fields to plotfile
   * - ``erf.dust.n_size_bins``
     - int
     - 3
     - Number of particle size bins
   * - ``erf.dust.particle_density``
     - Real
     - 2650.0
     - Bulk particle density [kg/m³]
   * - ``erf.dust.z0_dust``
     - Real
     - 0.01
     - Aerodynamic roughness length [m]
   * - ``erf.dust.silt_fraction``
     - Real
     - 0.10
     - Surface silt mass fraction [dimensionless]
   * - ``erf.dust.threshold_A_coeff``
     - Real
     - 0.0123
     - Bagnold threshold coefficient A [dimensionless]
   * - ``erf.dust.crust_index``
     - Real
     - 0.0
     - Surface crust strength index [0,1]
   * - ``erf.dust.grid_ratio``
     - int
     - 1
     - Refinement factor of dust grid relative to atmospheric grid [dimensionless]
   * - ``erf.dust.phreeqc_update_interval_s``
     - Real
     - 86400.0
     - Interval for re-reading PHREEQC files [s]
   * - ``erf.dust.surface_map_file``
     - string
     - ""
     - Path to NetCDF surface property map
   * - ``erf.dust.phreeqc_output_file``
     - string
     - ""
     - Path to PHREEQC output NetCDF
   * - ``erf.dust.soil_type_file``
     - string
     - ""
     - Path to soil type raster (.asc or .nc)
   * - ``erf.dust.silt_fraction_file``
     - string
     - ""
     - Path to silt fraction raster. Empty = uniform silt_fraction.
   * - ``erf.dust.crust_index_file``
     - string
     - ""
     - Path to crust index raster. Empty = uniform crust_index.
   * - ``erf.dust.moisture_flag_file``
     - string
     - ""
     - Path to moisture inhibition flag raster. Empty = uniform 0.
   * - ``erf.dust.suppression_file``
     - string
     - ""
     - Path to suppression coverage raster. Empty = uniform 0.

2D Dust Grid
------------

The dust emission layer operates on a 2D horizontal slab (k-extent = 1)
that covers the same physical domain as the ERF level-0 atmospheric grid.
The dust grid can be refined relative to the atmospheric grid by the factor
``erf.dust.grid_ratio``. With ``grid_ratio = 1`` (default), one dust cell
corresponds to one atmospheric cell. With ``grid_ratio = N``, each
atmospheric cell is subdivided into N² dust cells, allowing sub-atmospheric-
cell resolution of surface emission heterogeneity (e.g., pit walls, road
segments, pond boundaries).

The MPI decomposition of the dust grid follows the atmospheric
``DistributionMapping``, so that each MPI rank owns the same horizontal
tiles in both grids. This avoids inter-rank communication during wind
extraction (Phase 9) and atmospheric flux injection (Phase 10).

All atmospheric box sizes in x and y must be divisible by ``grid_ratio``.
If this constraint is not satisfied, the prerequisite check in
``verify_dust_prerequisites`` aborts the simulation with a diagnostic
message.

The prerequisite verification function (``verify_dust_prerequisites``)
checks eight conditions before the DustLayer is allocated:

1. ``SurfaceLayer`` is initialized (``zlo.type = "surface_layer"`` required).
2. Roughness length ``erf.most.z0`` is set.
3. No z-direction MPI decomposition (``amrex.max_grid_size_z = domain_nz``).
4. ``erf.dust.grid_ratio >= 1``.
5. All atmospheric box x,y sizes are divisible by ``grid_ratio``.
6. ``DistributionMapping`` size matches ``BoxArray`` size.
7. Domain z-index starts at 0.
8. Domain physical height exceeds the MRF reference height ``erf.most.zref``.

Dust State Variables
--------------------

The following ``MultiFab``\s are allocated on the 2D dust grid in Phase 2.
All fields have 1 ghost cell in x and y (z ghost cells are not used). The
number of components in ``dust_emission_flux`` equals ``erf.dust.n_size_bins``.

.. list-table::
   :header-rows: 1
   :widths: 30 10 60

   * - Field
     - ncomp
     - Description
   * - ``dust_ustar_t``
     - 1
     - Threshold friction velocity :math:`u^*_t` [m/s]. Initialized to the
       Bagnold value for bin 0 (coarsest bin). Updated by surface chemistry
       and moisture in Phases 4 and 5.
   * - ``dust_soil_type``
     - 1
     - Soil type index stored as Real. Codes 1--16 follow STATSGO. Codes
       :math:`\geq 100` are reserved for mine-specific surface types
       (assigned in Phase 3).
   * - ``dust_silt_fraction``
     - 1
     - Surface silt mass fraction [--]. Initialized from
       ``erf.dust.silt_fraction``. Updated by surface reader in Phase 3.
   * - ``dust_crust_index``
     - 1
     - Mineral crust strength index [0, 1]. Initialized from
       ``erf.dust.crust_index``. Updated by PHREEQC reader in Phase 4.
   * - ``dust_moisture_flag``
     - 1
     - Surface moisture inhibition flag [0, 1]. Set to 0 at initialization
       (dry surface). Updated from ERF surface layer in Phase 5.
   * - ``dust_suppression``
     - 1
     - Suppression agent coverage fraction [--]. Set to 0 at initialization.
       Updated by suppression model in Phase 8.
   * - ``dust_emission_flux``
     - ``n_size_bins``
     - Vertical dust emission flux per size bin [kg/m²/s]. Set to 0 at
       initialization. Computed by emission model in Phase 6.

Surface Property Maps
---------------------

Phase 3 adds support for reading spatial heterogeneity in surface properties
from raster files. Maps are read during ``DustLayer::initialize`` and
interpolated onto the 2D dust grid via bilinear interpolation.

Supported File Formats
~~~~~~~~~~~~~~~~~~~~~~

**ESRI ASCII Grid** (.asc)

- Simple ASCII format with 6-line header followed by data rows.
- Header format (case-insensitive keys):

  .. code-block:: text

     ncols <integer>
     nrows <integer>
     xllcorner <real>
     yllcorner <real>
     cellsize <real>
     nodata_value <real>

- Data rows are arranged northernmost-to-southernmost (row 0 is at y_max).
  The reader automatically reverses rows to match the ERF domain convention
  (row 0 is at y_min), following the pattern in ``ERF_FuelMap.H``.
- Reference: `ESRI ASCII Raster Format
  <https://desktop.arcgis.com/en/arcmap/latest/manage-data/raster-and-images/
  esri-ascii-raster-format.htm>`_

**NetCDF** (.nc)

- Requires ``ERF_ENABLE_NETCDF=ON`` at CMake configure time.
- Implementation deferred to Phase 4; currently aborts if requested.
- Variable naming convention to be determined.

Row Reversal Convention
~~~~~~~~~~~~~~~~~~~~~~~

The ESRI ASCII format stores the northernmost row first (row 0 in file = y_max).
The reader reverses rows during I/O so that row 0 in memory = y_min, consistent
with the ERF domain layout and the ``ERF_FuelMap.H`` convention.

Surface Type Codes
~~~~~~~~~~~~~~~~~~~

The ``soil_type`` field uses integer codes to identify surface categories:

+----------+-------------------------------------------------------------+
| Code     | Category                                                    |
+==========+=============================================================+
| 0        | Undefined (default, no emissions)                           |
+----------+-------------------------------------------------------------+
| 1--16    | STATSGO natural soil types. Reference:                      |
|          | `STATSGO2 Data                                              |
|          | <https://www.nrcs.usda.gov/resources/data-and-reports/      |
|          | ssurgo/statsgo2-data>`_                                    |
+----------+-------------------------------------------------------------+
| 100      | Mine tailings (general)                                     |
+----------+-------------------------------------------------------------+
| 101      | Lithium brine pond                                          |
+----------+-------------------------------------------------------------+
| 102      | Rare-earth-elements (REE) tailings                          |
+----------+-------------------------------------------------------------+
| 103      | Copper tailings                                             |
+----------+-------------------------------------------------------------+
| 104      | Unpaved haul road. Reference: `EPA AP-42 Chapter 13.2.2     |
|          | <https://www.epa.gov/air-emissions-factors-and-             |
|          | quantification/ap-42-compilation-air-emission-factors>`_   |
+----------+-------------------------------------------------------------+

Surface Property Map Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Maps are specified via five optional file path parameters, all under the
``erf.dust`` ParmParse prefix. If a path is empty (default), the corresponding
MultiFab retains its uniform default value set during initialization.

.. list-table::
   :header-rows: 1
   :widths: 25 10 65

   * - Parameter
    - Type
    - Description
   * - ``erf.dust.soil_type_file``
    - string
    - Path to soil type raster (.asc or .nc). Code values as above.
   * - ``erf.dust.silt_fraction_file``
    - string
    - Path to silt fraction raster. Values should be in [0,1]. Empty =
      uniform ``erf.dust.silt_fraction``.
   * - ``erf.dust.crust_index_file``
    - string
    - Path to crust strength index raster. Values in [0,1] where 0 = uncrusted,
      1 = fully crusted. Empty = uniform ``erf.dust.crust_index``.
   * - ``erf.dust.moisture_flag_file``
    - string
    - Path to surface moisture inhibition flag raster. Values in [0,1].
      Empty = uniform 0.0 (dry surface).
   * - ``erf.dust.suppression_file``
    - string
    - Path to suppression agent coverage fraction raster. Values in [0,1].
      Empty = uniform 0.0 (no suppression).

I/O Pattern
~~~~~~~~~~~

Reading follows the ``ERF_FireTerrainReader.cpp`` and ``ERF_FuelMap.H`` pattern:

1. **Rank-0 I/O**: The MPI rank with ``ParallelDescriptor::IOProcessor() == true``
   opens and reads the file.
2. **Broadcast**: Grid dimensions (``ncols``, ``nrows``, corner coordinates,
   cellsize) are broadcast to all ranks via ``ParallelDescriptor::Bcast``.
3. **GPU Copy**: Data is copied to device via ``Gpu::copy(Gpu::hostToDevice, ...)``.
4. **Bilinear Interpolation**: A GPU kernel loops over MultiFab tiles and
   computes values at dust-grid cell centers using bilinear interpolation.

PHREEQC Geochemical Coupling
----------------------------

PHREEQC (Parkhurst & Appelo 2013) is used offline to compute geochemical
speciation, mineral crust strength, salt efflorescence state, and toxic
metal composition of mine tailings surfaces. ERF-Dust reads PHREEQC
output files at prescribed intervals and updates the corresponding dust
grid ``MultiFab``\s.

Reference:
  Parkhurst, D.L., & Appelo, C.A.J. (2013). PHREEQC version 3.
  USGS Techniques and Methods, book 6, chap. A43.
  https://pubs.usgs.gov/tm/06/a43/

Timescale Justification
~~~~~~~~~~~~~~~~~~~~~~~

Geochemical processes (crust formation, efflorescence, salt precipitation)
evolve over days to weeks. The ERF atmospheric timestep is seconds. Periodic
file-based updates at intervals of hours to days are sufficient to capture
the evolution of surface chemistry. No runtime PHREEQC calls are made inside
the ERF timestepping loop.

Threshold Friction Velocity Update
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After reading ``crust_index`` and ``efflorescence`` from PHREEQC output,
``update_ustar_t_from_chemistry`` applies:

.. math::

   u^*_{t,\mathrm{new}} = u^*_{t,\mathrm{base}}
     \times (1 - \alpha_c \cdot C_I)
     \times (1 - \alpha_e \cdot E_f)

where :math:`C_I` is the crust index [0,1], :math:`E_f` is the efflorescence
fraction [0,1], and :math:`\alpha_c`, :math:`\alpha_e` are reduction
coefficients set by ``erf.dust.alpha_crust`` and ``erf.dust.alpha_efflor``.
The result is clamped to ``PhreeqcDustConst::USTAR_T_MIN`` = 0.05 m/s.

Reference:
  Marticorena, B., & Bergametti, G. (1995). J. Geophys. Res., 100, 16415.
  https://doi.org/10.1029/95JD00690

PHREEQC ParmParse Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Parameter
     - Description
   * - ``erf.dust.phreeqc_output_file``
     - Path to PHREEQC output file (.csv or .nc). Empty = no update.
   * - ``erf.dust.phreeqc_update_interval_s``
     - Interval between PHREEQC file reads [s]. Default = 86400 (1 day).
   * - ``erf.dust.alpha_crust``
     - Crust reduction coefficient for u*_t [-]. Default = 0.5.
   * - ``erf.dust.alpha_efflor``
     - Efflorescence reduction coefficient for u*_t [-]. Default = 0.3.
   * - ``erf.dust.phreeqc_crust_var``
     - CSV column / NetCDF variable name for crust index.
   * - ``erf.dust.phreeqc_silt_var``
     - CSV column / NetCDF variable name for silt fraction.
   * - ``erf.dust.phreeqc_efflor_var``
     - CSV column / NetCDF variable name for efflorescence fraction.
   * - ``erf.dust.phreeqc_supp_var``
     - CSV column / NetCDF variable name for suppression modifier.
   * - ``erf.dust.phreeqc_metal_var``
     - CSV column / NetCDF variable name for toxic metal mass fraction (bin 0).

Threshold Friction Velocity Computation (Phase 5)
---------------------------------------------------

Phase 5 implements the ``ERF_DustThreshold.H`` module that computes the per-cell
threshold friction velocity :math:`u^*_t` incorporating chemistry, moisture,
and suppression modifiers. The computation chain is applied after PHREEQC
updates and surface moisture changes.

Base Threshold (Bagnold 1941)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The base threshold friction velocity for emission onset is:

.. math::

   u^*_{t,\mathrm{base}} = A \sqrt{\frac{\rho_p \, g \, d}{\rho_a}}

where:
   - :math:`A` = Bagnold threshold coefficient (default 0.0123) [-]
   - :math:`\rho_p` = particle density (default 2650 kg/m³) [kg/m³]
   - :math:`g` = gravitational acceleration (9.81) [m/s²]
   - :math:`d` = particle diameter [m]
   - :math:`\rho_a` = air density (1.225 at sea level) [kg/m³]

References:
  Bagnold, R. A. (1941). *The Physics of Blown Sand and Desert Dunes*.
    Methuen, London.

Modifier Chain
~~~~~~~~~~~~~~

The base threshold is adjusted by three dimensionless factors:

**Chemistry (PHREEQC, Phase 4):**

.. math::

   f_{\mathrm{chem}} = (1 - \alpha_c \cdot C_I) \times (1 - \alpha_e \cdot E_f)

where :math:`C_I` is the crust strength index [0,1], :math:`E_f` is
efflorescence fraction [0,1], and :math:`\alpha_c`, :math:`\alpha_e` are
reduction coefficients (defaults 0.5 and 0.3). This reduction is justified
by Marticorena & Bergametti (1995): crusty and efflorescent surfaces
present higher cohesion and thus higher thresholds.

**Moisture Inhibition (Phase 5):**

.. math::

   f_{\mathrm{moist}} = 1 + \beta_m \cdot w

where :math:`w` is surface moisture [0,1] (clamped) and :math:`\beta_m`
is the moisture inhibition coefficient (default 4.0, dimensionless).
This formulation follows Fecan et al. (1999): soil moisture raises
capillary cohesion, suppressing emission. The additive structure ensures
:math:`f_{\mathrm{moist}} \geq 1`, consistently increasing :math:`u^*_t`
when moisture is present.

**Suppression Elevation (Phase 8):**

.. math::

   f_{\mathrm{supp}} = 1 + \gamma_s \cdot s

where :math:`s` is suppression agent coverage [0,1] (clamped) and
:math:`\gamma_s` is the suppression elevation coefficient (default 6.0,
dimensionless). This structure ensures :math:`f_{\mathrm{supp}} \geq 1`,
consistently raising :math:`u^*_t` to inhibit emission when suppression
is applied.

Final Threshold
~~~~~~~~~~~~~~~

The final per-cell threshold is computed as:

.. math::

   u^*_t = u^*_{t,\mathrm{base}} \times \frac{f_{\mathrm{chem}}}{f_{\mathrm{moist}} \times f_{\mathrm{supp}}}

clamped to the valid range [:math:`u^*_{t,\mathrm{min}}`, :math:`u^*_{t,\mathrm{max}}`],
where:
   - :math:`u^*_{t,\mathrm{min}}` = 0.05 m/s (floor; very fine particles)
   - :math:`u^*_{t,\mathrm{max}}` = 5.0 m/s (ceiling; very coarse or cohesive)

The clamp prevents unphysical regimes and ensures numerical stability.

Implementation
~~~~~~~~~~~~~~~

The computation is implemented as a header-only GPU kernel in
``Source/Dust/ERF_DustThreshold.H``. The module provides:

- ``compute_ustar_t_bagnold``: Device function for base Bagnold computation.
- ``compute_moisture_inhibition``: Device function for moisture factor.
- ``compute_suppression_factor``: Device function for suppression factor.
- ``compute_ustar_t_full``: Device function combining all modifiers and clamping.
- ``recompute_dust_ustar_t``: MultiFab GPU kernel called at each timestep.

The kernel is invoked in ``DustLayer::initialize`` (after surface maps
are populated) and in ``DustLayer::advance`` (after PHREEQC updates and
moisture changes). This ensures :math:`u^*_t` is always consistent with
the current state of surface chemistry, moisture, and suppression.

Phase 5 References
~~~~~~~~~~~~~~~~~~~

  Bagnold, R. A. (1941). *The Physics of Blown Sand and Desert Dunes*.
    Methuen, London.

  Marticorena, B., & Bergametti, G. (1995). Modeling the atmospheric
    dust cycle: 1. Design of a soil-derived dust emission scheme.
    *J. Geophys. Res.*, 100, 16415-16430.
    https://doi.org/10.1029/95JD00690

  Fecan, F., Marticorena, B., & Bergametti, G. (1999). Parametrization
    of the increase of aeolian erosion threshold wind friction velocity
    due to soil moisture for arid and semi-arid areas.
    *Ann. Geophys.*, 17, 149-157.
    https://doi.org/10.1007/s00585-999-0149-7

  Shao, Y., & Lu, H. (2000). A simple expression for wind erosion
    threshold friction velocity. *J. Geophys. Res.*, 105, 22437-22443.
    https://doi.org/10.1029/2000JD900304

Saltation and Vertical Emission Flux (Phase 6)
-----------------------------------------------

Phase 6 implements the Marticorena & Bergametti (1995) saltation bombardment
model for dust emission. Once the friction velocity :math:`u^*` (placeholder
in Phase 6, extracted from MRF in Phase 9) exceeds the threshold :math:`u^*_t`
computed in Phase 5, the horizontal saltation flux and vertical emission per
size bin are computed.

Horizontal Saltation Flux
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Owen (1964) form adopted by Marticorena & Bergametti (1995):

.. math::

  Q_s = C_s \frac{\rho_a}{g} {u^*}^3
        \left(1 - \frac{u^*_t}{u^*}\right)
        \left(1 + \frac{u^*_t}{u^*}\right)^2

where:
  - :math:`C_s` = saltation efficiency coefficient (default 2.61) [-]
  - :math:`\rho_a` = air density [kg/m³]
  - :math:`g` = gravitational acceleration (9.81) [m/s²]
  - :math:`u^*` = friction velocity [m/s]
  - :math:`u^*_t` = threshold friction velocity [m/s]
  - :math:`Q_s = 0` when :math:`u^* \leq u^*_t`.

Units: [kg m⁻¹ s⁻¹]

References:
  Owen, P. R. (1964). Saltation of uniform grains in air.
   *J. Fluid Mech.*, 20, 225-242.
   https://doi.org/10.1017/S0022112064001173

  Marticorena, B., & Bergametti, G. (1995). Modeling the atmospheric
   dust cycle: 1. Design of a soil-derived dust emission scheme.
   *J. Geophys. Res.*, 100, 16415-16430.
   https://doi.org/10.1029/95JD00690

Vertical Emission per Size Bin
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

  F_i = \alpha_i \cdot f_\text{silt} \cdot Q_s

where:
  - :math:`\alpha_i` = sandblasting efficiency [m⁻¹]
  - :math:`f_\text{silt}` = surface silt mass fraction [0,1]
  - :math:`Q_s` = horizontal saltation flux [kg m⁻¹ s⁻¹]

The sandblasting efficiency is computed empirically from silt fraction:

.. math::

  \log_{10}(\alpha_i) = 0.134 \cdot f_\text{clay} - 6.0

For mine tailings surfaces, clay fraction is approximated as:

.. math::

  f_\text{clay} = 0.2 \cdot f_\text{silt}

Emission flux is clamped to a maximum of 10⁻³ kg/(m² s) to prevent
physically unreasonable values that could drive numerical instabilities.

In Phase 6, all size bins use the same sandblasting efficiency.
Phase 11 will apply per-bin size-dependent corrections.

Units: [kg m⁻² s⁻¹]

Placeholder Friction Velocity in Phase 6
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In Phase 6, the friction velocity :math:`u^*` is a uniform prescribed
value from the parameter ``erf.dust.test_ustar`` (default 0, no emission).
This placeholder allows testing of the emission physics before Phase 9
implements actual :math:`u^*` extraction from the MRF PBL surface layer.

To enable emission in test cases, set ``erf.dust.test_ustar`` to a value
exceeding the local :math:`u^*_t` (typically 0.3 m/s or higher).

Phase 6 Implementation
~~~~~~~~~~~~~~~~~~~~~~~

The computation is implemented as a GPU kernel in ``Source/Dust/ERF_DustEmission.H``
(header-only device functions) and ``Source/Dust/ERF_DustEmission.cpp``
(MultiFab kernel). The module provides:

- ``compute_saltation_flux``: Device function for horizontal saltation flux.
- ``compute_sandblasting_efficiency``: Device function for alpha from silt fraction.
- ``compute_vertical_emission_flux``: Device function for vertical emission per bin.
- ``compute_dust_emission_flux``: MultiFab GPU kernel called each timestep.

The kernel is invoked in ``DustLayer::advance`` after threshold recomputation,
filling the ``dust_emission_flux`` MultiFab (one component per size bin) on
the 2D dust grid.

Phase 6 References
~~~~~~~~~~~~~~~~~~~

  Marticorena, B., & Bergametti, G. (1995). Modeling the atmospheric
   dust cycle: 1. Design of a soil-derived dust emission scheme.
   *J. Geophys. Res.*, 100, 16415-16430.
   https://doi.org/10.1029/95JD00690

  Owen, P. R. (1964). Saltation of uniform grains in air.
   *J. Fluid Mech.*, 20, 225-242.
   https://doi.org/10.1017/S0022112064001173

  Shao, Y., & Lu, H. (2000). A simple expression for wind erosion
   threshold friction velocity. *J. Geophys. Res.*, 105, 22437-22443.
   https://doi.org/10.1029/2000JD900304

  Bagnold, R. A. (1941). *The Physics of Blown Sand and Desert Dunes*.
   Methuen, London.

Blasting Schedule Emission (Phase 7)
------------------------------------

Mine blasting events generate sudden high-mass dust injection that
cannot be captured by a steady-state saltation model. Phase 7 implements
a timed blast schedule that injects dust mass into ``dust_emission_flux``
at specified times and locations, following the same event dispatch
pattern as ``ERF_IgnitionSchedule.H`` in ``Source/Fire/``.

CSV File Format
~~~~~~~~~~~~~~~

The blast schedule CSV has one event per row (lines starting with ``#``
or ``!`` are comments):

.. code-block:: text

  time_s  cx  cy  radius  mass_kg_m2  [mineral_type]  [priority]

where:

- ``time_s`` — simulation time of blast [s]
- ``cx``, ``cy`` — blast center in physical coordinates [m]
- ``radius`` — affected radius [m]
- ``mass_kg_m2`` — dust mass per unit area injected [kg/m²]
- ``mineral_type`` — 0=quartz tailings (default), 1=Li brine,
  2=REE tailings, 3=Cu tailings
- ``priority`` — higher integer fires first within the same timestep

Injection Formula
~~~~~~~~~~~~~~~~~

For each dust grid cell within ``radius`` of the blast center:

.. math::

  F_\text{blast} = \frac{m \cdot r_b}{dt \cdot N_\text{bins}}

where :math:`m` = ``mass_kg_m2``, :math:`r_b` = ``blast_reactivity``
(default 2.0, representing higher erodibility of fresh blasted surfaces),
:math:`dt` = timestep duration [s], :math:`N_\text{bins}` = number of
size bins. The result is added to ``dust_emission_flux`` and clamped to
``DustBlastConst::BLAST_FLUX_MAX`` = 10⁻² kg/m²/s.

Implementation
~~~~~~~~~~~~~~

``ERF_DustBlastSchedule.H`` is a header-only module that adapts
``ERF_IgnitionSchedule.H``: rank-0 reads the CSV, broadcasts via
``MPI_Bcast``, and applies events per timestep via a host-side event
loop + GPU ``ParallelFor`` per event.

ParmParse Parameters
~~~~~~~~~~~~~~~~~~~~

.. list-table::
  :header-rows: 1
  :widths: 40 60

  * - Parameter
    - Description
  * - ``erf.dust.blast_schedule_file``
    - Path to blast schedule CSV. Empty = no blast events.
  * - ``erf.dust.blast_reactivity``
    - Reactivity multiplier for fresh blast surface [-]. Default = 2.0.

Phase 7 References
~~~~~~~~~~~~~~~~~~

  U.S. EPA AP-42, Chapter 13.2.3 (Heavy Construction Operations).
   https://www.epa.gov/air-emissions-factors-and-quantification/ap-42-compilation-air-emission-factors

  MSHA 30 CFR Part 70 (Worker dust exposure limits).
   https://www.ecfr.gov/current/title-30/chapter-I/subchapter-O/part-70

Suppression Agent Degradation (Phase 8)
---------------------------------------

Mine haul roads and tailings surfaces are treated with dust suppression
agents (water, MgCl₂, lignin sulfonate) to reduce PM₁₀ and PM₂.₅ emissions.
Phase 8 tracks the time-decay of suppression agent coverage and produces
a per-cell re-treatment flag when coverage drops below a threshold.

Decay Model
~~~~~~~~~~~

Suppression coverage :math:`C` [0,1] decays each timestep as:

.. math::

   C_\text{new} = C_\text{old} \cdot \exp\!\left(-\frac{\Delta t}{\tau_\text{eff}}\right)

where the effective half-life :math:`\tau_\text{eff}` [s] is:

.. math::

   \tau_\text{eff} = \frac{\tau_\text{base}}{f_T \cdot f_\text{wind}}

with:

.. math::

   f_T    &= \exp(k_T (T_\text{surf} - T_\text{ref})) \\
   f_\text{wind} &= 1 + k_\text{wind} \cdot u_{10}

Constants: :math:`\tau_\text{base}` = ``erf.dust.supp_tau_base_s`` (default
3600 s for water), :math:`T_\text{ref}` = 293.15 K, :math:`k_T` = 0.05 K⁻¹,
:math:`k_\text{wind}` = 0.1 s/m.

Coverage values below ``DustSuppressionConst::COVERAGE_FLOOR`` = 0.01 are
set to zero. Cells where coverage drops below
``DustSuppressionConst::RETREAT_THRESHOLD`` = 0.2 receive a re-treatment
flag value of 1.0 in ``dust_retreat_flag``.

Structural Reference
~~~~~~~~~~~~~~~~~~~~

``ERF_DustSuppression.H`` follows the pattern of ``ERF_FireBreak.H`` in
``Source/Fire/``: a local 2D property modification on the surface grid via
a GPU ``ParallelFor`` kernel.

References:
  Gillies, J.A., et al. (2005). J. Air Waste Manage. Assoc., 55, 1168.
  https://doi.org/10.1080/10473289.2005.10464701

  U.S. EPA AP-42, Chapter 13.2.2 (Unpaved Roads).
  https://www.epa.gov/air-emissions-factors-and-quantification/ap-42-compilation-air-emission-factors

  MSHA 30 CFR Part 70. https://www.ecfr.gov/current/title-30/chapter-I/subchapter-O/part-70

ParmParse Parameters
~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Parameter
     - Description
   * - ``erf.dust.supp_tau_base_s``
     - Base suppression half-life [s]. Default = 3600 (water agent).
   * - ``erf.dust.test_surf_temp_K``
     - Surface temperature placeholder [K]. Phase 9 uses MRF T_sfc.
   * - ``erf.dust.test_wind_speed``
     - 10 m wind speed placeholder [m/s]. Phase 9 uses MRF wind.

Development Phases
------------------

Development is divided into 24 phases. Phase completion status is recorded
in ``Source/Dust/DUST_DEVELOPMENT.md``. Source code comments are restricted
to technical descriptions of algorithms and data structures.

Phases 1-8 construct the 2D grid infrastructure and PHREEQC file coupling.
Phases 9-16 implement the two-way coupling with the ERF 3D solver.
Phases 17-24 produce regulatory output aligned with EPA NAAQS, MSHA, and
DOE Critical Materials Assessment requirements.

References
----------

.. [Bagnold1941] Bagnold, R. A. (1941). *The Physics of Blown Sand and
   Desert Dunes*. Methuen, London.

.. [Hong1996] Hong, S.-Y., and H.-L. Pan (1996). Nonlocal boundary layer
   vertical diffusion in a medium-range forecast model.
   *Mon. Wea. Rev.*, 124, 2322-2339.
   https://doi.org/10.1175/1520-0493(1996)124<2322:NBLVDI>2.0.CO;2

.. [Marticorena1995] Marticorena, B., and G. Bergametti (1995). Modeling
   the atmospheric dust cycle: 1. Design of a soil-derived dust emission
   scheme. *J. Geophys. Res.*, 100, 16415-16430.
   https://doi.org/10.1029/95JD00690

.. [Shao2000] Shao, Y., and H. Lu (2000). A simple expression for wind
   erosion threshold friction velocity. *J. Geophys. Res.*, 105,
   22437-22443.
   https://doi.org/10.1029/2000JD900304

.. [Parkhurst2013] Parkhurst, D. L., and C. A. J. Appelo (2013).
   Description of input and examples for PHREEQC version 3.
   *U.S. Geological Survey Techniques and Methods*, book 6, chap. A43.
   https://pubs.usgs.gov/tm/06/a43/

.. [EPA_AP42] U.S. EPA (1998). AP-42, Compilation of Air Emission Factors,
   Chapter 13.2.2: Unpaved Roads.
   https://www.epa.gov/air-emissions-factors-and-quantification/ap-42-compilation-air-emission-factors

.. [MSHA2016] MSHA (2016). Lowering Miners' Exposure to Respirable Coal
   Mine Dust, Including Continuous Personal Dust Monitors.
   *Federal Register*, 81, 16285.
   https://www.federalregister.gov/documents/2014/05/01/2014-09084/

.. [EO14156] Executive Order 14156 (2025). Declaring a National Energy
   Emergency. *Federal Register*, 90, 8433.

.. [EO14157] Executive Order 14157 (2025). Unleashing American Energy.
   *Federal Register*, 90, 8437.

.. [DOE_CMA] U.S. Department of Energy (2023). Critical Materials
   Assessment. Office of Energy Efficiency and Renewable Energy.
   https://www.energy.gov/eere/vehicles/articles/doe-critical-materials-assessment
