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

Wind and Surface Field Extraction (Phase 9)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Phase 9 extracts atmospheric wind and surface fields from the ERF 3D solver
to the 2D dust grid each timestep. The MRF PBL parameterization provides:

- Friction velocity :math:`u^* = \sqrt{\tau_0 / \rho_a}` [m/s]
- Surface temperature :math:`T_{\text{surf}}` [K]
- PBL height :math:`H_{\text{PBL}}` [m]

These are interpolated from level 0 (first vertical cell) of the 3D grid
onto the dust grid using inverse-distance weighting on the horizontal mesh.
Wind speed at reference height :math:`z_{\text{ref}}` (default 10 m) is
extracted via cubic Hermite interpolation along the vertical to support
friction velocity diagnostics and emissions testing.

Structural Reference
  ``ERF_DustWindExtract.H/cpp`` adapts ``fill_fire_wind_from_interpolation``
  from ``Source/Fire/``. Uses ``SurfaceLayer`` accessors:
  ``get_u_star(0)``, ``get_t_surf(0)``, ``get_pblh(0)``.

Aerosol Injection into ERF-Atm (Phase 10)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Phase 10 implements one-way aerosol injection: the 2D dust emission flux
(computed in Phases 5-9) is coarsened to the 3D atmospheric grid and injected
as a mass source into a passive scalar slot.

Scalar Slot Design
  Phase 10 uses **one passive scalar slot** (``RhoAdv_comp``, index 4 in the
  conserved state) to carry total dust mass density [kg/m³]. Per-bin composition
  is tracked on the 2D dust grid and is not needed in 3D until Phase 11.
  This minimizes advection overhead. A new parameter ``transport_bins_separately``
  (default ``false``) controls scalar expansion in Phase 11.

Injection Formula
  At the surface (k=0), dust mass receives a source tendency:

  .. math::

    \frac{\partial (\rho_\text{dust})}{\partial t}\bigg|_{k=0}
      = \frac{F_\text{dust}}{\Delta z_{k=0}}

  where :math:`F_\text{dust}` [kg/m²/s] is the total emission flux coarsened
  from the dust grid to the atmosphere grid, and :math:`\Delta z_{k=0}` [m]
  is the depth of the first atmospheric cell. All levels :math:`k > 0` receive
  zero tendency from this injection; transport away from the surface is handled
  by advection and diffusion in the dycore.

Coupling Strength
  The parameter ``atm_feedback`` [0,1] scales the injection strength (default 1.0).
  Setting ``atm_feedback = 0`` disables injection for diagnostic tests.

Temporal Lag
  The injection uses a one-step explicit lag: emission flux from time step
  :math:`n` is injected at the start of step :math:`n+1`, before
  ``advance_dycore()``. This is consistent with the fire module pattern
  in ``ERF_FireAtmCoupling.H``.

Coarsening
  If the dust grid is refined (``grid_ratio > 1``), the flux is coarsened
  to the atmospheric mesh via ``amrex::average_down``. For ``grid_ratio = 1``,
  coarsening is a direct copy.

Structural Reference
  ``ERF_DustAtmCoupling.H/cpp`` adapts ``coarsen_fire_flux_to_atm`` and
  ``apply_fire_tendency_to_cc_source`` from ``Source/Fire/``. The header
  file documents the WRF-Fire coupling approach; dust injection uses the
  same pattern but does not apply a vertical profile (injection is
  surface-only).

References
  Mandel, J., et al. (2011). An adaptive time stepping algorithm for
  coupled atmosphere–wildfire simulation.
  *Math. Comput. Simul.*, 79 (3), 584–606.
  https://doi.org/10.1016/j.matcom.2008.03.015

Stokes Settling Velocity (Phase 11)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Phase 11 implements vertical settling of dust particles via first-order upwind
advection, applying the Stokes law for particle dynamics. Particles settle out
of the domain at the lower boundary (deposition is added in Phase 12).

Stokes Law Settling Velocity
  For a spherical particle in Stokes flow (valid for diameters < 100 µm):

  .. math::

    v_s = \frac{(\rho_p - \rho_a) g d^2}{18 \mu_a}

  where :math:`\rho_p` [kg/m³] is particle density, :math:`\rho_a` [kg/m³]
  is local air density, :math:`g` = 9.81 m/s² is gravitational acceleration,
  :math:`d` [m] is particle diameter, and :math:`\mu_a` [Pa·s] is dynamic
  viscosity of air.

Cunningham Slip Correction
  For sub-micron particles, molecular slip must be accounted for:

  .. math::

    C_c = 1 + \frac{2\lambda}{d} \left( A_1 + A_2 \exp\left(-\frac{A_3 d}{\lambda}\right) \right)

  where :math:`\lambda` = 0.066 µm is the mean free path at STP,
  :math:`A_1` = 1.257, :math:`A_2` = 0.400, :math:`A_3` = 0.55
  (Allen & Raabe 1985). For diameters > 1 µm the correction is < 15%.
  The corrected settling velocity is :math:`v_s^{\text{corr}} = v_s C_c`,
  clamped to [0.0, 1.0] m/s.

Vertical Tendency (First-Order Upwind Flux Divergence)
  The settling tendency is applied as a divergence of the downward mass flux:

  .. math::

    \frac{\partial (\rho_\text{dust})}{\partial t}\bigg|_{\text{settling}}
      = -\frac{1}{\Delta z} \left( F_\text{down}(k) - F_\text{down}(k-1) \right)

  where the downward flux at the top face of cell :math:`k` is:

  .. math::

    F_\text{down}(k) = v_s(k) \cdot \rho_\text{dust}(k)

  and :math:`v_s(k)` is evaluated using local air density :math:`\rho_a(k)`
  from the previous timestep. At the lower boundary (k=0), the flux leaving
  the domain is held at zero in Phase 11 (Phase 12 adds dry deposition).

Multi-Bin Scalar Expansion
  When ``transport_bins_separately = false`` (default, lower cost),
  all particles are advected as a single 3D scalar using the diameter
  of the first bin. When ``transport_bins_separately = true``, each
  particle size bin receives its own 3D scalar slot, allowing accurate
  representation of size-dependent settling and transport.

  The parameter ``bin_diameters`` [m] specifies the representative diameter
  for each bin (default: 7.0 µm, 2.5 µm, 50.0 µm for PM10, PM2.5, and
  coarse mineral dust respectively). Phase 10 injection is expanded to
  populate all active bins in Phase 11 when this flag is enabled.

Structural Reference
  ``ERF_DustSettling.H`` contains device functions for Cunningham correction
  and settling velocity computation, plus an inline MultiFab kernel that
  applies the first-order upwind tendency. It follows the pattern of
  ``ERF_FireAtmCoupling.H`` in adapting ``apply_fire_tendency_to_cc_source``.

References
  Stokes, G. G. (1851). On the effect of the internal friction of fluids
  on the motion of pendulums. *Trans. Cambridge Phil. Soc.*, 9, 8–106.

  Allen, M. D., and O. G. Raabe (1985). Slip correction measurements of
  spherical aerosol particles at known Stokes numbers.
  *Aerosol Sci. Tech.*, 4, 269–286.
  https://doi.org/10.1080/02786828508959055

  Seinfeld, J. H., and S. N. Pandis (2006).
  *Atmospheric Chemistry and Physics: From Air Pollution to Climate Change*,
  2nd ed., Chapter 9. Wiley.
  ISBN 978-0-471-72018-8.

  Ginoux, P., M. Chin, I. Tegen, J. Prospero, B. Holben, O. Dubovik,
  and S.-J. Lin (2001). Sources and distributions of dust aerosols
  simulated with the GOCART model. *J. Geophys. Res.*, 106, 20,255–20,273.
  https://doi.org/10.1029/2000JD901323

Dry Deposition Lower Boundary Condition (Phase 12)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Phase 12 implements physically correct dry deposition at the lower boundary,
replacing the zero-flux condition in Phase 11. Deposited mass is accumulated
on a 2D dust grid for Phase 18 (MSHA worker exposure) and Phase 21 (PHREEQC
feedback file output).

Zhang et al. (2001) Resistance Scheme
  Dry deposition velocity combines gravitational settling and surface
  collection efficiency:

  .. math::

    v_d = v_s + \frac{1}{r_a + r_s + r_a r_s v_s}

  where :math:`v_s` [m/s] is the Stokes settling velocity (Phase 11),
  :math:`r_a` [s/m] is aerodynamic resistance, and :math:`r_s` [s/m] is
  surface resistance.

Aerodynamic and Surface Resistances
  The aerodynamic resistance depends on friction velocity and von Kármán
  constant:

  .. math::

    r_a = \frac{1}{u_* \kappa}

  where :math:`u_*` [m/s] is friction velocity (extracted from SurfaceLayer
  in Phase 9) and :math:`\kappa` = 0.4. The surface resistance accounts
  for collection efficiency:

  .. math::

    r_s = \frac{1}{u_* E_0}

  where :math:`E_0` [-] is the collection efficiency. Values from Zhang et al.
  (2001):

  * Bare mine surface: :math:`E_0` = 3.0×10⁻³ (default, ``deposition_E0`` parameter)
  * Paved road: :math:`E_0` = 1.0×10⁻⁴
  * Vegetated buffer: :math:`E_0` = 1.0×10⁻²

Lower Boundary Flux and Accumulator
  At the lowest grid cell (k=0), the deposition flux is:

  .. math::

    F_{\text{dep}} = v_d \cdot \rho_{\text{dust}}(k=0) \quad [\text{kg/m}^2/\text{s}]

  This flux is applied as a source term:

  .. math::

    \frac{\partial (\rho_\text{dust})}{\partial t}\bigg|_{k=0,\,\text{dep}}
      = -\frac{v_d \cdot \rho_\text{dust}(k=0)}{\Delta z_{k=0}}

  Each timestep, the deposition flux is accumulated on the 2D dust grid
  in ``dust_deposition_rate`` [kg/m²], initialized to zero and never reset:

  .. math::

    M_{\text{dep}}(t) = M_{\text{dep}}(t-1) + v_d(t) \cdot \rho_\text{dust}(k=0, t) \cdot \Delta t

  This running total is used by Phase 18 to track cumulative exposure and
  by Phase 21 to compute geochemical feedback.

Test Configuration Requirements
  All Phase 9–12 tests use the Moeng-Sullivan neutral ABL with:

  * Domain: 3000 m × 3000 m × 1024 m
  * Grid: 8 × 8 × 64 cells
  * Geostrophic wind: 15.0 m/s → u* ≈ 0.56 m/s from MRF
  * Reference height: ``erf.most.zref`` = ``erf.dust.zref`` = 24.0 m
  * ``erf.transport_scalar = true`` required for scalar advection (PR #145 fix)
  * Sounding file: ``sounding_neutral_abl`` copied to each test directory

Structural Reference
  ``ERF_DustDeposition.H`` contains device functions for aerodynamic and
  surface resistance computation, plus inline MultiFab kernels for the
  deposition tendency and accumulator update. The module follows the pattern
  of ``ERF_FireAtmCoupling.H`` and ``ERF_DustSettling.H``.

References
  Seinfeld, J. H., and S. N. Pandis (2006).
  *Atmospheric Chemistry and Physics: From Air Pollution to Climate Change*,
  2nd ed., Chapter 19. Wiley.
  ISBN 978-0-471-72018-8.

  Zhang, L., R. Gong, P. A. Dawson, and C. P. J. Gong (2001). A parametrization
  of dry particle deposition velocity for modeling aerosol transport and
  deposition. *Atmos. Environ.*, 35, 549–560.
  https://doi.org/10.1016/S1352-2310(00)00326-5

  Slinn, W. G. N. (1982). Predictions for particle deposition to vegetative
  canopies. *Atmos. Environ.*, 16, 1785–1794.
  https://doi.org/10.1016/0004-6981(82)90271-2

Two-Way Coupling: Return Fields from ERF-Atm to ERF-Dust (Phase 13)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Phase 13 implements the return path, mapping three diagnostic and feedback
fields from the 3D atmospheric conserved state back to the 2D dust surface
grid. These fields enable loading feedback on the threshold and dynamic
moisture inhibition, with the deposition accumulator confirming Phase 12
correctness.

Near-Surface Dust Concentration
  The dust mass density at the lowest atmospheric cell (k=klo) is extracted
  from the conserved state and mapped to the 2D dust grid as ``dust_conc_sfc``
  [kg/m³]:

  .. math::

   C_{\text{sfc}}(i,j) = \max(\rho_\text{dust}(i/C, j/C, k=\text{klo}), 0)

  where :math:`C` is the grid refinement ratio (dimensionless). This field
  feeds the Shao (2001) loading feedback correction to the threshold friction
  velocity.

Shao (2001) Loading Feedback
  When ``loading_feedback_coeff`` > 0, the effective threshold is increased by
  near-surface dust concentration:

  .. math::

   u_{*,t}^{\text{eff}} = u_{*,t} \cdot (1 + \alpha_L \cdot C_{\text{sfc}})

  where :math:`\alpha_L` = ``loading_feedback_coeff`` [m³/kg]. This feedback
  models the suppression of saltation when sediment concentration builds up
  near the surface, reducing the potential for further emission. Setting
  :math:`\alpha_L = 0` (default) disables the feedback. Typical values from
  observations are in the range 10³–10⁷ m³/kg.

  Reference: Shao, Y. (2001). *J. Geophys. Res.*, 106, 20599–20610.
  https://doi.org/10.1029/2001JD900171

Surface Moisture Flux
  The vertical moisture flux at the bottom z-face (k=klo, face-staggered z)
  is extracted from ``Q1fx3`` and mapped to the 2D dust grid as
  ``dust_surf_moist`` [kg/m²/s]:

  .. math::

   q_{\text{flux}}(i,j) = Q1fx3(i/C, j/C, k=\text{klo})

  This field is null-safe: when the atmospheric model does not use a moisture
  scheme (``moisture_type == None``, all current dust tests), ``Q1fx3`` is null
  and ``dust_surf_moist`` is zeroed. The Phase 5 static moisture flag is used
  unchanged in this case. When a moisture scheme is active and
  ``use_dynamic_moisture = true``, the flux is used to compute the Fecan et al.
  (1999) dynamic moisture inhibition factor.

Fecan (1999) Dynamic Moisture Inhibition
  When ``use_dynamic_moisture = true`` and a moisture scheme is active,
  gravimetric water content :math:`w` [kg/kg] is derived from the surface flux:

  .. math::

   w = \frac{q_{\text{flux}}}{L_v \rho_a}

  where :math:`L_v = 2.501 \times 10^6` J/kg is the latent heat of vaporization
  and :math:`\rho_a = 1.225` kg/m³ is reference air density. The Fecan inhibition
  factor is:

  .. math::

   f_{\text{moist}} = \sqrt{1 + a_f \max(w - w', 0)}

  where :math:`a_f = 1.21` is the Fecan coefficient and :math:`w' = 0.003` kg/kg
  is the residual moisture threshold. The threshold is then:

  .. math::

   u_{*,t}^{\text{final}} = u_{*,t}^{\text{eff}} \cdot f_{\text{moist}}

  When ``Q1fx3`` is null (no moisture scheme), :math:`w = 0`, so
  :math:`f_{\text{moist}} = 1.0` and no inhibition is applied. This is correct
  for dry conditions.

  Reference: Fecan, F., B. Marticorena, and G. Bergametti (1999). Parametrization
  of the increase of the aeolian erosion threshold wind friction velocity due to
  soil moisture for arid and semi-arid areas. *Ann. Geophys.*, 17, 149–157.
  https://doi.org/10.1007/s00585-999-0149-7

Deposition Accumulator Diagnostic
  Phase 13 adds a domain-integral print each timestep confirming that the Phase 12
  ``dust_deposition_rate`` accumulator is active. The print appears in the
  ``[DUST DEBUG] Phase 13`` output as ``dep_total`` [kg/m²].

Coarsening Pattern
  All three fields use the same coarsening map as Phase 9 wind extraction and
  Phase 10 emission injection:

  .. code-block:: cpp

   dust_field(i,j,0) = atm_field(i/C, j/C, k_ref)

  This pattern ensures consistency between forward (emission injection) and
  backward (concentration/flux extraction) coupling on the two grids.

Implementation
  The extraction functions are header-only in ``ERF_DustAtmReturn.H``:

  * ``fill_dust_conc_from_atm``: Extract RhoDust at k=klo from conserved state.
  * ``fill_dust_moist_from_atm``: Extract Q1fx3 at z-face klo. Null-safe.

  Called from ``ERF_TI_slow_rhs_post.H`` after ``erf_slow_rhs_post`` completes,
  ensuring that S_new conserved state and Q1fx3 flux are updated before extraction.
  The Phase 5 threshold recomputation then applies the loading and moisture
  feedbacks before Phase 6 computes the emission flux for the next timestep.

Parameters
  * ``loading_feedback_coeff`` [m³/kg] (default: 0.0): Shao loading feedback
   coefficient. 0 = disabled. Positive values enable feedback.
  * ``use_dynamic_moisture`` [boolean] (default: false): Enable Fecan dynamic
   moisture inhibition. Safe to set true with no moisture scheme: Q1fx3 will
   be null and function returns immediately.

Test Configuration
  The ``DustAtmReturn`` test uses the neutral ABL configuration and verifies:

  * ``conc_sfc_max > 0`` after step 1 (dust injected by Phase 10 appears in 3D scalar)
  * ``conc_sfc_sum`` changes with step count
  * ``moisture_path_active = 0`` (Q1fx3 null, dry neutral ABL)
  * ``moist_flux_max = 0`` (no moisture flux in dry case)
  * ``dep_total`` increases each step (Phase 12 accumulator confirmed active)
  * Run exits with code 0 at ``max_step = 5``

  To test loading feedback: set ``loading_feedback_coeff = 1e6`` and verify that
  ``emission_flux_max`` decreases as ``conc_sfc_max`` grows over successive steps.

  To test the moisture path: enable a moisture scheme (``erf.moisture_type != None``)
  and set ``use_dynamic_moisture = true``. Q1fx3 will be non-null,
  ``moisture_path_active = 1`` will appear in debug output, and the Fecan
  inhibition factor will reduce u*_t in cells with surface moisture flux.

MRF Nonlocal Countergradient Extension to Dust Scalar (Phase 14)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Phase 14 extends the Hong & Pan (1996) MRF boundary layer scheme's nonlocal
countergradient transport framework to the dust scalar component, ensuring that
dust particles are diffused vertically with the same eddy diffusivity profile
:math:`K_h(z)` used for heat and moisture. Since dust has no prescribed surface
flux gradient (unlike heat and moisture), the countergradient term for dust
vanishes.

The MRF Nonlocal Diffusion Profile
  In the Hong & Pan (1996) MRF scheme, the vertical eddy viscosity in the PBL
  is computed as:

  .. math::

    K_h(z) = w_* \kappa h \frac{z}{h} \left(1 - \frac{z}{h}\right)^2

  where:

  * :math:`w_*` = convective velocity scale [m/s] (related to friction velocity :math:`u_*`)
  * :math:`\kappa = 0.4` = von Kármán constant
  * :math:`h` = PBL height [m]
  * :math:`z` = height above ground [m]

  This profile peaks at :math:`z/h = 1/3`, where
  :math:`K_{h,\text{peak}} \approx 0.037 \, w_* \, h`. For neutral boundary layers
  with u* ~ 0.56 m/s and h ~ 500 m (typical MRF diagnostic), the peak diffusivity
  is ~27.7 m²/s.

Scalar Diffusivity for Dust
  Phase 14 sets ``EddyDiff::Scalar_v`` (vertical scalar diffusivity) equal to
  ``EddyDiff::Theta_v`` (vertical heat diffusivity) in ``ComputeDiffusivityMRF.cpp``.
  The dust scalar component is then automatically diffused using this profile when
  ``erf.transport_scalar = true``, via the existing ``DiffusionSrcForState_N/T/S``
  routines which reference ``EddyDiff::Scalar_v`` for any passive scalar in the
  ``RhoScalar_comp`` range.

Zero Countergradient for Dust
  The MRF scheme includes countergradient corrections for heat
  (``HGAMT_v``, Hong & Pan 1996 Eq. 3) and moisture (``HGAMQ_v``):

  .. math::

    K_h \frac{\partial \theta}{\partial z} = K_h \left(\frac{\partial \overline{\theta}}{\partial z} - \gamma_\theta\right)

  where :math:`\gamma_\theta` represents the nonlocal countergradient transport.
  For dust, :math:`\gamma_{\text{dust}} = 0` because there is no prescribed surface
  dust flux analogous to the sensible heat flux :math:`\overline{w''\theta''}_{\text{sfc}}`.
  The diffusion is therefore purely gradient-driven, consistent with passive
  aerosol tracer transport in the literature.

Turbulent Schmidt Number Scaling
  By default, dust diffusivity is identical to heat diffusivity, corresponding to
  a turbulent Schmidt number :math:`Sc_t` equal to the turbulent Prandtl number
  :math:`Pr_t`. This is the standard assumption for passive aerosol tracers at large
  Reynolds numbers (Seinfeld & Pandis 2006, Ch. 16).

  For fine-tuning, Phase 14 introduces the parameter ``erf.dust.mrf_Sc_t``
  [dimensionless]:

  .. math::

    K_h^{\text{dust}} = K_h^{\text{heat}} \times \frac{Pr_t}{Sc_t}

  * Default: ``Sc_t = 0`` or ``Sc_t = Pr_t`` → no scaling (dust diffusivity = heat diffusivity)
  * Example: ``Sc_t = 2.0`` → dust diffusivity = 0.5 × heat diffusivity (halved relative to heat)

  Reference: Seinfeld, J.H., and S.N. Pandis (2006). *Atmospheric Chemistry and
  Physics: From Air Pollution to Climate Change*, 2nd ed., Ch. 16. Wiley.
  ISBN 978-0-471-72018-8.

Diagnostics
  Phase 14 adds two debug outputs (when ``erf.dust.dust_debug = true``):

  1. In ``extract_atm_return_fields``: Confirmation that gamma_dust = 0 (no countergradient)
  2. In ``ComputeTurbulentViscosity``: Print of :math:`K_{h,\text{max}}^{\text{Scalar}}` and
    :math:`K_{h,\text{max}}^{\text{Heat}}` for verification of Sc_t scaling

Parameters
  * ``erf.dust.mrf_Sc_t`` [dimensionless] (default: 0 → use Pr_t): Turbulent Schmidt
   number for dust in MRF PBL scheme. Controls ``EddyDiff::Scalar_v`` relative to
   ``EddyDiff::Theta_v``.

Test Configuration
  The ``DustMRFDiffusion`` test uses the neutral ABL configuration and verifies:

  * Phase 14 debug output shows ``Scalar_v_max`` and ``Theta_v_max`` each step
  * With ``Sc_t = 0`` (default): ``Scalar_v_max ≈ Theta_v_max`` (within 1%)
  * With ``Sc_t = 2.0``: ``Scalar_v_max ≈ 0.5 × Theta_v_max``
  * Dust concentration profile exhibits vertical diffusion (reduced gradient)
  * Phase 13 and 14 debug prints confirm MRF diffusion is active
  * Run exits with code 0 at ``max_step = 5``

Implementation
  Phase 14 modifications are concentrated in two files:

  * ``Source/PBL/ERF_ComputeDiffusivityMRF.cpp``: Set ``EddyDiff::Scalar_v`` and
   apply Sc_t scaling
  * ``Source/DataStructs/ERF_TurbStruct.H``: Add ``dust_mrf_Sc_t`` parameter and
   parser

  No diffusion calling code changes are required; ``DiffusionSrcForState_N/T/S``
  routines already reference ``EddyDiff::Scalar_v`` for dust via the existing
  ``is_valid_slow_var`` component index logic.

References
  .. [Hong1996MRF] (already defined above)
  .. [Seinfeld2006] Seinfeld, J.H., and S.N. Pandis (2006). *Atmospheric Chemistry
    and Physics: From Air Pollution to Climate Change*, 2nd ed., Ch. 16. Wiley.
    ISBN 978-0-471-72018-8.
  .. [Ginoux2001] Ginoux, P., M. Chin, I. N. Tegen, J. M. Prospero, B. Holben,
    O. Dubovik, and S.-J. Lin (2001). Sources and distributions of dust aerosols
    simulated with the GOCART model. *J. Geophys. Res.*, 106, 20255–20273.
    https://doi.org/10.1029/2000JD901323

Dust Diagnostics and Plotfile Output (Phase 16)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Phase 16 implements diagnostic output for the dust module, writing two types of files:

1. **Dust Plotfile** on native grid: 8-component 2D AMReX plotfile on the dust grid
   (independent of ERF AMR hierarchy)
2. **CSV time series**: Per-step scalar diagnostics to ``dust_diag.dat``

Plotfile Format
  The dust plotfile is written using the `VisMF::Write` + `WriteGenericPlotfileHeader`
  pattern, following the fire module (``ERF_FirePlotfile.cpp``). This approach is
  necessary because the dust grid has its own ``BoxArray`` and ``Geometry``
  (``m_dg.ba``, ``m_dg.geom``) that are independent of the ERF AMR level hierarchy.
  The standard ``WriteMultiLevelPlotfile`` cannot be used for this reason.

  **Output structure**:

  * ``plotfile_name/Header``: Standard AMReX ASCII header with variable names, geometry,
    and time metadata
  * ``plotfile_name/Level_0/Cell_H``, ``Cell_000000``, etc.: Binary data from ``VisMF::Write``
  * ``plotfile_name/DustMetadata.json``: Sidecar JSON with ``format_version``, ``time``,
    ``step``, ``grid_ratio``, ``n_variables``

  **Variables written per plotfile** (8 components):

  1. ``dust_emission_flux`` [kg/m²/s] — Total saltation + vertical emission flux.
     Marticorena & Bergametti (1995).
  2. ``dust_ustar_in`` [m/s] — Friction velocity u* extracted from MRF surface layer.
  3. ``dust_ustar_t`` [m/s] — Threshold friction velocity. Shao & Lu (2000).
  4. ``dust_deposition_rate`` [kg/m²] — Accumulated dry deposition. Zhang et al. (2001).
  5. ``dust_conc_sfc`` [kg/m³] — Near-surface dust concentration (loading feedback).
     Shao (2001).
  6. ``dust_surf_moist`` [kg/m²/s] — Surface moisture flux. Null-safe when no moisture
     scheme. Fecan et al. (1999).
  7. ``dust_suppression`` [0–1] — Suppression agent coverage. Gillies et al. (2005).
  8. ``dust_retreat_flag`` [0/1] — Re-treatment flag.

CSV Time Series
  The ``dust_diag.dat`` file follows the fire module CSV append pattern
  (``ERF_FireStatsOutput.H``): header written on first call, append one row per step.
  Only the IOProcessor (rank 0) writes; all ranks compute reductions in parallel.

  **Column definitions**:

  * ``step`` — Current timestep number
  * ``time_s`` — Simulation time [s]
  * ``emission_total_kg_s`` — Domain-integrated emission flux [kg/s]
  * ``deposition_total_kg_m2`` — Domain-accumulated deposition [kg/m²]
  * ``ustar_max_m_s`` — Maximum u* [m/s] (used for threshold checks)
  * ``flux_max_kg_m2_s`` — Maximum emission flux [kg/m²/s]
  * ``conc_sfc_max_kg_m3`` — Maximum near-surface concentration [kg/m³]

Write Timing and Final-Step Guard
  Plotfiles and diagnostics are written at:

  * Every ``dust_plot_int`` steps (if ``dust_plot_int > 0``)
  * At the final simulation step, even if the step number is not a multiple of
    ``dust_plot_int``

  The final-step guard uses ``m_last_dust_plot_step``, initialized to -1 and updated
  each time a plotfile is written. In ``WriteAtFinalTime``, the condition
  ``is_final=true && nstep > m_last_dust_plot_step`` ensures a write at simulation
  end regardless of ``dust_plot_int``.

Parameters
  * ``erf.dust.dust_plot_int`` [steps] (default: -1) — Interval between plotfile writes.
    Disabled when ``-1``. Set to 0 to write only at final time.
  * ``erf.dust.dust_plot_prefix`` [string] (default: ``"plt_dust_"``) — Prefix for
    plotfile directory names (format: ``{prefix}{step:05d}``)
  * ``erf.dust.dust_diag_file`` [string] (default: ``"dust_diag.dat"``) — Path to CSV
    time series file (relative to run directory)

Test Configuration
  The ``DustOutput`` regtest uses ``max_step = 7`` and ``dust_plot_int = 3``:

  * Plotfiles created at steps 0 (InitData_post), 3, 6 (multiples of 3), and 7 (final)
  * Each plotfile contains Header, Level_0/Cell* binary, and DustMetadata.json
  * CSV file has header row plus 8 data rows (steps 0–7)
  * Emissions increase from 0 at step 0 to max by step 1
  * Deposition increases monotonically with time
  * Phase 16 debug output confirms plotfile writes at correct steps

Disabling Plotfile Output
  Set ``erf.dust.dust_plot_int = -1`` to disable plotfile creation (default state).
  CSV diagnostics continue to be written every step regardless of ``dust_plot_int``.

Implementation
  Phase 16 introduces three new files:

  * ``Source/Dust/ERF_DustPlotfileCatalog.H`` — Header-only; provides
    ``dust_plotfile_var_names()`` and ``dust_plotfile_ncomp()``
  * ``Source/Dust/ERF_DustPlotfile.H`` — Forward declaration for ``WriteDustPlotfile``
  * ``Source/Dust/ERF_DustPlotfile.cpp`` — Implements collective VisMF::Write and
    IOProcessor-only Header/JSON writes
  * ``Source/Dust/ERF_DustStatsOutput.H`` — Header-only; provides
    ``write_dust_stats_header`` and ``append_dust_stats``

  ``DustLayer::write_output`` is called from:

  * ``ERF::InitData_post`` (initial state)
  * ``ERF::WriteAtIntermediateTime`` (every timestep)
  * ``ERF::WriteAtFinalTime`` (at end of run)

References
  .. [Zhang2001] Zhang, L., S. L. Gong, J. Padro, and L. Barrie (2001). A size-segregated
    particle dry deposition scheme for an atmospheric aerosol module.
    *Atmos. Environ.*, 35, 549–560.
    https://doi.org/10.1016/S1352-2310(00)00326-5
  .. [Fecan1999] Fecan, F., B. Marticorena, and G. Bergametti (1999). Parametrization
    of the increase of the aeolian erosion threshold wind friction velocity due to
    soil moisture for arid and semi-arid areas.
    *Ann. Geophys.*, 17, 149–157.
    https://doi.org/10.1007/s00585-999-0149-7
  .. [Gillies2005] Gillies, J. A., W. G. Nickling, and G. S. King (2005). Shear stress
    partitioning in large-scale surface roughness for shortgrass rangeland.
    *J. Air Waste Manag. Assoc.*, 55, 539–549.
    https://doi.org/10.1080/10473289.2005.10464701

EPA NAAQS PM2.5/PM10 Compliance Module
---------------------------------------

Phase 17 adds regulatory compliance diagnostics for EPA National Ambient Air Quality
Standards (NAAQS) particulate matter (PM) standards. The module computes PM2.5 and
PM10 mass concentrations from the dust scalar field and tracks 24-hour running averages
and exceedance flags at each model grid cell.

**Regulatory Background**

The U.S. EPA NAAQS sets national air quality standards for particulate matter at ground
level [EPA2024]_. Phase 17 focuses on 24-hour averages:

* **PM2.5 24-hour**: 35 µg/m³ (98th percentile of 24-hour averages over 3 years). 40 CFR Part 50.18.
* **PM10 24-hour**: 150 µg/m³ (not to be exceeded more than once per year on average over 3 years). 40 CFR Part 50.6.

Mine sites are regulated under 40 CFR Part 60 Subpart OOO (nonmetallic mineral processing)
and Subpart LL (metallic mineral processing) with fugitive dust provisions.

**Size Fractionation**

Dust size bins are classified into PM classes based on aerodynamic diameter following
EPA convention [Watson1994]_:

* **PM2.5**: aerodynamic diameter ≤ 2.5 µm
* **PM10**: aerodynamic diameter ≤ 10 µm (includes PM2.5)

A bin with diameter :math:`d` [m] contributes its full mass to the PM class if
:math:`d` ≤ threshold. No partial fractional attribution is applied within a bin at
this phase.

**Running Time Averages**

At each grid cell, 24-hour arithmetic mean PM concentrations are maintained using an
exponential moving average:

.. math::

  C_{\text{avg}}(t) = C_{\text{avg}}(t - \Delta t) \cdot \frac{T_{\text{avg}} - \Delta t}{T_{\text{avg}}} + C_{\text{now}} \cdot \frac{\Delta t}{T_{\text{avg}}}

where :math:`T_{\text{avg}} = 86400` s (24 hours), :math:`\Delta t` is the timestep, and
:math:`C_{\text{now}}` is the instantaneous concentration.

**Exceedance Flags**

At each grid cell, binary flags [0/1] indicate whether the 24-hour PM average exceeds
the NAAQS threshold:

* ``dust_pm25_exceed`` = 1 where PM2.5 24h > 35 µg/m³
* ``dust_pm10_exceed`` = 1 where PM10 24h > 150 µg/m³

**CSV Output**

Per-step domain statistics are written to ``dust_naaqs.csv`` with columns:

* ``step`` — Step number
* ``time_s`` — Simulation time [s]
* ``PM25_max_ug_m3`` — Instantaneous PM2.5 maximum [µg/m³]
* ``PM25_24h_max_ug_m3`` — 24-hour running average PM2.5 maximum [µg/m³]
* ``PM10_max_ug_m3`` — Instantaneous PM10 maximum [µg/m³]
* ``PM10_24h_max_ug_m3`` — 24-hour running average PM10 maximum [µg/m³]
* ``n_cells_PM25_exceed`` — Number of cells exceeding PM2.5 threshold
* ``n_cells_PM10_exceed`` — Number of cells exceeding PM10 threshold

**Plotfile Integration**

Phase 17 adds 6 new variables to the 2D dust plotfile (total 14 fields):

* ``dust_pm25_ug_m3`` — Instantaneous PM2.5 [µg/m³]
* ``dust_pm10_ug_m3`` — Instantaneous PM10 [µg/m³]
* ``dust_pm25_24h_ug_m3`` — 24-hour running average PM2.5 [µg/m³]
* ``dust_pm10_24h_ug_m3`` — 24-hour running average PM10 [µg/m³]
* ``dust_pm25_exceed`` — PM2.5 exceedance flag [0/1]
* ``dust_pm10_exceed`` — PM10 exceedance flag [0/1]

**Transport Mode Considerations**

When ``erf.dust.transport_bins_separately = false`` (default):
  A single 3D scalar holds total dust mass. PM attribution uses bin 0 diameter.

When ``erf.dust.transport_bins_separately = true``:
  Each bin has its own 3D scalar. PM attribution is per-bin, giving separate mass
  contributions for PM2.5 and PM10 from each size bin.

**Implementation Limitation**

Bins straddling the 2.5 µm or 10 µm threshold are classified by their nominal diameter.
No size distribution or partial bin attribution is applied. For higher fidelity, future
phases may use a size-resolved bin moment or lognormal distribution within each bin.

Parameters
  * ``erf.dust.dust_naaqs_file`` [string] (default: ``"dust_naaqs.csv"``) — Path to NAAQS
   compliance CSV output (relative to run directory)

Implementation
  Phase 17 introduces two new files:

  * ``Source/Dust/ERF_DustPM.H`` — Header-only; provides size fractionation functions
  * ``Source/Dust/ERF_DustNAAQSOutput.H`` — Header-only; provides CSV writer

  ``DustLayer::compute_naaqs_diagnostics`` is called from ``DustLayer::advance()``
  each timestep, after ``extract_atm_return_fields`` has updated ``dust_conc_sfc``.

MSHA Worker Exposure Tracking Module
-------------------------------------

Phase 18 introduces a time-weighted average (TWA) dose accumulator for MSHA
(Mine Safety and Health Administration) worker exposure protection. The module
computes an 8-hour rolling TWA of respirable dust concentration (PM10) and
detects exceedances of the MSHA Permissible Exposure Limit (PEL) set at 5 mg/m³
under 30 CFR Part 56 (Surface Mines) [30CFR56]_.

**Physical Model**

The module processes the instantaneous PM10 concentration (computed by Phase 17)
to track cumulative worker exposure over an 8-hour shift:

.. math::

    C_\mathrm{resp}(t) = \mathrm{PM10}(t) \times 10^{-3} \, [\mathrm{mg/m^3}]

    \mathrm{dose}(t) = \int_0^t C_\mathrm{resp}(\tau) \, \mathrm{d}\tau \, [\mathrm{mg/m^3/h}]

    \mathrm{TWA}_{8h}(t) = \frac{\mathrm{dose}(t)}{8.0} \, [\mathrm{mg/m^3}]

    \mathrm{exceed}(t) = \begin{cases} 1 & \text{if } \mathrm{TWA}_{8h}(t) > 5 \, \mathrm{mg/m^3} \\
                                      0 & \text{otherwise} \end{cases}

Shift reset: When a new shift begins (detected by monitoring
:math:`\lfloor t / T_\mathrm{shift} \rfloor`), the dose accumulator is zeroed
and the peak TWA during the ending shift is recorded.

**Output Files**

Phase 18 produces three per-site CSV files in the run directory:

1. ``msha_exposure.csv`` — Per-timestep exposure snapshot
   * Columns: ``step, time_s, TWA_max_mg_m3, n_cells_exceed, dose_max_mg_m3_h``
   * TWA_max is the maximum TWA over all dust grid cells
   * n_cells_exceed counts how many cells exceed the PEL
   * dose_max is the maximum cumulative dose over all cells

2. ``msha_shift_summary.csv`` — End-of-shift summary
   * Columns: ``shift_number, shift_end_time_s, TWA_peak_mg_m3, n_cells_exceed``
   * Written once per shift when the shift boundary is crossed
   * TWA_peak is the maximum TWA at the end of that shift

3. ``msha_receptor_<name>.csv`` — Named receptor point samples (one file per receptor)
   * Columns: ``step, time_s, PM10_ug_m3``
   * Samples PM10 at fixed points (x, y) specified in input
   * Uses nearest-cell-only (no interpolation)
   * Allows fence-line and downwind monitoring stations

**Parameters**

Configurable MSHA parameters in input files via ``erf.dust.*`` prefix:

* ``msha_pel_mg_m3`` [Real] (default: ``5.0``) — PEL threshold [mg/m³]
* ``msha_shift_duration_s`` [Real] (default: ``28800.0`` = 8 hours) — Shift length [s]
* ``msha_exposure_file`` [string] (default: ``"msha_exposure.csv"``) — Per-step CSV
* ``msha_shift_file`` [string] (default: ``"msha_shift_summary.csv"``) — Shift summary CSV
* ``msha_receptor_names`` [array of strings] (default: empty) — Receptor point names
* ``msha_receptor_x`` [array of Real] (default: empty) — Receptor x-coordinates [m]
* ``msha_receptor_y`` [array of Real] (default: empty) — Receptor y-coordinates [m]

Example input:

.. code-block:: bash

   erf.dust.msha_pel_mg_m3 = 5.0              # 30 CFR 56.5001
   erf.dust.msha_shift_duration_s = 28800.0   # 8 hours
   erf.dust.msha_exposure_file = "msha_exposure.csv"
   erf.dust.msha_shift_file = "msha_shift_summary.csv"
   erf.dust.msha_receptor_names = "fence_N" "fence_S" "downwind"
   erf.dust.msha_receptor_x = 500.0 500.0 1500.0
   erf.dust.msha_receptor_y = 2900.0 100.0 1500.0

**Plotfile Integration**

Phase 18 adds 4 new variables to the 2D dust plotfile (total 18 fields):

* ``dust_msha_dose_mg_m3_h`` — Cumulative dose [mg/m³/h]
* ``dust_msha_twa_mg_m3`` — 8-hour TWA [mg/m³]
* ``dust_msha_exceed`` — PEL exceedance [0/1]
* ``dust_msha_shift_twa`` — TWA at last shift end [mg/m³]

**Implementation Limitation**

TWA is computed over the domain-averaged dust grid cells. This is the
instantaneous concentration at fixed monitoring points, not the personal
(time-weighted) exposure of a mobile worker. For full OSHA compliance,
measurements should use personal dust samplers on workers. This module
provides spatial diagnostics (which grid cells are above PEL) and
static-point time series (how concentration evolves at fenceline stations).

Implementation
  Phase 18 introduces two new files:

  * ``Source/Dust/ERF_DustMSHA.H`` — Header-only; provides dose and exceedance functions
  * ``Source/Dust/ERF_DustMSHAOutput.H`` — Header-only; provides CSV writers for exposure, shift summary, and receptors

  ``DustLayer::compute_msha_exposure`` is called from ``DustLayer::advance()``
  each timestep, after ``compute_naaqs_diagnostics``.

+Lagrangian Super-Particle Source-Receptor Attribution (Phase 19)
+===============================================================
+
+Phase 19 implements Lagrangian super-particles for source-receptor attribution
+tracking. The module tracks individual dust particles from their emission point
+to deposition, accumulating deposition mass as a function of source location.
+
+**Overview**
+
+A super-particle represents an ensemble of dust grains emitted from a single
+grid cell. Each particle carries:
+
+* Mass in [kg]
+* Stokes settling velocity [m/s]
+* Release time [s]
+* Source grid cell indices (i, j) as Real for sub-cell accuracy
+
+Particles are advected by nearest-cell interpolation of ERF face-staggered
+velocities (xvel, yvel, zvel) and subject to gravitational settling. When
+a particle descends below z = z_lo + 0.5*dz_atm, it is deposited and its
+mass is accumulated in ``dust_source_map``, a 2D field on the dust grid.
+
+**Mathematics**
+
+Position update (Euler step):
+
+.. math::
+
+   \mathbf{x}_{n+1} = \mathbf{x}_n + (\mathbf{u}, \mathbf{v}, w_{adj}) \Delta t
+
+where the adjusted vertical velocity is:
+
+.. math::
+
+   w_{adj} = w - v_s
+
+The settling velocity is computed from Stokes law with Cunningham slip correction
+(see Phase 11, ``compute_stokes_settling`` in ``ERF_DustSettling.H``):
+
+.. math::
+
+   v_s = \frac{(\rho_p - \rho_a) g d^2 C_c}{18 \mu_a}
+
+**Velocity Interpolation**
+
+Velocities are obtained at the particle position via nearest-cell averaging.
+For each component:
+
+.. math::
+
+   u = 0.5 \left( u_{i,j,k} + u_{i+1,j,k} \right)
+   v = 0.5 \left( v_{i,j,k} + v_{i,j+1,k} \right)
+   w = 0.5 \left( w_{i,j,k} + w_{i,j,k+1} \right)
+
+This pattern is identical to ``AdvectWithFlow`` in ``ERFPCEvolve.cpp``.
+
+**Deposition and Source Map**
+
+When z < z_lo + 0.5*dz_atm, the particle is removed and its mass is added
+to ``dust_source_map`` at the original source cell (src_i, src_j):
+
+.. math::
+
+   S(i,j) += m
+
+The source map is never reset during a simulation, so it accumulates total
+deposition mass per source cell over all timesteps. The field has units [kg/m²].
+
+**Compile Guard**
+
+All particle code is guarded by:
+
+.. code-block:: cpp
+
+   #if defined(ERF_USE_DUST) && defined(ERF_USE_PARTICLES)
+
+When ``ERF_USE_PARTICLES`` is not defined:
+
+* No particle container is created
+* ``dust_source_map`` is zeroed in plotfile output
+* Plotfile still contains 19 components (last one is zero)
+
+**Parameters**
+
+Phase 19 adds two parameters to the dust configuration:
+
+* ``erf.dust.enable_particles`` — Boolean; default false. When true, particles
+  are released and tracked.
+* ``erf.dust.particle_release_interval`` — Integer; default 1. Particles are
+  released every N steps. Set > 1 to reduce computational cost.
+
+**Plotfile Integration**
+
+Phase 19 adds 1 new variable to the 2D dust plotfile (total 19 fields):
+
+* ``dust_source_map`` — Cumulative source deposition [kg/m²]
+
+The ``dust_source_map`` field is only updated when particles are tracked
+(``enable_particles`` = true and ``ERF_USE_PARTICLES`` defined). Otherwise,
+this field is zero in plotfile output.
+
+**Implementation**
+
+Phase 19 introduces three new files:
+
+* ``Source/Particles/ERF_DustPC.H`` — Header; declares the ``ERFDustPC`` container
+* ``Source/Particles/ERF_DustPC.cpp`` — Implementation; ``ReleaseParticles`` and ``AdvanceParticles`` methods
+
+``DustLayer::advance_particles`` is called from ``DustLayer::advance()`` after
+``compute_msha_exposure``, but only when ``enable_particles`` = true.
+
+**References**
+
+* AMReX Particles: https://amrex-codes.github.io/amrex/docs_html/Particles.html
+* Phase 11 (Stokes settling): ``ERF_DustSettling.H`` ``compute_stokes_settling``
+* Phase 9 pattern (nearest-cell advection): ``ERFPCEvolve.cpp`` ``AdvectWithFlow``
+* Shao, Y. (2008). *Physics and Modelling of Wind Erosion*. Springer.
+* Seinfeld, J. H., and S. N. Pandis (2006). *Atmospheric Chemistry and Physics*,
+  2nd ed., Ch. 9. Wiley. ISBN 978-0-471-72018-8.
+
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

.. [EPA2024] U.S. EPA (2024). Review of the NAAQS for Particulate Matter.
   https://www.epa.gov/pm-pollution/review-national-ambient-air-quality-standards-naaqs-particulate-matter-pm

.. [Watson1994] Watson, J. G., J. C. Chow, L. C. Pritchett, W. R. Pierson,
   W. A. Frazier, and R. G. Egami (1994). Receptor modeling application
   framework for particle source apportionment.
   *Atmos. Environ.*, 28, 2493-2509.
   https://doi.org/10.1016/1352-2310(94)90400-6

.. [Parkhurst2013] Parkhurst, D. L., and C. A. J. Appelo (2013).
   Description of input and examples for PHREEQC version 3.
   *U.S. Geological Survey Techniques and Methods*, book 6, chap. A43.
   https://pubs.usgs.gov/tm/06/a43/

.. [EPA_AP42] U.S. EPA (1998). AP-42, Compilation of Air Emission Factors,
   Chapter 13.2.2: Unpaved Roads.
   https://www.epa.gov/air-emissions-factors-and-quantification/ap-42-compilation-air-emission-factors

.. [30CFR56] Title 30, Part 56. Safety and Health Regulations for Surface
   Mining. U.S. Mine Safety and Health Administration (MSHA), Department of Labor.
   https://www.ecfr.gov/current/title-30/chapter-I/subchapter-K/part-56

.. [30CFR57] Title 30, Part 57. Safety and Health Regulations for
   Underground Mining. U.S. Mine Safety and Health Administration (MSHA),
   Department of Labor.
   https://www.ecfr.gov/current/title-30/chapter-I/subchapter-K/part-57

.. [Stayner1998] Stayner, L. T., N. J. Dankovic, and R. A. Lemen (1998).
   Occupational exposure to chrysotile asbestos and cancer risk: A review of
   the amphibole hypothesis. *Am. J. Respir. Crit. Care Med.*, 157, 69-89.
   https://doi.org/10.1164/ajrccm.157.1.9703022

.. [Attfield1992] Attfield, M. D., and K. J. Morring (1992). An investigation
   into the relationship between coal rank and coal worker's pneumoconiosis.
   *Am. J. Ind. Med.*, 22, 417-428.
   https://doi.org/10.1002/ajim.4700220311

.. [Msha2016] MSHA (2016). Lowering Miners' Exposure to Respirable Coal
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
