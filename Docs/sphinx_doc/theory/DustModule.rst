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
