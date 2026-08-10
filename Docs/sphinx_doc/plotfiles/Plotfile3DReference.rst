.. _sec:Plotfile3DReference:

3D plotfile variable reference
==============================

Use ``erf.plot_vars_1`` and ``erf.plot_vars_2`` to select cell-centered and
derived fields for the two standard 3D plotfile streams. ERF writes each stream
on all active AMR levels. Some fields require a compile-time option, a physics
package, or stored diagnostic state; the tables below state those restrictions.

Time-averaged velocity fields require ``erf.time_avg_vel = true``. Set
``erf.plot_face_vels = true`` to write the native staggered velocity components
in separate face-centered outputs associated with the configured 3D stream.

3D output variables
-------------------

The fixed core names are selected from the input list and emitted in the order
defined by ERF. Some names require a compile-time option, a physics package, or
stored diagnostic state; those restrictions are noted below.

The default subvolume inventory is documented on :ref:`sec:Plotfiles`.

+-----------------------------+------------------+
| Variable                    | Definition       |
|                             |                  |
+=============================+==================+
| **x_velocity**              | Velocity in x    |
|                             | direction        |
|                             | [m/s]            |
+-----------------------------+------------------+
| **y_velocity**              | Velocity in y    |
|                             | direction        |
|                             | [m/s]            |
+-----------------------------+------------------+
| **z_velocity**              | Velocity in z    |
|                             | direction        |
|                             | [m/s]            |
+-----------------------------+------------------+
| **density**                 | Dry density      |
|                             | [kg/m^3]         |
|                             |                  |
+-----------------------------+------------------+
| **moist_density**           | Dry-air density  |
|                             | plus vapor and   |
|                             | non-precipitating|
|                             | condensate       |
|                             | density [kg/m^3] |
|                             |                  |
+-----------------------------+------------------+
| **dens_hse**                | Hydrostatic      |
|                             | density          |
|                             | [kg/m^3]         |
+-----------------------------+------------------+
| **pert_dens**               | Perturbational   |
|                             | density          |
|                             | [kg/m^3]         |
+-----------------------------+------------------+
| **pressure**                | Total pressure   |
|                             | [Pa]             |
|                             |                  |
+-----------------------------+------------------+
| **pres_hse**                | Hydrostatic      |
|                             | pressure         |
|                             | [Pa]             |
+-----------------------------+------------------+
| **theta_hse**               | Hydrostatic      |
|                             | potential        |
|                             | temperature [K]  |
+-----------------------------+------------------+
| **qv_hse**                  | Base-state water |
|                             | vapor mixing     |
|                             | ratio [kg/kg]    |
+-----------------------------+------------------+
| **pert_pres**               | Perturbational   |
|                             | pressure         |
|                             | [Pa]             |
+-----------------------------+------------------+
| **pres_hse_x**              | Derivative of    |
|                             | hydrostatic      |
|                             | pressure in x    |
|                             | [Pa/m]           |
+-----------------------------+------------------+
| **pres_hse_y**              | Derivative of    |
|                             | hydrostatic      |
|                             | pressure in y    |
|                             | [Pa/m]           |
+-----------------------------+------------------+
| **dpdx**                    | Pressure gradient|
|                             | in x direction   |
|                             | [Pa/m]           |
+-----------------------------+------------------+
| **dpdy**                    | Pressure gradient|
|                             | in y direction   |
|                             | [Pa/m]           |
+-----------------------------+------------------+
| **dpdz**                    | Pressure gradient|
|                             | in z direction   |
|                             | [Pa/m]           |
+-----------------------------+------------------+
| **temp**                    | Temperature      |
|                             | [K]              |
|                             |                  |
+-----------------------------+------------------+
| **theta**                   | Potential        |
|                             | temperature [K]  |
|                             |                  |
+-----------------------------+------------------+
| **buoyancy**                | Buoyancy term    |
|                             | used by the      |
|                             | z-momentum       |
|                             | equation         |
+-----------------------------+------------------+
| **eq_pot_temp**             | Equivalent       |
|                             | potential        |
|                             | temperature [K]  |
+-----------------------------+------------------+
| **VPD**                     | Vapor pressure   |
|                             | deficit [kPa]    |
|                             |                  |
+-----------------------------+------------------+
| **rhotheta**                | Density * theta  |
|                             | [kg K/m^3]       |
|                             |                  |
+-----------------------------+------------------+
| **KE**                      | Prognostic       |
|                             | turbulent kinetic|
|                             | energy when      |
|                             | populated by the |
|                             | active closure,  |
|                             | including native |
|                             | SHOC [m^2/s^2]   |
+-----------------------------+------------------+
| **rhoKE**                   | Density * KE     |
|                             | [kg/(m s^2)]     |
|                             |                  |
+-----------------------------+------------------+
| **scalar**                  | Scalar magnitude |
|                             | [problem-dep.]   |
|                             |                  |
+-----------------------------+------------------+
| **reflectivity**            | reflectivity     |
|                             | cell-by-cell     |
|                             | [dBZ]            |
+-----------------------------+------------------+
| **max_reflectivity**        | max of           |
|                             | reflectivity     |
|                             | over a column    |
|                             | [dBZ]            |
+-----------------------------+------------------+
| **precipitable**            | precipitable     |
|                             | water (integral  |
|                             | over column)     |
|                             | [kg/m^2]         |
+-----------------------------+------------------+
| **mucape**                  | most unstable    |
|                             | CAPE over a      |
|                             | column [J/kg]    |
+-----------------------------+------------------+
| **vorticity_x**             | x-component of   |
|                             | vorticity [1/s]  |
|                             |                  |
+-----------------------------+------------------+
| **vorticity_y**             | y-component of   |
|                             | vorticity [1/s]  |
|                             |                  |
+-----------------------------+------------------+
| **vorticity_z**             | z-component of   |
|                             | vorticity [1/s]  |
|                             |                  |
+-----------------------------+------------------+
| **local_helicity**          | helicity         |
|                             | cell-by-cell     |
|                             | [m/s^2]          |
+-----------------------------+------------------+
| **helicity**                | helicity         |
|                             | (integral over   |
|                             | column)          |
|                             | [m^2/s^2]        |
+-----------------------------+------------------+
| **magvel**                  | magnitude of     |
|                             | velocity [m/s]   |
|                             |                  |
+-----------------------------+------------------+
| **divU**                    | divergence of    |
|                             | velocity [1/s]   |
|                             |                  |
+-----------------------------+------------------+
| **u_t_avg**                 | time average of  |
|                             | x-component of   |
|                             | velocity [m/s]   |
+-----------------------------+------------------+
| **v_t_avg**                 | time average of  |
|                             | y-component of   |
|                             | velocity [m/s]   |
+-----------------------------+------------------+
| **w_t_avg**                 | time average of  |
|                             | z-component of   |
|                             | velocity [m/s]   |
+-----------------------------+------------------+
| **umag_t_avg**              | time average of  |
|                             | velocity mag     |
|                             | [m/s]            |
+-----------------------------+------------------+
| **rhoadv_0**                | Conserved scalar |
|                             | [problem-dep.]   |
|                             |                  |
+-----------------------------+------------------+
| **soundspeed**              | Sound speed      |
|                             | [m/s]            |
|                             |                  |
+-----------------------------+------------------+
| **z_phys**                  | Terrain height   |
|                             | [m]              |
|                             |                  |
+-----------------------------+------------------+
| **detJ**                    | Jacobian         |
|                             | determinant [1]  |
|                             |                  |
+-----------------------------+------------------+
| **mapfac**                  | Map scale factor |
|                             | [1]              |
|                             |                  |
+-----------------------------+------------------+
| **lat_m**                   | Latitude at mass |
|                             | points           |
|                             | [deg]            |
+-----------------------------+------------------+
| **lon_m**                   | Longitude at     |
|                             | mass points      |
|                             | [deg]            |
+-----------------------------+------------------+
| **nut**                     | Eddy viscosity,  |
|                             | nu_t [m^2/s]     |
+-----------------------------+------------------+
| **Kmv**                     | Vertical         |
|                             | Eddy Diffusivity |
|                             | of Momentum      |
|                             | [kg/(m s)]       |
+-----------------------------+------------------+
| **Kmh**                     | Horizontal       |
|                             | Eddy Diffusivity |
|                             | of Momentum      |
|                             | (Note: For LES,  |
|                             | this is the      |
|                             | _dynamic_ eddy   |
|                             | viscosity, mu_t  |
|                             | = rho * nu_t     |
|                             | and Kmh==Kmv)    |
|                             | [kg/(m s)]       |
+-----------------------------+------------------+
| **Khv**                     | Vertical         |
|                             | Eddy Diffusivity |
|                             | of Heat          |
|                             | [kg/(m s)]       |
+-----------------------------+------------------+
| **Khh**                     | Horizontal       |
|                             | Eddy Diffusivity |
|                             | of Heat          |
|                             | [kg/(m s)]       |
+-----------------------------+------------------+
| **Lturb**                   | Turbulence       |
|                             | length scale     |
|                             | with             |
|                             | ``use_kturb``    |
|                             | [m]              |
+-----------------------------+------------------+
| **walldist**                | Wall distance    |
|                             | for RANS models  |
|                             | only [m]         |
+-----------------------------+------------------+
| **diss**                    | Subfilter-scale  |
|                             | dissipation      |
|                             | with diffusion / |
|                             | turbulence       |
|                             | [kg/(m s^3)]     |
+-----------------------------+------------------+
| **qt**                      | Total water      |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **qn**                      | Nonprecipitating |
|                             | water (qv + qc + |
|                             | qi)              |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **qp**                      | Precipitating    |
|                             | water (rain +    |
|                             | snow + graupel)  |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **qc**                      | Cloud water      |
|                             | mixing ratio     |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **qi**                      | Cloud ice        |
|                             | mixing ratio     |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **qrain**                   | Rain-water       |
|                             | mixing ratio     |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **qsnow**                   | Snow mixing      |
|                             | ratio [kg/kg]    |
+-----------------------------+------------------+
| **qgraup**                  | Graupel mixing   |
|                             | ratio [kg/kg]    |
+-----------------------------+------------------+
| **qv**                      | Water vapor      |
|                             | mixing ratio     |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **qsat**                    | Saturation water |
|                             | vapor mixing     |
|                             | ratio [kg/kg]    |
+-----------------------------+------------------+
| **rain_accum**              | Accumulated rain |
|                             | amount with      |
|                             | precipitating    |
|                             | moisture models  |
|                             | [mm]             |
+-----------------------------+------------------+
| **snow_accum**              | Accumulated snow |
|                             | amount with SAM  |
|                             | or Morrison      |
|                             | microphysics     |
|                             | [mm]             |
+-----------------------------+------------------+
| **graup_accum**             | Accumulated      |
|                             | graupel amount   |
|                             | with SAM or      |
|                             | Morrison         |
|                             | microphysics     |
|                             | [mm]             |
+-----------------------------+------------------+
| **rel_humidity**            | Relative         |
|                             | humidity;        |
|                             | currently filled |
|                             | only for         |
|                             | SuperDroplets    |
|                             | [1]              |
+-----------------------------+------------------+
| **condensation_rate**       | Condensation     |
|                             | rate with        |
|                             | SuperDroplets    |
|                             | only             |
|                             | [kg/kg/s]        |
+-----------------------------+------------------+
| **terrain_IB_mask**         | Immersed-boundary|
|                             | terrain/building |
|                             | mask; available  |
|                             | for immersed     |
|                             | forcing terrain  |
|                             | or buildings     |
|                             | [1]              |
+-----------------------------+------------------+
| **volfrac**                 | EB / immersed    |
|                             | boundary volume  |
|                             | fraction; unity  |
|                             | elsewhere        |
|                             | [1]              |
+-----------------------------+------------------+
| **qsrc_sw**                 | Shortwave        |
|                             | radiative        |
|                             | heating source   |
|                             | term with        |
|                             | radiation        |
|                             | [K/s]            |
+-----------------------------+------------------+
| **qsrc_lw**                 | Longwave         |
|                             | radiative        |
|                             | heating source   |
|                             | term with        |
|                             | radiation        |
|                             | [K/s]            |
+-----------------------------+------------------+
| **tracer_particles_count**  | Tracer particle  |
|                             | count per cell   |
|                             | requires         |
|                             | ERF_USE_PARTICLES|
|                             | to be defined    |
|                             | [count]          |
+-----------------------------+------------------+

The ``qrain``, ``qsnow``, and ``qgraup`` rows are available when the active
moisture scheme provides the corresponding rain, snow, or graupel component.

Fixed-field capability matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fixed moisture names are selected only when both the active scheme gives the
source its documented physical meaning and the corresponding conserved or
auxiliary component exists. A reduced scheme can retain a fixed-width state
vector without making every slot a valid diagnostic. Unsupported requests are
warned about and omitted before plotfile allocation.

The following matrix summarizes the fixed moisture capabilities. ``qv`` is
vapor, ``qc`` is cloud liquid, ``qi`` is cloud ice, and ``qr``, ``qs``, and
``qg`` are rain, snow, and graupel. ``N`` denotes the corresponding number
concentration and ``A`` the corresponding accumulation field.

.. list-table::
   :header-rows: 1
   :widths: 27 31 25 22

   * - Moisture type
     - Mass fields
     - Number fields
     - Accumulations
   * - ``None``
     - none
     - none
     - none
   * - ``MoistNoCondensation``
     - qv
     - none
     - none
   * - ``SatAdj``
     - qv, qc
     - none
     - none
   * - ``Kessler_NoRain``
     - qv, qc
     - none
     - none
   * - ``Kessler``
     - qv, qc, qr
     - none
     - rain
   * - ``SAM_NoPrecip_NoIce``
     - qv, qc
     - none
     - none
   * - ``SAM_NoIce``
     - qv, qc, qr
     - none
     - rain
   * - ``SAM``
     - qv, qc, qi, qr, qs, qg
     - none
     - rain, snow, graupel
   * - ``Morrison_NoIce``
     - qv, qc, qr
     - Nc, Nr
     - rain
   * - ``Morrison``
     - qv, qc, qi, qr, qs, qg
     - Nc, Ni, Nr, Ns, Ng
     - rain, snow, graupel
   * - ``WSM6``
     - qv, qc, qi, qr, qs, qg
     - none
     - rain, snow, graupel
   * - ``SuperDroplets``
     - qv, qc, qr
     - none
     - rain

For ``SuperDroplets``, water vapor, cloud water, and rain water occupy the
fixed ``RhoQ1``, ``RhoQ2``, and ``RhoQ3`` conserved components. The fixed
``qv``, ``qc``, ``qrain``, ``qt``, ``qn``, and ``qp`` diagnostics are therefore
available. For this scheme, ``qt = qv + qc + qrain``, ``qn = qv + qc``, and
``qp = qrain``. Provider-generated names such as ``qv_<species>``,
``qc_<species>``, and ``qt_<species>`` apply only to additional non-water
condensable species. The fixed ``rain_accum``, ``rel_humidity``, and
``condensation_rate`` names are available only when their documented
auxiliary qmoist storage is present.
The aggregate ``qt``, ``qn``, ``qp``, ``moist_density``, ``qsat``, and
``precipitable`` fields additionally require their source-state bounds.

In general, an aggregate fixed field is selected only when the complete
inclusive q-component range used by its writer is present in the actual
conserved state. SuperDroplets satisfies this rule for its q1:q3 water state.

For the fixed mixed-phase layouts, ``qi`` is read from ``RhoQ3`` and
``qrain`` from ``RhoQ4``. Warm-rain layouts have no ice component and place
``qrain`` in ``RhoQ3``. Selection checks the exact source component used by
the writer.

Optional storage restrictions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following fixed diagnostics are selected only if the storage exists on
every AMR level in the plotfile:

* ``u_t_avg``, ``v_t_avg``, ``w_t_avg``, and ``umag_t_avg`` require
  ``erf.time_avg_vel = true``. If no samples have been accumulated yet, the
  output value is defined as zero rather than dividing by zero.
* ``qsrc_sw`` and ``qsrc_lw`` require a non-``None`` radiation choice.
* ``nut``, ``Kmv``, ``Kmh``, ``Khv``, ``Khh``, and ``Lturb`` require
  ``use_kturb = true`` at every AMR level.
* ``diss`` requires molecular diffusion or ``use_kturb`` at every level.
* ``walldist`` requires a non-``None`` RANS choice at every level.

These restrictions apply to fixed names only. Dynamic names supplied by an
active microphysics or particle provider remain governed by that provider's
own plot-variable contract.

Fixed conserved-state fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The conserved-state inventory also includes the moisture-density components
below. The active moisture model determines which components are selectable.

.. code-block:: text

   rhoQ1 rhoQ2 rhoQ3 rhoQ4 rhoQ5 rhoQ6
   rhoQ7 rhoQ8 rhoQ9 rhoQ10 rhoQ11

These are conserved moisture or species densities. The component names are
fixed; the active moisture model controls which components are available.
ERF checks the active conserved-state and microphysics sizes before selecting a
requested component. Unsupported fixed names are omitted and reported as
unavailable; they do not reserve a plotfile component. Provider-supplied names
remain dynamic and are appended only when the active provider exposes them.

Velocity output behavior
~~~~~~~~~~~~~~~~~~~~~~~~

Requesting any one of ``x_velocity``, ``y_velocity``, or ``z_velocity`` selects
all three cell-centered velocity components. When ``erf.plot_face_vels = true``,
native AMReX output also writes three separate face-centered plotfiles with the
components ``x_velocity_stag``, ``y_velocity_stag``, and ``z_velocity_stag``.

For stream 1, the face files use ``<plot_file_1>U<step>``,
``<plot_file_1>V<step>``, and ``<plot_file_1>W<step>``. Stream 2 uses the
corresponding ``<plot_file_2>U<step>``, ``V<step>``, and ``W<step>`` prefixes.
These face-velocity filenames always use level-0 step numbering, including when
``erf.use_real_time_in_pltname = true``. The face outputs are native AMReX
artifacts; the NetCDF writer does not enter this face-output path.

Supplemental fixed derived fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The main table contains the ordinary fixed fields. The following fixed names
are conditional or specialized and remain part of the same source-defined
inventory. A requested name is still subject to the runtime selection
conditions in ``setPlotVariables``.

The following fixed fields are supplied by the native SHOC driver. The
native-only fields are written as ``-999`` when native SHOC diagnostics are not
available.

.. list-table::
   :header-rows: 1
   :widths: 24 20 56

   * - Variable
     - Units
     - Native SHOC meaning
   * - ``pblh``
     - m
     - Native SHOC ``pblh`` is reported in metres above local ground (AGL).
   * - ``shoc_cldfrac``
     - 1
     - Subgrid cloud fraction diagnosed by the native SHOC PDF.
   * - ``shoc_ql``
     - kg/kg
     - Cloud-liquid mixing ratio diagnosed by the native SHOC PDF.
   * - ``shoc_ql2``
     - (kg/kg)^2
     - Variance of the diagnosed cloud-liquid mixing ratio.
   * - ``shoc_cond``
     - kg/kg/s
     - Positive PDF condensate-change rate. It is zero unless ``erf.shoc.extra_shoc_diags = true``.
   * - ``wqls_sec``
     - (kg/kg) m/s
     - Turbulent vertical flux of the diagnosed SHOC cloud-liquid quantity.
   * - ``wthv_sec``
     - K m/s
     - Turbulent vertical flux of virtual potential temperature.
   * - ``w_sec``
     - m^2/s^2
     - Vertical-velocity variance.
   * - ``thl_sec``
     - K^2
     - Liquid-water potential-temperature variance.
   * - ``qw_sec``
     - (kg/kg)^2
     - Total-water variance.
   * - ``qwthl_sec``
     - K kg/kg
     - Covariance of total water and liquid-water potential temperature.
   * - ``wthl_sec``
     - K m/s
     - Turbulent vertical flux of liquid-water potential temperature.
   * - ``wqw_sec``
     - (kg/kg) m/s
     - Turbulent vertical flux of total water.
   * - ``w3``
     - m^3/s^3
     - Third central moment of vertical velocity.
   * - ``brunt``
     - s^-2
     - Squared Brunt-Vaisala frequency, :math:`N^2`.
   * - ``isotropy``
     - s
     - Native SHOC isotropy timescale.
   * - ``shear_prod``
     - m^2/s^3
     - Shear-production contribution to the native SHOC TKE budget.
   * - ``buoy_prod``
     - m^2/s^3
     - Buoyancy production or destruction contribution to the native SHOC TKE budget.
   * - ``diss_tke``
     - m^2/s^3
     - Dissipation contribution diagnosed for the native SHOC TKE budget.

When native SHOC diagnostics are available, the standard turbulence fields use
the native SHOC coefficients: ``Kmv`` is :math:`\rho K_m`, ``Khv`` is
:math:`\rho K_h`, ``Lturb`` is the native SHOC mixing length, and ``nut`` is
the kinematic momentum diffusivity ``Kmv / density``.

* ``nc``, ``ni``, ``nr``, ``ns``, and ``ng`` are moisture number
  concentrations. They are available when the active moisture model provides
  the corresponding conserved component; the output units follow that model's
  number-concentration convention.
* ``xvel_err``, ``yvel_err``, ``zvel_err``, and ``pp_err`` are error diagnostics
  available when ``ERF_COMPUTE_ERROR`` is enabled. Their units follow the
  corresponding velocity or pressure quantity.

The fixed conserved-state inventory is ``density``, ``rhotheta``, ``rhoKE``,
``rhoadv_0``, and ``rhoQ1`` through ``rhoQ11``. The ``rhoQ`` components are
active only to the extent that the selected moisture model provides them.

Wind-farm-only variables
~~~~~~~~~~~~~~~~~~~~~~~~

The following quantities are available only in builds with
``ERF_USE_WINDFARM`` enabled.

+-----------------------------+------------------+
| Variable                    | Definition       |
+=============================+==================+
| **num_turb**                | Number of wind   |
|                             | turbines in cell |
|                             | for Fitch, EWP,  |
|                             | SimpleAD, and    |
|                             | GeneralAD        |
|                             | [count]          |
+-----------------------------+------------------+
| **SMark0**                  | Windfarm marker  |
|                             | component 0 for  |
|                             | Fitch, EWP,      |
|                             | SimpleAD, and    |
|                             | GeneralAD        |
|                             | [1]              |
+-----------------------------+------------------+
| **SMark1**                  | Windfarm marker  |
|                             | component 1 for  |
|                             | SimpleAD and     |
|                             | GeneralAD        |
|                             | [1]              |
+-----------------------------+------------------+

Morrison microphysics variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When using Morrison two-moment microphysics, additional diagnostic variables
are available for output. These variables provide detailed information about
cloud and precipitation processes. To enable Morrison output, include any of
the variables below in your **erf.plot_vars_1** or **erf.plot_vars_2** list.

**Thermodynamic State Variables:**

+-----------------------------+------------------+
| Variable                    | Definition       |
+=============================+==================+
| **micro_rho**               | Air density      |
|                             | [kg/m^3]         |
+-----------------------------+------------------+
| **micro_theta**             | Potential        |
|                             | temperature [K]  |
+-----------------------------+------------------+
| **micro_temp**              | Absolute         |
|                             | temperature [K]  |
+-----------------------------+------------------+
| **micro_pres**              | Pressure [Pa]    |
|                             |                  |
+-----------------------------+------------------+

**Non-Precipitating Moisture Variables (mixing ratios in kg/kg):**

+-----------------------------+------------------+
| Variable                    | Definition       |
+=============================+==================+
| **micro_qv**                | Water vapor      |
|                             | mixing ratio     |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **micro_qc**                | Cloud liquid     |
|                             | water mixing     |
|                             | ratio            |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **micro_qi**                | Cloud ice        |
|                             | mixing ratio     |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **micro_qn**                | Total cloud      |
|                             | condensate       |
|                             | (qc + qi)        |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **micro_qt**                | Total water      |
|                             | mixing ratio     |
|                             | (qv + qn)        |
|                             | [kg/kg]          |
+-----------------------------+------------------+

**Precipitating Hydrometeor Variables (mixing ratios in kg/kg):**

+-----------------------------+------------------+
| Variable                    | Definition       |
+=============================+==================+
| **micro_qp**                | Total            |
|                             | precipitation    |
|                             | (qrain + qsnow + |
|                             | qgraup)          |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **micro_qrain**             | Rain water       |
|                             | mixing ratio     |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **micro_qsnow**             | Snow mixing      |
|                             | ratio            |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **micro_qgraup**            | Graupel mixing   |
|                             | ratio            |
|                             | [kg/kg]          |
+-----------------------------+------------------+

**Number Concentrations (1/kg):**

+-----------------------------+------------------+
| Variable                    | Definition       |
+=============================+==================+
| **micro_nc**                | Cloud droplet    |
|                             | number           |
|                             | concentration    |
|                             | [1/kg]           |
+-----------------------------+------------------+
| **micro_nr**                | Rain drop number |
|                             | concentration    |
|                             | [1/kg]           |
+-----------------------------+------------------+
| **micro_ni**                | Cloud ice number |
|                             | concentration    |
|                             | [1/kg]           |
+-----------------------------+------------------+
| **micro_ns**                | Snow number      |
|                             | concentration    |
|                             | [1/kg]           |
+-----------------------------+------------------+
| **micro_ng**                | Graupel number   |
|                             | concentration    |
|                             | [1/kg]           |
+-----------------------------+------------------+

**Dynamical Variables:**

+-----------------------------+------------------+
| Variable                    | Definition       |
+=============================+==================+
| **micro_omega**             | Grid-scale       |
|                             | vertical         |
|                             | velocity [m/s]   |
|                             | used as input to |
|                             | Morrison scheme  |
+-----------------------------+------------------+

**Example Usage:**

To output Morrison diagnostic variables, add them to your plot variables list:

.. code-block:: text

   erf.plot_vars_1 = density theta qv micro_qc micro_qrain micro_nc micro_nr

This will output the base ERF variables (density, theta, qv) along with Morrison
cloud water, rain water, cloud droplet number concentration, and rain drop number
concentration.

Scheme-provided dynamic fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After the fixed names are selected, ERF asks the active microphysics provider for
additional plot names. There is no universal static list for these fields.

The current providers expose the following families:

* Morrison exposes the 19 ``micro_*`` names listed in the section above.
* SuperDroplets generates ``qv_<species>``, ``qc_<species>``,
  ``qt_<species>``, ``sat_ratio_<species>``, and ``accum_<species>`` names for
  configured species, plus ``accum_<aerosol>`` names for configured aerosols.
* Providers that inherit the empty ``NullMoist`` plot-name implementation add
  no scheme-specific names.

Particle-provided fields
~~~~~~~~~~~~~~~~~~~~~~~~

Particle builds can add Eulerian mesh fields after particle containers are
initialized. ERF prefixes each provider name with the particle-container name.
The default particle container provides ``<container>_mass_density``. The
SuperDroplet provider additionally generates mass-flux, number-density,
species-density, species-flux, aerosol-density, and aerosol-flux families from
the configured species and aerosols. Particle count names use the form
``<container>_count``. These fields are dynamic; their exact names depend on the
particle configuration. Every configured container can be selected, including
one that has not yet been allocated; its mesh count is then initialized to zero.

Before constructing a 3D plotfile in a Lagrangian-microphysics run with
two-way AMR coupling and at least two AMR levels, ERF averages fine-level
microphysics storage, ``RhoTheta``, and
the active moist conserved components onto covered coarse cells. This is a
consistency update so coarse output reflects fine-level deposits; it is not a
physical time tendency. Other 3D output assembly is diagnostic construction.
