.. role:: cpp(code)
  :language: c++

.. _sec:Plotfiles:

*********
Plotfiles
*********
.. toctree::
   :maxdepth: 1

There are three plotfile output paths in ERF.

The standard plotfile writes 3D data on all active AMR levels. The user selects
the variables with the plotfile variable lists.

The 2D plotfile writes a pseudo-2D slab. It stores fields that depend on
horizontal position, such as map factors, latitude, longitude, surface-layer
diagnostics, selected surface fluxes, surface pressure, and column-integrated
water vapor.

The subvolume plotfile writes 3D data from one selected region of the domain.
Primary native 2D and 3D AMReX plotfiles contain the same execution and
ancestry record in ``job_info``, with a distinct artifact UUID for each
output. NetCDF output does not yet contain this metadata. See
:ref:`sec:Provenance`.

Controlling PlotFile Generation
===============================

Plotfiles can be written very efficiently in parallel in a native AMReX format.
They can also be written in NetCDF. It is possible to output plotfiles in the
same or separate formats at two distinct frequencies.

The computational cost associated with writing plotfiles
in the AMReX native format is typically negligible relative to the overall cost of the simulation;
in a recent performance study the cost of writing a plotfile was roughly a percent or two
of the cost of a single timestep.

If NetCDF output is preferred, one suggestion is to write the plotfiles in the native AMReX
format for efficient I/O performance, then to convert the plotfiles to NetCDF files using
the executable you can build in Exec/Tools (using gmake, or with the ``ERF_ENABLE_TOOLS`` flag
if using cmake).

The following options in the inputs file control the generation of plotfiles.
Note that plotfiles can be written at two different frequencies; the names,
frequency and content of the two streams are controlled separately.

.. _list-of-parameters-9:

List of Parameters for Both 2D and 3D Plotfiles
-----------------------------------------------

+----------------------------------+------------------+-----------------------+------------+
| Parameter                        | Definition       | Acceptable            | Default    |
|                                  |                  | Values                |            |
+----------------------------------+------------------+-----------------------+------------+
| **erf.plotfile_type**            | AMReX or NETCDF  | "amrex" or            | "amrex"    |
|                                  | format           | "netcdf" or "NetCDF"  |            |
+----------------------------------+------------------+-----------------------+------------+
| **erf.use_real_time_in_pltname** | Use real time    | Boolean               | false      |
|                                  | instead of time  |                       |            |
|                                  | step for         |                       |            |
|                                  | plotfile names   |                       |            |
+----------------------------------+------------------+-----------------------+------------+
| **erf.file_name_digits**         | Number of digits | Integer               | 5          |
|                                  | to be appended   | :math:`> 0`           |            |
|                                  | to the plotfile  |                       |            |
|                                  | and checkpoint   |                       |            |
|                                  | file names if    |                       |            |
|                                  | using time step  |                       |            |
+----------------------------------+------------------+-----------------------+------------+

List of Parameters for 3D Plotfiles
-----------------------------------

.. list-table::
   :header-rows: 1
   :widths: 24 30 24 14

   * - Parameter
     - Definition
     - Acceptable Values
     - Default
   * - ``erf.plot_file_1``
     - Prefix for plotfiles at first output frequency.
     - String
     - ``plt_1_*``
   * - ``erf.plot_file_2``
     - Prefix for plotfiles at second output frequency.
     - String
     - ``plt_2_*``
   * - ``erf.plot_int_1``
     - Write plot files every this many level-0 time steps for the first stream.
     - Integer :math:`> 0`
     - -1
   * - ``erf.plot_int_2``
     - Write plot files every this many level-0 time steps for the second stream.
     - Integer :math:`> 0`
     - -1
   * - ``erf.plot_per_1``
     - Write plot files every this much simulation time for the first stream.
     - Real :math:`> 0`
     - -1.0
   * - ``erf.plot_per_2``
     - Write plot files every this much simulation time for the second stream.
     - Real :math:`> 0`
     - -1.0
   * - ``erf.plot_vars_1``
     - Variables to include in the first plotfile stream.
     - List of names
     - None
   * - ``erf.plot_vars_2``
     - Variables to include in the second plotfile stream.
     - List of names
     - None
   * - ``erf.plot_face_vels``
     - Output ``{prefix}U``, ``{prefix}V``, and ``{prefix}W`` with velocity
       components on the staggered grid.
     - Boolean
     - false

List of Parameters for 2D Plotfiles
-----------------------------------

+-----------------------------+--------------------------------------+-----------------------+--------------+
| Parameter                   | Definition                           | Acceptable Values     | Default      |
+=============================+======================================+=======================+==============+
| **erf.plot2d_file_1**       | Prefix for 2D plotfiles at the first | String                | ``plt2d_1_`` |
|                             | output frequency.                    |                       |              |
+-----------------------------+--------------------------------------+-----------------------+--------------+
| **erf.plot2d_file_2**       | Prefix for 2D plotfiles at the       | String                | ``plt2d_2_`` |
|                             | second output frequency.             |                       |              |
+-----------------------------+--------------------------------------+-----------------------+--------------+
| **erf.plot2d_int_1**        | Write 2D plotfiles every this many   | Integer :math:`> 0`   | -1           |
|                             | level-0 time steps for stream 1.     |                       |              |
+-----------------------------+--------------------------------------+-----------------------+--------------+
| **erf.plot2d_int_2**        | Write 2D plotfiles every this many   | Integer :math:`> 0`   | -1           |
|                             | level-0 time steps for stream 2.     |                       |              |
+-----------------------------+--------------------------------------+-----------------------+--------------+
| **erf.plot2d_per_1**        | Write 2D plotfiles every this much   | Real :math:`> 0`      | -1.0         |
|                             | simulation time for stream 1.        |                       |              |
+-----------------------------+--------------------------------------+-----------------------+--------------+
| **erf.plot2d_per_2**        | Write 2D plotfiles every this much   | Real :math:`> 0`      | -1.0         |
|                             | simulation time for stream 2.        |                       |              |
+-----------------------------+--------------------------------------+-----------------------+--------------+
| **erf.plot2d_vars_1**       | Variables to include in the first    | List of names         | None         |
|                             | 2D plotfile stream.                  |                       |              |
+-----------------------------+--------------------------------------+-----------------------+--------------+
| **erf.plot2d_vars_2**       | Variables to include in the second   | List of names         | None         |
|                             | 2D plotfile stream.                  |                       |              |
+-----------------------------+--------------------------------------+-----------------------+--------------+

Notes
-----

- ERF writes 2D plotfiles as one-cell-thick horizontal slabs.
- The two 2D streams are independent. Each stream has its own file prefix,
  write interval, write period, and variable list.
- Variables in a 2D plotfile appear in ERF's canonical order. The order in
  ``erf.plot2d_vars_1`` or ``erf.plot2d_vars_2`` does not change the component
  order in the file.
- If a requested 2D variable is not available, ERF skips it and prints a
  warning that names the input parameter and the skipped variable.
- NetCDF plotfile output requires an ERF build with NetCDF enabled.

List of Parameters for Subvolumes
-----------------------------------

+-----------------------------+-------------------+-----------------------+---------------+
| Parameter                   | Definition        | Acceptable            | Default       |
|                             |                   | Values                |               |
+=============================+===================+=======================+===============+
| **erf.subvol_file**         | prefix for        | String                | “*subvol*”    |
|                             | subvolume         |                       |               |
|                             | file names        |                       |               |
+-----------------------------+-------------------+-----------------------+---------------+
| **erf.subvol_int**          | how often (by     | Integer               | -1            |
|                             | level-0 time      | :math:`> 0`           |               |
|                             | steps) to write   |                       |               |
|                             | subvol files      |                       |               |
+-----------------------------+-------------------+-----------------------+---------------+
| **erf.subvol_per**          | how often (in     | Real                  | -1.0          |
|                             | simulation time)  | :math:`> 0`           |               |
|                             | to write          |                       |               |
|                             | subvol files      |                       |               |
+-----------------------------+-------------------+-----------------------+---------------+
| **erf.subvol.origin**       | lower left corner | Reals                 | None -- must  |
|                             | of region to be   |                       | be specified  |
|                             | output            |                       | if outputting |
|                             |                   |                       | subvolumes    |
+-----------------------------+-------------------+-----------------------+---------------+
| **erf.subvol.nxnynz**       | dimensions        | Integers              | None -- must  |
|                             | of region to be   |                       | be specified  |
|                             | output            |                       | if outputting |
|                             |                   |                       | subvolumes    |
+-----------------------------+-------------------+-----------------------+---------------+
| **erf.subvol.dxdydz**       | resolution        | Reals                 | None -- must  |
|                             | of region to be   |                       | be specified  |
|                             | output            |                       | if outputting |
|                             |                   |                       | subvolumes    |
+-----------------------------+-------------------+-----------------------+---------------+

.. _notes-5:

Notes
-----

-  The NetCDF option for writing plotfiles is only available if ERF has been built with USE_NETCDF enabled.

.. _examples-of-usage-8:

Examples of Usage
-----------------

-  **erf.plotfile_type** = *amrex*

-  **erf.plot_file_1** = *plt_run*

-  **erf.plot_int_1** = 10

   means that native plot files (actually directories) starting with the prefix
   “*plt_run*” will be generated every 10 level-0 time steps. If using
   amrex format, that directory names will be *plt_run00000*, *plt_run00010*,
   *plt_run00020*, etc.  If using NetCDF format, the names will have ".nc" appended.

   In addition, while the amrex plotfiles will contain data at all of the refinement
   levels,  NetCDF files are separated by level.

PlotFile Outputs
================

Plotfiles can include the quantities of several simulation parameters as output.
They are summarized in the list below. Note that temporally averaged quantities
(e.g., ``u_t_avg, v_t_avg, w_t_avg, umag_t_avg``) require the user to enable the
storage of the time averaged variables with ``erf.time_avg_vel = true``.
Some optional quantities are only available when the corresponding compile-time
option or physics package is enabled; those restrictions are noted in the table.

Subvolumes current default to plotting only the three velocity components but will
be generalized in future.

Output Options for 3D plotfiles
-------------------------------

+-----------------------------+------------------+
| Parameter                   | Definition       |
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
| **moist_density**           | Total density    |
|                             | [kg/m^3]         |
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
| **KE**                      | SGS turbulent    |
|                             | kinetic energy   |
|                             | (from Deardorff  |
|                             | or MYNN)         |
|                             | [m^2/s^2]        |
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

Windfarm-only 3D plotfile variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following quantities are available only in builds with
``ERF_USE_WINDFARM`` enabled.

+-----------------------------+------------------+
| Parameter                   | Definition       |
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

Morrison Microphysics Output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When using Morrison two-moment microphysics, additional diagnostic variables
are available for output. These variables provide detailed information about
cloud and precipitation processes. To enable Morrison output, include any of
the variables below in your **erf.plot_vars_1** or **erf.plot_vars_2** list.

**Thermodynamic State Variables:**

+-----------------------------+------------------+
| Parameter                   | Definition       |
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
| Parameter                   | Definition       |
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
| Parameter                   | Definition       |
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
| Parameter                   | Definition       |
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
| Parameter                   | Definition       |
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

Output Options for 2D Plotfiles
-------------------------------

.. toctree::
   :maxdepth: 1

   LandSurfaceDiagnostics.rst

ERF supports two independent 2-D plotfile streams. Each output level is a
one-cell-thick horizontal slab. Use ``erf.plot2d_vars_1`` and
``erf.plot2d_vars_2`` to request built-in diagnostics. Use
``erf.plot2d_level_sets_1`` and ``erf.plot2d_level_sets_2`` to request
generated sampled-level fields.

ERF writes selected fixed diagnostics in canonical catalog order, not input
order. It appends sampled-level fields in level-set, target-value, and field
order. ERF warns and skips an unknown or configuration-ineligible request name.
An accepted component follows its documented runtime value contract.

The unified land-surface diagnostics ``temperature_2m``,
``water_vapor_mixing_ratio_2m``, and ``near_surface_diagnostic_source`` are
defined in :ref:`sec:LandSurfaceDiagnostics`.

.. _sec:Plotfile2DBuiltInCatalog:

Built-In 2D Diagnostic Catalog
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Descriptor metadata** means the fixed catalog's name, canonical order, long
name, unit string, category, and missing-value policy. **Selectable** means
that ERF accepts a request name for the active configuration and creates an
output component. **Runtime value** means the value written in that component
at a cell and output time. A **provider sentinel** is an internal invalid
value translated at the public output boundary. A **categorical source code**
is an integer state whose zero value means no source; it is not a continuous
fill value. A descriptor's metadata policy does not guarantee that the name is
selectable in every configuration.

.. BEGIN ERF BUILT-IN 2D DIAGNOSTIC CATALOG

.. list-table::
   :header-rows: 1
   :widths: 30 17 20 38 80

   * - Variable
     - Category
     - Units
     - Metadata policy
     - Description
   * - ``z_surf``
     - ``Geometry``
     - ``m``
     - ``AlwaysAvailable``
     - Surface elevation
   * - ``landmask``
     - ``Geometry``
     - ``1``
     - ``AlwaysAvailable``
     - Land-sea mask
   * - ``mapfac``
     - ``Geometry``
     - ``1``
     - ``AlwaysAvailable``
     - Map factor at mass points
   * - ``lat_m``
     - ``Geometry``
     - ``deg``
     - ``FillZeroWhenUnavailable``
     - Latitude at unstaggered mass points
   * - ``lon_m``
     - ``Geometry``
     - ``deg``
     - ``FillZeroWhenUnavailable``
     - Longitude at unstaggered mass points
   * - ``u_star``
     - ``SurfaceLayer``
     - ``m/s``
     - ``FillMinus999WhenUnavailable``
     - Friction velocity from the surface layer
   * - ``w_star``
     - ``SurfaceLayer``
     - ``m/s``
     - ``FillMinus999WhenUnavailable``
     - Convective velocity scale from the surface layer
   * - ``t_star``
     - ``SurfaceLayer``
     - ``K``
     - ``FillMinus999WhenUnavailable``
     - Temperature scale from the surface layer
   * - ``q_star``
     - ``SurfaceLayer``
     - ``kg/kg``
     - ``FillMinus999WhenUnavailable``
     - Humidity scale from the surface layer
   * - ``Olen``
     - ``SurfaceLayer``
     - ``m``
     - ``FillMinus999WhenUnavailable``
     - Obukhov length from the surface layer
   * - ``pblh``
     - ``SurfaceLayer``
     - ``m``
     - ``FillMinus999WhenUnavailable``
     - Planetary boundary layer height from the active PBL diagnostic provider
   * - ``t_surf``
     - ``SurfaceLayer``
     - ``K``
     - ``FillMinus999WhenUnavailable``
     - Surface temperature from the surface layer
   * - ``q_surf``
     - ``SurfaceLayer``
     - ``kg/kg``
     - ``FillMinus999WhenUnavailable``
     - Surface humidity from the surface layer
   * - ``z0``
     - ``SurfaceLayer``
     - ``m``
     - ``FillMinus999WhenUnavailable``
     - Roughness height from the surface layer
   * - ``OLR``
     - ``Radiation``
     - ``W/m^2``
     - ``FillMinus999WhenUnavailable``
     - Outgoing longwave radiation at the model top
   * - ``sens_flux``
     - ``SurfaceFlux``
     - ``kg K m^-2 s^-1``
     - ``FillMinus999WhenUnavailable``
     - Surface sensible heat flux
   * - ``laten_flux``
     - ``SurfaceFlux``
     - ``kg m^-2 s^-1``
     - ``FillMinus999WhenUnavailable``
     - Surface moisture flux (legacy output name)
   * - ``surf_pres``
     - ``SurfaceState``
     - ``Pa``
     - ``AlwaysAvailable``
     - Surface pressure
   * - ``precip_total_accum``
     - ``Precipitation``
     - ``kg/m^2``
     - ``FillZeroWhenUnavailable``
     - Accumulated surface precipitation, liquid-water equivalent
   * - ``precip_rain_accum``
     - ``Precipitation``
     - ``kg/m^2``
     - ``FillZeroWhenUnavailable``
     - Accumulated surface rain precipitation, liquid-water equivalent
   * - ``precip_snow_accum``
     - ``Precipitation``
     - ``kg/m^2``
     - ``FillZeroWhenUnavailable``
     - Accumulated surface snow precipitation, liquid-water equivalent
   * - ``precip_graupel_accum``
     - ``Precipitation``
     - ``kg/m^2``
     - ``FillZeroWhenUnavailable``
     - Accumulated surface graupel precipitation, liquid-water equivalent
   * - ``precip_hail_accum``
     - ``Precipitation``
     - ``kg/m^2``
     - ``FillZeroWhenUnavailable``
     - Accumulated surface hail precipitation, liquid-water equivalent
   * - ``precip_frozen_accum``
     - ``Precipitation``
     - ``kg/m^2``
     - ``FillZeroWhenUnavailable``
     - Accumulated frozen surface precipitation, liquid-water equivalent
   * - ``integrated_qv``
     - ``ColumnIntegral``
     - ``kg/m^2``
     - ``FillZeroWhenUnavailable``
     - Column-integrated water vapor
   * - ``integrated_qc``
     - ``ColumnIntegral``
     - ``kg/m^2``
     - ``FillZeroWhenUnavailable``
     - Column-integrated cloud liquid water
   * - ``integrated_qi``
     - ``ColumnIntegral``
     - ``kg/m^2``
     - ``FillZeroWhenUnavailable``
     - Column-integrated cloud ice
   * - ``integrated_qr``
     - ``ColumnIntegral``
     - ``kg/m^2``
     - ``FillZeroWhenUnavailable``
     - Column-integrated rain water
   * - ``integrated_qs``
     - ``ColumnIntegral``
     - ``kg/m^2``
     - ``FillZeroWhenUnavailable``
     - Column-integrated snow
   * - ``integrated_qg``
     - ``ColumnIntegral``
     - ``kg/m^2``
     - ``FillZeroWhenUnavailable``
     - Column-integrated graupel
   * - ``surface_diagnostic_source``
     - ``SurfaceLayer``
     - ``1``
     - ``AlwaysAvailable``
     - Surface diagnostic source code
   * - ``sensible_heat_flux``
     - ``SurfaceFlux``
     - ``W m^-2``
     - ``FillMinus999WhenUnavailable``
     - Surface sensible heat flux
   * - ``latent_heat_flux``
     - ``SurfaceFlux``
     - ``W m^-2``
     - ``FillMinus999WhenUnavailable``
     - Surface latent heat flux
   * - ``shoc_u_star``
     - ``PBL``
     - ``m/s``
     - ``FillMinus999WhenUnavailable``
     - Native SHOC friction velocity diagnostic
   * - ``shoc_Olen``
     - ``PBL``
     - ``m``
     - ``FillMinus999WhenUnavailable``
     - Native SHOC Obukhov length diagnostic
   * - ``shoc_wthv_sfc``
     - ``PBL``
     - ``K m s^-1``
     - ``FillMinus999WhenUnavailable``
     - Native SHOC surface virtual potential temperature flux
   * - ``t_sfc``
     - ``LandSurface``
     - ``K``
     - ``FillMinus999WhenUnavailable``
     - Noah-MP radiative surface temperature
   * - ``sfc_emis``
     - ``LandSurface``
     - ``1``
     - ``FillMinus999WhenUnavailable``
     - Surface bulk emissivity
   * - ``sfc_alb_dir_vis``
     - ``LandSurface``
     - ``1``
     - ``FillMinus999WhenUnavailable``
     - Direct visible surface albedo
   * - ``sfc_alb_dir_nir``
     - ``LandSurface``
     - ``1``
     - ``FillMinus999WhenUnavailable``
     - Direct near-infrared surface albedo
   * - ``sfc_alb_dif_vis``
     - ``LandSurface``
     - ``1``
     - ``FillMinus999WhenUnavailable``
     - Diffuse visible surface albedo
   * - ``sfc_alb_dif_nir``
     - ``LandSurface``
     - ``1``
     - ``FillMinus999WhenUnavailable``
     - Diffuse near-infrared surface albedo
   * - ``cos_zenith_angle``
     - ``LandSurface``
     - ``1``
     - ``FillMinus999WhenUnavailable``
     - Cosine of the solar zenith angle supplied to Noah-MP
   * - ``sw_flux_dn``
     - ``LandSurface``
     - ``W m^-2``
     - ``FillMinus999WhenUnavailable``
     - Downwelling shortwave flux supplied to Noah-MP
   * - ``sw_flux_dn_dir_vis``
     - ``LandSurface``
     - ``W m^-2``
     - ``FillMinus999WhenUnavailable``
     - Direct visible downwelling shortwave flux
   * - ``sw_flux_dn_dir_nir``
     - ``LandSurface``
     - ``W m^-2``
     - ``FillMinus999WhenUnavailable``
     - Direct near-infrared downwelling shortwave flux
   * - ``sw_flux_dn_dif_vis``
     - ``LandSurface``
     - ``W m^-2``
     - ``FillMinus999WhenUnavailable``
     - Diffuse visible downwelling shortwave flux
   * - ``sw_flux_dn_dif_nir``
     - ``LandSurface``
     - ``W m^-2``
     - ``FillMinus999WhenUnavailable``
     - Diffuse near-infrared downwelling shortwave flux
   * - ``lw_flux_dn``
     - ``LandSurface``
     - ``W m^-2``
     - ``FillMinus999WhenUnavailable``
     - Downwelling longwave flux supplied to Noah-MP
   * - ``grdflx``
     - ``LandSurface``
     - ``W m^-2``
     - ``FillMinus999WhenUnavailable``
     - Ground or snow heat flux
   * - ``fira``
     - ``LandSurface``
     - ``W m^-2``
     - ``FillMinus999WhenUnavailable``
     - Total net longwave flux, positive toward the atmosphere
   * - ``sav``
     - ``LandSurface``
     - ``W m^-2``
     - ``FillMinus999WhenUnavailable``
     - Solar radiation absorbed by vegetation
   * - ``sag``
     - ``LandSurface``
     - ``W m^-2``
     - ``FillMinus999WhenUnavailable``
     - Solar radiation absorbed by the ground
   * - ``albedo``
     - ``LandSurface``
     - ``1``
     - ``FillMinus999WhenUnavailable``
     - Broadband surface albedo
   * - ``sfcrunoff``
     - ``LandSurface``
     - ``m``
     - ``FillMinus999WhenUnavailable``
     - Accumulated surface runoff
   * - ``udrunoff``
     - ``LandSurface``
     - ``m``
     - ``FillMinus999WhenUnavailable``
     - Accumulated subsurface runoff
   * - ``noahmp_temperature_2m_vegetated``
     - ``LandSurface``
     - ``K``
     - ``FillMinus999WhenUnavailable``
     - Noah-MP 2-m temperature over the vegetated fraction
   * - ``noahmp_temperature_2m_bare``
     - ``LandSurface``
     - ``K``
     - ``FillMinus999WhenUnavailable``
     - Noah-MP 2-m temperature over the bare fraction
   * - ``noahmp_water_vapor_mixing_ratio_2m_vegetated``
     - ``LandSurface``
     - ``kg kg^-1 dry air``
     - ``FillMinus999WhenUnavailable``
     - Noah-MP 2-m water-vapor mixing ratio over the vegetated fraction
   * - ``noahmp_water_vapor_mixing_ratio_2m_bare``
     - ``LandSurface``
     - ``kg kg^-1 dry air``
     - ``FillMinus999WhenUnavailable``
     - Noah-MP 2-m water-vapor mixing ratio over the bare fraction
   * - ``noahmp_vegetation_fraction``
     - ``LandSurface``
     - ``1``
     - ``FillMinus999WhenUnavailable``
     - Noah-MP vegetation fraction
   * - ``temperature_2m``
     - ``LandSurface``
     - ``K``
     - ``FillMinus999WhenUnavailable``
     - Physical air temperature 2 m above the local surface
   * - ``water_vapor_mixing_ratio_2m``
     - ``LandSurface``
     - ``kg kg^-1 dry air``
     - ``FillMinus999WhenUnavailable``
     - Water-vapor mixing ratio 2 m above the local surface per unit dry-air mass
   * - ``near_surface_diagnostic_source``
     - ``LandSurface``
     - ``1``
     - ``AlwaysAvailable``
     - Source code for the unified near-surface diagnostic bundle
.. END ERF BUILT-IN 2D DIAGNOSTIC CATALOG

.. note::

   For the ``landmask`` variable, land is ``1`` and sea is ``0``. Buildings are ``2`` when using ImmersedForcing.

If a requested fixed diagnostic is not selectable for the active configuration,
ERF warns and skips it. The descriptor table above records metadata; it does
not define request acceptance or guarantee a non-missing runtime value.

.. _sec:Plotfile2DSelectionRules:

Selection and runtime-value rules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The selection contract and the value written after selection are separate:

.. list-table::
   :header-rows: 1
   :widths: 30 35 55

   * - Name or family
     - Selectable when
     - Runtime value
   * - ``z_surf``, ``landmask``, ``mapfac``
     - Selectable: fixed geometry or state writer path.
     - Value: ERF writes the corresponding geometry or state value.
   * - ``lat_m``, ``lon_m``
     - Selectable: fixed request names.
     - Value: ERF writes the coordinate data, or zero when the coordinate source is absent.
   * - ``u_star``, ``w_star``, ``t_star``, ``q_star``, ``Olen``, ``t_surf``, ``q_surf``, ``z0``
     - Selectable: these names are not filtered by the catalog builder.
     - Value: SurfaceLayer scalar values; ``-999`` when the source pointer is absent.
   * - ``pblh``
     - Selectable: fixed request name.
     - Value: native SHOC PBL height when present, otherwise SurfaceLayer; ``-999`` if neither exists.
   * - ``OLR``
     - Selectable: fixed request name.
     - Value: radiation output; ``-999`` when the radiation source is absent.
   * - ``sens_flux``, ``laten_flux``
     - Selectable: fixed request names.
     - Value: legacy conservative surface flux outputs.
   * - ``surf_pres``
     - Selectable: fixed request name.
     - Value: pressure computed from the lowest atmospheric state.
   * - ``surface_diagnostic_source``
     - Selectable: fixed request name.
     - Value: categorical code ``0`` through ``6``; code ``0`` means no source, including when SurfaceLayer is absent.
   * - ``sensible_heat_flux``, ``latent_heat_flux``
     - Selectable: fixed request names.
     - Value: the selected conservative sources converted to energy-flux units; ``-999`` when unavailable.
   * - ``shoc_u_star``, ``shoc_Olen``, ``shoc_wthv_sfc``
     - Selectable: fixed request names.
     - Value: native SHOC diagnostics; ``-999`` when absent.
   * - ``precip_total_accum``, ``precip_frozen_accum``
     - Selectable: the active moisture scheme has any rain, snow, or graupel mass component.
     - Value: normalized liquid-water-equivalent accumulations; ``precip_frozen_accum`` is defined zero in a rain-only scheme.
   * - ``precip_rain_accum``
     - Selectable: the active moisture scheme has a rain mass component.
     - Value: normalized rain accumulation, or the derived rain value when a runtime total source is used without a direct rain source.
   * - ``precip_snow_accum``
     - Selectable: the active moisture scheme has a snow mass component.
     - Value: normalized snow accumulation.
   * - ``precip_graupel_accum``
     - Selectable: the active moisture scheme has a graupel mass component.
     - Value: normalized graupel accumulation.
   * - ``precip_hail_accum``
     - Selectable: never for current moisture schemes; it remains a fixed descriptor.
     - Value: zero through the fixed metadata policy unless a future scheme supplies a hail source.
   * - ``integrated_qv``
     - Selectable: always.
     - Value: column-integrated water vapor; zero when moisture is disabled.
   * - ``integrated_qc``, ``integrated_qi``, ``integrated_qr``, ``integrated_qs``, ``integrated_qg``
     - Selectable: the active moisture model exposes the corresponding conserved mass component.
     - Value: the corresponding column water path. An unsupported species is skipped during request validation.
   * - Fixed land-surface provider names from ``t_sfc`` through ``noahmp_vegetation_fraction``
     - Selectable: the active land-surface inventory contains the exact name.
     - Value: finite provider values pass through when valid; absent, nonfinite, or sentinel values become ``-999``. Valid zero and signed values pass through when the field permits them. Raw fields never use MOST fallback.
   * - ``temperature_2m``
     - Selectable: a Noah-MP or SurfaceLayer pathway exists.
     - Value: a coherent native or MOST value, or ``-999`` when no complete source satisfies the request.
   * - ``water_vapor_mixing_ratio_2m``
     - Selectable: moisture is enabled and a Noah-MP or SurfaceLayer pathway exists.
     - Value: a coherent native or MOST mixing ratio, or ``-999`` when no complete source satisfies the request.
   * - ``near_surface_diagnostic_source``
     - Selectable: a Noah-MP or SurfaceLayer pathway exists, including a dry source-only request.
     - Value: categorical code ``0`` through ``6``; its ``AlwaysAvailable`` metadata policy records ``missing_value: null``.

Precipitation selections follow the active moisture components. At runtime,

.. math::

   P_{frozen} = P_{snow} + P_{graupel} + P_{hail},

.. math::

   P_{total} =
   \begin{cases}
   P_{total}^{native}, & \text{when a native total source exists},\\
   P_{rain} + P_{frozen}, & \text{otherwise},
   \end{cases}

and, when rain is not supplied directly,

.. math::

   P_{rain} = \max\left(0, P_{total} - P_{frozen}\right).

All public precipitation accumulations use liquid-water-equivalent ``kg/m^2``
metadata units.

.. _sec:Plotfile2DDynamicSoil:

Dynamic soil diagnostic families
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Dynamic soil fields are not part of the fixed catalog. ERF exposes only names
reported by the active land-surface provider inventory, in that inventory's
runtime order. ``<layer>`` is a one-based layer index.

.. BEGIN ERF DYNAMIC SOIL DIAGNOSTIC FAMILIES

.. list-table::
   :header-rows: 1
   :widths: 28 42 20 35

   * - Name family
     - Meaning
     - Units
     - Metadata
   * - ``smois_<layer>``
     - Volumetric total soil moisture
     - ``m^3 m^-3``
     - Category ``LandSurface``; policy ``FillMinus999WhenUnavailable``; invalid public values become ``-999``.
   * - ``sh2o_<layer>``
     - Volumetric liquid soil water
     - ``m^3 m^-3``
     - Category ``LandSurface``; policy ``FillMinus999WhenUnavailable``; invalid public values become ``-999``.
   * - ``tslb_<layer>``
     - Soil temperature
     - ``K``
     - Category ``LandSurface``; policy ``FillMinus999WhenUnavailable``; invalid public values become ``-999``.

.. END ERF DYNAMIC SOIL DIAGNOSTIC FAMILIES

Flux Diagnostics
^^^^^^^^^^^^^^^^

``sens_flux`` and ``laten_flux`` are legacy ERF conservative scalar flux
outputs. ``sensible_heat_flux`` and ``latent_heat_flux`` convert the same
selected conservative flux sources to ``W m^-2`` using ``Cp_d`` and ``L_v``,
respectively.

For non-SHOC configurations and native SHOC host-diffusion mode, these outputs
use the host vertical surface flux arrays. In native SHOC ``state_update`` mode,
SHOC consumes those surface fluxes before the host diffusion path clears the
overlapping arrays. In that mode, the 2D flux diagnostics use SHOC's preserved
consumed-flux snapshots, component by component. If the corresponding host flux
field was unavailable before SHOC consumed it, ERF writes ``-999`` rather than a
zero SHOC snapshot.

The sign convention follows ERF's lower-boundary flux convention. No
source-specific Noah-MP, MOST, or SHOC conversion is applied in the 2D
diagnostic layer. Noah-MP LSM kinematic-to-conservative conversion happens
upstream in SurfaceLayer. Native SHOC converts the conservative host fluxes to
kinematic column inputs internally, then the diagnostic snapshot preserves the
consumed flux in conservative ERF output units.

In native SHOC ``state_update`` mode, SHOC is the transport owner. The
``surface_diagnostic_source`` field still describes the upstream surface-flux
source path used before SHOC consumes it.

2D AMReX Metadata Sidecar
^^^^^^^^^^^^^^^^^^^^^^^^^

Native AMReX 2-D plotfiles write a JSON metadata sidecar named
``2DMetadata.json`` in the plotfile directory. The sidecar lists selected
outputs in component order. Fixed outputs use descriptor metadata. Sampled-
level outputs add source-field and vertical-coordinate metadata.

The metadata policy maps to JSON as follows:

* ``AlwaysAvailable`` records ``missing_value: null``.
* ``FillZeroWhenUnavailable`` records ``missing_value: 0``.
* ``FillMinus999WhenUnavailable`` records ``missing_value: -999``.

The sidecar uses the same 2D diagnostic catalog that defines the built-in
plotfile variables. Sampled-level outputs carry their own metadata record. The
sidecar does not record the runtime source chosen at each cell. Source choice is
carried by categorical output fields. It does not change field values.

NetCDF 2-D output uses the same variable names but does not write this JSON
sidecar or sampled-level metadata attributes. The sidecar format version is
``2``.

Example built-in metadata:

.. code-block:: json

   {
     "format_version": 2,
     "kind": "ERF 2D plotfile metadata",
     "n_variables": 2,
     "variables": [
       {
         "component_index": 0,
         "name": "z_surf",
         "long_name": "Surface elevation",
         "units": "m",
         "category": "Geometry",
         "missing_policy": "AlwaysAvailable",
         "missing_value": null
       },
       {
         "component_index": 1,
         "name": "latent_heat_flux",
         "long_name": "Surface latent heat flux",
         "units": "W m^-2",
         "category": "SurfaceFlux",
         "missing_policy": "FillMinus999WhenUnavailable",
         "missing_value": -999
       }
     ]
   }

Example sampled-level metadata:

.. code-block:: json

   {
     "format_version": 2,
     "kind": "ERF 2D plotfile metadata",
     "n_variables": 1,
     "variables": [
       {
         "component_index": 0,
         "name": "theta_p_850hPa",
         "long_name": "Potential temperature sampled on pressure levels",
         "units": "K",
         "category": "SampledLevel",
         "missing_policy": "FillMinus999WhenUnavailable",
         "missing_value": -999,
         "source_field": "theta",
         "vertical_coordinate": {
           "type": "pressure",
           "value": 850,
           "units": "hPa",
           "canonical_value": 85000,
           "canonical_units": "Pa",
           "interpolation": "linear"
         }
       }
     ]
   }

2D Sampled-Level Diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ERF can write 2D fields sampled on model-index, height, or pressure levels. This
mode lets one 2D output stream combine built-in diagnostics with fields on
user-defined vertical targets.

Select sampled level sets for each 2D output stream with:

.. code-block:: text

   erf.plot2d_level_sets_1 = upper_air bl_heights
   erf.plot2d_level_sets_2 = native_k

Define each level set under ``erf.plot2d.level_set.<name>.``.

Required keys:

.. code-block:: text

   erf.plot2d.level_set.<name>.coordinate = ...
   erf.plot2d.level_set.<name>.values = ...
   erf.plot2d.level_set.<name>.fields = ...

Optional keys:

.. code-block:: text

   erf.plot2d.level_set.<name>.units = ...
   erf.plot2d.level_set.<name>.interpolation = ...
   erf.plot2d.level_set.<name>.missing_value = ...

The default missing value is ``-999``. ERF does not extrapolate sampled-level
diagnostics. If a target lies outside the column, ERF writes the level set's
missing value.

Supported coordinates:

.. list-table::
   :header-rows: 1
   :widths: 20 48 16 16

   * - Coordinate
     - Meaning
     - Units
     - Interpolation
   * - ``model_index``
     - Cell-centered model level index.
     - ``1``
     - ``none``
   * - ``height_msl``
     - Cell-centered height above mean sea level.
     - ``m``
     - ``linear``
   * - ``height_agl``
     - Cell-centered height above local terrain.
     - ``m``
     - ``linear``
   * - ``pressure``
     - Cell-centered pressure.
     - ``Pa`` or ``hPa``
     - ``linear``

ERF converts pressure targets in ``hPa`` to canonical ``Pa`` values in metadata.
Model-index values must be integers. ERF rejects unsupported coordinates and
units during input parsing. Isentropic output is not enabled because it needs a
crossing policy for non-monotonic potential-temperature columns.

Supported fields:

.. list-table::
   :header-rows: 1
   :widths: 18 52 14 26

   * - Field
     - Meaning
     - Units
     - Availability
   * - ``rho``
     - Density.
     - ``kg m^-3``
     - Always available.
   * - ``theta``
     - Potential temperature.
     - ``K``
     - Always available.
   * - ``temp``
     - Temperature.
     - ``K``
     - Always available.
   * - ``pressure``
     - Pressure.
     - ``Pa``
     - Always available.
   * - ``height_msl``
     - Height above mean sea level.
     - ``m``
     - Always available.
   * - ``height_agl``
     - Height above local terrain.
     - ``m``
     - Always available.
   * - ``qv``
     - Water vapor mixing ratio.
     - ``kg kg^-1``
     - Available when the active moisture model has ``qv``.
   * - ``qc``
     - Cloud liquid water mixing ratio.
     - ``kg kg^-1``
     - Available when the active moisture model has ``qc``.
   * - ``qi``
     - Cloud ice mixing ratio.
     - ``kg kg^-1``
     - Available when the active moisture model has ``qi``.
   * - ``qr``
     - Rain mixing ratio.
     - ``kg kg^-1``
     - Available when the active moisture model has ``qr``.
   * - ``qs``
     - Snow mixing ratio.
     - ``kg kg^-1``
     - Available when the active moisture model has ``qs``.
   * - ``qg``
     - Graupel mixing ratio.
     - ``kg kg^-1``
     - Available when the active moisture model has ``qg``.
   * - ``u_east``
     - Eastward wind.
     - ``m/s``
     - Always available.
   * - ``v_north``
     - Northward wind.
     - ``m/s``
     - Always available.
   * - ``w``
     - Vertical wind.
     - ``m/s``
     - Always available.
   * - ``wind_speed``
     - Horizontal wind speed.
     - ``m/s``
     - Always available.
   * - ``wind_dir``
     - Meteorological wind direction.
     - ``degrees``
     - Always available. ERF writes the missing value for calm winds.

Use the canonical sampled species names ``qv``, ``qc``, ``qi``, ``qr``, ``qs``,
and ``qg``. Do not use 3D derived-variable aliases such as ``qrain``, ``qsnow``,
or ``qgraup`` in sampled-level output.

If a requested moisture field is unavailable for the active moisture model, ERF
warns and skips that field. If all requested fields in a level set are
unavailable, ERF aborts with an input error.

Pressure-level example:

.. code-block:: text

   erf.plot2d_level_sets_1 = upper_air

   erf.plot2d.level_set.upper_air.coordinate = pressure
   erf.plot2d.level_set.upper_air.units = hPa
   erf.plot2d.level_set.upper_air.values = 850 700 500
   erf.plot2d.level_set.upper_air.fields = theta temp qv qc

Height-above-ground example:

.. code-block:: text

   erf.plot2d_level_sets_1 = bl_heights

   erf.plot2d.level_set.bl_heights.coordinate = height_agl
   erf.plot2d.level_set.bl_heights.values = 100 500 1000
   erf.plot2d.level_set.bl_heights.fields = theta qv qc

Model-index example:

.. code-block:: text

   erf.plot2d_level_sets_1 = native_k

   erf.plot2d.level_set.native_k.coordinate = model_index
   erf.plot2d.level_set.native_k.values = 0 10 20
   erf.plot2d.level_set.native_k.fields = theta qv qc

Sampled-level output names follow this pattern:

.. code-block:: text

   <field>_<coordinate_tag>_<value_tag>

Examples:

.. code-block:: text

   theta_p_850hPa
   qv_p_700hPa
   qc_z_agl_500m
   theta_z_msl_1000m
   theta_k_10
   u_east_p_850hPa
   v_north_p_850hPa
   w_z_agl_500m
   wind_speed_k_10
   wind_speed_z_agl_500m
   wind_dir_k_10

Wind example:

.. code-block:: text

   erf.plot2d_level_sets_1 = upper_air_winds

   erf.plot2d.level_set.upper_air_winds.coordinate = pressure
   erf.plot2d.level_set.upper_air_winds.units = hPa
   erf.plot2d.level_set.upper_air_winds.values = 850 700
   erf.plot2d.level_set.upper_air_winds.fields = u_east v_north wind_speed wind_dir

ERF writes sampled-level variables after built-in variables in the same 2D
stream. Within sampled-level output, ERF writes variables in level-set order,
target-value order, and field order.

For linear sampling, ERF finds adjacent cell centers :math:`k_0` and
:math:`k_1` that bracket the target coordinate :math:`C_t`.

.. math::

   F_t = (1 - w) F_{k_0} + w F_{k_1}

with

.. math::

   w = \frac{C_t - C_{k_0}}{C_{k_1} - C_{k_0}}.

The bracket test works for increasing coordinates, such as height, and
decreasing coordinates, such as pressure. Model-index sampling copies the
requested level exactly.

Microphysical mass species are sampled as mixing ratios. ERF stores each
species as a conserved density :math:`\rho q_x`, so sampled output uses

.. math::

   q_x = \frac{\rho q_x}{\rho}.

ERF uses ``qv = 0`` for dry-run pressure and temperature calculations, but it
does not expose sampled ``qv`` unless the active moisture model has a water-vapor
state component.

Wind fields are cell-centered sampled-level diagnostics. ERF destaggers the
native face-centered velocity components to scalar cell centers before vertical
sampling. ``u_east`` and ``v_north`` are earth-relative horizontal winds. When
map-rotation coefficients are available, ERF rotates grid-relative horizontal
winds before output. Otherwise, ERF uses identity rotation. ``w`` is the
cell-centered vertical wind.

``wind_speed`` is the horizontal speed computed from ``u_east`` and
``v_north``. ``wind_dir`` is the meteorological wind direction in degrees
clockwise from north, indicating where the wind comes from. ERF writes the
level-set missing value for ``wind_dir`` when the horizontal wind is calm.
``wind_dir`` is derived from the interpolated vector; ERF does not vertically
interpolate wind direction as a scalar angle.

Sampled-level metadata records ``source_field`` and ``vertical_coordinate`` in
the native AMReX ``2DMetadata.json`` sidecar. NetCDF 2D output uses the same
sampled-level variable names, but it does not write sampled-level metadata
attributes.

Isentropic levels, vorticity, potential vorticity, and staggered velocity fields
need additional sampling rules. They are not part of this sampled-level output
mode.

2D Surface Precipitation Accumulation Diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ERF can write cumulative surface precipitation accumulations from the active
microphysics scheme in 2D plotfiles. These diagnostics report the liquid-water
equivalent mass that has reached the lower boundary since model start or the
most recent restart. Some schemes store explicit rain/snow/graupel species
accumulators, while others store a total accumulator plus frozen-species
subsets. ERF normalizes each available scheme-native source to ``kg/m^2``
before deriving the public 2-D fields. Request acceptance follows the active
rain, snow, and graupel mass components: total and frozen require any one of
those components; each species field requires its own component; hail is not
selectable for current schemes. Runtime mapping may still derive rain from a
total source without a direct rain source. The public fields use ``kg/m^2``
because downstream land-surface forcing consumes precipitation as mass over
area, even though ``1 kg/m^2`` is numerically equal to ``1 mm`` of liquid-water
equivalent.

+-----------------------------+--------------------------------------------------+----------+
| Field                       | Meaning                                          | Units    |
+=============================+==================================================+==========+
| ``precip_total_accum``      | Accumulated surface precipitation, liquid-water  | kg/m^2   |
|                             | equivalent                                       |          |
+-----------------------------+--------------------------------------------------+----------+
| ``precip_rain_accum``       | Accumulated surface rain precipitation,          | kg/m^2   |
|                             | liquid-water equivalent                          |          |
+-----------------------------+--------------------------------------------------+----------+
| ``precip_snow_accum``       | Accumulated surface snow precipitation,          | kg/m^2   |
|                             | liquid-water equivalent                          |          |
+-----------------------------+--------------------------------------------------+----------+
| ``precip_graupel_accum``    | Accumulated surface graupel precipitation,       | kg/m^2   |
|                             | liquid-water equivalent                          |          |
+-----------------------------+--------------------------------------------------+----------+
| ``precip_hail_accum``       | Accumulated surface hail precipitation,          | kg/m^2   |
|                             | liquid-water equivalent                          |          |
+-----------------------------+--------------------------------------------------+----------+
| ``precip_frozen_accum``     | Accumulated frozen surface precipitation,        | kg/m^2   |
|                             | liquid-water equivalent                          |          |
+-----------------------------+--------------------------------------------------+----------+

``precip_total_accum`` is the normalized native total when one exists, or the
sum of normalized rain and frozen accumulations otherwise.
``precip_frozen_accum`` is the sum of normalized snow, graupel, and hail
accumulations. Species absent from the active runtime source contribute zero to
derived totals. ``precip_hail_accum`` is a fixed descriptor reserved for a
future scheme that exposes a distinct hail accumulator.

To diagnose the frozen fraction over a coupling interval, use accumulation
differences rather than a ratio of cumulative values:

.. math::

   f_{frozen} =
   \begin{cases}
   \dfrac{\Delta P_{frozen}}{\Delta P_{total}}, & \Delta P_{total} > 0 \\
   0, & \text{otherwise}
   \end{cases}

with

.. math::

   \Delta P_{total} = P_{total}(t_1) - P_{total}(t_0),

and

.. math::

   \Delta P_{frozen} = P_{frozen}(t_1) - P_{frozen}(t_0).

2D Water-Path Diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^

ERF can write scheme-aware 2D water-path diagnostics for prognostic condensed
water species. These diagnostics integrate the conserved ``rho*q`` species
component over the model column. They are available only when the active
microphysics scheme exposes the corresponding conserved mass component.

For a condensed species :math:`q_x`, ERF writes

.. math::

   W_x(x,y) = \int_{z_b}^{z_t} \rho(x,y,z)\,q_x(x,y,z)\,dz,

where :math:`W_x` has units of :math:`\mathrm{kg\,m^{-2}}`.

The discrete diagnostic uses the same metric convention as ``integrated_qv``:

.. math::

   W_x(i,j) = \Delta z \sum_k (\rho q_x)_{i,j,k}

for constant-:math:`\Delta z` meshes, and

.. math::

   W_x(i,j) = \Delta z \sum_k (\rho q_x)_{i,j,k} J_{i,j,k}

when ERF uses the column metric factor :math:`J`.

+-------------------------------+-----------------------------------------------+----------+
| **integrated_qc**             | Cloud liquid water path                       | kg/m^2   |
+-------------------------------+-----------------------------------------------+----------+
| **integrated_qi**             | Cloud ice water path                          | kg/m^2   |
+-------------------------------+-----------------------------------------------+----------+
| **integrated_qr**             | Rain water path                               | kg/m^2   |
+-------------------------------+-----------------------------------------------+----------+
| **integrated_qs**             | Snow water path                               | kg/m^2   |
+-------------------------------+-----------------------------------------------+----------+
| **integrated_qg**             | Graupel water path                            | kg/m^2   |
+-------------------------------+-----------------------------------------------+----------+

A name is selectable only if the active moisture model has that species as a
conserved mass component. Unsupported species are skipped during request
validation. Two-moment number concentrations are not water mass paths and are
not included.

The metadata sidecar records these fields as ``ColumnIntegral`` diagnostics
with ``FillZeroWhenUnavailable`` missing-value policy. Water-path diagnostics
use the built-in metadata fields and add no water-path-specific metadata keys.

.. _sec:Plotfile2DSourceCodes:

Surface Diagnostic Source Codes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``surface_diagnostic_source`` is a cell-centered categorical diagnostic. It
reports the source path used by the SurfaceLayer scalar diagnostic path. It
does not report fractional land-cover contributions. The same code table
defines ``near_surface_diagnostic_source``; that field reports the request-
aware source selected for the unified 2-m output bundle.

If an input dataset contains fractional land information, this diagnostic still
reports the categorical source used by ERF's active SurfaceLayer scalar flux
path.

The fields do not describe fractional cover or complete staggered stress-face
provenance. Consumers must compare the numeric code with this table rather than
infer provenance from ``landmask``.

.. list-table::
   :header-rows: 1
   :widths: 10 32 58

   * - Value
     - Token
     - Meaning
   * - 0
     - ``missing``
     - No source was selected.
   * - 1
     - ``surface_layer_land``
     - SurfaceLayer supplied the land value.
   * - 2
     - ``lsm_land``
     - The land-surface model supplied the land value.
   * - 3
     - ``surface_layer_fallback``
     - A land-surface path existed, but its required result was invalid or incomplete; SurfaceLayer supplied the value.
   * - 4
     - ``surface_layer_sea``
     - SurfaceLayer supplied the value over water.
   * - 5
     - ``custom``
     - The custom surface pathway supplied the SurfaceLayer state.
   * - 6
     - ``rico``
     - The RICO pathway supplied the SurfaceLayer state.

Examples of Usage
-----------------

The following inputs write a 2D plotfile every 10 level-0 time steps:

.. code-block:: none

   erf.plot2d_file_1 = plt2d_
   erf.plot2d_int_1  = 10
   erf.plot2d_vars_1 = z_surf mapfac lat_m lon_m u_star surf_pres integrated_qv

The variable list may appear in any order. ERF writes the selected variables in
its canonical 2D plotfile order.
