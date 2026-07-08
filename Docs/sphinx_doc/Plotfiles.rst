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

ERF supports two 2D plotfile streams. Use ``erf.plot2d_vars_1`` and
``erf.plot2d_vars_2`` to request built-in 2D diagnostics. Use
``erf.plot2d_level_sets_1`` and ``erf.plot2d_level_sets_2`` to request
sampled-level diagnostics.

Built-in variables and sampled-level variables can appear in the same 2D output
stream. ERF writes built-in variables first. It writes selected built-in
variables in the catalog order shown below, not in input order. ERF then appends
sampled-level variables in level-set order, target-value order, and field order.

Built-In 2D Diagnostic Catalog
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 24 16 24 56

   * - Variable
     - Units
     - Availability
     - Description
   * - ``z_surf``
     - ``m``
     - Always available.
     - Surface elevation.
   * - ``landmask``
     - ``1``
     - Always available.
     - Land-sea mask. Land is ``1`` and sea is ``0``.
   * - ``mapfac``
     - ``1``
     - Always available.
     - Map factor at mass points.
   * - ``lat_m``
     - ``deg``
     - Available when latitude data are present.
     - Latitude at unstaggered mass points.
   * - ``lon_m``
     - ``deg``
     - Available when longitude data are present.
     - Longitude at unstaggered mass points.
   * - ``u_star``
     - ``m s^-1``
     - ``-999`` if unavailable.
     - Friction velocity from the surface layer.
   * - ``w_star``
     - ``m s^-1``
     - ``-999`` if unavailable.
     - Convective velocity scale from the surface layer.
   * - ``t_star``
     - ``K``
     - ``-999`` if unavailable.
     - Temperature scale from the surface layer.
   * - ``q_star``
     - ``kg kg^-1``
     - ``-999`` if unavailable.
     - Humidity scale from the surface layer.
   * - ``Olen``
     - ``m``
     - ``-999`` if unavailable.
     - Obukhov length from the surface layer.
   * - ``pblh``
     - ``m``
     - ``-999`` if unavailable.
     - Planetary boundary layer height. Native SHOC provides this value when
       available; otherwise ERF uses SurfaceLayer when present.
   * - ``t_surf``
     - ``K``
     - ``-999`` if unavailable.
     - Surface temperature from the surface layer.
   * - ``q_surf``
     - ``kg kg^-1``
     - ``-999`` if unavailable.
     - Surface humidity from the surface layer.
   * - ``z0``
     - ``m``
     - ``-999`` if unavailable.
     - Roughness height from the surface layer.
   * - ``OLR``
     - ``W m^-2``
     - ``-999`` if unavailable.
     - Outgoing longwave radiation at the model top.
   * - ``sens_flux``
     - ``kg K m^-2 s^-1``
     - ``-999`` if unavailable.
     - Conservative surface sensible heat flux.
   * - ``laten_flux``
     - ``kg m^-2 s^-1``
     - ``-999`` if unavailable.
     - Conservative surface moisture flux. This is a legacy output name.
   * - ``surf_pres``
     - ``Pa``
     - Always available.
     - Surface pressure.
   * - ``integrated_qv``
     - ``kg m^-2``
     - Zero when moisture is disabled.
     - Column-integrated water vapor.
   * - ``integrated_qc``
     - ``kg m^-2``
     - Available when the active moisture model has ``qc``.
     - Column-integrated cloud liquid water.
   * - ``integrated_qi``
     - ``kg m^-2``
     - Available when the active moisture model has ``qi``.
     - Column-integrated cloud ice.
   * - ``integrated_qr``
     - ``kg m^-2``
     - Available when the active moisture model has ``qr``.
     - Column-integrated rain water.
   * - ``integrated_qs``
     - ``kg m^-2``
     - Available when the active moisture model has ``qs``.
     - Column-integrated snow.
   * - ``integrated_qg``
     - ``kg m^-2``
     - Available when the active moisture model has ``qg``.
     - Column-integrated graupel.
   * - ``surface_diagnostic_source``
     - ``1``
     - ``-999`` when SurfaceLayer is not active.
     - Source code for the cell-centered SurfaceLayer scalar diagnostic path.
   * - ``sensible_heat_flux``
     - ``W m^-2``
     - ``-999`` if unavailable.
     - Surface sensible heat flux computed from the same conservative source as
       ``sens_flux``.
   * - ``latent_heat_flux``
     - ``W m^-2``
     - ``-999`` if unavailable.
     - Surface latent heat flux computed from the same conservative source as
       ``laten_flux``.
   * - ``shoc_u_star``
     - ``m s^-1``
     - ``-999`` if unavailable.
     - Native SHOC friction velocity diagnostic.
   * - ``shoc_Olen``
     - ``m``
     - ``-999`` if unavailable.
     - Native SHOC Obukhov length diagnostic.
   * - ``shoc_wthv_sfc``
     - ``K m s^-1``
     - ``-999`` if unavailable.
     - Native SHOC surface virtual potential temperature flux.

If a requested built-in diagnostic is not available for the active configuration,
ERF warns and skips that diagnostic.

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

Native AMReX 2D plotfiles include a JSON metadata sidecar named
``2DMetadata.json`` in the plotfile directory. The sidecar lists the selected
2D variables in output component order. For each variable it records the
component index, name, long name, units, diagnostic category, missing-value
policy, and documented missing value. Sampled-level outputs add a source field
and vertical-coordinate record.

The sidecar uses the same 2D diagnostic catalog that defines the built-in
plotfile variables. Sampled-level outputs carry their own metadata record. The
sidecar does not change field values. For example, in native SHOC
``state_update`` mode the flux fields may come from preserved SHOC-consumed
flux snapshots, as described above, but the sidecar records the public
diagnostic metadata for the selected variables.

Native AMReX 2D plotfiles write sampled-level metadata in ``2DMetadata.json``.
NetCDF 2D output uses the same sampled-level variable names, but does not write
sampled-level metadata attributes. The sidecar format version is ``2``.

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
before deriving the public 2D fields. When a scheme provides a total source
but no explicit rain source, ERF computes ``precip_rain_accum`` as
``precip_total_accum - precip_frozen_accum`` and clips small negative residuals
to zero. The public fields use ``kg/m^2`` because downstream land surface
forcing consumes precipitation as mass over area, even though ``1 kg/m^2`` is
numerically equal to ``1 mm`` of liquid water equivalent.

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

``precip_total_accum`` is the normalized total accumulation when the scheme
stores one, or the sum of the normalized species accumulations otherwise.
``precip_frozen_accum`` is the sum of the normalized frozen species
accumulations. Species that are not available in the active microphysics
scheme contribute zero to the derived totals. ``precip_hail_accum`` is reserved
for a future scheme that exposes a distinct hail accumulator.

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

A name is available only if the active moisture model has that species as a
conserved mass component. Two-moment number concentrations are not water mass
paths and are not included.

The metadata sidecar records these fields as ``ColumnIntegral`` diagnostics
with ``FillZeroWhenUnavailable`` missing-value policy. Water-path diagnostics
use the built-in metadata fields and add no water-path-specific metadata keys.

Surface Diagnostic Source Codes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``surface_diagnostic_source`` is a cell-centered categorical diagnostic. It
reports the source path used by the SurfaceLayer scalar diagnostic path. It
does not report fractional land-cover contributions.

If an input dataset contains fractional land information, this diagnostic still
reports the categorical source used by ERF's active SurfaceLayer scalar flux
path.

The diagnostic does not fully describe staggered stress-face provenance. ERF
may average adjacent LSM and non-LSM contributions when applying staggered
surface stresses.

+------+---------------------------------------------------------------+
| Code | Meaning                                                       |
+======+===============================================================+
| -999 | SurfaceLayer is inactive or the field is unavailable.         |
+------+---------------------------------------------------------------+
| 0    | Missing or unset.                                             |
+------+---------------------------------------------------------------+
| 1    | Non-LSM SurfaceLayer over land.                               |
+------+---------------------------------------------------------------+
| 2    | LSM over land.                                                |
+------+---------------------------------------------------------------+
| 3    | Non-LSM SurfaceLayer fallback where an LSM flux was undefined |
|      | or unavailable for a land cell.                               |
+------+---------------------------------------------------------------+
| 4    | Non-LSM SurfaceLayer over sea.                                |
+------+---------------------------------------------------------------+
| 5    | Custom prescribed surface-layer values.                       |
+------+---------------------------------------------------------------+
| 6    | RICO prescribed surface-layer values.                         |
+------+---------------------------------------------------------------+

Examples of Usage
-----------------

The following inputs write a 2D plotfile every 10 level-0 time steps:

.. code-block:: none

   erf.plot2d_file_1 = plt2d_
   erf.plot2d_int_1  = 10
   erf.plot2d_vars_1 = z_surf mapfac lat_m lon_m u_star surf_pres integrated_qv

The variable list may appear in any order. ERF writes the selected variables in
its canonical 2D plotfile order.
