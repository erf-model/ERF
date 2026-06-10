.. role:: cpp(code)
  :language: c++

.. _sec:Plotfiles:

*********
Plotfiles
*********
.. toctree::
   :maxdepth: 1

There are three types of plotfiles that can be written from ERF.

The first is the standard type of plotfile which includes 3D data on all levels for those variables
specified by the user in the inputs file.

The second is a pseudo-2D plotfile that contains data that is only defined as a function
of horizontal position, such as map factors, latitude and longitude.

The third type of plotfile contains
3D data on one level only and in a specified region of the domain.   We refer to this
latter capability as "Subvolumes" below.  The level at which the data is written
out is determined by the mesh spacing specified by the user.

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
|                                  | format           | "netcdf / "NetCDF"    |            |
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

+----------------------------------+------------------+-----------------------+------------+
| Parameter                        | Definition       | Acceptable            | Default    |
|                                  |                  | Values                |            |
+==================================+==================+=======================+============+
| **erf.plot_file_1**              | prefix for       | String                | “*plt_1_*” |
|                                  | plotfiles        |                       |            |
|                                  | at first freq.   |                       |            |
+----------------------------------+------------------+-----------------------+------------+
| **erf.plot_file_2**              | prefix for       | String                | “*plt_2_*” |
|                                  | plotfiles        |                       |            |
|                                  | at seoncd freq.  |                       |            |
+----------------------------------+------------------+-----------------------+------------+
| **erf.plot_int_1**               | how often (by    | Integer               | -1         |
|                                  | level-0 time     | :math:`> 0`           |            |
|                                  | steps) to write  |                       |            |
|                                  | plot files       |                       |            |
|                                  | at first freq.   |                       |            |
+----------------------------------+------------------+-----------------------+------------+
| **erf.plot_int_2**               | how often (by    | Integer               | -1         |
|                                  | level-0 time     | :math:`> 0`           |            |
|                                  | steps) to write  |                       |            |
|                                  | plot files       |                       |            |
|                                  | at second freq.  |                       |            |
+----------------------------------+------------------+-----------------------+------------+
| **erf.plot_per_1**               | how often (in    | Real                  | -1.0       |
|                                  | simulation time) | :math:`> 0`           |            |
|                                  | to write         |                       |            |
|                                  | plot files       |                       |            |
|                                  | at first freq.   |                       |            |
+----------------------------------+------------------+-----------------------+------------+
| **erf.plot_per_2**               | how often (in    | Real                  | -1.0       |
|                                  | simulation time) | :math:`> 0`           |            |
|                                  | to write         |                       |            |
|                                  | plot files       |                       |            |
|                                  | at second freq.  |                       |            |
+----------------------------------+------------------+-----------------------+------------+
| **erf.plot_vars_1**              | name of          | list of names         | None       |
|                                  | variables to     |                       |            |
|                                  | include in       |                       |            |
|                                  | plotfiles        |                       |            |
|                                  | at first freq.   |                       |            |
+----------------------------------+------------------+-----------------------+------------+
| **erf.plot_vars_2**              | name of          | list of names         | None       |
|                                  | variables to     |                       |            |
|                                  | include in       |                       |            |
|                                  | plotfiles        |                       |            |
|                                  | at second freq.  |                       |            |
+----------------------------------+------------------+-----------------------+------------+
| **erf.plot_face_vels**           | output plotfiles | Boolean               | false      |
|                                  | "{prefix}U",     |                       |            |
|                                  | "{prefix}V", and |                       |            |
|                                  | "{prefix}W"      |                       |            |
|                                  | with velocity    |                       |            |
|                                  | components on the|                       |            |
|                                  | staggered grid.  |                       |            |
+----------------------------------+------------------+-----------------------+------------+

List of Parameters for 2D Plotfiles
-----------------------------------

+-----------------------------+------------------+-----------------------+--------------+
| Parameter                   | Definition       | Acceptable            | Default      |
|                             |                  | Values                |              |
+=============================+==================+=======================+==============+
| **erf.plot2d_file_1**       | prefix for       | String                | “*plt2d_1_*” |
|                             | 2d plotfiles     |                       |              |
|                             | at first freq.   |                       |              |
+-----------------------------+------------------+-----------------------+--------------+
| **erf.plot2d_file_2**       | prefix for       | String                | “*plt2d_2_*” |
|                             | 2d plotfiles     |                       |              |
|                             | at seoncd freq.  |                       |              |
+-----------------------------+------------------+-----------------------+--------------+
| **erf.plot2d_int_1**        | how often (by    | Integer               | -1           |
|                             | level-0 time     | :math:`> 0`           |              |
|                             | steps) to write  |                       |              |
|                             | 2d plot files    |                       |              |
|                             | at first freq.   |                       |              |
+-----------------------------+------------------+-----------------------+--------------+
| **erf.plot2d_int_2**        | how often (by    | Integer               | -1           |
|                             | level-0 time     | :math:`> 0`           |              |
|                             | steps) to write  |                       |              |
|                             | 2d plot files    |                       |              |
|                             | at second freq.  |                       |              |
+-----------------------------+------------------+-----------------------+--------------+
| **erf.plot2d_per_1**        | how often (in    | Real                  | -1.0         |
|                             | simulation time) | :math:`> 0`           |              |
|                             | to write         |                       |              |
|                             | 2d plot files    |                       |              |
|                             | at first freq.   |                       |              |
+-----------------------------+------------------+-----------------------+--------------+
| **erf.plot2d_per_2**        | how often (in    | Real                  | -1.0         |
|                             | simulation time) | :math:`> 0`           |              |
|                             | to write         |                       |              |
|                             | 2d plot files    |                       |              |
|                             | at second freq.  |                       |              |
+-----------------------------+------------------+-----------------------+--------------+
| **erf.plot2d_vars_1**       | name of          | list of names         | None         |
|                             | variables to     |                       |              |
|                             | include in       |                       |              |
|                             | 2d plotfiles     |                       |              |
|                             | at first freq.   |                       |              |
+-----------------------------+------------------+-----------------------+--------------+
| **erf.plot2d_vars_2**       | name of          | list of names         | None         |
|                             | variables to     |                       |              |
|                             | include in       |                       |              |
|                             | p2d lotfiles     |                       |              |
|                             | at second freq.  |                       |              |
+-----------------------------+------------------+-----------------------+--------------+
| **erf.file_name_digits**    | Number of digits | Integer               | 5            |
|                             | to be appended   | :math:`> 0`           |              |
|                             | to the plot file |                       |              |
|                             | names            |                       |              |
+-----------------------------+------------------+-----------------------+--------------+

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

Output Options for 2D plotfiles
-------------------------------

+-------------------+----------------------------+
| Parameter         | Definition                 |
|                   |                            |
+===================+============================+
| **z_surf**        | Surface elevation [m]      |
+-------------------+----------------------------+
| **landmask**      | Land-sea mask              |
|                   | (land=1, sea=0)            |
|                   | [1]                        |
+-------------------+----------------------------+
| **mapfac**        | Map factors [1]            |
+-------------------+----------------------------+
| **lat_m**         | Latitude (at unstaggered   |
|                   | "mass" points) [deg]       |
+-------------------+----------------------------+
| **u_star**        | Friction velocity          |
|                   | (with SurfaceLayer only)   |
|                   | [m/s]                      |
+-------------------+----------------------------+
| **w_star**        | Convective velocity scale  |
|                   | (with SurfaceLayer only)   |
|                   | [m/s]                      |
+-------------------+----------------------------+
| **t_star**        | Temperature scale          |
|                   | (with SurfaceLayer only)   |
|                   | [K]                        |
+-------------------+----------------------------+
| **q_star**        | Humidity scale             |
|                   | (with SurfaceLayer only)   |
|                   | [kg/kg]                    |
+-------------------+----------------------------+
| **Olen**          | Obukhov length             |
|                   | (with SurfaceLayer only)   |
|                   | [m]                        |
+-------------------+----------------------------+
| **pblh**          | Diagnosed PBL height       |
|                   | (with SurfaceLayer only)   |
|                   | [m]                        |
+-------------------+----------------------------+
| **t_surf**        | Surface temperature        |
|                   | (with SurfaceLayer only)   |
|                   | [K]                        |
+-------------------+----------------------------+
| **q_surf**        | Surface humidity           |
|                   | (with SurfaceLayer only)   |
|                   | [kg/kg]                    |
+-------------------+----------------------------+
| **z0**            | Roughness height           |
|                   | (with SurfaceLayer only)   |
|                   | [m]                        |
+-------------------+----------------------------+
| **OLR**           | Outgoing long wavelength   |
|                   | radiation (with RRTMGP)    |
|                   | [W/m^2]                    |
+-------------------+----------------------------+
| **sens_flux**     | Sensible heat flux         |
|                   | (with SurfaceLayer only)   |
|                   | [kg m^-2 s^-1 K]           |
+-------------------+----------------------------+
| **laten_flux**    | Latent heat flux           |
|                   | (with SurfaceLayer only)   |
|                   | [kg m^-2 s^-1]             |
+-------------------+----------------------------+
| **surf_pres**     | Surface pressure           |
|                   | [Pa]                       |
+-------------------+----------------------------+
| **integrated_qv** | Column-integrated water    |
|                   | vapor; zero if moisture is |
|                   | disabled [kg/m^2]          |
+-------------------+----------------------------+

Examples of Usage
-----------------

In an input file, the user can select parameters to plot by supplying a space-delimited
list to **erf.plot_vars_1** or **erf.plot_vars_2**.

-  **erf.plot_vars_1** = *option1* *option2* *option3*
