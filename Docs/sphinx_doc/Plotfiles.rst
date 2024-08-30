.. role:: cpp(code)
  :language: c++

.. _sec:Plotfiles:

*********
Plotfiles
*********
.. toctree::
   :maxdepth: 1

Controlling PlotFile Generation
===============================

"Plotfiles" can be written very efficiently in parallel in a native AMReX format
or in HDF5.  They can also be written in NetCDF.

The following options in the inputs file control the generation of plotfiles.
Note that plotfiles can be written at two different frequencies; the names,
frequency and content of the two streams are controlled separately.

.. _list-of-parameters-9:

List of Parameters
------------------

+-----------------------------+------------------+-----------------------+------------+
| Parameter                   | Definition       | Acceptable            | Default    |
|                             |                  | Values                |            |
+=============================+==================+=======================+============+
| **erf.plotfile_type**       | AMReX, NETCDF    | "amrex" or            | "amrex"    |
|                             | or HDF5          | "netcdf / "NetCDF" or |            |
|                             |                  | "hdf5" / "HDF5"       |            |
+-----------------------------+------------------+-----------------------+------------+
| **erf.plot_file_1**         | prefix for       | String                | “*plt_1_*” |
|                             | plotfiles        |                       |            |
|                             | at first freq.   |                       |            |
+-----------------------------+------------------+-----------------------+------------+
| **erf.plot_file_2**         | prefix for       | String                | “*plt_2_*” |
|                             | plotfiles        |                       |            |
|                             | at seoncd freq.  |                       |            |
+-----------------------------+------------------+-----------------------+------------+
| **erf.plot_int_1**          | how often (by    | Integer               | -1         |
|                             | level-0 time     | :math:`> 0`           |            |
|                             | steps) to write  |                       |            |
|                             | plot files       |                       |            |
|                             | at first freq.   |                       |            |
+-----------------------------+------------------+-----------------------+------------+
| **erf.plot_int_2**          | how often (by    | Integer               | -1         |
|                             | level-0 time     | :math:`> 0`           |            |
|                             | steps) to write  |                       |            |
|                             | plot files       |                       |            |
|                             | at second freq.  |                       |            |
+-----------------------------+------------------+-----------------------+------------+
| **erf.plot_per_1**          | how often in     | Real                  | -1.0       |
|                             | simulation time  | :math:`> 0`           |            |
|                             | to write         |                       |            |
|                             | plot files       |                       |            |
|                             | at first freq.   |                       |            |
+-----------------------------+------------------+-----------------------+------------+
| **erf.plot_per_2**          | how often in     | Real                  | -1.0       |
|                             | simulation time  | :math:`> 0`           |            |
|                             | to write         |                       |            |
|                             | plot files       |                       |            |
|                             | at second freq.  |                       |            |
+-----------------------------+------------------+-----------------------+------------+
| **erf.plot_vars_1**         | name of          | list of names         | None       |
|                             | variables to     |                       |            |
|                             | include in       |                       |            |
|                             | plotfiles        |                       |            |
|                             | at first freq.   |                       |            |
+-----------------------------+------------------+-----------------------+------------+
| **erf.plot_vars_2**         | name of          | list of names         | None       |
|                             | variables to     |                       |            |
|                             | include in       |                       |            |
|                             | plotfiles        |                       |            |
|                             | at seoncd freq.  |                       |            |
+-----------------------------+------------------+-----------------------+------------+

.. _notes-5:

Notes
-----

-  The NeTCDF option is only available if ERF has been built with USE_NETCDF enabled.

.. _examples-of-usage-8:

Examples of Usage
-----------------

-  **erf.plotfile_type** = *amrex*

-  **erf.plot_file_1** = *plt_run*

-  **erf.plot_int_1** = 10

   means that native plot files (actually directories) starting with the prefix
   “*plt_run*” will be generated every 10 level-0 time steps. If using
   amrex format, that directory names will be *plt_run00000*, *plt_run00010*,
   *plt_run00020*, etc.  If using HDF5 format, the names will have ".h5"
   appended;  if using NetCDF format, the names will have ".nc" appended.

   In addition, while the amrex plotfiles will contain data at all of the refinement
   levels,  NetCDF files are separated by level.

PlotFile Outputs
================

Plotfiles can include the quantities of several simulation parameters as output.
They are summarized in the list below. Note that temporally averaged quantities
(e.g., ``u_t_avg, v_t_avg, w_t_avg, umag_t_avg``) require the user to enable the
storage of the time averaged variables with ``erf.time_avg_vel = true``.

Output Options
--------------

+-----------------------------+------------------+
| Parameter                   | Definition       |
|                             |                  |
+=============================+==================+
| **x_velocity**              | Velocity in x    |
|                             | direction        |
|                             |                  |
+-----------------------------+------------------+
| **y_velocity**              | Velocity in y    |
|                             | direction        |
|                             |                  |
+-----------------------------+------------------+
| **z_velocity**              | Velocity in z    |
|                             | direction        |
|                             |                  |
+-----------------------------+------------------+
| **density**                 | Total density    |
|                             |                  |
|                             |                  |
+-----------------------------+------------------+
| **dens_hse**                | Hydrostatic      |
|                             | density          |
|                             |                  |
+-----------------------------+------------------+
| **pert_dens**               | Perturbational   |
|                             | density          |
|                             |                  |
+-----------------------------+------------------+
| **pressure**                | Total pressure   |
|                             |                  |
|                             |                  |
+-----------------------------+------------------+
| **pres_hse**                | Hydrostatic      |
|                             | pressure         |
|                             |                  |
+-----------------------------+------------------+
| **pert_pres**               | Perturbational   |
|                             | pressure         |
|                             |                  |
+-----------------------------+------------------+
| **pres_hse_x**              | Derivative of    |
|                             | hydrostatic      |
|                             | pressure in x    |
+-----------------------------+------------------+
| **pres_hse_y**              | Derivative of    |
|                             | hydrostatic      |
|                             | pressure in y    |
+-----------------------------+------------------+
| **dpdx**                    | Pressure gradient|
|                             | in x direction   |
|                             |                  |
+-----------------------------+------------------+
| **dpdy**                    | Pressure gradient|
|                             | in y direction   |
|                             |                  |
+-----------------------------+------------------+
| **temp**                    | Temperature      |
|                             |                  |
|                             |                  |
+-----------------------------+------------------+
| **theta**                   | Potential        |
|                             | temperature      |
|                             |                  |
+-----------------------------+------------------+
| **rhotheta**                | Density * theta  |
|                             |                  |
|                             |                  |
+-----------------------------+------------------+
| **KE**                      | Kinetic energy   |
|                             |                  |
|                             |                  |
+-----------------------------+------------------+
| **QKE**                     | Turbulent        |
|                             | kinetic energy   |
|                             | * 2              |
+-----------------------------+------------------+
| **rhoKE**                   | Density * KE     |
|                             |                  |
|                             |                  |
+-----------------------------+------------------+
| **rhoQKE**                  | Density * QKE    |
|                             |                  |
|                             |                  |
+-----------------------------+------------------+
| **scalar**                  | Scalar magnitude |
|                             |                  |
|                             |                  |
+-----------------------------+------------------+
| **vorticity_x**             | x-component of   |
|                             | vorticity        |
|                             |                  |
+-----------------------------+------------------+
| **vorticity_y**             | y-component of   |
|                             | vorticity        |
|                             |                  |
+-----------------------------+------------------+
| **vorticity_z**             | z-component of   |
|                             | vorticity        |
|                             |                  |
+-----------------------------+------------------+
| **u_t_avg**                 | time average of  |
|                             | x-component of   |
|                             | velocity         |
+-----------------------------+------------------+
| **v_t_avg**                 | time average of  |
|                             | y-component of   |
|                             | velocity         |
+-----------------------------+------------------+
| **w_t_avg**                 | time average of  |
|                             | z-component of   |
|                             | velocity         |
+-----------------------------+------------------+
| **umag_t_avg**              | time average of  |
|                             | velocity mag     |
|                             |                  |
+-----------------------------+------------------+
| **rhoadv_0**                | Conserved scalar |
|                             |                  |
|                             |                  |
+-----------------------------+------------------+
| **soundspeed**              | Sound speed      |
|                             |                  |
|                             |                  |
+-----------------------------+------------------+
| **z_phys**                  | Terrain height   |
|                             |                  |
|                             |                  |
+-----------------------------+------------------+
| **detJ**                    | Jacobian         |
|                             | determinant      |
|                             |                  |
+-----------------------------+------------------+
| **mapfac**                  | Map scale factor |
|                             |                  |
|                             |                  |
+-----------------------------+------------------+
| **lat_m**                   | Latitude at mass |
|                             | points           |
|                             |                  |
+-----------------------------+------------------+
| **lon_m**                   | Longitude at     |
|                             | mass points      |
|                             |                  |
+-----------------------------+------------------+
| **Kmv**                     | Vertical         |
|                             | Eddy Diffusivity |
|                             | of Momentum      |
+-----------------------------+------------------+
| **Kmh**                     | Horizontal       |
|                             | Eddy Diffusivity |
|                             | of Momentum      |
+-----------------------------+------------------+
| **Khv**                     | Vertical         |
|                             | Eddy Diffusivity |
|                             | of Heat          |
+-----------------------------+------------------+
| **Khh**                     | Horizontal       |
|                             | Eddy Diffusivity |
|                             | of Heat          |
+-----------------------------+------------------+
| **qt**                      | Nonprecipitating |
|                             | water (qv + qc + |
|                             | qi)              |
+-----------------------------+------------------+
| **qp**                      | Precipitating    |
|                             | water (rain +    |
|                             | snow + graupel)  |
+-----------------------------+------------------+
| **qc**                      | Cloud water      |
|                             | mixing ratio     |
|                             |                  |
+-----------------------------+------------------+
| **qi**                      | Cloud ice        |
|                             | mixing ratio     |
|                             |                  |
+-----------------------------+------------------+
| **qv**                      | Water vapor      |
|                             | mixing ratio     |
|                             |                  |
+-----------------------------+------------------+
| **rhoQt**                   | Density * qt     |
|                             |                  |
|                             |                  |
+-----------------------------+------------------+
| **rhoQp**                   | Density * qp     |
|                             |                  |
|                             |                  |
+-----------------------------+------------------+

Examples of Usage
-----------------

In an input file, the user can select parameters to plot by supplying a space-delimited
list to **erf.plot_vars_1** or **erf.plot_vars_2**.

-  **erf.plot_vars_1** = *option1* *option2* *option3*

