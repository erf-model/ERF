.. role:: cpp(code)
  :language: c++

.. _sec:Plotfiles:

******
Plotfiles
******
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
|                             | at seoncd freq.  |                       |            |
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
