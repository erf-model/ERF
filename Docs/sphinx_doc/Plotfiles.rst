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

"Plotfiles" (which are really directories) can be written very efficiently
in parallel in a native AMReX format or in HDF5.  They can also be written
in NetCDF.

The following options in the inputs file control the generation of plotfiles.

.. _list-of-parameters-9:

List of Parameters
------------------

+-----------------------------+------------------+-----------------------+---------+
| Parameter                   | Definition       | Acceptable            | Default |
|                             |                  | Values                |         |
+=============================+==================+=======================+=========+
| **erf.plotfile_type**       | AMReX, NETCDF    | "amrex" or            | "amrex" |
|                             | or HDF5          | "netcdf / "NetCDF" or |         |
|                             |                  | "hdf5" / "HDF5"       |         |
+-----------------------------+------------------+-----------------------+---------+
| **erf.plot_file**           | prefix for       | String                | “*plt*” |
|                             | plotfiles        |                       |         |
+-----------------------------+------------------+-----------------------+---------+
| **erf.plot_int**            | how often (by    | Integer               | -1      |
|                             | level-0 time     | :math:`> 0`           |         |
|                             | steps) to write  |                       |         |
|                             | plot files       |                       |         |
+-----------------------------+------------------+-----------------------+---------+
| **erf.plot_vars**           | name of          | list of names         | None    |
|                             | variables to     |                       |         |
|                             | include in       |                       |         |
|                             | plotfiles        |                       |         |
+-----------------------------+------------------+-----------------------+---------+

.. _notes-5:

Notes
-----

-  The NeTCDF option is only available if ERF has been built with USE_NETCDF enabled.

.. _examples-of-usage-8:

Examples of Usage
-----------------

-  **erf.plotfile_type** = *amrex*

-  **erf.plot_file** = *plt_run*

-  **erf.plot_int** = 10

   means that native plot files (actually directories) starting with the prefix
   “*plt_run*” will be generated every 10 level-0 time steps. If using
   amrex format, that directory names will be *plt_run00000*, *plt_run00010*,
   *plt_run00020*, etc.  If using HDF5 format, the names will have ".h5"
   appended;  if using NetCDF format, the names will have ".nc" appended.

   In addition, while the amrex plotfiles will contain data at all of the refinement
   levels,  NetCDF files are separated by level.
