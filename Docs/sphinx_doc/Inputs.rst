.. role:: cpp(code)
  :language: c++

.. _sec:Inputs:

******
Inputs
******
.. toctree::
   :maxdepth: 1

The ERF executable reads run-time information from an “inputs” file which you put on the command line.
This section describes the inputs which can be specified either in the inputs file or on the command line.
If a value is specified on the command line, that value will override a value specified in the inputs file.

Problem Geometry
================

List of Parameters
------------------

+--------------------------+-----------------+-----------------+-------------+
| Parameter                | Definition      | Acceptable      | Default     |
|                          |                 | Values          |             |
+==========================+=================+=================+=============+
| **geometry.prob_lo**     | physical        | Real            | must be set |
|                          | location of low |                 |             |
|                          | corner of the   |                 |             |
|                          | domain          |                 |             |
+--------------------------+-----------------+-----------------+-------------+
| **geometry.prob_hi**     | physical        | Real            | must be set |
|                          | location of     |                 |             |
|                          | high corner of  |                 |             |
|                          | the domain      |                 |             |
+--------------------------+-----------------+-----------------+-------------+
| **geometry.is_periodic** | is the domain   | 0 if false, 1   | 0 0 0       |
|                          | periodic in     | if true         |             |
|                          | this direction  |                 |             |
+--------------------------+-----------------+-----------------+-------------+

Examples of Usage
-----------------

-  **geometry.prob_lo** = 0 0 0
   defines the low corner of the domain at (0,0,0) in physical space.

-  **geometry.prob_hi** = 1.e8 2.e8 2.e8
   defines the high corner of the domain at (1.e8,2.e8,2.e8) in
   physical space.

-  **geometry.is_periodic** = 0 1 0
   says the domain is periodic in the y-direction only.

Domain Boundary Conditions
==========================

.. _list-of-parameters-1:

List of Parameters
------------------

+---------------+---------------------------------+-------------------+----------------------------+
| Parameter     | Definition                      | Acceptable Values | Default                    |
+===============+=================================+===================+============================+
| **xlo.type**  | boundary type of xlo face       |                   | must be set if not periodic|
+---------------+---------------------------------+-------------------+----------------------------+
| **xhi.type**  | boundary type of xhi face       |                   | must be set if not periodic|
+---------------+---------------------------------+-------------------+----------------------------+
| **ylo.type**  | boundary type of ylo face       |                   | must be set if not periodic|
+---------------+---------------------------------+-------------------+----------------------------+
| **yhi.type**  | boundary type of yhi face       |                   | must be set if not periodic|
+---------------+---------------------------------+-------------------+----------------------------+
| **zlo.type**  | boundary type of zlo face       |                   | must be set if not periodic|
+---------------+---------------------------------+-------------------+----------------------------+
| **zhi.type**  | boundary type of zhi face       |                   | must be set if not periodic|
+---------------+---------------------------------+-------------------+----------------------------+


Resolution
==========

.. _list-of-parameters-2:

List of Parameters
------------------

+---------------------------+-----------------+-----------------+-------------+
| Parameter                 | Definition      | Acceptable      | Default     |
|                           |                 | Values          |             |
+===========================+=================+=================+=============+
| **amr.n_cell**            | number of cells | Integer > 0     | must be set |
|                           | in each         |                 |             |
|                           | direction at    |                 |             |
|                           | the coarsest    |                 |             |
|                           | level           |                 |             |
+---------------------------+-----------------+-----------------+-------------+
| **amr.max_level**         | number of       | Integer >= 0    | must be set |
|                           | levels of       |                 |             |
|                           | refinement      |                 |             |
|                           | above the       |                 |             |
|                           | coarsest level  |                 |             |
+---------------------------+-----------------+-----------------+-------------+
| **amr.ref_ratio**         | ratio of coarse | 2 / 3 / 4       | 2 for all   |
|                           | to fine grid    | (one per level) | levels      |
|                           | spacing between |                 |             |
|                           | subsequent      |                 |             |
|                           | levels          |                 |             |
+---------------------------+-----------------+-----------------+-------------+
| **amr.ref_ratio_vect**    | ratio of coarse | 3 integers      | 2 for all   |
|                           | to fine grid    | (one per dir)   | directions  |
|                           | spacing between | 2 / 3 / 4       |             |
|                           | subsequent      |                 |             |
|                           | levels          |                 |             |
+---------------------------+-----------------+-----------------+-------------+
| **amr.regrid_int**        | how often to    | Integer > 0     | must be set |
|                           | regrid          |                 |             |
+---------------------------+-----------------+-----------------+-------------+
| **amr.regrid_on_restart** | should we       | 0 or 1          | 0           |
|                           | regrid          |                 |             |
|                           | immediately     |                 |             |
|                           | after           |                 |             |
|                           | restarting      |                 |             |
+---------------------------+-----------------+-----------------+-------------+

Note: if **amr.max_level** = 0 then you do not need to set
**amr.ref_ratio** or **amr.regrid_int**.

.. _examples-of-usage-2:

Examples of Usage
-----------------

-  **amr.n_cell** = 32 64 64

   would define the domain to have 32 cells in the x-direction, 64 cells
   in the y-direction, and 64 cells in the z-direction *at the coarsest level*.

-  | **amr.max_level** = 2
   | would allow a maximum of 2 refined levels in addition to the coarse
     level. Note that these additional levels will only be created only
     if the tagging criteria are such that cells are flagged as needing
     refinement. The number of refined levels in a calculation must be
     :math:`\leq` **amr.max_level**, but can change in time and need not
     always be equal to **amr.max_level**.

-  | **amr.ref_ratio** = 2 3
   | would set factor 2 refinement between levels 0 and 1, and factor 3
     refinement between levels 1 and 2. Note that you must have at least
     **amr.max_level** values of **amr.ref_ratio** (Additional values
     may appear in that line and they will be ignored).

-  | **amr.ref_ratio_vect** = 2 4 3
   | would set factor {2 in x-dir, 4 in y-dir, 3 in z-dir} refinement between
     all adjacent levels.    Note that you must specify 3 values, one for
     each coordinate direction.

-  | **amr.regrid_int** = 2 2
   | tells the code to regrid every 2 steps. Thus in this example, new
     level-1 grids will be created every 2 level-0 time steps, and new
     level-2 grids will be created every 2 level-1 time steps.

Regridding
==========

Overview
--------

The user defines how to tag individual cells at a given level for refinement.
This list of tagged cells is sent to a grid generation routine, which uses the
Berger–Rigoutsos algorithm to create rectangular grids that contain the
tagged cells.

See :ref:`MeshRefinement` for more details on how to specify regions for
refinement.

.. _list-of-parameters-4:

List of Parameters
------------------

+----------------------------+----------------+----------------+----------------+
| Parameter                  | Definition     | Acceptable     | Default        |
|                            |                | Values         |                |
+============================+================+================+================+
| **amr.regrid_file**        | name of file   | text           | no file        |
|                            | from which to  |                |                |
|                            | read the grids |                |                |
+----------------------------+----------------+----------------+----------------+
| **amr.grid_eff**           | grid           | Real > 0, < 1  | 0.7            |
|                            | efficiency at  |                |                |
|                            | coarse level   |                |                |
|                            | at which grids |                |                |
|                            | are created    |                |                |
+----------------------------+----------------+----------------+----------------+
| **amr.n_error_buf**        | radius of      | Integer >= 0   | 1              |
|                            | additional     |                |                |
|                            | tagging around |                |                |
|                            | already tagged |                |                |
|                            | cells          |                |                |
+----------------------------+----------------+----------------+----------------+
| **amr.max_grid_size**      | maximum size   | Integer > 0    | 32             |
|                            | of a grid in   |                |                |
|                            | any direction  |                |                |
+----------------------------+----------------+----------------+----------------+
| **amr.max_grid_size**      | maximum size   | Integer        | 32             |
+----------------------------+----------------+----------------+----------------+
| **amr.blocking_factor**    | grid size must | Integer > 0    | 2              |
|                            | be a multiple  |                |                |
|                            | of this        |                |                |
+----------------------------+----------------+----------------+----------------+
| **amr.refine_grid_layout** | refine grids   | 0 if false, 1  | 1              |
|                            | more if # of   | if true        |                |
|                            | processors     |                |                |
|                            | :math:`>` # of |                |                |
|                            | grids          |                |                |
+----------------------------+----------------+----------------+----------------+

.. _notes-2:

Notes
-----

-  **amr.n_error_buf**, **amr.max_grid_size** and
   **amr.blocking_factor** can be read in as a single value which is
   assigned to every level, or as multiple values, one for each level

-  **amr.max_grid_size** at every level must be even

-  **amr.blocking_factor** at every level must be a power of 2

-  the domain size **amr.n_cell** must be a multiple of
   **amr.blocking_factor** at level 0

-  **amr.max_grid_size** must be a multiple of **amr.blocking_factor**
   at every level

.. _examples-of-usage-3:

Examples of Usage
-----------------

-  | **amr.regrid_file** = *fixed_grids*
   | In this case the list of grids at each fine level are contained in
     the file *fixed_grids*, which will be read during the gridding
     procedure. These grids must not violate the **amr.max_grid_size**
     criterion. The rest of the gridding procedure described below will
     not occur if **amr.regrid_file** is set.

-  | **amr.grid_eff** = 0.9
   | During the grid creation process, at least 90% of the cells in each
     grid at the level at which the grid creation occurs must be tagged
     cells. Note that this is applied at the coarsened level at which
     the grids are actually made, and before **amr.max_grid_size** is
     imposed.

-  | **amr.max_grid_size** = 64
   | The final grids will be no longer than 64 cells on a side at every
     level.

-  | **amr.max_grid_size** = 64 32 16
   | The final grids will be no longer than 64 cells on a side at level
     0, 32 cells on a side at level 1, and 16 cells on a side at level
     2.

-  | **amr.blocking_factor** = 32
   | The dimensions of all the final grids will be multiples of 32 at
     all levels.

-  | **amr.blocking_factor** = 32 16 8
   | The dimensions of all the final grids will be multiples of 32 at
     level 0, multiples of 16 at level 1, and multiples of 8 at level 2.

.. _subsec:grid-generation:

Gridding and Load Balancing
---------------------------

The ERF gridding and load balancing strategy is based on that in AMReX.
See the `Gridding`_ section of the AMReX documentation for details.

.. _`Gridding`: https://amrex-codes.github.io/amrex/docs_html/ManagingGridHierarchy_Chapter.html

Simulation Time
===============

.. _list-of-parameters-5:

List of Parameters
------------------

+-----------------+---------------------------+--------------+---------+
| Parameter       | Definition                | Acceptable   | Default |
|                 |                           | Values       |         |
+=================+===========================+==============+=========+
| **max_step**    | maximum number of level 0 | Integer >= 0 | -1      |
|                 | time steps                |              |         |
+-----------------+---------------------------+--------------+---------+
| **stop_time**   | final simulation          | Real >= 0    | -1.0    |
|                 | time                      |              |         |
+-----------------+---------------------------+--------------+---------+

.. _notes-3:

Notes
-----

To control the number of time steps, you can limit by the maximum number
of level-0 time steps (**max_step**), or the final simulation time
(**stop_time**), or both. The code will stop at whichever criterion
comes first. Note that if the code reaches **stop_time** then the final
time step will be shortened so as to end exactly at **stop_time**, not
pass it.

.. _examples-of-usage-4:

Examples of Usage
-----------------

-  **max_step** = 1000

-  **stop_time** = 1.0

will end the calculation when either the simulation time reaches 1.0 or
the number of level-0 steps taken equals 1000, whichever comes first.

Time Step
=========

.. _list-of-parameters-6:

List of Parameters
------------------

+---------------------+----------------+----------------+----------------+
| Parameter           | Definition     | Acceptable     | Default        |
|                     |                | Values         |                |
+=====================+================+================+================+
| **erf.cfl**         | CFL number for | Real > 0 and   | 0.8            |
|                     | hydro          | <= 1           |                |
|                     |                |                |                |
|                     |                |                |                |
+---------------------+----------------+----------------+----------------+
| **erf.init_shrink** | factor by      | Real > 0 and   | 1.0            |
|                     | which to       | <= 1           |                |
|                     | shrink the     |                |                |
|                     | initial time   |                |                |
|                     | step           |                |                |
+---------------------+----------------+----------------+----------------+
| **erf.change_max**  | factor by      | Real >= 1      | 1.1            |
|                     | which the time |                |                |
|                     | step can grow  |                |                |
|                     | in subsequent  |                |                |
|                     | steps          |                |                |
+---------------------+----------------+----------------+----------------+
| **erf.fixed_dt**    | level-0 time   | Real > 0       | unused if not  |
|                     | step           |                | set            |
|                     | regardless of  |                |                |
|                     | cfl or other   |                |                |
|                     | settings       |                |                |
+---------------------+----------------+----------------+----------------+
| **erf.initial_dt**  | initial        | Real > 0       | unused if not  |
|                     | level-0 time   |                | set            |
|                     | step           |                |                |
|                     | regardless of  |                |                |
|                     | other settings |                |                |
+---------------------+----------------+----------------+----------------+
| **erf.max_dt**      | time step      | Real > 0       | 1.e20          |
|                     | above which    |                |                |
|                     | calculation    |                |                |
|                     | will abort     |                |                |
+---------------------+----------------+----------------+----------------+
| **erf.dt_cutoff**   | time step      | Real > 0       | 0.0            |
|                     | below which    |                |                |
|                     | calculation    |                |                |
|                     | will abort     |                |                |
+---------------------+----------------+----------------+----------------+

.. _examples-of-usage-5:

Examples of Usage
-----------------

-  | **erf.cfl** = 0.9
   | defines the timestep as dt = cfl \* dx / (u+c)

-  | **erf.init_shrink** = 0.01
   | sets the initial time step to 1% of what it would be otherwise.

-  | **erf.change_max** = 1.1
   | allows the time step to increase by no more than 10% in this case.
     Note that the time step can shrink by any factor; this only
     controls the extent to which it can grow.

-  | **erf.fixed_dt** = 1.e-4
   | sets the level-0 time step to be 1.e-4 for the entire simulation,
     ignoring the other timestep controls. Note that if
     **erf.init_shrink** :math:`\neq 1` then the first time step will in
     fact be **erf.init_shrink** \* **erf.fixed_dt**.

-  | **erf.initial_dt** = 1.e-4
   | sets the *initial* level-0 time step to be 1.e-4 regardless of
     **erf.cfl** or **erf.fixed_dt**. The time step can grow in
     subsequent steps by a factor of **erf.change_max** each step.

-  | **erf.dt_cutoff** = 1.e-20
   | tells the code to abort if the time step ever gets below 1.e-20.
     This is a safety mechanism so that if things go nuts you don’t burn
     through your entire computer allocation because you don’t realize
     the code is misbehaving.

Restart Capability
==================

|  has a standard sort of checkpointing and restarting capability. In
  the inputs file, the following options control the generation of
  checkpoint files (which are really directories):

.. _list-of-parameters-8:

List of Parameters
------------------

+---------------------------------+----------------+----------------+----------------+
| Parameter                       | Definition     | Acceptable     | Default        |
|                                 |                | Values         |                |
+=================================+================+================+================+
| **amr.check_file**              | prefix for     | String         | “*chk*”        |
|                                 | restart files  |                |                |
+---------------------------------+----------------+----------------+----------------+
| **amr.check_int**               | how often (by  | Integer        | -1             |
|                                 | level-0 time   | :math:`> 0`    |                |
|                                 | steps) to      |                |                |
|                                 | write restart  |                |                |
|                                 | files          |                |                |
+---------------------------------+----------------+----------------+----------------+
| **amr.check_per**               | how often (by  | Real           | -1.0           |
|                                 | simulation     | :math:`> 0`    |                |
|                                 | time) to write |                |                |
|                                 | restart files  |                |                |
+---------------------------------+----------------+----------------+----------------+
| **amr.restart**                 | name of the    | String         | not used if    |
|                                 | file           |                | not set        |
|                                 | (directory)    |                |                |
|                                 | from which to  |                |                |
|                                 | restart        |                |                |
+---------------------------------+----------------+----------------+----------------+
| **amr.checkpoint_files_output** | should we      | 0 or 1         | 1              |
|                                 | write          |                |                |
|                                 | checkpoint     |                |                |
|                                 | files          |                |                |
+---------------------------------+----------------+----------------+----------------+
| **amr.check_nfiles**            | how parallel   | Integer        | 64             |
|                                 | is the writing | :math:`\geq 1` |                |
|                                 | of the         |                |                |
|                                 | checkpoint     |                |                |
|                                 | files          |                |                |
+---------------------------------+----------------+----------------+----------------+
| **amr.checkpoint_on_restart**   | should we      | 0 or 1         | 0              |
|                                 | write a        |                |                |
|                                 | checkpoint     |                |                |
|                                 | immediately    |                |                |
|                                 | after          |                |                |
|                                 | restarting     |                |                |
+---------------------------------+----------------+----------------+----------------+
| **erf.dump_old**                | do we store    | true / false   | false          |
|                                 | old data in    |                |                |
|                                 | checkpoints?   |                |                |
+---------------------------------+----------------+----------------+----------------+
| **erf.checkpoint_type**         | AMReX chk file | "amrex" or     | "amrex"        |
|                                 | type or NetCDF | "NetCDF"       |                |
|                                 |                |                |                |
+---------------------------------+----------------+----------------+----------------+

.. _notes-4:

Notes
-----

-  You should specify either **amr.check_int** or **amr.check_per**. Do
   not try to specify both.

-  Note that if **amr.check_per** is used then in order to hit that
   exact time the code may modify the time step slightly, which will
   change your results ever so slightly than if you didn’t set this
   flag.

-  Note that **amr.plotfile_on_restart** and
   **amr.checkpoint_on_restart** only take effect if
   **amr.regrid_on_restart** is in effect.

-  See the Software Section for more details on parallel I/O and the
   **amr.check_nfiles** parameter.

-  If you are doing a scaling study then set
   **amr.checkpoint_files_output** = 0 so you can test scaling of the
   algorithm without I/O.

-  If you compile with NetCDF capability, you may choose to dump a
   NetCDF format checkpoint file instead of the default.

.. _examples-of-usage-7:

Examples of Usage
-----------------

-  **amr.check_file** = *chk_run*

-  **amr.check_int** = 10

   means that restart files (really directories) starting with the
   prefix “*chk_run*” will be generated every 10 level-0 time steps. The
   directory names will be *chk_run00000*, *chk_run00010*,
   *chk_run00020*, etc.

If instead you specify

-  **amr.check_file** = *chk_run*

-  **amr.check_per** = 0.5

   then restart files (really directories) starting with the prefix
   “*chk_run*” will be generated every 0.1 units of simulation time. The
   directory names will be *chk_run00000*, *chk_run00043*,
   *chk_run00061*, etc, where :math:`t = 0.1` after 43 level-0 steps,
   :math:`t = 0.2` after 61 level-0 steps, etc.

To restart from *chk_run00061*,for example, then set

-  **amr.restart** = *chk_run00061*

.. _sec:PlotFiles:

Controlling PlotFile Generation
===============================

The main output from  is in the form of plotfiles (which are really
directories) written in a native AMReX format, which can be read by
Visit, ParaView and yt. The following options in the inputs file control the
generation of plotfiles

.. _list-of-parameters-9:

List of Parameters
------------------

+-----------------------------+------------------+------------------+---------+
| Parameter                   | Definition       | Acceptable       | Default |
|                             |                  | Values           |         |
+=============================+==================+==================+=========+
| **amr.plot_file**           | prefix for       | String           | “*plt*” |
|                             | plotfiles        |                  |         |
+-----------------------------+------------------+------------------+---------+
| **amr.plot_int**            | how often (by    | Integer          | -1      |
|                             | level-0 time     | :math:`> 0`      |         |
|                             | steps) to write  |                  |         |
|                             | plot files       |                  |         |
+-----------------------------+------------------+------------------+---------+
| **amr.plot_per**            | how often (by    | Real :math:`> 0` | -1.0    |
|                             | simulation time) |                  |         |
|                             | to write plot    |                  |         |
|                             | files            |                  |         |
+-----------------------------+------------------+------------------+---------+
| **amr.plot_vars**           | name of state    | ALL, NONE or     | ALL     |
|                             | variables to     | list             |         |
|                             | include in       |                  |         |
|                             | plotfiles        |                  |         |
+-----------------------------+------------------+------------------+---------+
| **amr.derive_plot_vars**    | name of derived  | ALL, NONE or     | NONE    |
|                             | variables to     | list             |         |
|                             | include in       |                  |         |
|                             | plotfiles        |                  |         |
+-----------------------------+------------------+------------------+---------+
| **amr.plot_files_output**   | should we write  | 0 or 1           | 1       |
|                             | plot files       |                  |         |
+-----------------------------+------------------+------------------+---------+
| **amr.plotfile_on_restart** | should we write  | 0 or 1           | 0       |
|                             | a plotfile       |                  |         |
|                             | immediately      |                  |         |
|                             | after restarting |                  |         |
+-----------------------------+------------------+------------------+---------+
| **amr.plot_nfiles**         | how parallel is  | Integer          | 64      |
|                             | the writing of   | :math:`\geq 1`   |         |
|                             | the plotfiles    |                  |         |
+-----------------------------+------------------+------------------+---------+
| **erf.plotfile_type**       | AMReX plt file   | "amrex" or       | "amrex" |
|                             | type or NetCDF   | "NetCDF"         |         |
|                             |                  |                  |         |
+-----------------------------+------------------+------------------+---------+

.. _notes-5:

Notes
-----

-  You should specify either **amr.plot_int** or **amr.plot_per**. Do
   not try to specify both.

-  Note that if **amr.plot_per** is used then in order to hit that exact
   time the code may modify the time step slightly, which will change
   your results ever so slightly than if you didn’t set this flag.

-  See the Software Section for more details on parallel I/O and the
   **amr.plot_nfiles** parameter.

-  If you are doing a scaling study then set **amr.plot_files_output** =
   0 so you can test scaling of the algorithm without I/O.

-  By default, plotfiles are written in double precision (NATIVE
   format). If you want to save space by writing them in single
   precision, set the fab.format flag to IEEE32.

-  If you compile with NetCDF capability, you may choose to dump a
   NetCDF format plot file instead of the default.

.. _examples-of-usage-8:

Examples of Usage
-----------------

-  **amr.plot_file** = *plt_run*

-  **amr.plot_int** = 10

   means that plot files (really directories) starting with the prefix
   “*plt_run*” will be generated every 10 level-0 time steps. The
   directory names will be *plt_run00000*, *plt_run00010*,
   *plt_run00020*, etc.

If instead you specify

-  **amr.plot_file** = *plt_run*

-  **amr.plot_per** = 0.5

   then restart files (really directories) starting with the prefix
   “plt_run” will be generated every 0.1 units of simulation time. The
   directory names will be *plt_run00000*, *plt_run00043*,
   *plt_run00061*, etc, where :math:`t = 0.1` after 43 level-0 steps,
   :math:`t = 0.2` after 61 level-0 steps, etc.

Screen Output
=============

.. _list-of-parameters-10:

List of Parameters
------------------

+----------------------------+------------------+----------------+----------------+
| Parameter                  | Definition       | Acceptable     | Default        |
|                            |                  | Values         |                |
+============================+==================+================+================+
| **amr.v**                  | verbosity of     | 0 or 1         | 0              |
|                            | Amr.cpp          |                |                |
+----------------------------+------------------+----------------+----------------+
| **erf.v**                  | verbosity of     | 0 or 1         | 0              |
|                            | ERF.cpp          |                |                |
+----------------------------+------------------+----------------+----------------+
| **erf.sum_interval**       | if               |                |                |
|                            | :math:`> 0,`     |                |                |
|                            | how often (in    |                |                |
|                            | level-0 time     |                |                |
|                            | steps)           |                |                |
|                            | to compute and   | Integer        | -1             |
|                            | print integral   |                |                |
|                            | quantities       |                |                |
+----------------------------+------------------+----------------+----------------+

.. _examples-of-usage-9:

Examples of Usage
-----------------

-  | **erf.sum_interval** = 2
   | if **erf.sum_interval** :math:`> 0` then the code computes and
     prints certain integral quantities, such as total mass, momentum
     and energy in the domain every **erf.sum_interval** level-0 steps.
     In this example the code will print these quantities every two
     coarse time steps. The print statements have the form
   | TIME= 1.91717746 MASS= 1.792410279e+34
   | for example. If this line is commented out then it will not compute
     and print these quanitities.


Diffusive Physics
=================

.. _list-of-parameters-12:

List of Parameters
------------------

+----------------------------------+--------------------+-------------------+-------------+
| Parameter                        | Definition         | Acceptable        | Default     |
|                                  |                    | Values            |             |
+==================================+====================+===================+=============+
| **erf.alpha_T**                  | Diffusion coeff.   | Real              | 0.0         |
|                                  | for temperature    |                   |             |
+----------------------------------+--------------------+-------------------+-------------+
| **erf.alpha_C**                  | Diffusion coeff.   | Real              | 0.0         |
|                                  | for scalar         |                   |             |
+----------------------------------+--------------------+-------------------+-------------+
| **erf.rho0_trans**               | Reference density  | Real              | 1.0         |
|                                  | to compute const.  |                   |             |
|                                  | rho*Alpha          |                   |             |
+----------------------------------+--------------------+-------------------+-------------+
| **erf.les_type**                 | Using an LES       | "None",           | "None"      |
|                                  | model, and if so,  | "Smagorinsky",    |             |
|                                  | which type?        | "Deardorff"       |             |
+----------------------------------+--------------------+-------------------+-------------+
| **erf.molec_diff_type**          | Using molecular    | "None",           | "Constant", |
|                                  | viscosity and      | "Constant"        | but "None"  |
|                                  | diffusivity?       |                   | for LES     |
+----------------------------------+--------------------+-------------------+-------------+
| **erf.dynamicViscosity**         | Viscous coeff. if  | Real              | 0.0         |
|                                  | DNS                |                   |             |
+----------------------------------+--------------------+-------------------+-------------+
| **erf.Cs**                       | Constant           | Real              | 0.0         |
|                                  | Smagorinsky coeff. |                   |             |
+----------------------------------+--------------------+-------------------+-------------+
| **erf.Pr_t**                     | Turbulent Prandtl  | Real              | 1.0         |
|                                  | Number             |                   |             |
+----------------------------------+--------------------+-------------------+-------------+
| **erf.Sc_t**                     | Turbulent Schmidt  | Real              | 1.0         |
|                                  | Number             |                   |             |
+----------------------------------+--------------------+-------------------+-------------+
| **erf.spatial_order**            |                    | 1 / 2 / 4         | 2           |
+----------------------------------+--------------------+-------------------+-------------+


Forcing Terms
=============

.. _list-of-parameters-13:

List of Parameters
------------------

+----------------------------------+-------------------+-------------------+-------------+
| Parameter                        | Definition        | Acceptable        | Default     |
|                                  |                   | Values            |             |
+==================================+===================+===================+=============+
| **erf.abl_driver_type**          | Type of external  | None,             | None        |
|                                  | forcing term      | PressureGradient  |             |
|                                  |                   | GeostrophicWind   |             |
+----------------------------------+-------------------+-------------------+-------------+
| **erf.abl_pressure_grad**        | Pressure gradient | 3 Reals           | (0.,0.,0.)  |
|                                  | forcing term      |                   |             |
|                                  | (only if          |                   |             |
|                                  | abl.driver_type = |                   |             |
|                                  | PressureGradient) |                   |             |
+----------------------------------+-------------------+-------------------+-------------+
| **erf.abl_geo_wind**             | Geostrophic       | 3 Reals           | (0.,0.,0.)  |
|                                  | forcing term      |                   |             |
|                                  | (only if          |                   |             |
|                                  | abl.driver_type = |                   |             |
|                                  | GeostrophicWind)  |                   |             |
+----------------------------------+-------------------+-------------------+-------------+
| **erf.use_gravity**              | Include gravity   | true / false      | false       |
|                                  | in momentum       |                   |             |
|                                  | update?  If true, |                   |             |
|                                  | there is buoyancy |                   |             |
+----------------------------------+-------------------+-------------------+-------------+
| **erf.use_coriolis**             | Include Coriolis  | true / false      | false       |
|                                  | forcing           |                   |             |
+----------------------------------+-------------------+-------------------+-------------+
| **erf.use_rayleigh_damping **    | Include explicit  | true / false      | false       |
|                                  | Rayleigh damping  |                   |             |
+----------------------------------+-------------------+-------------------+-------------+

