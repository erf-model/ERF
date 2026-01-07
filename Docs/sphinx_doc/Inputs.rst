.. role:: cpp(code)
  :language: c++

.. _sec:Inputs:

******
Inputs
******
.. toctree::
   :maxdepth: 1

The ERF executable reads run-time information from an inputs file which you name on the command line.
This section describes the inputs which can be specified either in the inputs file or on the command line.
A value specified on the command line will override a value specified in the inputs file.

Governing Equations
===================
+--------------------------+-----------------------------+---------------+-------------+
| Parameter                | Definition                  | Acceptable    | Default     |
|                          |                             | Values        |             |
+==========================+=============================+===============+=============+
| **erf.anelastic**        | if 1, solve the anelastic   | 0, 1          | 0           |
|                          | equations (instead of       |               |             |
|                          | the compressible equations);|               |             |
|                          | can be a single value (for  |               |             |
|                          | all amr levels) or a list of|               |             |
|                          | values (one per level)      |               |             |
+--------------------------+-----------------------------+---------------+-------------+
| **erf.buoyancy_type**    | Controls how buoyancy is    | 1, 2, 3, 4    | 1 (may be   |
|                          | calculated: 1=density       |               | auto-set to |
|                          | perturbation, 2/3=temp.     |               | 2 or 3 based|
|                          | perturbation, 4=potential   |               | on moisture |
|                          | temp. perturbation.         |               | and solver) |
|                          | See :ref:`Buoyancy` for     |               |             |
|                          | details                     |               |             |
+--------------------------+-----------------------------+---------------+-------------+
| **erf.use_fft**          | use FFT rather than         | true / false  | false       |
|                          | multigrid to solve the      |               |             |
|                          | the Poisson equations       |               |             |
+--------------------------+-----------------------------+---------------+-------------+
| **erf.mg_v**             | verbosity of the multigrid  | Integer >= 0  | 0           |
|                          | solver if used              |               |             |
|                          | the Poisson equations       |               |             |
+--------------------------+-----------------------------+---------------+-------------+

.. note::

   To use the FFT solver if running with anelastic, you must set USE_FFT = TRUE if using gmake
   or ERF_ENABLE_FFT if using cmake.

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
Domain boundary conditions in ERF may be broadly categorized as ``ideal`` or ``real`` where
the ideal BC types correspond those found in classic fluid solvers and real correspond to an
external data source that may be based upon observation data. The inputs for these types of BCs
are given in more detail in :ref:`sec:LateralBoundaryConditions`. We briefly note that the supported
ideal BC types are: ``inflow``, ``outflow``, ``slipwall``, ``noslipwall``, ``symmetry`` or ``MOST``.


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

+-----------------------------------+--------------------+-----------------+-------------+
| Parameter                         | Definition         | Acceptable      | Default     |
|                                   |                    | Values          |             |
+===================================+====================+=================+=============+
| **amr.n_cell**                    | number of cells    | Integer > 0     | must be set |
|                                   | in each            |                 |             |
|                                   | direction at       |                 |             |
|                                   | the coarsest       |                 |             |
|                                   | level              |                 |             |
+-----------------------------------+--------------------+-----------------+-------------+
| **amr.max_level**                 | number of          | Integer >= 0    | must be set |
|                                   | levels of          |                 |             |
|                                   | refinement         |                 |             |
|                                   | above the          |                 |             |
|                                   | coarsest level     |                 |             |
+-----------------------------------+--------------------+-----------------+-------------+
| **amr.ref_ratio**                 | ratio of coarse    | Integer >= 1    | 2 for all   |
|                                   | to fine grid       | (one per level) | levels      |
|                                   | spacing between    |                 |             |
|                                   | subsequent         |                 |             |
|                                   | levels             |                 |             |
+-----------------------------------+--------------------+-----------------+-------------+
| **amr.ref_ratio_vect**            | ratio of coarse    | 3 integers >= 1 | 2 for all   |
|                                   | to fine grid       | (one per dir)   | directions  |
|                                   | spacing between    |                 |             |
|                                   | subsequent         |                 |             |
|                                   | levels             |                 |             |
+-----------------------------------+--------------------+-----------------+-------------+
| **amr.regrid_int**                | how often to       | Integer > 0     | -1          |
|                                   | regrid             | (if negative,   |             |
|                                   |                    | no regridding)  |             |
+-----------------------------------+--------------------+-----------------+-------------+
| **amr.regrid_on_restart**         | should we          | 0 or 1          | 0           |
|                                   | regrid             |                 |             |
|                                   | immediately        |                 |             |
|                                   | after              |                 |             |
|                                   | restarting         |                 |             |
+-----------------------------------+--------------------+-----------------+-------------+
| **amr.regrid_level_0_on_restart** | should we          | true or false   | false       |
|                                   | regrid level       |                 |             |
|                                   | immediately        |                 |             |
|                                   | after              |                 |             |
|                                   | restarting         |                 |             |
+-----------------------------------+--------------------+-----------------+-------------+
| **amr.iterate_grids**             | do we iterate      | true, false     | true        |
|                                   | on the grids?      |                 |             |
|                                   |                    |                 |             |
+-----------------------------------+--------------------+-----------------+-------------+

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


Grid Stretching
===============

This automatically activates **erf.terrain_type = StaticFittedMesh**. By default, the
problem-specific terrain is initialized to be flat at an elevation of z=0.
These inputs are used to automatically generate the staggered z levels used to
calculate the grid metric transformation. Alternatively, arbitrary z levels may
be specified with the **erf.terrain_z_levels** parameter, which should vary
from 0 (at the surface) to the domain height.

.. _list-of-parameters-3:

List of Parameters
------------------

+-------------------------------+-----------------+-----------------+-------------+
| Parameter                     | Definition      | Acceptable      | Default     |
|                               |                 | Values          |             |
+===============================+=================+=================+=============+
| **erf.grid_stretching_ratio** | scaling factor  | Real > 1        | 0 (no grid  |
|                               | applied to      |                 | stretching) |
|                               | delta z at each |                 |             |
|                               | level           |                 |             |
+-------------------------------+-----------------+-----------------+-------------+
| **erf.initial_dz**            | vertical grid   | Real > 0        | must be set |
|                               | spacing for the |                 | if grid     |
|                               | cell above the  |                 | stretching  |
|                               | bottom surface  |                 | ratio is set|
+-------------------------------+-----------------+-----------------+-------------+
| **erf.terrain_z_levels**      | nominal         | List of Real    | NONE        |
|                               | staggered       |                 |             |
|                               | z levels        |                 |             |
+-------------------------------+-----------------+-----------------+-------------+

.. _notes-3:

Notes
-----

- If both **erf.terrain_z_levels** and **erf.grid_stretching_ratio** are
  specified, the simple grid stretching will be ignored.
- The number of input **erf.terrain_z_levels** must be equal **amr.n_cell** in
  the z direction + 1.

.. _examples-of-usage-3:

Examples of Usage
-----------------

-  **erf.grid_stretching_ratio** = 1.025

-  | **erf.initial_dz** = 5.0
   | the first cell center would be at z=2.5


Regridding
==========

Overview
--------

The user defines how to tag individual cells at a given level for refinement.
This list of tagged cells is sent to a grid generation routine, which uses the
Berger-Rigoutsos algorithm to create rectangular grids that contain the tagged cells.

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

+---------------------+---------------------------+--------------+---------+
| Parameter           | Definition                | Acceptable   | Default |
|                     |                           | Values       |         |
+=====================+===========================+==============+=========+
| **max_step**        | maximum number of level 0 | Integer >= 0 | -1      |
|                     | time steps                |              |         |
+---------------------+---------------------------+--------------+---------+
| **start_time**      | starting simulation       | Real >= 0    |  0.0    |
|                     | time                      |              |         |
+---------------------+---------------------------+--------------+---------+
| **stop_time**       | final simulation          | Real >= 0    | Very    |
|                     | time                      |              | Large   |
+---------------------+---------------------------+--------------+---------+
| **start_datetime**  | starting simulation       | String       | None    |
|                     | date/time                 | (see notes)  |         |
+---------------------+---------------------------+--------------+---------+
| **stop_datetime**   | final simulation          | String       | None    |
|                     | date/time                 | (see notes)  |         |
+---------------------+---------------------------+--------------+---------+

.. _notes-3:

Notes
-----

- To control the number of time steps, you can limit by the maximum number
  of level-0 time steps (**max_step**), or the final simulation time
  (**stop_time**), or both. The code will stop at whichever criterion
  comes first. Note that if the code reaches **stop_time** then the final
  time step will be shortened so as to end exactly at **stop_time**, not
  pass it.
- For real data cases, the start and stop times is the epoch time in seconds.
- **start_datetime** and **stop_datetime** are in UTC and use the following
  strftime format: "%Y-%m-%d %H:%M:%S", e.g., "2001-01-01 18:00:00".
  The start/stop datetime inputs have precedence over the time inputs.

.. _examples-of-usage-4:

Examples of Usage
-----------------

-  **max_step** = 1000

-  **stop_time** = 1.0

will end the calculation when either the simulation time reaches 1.0 or
the number of level-0 steps taken equals 1000, whichever comes first.

Time Step
=========

The solver timestep can be fixed by the user or computed dynamically at each timestep based on the user-specified CFL
number --- i.e., adaptive time stepping. For the compressible equations, the timestep calculation uses the acoustic CFL constraint.
We note that when using implicit substepping, the vertical mesh spacing does not appear in the time step calculation.
The number of acoustic sub-steps per timestep can also be specified by the user as a fixed value or by specifying the
number of substeps per RK stage.  For the anelastic equations, the timestep calculation uses the advective CFL constraint,
which means it is determined by the fluid speed rather than the sound speed and thus allows much larger timesteps.

.. _list-of-parameters-6:

List of Parameters
------------------

+----------------------------+----------------------+----------------+---------------------+
| Parameter                  | Definition           | Acceptable     | Default             |
|                            |                      | Values         |                     |
+============================+======================+================+=====================+
| **erf.substepping_type**   | Should we substep in | "Implicit" or  | "Implicit" if       |
|                            | each RK stage?       | "None"         | compressible,       |
|                            |                      |                | "None" if anelastic |
+----------------------------+----------------------+----------------+---------------------+
| **erf.vert_implicit**      | Do vertical implicit | Boolean        | false               |
|                            | solve for diffusion  |                |                     |
|                            | of u, v, and theta   |                |                     |
|                            | with default         |                |                     |
|                            | time-centering in    |                |                     |
|                            | each stage           |                |                     |
+----------------------------+----------------------+----------------+---------------------+
| **erf.vert_implicit_fac**  | How much implicit    | Real >= 0      | 0.0, 0.0, 0.0       |
|                            | vertical diffusion   | (explicit) and | (fully explicit)    |
|                            | to include in each   | <= 1 (implicit)|                     |
|                            | RK stage? Currently, |                |                     |
|                            | only applies to      |                |                     |
|                            | rho*theta component. |                |                     |
|                            |                      |                |                     |
|                            | Specify either one   |                |                     |
|                            | (the same for all    |                |                     |
|                            | stages) or three     |                |                     |
|                            | values, one per stage|                |                     |
+----------------------------+----------------------+----------------+---------------------+
| **erf.cfl**                | CFL number used to   | Real > 0 and   | 0.8                 |
|                            | compute level 0 dt   | <= 1           |                     |
+----------------------------+----------------------+----------------+---------------------+
| **erf.substepping_cfl**    | CFL number used to   | Real > 0 and   | 1.0                 |
|                            | compute the number   | <= 1           |                     |
|                            | of substeps          |                |                     |
+----------------------------+----------------------+----------------+---------------------+
| **erf.fixed_dt**           | set level 0 dt       | Real > 0       | unused if not       |
|                            | as this value        |                | set                 |
|                            | regardless of        |                |                     |
|                            | cfl or other         |                |                     |
|                            | settings             |                |                     |
+----------------------------+----------------------+----------------+---------------------+
| **erf.fixed_fast_dt**      | set fast dt          | Real > 0       |                     |
|                            | as this value        |                |                     |
+----------------------------+----------------------+----------------+---------------------+
| **erf.fixed_mri_dt_ratio** | set fast dt          | even int > 0   | only relevant if    |
|                            | as slow dt /         |                | substepping_type    |
|                            | this ratio           |                | is not None         |
+----------------------------+----------------------+----------------+---------------------+
| **erf.init_shrink**        | factor by which      | Real > 0 and   | 1.0                 |
|                            | to shrink the        | <= 1           |                     |
|                            | initial dt           |                |                     |
+----------------------------+----------------------+----------------+---------------------+
| **erf.change_max**         | factor by which      | Real >= 1      | 1.1                 |
|                            | dt can grow          |                |                     |
|                            | in subsequent        |                |                     |
|                            | steps                |                |                     |
+----------------------------+----------------------+----------------+---------------------+
| **erf.dt_max**             | maximum adaptive     | Real > 0       | 1e9                 |
|                            | timestep             |                |                     |
|                            | allowed by time      |                |                     |
|                            | stepping             |                |                     |
+----------------------------+----------------------+----------------+---------------------+
| **erf.dt_max_initial**     | maximum initial      | Real > 0       | 1.0                 |
|                            | timestep             |                |                     |
+----------------------------+----------------------+----------------+---------------------+
| **erf.dt_ref_ratio**       | ratio of coarse      | Integer >= 1   | same as             |
|                            | to fine grid         | (one per level)| maximum over        |
|                            | time steps between   |                | directions of       |
|                            | subsequent           |                | ref_ratio           |
|                            | levels               |                |                     |
+----------------------------+----------------------+----------------+---------------------+

Notes
-----------------

-  | If **erf.anelastic** is true then **substepping_type** is internally set to "None".

-  | The time step controls work somewhat differently depending on whether one is using
     acoustic substepping in time; this is determined by the value of **substepping_type**.

-  | If **erf.substepping_type = None** there is only one time step to be calculated,
     and **fixed_fast_dt** and **fixed_mri_dt_ratio** are not used.

   * | If **erf.fixed_dt** is also specified, the timestep will be set to **fixed_dt**.

   * | If **erf.fixed_dt** is not specified, the timestep will be computed using the CFL condition for compressible flow.
       If **erf.cfl** is specified, that CFL value will be used.  If not, the default value will be used.

-  | If **erf.substepping_type** is not **None** we must determine both the slow and fast timesteps.
   * | If **erf.fixed_dt** is specified, the slow timestep will be set to **fixed_dt**.

   * | If **erf.fixed_dt** is not set, the slow timestep will be computed using the CFL
       condition for incompressible/anelastic flow.  If **erf.cfl** is specified, that CFL value will be used.
       If not, the default value will be used.

   * | There are several consistency checks before the fast timestep is computed.  Specifically, if any
       of the following are true the code will abort while reading the inputs.

     * | If **erf.fixed_mri_dt_ratio** is specified but is not an even positive integer
     * | If **erf.fixed_dt** and **erf.fast_fixed_dt** are specified and the ratio of **fixed_dt** to **fast_fixed_dt**
         is not an even positive integer
     * | If **erf.fixed_dt** and **erf.fast_fixed_dt** and **erf.fixed_mri_dt_ratio** are all specified but are inconsistent

   * | Once the slow timestep is set and the inputs are allowed per the above criteria,
       the fast timestep is computed in one of several ways:

     * | If **erf.fixed_fast_dt** is specified, the fast timestep will be set to **fixed_fast_dt**.

     * | If **erf.fixed_mri_dt_ratio** is specified and **erf.fixed_fast_dt** is not specified,
         the fast timestep will be the slow timestep divided by **fixed_mri_dt_ratio.**

     * | If neither **erf.fixed_mri_dt_ratio** nor **erf.fixed_fast_dt** is specified, then the fast timestep
         will be computed using the CFL condition for compressible flow, then adjusted (reduced if necessary)
         as above so that the ratio of slow timestep to fine timestep is an even integer.
         If **erf.substepping_cfl** is specified, that CFL value will be used.  If not, the default value will be used.

.. _examples-of-usage-5:

Examples of Usage of Additional Parameters
-------------------------------------------

-  | **erf.init_shrink** = 0.01
   | sets the initial slow time step to 1% of what it would be otherwise.
     Note that if **erf.init_shrink** :math:`\neq 1` and **fixed_dt** is specified,
     then the first time step will in fact be **erf.init_shrink** \* **erf.fixed_dt**.

-  | **erf.change_max** = 1.1
   | allows the slow time step to increase by no more than 10% in this case.
     Note that the time step can shrink by any factor; this only
     controls the extent to which it can grow.

Restart Capability
==================

See :ref:`sec:Checkpoint` for how to control the checkpoint/restart capability.

PlotFiles
===============================

See :ref:`sec:Plotfiles` for how to control the types and frequency of plotfile
generation.


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
     and print these quantities.

Diagnostic Outputs
==================

If ``erf.v`` is set then one or more additional output files may be requested.
These include (1) a surface time history file, (2) a history of mean profiles,
(3) a history of vertical flux profiles (i.e., variances and covariances), and
(4) a history of modeled subgrid stresses. The number of specified output
filenames dictates the level of output. E.g., specifying 3 filenames will give
outputs (1), (2), and (3). Data files are only written if ``erf.profile_int >
0``. This output functionality has not been implemented for terrain.

.. _list-of-parameters-10a:


List of Parameters
------------------

+-------------------------------+------------------+----------------+----------------+
| Parameter                     | Definition       | Acceptable     | Default        |
|                               |                  | Values         |                |
+===============================+==================+================+================+
| **erf.data_log**              | Output           | Up to four     | NONE           |
|                               | filename(s)      | strings        |                |
+-------------------------------+------------------+----------------+----------------+
| **erf.der_data_log**          | Output           | Up to four     | NONE           |
|                               | filename(s) for  | strings        |                |
|                               | derived data     |                |                |
+-------------------------------+------------------+----------------+----------------+
| **erf.energy_data_log**       | Output           | Up to four     | NONE           |
|                               | filename(s) for  | strings        |                |
|                               | total energy     |                |                |
+-------------------------------+------------------+----------------+----------------+
| **erf.profile_int**           | Interval (number)| Integer        | -1             |
|                               | of steps between |                |                |
|                               | outputs          |                |                |
+-------------------------------+------------------+----------------+----------------+
| **erf.destag_profiles**       | Interpolate all  | Boolean        | true           |
|                               | outputs to       |                |                |
|                               | cell-center      |                |                |
|                               | heights          |                |                |
+-------------------------------+------------------+----------------+----------------+

By default, all profiles are planar-averaged quantities :math:`\langle\cdot\rangle`
that are destaggered by interpolating to cell centers where appropriate.
Setting ``erf.destag_profiles = false`` will
keep vertically staggered quantities on z faces -- quantities already at cell
centers or on x/y faces will remain at those locations. Note that all output
quantities -- whether cell-centered or face-centered -- will be output on the
staggered grid. The user should discard values at the highest z level
(corresponding to the z-dir of ``amr.n_cell`` + 1) for destaggered quantities.
Staggered quantities are indicated below.

The requested output files have the following columns:


* Surface time history

  #. Time (s)

  #. Friction velocity, :math:`u_*` (m/s)

  #. Surface-layer potential temperature scale, :math:`\theta_*` (K)

  #. Obukhov length, :math:`L` (m)

* Mean flow profiles

  #. Time (s)

  #. Height (m)

  #. X-velocity, :math:`\langle u \rangle` (m/s)

  #. Y-velocity, :math:`\langle v \rangle` (m/s)

  #. *Z-velocity*, :math:`\langle w \rangle` (m/s) -- *staggered*

  #. Dry air density, :math:`\langle \rho \rangle` (kg/m3)

  #. Total (moist) potential temperature, :math:`\langle \theta \rangle` (K)

  #. Turbulent kinetic energy (TKE), :math:`\langle k \rangle` (m2/s2) from the Deardorff subgrid model only

  #. Eddy diffusivity, :math:`\langle K \rangle` (kg/(m-s)) from the subgrid model

  #. Water vapor mixing ratio, :math:`\langle q_v \rangle` (kg/kg) from the microphysics model

  #. Cloud water mixing ratio, :math:`\langle q_c \rangle` (kg/kg) from the microphysics model

  #. Rain water mixing ratio, :math:`\langle q_r \rangle` (kg/kg) from the microphysics model, if available

  #. Ice mixing ratio, :math:`\langle q_i \rangle` (kg/kg) from the microphysics model, if available

  #. Snow mixing ratio, :math:`\langle q_s \rangle` (kg/kg) from the microphysics model, if available

  #. Graupel mixing ratio, :math:`\langle q_g \rangle` (kg/kg) from the microphysics model, if available

* Vertical flux profiles

  #. Time (s)

  #. Height (m)

  #. X-velocity variance, :math:`\langle u^\prime u^\prime \rangle` (m2/s2)

  #. X,Y-velocity covariance, :math:`\langle u^\prime v^\prime \rangle` (m2/s2)

  #. *X,Z-velocity covariance*, :math:`\langle u^\prime w^\prime \rangle` (m2/s2) -- *staggered*

  #. Y-velocity variance, :math:`\langle v^\prime v^\prime \rangle` (m2/s2)

  #. *Y,Z-velocity covariance*, :math:`\langle v^\prime w^\prime \rangle` (m2/s2) -- *staggered*

  #. *Z-velocity variance*, :math:`\langle w^\prime w^\prime \rangle` (m2/s2) -- *staggered*

  #. X-direction heat flux, :math:`\langle u^\prime \theta^\prime \rangle` (K m/s)

  #. Y-direction heat flux, :math:`\langle v^\prime \theta^\prime \rangle` (K m/s)

  #. *Z-direction heat flux*, :math:`\langle w^\prime \theta^\prime \rangle` (K m/s) -- *staggered*

  #. Temperature variance, :math:`\langle \theta^\prime \theta^\prime \rangle` (K m/s)

  #. X-direction turbulent transport of TKE, :math:`\langle u_i^\prime u_i^\prime u^\prime \rangle` (m3/s3)
     -- Note: :math:`u_i u_i = uu + vv + ww`

  #. Y-direction turbulent transport of TKE, :math:`\langle u_i^\prime u_i^\prime v^\prime \rangle` (m3/s3)

  #. *Z-direction turbulent transport of TKE*, :math:`\langle u_i^\prime u_i^\prime w^\prime \rangle` (m3/s3) -- *staggered*

  #. X-direction pressure transport of TKE, :math:`\langle p^\prime u^\prime \rangle` (m3/s3)

  #. Y-direction pressure transport of TKE, :math:`\langle p^\prime v^\prime \rangle` (m3/s3)

  #. *Z-direction pressure transport of TKE*, :math:`\langle p^\prime w^\prime \rangle` (m3/s3) -- *staggered*

  #. *Z-direction flux of water vapor*, :math:`\langle q_v^\prime w^\prime \rangle` (m/s) -- *staggered*

  #. *Z-direction flux of cloud water*, :math:`\langle q_c^\prime w^\prime \rangle` (m/s) -- *staggered*

  #. *Z-direction flux of rain water*, :math:`\langle q_r^\prime w^\prime \rangle` (m/s) -- *staggered*

  #. *Z-direction virtual temperature flux*, :math:`\langle w^\prime \theta_v^\prime \rangle` (K m/s) -- *staggered*

* Modeled subgrid-scale (SGS) profiles

  #. SGS stress tensor component, :math:`\tau_{11}` (m2/s2)

  #. SGS stress tensor component, :math:`\tau_{12}` (m2/s2)

  #. *SGS stress tensor component*, :math:`\tau_{13}` (m2/s2) -- *staggered*

  #. SGS stress tensor component, :math:`\tau_{22}` (m2/s2)

  #. *SGS stress tensor component*, :math:`\tau_{23}` (m2/s2) -- *staggered*

  #. SGS stress tensor component, :math:`\tau_{33}` (m2/s2)

  #. *SGS heat flux*, :math:`\tau_{\theta w}` (K m/s) -- *staggered*

  #. *SGS water vapor flux*, :math:`\tau_{q_v w}` (K m/s) -- *staggered*

  #. *SGS cloud water flux*, :math:`\tau_{q_c w}` (K m/s) -- *staggered*

  #. SGS turbulence dissipation, :math:`\epsilon` (m2/s3)


Data Sampling Outputs
==================

Data along query lines or planes may be output during the simulation if
``erf.do_line_sampling = true`` or  ``erf.do_plane_sampling = true``, respectively.
The potential temperature and wind-speed will be written to native ``plt_line/plane``
at the step frequency dictated by ``erf.sampler_interval = <int>``. For line sampling,
users must prescribe ``sample_line_lo`` and ``sample_line_hi`` inputs which are 3 integer
values corresponding to the (i,j,k) indices at the beginning and end of the line.
Additionally, users must specify ``sample_line_dir`` to prescribed the direction of
the line. The same inputs are used for the plane sampling except that ``sample_plane_lo/hi``
must be the physical locations of the plane corners. This output functionality has
not been implemented for terrain. By default, sampled line and plane data will have the
prefixes "plt_line" and "plt_plane", respectively. Names for sampled data may optionally
be provided with ``sample_line_name`` and/or ``sample_plane_name`` -- if provided, each
line and/or plane must be named.

Line and plane samples will be default be written to plotfiles, one plotfile per output
snapshot, with all output variables in the same file. Alternatively, line sampling has
the ``erf.line_sampling_text_output`` option, which writes one text file per output
variable, with all snapshots appended to the same file over time. This is similar to the
tslist output from WRF but output is provided only from the finest domain that contains
the entire requested sampling line; velocities are also destaggered.

The sampled variables can be selected with the ``erf.line_sampling_vars`` option and
includes a subset of the plotfile outputs: "density", "x_velocity", "y_velocity", "z_velocity",
"magvel", "theta", "qv", "qc", and "pressure". Velocities are output at cell centers only.
The water vapor mixing ratio "qv" will only output valid values if a moisture model is used.
Pressure is calculated from rho*theta and will account for moisture if qv is requested.

.. _list-of-parameters-10b:


List of Parameters
------------------

+-----------------------------------+------------------+----------------+----------------+
| Parameter                         | Definition       | Acceptable     | Default        |
|                                   |                  | Values         |                |
+===================================+==================+================+================+
| **erf.line_sampling_interval**,   | Output           | Integer        | -1             |
| **erf.plane_sampling_interval**   | frequency (steps)|                |                |
+-----------------------------------+------------------+----------------+----------------+
| **erf.line_sampling_per**,        | Output           | Real           | -1             |
| **erf.plane_sampling_per**        | frequency (time) |                |                |
+-----------------------------------+------------------+----------------+----------------+
| **erf.do_line_sampling**          | Flag to do line  | Boolean        | false          |
|                                   | sampling         |                |                |
|                                   |                  |                |                |
+-----------------------------------+------------------+----------------+----------------+
| **erf.do_plane_sampling**         | Flag to do plane | Boolean        | false          |
|                                   | sampling         |                |                |
|                                   |                  |                |                |
+-----------------------------------+------------------+----------------+----------------+
| **erf.sample_line_dir**           | Directionality   | Integer        | None           |
|                                   | of the line      |                |                |
+-----------------------------------+------------------+----------------+----------------+
| **erf.sample_plane_dir**          | Directionality   | Integer        | None           |
|                                   | of the plane     |                |                |
+-----------------------------------+------------------+----------------+----------------+
| **erf.sample_line_lo/hi**         | Bounding (i,j,k) | 3 Integers per | None           |
|                                   | on the line(s)   | line           |                |
+-----------------------------------+------------------+----------------+----------------+
| **erf.sample_plane_lo/hi**        | Bounding point   | 3 Reals per    | None           |
|                                   | on the plane(s)  | plane          |                |
+-----------------------------------+------------------+----------------+----------------+
| **erf.sample_line_name**          | Output prefix of | One string per | None           |
|                                   | each line        | sample line    |                |
+-----------------------------------+------------------+----------------+----------------+
| **erf.sample_plane_name**         | Output prefix of | One string per | None           |
|                                   | each plane       | sample plane   |                |
+-----------------------------------+------------------+----------------+----------------+
| **erf.line_sampling_text_output** | Write text files | Boolean        | false          |
|                                   | instead of AMReX |                |                |
|                                   | plotfiles        |                |                |
+-----------------------------------+------------------+----------------+----------------+
| **erf.line_sampling_vars**        | Specify sampled  | List of strings| theta, magvel  |
|                                   | variables        |                |                |
+-----------------------------------+------------------+----------------+----------------+

.. _examples-of-usage-10b:

Example of Usage
-----------------

::

   erf.do_line_sampling = true               # Do line sampling
   erf.line_sampling_interval = 1            # Write plt files every step
   erf.sample_line_lo   = 5 32  5 10   32 5  # Lo points for two lines
   erf.sample_line_hi   = 5 32 25 1000 32 5  # Hi points for two lines
   erf.sample_line_dir  = 2 0                # One line in z and one in x

   erf.do_plane_sampling = true              # Do plane sampling
   erf.plane_sampling_interval = 10          # Write plt files every 10 steps
   erf.sample_plane_lo   =  48.0  48.0  32.0 # Lo points for one plane
   erf.sample_plane_hi   = 320.0 320.0  32.0 # Hi points for one plane
   erf.sample_plane_dir  = 2                 # One plane with z normal


Advection Schemes
=================

.. _list-of-parameters-11:

List of Parameters
------------------

+----------------------------------+--------------------+---------------------+--------------+
| Parameter                        | Definition         | Acceptable          | Default      |
|                                  |                    | Values              |              |
+==================================+====================+=====================+==============+
| **erf.dycore_horiz_adv_type**    | Horizontal         | see below           | Upwind_3rd   |
|                                  | advection type     |                     |              |
|                                  | for dycore vars    |                     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.dycore_vert_adv_type**     | Vertical           | see below           | Upwind_3rd   |
|                                  | advection type     |                     |              |
|                                  | for dycore vars    |                     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.dryscal_horiz_adv_type**   | Horizontal         | see below           | Upwind_3rd   |
|                                  | advection type     |                     |              |
|                                  | for dry scalars    |                     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.dryscal_vert_adv_type**    | Vertical           | see below           | Upwind_3rd   |
|                                  | advection type     |                     |              |
|                                  | for dry scalars    |                     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.moistscal_horiz_adv_type** | Horizontal         | see below           | WENO3        |
|                                  | advection type     |                     |              |
|                                  | for moist scalars  |                     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.moistscal_vert_adv_type**  | Vertical           | see below           | WENO3        |
|                                  | advection type     |                     |              |
|                                  | for moist scalars  |                     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.use_efficient_advection**  | Use efficient      | true/false          | false        |
|                                  | advection scheme   |                     |              |
|                                  | for scalars        |                     |              |
+----------------------------------+--------------------+---------------------+--------------+

The allowed advection types for the dycore variables are:

* From WRF:
  "Centered_2nd",
  "Upwind_3rd",
  "Centered_4th",
  "Upwind_5th", and
  "Centered_6th"
* Blended schemes:
  "Blended_3rd4th" and
  "Blended_5th6th"
* WENO-JS schemes:
  "WENO3"
  "WENO5", and
  "WENO7"
* WENO-Z schemes:
  "WENOZ3"
  "WENOZ5", and
  "WENOZ7"

In addition to the above advection schemes for dycore variables, dry and moist
scalars also have the option of using "Upwind_3rd_SL" (slope-limited) and
"WENOMZQ3."

To use the blended schemes, ``erf.[vartype]_[horiz|vert]_upw_frac`` for the corresponding
``erf.[vartype]_[horiz|vert]_adv_type`` should be set to a scalar between 0 and 1, where 1 is fully
upwind and 0 is fully central.

Note: if using WENO schemes, the horizontal and vertical advection types must be set to
the same string.

The efficient advection schemes for dry and moist scalars exploit the substages of the
time advancing RK3 scheme by using lower order schemes in the first two substages and the
solver's choice of scheme in the final stage. Based on CPU-only runtimes on Perlmutter for
the scalar advection routine, the approximate computational savings for the scalar advection
schemes are as follows when using efficient advection option: roughly 30% for Centered_4th
and Centered_6th, 35% for Upwind_5th, roughly 45% for WENO5 and WENOZ5, and roughly 60% for
Upwind_3rd, WENO3, WENOZ3, and WENOMZQ3.


Diffusive Physics
=================

.. _list-of-parameters-12:

List of Parameters
------------------

+----------------------------------+--------------------+---------------------+--------------+
| Parameter                        | Definition         | Acceptable          | Default      |
|                                  |                    | Values              |              |
+==================================+====================+=====================+==============+
| **erf.alpha_T**                  | Diffusion coeff.   | Real                | 0.0          |
|                                  | for temperature    |                     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.alpha_C**                  | Diffusion coeff.   | Real                | 0.0          |
|                                  | for scalar         |                     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.rho0_trans**               | Reference density  | Real                | 1.0          |
|                                  | to compute const.  |                     |              |
|                                  | rho*Alpha          |                     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.les_type**                 | Using an LES       | "None",             | "None"       |
|                                  | model, and if so,  | "Smagorinsky",      |              |
|                                  | which type?        | "Deardorff"         |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.rans_type**                | Using a RANS       | "None" or "kEqn"    | "None"       |
|                                  | model, and if so,  |                     |              |
|                                  | which type?        |                     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.molec_diff_type**          | Using molecular    | "None",             | "None"       |
|                                  | viscosity and      | "Constant", or      |              |
|                                  | diffusivity?       | "ConstantAlpha"     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.dynamic_viscosity**        | Viscous coeff. if  | Real                | 0.0          |
|                                  | DNS                |                     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.Cs**                       | Constant           | Real                | 0.0          |
|                                  | Smagorinsky coeff. |                     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.use_moist_Ri_correction**  | Apply moist        | Boolean             | false        |
|                                  | Richardson number  |                     |              |
|                                  | limiter to the     |                     |              |
|                                  | Smagorinsky model  |                     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.Ck**                       | Constant           | Real                | 0.1          |
|                                  | Deardorff k coeff. |                     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.Ce**                       | Constant           | Real                | 0.93         |
|                                  | Deardorff epsilon  |                     |              |
|                                  | coeff.             |                     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.Ce_wall**                  | Constant           | Real                | 0.0          |
|                                  | Deardorff epsilon  |                     |              |
|                                  | coeff. at the wall;|                     |              |
|                                  | if > 0, then set   |                     |              |
|                                  | Ce to this at k=0  |                     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.sigma_k**                  | Constant Deardorff | Real                | 0.5          |
|                                  | coeff. in          |                     |              |
|                                  | downgradient       |                     |              |
|                                  | diffusion term     |                     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.theta_ref**                | Reference potential| Real                | 0.0          |
|                                  | temperature used   |                     |              |
|                                  | to characterize    |                     |              |
|                                  | stable             |                     |              |
|                                  | stratficiation;    |                     |              |
|                                  | constant if > 0,   |                     |              |
|                                  | otherwise the      |                     |              |
|                                  | instantaneous local|                     |              |
|                                  | value is used      |                     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.Pr_t**                     | Turbulent Prandtl  | Real                | 1.0          |
|                                  | Number             |                     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.Sc_t**                     | Turbulent Schmidt  | Real                | 1.0          |
|                                  | Number             |                     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.num_diff_coeff**           | Coefficient for    | Real                | 0.0          |
|                                  | 6th-order          | [0.0,  1.0]         |              |
|                                  | numerical          |                     |              |
|                                  | diffusion, set to 0|                     |              |
|                                  | to disable         |                     |              |
+----------------------------------+--------------------+---------------------+--------------+

Note: in the equations for the evolution of momentum, potential temperature and advected scalars, the
diffusion coefficients are written as :math:`\mu`, :math:`\rho \alpha_T` and :math:`\rho \alpha_C`, respectively.

If we set ``erf.molec_diff_type`` to ``Constant``, then

- ``erf.dynamic_viscosity`` is used as the value of :math:`\mu` in the momentum equation, and

- ``erf.alpha_T`` is multiplied by ``erf.rho0_trans`` to form the coefficient for potential temperature, and

- ``erf.alpha_C`` is multiplied by ``erf.rho0_trans`` to form the coefficient for an advected scalar.

If we set ``erf.molec_diff_type`` to ``ConstantAlpha``, then

- the dynamic viscosity in the momentum equation is assumed to have the form :math:`\mu = \rho \alpha_M`
  where :math:`\alpha_M` is a momentum diffusivity constant with units of kinematic viscosity, calculated as
  ``erf.dynamic_viscosity`` divided by ``erf.rho0_trans``;
  this diffusivity is multiplied by the instantaneous local density :math:`\rho` to form the coefficient in the momentum equation; and

- ``erf.alpha_T`` is multiplied by the instantaneous local density :math:`\rho` to form the coefficient for potential temperature, and

- ``erf.alpha_C`` is multiplied by the instantaneous local density :math:`\rho` to form the coefficient for an advected scalar.

Parameters for LES can either be set with one value that applies across all levels, or set with a number of values
equal to the number of levels, allowing unique values of the parameter to be set for each level.

PBL Scheme
==========

.. _list-of-parameters-13:

List of Parameters
------------------

+-----------------------------------------+--------------------+---------------------+-------------+
| Parameter                               | Definition         | Acceptable          | Default     |
|                                         |                    | Values              |             |
+=========================================+====================+=====================+=============+
| **erf.pbl_type**                        | Name of PBL Scheme | "None", "MYNN25",   | "None"      |
|                                         | to be used         | "MYNNEDMF", "MYJ",  |             |
|                                         |                    | "YSU", "MRF",       |             |
|                                         |                    | "SHOC"              |             |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_mynn_A1**                     | MYNN Constant A1   | Real                | 1.18        |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_mynn_A2**                     | MYNN Constant A2   | Real                | 0.665       |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_mynn_B1**                     | MYNN Constant B1   | Real                | 24.0        |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_mynn_B2**                     | MYNN Constant B2   | Real                | 15.0        |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_mynn_C1**                     | MYNN Constant C1   | Real                | 0.137       |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_mynn_C2**                     | MYNN Constant C1   | Real                | 0.75        |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_mynn_C3**                     | MYNN Constant C3   | Real                | 0.352       |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_mynn_C4**                     | MYNN Constant C4   | Real                | 0.0         |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_mynn_C5**                     | MYNN Constant C5   | Real                | 0.2         |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_mynn_SQfactor**               | MYNN ratio of      | Real                | 3.0         |
|                                         | stability functions|                     |             |
|                                         | SQ / SM            |                     |             |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_mynn_diffuse_moistvars**      | Diffuse moisture   | bool                | 0           |
|                                         | variables using    |                     |             |
|                                         | modeled eddy       |                     |             |
|                                         | diffusivity        |                     |             |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.advect_QKE**                      | Include advection  | bool                | 1           |
|                                         | terms in QKE eqn   |                     |             |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.diffuse_QKE_3D**                  | Include horizontal | bool                | 0           |
|                                         | turb. diffusion    |                     |             |
|                                         | terms in QKE eqn.  |                     |             |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_ysu_force_over_water**        | Treat whole domain | bool                | 0           |
|                                         | as over water for  |                     |             |
|                                         | YSU PBL scheme     |                     |             |
|                                         | regardless of      |                     |             |
|                                         | LSM/other inputs   |                     |             |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_ysu_land_Ribcr**              | Over land critical | Real                | 0.25        |
|                                         | Richardson number  |                     |             |
|                                         | for YSU PBL Scheme |                     |             |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_ysu_unst_Ribcr**              | Unstable critical  | Real                | 0.0         |
|                                         | Richardson number  |                     |             |
|                                         | for YSU PBL Scheme |                     |             |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_ysu_coriolis_freq**           | Coriolis frq. used | Real                | 1.0e-4      |
|                                         | for YSU PBL Scheme |                     |             |
|                                         | (1e-4 in WRF)      |                     |             |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_ysu_use_consistent_coriolis** | Ignore above param | Bool                | 0           |
|                                         | and use the value  |                     |             |
|                                         | from ERF coriolis  |                     |             |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_mrf_coriolis_freq**           | Coriolis frq. used | Real                | 1.0e-4      |
|                                         | for MRF PBL Scheme |                     |             |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_mrf__Ribcr**                  | Over land critical | Real                | 0.5         |
|                                         | Richardson number  |                     |             |
|                                         | for MRF PBL Scheme |                     |             |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_mrf_const_b**                 | Coefficient for the| Real                | 7.8         |
|                                         | countergradient    |                     |             |
|                                         | term               |                     |             |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_mrf_sf**                      | ratio of surface   | Real                | 0.1         |
|                                         | layer height to    |                     |             |
|                                         | boundary layer     |                     |             |
|                                         | height             |                     |             |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.mrf_moistvars**                   | Diffuse moisture   | bool                | 0           |
|                                         | variables using    |                     |             |
|                                         | modeled eddy       |                     |             |
|                                         | diffusivity        |                     |             |
+-----------------------------------------+--------------------+---------------------+-------------+

Note that both PBL schemes must be used in conjunction with a MOST boundary condition
at the surface (Zlo) boundary. The YSU scheme is work in progress currently.

If the PBL scheme is activated, it determines the turbulent diffusivity in the vertical
direction. If an LES model is also specified, it determines only the horizontal turbulent
diffusivity.

Right now, the QKE equation is solved if and only if the MYNN2.5 PBL model is selected. In that
transport equation, it is optional to advect QKE, and to apply LES diffusive transport for QKE
in the horizontal directions (the vertical component is always computed as part of the PBL
scheme).

Parameters for PBL schemes can either be set with one value that applies across all levels, or set with a number of values
equal to the number of levels, allowing unique values of the parameter to be set for each level.

Forcing Terms
=============

.. _list-of-parameters-14:

List of Parameters
------------------

+-------------------------------------+------------------------+-------------------+---------------------+
| Parameter                           | Definition             | Acceptable        | Default             |
|                                     |                        | Values            |                     |
+=====================================+========================+===================+=====================+
| **erf.abl_driver_type**             | Type of external       | None,             | None                |
|                                     | forcing term           | PressureGradient  |                     |
|                                     |                        | GeostrophicWind   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.abl_pressure_grad**           | Pressure gradient      | 3 Reals           | (0.,0.,0.)          |
|                                     | forcing term           |                   |                     |
|                                     | (only if               |                   |                     |
|                                     | abl.driver_type =      |                   |                     |
|                                     | PressureGradient)      |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.abl_geo_wind**                | Geostrophic            | 3 Reals           | (0.,0.,0.)          |
|                                     | forcing term           |                   |                     |
|                                     | (only if               |                   |                     |
|                                     | abl.driver_type =      |                   |                     |
|                                     | GeostrophicWind)       |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.abl_geo_wind_table**          | Path to text file      | String            | None                |
|                                     | containing a           |                   |                     |
|                                     | geostrophic wind       |                   |                     |
|                                     | profile                |                   |                     |
|                                     | (with z, Ug, and       |                   |                     |
|                                     |  Vg whitespace         |                   |                     |
|                                     |  delimited             |                   |                     |
|                                     |  columns)              |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.windfarm_type**               | Wind farm              | "None",           | "None"              |
|                                     | parameterization       | "Fitch", "EWP",   |                     |
|                                     |                        | "SimpleActuator", |                     |
|                                     |                        | "GeneralActuator" |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.const_massflux_u**            | Include a momentum     | Real              | 0.                  |
| **erf.const_massflux_v**            | source at each time,   |                   |                     |
|                                     | (e.g., representing a  |                   |                     |
|                                     | background driving     |                   |                     |
|                                     | pressure gradient),    |                   |                     |
|                                     | to obtain a desired    |                   |                     |
|                                     | mass flux with the     |                   |                     |
|                                     | specified bulk velocity|                   |                     |
|                                     | in x,y                 |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.const_massflux_layer_lo**     | Two heights defining   | Real              | None                |
| **erf.const_massflux_layer_hi**     | the layer over which   |                   |                     |
|                                     | the mass flux is       |                   |                     |
|                                     | integrated and compared|                   |                     |
|                                     | to the desired input(s)|                   |                     |
|                                     | specified above        |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.const_massflux_tau**          | Timescale over which   | Real              | None                |
|                                     | to adjust the          |                   |                     |
|                                     | background pressure    |                   |                     |
|                                     | gradient to match the  |                   |                     |
|                                     | specified mass flux    |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.use_gravity**                 | Include gravity        | true / false      | false               |
|                                     | in momentum            |                   |                     |
|                                     | update?  If true,      |                   |                     |
|                                     | there is buoyancy      |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.use_coriolis**                | Include Coriolis       | true / false      | false               |
|                                     | forcing                |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.variable_coriolis**           | Include Coriolis       | true / false      | false               |
|                                     | forcing that varies    |                   |                     |
|                                     | with latitude          |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.rotational_time_period**      | Used to calculate the  | Real              | 86400.0             |
|                                     | Coriolis frequency     |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.latitude**                    | Used to calculate the  | Real              | 90.0                |
|                                     | Coriolis frequency     |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.coriolis_3d**                 | Include z component in | true / false      | true                |
|                                     | the Coriolis forcing   |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.rayleigh_damping_type**       | Rayleigh damping       | "None",           | "SlowExplicit"      |
|                                     | type                   | "SlowExplicit",   |                     |
|                                     |                        | "FastExplicit"    |                     |
|                                     |                        | "FastImplicit"    |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.rayleigh_damp_U**             | Include explicit       | true / false      | false               |
|                                     | Rayleigh damping in    |                   |                     |
|                                     | the x-momentum equation|                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.rayleigh_damp_V**             | Include explicit       | true / false      | false               |
|                                     | Rayleigh damping in    |                   |                     |
|                                     | the y-momentum equation|                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.rayleigh_damp_W**             | Include                | true / false      | false               |
|                                     | Rayleigh damping in    |                   |                     |
|                                     | the z-momentum equation|                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.rayleigh_damp_T**             | Include explicit       | true / false      | false               |
|                                     | Rayleigh damping in    |                   |                     |
|                                     | the potential          |                   |                     |
|                                     | temperature equation   |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.rayleigh_dampcoef**           | Rayleigh damping       | Real              | 0.2                 |
|                                     | coefficient, an inverse|                   |                     |
|                                     | timescale              |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.rayleigh_zdamp**              | Rayleigh damping       | Real              | 500.0               |
|                                     | layer depth measured   |                   |                     |
|                                     | from the top of domain |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.nudging_from_input_sounding** | Add momentum source    | true / false      | false               |
|                                     | terms to nudge the     |                   |                     |
|                                     | solution towards the   |                   |                     |
|                                     | initial sounding       |                   |                     |
|                                     | profile                |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.input_sounding_file**         | Name(s) of the         | String(s)         | input_sounding_file |
|                                     | input sounding file(s) |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.forest_file**                 | Name(s) of the         | String            | None                |
|                                     | canopy forest file     |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.input_sounding_time**         | Time(s) of the         | Real(s)           | false               |
|                                     | input sounding file(s) |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.tau_nudging**                 | Time scale for         | Real              | 5.0                 |
|                                     | nudging                |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+

If ``erf.nudging_from_input_sounding`` is true, it is expected that at least one input sounding
file is available.  If there is only one, and no specification of time is made, it is assumed that
the one file corresponds to time = 0.0.   If the final time supplied in
``input_*_sounding_*_time``  is less than the final time in the calculation, the final sounding supplied
in ``input_*_sounding_*_file`` will be used for all times later than the final value in
in ``input_*_sounding_*_time``.

In addition, custom forcings or tendencies may be defined on a problem-specific
basis. This affords additional flexibility in defining the RHS source term as
a function of time and/or height. Implementation entails modifying problem
source code inside the Exec directory and overriding the ``update_*_sources()``
function(s).

+--------------------------------------------+-------------------+-------------------+-------------+
| Parameter                                  | Definition        | Acceptable        | Default     |
|                                            |                   | Values            |             |
+============================================+===================+===================+=============+
| **erf.custom_forcing_uses_primitive_vars** | User-defined      | true or false     | false       |
|                                            | source terms set  |                   |             |
|                                            | the tendency of   |                   |             |
|                                            | primitive         |                   |             |
|                                            | variables instead |                   |             |
|                                            | of conserved      |                   |             |
|                                            | quantities        |                   |             |
|                                            | (rho*prim_var)    |                   |             |
+--------------------------------------------+-------------------+-------------------+-------------+
| **erf.add_custom_rhotheta_forcing**        | Apply the         | true or false     | false       |
|                                            | user-defined      |                   |             |
|                                            | temperature source|                   |             |
|                                            | term              |                   |             |
+--------------------------------------------+-------------------+-------------------+-------------+
| **erf.add_custom_moisture_forcing**        | Apply the         | true or false     | false       |
|                                            | user-defined      |                   |             |
|                                            | qv source         |                   |             |
|                                            | term              |                   |             |
+--------------------------------------------+-------------------+-------------------+-------------+
| **erf.add_custom_w_subsidence**            | Apply the         | true or false     | false       |
|                                            | user-defined      |                   |             |
|                                            | vertical velocity |                   |             |
|                                            | profile for use in|                   |             |
|                                            | calculating       |                   |             |
|                                            | subsidence source |                   |             |
|                                            | terms             |                   |             |
+--------------------------------------------+-------------------+-------------------+-------------+
| **erf.add_custom_geostrophic_profile**     | Apply the         | true or false     | false       |
|                                            | user-defined      |                   |             |
|                                            | geostrophic wind  |                   |             |
|                                            | profile           |                   |             |
+--------------------------------------------+-------------------+-------------------+-------------+

Note that ``erf.add_custom_geostrophic_profile`` cannot be used in combination
with an ``erf.abl_geo_wind_table``.

- Wind farm parameterization requires ``USE_WINDFARM=TRUE`` at build time. See :ref:`WindFarmModels` for theory and examples.


Numerical Stability
===================

To enhance numerical stability, e.g., for operational runs, we provide some additional controls.

List of Parameters
------------------

+----------------------------+----------------------------------+-------------------+-------------+
| **erf.beta_s**             | Time off-centering coefficient   | Real              | 0.1         |
+----------------------------+----------------------------------+-------------------+-------------+
| **erf.w_damping**          | Enable vertical-velocity         | Bool              | false       |
|                            | damping                          |                   |             |
+----------------------------+----------------------------------+-------------------+-------------+
| **erf.w_damping_cfl**      | Critical vertical advective      | Real              | 1.0         |
|                            | CFL for w-damping to be applied, |                   |             |
|                            | if ``erf.w_damping`` is true     |                   |             |
+----------------------------+----------------------------------+-------------------+-------------+
| **erf.w_damping_const**    | w-damping coefficient (m/s2)     | Real              | -1          |
+----------------------------+----------------------------------+-------------------+-------------+
| **erf.w_damping_coeff**    | w-damping coefficient (-)        | Real              | -1          |
+----------------------------+----------------------------------+-------------------+-------------+

If ``erf.w_damping`` is true, then either ``erf.w_damping_const`` or ``erf.w_damping_coeff`` must
be specified. For WRF-like damping, set ``erf.w_damping_const = 0.3`` to give

.. math::

   f_d = - \rho\,\text{sgn}(w)\,(C-C_{crit}) \cdot \gamma,

where :math:`\gamma` is a constant damping coefficient with units of m/s2 and the advective Courant
number is :math:`C = |w|\Delta t / \Delta z`.

If ``erf.w_damping_coeff`` is specified instead, then

.. math::

   f_d = - \rho\,\text{sgn}(w)\,(|w|-w_{crit}) \cdot \frac{\alpha}{\Delta t},

where :math:`\alpha` is a dimensionless damping factor and :math:`w_{crit} = C_{crit} \Delta z / \Delta t`.
This is equivalent to:

.. math::

   f_d = - \rho\,\text{sgn}(w)\,(C-C_{crit}) \cdot \alpha \frac{\Delta z}{(\Delta t)^2}.

This approach gives a damping coefficient that is sensitive to the vertical grid spacing and
robustly damps towards the critical Courant number.


Initialization
==============

The initialization in ERF has two steps: creation of the background state and creation of initial perturbations from the background state.

The initialization strategy is determined at runtime by ``init_type``, which has six possible values.

See :ref:`sec:Initialization` for more detail about how to provide initial conditions for an ERF simulation.

In addition, there is a run-time option to project the initial velocity field to make it divergence-free.

List of Parameters
------------------

+----------------------------------+---------------------+--------------------+-----------------------+
| Parameter                        | Definition          | Acceptable         | Default               |
|                                  |                     | Values             |                       |
+==================================+=====================+====================+=======================+
| **erf.init_type**                | Initialization      | "None",            | "None"                |
|                                  | type                | "WRFInput",        |                       |
|                                  |                     | "Input_Sounding"   |                       |
|                                  |                     | "Metgrid"          |                       |
|                                  |                     | "NCFile"           |                       |
|                                  |                     | "Uniform"          |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.input_sounding_file**      | Path to WRF-style   |  String            | "input_sounding"      |
|                                  | input sounding      |                    |                       |
|                                  | file                |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.sounding_type**            | How to interpret    | "Ideal",           | Ideal                 |
|                                  | the sounding        | "Isentropic",      |                       |
|                                  | provided with       | "DryIsentropic",   |                       |
|                                  | init_type =         | "ConstantDensity"  |                       |
|                                  | "input_sounding"    |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.use_real_bcs**             | If init_type is     | true or false      | true if               |
|                                  | WRFInput or Metgrid |                    | if init_type          |
|                                  | do we want to use   |                    | is WRFInput or        |
|                                  | these bcs?          |                    | Metgrid;              |
|                                  |                     |                    | else false            |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.nc_init_file**             | NetCDF file with    |  String            | NONE                  |
|                                  | initial mesoscale   |                    |                       |
|                                  | data                |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.nc_bdy_file**              | NetCDF file with    |  String            | NONE                  |
|                                  | mesoscale data at   |                    |                       |
|                                  | lateral boundaries  |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.project_initial_velocity** | project initial     |  0 or 1            | 1 if anelastic;       |
|                                  | velocity?           |                    | 0 if compressible     |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.real_width**               | Lateral boundary    |  Integer           | 0                     |
|                                  | total width if      |                    |                       |
|                                  | use_real_bcs is     |                    |                       |
|                                  | true                |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.real_set_width**           | Lateral boundary    |  Integer           | 0                     |
|                                  | specified zone      |                    |                       |
|                                  | width if            |                    |                       |
|                                  | use_real_bcs is     |                    |                       |
|                                  | true                |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.real_extrap_w**            | First-order         | bool               | true                  |
|                                  | extrapolation of    |                    |                       |
|                                  | vertical velocities |                    |                       |
|                                  | on lateral          |                    |                       |
|                                  | boundaries (instead |                    |                       |
|                                  | of setting to 0) if |                    |                       |
|                                  | use_real_bcs is     |                    |                       |
|                                  | true                |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.metgrid_debug_quiescent**  | If init_type is     | true or false      | false                 |
|                                  | Metgrid, overwrite  |                    |                       |
|                                  | initial conditions  |                    |                       |
|                                  | and boundary        |                    |                       |
|                                  | conditions to be    |                    |                       |
|                                  | quiescent.          |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.metgrid_debug_isothermal** | If init_type is     | true or false      | false                 |
|                                  | Metgrid, overwrite  |                    |                       |
|                                  | theta to be 300 in  |                    |                       |
|                                  | initial conditions  |                    |                       |
|                                  | and boundary        |                    |                       |
|                                  | conditions.         |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.metgrid_debug_dry**        | If init_type is     | true or false      | false                 |
|                                  | Metgrid, overwrite  |                    |                       |
|                                  | qv to be dry in     |                    |                       |
|                                  | initial conditions  |                    |                       |
|                                  | and boundary        |                    |                       |
|                                  | conditions.         |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.metgrid_debug_msf**        | If init_type is     | true or false      | false                 |
|                                  | Metgrid, overwrite  |                    |                       |
|                                  | map scale factors   |                    |                       |
|                                  | to be 1.            |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.metgrid_debug_psfc**       | If init_type is     | true or false      | false                 |
|                                  | Metgrid, overwrite  |                    |                       |
|                                  | surface pressure    |                    |                       |
|                                  | to be 10**5.        |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.metgrid_interp_theta**     | If init_type is     | true or false      | false                 |
|                                  | Metgrid, calculate  |                    |                       |
|                                  | theta on origin     |                    |                       |
|                                  | model vertical      |                    |                       |
|                                  | levels and then     |                    |                       |
|                                  | interpolate onto    |                    |                       |
|                                  | the ERF vertical    |                    |                       |
|                                  | levels.             |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.metgrid_basic_linear**     | If init_type is     | true or false      | false                 |
|                                  | Metgrid, use        |                    |                       |
|                                  | linear vertical     |                    |                       |
|                                  | interpolation and   |                    |                       |
|                                  | no quality          |                    |                       |
|                                  | control?            |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.metgrid_use_below_sfc**    | If init_type is     | true or false      | true                  |
|                                  | Metgrid, use the    |                    |                       |
|                                  | origin data levels  |                    |                       |
|                                  | below the surface?  |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.metgrid_use_sfc**          | If init_type is     | true or false      | true                  |
|                                  | Metgrid, use the    |                    |                       |
|                                  | origin data level   |                    |                       |
|                                  | at the surface?     |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.metgrid_retain_sfc**       | If init_type is     | true or false      | false                 |
|                                  | Metgrid, assign     |                    |                       |
|                                  | the lowest level    |                    |                       |
|                                  | directly using the  |                    |                       |
|                                  | surface value from  |                    |                       |
|                                  | the origin data?    |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.metgrid_proximity**        | If init_type is     | Real               | 1000.                 |
|                                  | Metgrid, pressure   |                    |                       |
|                                  | differential for    |                    |                       |
|                                  | detecting origin    |                    |                       |
|                                  | levels that are     |                    |                       |
|                                  | problematically     |                    |                       |
|                                  | close together      |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.metgrid_order**            | If init_type is     | Integer            | 2                     |
|                                  | Metgrid, order of   |                    |                       |
|                                  | the Lagrange        |                    |                       |
|                                  | polynomial          |                    |                       |
|                                  | interpolation       |                    |                       |
|                                  | scheme for          |                    |                       |
|                                  | vertical            |                    |                       |
|                                  | interpolation       |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.metgrid_force_sfc_k**      | If init_type is     | Integer            | 0                     |
|                                  | Metgrid, force the  |                    |                       |
|                                  | origin data         |                    |                       |
|                                  | surface level to    |                    |                       |
|                                  | be included in the  |                    |                       |
|                                  | interpolation for   |                    |                       |
|                                  | this many ERF       |                    |                       |
|                                  | vertical levels     |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+

Notes
-----------------

If **erf.init_type = WRFInput** or **erf.init_type = NCFile**,
the problem is initialized with mesoscale data contained in a NetCDF file,
provided via ``erf.nc_init_file`` (e.g., "wrfinput_d01").

In addition, if **erf.use_real_bcs = true**, the lateral boundary conditions must be supplied in a NetCDF files
specified by ``erf.nc_bdy_file`` (e.g., "wrfbdy_d01").  (If **erf.use_real_bcs = false**, no file is read for the
boundary conditions so they must be specified in the inputs file.)

If **erf.use_real_bcs = true**,
the extent of the relaxation zone may be controlled with ``erf.real_width`` (corresponding to WRF's **spec_bdy_width**)
and ``erf.real_set_width`` (corresponding to WRF's **spec_zone**, typically set to 1), which corresponds to a relaxation zone with a
width of **real_width - real_set_width**.

Note that the **erf.project_initial_velocity** option is available for all **init_type** options.  If using the anelastic
formulation this will be true regardless of the input; if using the compressible formulation the default is false but
that can be over-written.

Map Scale Factors
=================

Map scale factors are always present in the evolution equations, but the values default to one
unless specified in the initialization when **erf.init_type = real**.

There is an option to test the map scale factors by setting  **erf.test_mapfactor = true**; this
arbitrarily sets the map factors to 0.5 in order to test the implementation.

Terrain
=======

ERF allows the use to specify whether terrain-fitted coordinates should be used by
setting **erf.terrain_type = StaticFittedMesh** (static mesh) or
**erf.terrain_type = MovingFittedMesh** (time-dependent mesh).
If using terrain-fitted coordinates, the user also has the option to specify one of three
methods for defining how the terrain-fitted coordinates given the topography:

- Basic Terrain Following (BTF):
    The influence of the terrain decreases linearly with height.
-  Smoothed Terrain Following (STF):
    Small-scale terrain structures are progressively smoothed out of the coordinate system as height increases.
- Sullivan Terrain Following (name TBD):
    The influence of the terrain decreases with the cube of height.

The user can also specify that terrain should be represented with an immersed forcing method
(with an optional wall model, see :ref:`Forcings` for more detail), or
with an embedded boundary / cut cell representation.

.. note:: The embedded boundary / cut cell representation is a work in progress and not ready for use!

The height at the surface nodes can be defined analytically or read from a text file
specified by ``erf.terrain_file_name``.  The ordering of data in the file is first
nx, ny (where nx is the number of values specified in the x-direction and
ny is the number of values specified in the y-direction; note that these need not match
the resolution of the problem).  Then nx x-values are read followed by ny y-values.   Finally,
the z-coordinates are read in the order z(x1,y1), z(x1,y2), z(x1,y3), ... z(x2,y1), ... z(nx,ny).
An example is given in Exec/ABL/erf_terrain_def.

List of Parameters
------------------

+-----------------------------+--------------------+--------------------+------------+
| Parameter                   | Definition         | Acceptable         | Default    |
|                             |                    | Values             |            |
+=============================+====================+====================+============+
| **erf.terrain_type**        | Is there terrain   | None               | None       |
|                             | and if so, how is  | StaticFittedMesh   |            |
|                             | it represented?    | MovingFittedMesh   |            |
|                             |                    | ImmersedForcing    |            |
|                             |                    | EB                 |            |
+-----------------------------+--------------------+--------------------+------------+
| **erf.terrain_smoothing**   | specify terrain    | 0,                 | 0          |
|                             | following          | 1,                 |            |
|                             |                    | 2                  |            |
+-----------------------------+--------------------+--------------------+------------+
| **erf.terrain_file_name**   | filename           | String             | NONE       |
+-----------------------------+--------------------+--------------------+------------+

Examples of Usage
-----------------

-  **erf.terrain_smoothing**  = 0
    BTF is used when generating the terrain following coordinate.

-  **erf.terrain_smoothing**  = 1
    STF is used when generating the terrain-fitted coordinate. Additionally,
    ``erf.terrain_gamma_m`` (default=0.5) may be used to set the minimum
    allowable fractional grid spacing. From Klemp 2011, MWR: "Values of 0.5-0.6
    seem to work best in 2D applications, while values about half this
    magnitude appear better for 3D real-terrain simulations."

-  **erf.terrain_smoothing**  = 2
    Sullivan TF is used when generating the terrain following coordinate.

-  When setting **erf.terrain_file_name**, the format of the file is expected to
    be (in raw text): the integer nx on the first line, the integer ny on the second line,
    then the nx values of the x-coordinate, then the ny values of
    the y-coordinate, then the (nx times ny) values of the z-coordinate associated
    with the (x,y) values we have just read in.  Note that the z-values are in the
    order z(x1,y1), z(x1,y2), z(x1,y3), ... which is contrary to standard Fortran ordering

Land Surface Model
==================

The land surface model provides energy and moisture fluxes at the lower boundary.

List of Parameters
------------------

+--------------------------------+----------------------------+--------------------+-------------+
| Parameter                      | Definition                 | Acceptable         | Default     |
|                                |                            | Values             |             |
+================================+============================+====================+=============+
| **erf.land_surface_model**     | Enables land surface       | "None",            | "None"      |
|                                | energy and moisture        | "NOAHMP",          |             |
|                                | fluxes                     | "MM5", "SLM"       |             |
+--------------------------------+----------------------------+--------------------+-------------+

.. note::

   Noah-MP requires ``USE_NOAHMP=TRUE`` at build time. See :ref:`CouplingToNoahMP` for details.

Coupling Type (Data Exchange)
==============================

For coupled simulations with AMR-Wind or WaveWatch3, this controls the directionality of data exchange.


List of Parameters
------------------

+--------------------------------+----------------------------+--------------------+-------------+
| Parameter                      | Definition                 | Acceptable         | Default     |
|                                |                            | Values             |             |
+================================+============================+====================+=============+
| **erf.coupling_type**          | Coupling mode for          | "OneWay",          | "None"      |
|                                | coupled simulations with   | "TwoWay"           |             |
|                                | AMR-Wind or WaveWatch3     |                    |             |
+--------------------------------+----------------------------+--------------------+-------------+

Notes
-----

- **TwoWay**: Bidirectional coupling - ERF both receives and sends data
- **OneWay**: ERF receives forcing but doesn't send data back

See :ref:`CouplingToAMRWind` and :ref:`CouplingToWW3` for more information.

Moisture
========

ERF has several different moisture models. The models that are currently implemented
are Eulerian models; however, ERF has the capability for Lagrangian models when
compiled with particles.

The following run-time options control how the full moisture model is used.

List of Parameters
------------------

+---------------------------------+--------------------------+-----------------------+------------+
| Parameter                       | Definition               | Acceptable            | Default    |
|                                 |                          | Values                |            |
+=================================+==========================+=======================+============+
| **erf.moisture_model**          | Name of moisture model   |  "None", "SAM",       | "None"     |
|                                 |                          |  "Kessler", "SatAdj"  |            |
|                                 |                          |  "Kessler_NoRain",    |            |
|                                 |                          |  "Morrison",          |            |
|                                 |                          |  "Morrison_NoIce",    |            |
|                                 |                          |  "SAM_NoPrecip_NoIce",|            |
|                                 |                          |  "SAM_NoIce", "P3"    |            |
+---------------------------------+--------------------------+-----------------------+------------+
| **erf.moisture_tight_coupling** | If true, advance         |  Bool                 | false      |
|                                 | microphysics after every |                       |            |
|                                 | slow step in the dycore; |                       |            |
|                                 | otherwise, update after  |                       |            |
|                                 | the dycore has been      |                       |            |
|                                 | advanced at each timestep|                       |            |
+---------------------------------+--------------------------+-----------------------+------------+

Radiation
=========

ERF allows for radiative heating computations with the RRTMGP library.
If building with cmake, the following flags must be enabled:
``-DERF_ENABLE_RRTMGP:BOOL=ON``, ``-DERF_ENABLE_NETCDF:BOOL=ON``, and ``-DERF_ENABLE_HDF5:BOOL=ON``;
see **ERF/Build/cmake_with_radiation.sh**.

If building with gmake, set ``USE_RRTMGP = TRUE`` and ``USE_NETCDF = TRUE`` in the GNUmakefile.

Notes
-----------------

-  | A rule of thumb for the radiation update frequency is 1 min per km of grid spacing (e.g., every 10 min on a 10-km grid)

-  | For idealized studies, constant latitude/longitude may be specified through **erf.rad_cons_lat**
   | and **erf.rad_cons_lon**.

List of Parameters
------------------

+--------------------------------+--------------------------+--------------------+-----------------------------------+
| Parameter                      | Definition               | Acceptable         | Default                           |
|                                |                          | Values             |                                   |
+================================+==========================+====================+===================================+
| **erf.radiation_model**        | Enables radiation        | "None",            | "None"                            |
|                                | transfer calculations    | "RRTMGP"           |                                   |
+--------------------------------+--------------------------+--------------------+-----------------------------------+
| **erf.rad_freq_in_steps**      | Number of steps between  |  int               |    1                              |
+--------------------------------+--------------------------+--------------------+-----------------------------------+
| **erf.rad_write_fluxes**       | Flag to write fluxes     |  Bool              |    false                          |
+--------------------------------+--------------------------+--------------------+-----------------------------------+
| **erf.rad_cons_lat**           | Constant latitude        |  Real              |    39.8                           |
+--------------------------------+--------------------------+--------------------+-----------------------------------+
| **erf.rad_cons_lon**           | Constant longitude       |  Real              |    -98.5                          |
+--------------------------------+--------------------------+--------------------+-----------------------------------+
| **erf.co2vmr**                 | CO2 volume mixing ratio  |  Real              | 388.717e-6                        |
+--------------------------------+--------------------------+--------------------+-----------------------------------+
| **erf.o3vmr**                  | O3 volume mixing ratio   |  Real              |   1.887e-7                        |
+--------------------------------+--------------------------+--------------------+-----------------------------------+
| **erf.n2ovmr**                 | N2O volume mixing ratio  |  Real              | 323.141e-9                        |
+--------------------------------+--------------------------+--------------------+-----------------------------------+
| **erf.covmr**                  | CO volume mixing ratio   |  Real              |   1.000e-7                        |
+--------------------------------+--------------------------+--------------------+-----------------------------------+
| **erf.ch4vmr**                 | CH4 volume mixing ratio  |  Real              |   1.807e-6                        |
+--------------------------------+--------------------------+--------------------+-----------------------------------+
| **erf.o2vmr**                  | O2 volume mixing ratio   |  Real              |   0.209                           |
+--------------------------------+--------------------------+--------------------+-----------------------------------+
| **erf.n2vmr**                  | N2 volume mixing ratio   |  Real              |   0.791                           |
+--------------------------------+--------------------------+--------------------+-----------------------------------+
| **erf.rrtmgp_file_path**       | path to NC files         |  String            |   "./"                            |
+--------------------------------+--------------------------+--------------------+-----------------------------------+
| **erf.rrtmgp_coeffs_sw**       | path to NC files         |  String            | rrtmgp-data-sw-g224-2018-12-04.nc |
+--------------------------------+--------------------------+--------------------+-----------------------------------+
| **erf.rrtmgp_coeffs_lw**       | path to NC files         |  String            | rrtmgp-data-lw-g256-2018-12-04.nc |
+--------------------------------+--------------------------+--------------------+-----------------------------------+
| **erf.rrtmgp_cloud_optics_sw** | path to NC files         |  String            | rrtmgp-cloud-optics-coeffs-sw.nc  |
+--------------------------------+--------------------------+--------------------+-----------------------------------+
| **erf.rrtmgp_cloud_optics_lw** | path to NC files         |  String            | rrtmgp-cloud-optics-coeffs-lw.nc  |
+--------------------------------+--------------------------+--------------------+-----------------------------------+

The lookup data may be downloaded as a package from `here <https://doi.org/10.22002/ppv8a-4q131>`_.

.. note::

   Using RRTMGP requires ``USE_RRTMGP=TRUE`` at build time. See :ref:`building` for build instructions.

Runtime Error Checking
======================

Through `AMReX functionality <https://amrex-codes.github.io/amrex/docs_html/Debugging.html>`_,
ERF supports options to raise errors when NaNs, division by zero, and overflow errors
are detected. These checks are activated at runtime using the input parameters below.

.. note::

   When running on Macs using the Apple-Clang compilers with optimization
   (``DEBUG = FALSE`` in the ``GNUmakefile``), these checks may lead to false positives
   due to optimizations performed by the compiler and the flags should be turned off.
   It is still possible to run with these error checks with Apple-Clang debug builds.

List of Parameters
------------------

+-----------------------------+---------------------------+-------------------+------------+
| Parameter                   | Definition                | Acceptable Values | Default    |
+=============================+===========================+===================+============+
| **erf.check_for_nans**      | Test solution for NaNs    |  int              | 0          |
+-----------------------------+---------------------------+-------------------+------------+
| **amrex.fpe_trap_invalid**  | Raise errors for NaNs     |  0 / 1            | 0          |
+-----------------------------+---------------------------+-------------------+------------+
| **amrex.fpe_trap_zero**     | Raise errors for divide   |  0 / 1            | 0          |
|                             | by zero                   |                   |            |
+-----------------------------+---------------------------+-------------------+------------+
| **amrex.fpe_trap_overflow** | Raise errors for overflow |  0 / 1            | 0          |
+-----------------------------+---------------------------+-------------------+------------+
