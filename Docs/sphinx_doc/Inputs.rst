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
| **amr.regrid_int**        | how often to    | Integer > 0     | -1          |
|                           | regrid          | (if negative,   |             |
|                           |                 | no regridding)  |             |
+---------------------------+-----------------+-----------------+-------------+
| **amr.regrid_on_restart** | should we       | 0 or 1          | 0           |
|                           | regrid          |                 |             |
|                           | immediately     |                 |             |
|                           | after           |                 |             |
|                           | restarting      |                 |             |
+---------------------------+-----------------+-----------------+-------------+
| **amr.iterate_grids**     | do we iterate   | true, false     | true        |
|                           | on the grids?   |                 |             |
|                           |                 |                 |             |
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


Grid Stretching
===============

This automatically activates **erf.use_terrain**. By default, the
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

+-----------------+---------------------------+--------------+---------+
| Parameter       | Definition                | Acceptable   | Default |
|                 |                           | Values       |         |
+=================+===========================+==============+=========+
| **max_step**    | maximum number of level 0 | Integer >= 0 | -1      |
|                 | time steps                |              |         |
+-----------------+---------------------------+--------------+---------+
| **start_time**  | starting simulation       | Real >= 0    |  0.0    |
|                 | time                      |              |         |
+-----------------+---------------------------+--------------+---------+
| **stop_time**   | final simulation          | Real >= 0    | Very    |
|                 | time                      |              | Large   |
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

+----------------------------+----------------------+----------------+-------------------+
| Parameter                  | Definition           | Acceptable     | Default           |
|                            |                      | Values         |                   |
+============================+======================+================+===================+
| **erf.no_substepping**     | Should we turn off   | int (0 or 1)   | 0                 |
|                            | substepping in time? |                |                   |
+----------------------------+----------------------+----------------+-------------------+
| **erf.cfl**                | CFL number for       | Real > 0 and   | 0.8               |
|                            | hydro                | <= 1           |                   |
|                            |                      |                |                   |
|                            |                      |                |                   |
+----------------------------+----------------------+----------------+-------------------+
| **erf.fixed_dt**           | set level 0 dt       | Real > 0       | unused if not     |
|                            | as this value        |                | set               |
|                            | regardless of        |                |                   |
|                            | cfl or other         |                |                   |
|                            | settings             |                |                   |
+----------------------------+----------------------+----------------+-------------------+
| **erf.fixed_fast_dt**      | set fast dt          | Real > 0       | only relevant     |
|                            | as this value        |                | if use_native_mri |
|                            |                      |                | is true           |
+----------------------------+----------------------+----------------+-------------------+
| **erf.fixed_mri_dt_ratio** | set fast dt          | even int > 0   | only relevant     |
|                            | as slow dt /         |                | if no_substepping |
|                            | this ratio           |                | is 0              |
+----------------------------+----------------------+----------------+-------------------+
| **erf.init_shrink**        | factor by which      | Real > 0 and   | 1.0               |
|                            | to shrink the        | <= 1           |                   |
|                            | initial dt           |                |                   |
+----------------------------+----------------------+----------------+-------------------+
| **erf.change_max**         | factor by which      | Real >= 1      | 1.1               |
|                            | dt can grow          |                |                   |
|                            | in subsequent        |                |                   |
|                            | steps                |                |                   |
+----------------------------+----------------------+----------------+-------------------+

Notes
-----------------

-  | The time step controls work somewhat differently depending on whether one is using
     acoustic substepping in time; this is determined by the value of **no_substepping**.

-  | If **erf.no_substepping = 1** there is only one time step to be calculated,
     and **fixed_fast_dt** and **fixed_mri_dt_ratio** are not used.

   * | If **erf.fixed_dt** is also specified, the timestep will be set to **fixed_dt**.

   * | If **erf.fixed_dt** is not specified, the timestep will be computed using the CFL condition for compressible flow.
       If **erf.cfl** is specified, that CFL value will be used.  If not, the default value will be used.

-  | If **erf.no_substepping = 0** we must determine both the slow and fast timesteps.
   * | If **erf.fixed_dt** is specified, the slow timestep will be set to **fixed_dt**.

   * | If **erf.fixed_dt** is not set, the slow timestep will be computed using the CFL
       condition for incompressible flow.  If **erf.cfl** is specified, that CFL value will be used.
       If not, the default value will be used.

   * | There are several consistency checks before the fast timestep is computed.  Specifically, if any
       of the following are true the code will abort while reading the inputs.

     * | If **erf.fixed_mri_dt_ratio** is specified but is not an even positive integer
     * | If **erf.fixed_dt** and **erf.fast_fixed_dt** are specified and the ratio of **fixed_dt** to **fast_fixed_dt**
         is not an even positive integer
     * | If **erf.fixed_dt** and **erf.fast_fixed_dt** and **erf.fixed_mri_dt_ratio** are all specified but are inconsitent

   * | Once the slow timestep is set and the inputs are allowed per the above criteria,
       the fast timestep is computed in one of several ways:

     * | If **erf.fixed_fast_dt** is specified, the fast timestep will be set to **fixed_fast_dt**.

     * | If **erf.fixed_mri_dt_ratio** is specified and **erf.fixed_fast_dt** is not specified,
         the fast timestep will be the slow timestep divided by **fixed_mri_dt_ratio.**

     * | If neither **erf.fixed_mri_dt_ratio** nor **erf.fixed_fast_dt** is specified, then the fast timestep
         will be computed using the CFL condition for compressible flow, then adjusted (reduced if necessary)
         as above so that the ratio of slow timestep to fine timestep is an even integer.
         If **erf.cfl** is specified, that CFL value will be used.  If not, the default value will be used.

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
| **erf.datalog**               | Output           | Up to four     | NONE           |
|                               | filename(s)      | strings        |                |
+-------------------------------+------------------+----------------+----------------+
| **erf.profile_int**           | Interval (number)| Integer        | -1             |
|                               | of steps between |                |                |
|                               | ouputs           |                |                |
+-------------------------------+------------------+----------------+----------------+
| **erf.interp_profiles_to_cc** | Interpolate all  | Boolean        | true           |
|                               | outputs to cell  |                |                |
|                               | centers          |                |                |
+-------------------------------+------------------+----------------+----------------+

By default, all profiles are planar-averaged quantities :math:`\langle\cdot\rangle`
interpolated to cell centers. Setting ``erf.interp_profiles_to_cc = false`` will
keep vertically staggered quantities on z faces (quantities already at cell
centers or on x/y faces will remain at those locations). Note that all output
quantities--whether cell-centered or face-centered--will be output on the
staggered grid. The user should discard the highest z level (corresponding to
the z-dir ``amr.n_cell`` + 1) for cell-centered quantities. Staggered quantities
are indicated below.

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

  #. Turbulent kinetic energy (TKE), :math:`\langle k \rangle` (m2/s2) for the subgrid model

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

* Modeled subgrid-scale (SGS) profiles

  #. SGS stress tensor component, :math:`\tau_{11}` (m2/s2)

  #. SGS stress tensor component, :math:`\tau_{12}` (m2/s2)

  #. *SGS stress tensor component*, :math:`\tau_{13}` (m2/s2) -- *staggered*

  #. SGS stress tensor component, :math:`\tau_{22}` (m2/s2)

  #. *SGS stress tensor component*, :math:`\tau_{23}` (m2/s2) -- *staggered*

  #. SGS stress tensor component, :math:`\tau_{33}` (m2/s2)

  #. SGS heat flux, :math:`\tau_{\theta w}` (K m/s)

  #. SGS turbulence dissipation, :math:`\epsilon` (m2/s3)


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

The allowed advection types for the dycore variables are
"Centered_2nd", "Upwind_3rd", "Blended_3rd4th", "Centered_4th", "Upwind_5th", "Blended_5th6th",
and "Centered_6th".

The allowed advection types for the dry and moist scalars are
"Centered_2nd", "Upwind_3rd", "Blended_3rd4th", "Centered_4th", "Upwind_5th", "Blended_5th6th",
"Centered_6th" and in addition,
"WENO3", "WENOZ3", "WENOMZQ3", "WENO5", and "WENOZ5."

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
| **erf.molec_diff_type**          | Using molecular    | "None",             | "None"       |
|                                  | viscosity and      | "Constant", or      |              |
|                                  | diffusivity?       | "ConstantAlpha"     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.dynamicViscosity**         | Viscous coeff. if  | Real                | 0.0          |
|                                  | DNS                |                     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.Cs**                       | Constant           | Real                | 0.0          |
|                                  | Smagorinsky coeff. |                     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.Pr_t**                     | Turbulent Prandtl  | Real                | 1.0          |
|                                  | Number             |                     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.Sc_t**                     | Turbulent Schmidt  | Real                | 1.0          |
|                                  | Number             |                     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.use_NumDiff**              | Use 6th order      | "true",             | "false"      |
|                                  | numerical diffusion| "false"             |              |
|                                  |                    |                     |              |
+----------------------------------+--------------------+---------------------+--------------+
| **erf.NumDiffCoeff**             | Coefficient for    | Real                | 0.0          |
|                                  | 6th order          | [0.0,  1.0]         |              |
|                                  | numerical diffusion|                     |              |
+----------------------------------+--------------------+---------------------+--------------+

Note: in the equations for the evolution of momentum, potential temperature and advected scalars, the
diffusion coefficients are written as :math:`\mu`, :math:`\rho \alpha_T` and :math:`\rho \alpha_C`, respectively.

If we set ``erf.molec_diff_type`` to ``Constant``, then

- ``erf.dynamicViscosity`` is used as the value of :math:`\mu` in the momentum equation, and

- ``erf.alpha_T`` is multiplied by ``erf.rho0_trans`` to form the coefficient for potential temperature, and

- ``erf.alpha_C`` is multiplied by ``erf.rho0_trans`` to form the coefficient for an advected scalar.

If we set ``erf.molec_diff_type`` to ``ConstantAlpha``, then

- the dynamic viscosity in the momentum equation is assumed to have the form :math:`\mu = \rho \alpha_M`
  where :math:`\alpha_M` is a momentum diffusivity constant with units of kinematic viscosity, calculated as
  ``erf.dynamicViscosity`` divided by ``erf.rho0_trans``;
  this diffusivity is multiplied by the instantaneous local density :math:`\rho` to form the coefficient in the momentum equation; and

- ``erf.alpha_T`` is multiplied by the instantaneous local density :math:`\rho` to form the coefficient for potential temperature, and

- ``erf.alpha_C`` is multiplied by the instantaneous local density :math:`\rho` to form the coefficient for an advected scalar.


PBL Scheme
==========

.. _list-of-parameters-13:

List of Parameters
------------------

+----------------------------------+--------------------+---------------------+-------------+
| Parameter                        | Definition         | Acceptable          | Default     |
|                                  |                    | Values              |             |
+==================================+====================+=====================+=============+
| **erf.pbl_type**                 | Name of PBL Scheme | "None", "MYNN2.5"   | "None"      |
|                                  | to be used         |                     |             |
+----------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_A1**                   | MYNN Constant A1   | Real                | 1.18        |
+----------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_A2**                   | MYNN Constant A2   | Real                | 0.665       |
+----------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_B1**                   | MYNN Constant B1   | Real                | 24.0        |
+----------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_B2**                   | MYNN Constant B2   | Real                | 15.0        |
+----------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_C1**                   | MYNN Constant C1   | Real                | 0.137       |
+----------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_C2**                   | MYNN Constant C1   | Real                | 0.75        |
+----------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_C3**                   | MYNN Constant C3   | Real                | 0.352       |
+----------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_C4**                   | MYNN Constant C4   | Real                | 0.0         |
+----------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_C5**                   | MYNN Constant C5   | Real                | 0.2         |
+----------------------------------+--------------------+---------------------+-------------+
| **erf.advect_QKE**               | Include advection  | bool                | 1           |
|                                  | terms in QKE eqn   |                     |             |
+----------------------------------+--------------------+---------------------+-------------+
| **erf.diffuse_QKE_3D**           | Include horizontal | bool                | 0           |
|                                  | turb. diffusion    |                     |             |
|                                  | terms in QKE eqn.  |                     |             |
+----------------------------------+--------------------+---------------------+-------------+

Note that the MYNN2.5 scheme must be used in conjunction with a MOST boundary condition
at the surface (Zlo) boundary.

If the PBL scheme is activated, it determines the turbulent diffusivity in the vertical
direction. If an LES model is also specified, it determines only the horizontal turbulent
diffusivity.

Right now, the QKE equation is solved if and only if the MYNN2.5 PBL model is selected. In that
transport equation, it is optional to advect QKE, and to apply LES diffusive transport for QKE
in the horizontal directions (the vertical component is always computed as part of the PBL
scheme).

Forcing Terms
=============

.. _list-of-parameters-14:

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
| **erf.use_rayleigh_damping**     | Include explicit  | true / false      | false       |
|                                  | Rayleigh damping  |                   |             |
+----------------------------------+-------------------+-------------------+-------------+

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

Initialization
==============

ERF can be initialized in different ways. These are listed below:

- Custom initialization:
    Several problems under **Exec** are initialized in a custom manner. The state and velocity components are specific to the problem. These problems are meant for demonstration and do not include any terrain or map scale factors.
- Initialization using a NetCDF file:
    Problems in ERF can be initialized using a NetCDF file containing the mesoscale data. The state and velocity components of the ERF domain are ingested from the mesoscale data. This is a more realistic problem with real atmospheric data used for initialization. The typical filename used for initialization is ``wrfinput_d01``, which is the outcome of running ``ideal.exe`` or ``real.exe`` of the WPS/WRF system.  These problems are run with both terrain and map scale factors.
- Initialization using an ``input_sounding`` file:
    Problems in ERF can be initialized using an ``input_sounding`` file containing the vertical profile. This file has the same format as used by ``ideal.exe`` executable in WRF. Using this option for initialization, running ``ideal.exe`` and reading from the resulting ``wrfinput_d01`` file are not needed. This option is used for initializing ERF domain to a horizontally homogeneous mesoscale state and does not include terrain or map scale factors.

In addition, there is a run-time option to project the initial velocity field to make it divergence-free.  To take
advantage of this option, the code must be built with ``USE_POISSON_SOLVE = TRUE`` in the GNUmakefile if using gmake, or with
``-DERF_ENABLE_POISSON_SOLVE:BOOL=ON`` in the cmake.sh file if using cmake.

List of Parameters
------------------

+----------------------------------+-------------------+--------------------+------------------+
| Parameter                        | Definition        | Acceptable         | Default          |
|                                  |                   | Values             |                  |
+==================================+===================+====================+==================+
| **erf.init_type**                | Initialization    | “custom”,          | “*custom*”       |
|                                  | type              | “ideal”,           |                  |
|                                  |                   | "real",            |                  |
|                                  |                   |"input_sounding"    |                  |
+----------------------------------+-------------------+--------------------+------------------+
| **erf.input_sounding_file**      | Path to WRF-style |  String            | "input_sounding" |
|                                  | input sounding    |                    |                  |
|                                  | file              |                    |                  |
+----------------------------------+-------------------+--------------------+------------------+
| **erf.init_sounding_ideal**      | Perform           |  true or false     | false            |
|                                  | initialization    |                    |                  |
|                                  | like WRF's        |                    |                  |
|                                  | ideal.exe         |                    |                  |
+----------------------------------+-------------------+--------------------+------------------+
| **erf.use_real_bcs**             | If init_type is   | true or false      | true if          |
|                                  | real or metgrid,  |                    | if init_type     |
|                                  | do we want to use |                    | is real or       |
|                                  | these bcs?        |                    | metgrid;         |
|                                  |                   |                    | else false       |
+----------------------------------+-------------------+--------------------+------------------+
| **erf.nc_init_file**             | NetCDF file with  |  String            | NONE             |
|                                  | initial mesoscale |                    |                  |
|                                  | data              |                    |                  |
+----------------------------------+-------------------+--------------------+------------------+
| **erf.nc_bdy_file**              | NetCDF file with  |  String            | NONE             |
|                                  | mesoscale data at |                    |                  |
|                                  | lateral boundaries|                    |                  |
+----------------------------------+-------------------+--------------------+------------------+
| **erf.project_initial_velocity** | project initial   |  Integer           | 1                |
|                                  | velocity?         |                    |                  |
+----------------------------------+-------------------+--------------------+------------------+

Notes
-----------------

If **erf.init_type = ideal**, the problem is initialized with mesoscale data contained in a NetCDF file, provided via ``erf.nc_init_file``. The mesoscale data are horizontally homogeneous, i.e., there is variation only in the vertical direction.

If **erf.init_type = real**, the problem is initialized with mesoscale data contained in a NetCDF file,
provided via ``erf.nc_init_file``. The mesoscale data are realistic with variation in all three directions.
In addition, the lateral boundary conditions must be supplied in a NetCDF files specified by **erf.nc_bdy_file = wrfbdy_d01**

If **erf.init_type = input_sounding**, a WRF-style input sounding is read from
``erf.input_sounding_file``. This text file includes any set of levels that
goes at least up to the model top height. The first line includes the surface
pressure [hPa], potential temperature [K], and water vapor mixing ratio [g/kg].
Each subsequent line has five input values: height [m above sea level], dry
potential temperature [K], vapor mixing ratio [g/kg], x-direction wind
component [m/s], and y-direction wind component [m/s]. Please pay attention to
the units of pressure and mixing ratio. If **erf.init_sounding_ideal = true**,
then moist and dry conditions throughout the air column are determined by
integrating the hydrostatic equation from the surface.

If **erf.init_type = custom** or **erf.init_type = input_sounding**, ``erf.nc_init_file`` and ``erf.nc_bdy_file`` do not need to be set.

Setting **erf.project_initial_velocity = 1** will have no effect if the code is not built with **ERF_USE_POISSON_SOLVE** defined.

Map Scale Factors
=================

Map scale factors are always present in the evolution equations, but the values default to one
unless specified in the initialization when **erf.init_type = real**.

There is an option to test the map scale factors by setting  **erf.test_mapfactor = true**; this
arbitrarily sets the map factors to 0.5 in order to test the implementation.

Terrain
=======

ERF allows the use to specify whether terrain-fitted coordinates should be used by
setting **erf.use_terrain** (default false).
If terrain-fitted coordinates are chosen, they are defined to be static (default)
or moving by setting **erf.terrain_type.**
If using terrain, the user also has the option to specify one of three
methods for defining how the terrain-fitted coordinates given the topography:

- Basic Terrain Following (BTF):
    The influence of the terrain decreases linearly with height.
- Smoothed Terrain Following (STF):
    Small-scale terrain structures are progressively smoothed out of the coordinate system as height increases.
- Sullivan Terrain Following (name TBD):
    The influence of the terrain decreases with the cube of height.

List of Parameters
------------------

+-----------------------------+--------------------+--------------------+------------+
| Parameter                   | Definition         | Acceptable         | Default    |
|                             |                    | Values             |            |
+=============================+====================+====================+============+
| **erf.use_terrain**         | use terrain-fitted |  true / false      | false      |
|                             | coordinates?       |                    |            |
+-----------------------------+--------------------+--------------------+------------+
| **erf.terrain_type**        | static or moving?  |  Static / Moving   | Static     |
+-----------------------------+--------------------+--------------------+------------+
| **erf.terrain_smoothing**   | specify terrain    | 0,                 | 0          |
|                             | following          | 1,                 |            |
|                             |                    | 2                  |            |
+-----------------------------+--------------------+--------------------+------------+


Examples of Usage
-----------------

-  **erf.terrain_smoothing**  = 0
    BTF is used when generating the terrain following coordinate.

-  **erf.terrain_smoothing**  = 1
    STF is used when generating the terrain following coordinate.

-  **erf.terrain_smoothing**  = 2
    Sullivan TF is used when generating the terrain following coordinate.

Moisture
========

ERF has several different moisture models.
The following run-time options control how the full moisture model is used.

List of Parameters
------------------

+-----------------------------+--------------------------+--------------------+------------+
| Parameter                   | Definition               | Acceptable         | Default    |
|                             |                          | Values             |            |
+=============================+==========================+====================+============+
| **erf.moisture_model**      | Name of moisture model   |  "SAM", "Kessler", | "Null"     |
|                             |                          |  "FastEddy"        |            |
+-----------------------------+--------------------------+--------------------+------------+
| **erf.do_cloud**            | use basic moisture model |  true / false      | true       |
+-----------------------------+--------------------------+--------------------+------------+
| **erf.do_precip**           | include precipitation    |  true / false      | true       |
|                             | in treatment of moisture |                    |            |
+-----------------------------+--------------------------+--------------------+------------+

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
| **amrex.fpe_trap_invalid**  | Raise errors for NaNs     |  0 / 1            | 0          |
+-----------------------------+---------------------------+-------------------+------------+
| **amrex.fpe_trap_zero**     | Raise errors for divide   |  0 / 1            | 0          |
|                             | by zero                   |                   |            |
+-----------------------------+---------------------------+-------------------+------------+
| **amrex.fpe_trap_overflow** | Raise errors for overflow |  0 / 1            | 0          |
+-----------------------------+---------------------------+-------------------+------------+
