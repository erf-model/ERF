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
| **erf.use_fft**          | use FFT rather than         | Boolean       | false       |
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
| **geometry.prob_lo**     | physical        | Real            | 0 0 0       |
|                          | location of low |                 |             |
|                          | corner of the   |                 |             |
|                          | domain          |                 |             |
+--------------------------+-----------------+-----------------+-------------+
| **geometry.prob_hi**     | physical        | Real            | must be set |
|                          | location of     |                 | unless      |
|                          | high corner of  |                 | prob_extent |
|                          | the domain      |                 | is set      |
+--------------------------+-----------------+-----------------+-------------+
| **geometry.prob_extent** | physical size   | Real            | must be set |
|                          | of the domain   |                 | unless      |
|                          | in each         |                 | prob_hi     |
|                          | direction       |                 | is set      |
+--------------------------+-----------------+-----------------+-------------+
| **geometry.is_periodic** | is the domain   | 0 if false, 1   | 0 0 0       |
|                          | periodic in     | if true         |             |
|                          | this direction  |                 |             |
+--------------------------+-----------------+-----------------+-------------+

.. note::

   The extent of the domain may be specified either by **geometry.prob_hi** or by
   **geometry.prob_extent**, where ``prob_hi = prob_lo + prob_extent``.  Exactly one
   of the two must be given; specifying both is an error.  Since **geometry.prob_lo**
   defaults to 0 in every direction, a domain whose low corner is at the origin can be
   defined by setting **geometry.prob_extent** alone, i.e. without setting either
   **geometry.prob_lo** or **geometry.prob_hi**.

   These parameters may also be given with the ``erf`` prefix, i.e. as **erf.prob_lo**,
   **erf.prob_hi**, **erf.prob_extent** and **erf.is_periodic**, but the ``erf``- and
   ``geometry``-prefixed forms may not be mixed.  Note that if the ``erf`` prefix is used,
   then **erf.prob_lo** must be set explicitly along with **erf.prob_hi** or
   **erf.prob_extent**.

Examples of Usage
-----------------

-  **geometry.prob_lo** = 0 0 0
   defines the low corner of the domain at (0,0,0) in physical space.

-  **geometry.prob_hi** = 1.e8 2.e8 2.e8
   defines the high corner of the domain at (1.e8,2.e8,2.e8) in
   physical space.

-  **geometry.prob_extent** = 1.e8 2.e8 2.e8
   defines a domain of size (1.e8,2.e8,2.e8) in physical space; combined with
   the default **geometry.prob_lo** = 0 0 0, this is equivalent to the two
   examples above.

-  **geometry.is_periodic** = 0 1 0
   says the domain is periodic in the y-direction only.

Domain Boundary Conditions
==========================
Domain boundary conditions in ERF may be broadly categorized as ``ideal`` or ``real`` where
the ideal BC types correspond those found in classic fluid solvers and real correspond to an
external data source that may be based upon observation data. The inputs for these types of BCs
are given in more detail in :ref:`sec:LateralBoundaryConditions`. We briefly note that the supported
ideal BC types are ``inflow``, ``inflow_outflow``, ``outflow``, ``ho_outflow``, ``open``,
``slipwall``, ``noslipwall``, ``symmetry`` and ``surface_layer``. The type strings are matched
case-insensitively, so e.g. ``NoSlipWall`` and ``noslipwall`` are equivalent.


.. _list-of-parameters-1:

List of Parameters
------------------

+-----------------+-----------------------------------+---------------------------------------------------------------------------+--------------------------------+
| Parameter       | Definition                        | Acceptable Values                                                         | Default                        |
+=================+===================================+===========================================================================+================================+
| **xlo.type**    | boundary type of the xlo face     | inflow, outflow, inflow_outflow, open, ho_outflow, slipwall, noslipwall,  | must be set if the direction   |
|                 |                                   | symmetry, surface_layer                                                   | is not periodic                |
+-----------------+-----------------------------------+---------------------------------------------------------------------------+--------------------------------+
| **xhi.type**    | boundary type of the xhi face     | inflow, outflow, inflow_outflow, open, ho_outflow, slipwall, noslipwall,  | must be set if the direction   |
|                 |                                   | symmetry, surface_layer                                                   | is not periodic                |
+-----------------+-----------------------------------+---------------------------------------------------------------------------+--------------------------------+
| **ylo.type**    | boundary type of the ylo face     | inflow, outflow, inflow_outflow, open, ho_outflow, slipwall, noslipwall,  | must be set if the direction   |
|                 |                                   | symmetry, surface_layer                                                   | is not periodic                |
+-----------------+-----------------------------------+---------------------------------------------------------------------------+--------------------------------+
| **yhi.type**    | boundary type of the yhi face     | inflow, outflow, inflow_outflow, open, ho_outflow, slipwall, noslipwall,  | must be set if the direction   |
|                 |                                   | symmetry, surface_layer                                                   | is not periodic                |
+-----------------+-----------------------------------+---------------------------------------------------------------------------+--------------------------------+
| **zlo.type**    | boundary type of the zlo face     | inflow, outflow, inflow_outflow, open, ho_outflow, slipwall, noslipwall,  | must be set if the direction   |
|                 |                                   | symmetry, surface_layer                                                   | is not periodic                |
+-----------------+-----------------------------------+---------------------------------------------------------------------------+--------------------------------+
| **zhi.type**    | boundary type of the zhi face     | inflow, outflow, inflow_outflow, open, ho_outflow, slipwall, noslipwall,  | must be set if the direction   |
|                 |                                   | symmetry, surface_layer                                                   | is not periodic                |
+-----------------+-----------------------------------+---------------------------------------------------------------------------+--------------------------------+


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
|                                   | subsequent levels, |                 |             |
|                                   | in every direction |                 |             |
+-----------------------------------+--------------------+-----------------+-------------+
| **amr.ref_ratio_vect**            | ratio of coarse    | 3 integers >= 1 | 2 for all   |
|                                   | to fine grid       | per level, one  | levels and  |
|                                   | spacing between    | per coordinate  | directions  |
|                                   | subsequent levels, | direction       |             |
|                                   | per direction      |                 |             |
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
| **amr.regrid_level_0_on_restart** | should we          | Boolean         | false       |
|                                   | regrid level       |                 |             |
|                                   | immediately        |                 |             |
|                                   | after              |                 |             |
|                                   | restarting         |                 |             |
+-----------------------------------+--------------------+-----------------+-------------+
| **amr.iterate_grids**             | do we iterate      | Boolean         | true        |
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
   | with **amr.max_level** = 1, this would set factor {2 in x-dir, 4 in y-dir,
     3 in z-dir} refinement between levels 0 and 1.  Note that you must specify
     3 values for every refined level, i.e. 3 :math:`\times` **amr.max_level**
     values in all, ordered (x,y,z) for level 0, then (x,y,z) for level 1, and
     so on.  **amr.ref_ratio** and **amr.ref_ratio_vect** may not both be given.

-  | **amr.ref_ratio_vect** = 2 2 1
   | refines by a factor of 2 in x and y but not at all in z, so the fine level
     has the same vertical grid spacing as the coarse level.  Note that
     **amr.ref_ratio** cannot express this, since it applies a single ratio in
     all three directions; see :ref:`MeshRefinement` for the distinction between
     refining in z and creating a finer level over a region.

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

.. _notes-1:

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

+------------------------------+----------------------------------+-----------------+---------------------+
| Parameter                    | Definition                       | Acceptable      | Default             |
|                              |                                  | Values          |                     |
+==============================+==================================+=================+=====================+
| **amr.regrid_file**          | name of file from which to read  | text            | no file             |
|                              | the grids                        |                 |                     |
+------------------------------+----------------------------------+-----------------+---------------------+
| **amr.grid_eff**             | grid efficiency at coarse level  | Real > 0, < 1   | 0.7                 |
|                              | at which grids are created       |                 |                     |
+------------------------------+----------------------------------+-----------------+---------------------+
| **amr.n_error_buf**          | radius of additional tagging     | Integer >= 0    | 0                   |
|                              | around already tagged cells      |                 |                     |
+------------------------------+----------------------------------+-----------------+---------------------+
| **amr.n_error_buf_x**        | radius of additional tagging in  | Integer >= 0    | **n_error_buf**     |
|                              | the x-direction                  |                 |                     |
+------------------------------+----------------------------------+-----------------+---------------------+
| **amr.n_error_buf_y**        | radius of additional tagging in  | Integer >= 0    | **n_error_buf**     |
|                              | the y-direction                  |                 |                     |
+------------------------------+----------------------------------+-----------------+---------------------+
| **amr.n_error_buf_z**        | radius of additional tagging in  | Integer >= 0    | **n_error_buf**     |
|                              | the z-direction                  |                 |                     |
+------------------------------+----------------------------------+-----------------+---------------------+
| **amr.n_proper**             | minimum number of coarse cells   | Integer > 0     | 2                   |
|                              | between coarse-fine boundaries   |                 |                     |
+------------------------------+----------------------------------+-----------------+---------------------+
| **amr.max_grid_size**        | maximum size of a grid in every  | Integer > 0     | 2048                |
|                              | direction                        |                 |                     |
+------------------------------+----------------------------------+-----------------+---------------------+
| **amr.max_grid_size_x**      | maximum size of a grid in the    | Integer > 0     | **max_grid_size**   |
|                              | x-direction                      |                 |                     |
+------------------------------+----------------------------------+-----------------+---------------------+
| **amr.max_grid_size_y**      | maximum size of a grid in the    | Integer > 0     | **max_grid_size**   |
|                              | y-direction                      |                 |                     |
+------------------------------+----------------------------------+-----------------+---------------------+
| **amr.max_grid_size_z**      | maximum size of a grid in the    | Integer > 0     | **max_grid_size**   |
|                              | z-direction                      |                 |                     |
+------------------------------+----------------------------------+-----------------+---------------------+
| **amr.blocking_factor**      | grid size must be a multiple of  | Integer > 0     | 1                   |
|                              | this in every direction          |                 |                     |
+------------------------------+----------------------------------+-----------------+---------------------+
| **amr.blocking_factor_x**    | grid size in x must be a         | Integer > 0     | **blocking_factor** |
|                              | multiple of this                 |                 |                     |
+------------------------------+----------------------------------+-----------------+---------------------+
| **amr.blocking_factor_y**    | grid size in y must be a         | Integer > 0     | **blocking_factor** |
|                              | multiple of this                 |                 |                     |
+------------------------------+----------------------------------+-----------------+---------------------+
| **amr.blocking_factor_z**    | grid size in z must be a         | Integer > 0     | **blocking_factor** |
|                              | multiple of this                 |                 |                     |
+------------------------------+----------------------------------+-----------------+---------------------+
| **amr.refine_grid_layout**   | chop grids further if # of       | 0 if false, 1   | 1                   |
|                              | processors > # of grids          | if true         |                     |
+------------------------------+----------------------------------+-----------------+---------------------+
| **amr.refine_grid_layout_x** | chop in x when refining the grid | 0 if false, 1   | 1                   |
|                              | layout                           | if true         |                     |
+------------------------------+----------------------------------+-----------------+---------------------+
| **amr.refine_grid_layout_y** | chop in y when refining the grid | 0 if false, 1   | 1                   |
|                              | layout                           | if true         |                     |
+------------------------------+----------------------------------+-----------------+---------------------+
| **amr.refine_grid_layout_z** | chop in z when refining the grid | 0 if false, 1   | 0                   |
|                              | layout                           | if true         |                     |
+------------------------------+----------------------------------+-----------------+---------------------+

.. _notes-2:

Notes
-----

-  ERF changes several of the AMReX gridding defaults; the values in the table
   above are the ERF defaults, and they are set in ``add_par`` in
   ``Source/main.cpp``.  In particular **amr.max_grid_size** defaults to a very
   large value (so that grids are chopped only when there are more processors
   than grids), **amr.blocking_factor** defaults to 1, and
   **amr.refine_grid_layout_z** defaults to 0 (the AMReX default is 1).

-  **amr.n_error_buf**, **amr.max_grid_size** and
   **amr.blocking_factor** can be read in as a single value which is
   assigned to every level, or as multiple values, one for each level

-  **amr.n_error_buf**, **amr.max_grid_size** and **amr.blocking_factor** apply
   to all coordinate directions; the per-direction forms
   **amr.n_error_buf_x/_y/_z**, **amr.max_grid_size_x/_y/_z** and
   **amr.blocking_factor_x/_y/_z** override them in the specified direction
   only, and may also be specified per level

-  **amr.max_grid_size** at every level must be even

-  **amr.blocking_factor** at every level must be a power of 2

-  the domain size **amr.n_cell** must be a multiple of
   **amr.blocking_factor** at level 0

-  **amr.max_grid_size** must be a multiple of **amr.blocking_factor**
   at every level

-  **amr.refine_grid_layout** is a single flag that sets all three of
   **amr.refine_grid_layout_x/_y/_z**; if a per-direction flag is also
   specified it takes precedence in that direction

.. _examples-of-usage-4:

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

-  | **amr.max_grid_size_x** = 64
   | **amr.max_grid_size_y** = 64
   | The final grids will be no longer than 64 cells in the x and y
     directions, but will not be limited in the z direction.  This is the
     way to limit the horizontal box size without decomposing the grids
     vertically.

-  | **amr.refine_grid_layout_z** = 1
   | Allow the grids to be chopped in the vertical direction as well as
     horizontally when there are more processors than grids.  This is *not*
     the ERF default.

.. _subsec:no-vertical-decomposition:

Avoiding Decomposition in the Vertical Direction
------------------------------------------------

Many ERF workflows want each box to span the full vertical extent of the domain
or of the refined region -- for example because a physics package operates on
entire columns.  (Native SHOC requires this and will abort if any box is split
in z.)  ERF is set up so that this is the default behavior, in both places where
grids are created:

-  **When the level 0 grids are created**, ERF decomposes the domain across the
   processors itself (see ``ERFPostProcessBaseGrids``).  It decomposes in the
   vertical direction only if **amr.max_grid_size_z** is smaller than
   **amr.n_cell** in the z direction.  Since **amr.max_grid_size** defaults to a
   very large value, no vertical decomposition is done unless it is requested.

-  **When grids are chopped to give at least one grid per processor**, either
   during regridding at a finer level or when regridding level 0 on restart,
   AMReX only chops in the directions for which the corresponding
   **amr.refine_grid_layout_x/_y/_z** flag is set.  Because ERF defaults
   **amr.refine_grid_layout_z** to 0, this load-balancing step never splits a
   box in the vertical direction.

The usual way to *accidentally* introduce a vertical decomposition is to set
**amr.max_grid_size** as a single value, since that limits the box size in all
three directions.  To limit the box size horizontally only, use the
per-direction forms, e.g.

::

     amr.max_grid_size_x = 64
     amr.max_grid_size_y = 64

and leave **amr.max_grid_size_z** at its (large) default.  The same holds for
**amr.blocking_factor** versus **amr.blocking_factor_x/_y**.

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

.. _examples-of-usage-5:

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
We note that when using implicit substepping, the vertical grid spacing does not appear in the time step calculation.
The number of acoustic sub-steps per timestep can also be specified by the user as a fixed value or by specifying the
number of substeps per RK stage.  For the anelastic equations, the timestep calculation uses the advective CFL constraint,
which means it is determined by the fluid speed rather than the sound speed and thus allows much larger timesteps.

.. _list-of-parameters-6:

List of Parameters
------------------

+-------------------------------------+----------------------+----------------+---------------------+
| Parameter                           | Definition           | Acceptable     | Default             |
|                                     |                      | Values         |                     |
+=====================================+======================+================+=====================+
| **erf.substepping_type**            | Should we substep in | "Implicit" or  | "Implicit" if       |
|                                     | each RK stage?       | "None"         | compressible,       |
|                                     |                      |                | "None" if anelastic |
+-------------------------------------+----------------------+----------------+---------------------+
| **erf.vert_implicit**               | Do vertical implicit | Boolean        | true if compressible|
|                                     | solve for diffusion  |                | false if anelastic  |
|                                     | of u, v, theta, KE,  |                |                     |
|                                     | and qv with default  |                |                     |
|                                     | time-centering in    |                |                     |
|                                     | each stage           |                |                     |
+-------------------------------------+----------------------+----------------+---------------------+
| **erf.vert_implicit_fac**           | How much implicit    | Real >= 0      | 1.0, 1.0, 0.0       |
|                                     | vertical diffusion   | (explicit) and |                     |
|                                     | to include in each   | <= 1 (implicit)|                     |
|                                     | RK stage?            |                |                     |
|                                     |                      |                |                     |
|                                     | Specify either one   |                |                     |
|                                     | (the same for all    |                |                     |
|                                     | stages) or three     |                |                     |
|                                     | values, one per stage|                |                     |
+-------------------------------------+----------------------+----------------+---------------------+
| **erf.implicit_thermal_diffusion**  | Do vertical implicit | Boolean        | true                |
|                                     | for theta?           |                |                     |
+-------------------------------------+----------------------+----------------+---------------------+
| **erf.implicit_moisture_diffusion** | Do vertical implicit | Boolean        | true                |
|                                     | for Qv?              |                |                     |
+-------------------------------------+----------------------+----------------+---------------------+
| **erf.implicit_ke_diffusion**       | Do vertical implicit | Boolean        | true                |
|                                     | for KE?              |                |                     |
+-------------------------------------+----------------------+----------------+---------------------+
| **erf.implicit_momentum_diffusion** | Do vertical implicit | Boolean        | true                |
|                                     | for U & V?           |                |                     |
+-------------------------------------+----------------------+----------------+---------------------+
| **erf.cfl**                         | CFL number used to   | Real > 0 and   | 0.8                 |
|                                     | compute level 0 dt   | <= 1           |                     |
+-------------------------------------+----------------------+----------------+---------------------+
| **erf.substepping_cfl**             | CFL number used to   | Real > 0 and   | 1.0                 |
|                                     | compute the number   | <= 1           |                     |
|                                     | of substeps          |                |                     |
+-------------------------------------+----------------------+----------------+---------------------+
| **erf.fixed_dt**                    | set level 0 dt       | Real > 0       | unused if not       |
|                                     | as this value        |                | set                 |
|                                     | regardless of        |                |                     |
|                                     | cfl or other         |                |                     |
|                                     | settings             |                |                     |
+-------------------------------------+----------------------+----------------+---------------------+
| **erf.fixed_fast_dt**               | set fast dt          | Real > 0       | unused if not       |
|                                     | as this value        |                | set                 |
+-------------------------------------+----------------------+----------------+---------------------+
| **erf.fixed_mri_dt_ratio**          | set fast dt          | even int > 0   | only relevant if    |
|                                     | as slow dt /         |                | substepping_type    |
|                                     | this ratio           |                | is not None         |
+-------------------------------------+----------------------+----------------+---------------------+
| **erf.init_shrink**                 | factor by which      | Real > 0 and   | 1.0                 |
|                                     | to shrink the        | <= 1           |                     |
|                                     | initial dt           |                |                     |
+-------------------------------------+----------------------+----------------+---------------------+
| **erf.change_max**                  | factor by which      | Real >= 1      | 1.1                 |
|                                     | dt can grow          |                |                     |
|                                     | in subsequent        |                |                     |
|                                     | steps                |                |                     |
+-------------------------------------+----------------------+----------------+---------------------+
| **erf.dt_max**                      | maximum adaptive     | Real > 0       | 1e9                 |
|                                     | timestep             |                |                     |
|                                     | allowed by time      |                |                     |
|                                     | stepping             |                |                     |
+-------------------------------------+----------------------+----------------+---------------------+
| **erf.dt_max_initial**              | maximum initial      | Real > 0       | 1.0                 |
|                                     | timestep             |                |                     |
+-------------------------------------+----------------------+----------------+---------------------+
| **erf.dt_ref_ratio**                | ratio of coarse      | Integer >= 1   | same as             |
|                                     | to fine grid         | (one per level)| maximum over        |
|                                     | time steps between   |                | directions of       |
|                                     | subsequent           |                | ref_ratio           |
|                                     | levels               |                |                     |
+-------------------------------------+----------------------+----------------+---------------------+

Notes
-----------------

-  | If **erf.anelastic** is true then **substepping_type** is internally set to "None".

-  | The time step controls work somewhat differently depending on whether one is using
     acoustic substepping in time; this is determined by the value of **substepping_type**.

-  | If **erf.substepping_type = None** there is only one time step to be calculated,
     and **erf.fixed_fast_dt** and **erf.fixed_mri_dt_ratio** are not used.

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

     * | If **erf.fixed_dt** and **erf.fixed_fast_dt** are specified and the ratio of **fixed_dt** to **fixed_fast_dt**
         is not an even positive integer

     * | If **erf.fixed_dt** and **erf.fixed_fast_dt** and **erf.fixed_mri_dt_ratio** are all specified but are inconsistent

   * | Once the slow timestep is set and the inputs are allowed per the above criteria,
       the fast timestep is computed in one of several ways:

     * | If **erf.fixed_fast_dt** is specified, the fast timestep will be set to **fixed_fast_dt**.

     * | If **erf.fixed_mri_dt_ratio** is specified and **erf.fixed_fast_dt** is not specified,
         the fast timestep will be the slow timestep divided by **fixed_mri_dt_ratio.**

     * | If neither **erf.fixed_mri_dt_ratio** nor **erf.fixed_fast_dt** is specified, then the fast timestep
         will be computed using the CFL condition for compressible flow, then adjusted (reduced if necessary)
         as above so that the ratio of slow timestep to fine timestep is an even integer.
         If **erf.substepping_cfl** is specified, that CFL value will be used.  If not, the default value will be used.

.. _examples-of-usage-6:

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

Additional Checkpoint Controls
------------------------------

These are AMReX-level checkpoint controls that may be specified either with the
``erf`` prefix or the ``amr`` prefix.

List of Parameters
~~~~~~~~~~~~~~~~~~

+--------------------+------------------------------+---------------------+---------+
| Parameter          | Definition                   | Acceptable Values   | Default |
+====================+==============================+=====================+=========+
| **amr.check_int**  | Checkpoint interval (steps)  | Integer >= 0        | -1      |
+--------------------+------------------------------+---------------------+---------+
| **amr.check_per**  | Checkpoint interval (time)   | Real >= 0           | -1.0    |
+--------------------+------------------------------+---------------------+---------+

Notes
~~~~~

- If both are specified, **amr.check_int** takes precedence and **amr.check_per**
  is ignored.
- These are in addition to the main checkpoint controls described in
  :ref:`sec:Checkpoint`.

PlotFiles
===============================

See :ref:`sec:Plotfiles` for how to control the types and frequency of plotfile
generation.

Boundary Files
==============

+----------------------+------------------------------+-------------------+----------------------+
| Parameter            | Definition                   | Acceptable Values | Default              |
+======================+==============================+===================+======================+
| **erf.write_erfbdy** | Write AMReX-native format    | Boolean           | true for non-restart |
|                      | boundary file for real-data  |                   | real data cases,     |
|                      | cases only                   |                   | otherwise false      |
+----------------------+------------------------------+-------------------+----------------------+
| **erf.erfbdy_file**  | Name of the boundary file    | String            | "erfbdy"             |
+----------------------+------------------------------+-------------------+----------------------+

Additional Plotfile Controls
----------------------------

List of Parameters
~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 32 40 18 10

   * - Parameter
     - Definition
     - Acceptable values
     - Default
   * - **erf.expand_plotvars_to_unif_rr**
     - Expand plot variables to a uniform refinement ratio for mixed-refinement cases
     - Boolean
     - false
   * - **erf.compute_mean_vars**
     - Accumulate interval first and second moments for resolved plot diagnostics
     - Boolean
     - false
   * - **erf.mean_vars_reset_mode**
     - Reset moments after a complete 3-D plotfile batch or once at a specified time
     - ``plotfile`` or ``time``
     - ``plotfile``
   * - **erf.mean_vars_reset_time**
     - Reset time when reset mode is ``time``
     - Real [s], >= 0
     - -1.0

When ``erf.compute_mean_vars = true``, ERF makes the interval fields
``u_mean``, ``v_mean``, ``w_mean``, ``theta_mean``, the corresponding second
moments, resolved variances/covariances, and resolved ``tke_resolved`` available to
``erf.plot_vars_1`` and ``erf.plot_vars_2``. See
:ref:`sec:Plotfile3DReference` for definitions. In ``time`` reset mode, the
accumulator is reset once, at the first step whose start time is at or beyond
``erf.mean_vars_reset_time``.
In ``plotfile`` mode, a reset occurs only after all 3-D streams due at the
same output event have been written and at least one emitted stream contains
an interval diagnostic; a density-only stream does not close the shared
averaging window.
In either mode, regridding a level rebuilds that level's accumulator, so the
averaging window restarts on the regridded level while the other levels keep
accumulating. A plotfile written shortly after a regrid can therefore mix
per-level averaging windows; see :ref:`sec:Plotfile3DReference`.


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
These are **useful for idealized, horizontally homogeneous simulation domains**
and include (1) a surface time history file, (2) a history of mean profiles,
(3) a history of vertical flux profiles (i.e., variances and covariances),
(4) a history of modeled subgrid stresses, and (5) a history of large scale
forcing and nudging data. The profiles are calculated from planar averages.

The number of output filenames specified through ``erf.data_log`` dictates the
level of output. E.g., specifying 3 filenames will give outputs (1), (2), and (3);
specifying 5 filenames will also include output (5). Data files are only written
if ``erf.profile_int > 0``.

This output functionality has not been implemented for terrain.
For **real simulation domains**, users should use 2-D and 3-D :ref:`sec:Plotfiles`.

.. _list-of-parameters-10a:


List of Parameters
------------------

+-------------------------------+------------------+----------------+----------------+
| Parameter                     | Definition       | Acceptable     | Default        |
|                               |                  | Values         |                |
+===============================+==================+================+================+
| **erf.data_log**              | Output           | Up to five     | NONE           |
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

All profiles are planar-averaged quantities :math:`\langle\cdot\rangle`
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

* Large scale forcing and nudging profiles

  #. Time (s)

  #. Height (m)

  #. Total applied temperature tendency, :math:`\partial \theta / \partial t` (K/s)

  #. Total applied water vapor tendency, :math:`\partial q_v / \partial t` (kg/kg/s)

  #. Large-scale subsidence velocity, :math:`w_{sub}` (m/s)

  #. Horizontal temperature tendency, :math:`\partial \theta / \partial t` (K/s)

  #. Horizontal water vapor tendency, :math:`\partial q_v / \partial t` (kg/kg/s)

  #. Vertical temperature tendency, :math:`\partial \theta / \partial t` (K/s)

  #. Vertical water vapor tendency, :math:`\partial q_v / \partial t` (kg/kg/s)

  #. Vertical cloud water tendency, :math:`\partial q_c / \partial t` (kg/kg/s)

  #. Temperature nudging tendency, :math:`\theta_{nudge}` (K/s)

  #. Water vapor nudging tendency, :math:`q_{v,nudge}` (kg/kg/s)

  #. X-velocity nudging tendency, :math:`u_{nudge}` (m/s2)

  #. Y-velocity nudging tendency, :math:`v_{nudge}` (m/s2)

.. note::

   The vertical tendency columns are only populated for interior cells between
   the bottom and top boundaries; they are zero at the domain edges.


Data Sampling Outputs
=====================

.. note::
   **For WRF users**: The ERF analog of **wrfout** and **auxhist** history
   files are AMReX-format :ref:`sec:Plotfiles`. Unlike WRF, which combines
   2-D and 3-D field quantities into a single history file, ERF allows separate
   outputs at different intervals.

   To obtain history files in netcdf format, it is recommended to convert from
   the native AMReX output using postprocessing tools provided in Exec/Tools if
   using gmake, or with the ``ERF_ENABLE_TOOLS`` flag if using cmake.

   The ERF analog of **tslist** output is the line sampling described in this
   section.

Data along query lines or planes may be output during the simulation if
``erf.do_line_sampling = true`` or  ``erf.do_plane_sampling = true``, respectively.
The potential temperature and wind-speed will be written to native ``plt_line/plane``
at the step frequency dictated by ``erf.line_sampling_interval = <int>`` or
``erf.plane_sampling_interval = <int>``. For line sampling,
users must prescribe either ``sample_line_lo`` and ``sample_line_hi`` inputs, which are
3 integer values corresponding to the (i,j,k) indices at the beginning and end of the line,
or ``sample_line_lo_real`` and ``sample_line_hi_real``, which are 3 real values corresponding
to the physical locations of the beginning and end of the line. These two ways of specifying
the line are mutually exclusive. Additionally, users must specify ``sample_line_dir`` to
prescribed the direction of the line. The same inputs are used for the plane sampling except
that ``sample_plane_lo/hi``
must be the physical locations of the plane corners. This output functionality has
not been implemented for terrain. By default, sampled line and plane data will have the
prefixes "plt_line" and "plt_plane", respectively. Names for sampled data may optionally
be provided with ``sample_line_name`` and/or ``sample_plane_name`` -- if provided, each
line and/or plane must be named.

Plane samples are written as true multi-level AMReX plotfiles: every refinement level
whose grids intersect the plane is written (``Level_0``, ``Level_1``, ...), so a plane
that cuts through static refinement patches retains the finer in-plane resolution there.
By default all intersecting levels are written; ``erf.plane_sampling_max_level = <int>``
caps the finest level (``0`` forces level-0-only output). The slice-normal direction is
resolved natively on each level by replicating the sampled plane across the level's cells,
so the resulting dataset has an isotropic refinement ratio and loads cleanly in yt/amrvis.

Line and plane samples will be default be written to plotfiles, one plotfile per output
snapshot, with all output variables in the same file. Alternatively, line sampling has
the ``erf.line_sampling_text_output`` option, which writes one text file per output
variable, with all snapshots appended to the same file over time. This is similar to the
tslist output from WRF but output is provided only from the finest domain that contains
the entire requested sampling line; velocities are also destaggered.

The sampled variables can be selected with the ``erf.line_sampling_vars`` and
``erf.plane_sampling_vars`` options and include a subset of the plotfile outputs:
``density``, ``x_velocity``, ``y_velocity``, ``z_velocity``, ``magvel``,
``theta``, ``qv``, ``qc``, and ``pressure``. Line sampling additionally supports
``sgs_tke``, ``sgs_tau13``, ``sgs_tau23``, and ``sgs_hfx3``. The SGS stress and
heat-flux values are averaged from their native edge/face locations to cell
centers; if their diffusion storage is unavailable, ERF writes zero. Velocities
are output at cell centers only. The water vapor mixing ratio ``qv`` will only
output valid values if a moisture model is used. Pressure is calculated from
rho*theta and will account for moisture if ``qv`` is requested.
``sgs_tke`` is zero unless the selected level has a prognostic TKE closure
(``Deardorff``, k-equation RANS, MYJ/MYNN, or a SHOC-family PBL); Smagorinsky
and ``les_type = None`` do not provide prognostic TKE and their unused state
slot is never sampled.

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
| **erf.sample_line_lo_real/**      | Bounding point   | 3 Reals per    | None           |
| **hi_real**                       | on the line(s)   | line           |                |
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
|                                   | variables; SGS   |                |                |
|                                   | fields are line- |                |                |
|                                   | only             |                |                |
+-----------------------------------+------------------+----------------+----------------+
| **erf.plane_sampling_vars**       | Specify sampled  | List of strings| theta, magvel  |
|                                   | variables        |                |                |
+-----------------------------------+------------------+----------------+----------------+
| **erf.plane_sampling_max_level**  | Cap on finest    | Integer        | -1 (all        |
|                                   | level written to |                | intersecting   |
|                                   | the plane        |                | levels)        |
|                                   | plotfile         |                |                |
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
| **erf.use_efficient_advection**  | Use efficient      | Boolean             | false        |
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
|                                         |                    | "YSU", "YSUNew",    |             |
|                                         |                    | "MRF",              |             |
|                                         |                    | "NATIVE_SHOC",      |             |
|                                         |                    | "EAMXX_SHOC";       |             |
|                                         |                    | legacy "SHOC" is    |             |
|                                         |                    | deprecated and maps |             |
|                                         |                    | to "EAMXX_SHOC"     |             |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.tke_min**                         | Minimum initial    | Real                | 1.e-6       |
|                                         | TKE used when a    |                     |             |
|                                         | prognostic-TKE     |                     |             |
|                                         | closure is active  |                     |             |
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
| **erf.pbl_mynn_C2**                     | MYNN Constant C2   | Real                | 0.75        |
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
| **erf.pbl_mynn_diffuse_moistvars**      | Diffuse moisture   | Boolean             | false       |
|                                         | variables using    |                     |             |
|                                         | modeled eddy       |                     |             |
|                                         | diffusivity        |                     |             |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.advect_QKE**                      | Include advection  | Boolean             | true        |
|                                         | terms in QKE eqn   |                     |             |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.diffuse_QKE_3D**                  | Include horizontal | Boolean             | false       |
|                                         | turb. diffusion    |                     |             |
|                                         | terms in QKE eqn.  |                     |             |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_ysu_force_over_water**        | Treat whole domain | Boolean             | false       |
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
| **erf.pbl_ysu_use_consistent_coriolis** | Ignore above param | Boolean             | false       |
|                                         | and use the value  |                     |             |
|                                         | from ERF coriolis  |                     |             |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_mrf_coriolis_freq**           | Coriolis frq. used | Real                | 1.0e-4      |
|                                         | for MRF PBL Scheme |                     |             |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_mrf_Ribcr**                   | Over land critical | Real                | 0.5         |
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
| **erf.mrf_moistvars**                   | Diffuse moisture   | Boolean             | false       |
|                                         | variables using    |                     |             |
|                                         | modeled eddy       |                     |             |
|                                         | diffusivity        |                     |             |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.enable_mrf_countergradient**      | Enable             | Boolean             | false       |
|                                         | countergradient    |                     |             |
|                                         | correction terms   |                     |             |
|                                         | in MRF PBL scheme  |                     |             |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.pbl_mrf_highres_bounds**          | Enable alternative | Boolean             | false       |
|                                         | high-resolution    |                     |             |
|                                         | grid-dependent     |                     |             |
|                                         | diffusivity bounds |                     |             |
+-----------------------------------------+--------------------+---------------------+-------------+
| **erf.enable_mrf_unbounded_vpert**      | Enable physically  | Boolean             | false       |
|                                         | superior           |                     |             |
|                                         | unbounded VPERT    |                     |             |
|                                         | in MRF PBL scheme  |                     |             |
+-----------------------------------------+--------------------+---------------------+-------------+

Note that both PBL schemes must be used in conjunction with a MOST boundary condition
at the surface (Zlo) boundary. The YSU scheme is work in progress currently.

If the PBL scheme is activated, it determines the turbulent diffusivity in the vertical
direction. If an LES model is also specified, it determines only the horizontal turbulent
diffusivity.

ERF uses prognostic TKE for Deardorff LES, k-equation RANS, MYJ, MYNN2.5,
MYNN-EDMF, and SHOC. When one of these closures is active, ERF initializes
TKE at startup to ``erf.tke_min`` (default ``1.e-6 m^2/s^2``) before any
optional problem-specific TKE perturbations or surface-layer-based TKE
initialization are applied.

Parameters for PBL schemes can either be set with one value that applies across all levels, or set with a number of values
equal to the number of levels, allowing unique values of the parameter to be set for each level.

Forest Canopy
=============

Forest drag can be configured either with the existing discrete-patch text file
or with gridded NetCDF fields. The two modes are mutually exclusive. A complete
file configuration enables forest drag automatically; if ``erf.do_forest_drag``
is specified explicitly, its value must agree with the presence or absence of
that configuration. See :ref:`sec:Forest` for the LAD profiles and the scoped
fixed-leaf-temperature exchange equations.

The gridded mode requires NetCDF support and both ``forest_lai_file`` and
``forest_height_file``. Supply exactly one of ``forest_cd_file`` and the
constant ``forest_cd``. The LAI, height, and optional drag-coefficient files
must share the same horizontal grid. ERF reads variables named ``LAI``,
``height``, and ``cd``; each may be two-dimensional or have a leading time
dimension, in which case the first time record is used. Coordinates must be
projected Cartesian ``x``/``y``; raw ``lon``/``lat`` coordinates are rejected
because ERF does not silently reinterpret degrees as meters. Each horizontal direction must contain
at least two finite, strictly increasing, uniformly spaced coordinates. All
forest fields must match the LAI dimensions, spacing, and horizontal origin;
ERF rejects misregistered fields before interpolation. A coordinate pair is
required; ERF does not assume a unit-spacing grid when metadata is absent.

.. list-table:: Forest canopy parameters
   :header-rows: 1
   :widths: 28 45 17 12

   * - Parameter
     - Definition
     - Acceptable values
     - Default
   * - **erf.do_forest_drag**
     - Optional consistency switch for forest drag
     - Boolean
     - false
   * - **erf.forest_file**
     - Discrete-patch canopy text file
     - String
     - None
   * - **erf.forest_lai_file**
     - Gridded NetCDF file containing ``LAI``
     - String
     - None
   * - **erf.forest_height_file**
     - Gridded NetCDF file containing canopy ``height`` [m]
     - String
     - None
   * - **erf.forest_cd_file**
     - Gridded NetCDF file containing ``cd``
     - String
     - None
   * - **erf.forest_cd**
     - Spatially constant drag coefficient; mutually exclusive with ``forest_cd_file``
     - Real
     - None
   * - **erf.forest_tree_type**
     - Vertical LAD profile: uniform or Lalic--Mihailovic
     - 1 or 2
     - 1
   * - **erf.forest_laimax**
     - Fraction of canopy height at maximum LAD for tree type 2; must satisfy ``0 <= laimax < 1``
     - Real
     - 0.8
   * - **erf.forest_substep**
     - Apply canopy momentum source during fast substeps
     - Boolean
     - false
   * - **erf.forest_biophysics**
     - Master switch for the scoped canopy exchange path
     - Boolean
     - false
   * - **erf.forest_biophysics_heat**
     - Apply fixed-leaf-temperature sensible and latent exchange
     - Boolean
     - false
   * - **erf.forest_leaf_theta_fixed**
     - Prescribed leaf potential temperature [K]; required and positive when heat exchange is enabled
     - Real [K]
     - None

``erf.forest_biophysics_heat = true`` requires
``erf.forest_biophysics = true`` and a positive
``erf.forest_leaf_theta_fixed``.

For example, a gridded canopy with fixed leaf potential temperature can be
configured as follows:

::

   erf.do_forest_drag          = true
   erf.forest_lai_file         = canopy_lai.nc
   erf.forest_height_file      = canopy_height.nc
   erf.forest_cd               = 0.15
   erf.forest_tree_type        = 2
   erf.forest_laimax           = 0.6
   erf.forest_biophysics       = true
   erf.forest_biophysics_heat  = true
   erf.forest_leaf_theta_fixed = 301.0

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
|                                     | erf.abl_driver_type =  |                   |                     |
|                                     | PressureGradient)      |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.abl_geo_wind**                | Geostrophic            | 3 Reals           | (0.,0.,0.)          |
|                                     | forcing term           |                   |                     |
|                                     | (only if               |                   |                     |
|                                     | erf.abl_driver_type =  |                   |                     |
|                                     | GeostrophicWind)       |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.abl_geo_forcing**             | Constant body force    | 3 Reals           | (0.,0.,0.)          |
|                                     | applied to momentum    |                   |                     |
|                                     | equations              |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.abl_geo_wind_table**          | Path to text file      | String            | None                |
|                                     | containing a           |                   |                     |
|                                     | geostrophic wind       |                   |                     |
|                                     | profile                |                   |                     |
|                                     | (with z, Ug, and       |                   |                     |
|                                     | Vg whitespace          |                   |                     |
|                                     | delimited              |                   |                     |
|                                     | columns)               |                   |                     |
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
| **erf.use_gravity**                 | Include gravity        | Boolean           | false               |
|                                     | in momentum            |                   |                     |
|                                     | update?  If true,      |                   |                     |
|                                     | there is buoyancy      |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.use_coriolis**                | Include Coriolis       | Boolean           | false               |
|                                     | forcing                |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.variable_coriolis**           | Include Coriolis       | Boolean           | false               |
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
| **erf.rayleigh_damping_type**       | Rayleigh damping       | "SlowExplicit",   | "SlowExplicit"      |
|                                     | type. Leave all        | "FastExplicit",   |                     |
|                                     | rayleigh_damp_* flags  | "FastImplicit"    |                     |
|                                     | false to disable       |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.rayleigh_damp_U**             | Include explicit       | Boolean           | false               |
|                                     | Rayleigh damping in    |                   |                     |
|                                     | the x-momentum equation|                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.rayleigh_damp_V**             | Include explicit       | Boolean           | false               |
|                                     | Rayleigh damping in    |                   |                     |
|                                     | the y-momentum equation|                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.rayleigh_damp_W**             | Include                | Boolean           | false               |
|                                     | Rayleigh damping in    |                   |                     |
|                                     | the z-momentum equation|                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.rayleigh_damp_T**             | Include explicit       | Boolean           | false               |
|                                     | Rayleigh damping in    |                   |                     |
|                                     | the potential          |                   |                     |
|                                     | temperature equation   |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.rayleigh_dampcoef**           | Inverse damping        | Real [1/s]        | 0.2                 |
|                                     | timescale multiplying  |                   |                     |
|                                     | the vertical damping   |                   |                     |
|                                     | weight                 |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.rayleigh_zdamp**              | Depth of upper damping | Real [m]          | 500.0               |
|                                     | layer below model top  |                   |                     |
|                                     | where sine-squared     |                   |                     |
|                                     | ramp is nonzero        |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.nudging_from_input_sounding** | Add momentum source    | Boolean           | false               |
|                                     | terms to nudge the     |                   |                     |
|                                     | solution towards the   |                   |                     |
|                                     | initial sounding       |                   |                     |
|                                     | profile                |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.nudging_t_z1**                | Bottom of the height   | Real [m]          | 0.0                 |
|                                     | range over which       |                   |                     |
|                                     | potential temperature  |                   |                     |
|                                     | is nudged              |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.nudging_t_z2**                | Top of the same        | Real [m]          | 10000.0             |
|                                     | range                  |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.nudging_q_z1**                | Bottom of the height   | Real [m]          | 0.0                 |
|                                     | range over which       |                   |                     |
|                                     | water vapor is         |                   |                     |
|                                     | nudged                 |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.nudging_q_z2**                | Top of the same        | Real [m]          | 10000.0             |
|                                     | range                  |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.large_scale_forcing**         | Apply time-varying     | Boolean           | false               |
|                                     | large-scale tendencies |                   |                     |
|                                     | and subsidence read    |                   |                     |
|                                     | from a forcing file    |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.large_scale_forcing_file**    | Name of the            | String            | None                |
|                                     | large-scale forcing    |                   |                     |
|                                     | file                   |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.forcing_timescale**           | Relaxation time scale  | Real [s]          | 0.0                 |
|                                     | for the u and v        |                   |                     |
|                                     | large-scale nudging;   |                   |                     |
|                                     | 0 disables it          |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.input_sounding_file**         | Name(s) of the         | String(s)         | input_sounding      |
|                                     | input sounding file(s) |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.input_sounding_time**         | Time(s) of the         | Real(s)           | 0.0                 |
|                                     | input sounding file(s) |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.tau_nudging**                 | Time scale for         | Real              | 5.0                 |
|                                     | nudging                |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.bdy_nudge_factor**            | Sets real bc nudging   | Real              | 10.0                |
|                                     | strength as 1/(VAL*dt) |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.use_wrf_bdy_density**         | Use WRF-reconstructed  | Boolean           | true                |
|                                     | dry-air density for    |                   |                     |
|                                     | real WRF boundaries    |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.use_wrf_bdy_qc_qi**           | Ingest WRF ``QCLOUD``  | Boolean           | false               |
|                                     | and active ``QICE`` at |                   |                     |
|                                     | real boundaries        |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.bdy_rho_nudge_factor**        | Density Davies factor; | Real              | -1.0                |
|                                     | non-positive uses      |                   |                     |
|                                     | ``bdy_nudge_factor``   |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+
| **erf.bdy_moist_nudge_type**        | Which strategy for     | int 0,1,2 or 3    | 1                   |
|                                     | nudging of moist vars  |                   |                     |
+-------------------------------------+------------------------+-------------------+---------------------+

For WRF real boundaries, moisture nudging type 0 relaxes only ``qv`` toward
``QVAPOR``. Type 1 applies a ``qv`` tendency based on ERF ``qv+qc+qi`` versus
WRF ``QVAPOR``. Type 2 relaxes ``qv`` toward ``QVAPOR`` and active ``qc``/``qi``
toward zero with the existing latent-heat terms. Type 3 directly relaxes active
``qv``, ``qc``, and ``qi`` toward ``QVAPOR``, ``QCLOUD``, and ``QICE`` without
an added latent-heat term. Type 3 requires ``erf.use_wrf_bdy_qc_qi = true``.

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
| **erf.custom_forcing_uses_primitive_vars** | User-defined      | Boolean           | false       |
|                                            | source terms set  |                   |             |
|                                            | the tendency of   |                   |             |
|                                            | primitive         |                   |             |
|                                            | variables instead |                   |             |
|                                            | of conserved      |                   |             |
|                                            | quantities        |                   |             |
|                                            | (rho*prim_var)    |                   |             |
+--------------------------------------------+-------------------+-------------------+-------------+
| **erf.add_custom_rhotheta_forcing**        | Apply the         | Boolean           | false       |
|                                            | user-defined      |                   |             |
|                                            | temperature source|                   |             |
|                                            | term              |                   |             |
+--------------------------------------------+-------------------+-------------------+-------------+
| **erf.add_custom_moisture_forcing**        | Apply the         | Boolean           | false       |
|                                            | user-defined      |                   |             |
|                                            | qv source         |                   |             |
|                                            | term              |                   |             |
+--------------------------------------------+-------------------+-------------------+-------------+
| **erf.add_custom_w_subsidence**            | Apply the         | Boolean           | false       |
|                                            | user-defined      |                   |             |
|                                            | vertical velocity |                   |             |
|                                            | profile for use in|                   |             |
|                                            | calculating       |                   |             |
|                                            | subsidence source |                   |             |
|                                            | terms             |                   |             |
+--------------------------------------------+-------------------+-------------------+-------------+
| **erf.add_do_theta_advection**             | When using custom | Boolean           | true        |
|                                            | w subsidence,     |                   |             |
|                                            | apply the         |                   |             |
|                                            | subsidence source |                   |             |
|                                            | term to the       |                   |             |
|                                            | (rho*theta)       |                   |             |
|                                            | equation          |                   |             |
+--------------------------------------------+-------------------+-------------------+-------------+
| **erf.add_do_mom_advection**               | When using custom | Boolean           | true        |
|                                            | w subsidence,     |                   |             |
|                                            | apply the         |                   |             |
|                                            | subsidence source |                   |             |
|                                            | terms to the      |                   |             |
|                                            | momentum          |                   |             |
|                                            | equations         |                   |             |
+--------------------------------------------+-------------------+-------------------+-------------+
| **erf.add_custom_geostrophic_profile**     | Apply the         | Boolean           | false       |
|                                            | user-defined      |                   |             |
|                                            | geostrophic wind  |                   |             |
|                                            | profile           |                   |             |
+--------------------------------------------+-------------------+-------------------+-------------+

Note that ``erf.add_custom_geostrophic_profile`` cannot be used in combination
with an ``erf.abl_geo_wind_table``.

- Wind farm parameterization requires ``USE_WINDFARM=TRUE`` (gmake)
  or ``-DERF_ENABLE_WINDFARM`` (cmake) at build time.
  See :ref:`sec:WindFarmModels` for theory and examples.


Large Scale Forcing
--------------------

Large-scale forcing can be used to prescribe time-varying vertical profiles for
temperature, water vapor, horizontal wind, and vertical subsidence.
To enable it, set::

  erf.large_scale_forcing      = true
  erf.large_scale_forcing_file = lsf.txt
  erf.forcing_timescale        = 3600.0

where ``erf.large_scale_forcing_file`` is the path to a text file containing the
large scale forcing data and ``erf.forcing_timescale`` is the relaxation timescale
(in seconds) for nudging velocities. The file contains one or more time blocks, where each time
block is denoted with a header containing the day, number of levels in the block, and reference
pressure for that time. Each row of data in the block contains ``z, p, tls, qls, uls, vls, wls``
corresponding to level height, pressure, large scale temperature and moisture tendency, large scale
zonal and meridional velocity, and large scale vertical subsidence velocity.

- Either a height grid ``z[m]`` or pressure grid ``p[mb]`` can be provided, where one can be a negative number
  to indicate the other grid type should be used (e.g. if the ``z`` values are negative, the file is
  interpreted as pressure-level data instead of height-level data). If both grids are provided, then the ``z`` grid is used.

- The large scale forcing grid is interpolated to the ERF grid, similar to the input sounding.

- ``erf.forcing_timescale`` is the relaxation time scale used when nudging the
  horizontal velocities toward the ``uls`` and ``vls`` profiles. The
  temperature, moisture, and subsidence tendencies are applied directly and
  are not scaled by this option.

- Times are converted to elapsed seconds relative to the first timestamp in the file.

- Temperature and moisture nudging can be restricted to a custom interval instead of entire domain with the
  ``erf.nudging_t_z1`` and ``erf.nudging_t_z2``, and ``erf.nudging_q_z1`` and ``erf.nudging_q_z2`` options, respectively.

- **NOTE:** When both large scale forcing and ``erf.nudging_from_input_sounding`` are used, the u and v velocities are nudged
  according to the observed large scale velocity ``uls`` and ``vls`` instead of the input sounding, using
  ``erf.forcing_timescale``. Temperature and moisture are nudged from input sounding as normal, using the input
  sounding ``tau_nudging`` timescale.

Example file format for large scale forcing:

.. code-block:: text

   z[m] p[mb]  tls[K/s] qls[kg/kg/s] uls  vls  wls[m/s]
   0.0,   4, 1000.0  day, levels, pres0
      0.0  1000.00  0.000e+00  0.000e+00  0.0  0.0  0.000e+00
    100.0   990.00  0.000e+00  0.000e+00  0.5  0.2  0.000e+00
   1000.0   900.00 -5.000e-06  0.000e+00  2.0  1.0  0.000e+00
   5000.0   500.00 -1.000e-05  0.000e+00  6.0  3.0  0.000e+00
   0.5,   4, 1000.0  day, levels, pres0
      0.0  1000.00  0.000e+00  0.000e+00  0.2  0.1  0.000e+00
    100.0   990.00  0.000e+00  0.000e+00  0.7  0.3  0.000e+00
   1000.0   900.00 -2.500e-06  0.000e+00  2.5  1.2  0.000e+00
   5000.0   500.00 -8.000e-06  0.000e+00  6.5  3.2  0.000e+00


Boundary Plane I/O (Coupling Support)
=====================================

These options control writing and reading boundary-plane data for coupling
workflows (e.g., AMR-Wind). See :ref:`CouplingToAMRWind` for examples.

List of Parameters
------------------

+--------------------------------------+---------------------------------+---------------------+---------+
| Parameter                            | Definition                      | Acceptable Values   | Default |
+======================================+=================================+=====================+=========+
| **erf.output_bndry_planes**          | Enable boundary-plane output    | 0 or 1              | 0       |
+--------------------------------------+---------------------------------+---------------------+---------+
| **erf.input_bndry_planes**           | Enable boundary-plane input     | 0 or 1              | 0       |
+--------------------------------------+---------------------------------+---------------------+---------+
| **erf.bndry_output_planes_interval** | Output interval (steps)         | Integer >= 0        | -1      |
+--------------------------------------+---------------------------------+---------------------+---------+
| **erf.bndry_output_planes_per**      | Output interval (time)          | Real >= 0           | -1.0    |
+--------------------------------------+---------------------------------+---------------------+---------+
| **erf.bndry_output_start_time**      | Start time for output           | Real >= 0           | 0.0     |
+--------------------------------------+---------------------------------+---------------------+---------+
| **erf.bndry_output_planes_file**     | Output directory                | String              | None    |
+--------------------------------------+---------------------------------+---------------------+---------+
| **erf.bndry_output_box_lo**          | Lower-left (x,y) of output box  | 2 Reals             | None    |
+--------------------------------------+---------------------------------+---------------------+---------+
| **erf.bndry_output_box_hi**          | Upper-right (x,y) of output box | 2 Reals             | None    |
+--------------------------------------+---------------------------------+---------------------+---------+
| **erf.bndry_output_var_names**       | Variables to write              | List of strings     | All     |
+--------------------------------------+---------------------------------+---------------------+---------+
| **erf.bndry_file**                   | Input boundary-plane directory  | String              | None    |
+--------------------------------------+---------------------------------+---------------------+---------+
| **erf.bndry_input_var_names**        | Variables to read               | List of strings     | All     |
+--------------------------------------+---------------------------------+---------------------+---------+

Notes
-----

- If both interval controls are set, output occurs when either criterion is met.
- Output is written for the finest level that fully contains the requested box.

.. _sec:NumericalStability:

Numerical Stability
===================

To enhance numerical stability, e.g., for operational runs, we provide some additional controls.

List of Parameters
------------------

+----------------------------+----------------------------------+-------------------+-------------+
| Parameter                  | Definition                       | Acceptable Values | Default     |
+============================+==================================+===================+=============+
| **erf.beta_s**             | Time off-centering coefficient   | Real              | 0.1         |
+----------------------------+----------------------------------+-------------------+-------------+
| **erf.w_damping**          | Enable vertical-velocity         | Boolean           | false       |
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

The initialization strategy is determined at runtime by ``init_type``, which has ten possible values.

See :ref:`sec:Initialization` for more detail about how to provide initial conditions for an ERF simulation.

In addition, there is a run-time option to project the initial velocity field to make it divergence-free.

List of Parameters
------------------

+----------------------------------+---------------------+--------------------+-----------------------+
| Parameter                        | Definition          | Acceptable         | Default               |
|                                  |                     | Values             |                       |
+==================================+=====================+====================+=======================+
| **erf.init_type**                | Initialization      | "None",            | "None"                |
|                                  | type                | "Uniform"          |                       |
|                                  |                     | "ConstantDensity", |                       |
|                                  |                     | "Isentropic",      |                       |
|                                  |                     | "MoistBaseState",  |                       |
|                                  |                     | "Input_Sounding"   |                       |
|                                  |                     | "Hindcast",        |                       |
|                                  |                     | "WRFInput",        |                       |
|                                  |                     | "Metgrid"          |                       |
|                                  |                     | "NCFile"           |                       |
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
| **erf.use_real_bcs**             | If init_type is     | Boolean            | true if               |
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
| **erf.avg_grid_faces_to_nodes**  | avg z heights       |  Boolean           | false                 |
|                                  | from z-face data in |                    |                       |
|                                  | wrfinput/metgrid or |                    |                       |
|                                  | make our own grid?  |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.rebalance_wrf_input**      | rebalance state     |  Boolean           | true                  |
|                                  | from wrfinput and   |                    |                       |
|                                  | wrfbdy?             |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.real_extrap_w**            | First-order         |  Boolean           | true                  |
|                                  | extrapolation of    |                    |                       |
|                                  | vertical velocities |                    |                       |
|                                  | on lateral          |                    |                       |
|                                  | boundaries (instead |                    |                       |
|                                  | of setting to 0) if |                    |                       |
|                                  | use_real_bcs is     |                    |                       |
|                                  | true                |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.metgrid_debug_quiescent**  | If init_type is     |  Boolean           | false                 |
|                                  | Metgrid, overwrite  |                    |                       |
|                                  | initial conditions  |                    |                       |
|                                  | and boundary        |                    |                       |
|                                  | conditions to be    |                    |                       |
|                                  | quiescent.          |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.metgrid_debug_isothermal** | If init_type is     |  Boolean           | false                 |
|                                  | Metgrid, overwrite  |                    |                       |
|                                  | theta to be 300 in  |                    |                       |
|                                  | initial conditions  |                    |                       |
|                                  | and boundary        |                    |                       |
|                                  | conditions.         |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.metgrid_debug_dry**        | If init_type is     |  Boolean           | false                 |
|                                  | Metgrid, overwrite  |                    |                       |
|                                  | qv to be dry in     |                    |                       |
|                                  | initial conditions  |                    |                       |
|                                  | and boundary        |                    |                       |
|                                  | conditions.         |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.metgrid_debug_msf**        | If init_type is     |  Boolean           | false                 |
|                                  | Metgrid, overwrite  |                    |                       |
|                                  | map scale factors   |                    |                       |
|                                  | to be 1.            |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.metgrid_debug_psfc**       | If init_type is     |  Boolean           | false                 |
|                                  | Metgrid, overwrite  |                    |                       |
|                                  | surface pressure    |                    |                       |
|                                  | to be 10**5.        |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.metgrid_interp_theta**     | If init_type is     |  Boolean           | false                 |
|                                  | Metgrid, calculate  |                    |                       |
|                                  | theta on origin     |                    |                       |
|                                  | model vertical      |                    |                       |
|                                  | levels and then     |                    |                       |
|                                  | interpolate onto    |                    |                       |
|                                  | the ERF vertical    |                    |                       |
|                                  | levels.             |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.metgrid_basic_linear**     | If init_type is     |  Boolean           | false                 |
|                                  | Metgrid, use        |                    |                       |
|                                  | linear vertical     |                    |                       |
|                                  | interpolation and   |                    |                       |
|                                  | no quality          |                    |                       |
|                                  | control?            |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.metgrid_use_below_sfc**    | If init_type is     |  Boolean           | true                  |
|                                  | Metgrid, use the    |                    |                       |
|                                  | origin data levels  |                    |                       |
|                                  | below the surface?  |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.metgrid_use_sfc**          | If init_type is     |  Boolean           | true                  |
|                                  | Metgrid, use the    |                    |                       |
|                                  | origin data level   |                    |                       |
|                                  | at the surface?     |                    |                       |
+----------------------------------+---------------------+--------------------+-----------------------+
| **erf.metgrid_retain_sfc**       | If init_type is     |  Boolean           | false                 |
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
the extent of the relaxation zone may be controlled with ``erf.real_width`` (corresponding to WRF's **spec_bdy_width**).

Note that the **erf.project_initial_velocity** option is available for all **init_type** options.  If using the anelastic
formulation this will be true regardless of the input; if using the compressible formulation the default is false but
that can be over-written.

Map Scale Factors
=================

Map scale factors are always present in the evolution equations, but the values default to one
unless specified in the initialization when **erf.init_type = WRFInput** or **erf.init_type = Metgrid**.

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

The height at the surface nodes can be defined analytically, read from a text file
specified by ``erf.terrain_file_name``, or read from a NetCDF file specified by
``erf.terrain_file_name_nc``. The NetCDF option requires NetCDF support and takes
precedence when both inputs are present. The explicit Cartesian format uses a
height variable named ``height``, ``z``, ``terrain``, ``HGT_M``, or ``AGL`` with
dimensions ``(y,x)`` or ``(time,y,x)`` (time must be the leading dimension).
It must also provide one-dimensional, finite, strictly increasing, uniformly
spaced ``x`` and ``y`` coordinate variables matching the height dimensions;
ERF reads the first time record when one is present.

A genuine WPS ``geo_em`` file may instead provide
``HGT_M(Time,south_north,west_east)`` without Cartesian ``x``/``y`` variables.
It must provide matching first-time-record ``XLAT_M`` and ``XLONG_M`` fields and
the WPS global attributes ``MAP_PROJ``, ``CEN_LAT``, ``CEN_LON``, ``STAND_LON``,
``DX``, ``DY``, ``WEST-EAST_GRID_DIMENSION``, and
``SOUTH-NORTH_GRID_DIMENSION``. ERF supports Lambert (1), polar (2), and
Mercator (3) WPS projections; Lambert additionally requires ``TRUELAT1`` and
``TRUELAT2``, while the other two require ``TRUELAT1``. The latitude/longitude
projection (6) is rejected, because WPS writes ``DX``/``DY`` in degrees for it
while ERF consumes them as metres. ``HGT_M``, ``XLAT_M``, and ``XLONG_M`` must
have the exact ``(Time,south_north,west_east)`` order, and the first time
record is selected.
The WPS grid dimensions, ``DX``/``DY``, and domain extent must match the ERF
level-0 Cartesian geometry.

WPS fields are mass-point fields. ERF maps their logical mass-point locations
to terrain nodes using ``x_mass = x_lo + (i+1/2) DX`` and
``y_mass = y_lo + (j+1/2) DY``. At the four outer node edges, where no WPS mass
point exists, ERF holds the nearest mass-point value (an explicit one-sided
mass-to-node conversion); this is not general extrapolation. For both formats,
bilinear interpolation includes the final coordinate. Explicit Cartesian
terrain nodes genuinely outside the file grid are assigned zero, and malformed
or nonfinite WPS metadata, coordinates, or heights are rejected.

For the text format, the ordering of data in the file is first
nx, ny (where nx is the number of values specified in the x-direction and
ny is the number of values specified in the y-direction; note that these need not match
the resolution of the problem).  Then nx x-values are read followed by ny y-values.   Finally,
the z-coordinates are read in the order z(x1,y1), z(x1,y2), z(x1,y3), ... z(x2,y1), ... z(nx,ny).
An example is given in Exec/CanonicalTests/ABL/erf_terrain_def.

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
| **erf.terrain_file_name_nc**| NetCDF terrain     | String             | NONE       |
|                             | filename           |                    |            |
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
|                                | fluxes                     | "SLM"              |             |
+--------------------------------+----------------------------+--------------------+-------------+
| **erf.use_coupled_sst**        | Expect sea-surface         | true / false       | false       |
|                                | temperature from an        |                    |             |
|                                | external ocean coupler     |                    |             |
+--------------------------------+----------------------------+--------------------+-------------+

.. note::

   Noah-MP requires ``USE_NOAHMP=TRUE`` at build time. See :ref:`CouplingToNoahMP` for details.

.. note::

   ``erf.use_coupled_sst`` is independent of ``erf.land_surface_model``. Coupled SST
   is applied through the lower-boundary path alongside ``wrflowinp`` SST/TSK, and
   only on the water cells the coupler actually covers: land keeps its land surface
   model, and water the ocean grid does not reach keeps the ``wrflowinp`` value. A
   coupled run can therefore also run Noah-MP.

   This replaces ``erf.land_surface_model = OceanSurf``, which delivered coupled SST
   as a land surface model and so occupied the one land surface model slot. That
   value has been removed and now aborts with a message pointing here.

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
|                                 |                          |  "SAM_NoIce", "WSM6", |            |
|                                 |                          |  "WDM6", "P3",        |            |
|                                 |                          |  "MoistNoCondensation"|            |
+---------------------------------+--------------------------+-----------------------+------------+
| **erf.moisture_tight_coupling** | If true, advance         |  Boolean              | false      |
|                                 | microphysics after every |                       |            |
|                                 | slow step in the dycore; |                       |            |
|                                 | otherwise, update after  |                       |            |
|                                 | the dycore has been      |                       |            |
|                                 | advanced at each timestep|                       |            |
+---------------------------------+--------------------------+-----------------------+------------+
| **erf.morrison_ndcnst**         | Constant cloud-droplet   | Positive real         | 250        |
|                                 | number concentration for | (cm^-3)               | (cm^-3)    |
|                                 | Morrison when constant   |                       |            |
|                                 | droplet number is active |                       |            |
+---------------------------------+--------------------------+-----------------------+------------+
| **wdm6.hail_opt**               | Graupel/hail regime      | 0, 1                  | 0          |
|                                 | selector for WDM6:       |                       |            |
|                                 | 0 = graupel regime,      |                       |            |
|                                 | 1 = hail regime          |                       |            |
+---------------------------------+--------------------------+-----------------------+------------+
| **wdm6.ccn0**                   | Background CCN number    | Positive real         | 100.0e6    |
|                                 | concentration for WDM6   | (m^-3)                | (m^-3)     |
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

+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| Parameter                           | Definition                             | Acceptable Values | Default / Notes                   |
+=====================================+========================================+===================+===================================+
| **erf.radiation_model**             | Enable radiation model                 | "None", "RRTMGP", | "None". "Simple" is a prescribed  |
|                                     |                                        | "Simple"          | cooling/heating profile with no   |
|                                     |                                        |                   | radiative transfer.               |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.rad_nvar**                    | Size of block memory allocation        | Integer > 0       | 12                                |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.rad_t_sfc**                   | Surface temperature if no LSM          | Real              | Must be set without LSM           |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.rad_freq_in_steps**           | Radiation update frequency (steps)     | Integer >= 1      | 1                                 |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.rad_ncol_chunk**              | Columns per RRTMGP kernel launch.      | Integer >= 1      | 5000. Lower values reduce peak    |
|                                     | Controls peak GPU memory by processing |                   | GPU memory; higher values reduce  |
|                                     | radiation in batches of this size.     |                   | kernel launch overhead.           |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.rad_write_fluxes**            | Write radiation fluxes to plotfiles    | Boolean           | false                             |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.rad_do_subcol_sampling**      | Enable MCICA subcolumn sampling        | Boolean           | true                              |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.rad_orbital_year**            | Fixed orbital year for zenith calcs    | Integer           | < 0 uses timestamp year           |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.rad_orbital_eccentricity**    | Override orbital eccentricity          | Real              | < 0 uses computed value           |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.rad_orbital_obliquity**       | Override orbital obliquity             | Real              | < 0 uses computed value           |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.rad_orbital_mvelp**           | Override mean longitude of perihelion  | Real              | < 0 uses computed value           |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.rad_cons_lat**                | Constant latitude for idealized cases  | Real              | 39.809860                         |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.rad_cons_lon**                | Constant longitude for idealized cases | Real              | -98.555183                        |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.fixed_total_solar_irradiance**| Fixed total solar irradiance (TOA)     | Real              | < 0 disables                      |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.fixed_solar_zenith_angle**    | Fixed solar zenith (passed as ``mu0``) | Real              | <= 0 disables                     |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.co2vmr**                      | CO2 volume mixing ratio                | Real              | 388.717e-6                        |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.o3vmr**                       | O3 volume mixing ratio profile         | Real or list      | From dataset if unset             |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.n2ovmr**                      | N2O volume mixing ratio                | Real              | 323.141e-9                        |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.covmr**                       | CO volume mixing ratio                 | Real              | 1.0e-7                            |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.ch4vmr**                      | CH4 volume mixing ratio                | Real              | 1807.851e-9                       |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.o2vmr**                       | O2 volume mixing ratio                 | Real              | 0.209448                          |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.n2vmr**                       | N2 volume mixing ratio                 | Real              | 0.7906                            |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.rad_do_aerosol**              | Enable aerosol forcing in radiation    | Boolean           | true                              |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.rad_extra_clnclrsky_diag**    | Extra clean and clear-sky diagnostics  | Boolean           | false                             |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.rad_extra_clnsky_diag**       | Extra clean-sky diagnostics            | Boolean           | false                             |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.rrtmgp_file_path**            | Path to RRTMGP data files              | String            | "."                               |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.rrtmgp_coeffs_sw**            | Shortwave k-distribution file          | String            | rrtmgp-data-sw-g224-2018-12-04.nc |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.rrtmgp_coeffs_lw**            | Longwave k-distribution file           | String            | rrtmgp-data-lw-g256-2018-12-04.nc |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.rrtmgp_cloud_optics_sw**      | Shortwave cloud optics file            | String            | rrtmgp-cloud-optics-coeffs-sw.nc  |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+
| **erf.rrtmgp_cloud_optics_lw**      | Longwave cloud optics file             | String            | rrtmgp-cloud-optics-coeffs-lw.nc  |
+-------------------------------------+----------------------------------------+-------------------+-----------------------------------+

Notes
=====

- ``erf.fixed_solar_zenith_angle`` is passed directly as ``mu0`` (cosine of the zenith angle).
- If ``erf.rad_orbital_year`` is negative, the orbital year is taken from the simulation timestamp.
  This is keyed off the value itself, not off whether the input was specified.
- When orbital parameters are negative, values are computed from the orbital year.  Each of
  ``erf.rad_orbital_eccentricity``, ``erf.rad_orbital_obliquity`` and ``erf.rad_orbital_mvelp``
  is honored independently, so specifying one does not change how the others are obtained.

The lookup data may be downloaded as a package from `here <https://doi.org/10.22002/ppv8a-4q131>`_.
The k-distribution files are also shipped with the RRTMGP submodule in
``Submodules/RRTMGP/rrtmgp/data``, and the cloud optics files in
``Submodules/RRTMGP/extensions/cloud_optics``.

All four files are expected in the single directory named by ``erf.rrtmgp_file_path``, which
is resolved relative to the run directory (not to the location of the inputs file). ERF checks
for them while constructing the radiation interface and aborts with a message listing the
directory and any missing file, so a mis-set path is reported at startup rather than as a
netCDF "No such file or directory" error.

.. note::

   Using RRTMGP requires ``USE_RRTMGP=TRUE`` at build time. See :ref:`sec:building` for build instructions.

SHOC
======================

ERF provides two SHOC implementations. Select the ERF-native implementation with
``erf.pbl_type = NATIVE_SHOC``. Select the optional EAMxx implementation with
``erf.pbl_type = EAMXX_SHOC``. The legacy value ``SHOC`` is deprecated and is
interpreted as ``EAMXX_SHOC``; it does not select native SHOC.

Both SHOC implementations use the ``erf.shoc`` namespace for closure tuning.
Native SHOC also has transport and native-only numerical controls. EAMxx SHOC
has a small number of interface-only controls. See :ref:`SHOC` for the model
description, coupling restrictions, transport behavior, and diagnostic output.

Shared SHOC closure inputs
--------------------------

The following option names are read by both native SHOC and EAMxx SHOC. The two
implementations are separate code paths, so identical input values do not imply
bit-for-bit identical calculations.

.. list-table::
   :header-rows: 1
   :widths: 36 44 18 12

   * - Parameter
     - Definition
     - Acceptable Values
     - Default
   * - **erf.shoc.lambda_low**
     - Lower bound for the SHOC stability-correction coefficient
     - Real > 0
     - 0.001
   * - **erf.shoc.lambda_high**
     - Upper bound for the SHOC stability-correction coefficient
     - Real >= ``lambda_low``
     - 0.04
   * - **erf.shoc.lambda_slope**
     - Slope used to vary the stability-correction coefficient
     - Real
     - 2.65
   * - **erf.shoc.lambda_thresh**
     - Threshold used in the stability-correction calculation
     - Real
     - 0.02
   * - **erf.shoc.thl2tune**
     - Liquid-water potential-temperature variance tuning factor
     - Real
     - 1.0
   * - **erf.shoc.qw2tune**
     - Total-water variance tuning factor
     - Real
     - 1.0
   * - **erf.shoc.qwthl2tune**
     - Total-water / liquid-water-potential-temperature covariance tuning factor
     - Real
     - 1.0
   * - **erf.shoc.w2tune**
     - Vertical-velocity variance tuning factor
     - Real
     - 1.0
   * - **erf.shoc.length_fac**
     - Turbulent length-scale tuning factor
     - Real > 0
     - 0.5
   * - **erf.shoc.c_diag_3rd_mom**
     - Diagnostic third-moment closure coefficient
     - Real
     - 7.0
   * - **erf.shoc.coeff_kh**
     - Scalar / heat diffusivity tuning coefficient
     - Real >= 0
     - 0.1
   * - **erf.shoc.coeff_km**
     - Momentum diffusivity tuning coefficient
     - Real >= 0
     - 0.1
   * - **erf.shoc.shoc_1p5tke**
     - Use the reduced 1.5-order TKE form instead of the higher-order SHOC moment form
     - true / false
     - false
   * - **erf.shoc.extra_shoc_diags**
     - Enable additional SHOC diagnostic calculations; plot variables must still be requested explicitly
     - true / false
     - false

Native SHOC inputs
------------------

The following controls apply to ``erf.pbl_type = NATIVE_SHOC``.

.. list-table::
   :header-rows: 1
   :widths: 36 52 24 14

   * - Parameter
     - Definition
     - Acceptable Values
     - Default
   * - **erf.shoc.transport_mode**
     - Selects how native SHOC transports thermodynamic variables, moisture/cloud state, and TKE
     - ``state_update``, ``host_diffusion``
     - ``state_update``
   * - **erf.shoc.momentum_transport**
     - Selects how native SHOC applies horizontal-momentum transport
     - ``none``, ``state_update``, ``host_diffusion``
     - ``host_diffusion``
   * - **erf.shoc.top_taper_depth**
     - Depth below the model top over which native SHOC smoothly tapers TKE/diffusivity and higher-order turbulence quantities [m]; zero disables the taper
     - Real >= 0
     - 0.0
   * - **erf.shoc.top_taper_min_factor**
     - Multiplicative factor at the model top when the native top taper is enabled; the factor reaches 1 below the taper layer
     - Real in [0, 1]
     - 0.0
   * - **erf.shoc.signed_tke_production**
     - If false, clip the sum of native SHOC shear and buoyancy production at zero before dissipation; if true, retain the signed production sum
     - true / false
     - false
   * - **erf.shoc.debug_summary**
     - Print a synchronized min/max summary after each native SHOC advance; intended for debugging rather than routine production
     - true / false
     - false

For moist native SHOC runs, use ``erf.shoc.transport_mode = state_update``.
The full ``host_diffusion`` transport mode is currently supported only when
``erf.moisture_model = None`` because native SHOC does not own cloud
macrophysics in this mode while SHOC-family microphysics condensation is
suppressed. It also requires ``erf.shoc.momentum_transport = host_diffusion``.
Native ``state_update`` rejects moisture layouts containing cloud-water or
cloud-ice number concentrations because a number closure has not yet been
implemented.

Native SHOC ``pblh`` is reported in metres above local ground (AGL).

The old value ``erf.shoc.transport_mode = tendencies`` is no longer accepted.
Use ``state_update`` instead.

EAMxx SHOC interface inputs
---------------------------

The following controls have a verified runtime effect in the optional
``EAMXX_SHOC`` interface and are not native-SHOC controls.

.. list-table::
   :header-rows: 1
   :widths: 36 52 18 12

   * - Parameter
     - Definition
     - Acceptable Values
     - Default
   * - **erf.shoc.apply_tms**
     - Apply the turbulent mountain stress contribution in the EAMxx SHOC interface
     - true / false
     - false
   * - **erf.shoc.check_flux_state**
     - Run the EAMxx SHOC interface flux/state consistency check
     - true / false
     - false

Only options documented above are part of the supported SHOC user-facing
runtime contract. Some historical or development ``erf.shoc`` keys may still
be parsed internally but are intentionally not documented as active features.

Minimal native SHOC setup
-------------------------

A native SHOC run requires a surface-layer lower boundary:

.. code-block:: text

   zlo.type = "surface_layer"
   erf.pbl_type = NATIVE_SHOC

The native transport defaults are:

.. code-block:: text

   erf.shoc.transport_mode = state_update
   erf.shoc.momentum_transport = host_diffusion

These two transport lines may be omitted when the defaults are desired.

Minimal EAMxx SHOC setup
------------------------

The EAMxx SHOC path requires an executable built with the CMake option
``ERF_ENABLE_EAMXX_SHOC=ON``. For example:

.. code-block:: bash

   cmake -DERF_ENABLE_EAMXX_SHOC=ON ...

The corresponding runtime selection is:

.. code-block:: text

   zlo.type = "surface_layer"
   erf.pbl_type = EAMXX_SHOC

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

Reproducibility
===============

List of Parameters
------------------

+--------------------------+-------------------------------+---------------------+---------+
| Parameter                | Definition                    | Acceptable Values   | Default |
+==========================+===============================+=====================+=========+
| **erf.fix_random_seed**  | Use a fixed random seed for   | 0 or 1              | 0       |
|                          | reproducible runs             |                     |         |
+--------------------------+-------------------------------+---------------------+---------+
| **erf.random_seed**      | Seed for the random number    | Integer >= 0        | -1      |
|                          | generator, offset by MPI      |                     |         |
|                          | rank; ignored if              |                     |         |
|                          | ``erf.fix_random_seed`` is 1  |                     |         |
+--------------------------+-------------------------------+---------------------+---------+

Embedded Boundary (EB) Tuning
=============================

List of Parameters
------------------

+------------------------+-----------------------------------------+---------------------+-----------+
| Parameter              | Definition                              | Acceptable Values   | Default   |
+========================+=========================================+=====================+===========+
| **eb2.small_volfrac**  | Volume-fraction threshold used          | Real > 0            | 1.0e-14   |
|                        | to treat cells as effectively empty     |                     |           |
+------------------------+-----------------------------------------+---------------------+-----------+

Particles
=========

List of Parameters
------------------

+--------------------------+------------------------------+---------------------+---------+
| Parameter                | Definition                   | Acceptable Values   | Default |
+==========================+==============================+=====================+=========+
| **particles.disable_plt**| Disable particle plotfile    | Boolean             | false   |
|                          | output                       |                     |         |
+--------------------------+------------------------------+---------------------+---------+

.. _sec:EnsembleInitialization:

Ensemble Initialization
=======================

When ``erf.is_init_for_ensemble`` is true, the initial state is perturbed with a
spatially correlated random field so that each member of an ensemble starts from
a slightly different state.  In that case all of the ``ensemble.*`` parameters
below must be set.

List of Parameters
------------------

+-------------------------------------------+--------------------------------+---------------------+-------------+
| Parameter                                 | Definition                     | Acceptable Values   | Default     |
+===========================================+================================+=====================+=============+
| **erf.is_init_for_ensemble**              | Perturb the initial state for  | Boolean             | false       |
|                                           | an ensemble member             |                     |             |
+-------------------------------------------+--------------------------------+---------------------+-------------+
| **ensemble.n_members**                    | Number of ensemble members     | Integer >= 2        | must be set |
+-------------------------------------------+--------------------------------+---------------------+-------------+
| **ensemble.coarse_bckgnd_data_file**      | Coarse background data file    | String              | must be set |
|                                           | interpolated onto the fine     |                     |             |
|                                           | mesh to form the background    |                     |             |
|                                           | state                          |                     |             |
+-------------------------------------------+--------------------------------+---------------------+-------------+
| **ensemble.ens_pert_amplitude**           | Amplitude of the perturbation  | Real > 0            | must be set |
|                                           | added to the background state  |                     |             |
+-------------------------------------------+--------------------------------+---------------------+-------------+
| **ensemble.ens_pert_correlated_radius**   | Spatial correlation radius of  | Real > 0            | must be set |
|                                           | the perturbation               |                     |             |
+-------------------------------------------+--------------------------------+---------------------+-------------+

.. _sec:SolverChoice:

Solver Options in ``SolverChoice``
==================================

The parameters in this section are the complete set of inputs read by
``SolverChoice`` in ``Source/DataStructs/ERF_DataStruct.H``, which holds most of
the algorithmic and physics-selection options for a run.  Unless noted
otherwise, every parameter below uses the ``erf.`` prefix.  Many of these
parameters also appear in the topic-oriented sections above, and in the theory
chapters, where their meaning is discussed in more detail; this section is
intended as a single reference list.

.. note::

   Several of the parameters below may be specified either as a single value,
   which is then used at every AMR level, or as a space-separated list with one
   value per level, e.g. ``erf.anelastic = 1 0 0``.  These are identified as
   "per-level" in the tables below.  For a list, the number of values given must
   be at least ``amr.max_level + 1``.

Equation Set
------------

+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| Parameter                                      | Definition                                             | Acceptable Values              | Default                |
+================================================+========================================================+================================+========================+
| **erf.anelastic**                              | if 1, solve the anelastic equations rather than the    | 0, 1                           | 0                      |
|                                                | fully compressible equations (per-level).  Setting     |                                |                        |
|                                                | this to 1 at a level also forces                       |                                |                        |
|                                                | ``project_initial_velocity`` = 1, ``fixed_density`` =  |                                |                        |
|                                                | 1, ``buoyancy_type`` = 3 and ``substepping_type`` =    |                                |                        |
|                                                | None at that level                                     |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.fixed_density**                          | if 1, hold the density fixed in time (per-level)       | 0, 1                           | 0 (1 if anelastic)     |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.project_initial_velocity**               | if 1, project the initial velocity field to make it    | 0, 1                           | 0 (1 if anelastic)     |
|                                                | satisfy the anelastic constraint (per-level)           |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.buoyancy_type**                          | which expression to use for the buoyancy term          | 1, 2, 3, 4                     | 1 (3 if anelastic)     |
|                                                | (per-level); see :ref:`Buoyancy`                       |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.use_gravity**                            | include the gravitational source term; primarily an    | Boolean                        | false                  |
|                                                | option for unit testing.  Forced to true if            |                                |                        |
|                                                | ``init_type`` is ``WRFInput`` or ``Metgrid``           |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.c_p**                                    | specific heat at constant pressure for dry air         | Real > 0                       | 1004.5                 |
|                                                | [J/(kg-K)]                                             |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.gradp_type**                             | which horizontal pressure gradient formulation to use  | 0, 1                           | 0                      |
|                                                | with terrain-fitted coordinates: 0 for dp/dx with a    |                                |                        |
|                                                | dp/dz correction, 1 for the gradient of the vertically |                                |                        |
|                                                | interpolated pressure (Klemp 2011)                     |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.use_pert_pres_gradient**                 | if true, the lateral pressure gradient in the momentum | Boolean                        | true                   |
|                                                | equation uses horizontal derivatives of the            |                                |                        |
|                                                | perturbational pressure; if false, of the full         |                                |                        |
|                                                | pressure                                               |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.transport_scalar**                       | transport the passive scalar component                 | Boolean                        | true                   |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.use_lagged_delta_rt**                    | use the lagged (rho*theta) increment in the fast       | Boolean                        | true                   |
|                                                | integrator; may only be set to false when              |                                |                        |
|                                                | ``terrain_type`` = ``MovingFittedMesh``                |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.coupling_type**                          | how data is exchanged between AMR levels.  Forced to   | OneWay, TwoWay                 | TwoWay                 |
|                                                | ``OneWay`` if some levels are anelastic and others     |                                |                        |
|                                                | compressible                                           |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+

Acoustic Substepping and the Poisson Solve
------------------------------------------

+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| Parameter                                      | Definition                                             | Acceptable Values              | Default                |
+================================================+========================================================+================================+========================+
| **erf.substepping_type**                       | acoustic substepping strategy for the compressible     | None, Implicit                 | Implicit               |
|                                                | solver (per-level).  Forced to ``None`` at any level   |                                |                        |
|                                                | where ``anelastic`` = 1                                |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.substepping_diag**                       | print additional CFL diagnostics for the compressible  | Boolean                        | false                  |
|                                                | solver with substepping                                |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.beta_s**                                 | time off-centering coefficient used in the acoustic    | Real                           | 0.1                    |
|                                                | substepping; positive values weight the update towards |                                |                        |
|                                                | the new time                                           |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.force_stage1_single_substep**            | if 1, take a single substep in the first Runge-Kutta   | 0, 1                           | 1                      |
|                                                | stage                                                  |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.ncorr**                                  | number of correction iterations in the projection      | Integer >= 1                   | 1                      |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.poisson_abstol**                         | absolute tolerance for the Poisson solve; raised to    | Real > 0                       | 1.e-8                  |
|                                                | 1.e-6 in single precision builds                       |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.poisson_reltol**                         | relative tolerance for the Poisson solve; raised to    | Real > 0                       | 1.e-8                  |
|                                                | 1.e-6 in single precision builds                       |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+

Initialization, Terrain and Vertical Mesh
-----------------------------------------

+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| Parameter                                      | Definition                                             | Acceptable Values              | Default                |
+================================================+========================================================+================================+========================+
| **erf.init_type**                              | source of the initial data; must be set.  See          | Input_Sounding, NCFile,        | must be set            |
|                                                | :ref:`sec:Initialization`                              | WRFInput, Metgrid, Uniform,    |                        |
|                                                |                                                        | ConstantDensity,               |                        |
|                                                |                                                        | ConstantDensityLinearTheta,    |                        |
|                                                |                                                        | Isentropic, MoistBaseState,    |                        |
|                                                |                                                        | HindCast                       |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.sounding_type**                          | how the profiles in an input sounding file are         | ConstantDensity, Ideal,        | Ideal                  |
|                                                | interpreted; read only when ``init_type`` =            | Isentropic, DryIsentropic      |                        |
|                                                | ``Input_Sounding``                                     |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.nc_bdy_file**                            | name of the NetCDF (``wrfbdy``) lateral boundary file; | String                         | None                   |
|                                                | read only when ``init_type`` = ``WRFInput``.  If not   |                                |                        |
|                                                | given, ``use_real_bcs`` is set to false                |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.use_real_bcs**                           | use the real-data lateral boundary conditions from the | Boolean                        | set by ``init_type``   |
|                                                | boundary file.  Defaults to true when ``init_type`` is |                                |                        |
|                                                | ``WRFInput`` or ``Metgrid``, false otherwise, and may  |                                |                        |
|                                                | only be used to turn the real boundary conditions off  |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.terrain_type**                           | how terrain and immersed boundaries are represented.   | None, StaticFittedMesh,        | None                   |
|                                                | Only ``StaticFittedMesh`` is allowed with              | MovingFittedMesh, EB,          |                        |
|                                                | ``init_type`` = ``WRFInput`` or ``Metgrid``            | ImmersedForcing                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.buildings_type**                         | how buildings are represented; see :ref:`Forcings`     | None, ImmersedForcing          | None                   |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.grid_stretching_ratio**                  | ratio by which each vertical cell size is increased    | 0 or Real >= 1                 | 0 (no stretching)      |
|                                                | relative to the one below it.  If set, the mesh type   |                                |                        |
|                                                | becomes ``StretchedDz`` and ``erf.initial_dz`` must    |                                |                        |
|                                                | also be given                                          |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.zsurface**                               | nominal surface height used when generating a          | Real                           | 0.0                    |
|                                                | stretched vertical mesh; read only when                |                                |                        |
|                                                | ``grid_stretching_ratio`` >= 1                         |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.initial_dz**                             | vertical cell size of the lowest cell when generating  | Real > 0                       | must be set            |
|                                                | a stretched vertical mesh; required when               |                                |                        |
|                                                | ``grid_stretching_ratio`` >= 1                         |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.terrain_z_levels**                       | explicit list of nominal (terrain-free) staggered z    | List of Reals                  | None                   |
|                                                | levels; specifying these makes the mesh type           |                                |                        |
|                                                | ``StretchedDz``.  See :ref:`sec:Meshing`               |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.avg_grid_faces_to_nodes**                | average the ``wrfinput`` / ``met_em`` heights onto the | Boolean                        | false                  |
|                                                | nodes rather than reconstructing nodal heights whose   |                                |                        |
|                                                | four-node average reproduces them                      |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.rebalance_wrf_input**                    | rebalance (hydrostatically re-integrate) the state     | Boolean                        | true                   |
|                                                | read from ``wrfinput`` and ``wrfbdy``.  Forced to true |                                |                        |
|                                                | if ``avg_grid_faces_to_nodes`` is false                |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+

Physics Model Selection
-----------------------

+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| Parameter                                      | Definition                                             | Acceptable Values              | Default                |
+================================================+========================================================+================================+========================+
| **erf.moisture_model**                         | which moisture / microphysics model to use             | None, SatAdj, Kessler,         | None                   |
|                                                |                                                        | Kessler_NoRain, SAM,           |                        |
|                                                |                                                        | SAM_NoIce, SAM_NoPrecip_NoIce, |                        |
|                                                |                                                        | Morrison, Morrison_NoIce,      |                        |
|                                                |                                                        | WSM6, WDM6, SuperDroplets,     |                        |
|                                                |                                                        | MoistNoCondensation            |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.moisture_tight_coupling**                | tightly couple the moisture update to the dynamics;    | Boolean                        | false                  |
|                                                | read only when a moisture model is active              |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.land_surface_model**                     | which land surface model to use                        | None, SLM, NOAHMP              | None                   |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.use_coupled_sst**                        | an external ocean coupler supplies the sea surface     | Boolean                        | false                  |
|                                                | temperature for some or all water cells; this must be  |                                |                        |
|                                                | set at problem definition time because the surface     |                                |                        |
|                                                | layer uses it to select its temperature boundary       |                                |                        |
|                                                | treatment                                              |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.is_land**                                | whether each AMR level is treated as land (1) or water | 0, 1                           | 1                      |
|                                                | (0) (per-level)                                        |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.radiation_model**                        | which radiation model to use; ``RRTMGP`` requires that | None, RRTMGP, Simple           | None                   |
|                                                | ERF was built with RRTMGP enabled                      |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.four_stream_radiation**                  | use the four-stream radiation approximation            | Boolean                        | false                  |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+

Implicit Vertical Diffusion
---------------------------

These parameters control the time-centering of the vertical differences in the
diffusive terms.  The implicit solve is turned off automatically if there is no
molecular or turbulent diffusion, if the level is anelastic, if ``terrain_type``
= ``EB``, or if the active PBL scheme owns *all* of the vertical diffusion
itself.  That last case covers ``erf.pbl_type = EAMXX_SHOC`` (which always owns
both scalar and momentum diffusion) and ``erf.pbl_type = NATIVE_SHOC`` when
``erf.shoc.transport_mode = state_update`` *and*
``erf.shoc.momentum_transport`` is ``state_update`` or ``none``.

When a SHOC-family scheme is active, the ``erf.implicit_*_diffusion`` flags are
additionally restricted to the components that SHOC hands back to the host:
``erf.shoc.momentum_transport = host_diffusion`` (the native default) keeps
``erf.implicit_momentum_diffusion``, and ``erf.shoc.transport_mode =
host_diffusion`` keeps ``erf.implicit_thermal_diffusion``,
``erf.implicit_moisture_diffusion``, and ``erf.implicit_ke_diffusion``.
Requesting a component that SHOC owns is ignored, with a message in the run
log.  ``erf.vert_implicit`` and ``erf.vert_implicit_fac`` are always honored,
so the solve can still be disabled entirely from the inputs file.

+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| Parameter                                      | Definition                                             | Acceptable Values              | Default                |
+================================================+========================================================+================================+========================+
| **erf.vert_implicit**                          | if false, turn off the vertical implicit solve at all  | Boolean                        | true                   |
|                                                | levels                                                 |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.vert_implicit_fac**                      | time-centering factor for the vertical diffusive       | 1 or 3 Reals in [0,1]          | 1.0 1.0 0.0            |
|                                                | terms, where 0 is fully explicit and 1 is fully        |                                |                        |
|                                                | implicit.  Specify either one value used in all        |                                |                        |
|                                                | Runge-Kutta stages, or three values, one per stage     |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.implicit_thermal_diffusion**             | include the implicit contribution to vertical thermal  | Boolean                        | true                   |
|                                                | diffusion                                              |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.implicit_moisture_diffusion**            | include the implicit contribution to vertical moisture | Boolean                        | true                   |
|                                                | diffusion                                              |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.implicit_ke_diffusion**                  | include the implicit contribution to vertical          | Boolean                        | true                   |
|                                                | turbulent kinetic energy diffusion                     |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.implicit_momentum_diffusion**            | include the implicit contributions in tau13 and tau23  | Boolean                        | true                   |
|                                                | (and tau33 if built with ``ERF_IMPLICIT_W``) used to   |                                |                        |
|                                                | correct the momenta                                    |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.implicit_before_substep**                | if true, the vertical implicit diffusive solve is done | Boolean                        | true                   |
|                                                | before the acoustic substepping rather than after.     |                                |                        |
|                                                | Forced to true if any level has ``substepping_type`` = |                                |                        |
|                                                | ``None``                                               |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+

Coriolis and Large-Scale ABL Forcing
------------------------------------

+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| Parameter                                      | Definition                                             | Acceptable Values              | Default                |
+================================================+========================================================+================================+========================+
| **erf.use_coriolis**                           | include the Coriolis forcing                           | Boolean                        | false                  |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.coriolis_3d**                            | include the vertical (cosine of latitude) Coriolis     | Boolean                        | true                   |
|                                                | terms in addition to the traditional terms             |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.variable_coriolis**                      | allow the Coriolis frequency to vary spatially rather  | Boolean                        | false                  |
|                                                | than using a single latitude                           |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.rotational_time_period**                 | rotational period of the planet [s], used to form the  | Real > 0                       | 86400.0                |
|                                                | Coriolis factor; read only if ``use_coriolis`` is true |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.latitude**                               | latitude [degrees] at which the constant-latitude      | Real                           | 90.0                   |
|                                                | Coriolis frequency is evaluated; read only if          |                                |                        |
|                                                | ``use_coriolis`` is true                               |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.abl_driver_type**                        | type of large-scale forcing used to drive the boundary | None, PressureGradient,        | None                   |
|                                                | layer                                                  | GeostrophicWind                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.abl_pressure_grad**                      | imposed horizontal pressure gradient vector, used when | 3 Reals                        | 0.0 0.0 0.0            |
|                                                | ``abl_driver_type`` = ``PressureGradient``             |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.abl_geo_forcing**                        | geostrophic forcing vector applied directly, rather    | 3 Reals                        | 0.0 0.0 0.0            |
|                                                | than being computed from a geostrophic wind            |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.abl_geo_wind**                           | geostrophic wind vector from which the geostrophic     | 3 Reals                        | 0.0 0.0 0.0            |
|                                                | forcing is computed; read only when ``use_coriolis``   |                                |                        |
|                                                | is true and ``abl_driver_type`` = ``GeostrophicWind``, |                                |                        |
|                                                | and ignored if ``abl_geo_wind_table`` is given         |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.abl_geo_wind_table**                     | file containing a height-varying geostrophic wind      | String                         | None                   |
|                                                | profile; may not be combined with                      |                                |                        |
|                                                | ``add_custom_geostrophic_profile``                     |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.const_massflux_u**                       | target mass flux in the x direction; requires periodic | Real                           | 0.0 (not used)         |
|                                                | boundaries in x                                        |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.const_massflux_v**                       | target mass flux in the y direction; requires periodic | Real                           | 0.0 (not used)         |
|                                                | boundaries in y                                        |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.const_massflux_tau**                     | relaxation time scale used by the constant-mass-flux   | Real > 0                       | 1.0                    |
|                                                | forcing                                                |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.const_massflux_layer_lo**                | height of the bottom of the layer over which the       | Real                           | bottom of domain       |
|                                                | constant mass flux is enforced                         |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.const_massflux_layer_hi**                | height of the top of the layer over which the constant | Real                           | top of domain          |
|                                                | mass flux is enforced                                  |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+

Custom, Nudging and Numerical Diffusion Forcing
-----------------------------------------------

+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| Parameter                                      | Definition                                             | Acceptable Values              | Default                |
+================================================+========================================================+================================+========================+
| **erf.add_custom_rhotheta_forcing**            | add the problem-specific (rho*theta) source term       | Boolean                        | false                  |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.add_custom_moisture_forcing**            | add the problem-specific moisture source term          | Boolean                        | false                  |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.add_custom_w_subsidence**                | add the problem-specific vertical subsidence velocity  | Boolean                        | false                  |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.add_do_theta_advection**                 | apply the custom subsidence to the (rho*theta)         | Boolean                        | true                   |
|                                                | equation; only used with ``add_custom_w_subsidence``   |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.add_do_mom_advection**                   | apply the custom subsidence to the momentum equations; | Boolean                        | true                   |
|                                                | only used with ``add_custom_w_subsidence``             |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.add_custom_geostrophic_profile**         | use the problem-specific height-varying geostrophic    | Boolean                        | false                  |
|                                                | wind profile; may not be combined with                 |                                |                        |
|                                                | ``abl_geo_wind_table``                                 |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.custom_forcing_uses_primitive_vars**     | interpret the custom forcing terms as sources for the  | Boolean                        | false                  |
|                                                | primitive variables rather than the conserved          |                                |                        |
|                                                | variables                                              |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.spatial_rhotheta_forcing**               | the custom (rho*theta) forcing varies horizontally as  | Boolean                        | false                  |
|                                                | well as vertically                                     |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.spatial_moisture_forcing**               | the custom moisture forcing varies horizontally as     | Boolean                        | false                  |
|                                                | well as vertically                                     |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.nudging_from_input_sounding**            | nudge the solution towards the time-varying profiles   | Boolean                        | false                  |
|                                                | given in the input sounding files                      |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.nudging_q_z1**                           | height below which moisture is not nudged              | Real                           | 0.0                    |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.nudging_q_z2**                           | height above which moisture is not nudged              | Real                           | 10000.0                |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.nudging_t_z1**                           | height below which temperature is not nudged           | Real                           | 0.0                    |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.nudging_t_z2**                           | height above which temperature is not nudged           | Real                           | 10000.0                |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.large_scale_forcing**                    | apply the large-scale forcing read from a forcing file | Boolean                        | false                  |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.num_diff_coeff**                         | coefficient of the sixth-order numerical diffusion; 0  | Real in [0,1]                  | 0.0                    |
|                                                | turns it off                                           |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+

Immersed Forcing and Canopy Source Terms
----------------------------------------

See :ref:`Forcings` for the formulation of the immersed forcing terms and for
guidance on choosing the drag coefficients, which have different defaults for
the compressible and anelastic solvers.

+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| Parameter                                      | Definition                                             | Acceptable Values              | Default                |
+================================================+========================================================+================================+========================+
| **erf.immersed_forcing_substep**               | apply the immersed forcing source terms during the     | Boolean                        | true if compressible,  |
|                                                | acoustic substeps only                                 |                                | false if anelastic     |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.forest_substep**                         | apply the forest canopy source terms during the        | Boolean                        | false                  |
|                                                | acoustic substeps only                                 |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.if_Cd_momentum**                         | immersed forcing drag coefficient for momentum         | Real > 0                       | 500.0 if compressible, |
|                                                |                                                        |                                | 50.0 if anelastic      |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.if_Cd_scalar**                           | immersed forcing drag coefficient for scalars          | Real > 0                       | 50.0 if compressible,  |
|                                                |                                                        |                                | 5.0 if anelastic       |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.if_implicit_drag**                       | use a point-implicit (linearly implicit) form of the   | Boolean                        | false                  |
|                                                | immersed forcing drag rather than an explicit          |                                |                        |
|                                                | forward-Euler source                                   |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.if_z0**                                  | roughness length [m] used by the immersed forcing wall | Real > 0                       | 0.1                    |
|                                                | model                                                  |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.if_surf_temp_flux**                      | surface temperature flux [K m/s] imposed at immersed   | Real                           | 1.e-8                  |
|                                                | surfaces; only one of ``if_surf_temp_flux``,           |                                |                        |
|                                                | ``if_init_surf_temp`` and ``if_Olen`` may be set       |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.if_init_surf_temp**                      | surface temperature [K] imposed at immersed surfaces;  | Real                           | 0.0                    |
|                                                | only one of ``if_surf_temp_flux``,                     |                                |                        |
|                                                | ``if_init_surf_temp`` and ``if_Olen`` may be set       |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.if_surf_heating_rate**                   | rate of change [K/hr] of the immersed surface          | Real                           | 0.0                    |
|                                                | temperature (converted internally to K/s)              |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.if_Olen**                                | Obukhov length [m] imposed at immersed surfaces; only  | Real                           | 1.e-8                  |
|                                                | one of ``if_surf_temp_flux``, ``if_init_surf_temp``    |                                |                        |
|                                                | and ``if_Olen`` may be set                             |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.if_use_most**                            | use the Monin-Obukhov similarity theory wall model at  | Boolean                        | false                  |
|                                                | immersed surfaces                                      |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.if_stability_correction**                | include the stability corrections in the immersed      | Boolean                        | false                  |
|                                                | forcing similarity functions; use with caution for     |                                |                        |
|                                                | horizontal walls                                       |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.if_ws_floor**                            | lower bound [m/s] on the wind speed used by the        | Real > 0                       | 0.001                  |
|                                                | immersed forcing wall model                            |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.if_damp_alpha**                          | damping coefficient used in the immersed forcing wall  | Real                           | 0.5                    |
|                                                | model                                                  |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.use_rotate_surface_flux**                | rotate the MOST surface stresses to be normal to the   | Boolean                        | false                  |
|                                                | terrain surface; requires that terrain be used         |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+

Lateral Boundary Nudging for Real-Data Runs
-------------------------------------------

+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| Parameter                                      | Definition                                             | Acceptable Values              | Default                |
+================================================+========================================================+================================+========================+
| **erf.bdy_nudge_factor**                       | strength of the Davies relaxation in the lateral       | Real > 0                       | 10.0                   |
|                                                | boundary region; the nudging rate is                   |                                |                        |
|                                                | 1/(``bdy_nudge_factor`` * dt)                          |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.bdy_rho_nudge_factor**                   | Davies factor used for density specifically; a         | Real                           | -1.0 (use              |
|                                                | non-positive value means ``bdy_nudge_factor`` is used  |                                | ``bdy_nudge_factor``)  |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.use_wrf_bdy_density**                    | use the dry-air density reconstructed from the         | Boolean                        | true                   |
|                                                | ``wrfbdy`` data.  Ignored (set to false) unless real   |                                |                        |
|                                                | boundary conditions are used with ``init_type`` =      |                                |                        |
|                                                | ``WRFInput``                                           |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.use_wrf_bdy_qc_qi**                      | ingest cloud water and cloud ice from ``wrfinput`` and | Boolean                        | false                  |
|                                                | ``wrfbdy``; requires an active moisture model and real |                                |                        |
|                                                | boundary conditions with ``init_type`` = ``WRFInput``  |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.bdy_moist_nudge_type**                   | which approach is used to nudge the moist variables in | 0, 1, 2, 3                     | 1                      |
|                                                | the boundary region; 3 requires ``use_wrf_bdy_qc_qi``  |                                |                        |
|                                                | = true                                                 |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+

Inflow Perturbation
-------------------

+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| Parameter                                      | Definition                                             | Acceptable Values              | Default                |
+================================================+========================================================+================================+========================+
| **erf.perturbation_type**                      | type of inflow turbulence generation used at a level   | None, Source, Direct, CPM,     | None                   |
|                                                | (per-level).  See                                      | CPM_W                          |                        |
|                                                | :ref:`sec:InflowTurbulenceGeneration`                  |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+

Wind Farm Models
----------------

See :ref:`sec:WindFarmModels` for a description of the wind farm
parametrizations and the format of the tables below.

+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| Parameter                                      | Definition                                             | Acceptable Values              | Default                |
+================================================+========================================================+================================+========================+
| **erf.windfarm_type**                          | which wind farm parametrization to use                 | None, Fitch, EWP, SimpleAD,    | None                   |
|                                                |                                                        | GeneralAD                      |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.windfarm_loc_type**                      | coordinate system in which the turbine locations are   | None, lat_lon, x_y             | None                   |
|                                                | given                                                  |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.windfarm_loc_table**                     | file listing the turbine locations                     | String                         | None                   |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.windfarm_spec_table**                    | file containing the turbine specifications             | String                         | None                   |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.windfarm_spec_table_extra**              | file containing additional turbine specifications      | String                         | None                   |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.windfarm_blade_table**                   | file containing the blade geometry, used by the        | String                         | None                   |
|                                                | generalized actuator disk                              |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.windfarm_airfoil_tables**                | directory containing the airfoil tables, used by the   | String                         | None                   |
|                                                | generalized actuator disk                              |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.sampling_distance_by_D**                 | distance upstream of the turbine, as a multiple of the | Real > 0                       | -1.0 (must be set)     |
|                                                | rotor diameter, at which the incoming free stream      |                                |                        |
|                                                | velocity is sampled; required with ``windfarm_type`` = |                                |                        |
|                                                | ``SimpleAD``                                           |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.turb_disk_angle_from_x**                 | angle [degrees] of the face of the turbine disk from   | Real > 0                       | -1.0 (must be set)     |
|                                                | the x axis; a turbine facing flow in the x direction   |                                |                        |
|                                                | has a value of 90.  Required with ``windfarm_type`` =  |                                |                        |
|                                                | ``SimpleAD`` or ``GeneralAD``                          |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.windfarm_x_shift**                       | distance by which the bounding box of the wind farm is | Real > 0                       | -1.0 (must be set)     |
|                                                | shifted from the x axis; required with                 |                                |                        |
|                                                | ``windfarm_loc_type`` = ``lat_lon``                    |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.windfarm_y_shift**                       | distance by which the bounding box of the wind farm is | Real > 0                       | -1.0 (must be set)     |
|                                                | shifted from the y axis; required with                 |                                |                        |
|                                                | ``windfarm_loc_type`` = ``lat_lon``                    |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+

Hindcast Forcing and Hurricane Tracking
---------------------------------------

See :ref:`sec:HindCast` for a description of the hindcast forcing capability
and the expected format of the boundary and surface data.

+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| Parameter                                      | Definition                                             | Acceptable Values              | Default                |
+================================================+========================================================+================================+========================+
| **erf.hindcast_lateral_forcing**               | drive the lateral boundaries from hindcast data        | Boolean                        | false                  |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.hindcast_boundary_data_dir**             | directory containing the hindcast lateral boundary     | String                         | None (must be set)     |
|                                                | data; required if ``hindcast_lateral_forcing`` is true |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.hindcast_data_interval_in_hrs**          | time interval [hr] between successive hindcast data    | Real > 0                       | -1.0 (must be set)     |
|                                                | files; required if ``hindcast_lateral_forcing`` is     |                                |                        |
|                                                | true                                                   |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.hindcast_lateral_sponge_strength**       | strength of the lateral sponge used with hindcast      | Real > 0                       | -1.0 (must be set)     |
|                                                | forcing; required if ``hindcast_lateral_forcing`` is   |                                |                        |
|                                                | true                                                   |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.hindcast_lateral_sponge_length**         | width of the lateral sponge region used with hindcast  | Real > 0                       | -1.0 (must be set)     |
|                                                | forcing; required if ``hindcast_lateral_forcing`` is   |                                |                        |
|                                                | true                                                   |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.hindcast_zhi_sponge_damping**            | apply sponge damping near the upper boundary with      | Boolean                        | false                  |
|                                                | hindcast forcing                                       |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.hindcast_zhi_sponge_length**             | depth of the upper sponge region; required if          | Real > 0                       | -1.0 (must be set)     |
|                                                | ``hindcast_zhi_sponge_damping`` is true                |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.hindcast_zhi_sponge_strength**           | strength of the upper sponge; required if              | Real > 0                       | -1.0 (must be set)     |
|                                                | ``hindcast_zhi_sponge_damping`` is true                |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.hindcast_surface_bcs**                   | set the surface boundary conditions from hindcast data | Boolean                        | false                  |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.hindcast_surface_data_dir**              | directory containing the hindcast surface data; read   | String                         | None                   |
|                                                | only if ``hindcast_surface_bcs`` is true               |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.io_hurricane_eye_tracker**               | write out files tracking the eye of a hurricane        | Boolean                        | false                  |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.hurricane_eye_latitude**                 | approximate latitude of the hurricane eye in the       | Real                           | must be set            |
|                                                | initial condition; required if                         |                                |                        |
|                                                | ``io_hurricane_eye_tracker`` is true                   |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.hurricane_eye_longitude**                | approximate longitude of the hurricane eye in the      | Real                           | must be set            |
|                                                | initial condition; required if                         |                                |                        |
|                                                | ``io_hurricane_eye_tracker`` is true                   |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+

Diagnostics and Testing
-----------------------

+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| Parameter                                      | Definition                                             | Acceptable Values              | Default                |
+================================================+========================================================+================================+========================+
| **erf.Ave_Plane**                              | index of the direction normal to the planes used when  | 0, 1, 2                        | 2                      |
|                                                | computing horizontal averages                          |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.time_avg_vel**                           | accumulate and output time-averaged velocity fields    | Boolean                        | false                  |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.test_mapfactor**                         | set the map scale factors to 0.5 instead of 1 for      | Boolean                        | false                  |
|                                                | testing; see :ref:`sec:MapFactors`                     |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+

Ensemble Initialization
-----------------------

``erf.is_init_for_ensemble`` is read here, and when it is true the remaining
``ensemble.*`` parameters listed in :ref:`sec:EnsembleInitialization` are read
as well.

Deprecated and Removed Inputs
-----------------------------

``SolverChoice`` explicitly checks for the following inputs so that runs using
older inputs files fail with an informative message rather than silently
ignoring them.

+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| Parameter                                      | Definition                                             | Acceptable Values              | Default                |
+================================================+========================================================+================================+========================+
| **erf.use_terrain**                            | removed; set ``erf.terrain_type`` instead              | none -- input removed          | aborts if present      |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.use_moist_background**                   | removed; set ``erf.init_type`` = ``MoistBaseState``    | none -- input removed          | aborts if present      |
|                                                | instead                                                |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.moisture_type**                          | renamed; set ``erf.moisture_model`` instead            | none -- input renamed          | aborts if present      |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.no_substepping**                         | removed; set ``erf.substepping_type`` instead          | none -- input removed          | aborts if present      |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.init_type**                              | the values ``Real`` and ``Ideal`` were replaced by     | ``WRFInput`` in place of       | aborts if ``Real`` or  |
|                                                | ``WRFInput``                                           | either retired value           | ``Ideal`` is used      |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.terrain_type**                           | the values ``Static`` and ``Moving`` were replaced by  | ``StaticFittedMesh``,          | warns if ``Static`` or |
|                                                | ``StaticFittedMesh`` and ``MovingFittedMesh``; the old | ``MovingFittedMesh``           | ``Moving`` is used     |
|                                                | names still work                                       |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.land_surface_model**                     | the value ``OceanSurf`` was removed; set               | ``None``, ``SLM``, ``NOAHMP``  | aborts if              |
|                                                | ``erf.use_coupled_sst`` = 1 and use                    |                                | ``OceanSurf`` is used  |
|                                                | ``erf.land_surface_model`` for the land surface itself |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+

.. _sec:DampingChoice:

Damping Options in ``DampingChoice``
====================================

The parameters in this section are the complete set of inputs read by
``DampingChoice`` in ``Source/DataStructs/ERF_DampingStruct.H``, which holds the
Rayleigh damping and vertical-velocity damping options.  All of them use the
``erf.`` prefix.

Rayleigh Damping
----------------

Rayleigh damping relaxes the solution towards a reference state within a layer
of depth ``erf.rayleigh_zdamp`` below the model top, using a sine-squared ramp
that is zero at the bottom of that layer.  It is off unless at least one of the
``erf.rayleigh_damp_*`` flags is set.  See :ref:`Forcings` for the form of the
damping terms.

+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| Parameter                                      | Definition                                             | Acceptable Values              | Default                |
+================================================+========================================================+================================+========================+
| **erf.rayleigh_damping_type**                  | how the Rayleigh damping terms are integrated in time: | SlowExplicit, FastExplicit,    | SlowExplicit           |
|                                                | ``SlowExplicit`` applies them explicitly in the slow   | FastImplicit                   |                        |
|                                                | (Runge-Kutta) update, while ``FastExplicit`` and       |                                |                        |
|                                                | ``FastImplicit`` apply them explicitly or implicitly   |                                |                        |
|                                                | within the acoustic substeps.  Case sensitive; an      |                                |                        |
|                                                | unrecognized value aborts                              |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.rayleigh_damp_U**                        | apply Rayleigh damping in the x-momentum equation      | Boolean                        | false                  |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.rayleigh_damp_V**                        | apply Rayleigh damping in the y-momentum equation      | Boolean                        | false                  |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.rayleigh_damp_W**                        | apply Rayleigh damping in the z-momentum equation      | Boolean                        | false                  |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.rayleigh_damp_T**                        | apply Rayleigh damping in the potential temperature    | Boolean                        | false                  |
|                                                | equation                                               |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.rayleigh_dampcoef**                      | inverse damping time scale [1/s] multiplying the       | Real > 0                       | 0.2                    |
|                                                | vertical damping weight                                |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.rayleigh_zdamp**                         | depth [m] of the damping layer below the model top     | Real > 0                       | 500.0                  |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+

Vertical-Velocity Damping
-------------------------

Vertical-velocity damping adds a source term that opposes vertical motion
wherever the vertical advective Courant number exceeds ``erf.w_damping_cfl``.
If ``erf.w_damping`` is true, at least one of ``erf.w_damping_const`` or
``erf.w_damping_coeff`` must be given a positive value, otherwise the run
aborts; if both are positive, ``erf.w_damping_const`` is used.  See
:ref:`sec:NumericalStability` for the two damping formulations.

+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| Parameter                                      | Definition                                             | Acceptable Values              | Default                |
+================================================+========================================================+================================+========================+
| **erf.w_damping**                              | enable vertical-velocity damping                       | Boolean                        | false                  |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.w_damping_cfl**                          | critical vertical advective Courant number above which | Real > 0                       | 1.0                    |
|                                                | the damping is applied                                 |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.w_damping_const**                        | constant damping coefficient [m/s^2], as used by WRF   | Real                           | -1.0 (not used)        |
|                                                | (which uses 0.3); a positive value selects this        |                                |                        |
|                                                | formulation                                            |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| **erf.w_damping_coeff**                        | dimensionless damping factor giving a grid-dependent   | Real                           | -1.0 (not used)        |
|                                                | coefficient proportional to dz/dt^2; a positive value  |                                |                        |
|                                                | selects this formulation                               |                                |                        |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+

Deprecated Damping Inputs
-------------------------

+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
| Parameter                                      | Definition                                             | Acceptable Values              | Default                |
+================================================+========================================================+================================+========================+
| **erf.rayleigh_damp_substep**                  | removed; set ``erf.rayleigh_damping_type`` instead     | none -- input removed          | aborts if present      |
+------------------------------------------------+--------------------------------------------------------+--------------------------------+------------------------+
