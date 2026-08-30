
 .. role:: cpp(code)
    :language: c++

 .. _MeshRefinement:

Mesh Refinement
===============

ERF also allows both dynamic and static mesh refinement with sub-cycling in time at finer levels of refinement.
Arbitrary integer refinement ratios are allowed although typically ratios of 2, 3 or 4 are used; refinement can also be anisotropic,
allowing refinement in one coordinate direction but not another.

We utilize two-way coupling, in which the coarse solution is used to provide boundary conditions for the fine solution
and the fine solution is averaged down onto the coarser level.  In addition, we reflux all advected scalars to ensure conservation.
For coarse-to-fine communication we provide "ghost cell" data for cell-centered data and tangential momentum components to
the fine level by interpolating in space and time outside the region covered by the fine level.
We also interpolate the normal momentum from coarse to fine on the coarse-fine interface itself;
this ensures mass conservation since the normal momentum is in fact the flux for the density field.
In order to ensure that the fine momentum on the coarse-fine boundary stays consistent with the
interpolated coarse values throughout a fine timestep, we also interpolate the source term for
the normal momentum on the coarse-fine interface.
When using the anelastic approximation, this ensures that the computation of the updates to the
fine momentum do not use any values of the perturbational pressure from the coarser level since
the perturbational pressure is not synchronized between levels.

Note that any tagged region will be covered by one or more boxes.  The user may
specify the refinement criteria and/or region to be covered, but not the decomposition of the region into
individual grids.

See the `Gridding`_ section of the AMReX documentation for details of how individual grids are created.

By default the boxes covering a tagged region span the full vertical extent of
that region, i.e. the grids are not decomposed in the z direction; see
:ref:`subsec:no-vertical-decomposition` for the parameters that control this.

.. _`Gridding`: https://amrex-codes.github.io/amrex/docs_html/ManagingGridHierarchy_Chapter.html

.. note::

   Three different things are easily conflated when talking about the vertical
   direction, and they are controlled by three different sets of parameters.

   **Refining in z** means changing the grid spacing: the finer level has
   :math:`dz_{fine} = dz_{coarse} / r_z`, where :math:`r_z` is the refinement
   ratio in the z direction.  Creating a level greater than 0 over some region
   does *not* by itself imply that dz changes there.  We nevertheless say that
   such a region is "refined", so the word is used for both ideas; whether the
   vertical spacing actually changes is a separate question, decided only by the
   refinement ratio.  ``amr.ref_ratio`` sets one ratio per level which is then
   applied in all three directions, so using it always refines in z as well.  To
   create finer levels without changing dz, use ``amr.ref_ratio_vect`` with 1 in
   the third slot, for example

   ::

          amr.max_level = 1
          amr.ref_ratio_vect = 2 2 1

   which refines by a factor of 2 horizontally while leaving the vertical grid
   spacing of the fine level equal to that of the coarse level.  Initialization
   from WRF input files supports either choice: a wrfinput file for a nest carries
   the same eta levels as its parent, so if the refinement ratio in z is greater
   than 1 then ERF interpolates the file onto the finer vertical grid as it reads
   it.  See :ref:`sec:nested-wrfinput`.

   **The vertical extent of the refined region** -- how much of the depth of the
   domain the finer level covers -- is a property of the region being tagged, not
   of the refinement ratio.  It is set by the refinement box for static
   refinement (see :ref:`subsec:full-depth-refinement`) and by the tagging
   criterion, optionally buffered in z, for dynamic refinement (see
   :ref:`subsec:dynamic-full-depth`).

   **The vertical decomposition of that region into individual grids** -- whether
   the region is chopped in z into several boxes so that they can be distributed
   across processors -- is yet another question, controlled by
   ``amr.max_grid_size_z`` and ``amr.refine_grid_layout_z``; see
   :ref:`subsec:no-vertical-decomposition`.

   A level can therefore be created with the same dz as its parent, spanning the
   full depth of the domain, in boxes that are not split vertically -- and each of
   those three statements is arranged by a different parameter.

Static Mesh Refinement
----------------------

For static refinement, we can control the placement of grids by specifying
the low and high extents (in physical space or index space) of each box.

The following example demonstrates how to tag regions for static refinement.
In this first example, all cells in the region ((1200,1400,0.)(3000,2400,2048.))
and in the region ((3200,3400,0.)(4000,4000,2048.)) are tagged for
one level of refinement.

::

          amr.max_level = 1
          amr.ref_ratio = 2

          erf.refinement_indicators = box1 box2

          erf.box1.in_box_lo = 1200 1400 0
          erf.box1.in_box_hi = 3000 2400 2048
          erf.box1.max_level = 1

          erf.box2.in_box_lo = 3200 3400 0
          erf.box2.in_box_hi = 4000 4000 2048
          erf.box2.max_level = 1

It is possible to refine the box by multiple levels, and to have different levels of refinement for different boxes. For example, if we set

::

          amr.max_level = 2
          amr.ref_ratio = 2 2

          erf.refinement_indicators = box1 box2

          erf.box1.in_box_lo = 1200 1400 0
          erf.box1.in_box_hi = 3000 2400 2048
          erf.box1.max_level = 1

          erf.box2.in_box_lo = 3200 3400 0
          erf.box2.in_box_hi = 4000 4000 2048
          erf.box2.max_level = 2

box1 is refined up to level 1, but box2 is refined up to level 2. If no max_level is specified for a
refinement indicator, the region is refined up to the maximum level specified in ``amr.max_level``.
Error messages are generated to help users adjust the refinements to create a valid box at each level.


We note that instead of specifying the physical extent enclosed, we can instead specify the indices of
the bounding box of the refined region in the index space of that fine level.
To do this we use
``in_box_lo_indices`` and ``in_box_hi_indices`` instead of ``in_box_lo`` and ``in_box_hi``.
If we want to refine the inner region (spanning half the width in each direction) by one level of
factor 2 refinement, and the domain has 32x64x8 cells at level 0 covering the domain, then we would set

::

          amr.max_level = 1
          amr.ref_ratio = 2

          erf.refinement_indicators = box1

          erf.box1.in_box_lo_indices = 16 32  4
          erf.box1.in_box_hi_indices = 47 95 11
          erf.box1.max_level = 1

There is also an option to specify the indices of the bounding box of the refined region in the index space of the coarser level, using
``in_box_lo_indices_crse`` and ``in_box_hi_indices_crse``.  This is useful when the user has a particular region in mind that they want to refine,
and they know the indices of that region on the coarser level but not on the finer level.  In this case, the code will automatically adjust the
indices to create a valid box at the finer level.

::

          amr.max_level = 1
          amr.ref_ratio = 2

          erf.refinement_indicators = box1

          erf.box1.in_box_lo_indices_crse = 16 32  4
          erf.box1.in_box_hi_indices_crse = 47 95 11
          erf.box1.max_level = 1


The lo_indices should be divisible by the refinement ratio, and the hi_indices should be one less than a number divisible by the refinement ratio.
There are no such requirements for the coarse level indices, since the code will adjust them as needed to create a valid box at the finer level.

.. _subsec:full-depth-refinement:

Refining the Full Depth of the Domain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is very common to want a refined region that covers the full depth of the
domain, either because a physics package needs entire columns or simply because
there is no reason to limit the refinement vertically.  Rather than spelling out
the vertical extent, the user may specify **only the x and y values** for any of
the three forms of the box specification, in which case ERF fills in the vertical
extent with the full extent of the domain:

-  for ``in_box_lo`` / ``in_box_hi``, the vertical extent is set to
   ``geometry.prob_lo`` and ``geometry.prob_hi`` in z

-  for ``in_box_lo_indices`` / ``in_box_hi_indices``, it is set to the full range
   of cells in z in the domain at the level being created

-  for ``in_box_lo_indices_crse`` / ``in_box_hi_indices_crse``, it is set to the
   full range of cells in z in the domain at the next coarser level

So the first example above could equivalently have been written as

::

          amr.max_level = 1
          amr.ref_ratio = 2

          erf.refinement_indicators = box1 box2

          erf.box1.in_box_lo = 1200 1400
          erf.box1.in_box_hi = 3000 2400
          erf.box1.max_level = 1

          erf.box2.in_box_lo = 3200 3400
          erf.box2.in_box_hi = 4000 4000
          erf.box2.max_level = 1

if the domain runs from 0 to 2048 in z.  This form has the advantage that the
inputs file does not have to be edited if the vertical extent of the domain or
the number of cells in z is changed.  Each of ``in_box_lo`` and ``in_box_hi``
(and their index-space counterparts) should be given either two or three values,
and the lo and hi specifications for a given indicator must have the same number
of values.

Some options **require** a refinement region that spans the full depth of the
domain.  The PBL models -- MYJ, MYNN2.5, MYNN-EDMF, YSU and MRF -- compute the
height of the boundary layer by working up each column, so every box must go from
the bottom of the domain to the top.  If a refinement box is specified that does
not, ERF prints a message telling the user how to correct the inputs file, then
aborts.  Using the two-value form above is the simplest way to guarantee this.
Similarly, the SHOC PBL model requires that no box be split in the vertical
direction; see :ref:`subsec:no-vertical-decomposition`.

A refined region that deliberately covers only part of the depth is common in
terrain runs; see :ref:`subsec:partial-depth-refinement`.

The two-value form applies to statically specified boxes.  For dynamically
created grids there is a separate, blunter mechanism based on the vertical
tagging buffer; see :ref:`subsec:dynamic-full-depth`.

If the vertical extent *is* given explicitly, note that the number of cells in the
z direction should be divisible by the refinement ratio in the vertical direction,
and that the specified indices are snapped outward to the nearest indices aligned
with the refinement ratio (a message is printed whenever this snapping changes the
box).

.. _subsec:partial-depth-refinement:

Refining Only Part of the Depth with Terrain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A nested LES region rarely needs to be refined all the way to the model top, so a
refinement box whose vertical extent stops well below ``geometry.prob_hi`` in z is a
natural thing to want.  This is supported with terrain-fitted coordinates for both of the
initialization pathways that build a stratified base state:

-  ``erf.init_type = input_sounding``, with the terrain read from a text file via
   ``erf.terrain_file_name``; and

-  ``erf.init_type = WRFInput``, whether each level has its own wrfinput file or only
   level 0 does.

For example, with a domain 17500 m deep::

          amr.max_level      = 1
          amr.ref_ratio_vect = 3 3 1

          erf.refinement_indicators = box1
          erf.box1.max_level = 1
          erf.box1.in_box_lo =  70000.   90000.     0.
          erf.box1.in_box_hi = 110000.  150000.  6000.

Two things to watch, both specific to how the box is specified:

-  All three components must be given.  The two-value form described above fills the
   vertical extent with the full depth of the domain, which is the opposite of what is
   wanted here.

-  The z values are located on the nominal (undeformed) vertical grid -- that is, on
   ``zlevels_stag`` -- and not on the terrain-following heights, so the box may not reach
   the height you expect.  This matters most with ``erf.init_type = WRFInput``, where the
   nominal grid is uniform while the heights the wrfinput file supplies are strongly
   compressed near the surface.  Check the box ERF reports
   (``Saving in 'boxes at level'``) against the ``z_phys`` field in the plotfile.

With ``erf.init_type = WRFInput``, if a file is given for the finer level then the box
must be contained in the region that file covers; ERF checks this and aborts with both
boxes printed if it is not.

The base state is built independently on each level, against that level's own heights,
and its vertical integration is causal upward -- so a level that stops partway up gets
exactly the base state it would have had if it had been refined to the model top.  See
:ref:`subsec:base-state-refined` for that, and :ref:`sec:nested-wrfinput` for the
wrfinput-specific details.

The restrictions are the same for both pathways: the PBL models and the column-integral
derived quantities need entire columns and cannot be used on a level whose grids do not
span the domain in z.  See :ref:`sec:nested-wrfinput`.

Moving Refinement Regions
~~~~~~~~~~~~~~~~~~~~~~~~~

We can also adapt this static refinement paradigm to specify rectangular regions whose locations
are a prescribed function of time.   Following the WRF paradigm for moving nested grids, we specify
motion in terms of the speed of the grid in separate time intervals.  Modifying the first example above,
we could have the following, where the first box would be static but the second box would change in time.
In particular, the box would move east at 1 m/s for the first 10 minutes, be at rest from 10 minutes to 20 minutes,
then move north at 1 m/2 from 20 minutes to 80 minutes after the start time.

::

          amr.max_level = 1
          amr.ref_ratio = 2

          erf.refinement_indicators = box1 box2

          erf.box1.in_box_lo = 1200 1400 0
          erf.box1.in_box_hi = 3000 2400 2048
          erf.box1.max_level = 1

          erf.box2.in_box_lo = 3200 3400 0
          erf.box2.in_box_hi = 4000 4000 2048
          erf.box2.max_level = 1
          erf.box2.move_start_time = 0.   1200.
          erf.box2.move_stop_time  = 600. 4800.
          erf.box2.move_speed_x    = 1. 0.
          erf.box2.move_speed_y    = 0. 1.

Dynamic Mesh Refinement
-----------------------

Dynamically created tagging functions are based on runtime data specified in the inputs file.
These dynamically generated functions test on either state variables or derived variables
defined in ERF_derive.cpp and included in the derive_lst in Setup.cpp.
(We note that static refinement can also be achieved by using the refinement criteria as specified below
but setting ``erf.regrid_int`` to a number greater than the total number of steps that will be taken.)

Available tests include

-  “greater\_than”: :math:`field >= threshold`

-  “less\_than”: :math:`field <= threshold`

-  “adjacent\_difference\_greater”: :math:`max( | \text{difference between any nearest-neighbor cell} | ) >= threshold`

This example adds three user-named criteria –
hi\_rho: cells with density greater than 1 on level 0, and greater than 2 on level 1 and higher;
lo\_theta: cells with theta less than 300 that are inside the region ((.25,.25,.25)(.75,.75,.75));
and adv_diff: cells having a difference in the scalar of 0.01 or more from that of any immediate neighbor.
The first will trigger up to AMR level 3, the second only to level 1, and the third to level 2.
The third will be active only when the problem time is between 0.001 and 0.002 seconds.

Note that density and rhoadv_0 are the names of state variables, whereas theta is the name of a derived variable,
computed by dividing the variable named rhotheta by the variable named density.

::

          erf.refinement_indicators = hi_rho lo_theta advdiff

          erf.hi_rho.max_level = 3
          erf.hi_rho.value_greater = 1. 2.
          erf.hi_rho.field_name = density

          erf.lo_theta.max_level = 1
          erf.lo_theta.value_less = 300
          erf.lo_theta.field_name = rhotheta
          erf.lo_theta.in_box_lo = .25 .25 .25
          erf.lo_theta.in_box_hi = .75 .75 .75

          erf.advdiff.max_level = 2
          erf.advdiff.adjacent_difference_greater = 0.01
          erf.advdiff.field_name = rhoadv_0
          erf.advdiff.start_time = 0.001
          erf.advdiff.end_time = 0.002

.. _subsec:dynamic-full-depth:

Full-Depth Refinement with Dynamic Tagging
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Note that when ``in_box_lo`` / ``in_box_hi`` are combined with a field-based
criterion, as in ``lo_theta`` above, the box only restricts the region in which
cells may be tagged; it does not force the resulting fine grids to fill that
region.  Left to itself, the vertical extent of a dynamically created grid simply
follows the cells that satisfy the criterion.

The full depth of the domain can nonetheless be enforced by using the vertical
buffer ``amr.n_error_buf_z``.  Before the grids are generated, the set of tagged
cells is grown by ``amr.n_error_buf`` cells in each direction, and that buffer
may be set per direction.  If the vertical buffer is at least as large as the
number of cells in the z direction at the level being tagged, then every tagged
cell is grown into a full column and the resulting boxes reach from the bottom of
the domain to the top.  For a domain with 64 cells in the vertical, for example,

::

          amr.n_error_buf_x = 2
          amr.n_error_buf_y = 2
          amr.n_error_buf_z = 64

This is not the most efficient way to obtain full-depth refinement -- the tag
arrays are allocated with ``n_error_buf`` ghost cells and the buffering is redone
at every regrid, so a large vertical buffer costs both memory and time -- but it
requires no advance knowledge of where the criterion will trigger, which is
exactly the situation in which dynamic refinement is used.

This technique applies only to dynamic refinement: ERF aborts with
``Don't use n_error_buf > 0 when setting the box explicitly`` if a nonzero
``n_error_buf`` is combined with an explicitly specified refinement box.  For
static refinement, use the two-value form described in
:ref:`subsec:full-depth-refinement` instead.
