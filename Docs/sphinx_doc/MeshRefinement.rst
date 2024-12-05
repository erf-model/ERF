
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
For coarse-to-fine communication we provide ``ghost cell'' data for cell-centered data and tangential momentum components to
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

.. _`Gridding`: https://amrex-codes.github.io/amrex/docs_html/ManagingGridHierarchy_Chapter.html

Static Mesh Refinement
----------------------

For static refinement, we can control the placement of grids by specifying
the low and high extents (in physical space or index space) of each box.

The following example demonstrates how to tag regions for static refinement.
In this first example, all cells in the region ((.15,.25,0.)(.35,.45,1.))
and in the region ((.65,.75,0.0)(.85,.95,1.0)) are tagged for
one level of refinement.

::

          amr.max_level = 1
          amr.ref_ratio = 2

          erf.refinement_indicators = box1 box2

          erf.box1.in_box_lo = .15 .25 0.0
          erf.box1.in_box_hi = .35 .45 1.0

          erf.box2.in_box_lo = .65 .75 0.0
          erf.box2.in_box_hi = .85 .95 1.0

In the example below, we refine the region ((.15,.25,0.)(.35,.45,.5))
by two levels of factor 3 refinement. In this case, the refined region at level 1 will
be sufficient to enclose the refined region at level 2.

::

          amr.max_level = 2
          amr.ref_ratio = 3 3

          erf.refinement_indicators = box1

          erf.box1.in_box_lo = .15 .25 0.0
          erf.box1.in_box_hi = .35 .45 1.0

And in this final example, the region ((.15,.25,0.)(.35,.45,1.))
will be refined by two levels of factor 3, but the larger region, ((.05,.05,0.)(.75,.75,1.))
will be refined by a single factor 3 refinement.

::

          amr.max_level = 2
          amr.ref_ratio = 3 3

          erf.refinement_indicators = box1 box2

          erf.box1.in_box_lo = .15 .25 0.0
          erf.box1.in_box_hi = .35 .45 1.0

          erf.box2.max_level = 1
          erf.box2.in_box_lo = .05 .05 0.0
          erf.box2.in_box_hi = .75 .75 1.0


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
