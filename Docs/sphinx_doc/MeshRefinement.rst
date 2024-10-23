
 .. role:: cpp(code)
    :language: c++

 .. _MeshRefinement:

Mesh Refinement
===============

ERF allows both static and dynamic mesh refinement, as well as the choice of one-way or two-way coupling.
Dynamic refinement is currently only allowed when terrain is not being used.

The refinement ratio is specified by the user at runtime. Refinement is not allowed in the vertical.

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

Coupling Types
--------------

ERF supports one-way and two-way coupling between levels; this is a run-time input

::

      erf.coupling_type = "OneWay" or "TwoWay"

By one-way coupling, we mean that between each pair of refinement levels,
the coarse level communicates data to the fine level to serve as boundary conditions
for the time advance of the fine solution. For cell-centered quantities,
and face-baced normal momenta on the coarse-fine interface, the coarse data is conservatively
interpolated to the fine level.

The interpolated data is utilized to specify ghost cell data (outside of the valid fine region)
as well as specified data inside the lateral boundaries of the fine region.
See :ref:`sec:LateralBoundaryConditions` for the details of how the relaxation works; when
used in the context of mesh refinement we fill the specified values by interpolation from the
coarser level rather than reading from the external file. For coarse/fine boundaries,
a user may specify the total width of the interior specified (Dirichlet) and relaxation region with
``erf.cf_width = <Int>`` (yellow + blue)
and analogously the width of the interior specified (Dirichlet) region may be specified with
``erf.cf_set_width = <Int>`` (yellow).

Setting ``erf.cf_set_width = 0`` designates that we interpolate the momenta
at faces only on the coarse-fine boundary itself; no interior cell-centered data, or momenta
inside the fine region, are filled from the coarser level.

By two-way coupling, we mean that in additional to interpolating data from the coarser level
to supply boundary conditions for the fine regions,
the fine level also communicates data back to the coarse level in two ways:

- The fine cell-centered data are conservatively averaged onto the coarse mesh covered by fine mesh.

- The fine momenta are conservatively averaged onto the coarse faces covered by fine mesh.

- A "reflux" operation is performed for all cell-centered data; this updates values on the coarser
  level outside of regions covered by the finer level.

We note that when one-way coupling is used, quantities which are advanced in conservation form
potentially violate global conservation.  Two-way coupling ensures conservation of mass, and of the advective contribution
to all scalar updates, but does not account for loss of conservation due to diffusive or source terms.
