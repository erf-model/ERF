
 .. role:: cpp(code)
    :language: c++

 .. role:: fortran(code)
    :language: fortran

 .. _MeshRefinement:

Mesh Refinement
===============

ERF allows both static and dynamic mesh refinement.

Note that any tagged region will be covered by one or more boxes.  The user may
specify the refinement criteria and/or region to be covered, but not the decomposition of the region into
individual grids.

See the `Gridding`_ section of the AMReX documentation for details of how individual grids are created.

.. _`Gridding`: https://amrex-codes.github.io/amrex/docs_html/ManagingGridHierarchy_Chapter.html

Static Mesh Refinement
----------------------

For static refinement, we control the placement of grids by specifying
the low and high extents (in physical space) of each box in the lateral
directions.   ERF enforces that all refinement spans the entire vertical direction.

The following example demonstrates how to tag regions for static refinement.
In this first example, all cells in the region ((.15,.25,prob_lo_z)(.35,.45,prob_hi_z))
and in the region ((.65,.75,prob_lo_z)(.85,.95,prob_hi_z)) are tagged for
one level of refinement, where prob_lo_z and prob_hi_z are the vertical extents of the domain:

::

          amr.max_level = 1
          amr.ref_ratio = 2

          amr.refinement_indicators = box1 box2

          amr.box1.in_box_lo = .15 .25
          amr.box1.in_box_hi = .35 .45

          amr.box2.in_box_lo = .65 .75
          amr.box2.in_box_hi = .85 .95

In the example below, we refine the region ((.15,.25,prob_lo_z)(.35,.45,prob_hi_z))
by two levels of factor 3 refinement. In this case, the refined region at level 1 will
be sufficient to enclose the refined region at level 2.

::

          amr.max_level = 2
          amr.ref_ratio = 3 3

          amr.refinement_indicators = box1

          amr.box1.in_box_lo = .15 .25
          amr.box1.in_box_hi = .35 .45

And in this final example, the region ((.15,.25,prob_lo_z)(.35,.45,prob_hi_z))
will be refined by two levels of factor 3, but the larger region, ((.05,.05,prob_lo_z)(.75,.75,prob_hi_z))
will be refined by a single factor 3 refinement.

::

          amr.max_level = 2
          amr.ref_ratio = 3 3

          amr.refinement_indicators = box1 box2

          amr.box1.in_box_lo = .15 .25
          amr.box1.in_box_hi = .35 .45

          amr.box2.in_box_lo = .05 .05
          amr.box2.in_box_hi = .75 .75
          amr.box2.max_level = 1


Dynamic Mesh Refinement
-----------------------

Dynamically created tagging functions are based on runtime data specified in the inputs file.
These dynamically generated functions test on either state variables or derived variables
defined in ERF_derive.cpp and included in the derive_lst in Setup.cpp.

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

          amr.refinement_indicators = hi_rho lo_theta advdiff

          amr.hi_rho.max_level = 3
          amr.hi_rho.value_greater = 1. 2.
          amr.hi_rho.field_name = density

          amr.lo_theta.max_level = 1
          amr.lo_theta.value_less = 300
          amr.lo_theta.field_name = rhotheta
          amr.lo_theta.in_box_lo = .25 .25 .25
          amr.lo_theta.in_box_hi = .75 .75 .75

          amr.advdiff.max_level = 2
          amr.advdiff.adjacent_difference_greater = 0.01
          amr.advdiff.field_name = rhoadv_0
          amr.advdiff.start_time = 0.001
          amr.advdiff.end_time = 0.002

Coupling Types
--------------

ERF supports both one-way and two-coupling; this is a run-time input

::

      erf.coupling_type = "OneWay" or "TwoWay"

By one-way coupling, we mean that between each pair of refinement levels,
the coarse mesh communicates Dirichlet data to the fine mesh in the form of ghost cell
data (outside of the valid fine region) for cell-centered quantities, and face-baced normal
momenta on the coarse-fine interface.  In both of these, the coarse data is conservatively
interpolated to the fine mesh.

By two-way coupling, we mean that in additional to the one-way coupling operations, the fine mesh
communicates data back to the coarse mesh in two ways:

- The fine cell-centered data is conservatively averaged onto the coarse mesh covered by fine mesh.

- A "reflux" operation is performed for all cell-centered data.

Because the normal momentum at the fine level is itself interpolated from the coarse, the
difference between fluxes -- along the coarse-fine interfaces -- used to update the coarse data and fluxes
used to update the fine data is due to the difference in the averaging of the advected quantity to the face
where the flux is defined.

We note that both coupling schemes are conservative for mass because the fluxes for the continuity
equation are the momenta themselves, which are interpolated on faces at the coarse-fine interface.  Other advected
quantities which are advanced in conservation form will lose conservation with one-way coupling.
Two-way coupling is conservative for these scalars as long as the refluxing operation is included with the
averaging down.
