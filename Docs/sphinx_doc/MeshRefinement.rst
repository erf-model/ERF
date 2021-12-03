
 .. role:: cpp(code)
    :language: c++

 .. role:: fortran(code)
    :language: fortran

 .. _MeshRefinement:

Mesh Refinement
===============

Note: Mesh refinement is a WIP -- for now the documentation is aspirational.

Grid Creation
-------------

ERF allows both static and dynamic mesh refinement.  For the static refinement, we currently control
the placement of grids using

+--------------------------+------------------+-----------------+-------------+
| Parameter                | Definition       | Acceptable      | Default     |
|                          |                  | Values          |             |
+==========================+==================+=================+=============+
| **tagging.tag_region**   | are we going to  | true / false    | false       |
|                          | specify a tagged |                 |             |
|                          | region for       |                 |             |
|                          | refinement       |                 |             |
+--------------------------+------------------+-----------------+-------------+
| **tagging.region_lo**    | low corner of    | 3 Reals         | None        |
|                          | physical         |                 |             |
|                          | location for     |                 |             |
|                          | refinement if    |                 |             |
|                          | tag_region true  |                 |             |
+--------------------------+------------------+-----------------+-------------+
| **tagging.region_hi**    | high corner of   | 3 Reals         | None        |
|                          | physical         |                 |             |
|                          | location for     |                 |             |
|                          | refinement if    |                 |             |
|                          | tag_region true  |                 |             |
+--------------------------+------------------+-----------------+-------------+

Note that the tagged region will be covered by one or more boxes.  The user may
specify the region to be covered but not the decompostion of the region into
individual grids.

See the `Gridding`_ section of the AMReX documentation for details of how individual grids are created.

.. _`Gridding`: https://amrex-codes.github.io/amrex/docs_html/ManagingGridHierarchy_Chapter.html

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
