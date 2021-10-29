
 .. role:: cpp(code)
    :language: c++

 .. role:: fortran(code)
    :language: fortran

 .. _MeshRefinement:

Mesh Refinement
---------------

Note: Mesh refinement is a WIP -- for now the documentation is aspirational.

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
