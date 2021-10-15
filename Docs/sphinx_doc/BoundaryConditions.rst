
 .. role:: cpp(code)
    :language: c++

 .. role:: fortran(code)
    :language: fortran

 .. _BCs:

Boundary Conditions
-------------------

ERF manages boundary conditions in a form consistent with many AMReX codes. Ghost cell data are updated over an AMR level during a ``FillPatch`` operation and fluxes are then computed over the entire box without specifically recognizing boundary cells. The Fortran routine ``pc_hypfill`` in ``bc_fill_nd.F90`` is called to set state data at physical boundaries for this purpose.  A generic boundary filler function, ``filcc_nd``, is supplied to fill standard boundary condition types that do not require user input, including:

* *Interior* - Copy-in-intersect in index space (same as periodic boundary conditions). Periodic boundaries are set in the ERF inputs file
* *Symmetry* - All conserved quantities and the tangential momentum component are reflected from interior cells without
  sign change (REFLECT_EVEN) while the normal component is reflected with a sign change (REFLECT_ODD)
* *NoSlipWall* - REFLECT_EVEN is applied to all conserved quantities except for both tangential and normal momentum components which are updated
  using REFLECT_ODD
* *SlipWall*  - SlipWall is identical to Symmetry
* *FOExtrap* - First-order extrapolation: the value in the ghost-cells are a copy of the last interior cell.

More complex boundary conditions require user input that is prescribed explicitly.  In the code, all types are formally handled in ``pc_hypfill``; ``filcc_nd`` is called first to handle all the above types.  Boundaries identified as ``UserBC`` in the inputs will be tagged as ``EXT_DIR`` in ``pc_hypfill`` and will be ignored by ``filcc_nd``.  Users will then fill the Dirichlet boundary values, typically by calling the helper function, ``bcnormal``. The indirection here is not required, but is recommended for reasons discussed below.

