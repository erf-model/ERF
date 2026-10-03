
 .. role:: cpp(code)
    :language: c++

.. _subsec:LinearSolvers:

Linear Solvers
==============

Evolving the anelastic equation set requires solution of a Poisson equation in which we solve for the update to the perturbational pressure at cell centers.
ERF uses several solver options available through AMReX: geometric multigrid, Fast Fourier Transforms (FFTs) and preconditioned GMRES.
For simulations without terrain-fitted coordinates or grid stretching, one of the FFT options is generally the fastest solver,
followed by multigrid.  We note that the multigrid solver has the option to "ignore" a coordinate direction
if the domain is only one cell wide in that direction; this allows for efficient solution of effectively 2D problems.
Multigrid can also be used when the union of grids at a level is not in itself rectangular; the FFT solvers do not work in that general case.

For simulations using grid stretching in the vertical but a flat lower boundary, we must use the hybrid FFT solver in which
we perform 2D transforms only in the lateral directions and couple the solution in the vertical direction with a tridiagonal solve.
In both these cases we use a 7-point stencil.

To solve the Poisson equation using terrain-fitted coordinates with general terrain, the stencil has variable
coefficients and couples each cell to its neighbours in the vertical planes through the metric cross terms
(15 points: the 7-point part scaled by the face areas and the Jacobian, plus the :math:`x`-:math:`z` and
:math:`y`-:math:`z` corner terms through :math:`h_\xi` and :math:`h_\eta`).
Three solvers are available for it, chosen with ``erf.terrain_poisson_solver``:

* ``gmres_fft`` (the default): GMRES preconditioned by the hybrid FFT solve of the flat problem.
  It requires the FFT build and a rectangular union of grids on the level.
* ``mlmg``: geometric multigrid on the full terrain stencil (the ``MLTerrainPoisson`` operator below).
  It works in every build and on any union of grids, and is the only option for a refined level whose
  boxes do not form a rectangle.
* ``gmres_mlmg``: GMRES on the multigrid operator (``amrex::GMRESMLMG``) with
  ``erf.terrain_mlmg_precond_iters`` V-cycles as the preconditioner.  The V-cycles use a fixed number of
  smoothing sweeps as their bottom solve, so that the preconditioner is a fixed linear operator: a bottom
  solve iterated to a tolerance is not one, and GMRES then stops on a residual estimate that has drifted
  from the true residual (observed on the ridge case below: a reported convergence with a divergence of
  1e-5 left).  The GMRES driver of ``gmres_fft`` can take the same V-cycles through
  ``TerrainPoisson::setPrecondFunction``, which is unit-tested and left for a cheaper preconditioner.

All three solve the same discrete system to the tolerances ``erf.poisson_reltol`` / ``erf.poisson_abstol``
and compute the face fluxes of the velocity correction from the same stencil, so the projected velocity is
discretely divergence-free for the terrain operator whichever solver is used.

.. note::

   The FFT solver / preconditioner can only be used when the union of grids at a level is itself rectangular.

Multigrid for terrain-fitted coordinates
----------------------------------------

``MLTerrainPoisson`` (``Source/LinearSolvers/ERF_MLTerrainPoisson.H``) is an ``amrex::MLCellLinOp`` whose
matrix-vector product is the terrain stencil of ``ERF_TerrainPoisson_3D_K.H``, the one the GMRES operator applies,
so the two are identical to the bit on the finest level.  It differs from the standard AMReX cell-centred operators
in the following ways.

**Coarse levels.**  Every multigrid level carries its own nodal terrain surface, obtained by keeping every other node
of the finer level in each coarsened direction (the coarse surface is exactly the fine surface at the coarse nodes),
and the face areas :math:`a_x, a_y, a_z` and the Jacobian :math:`J` recomputed from it with the same routines ERF uses on
the finest level.  The coarse operators therefore discretize the same terrain on a coarser mesh rather than averaging
the fine metric coefficients.  Ghost nodes of a coarse surface outside the domain follow ERF's rules for the finest
level: the nearest node inside the domain laterally, linear extrapolation below the surface and above the top.
When map factors scale the face areas, the scale is carried to the coarse levels by subsampling.

**Boundary conditions.**  Periodic, homogeneous Neumann (every boundary type except ``Outflow`` and ``Open``) and
homogeneous Dirichlet (``Outflow``/``Open``, filled at second order with ``ghost = -interior`` as in the GMRES operator).
On a refined level the boundary of the union of boxes is a coarse/fine boundary where the momentum is fixed by the
coarser level, so the correction flux through it is set to zero exactly, as through the surface.  (A Neumann mirror
of the ghost cells is not enough there: on sloping terrain it leaves the cross-term part of the flux, and an all-Neumann
problem on the union is then inconsistent and no solver converges on it.)  Because the stencil reads the edge ghost
cells in the :math:`x`-:math:`z` and :math:`y`-:math:`z` planes, the operator fills its ghost cells itself from a
precomputed source map, on the cells no box of the level covers.

**Smoother.**  Red-black (by the parity of the column index :math:`i+j`) line relaxation in :math:`z`: every column
segment of a box is relaxed with the exact tridiagonal block of the operator.  The stencil never couples cells whose
:math:`i` and :math:`j` both differ, so the sweep is an exact block Gauss-Seidel.  The tridiagonal coefficients,
including every boundary fold, are read off the stencil itself by applying it to the unit vector of each cell as the
ghost fill sees it, so they cannot drift from the operator.

**Singular problems.**  The operator is :math:`J^{-1}` times a conservative divergence, so its left null vector is
:math:`J` rather than the constant: when every boundary is periodic or Neumann the solvability offset subtracted from
the right-hand side is the :math:`J`-weighted mean, not the plain mean.

**Hierarchy.**  The multigrid coarsens every direction by two as long as all of them allow it and then the remaining
ones alone (AMReX semicoarsening), so a domain with few vertical cells still gets a deep hierarchy and a small
bottom problem; the column smoother suits the flat cells this produces.  AMReX can only coarsen as far as every box
of the level allows: a box 150 cells wide allows one coarsening although the 600-cell domain it tiles allows three,
and the two-level hierarchy that results leaves a bottom problem too large for the bottom solver (the ridge case
below stalled at a relative residual of 1e-5 that way).  When the level's boxes allow fewer coarsenings than the
domain in some direction, the solve is therefore re-gridded onto boxes whose extents are multiples of the full
coarsening ratio, the metrics and right-hand side are copied over and the solution and fluxes copied back; the
operator does not depend on the layout, so only the solver's round-off changes.  ``erf.mg_v = 1`` reports when this
happens.  The operator, with its coarse metrics, ghost source maps and column coefficients, is built once per level
and layout and kept between projections of a static terrain (it is rebuilt when the level is remade, and never kept
for a moving terrain).

**Cost.**  Measured on one node (two MPI ranks, Release build, tolerances 1e-8), per projection: on the 2D
Witch-of-Agnesi ridge of the RANS suite (128 x 64 cells) the FFT-preconditioned GMRES takes 7 iterations,
the multigrid 14 V-cycles and the multigrid-preconditioned GMRES 12 iterations, at 1 s, 3 s and 3 s for 50 steps;
on a 2D ridge of slope 0.63 at 4 m (600 x 224 cells, Inflow/Outflow, terrain smoothing) 19 GMRES iterations against
7 V-cycles and 13 preconditioned GMRES iterations, at 52 s, 95 s and 140 s for 100 steps; on Askervein
(300 x 300 x 18 cells, stretched) 7 GMRES iterations against 14 V-cycles and 43 preconditioned GMRES iterations, at
46 s, 390 s and 780 s for 20 steps.
The multigrid takes fewer iterations but each V-cycle costs about nine applications of the 15-point stencil (the
column relaxation re-evaluates the operator in every sweep: 60 % of the solve time on the ridge, the operator
applications of the bottom solve another 25 %), while the FFT preconditioner is nearly exact for a terrain-following
mesh and needs some twenty applications in all; on the ridge a projection costs 0.11 s with the FFT, 0.31 s with the
multigrid and 0.64 s with the multigrid-preconditioned GMRES.  Keeping the operator between projections changes none of
this (its setup was below 5 % of the time), and one smoothing sweep instead of two needs twice the cycles and is slower.
A domain whose extents have few factors of two (300 = 4 x 75, 18 = 2 x 9) also gives a shallow hierarchy with a large
bottom problem.  The FFT-preconditioned GMRES therefore stays the default; the multigrid options are for builds without
FFT, for refined levels whose boxes do not form a rectangle, and as a reference solver.  Precomputing the metric
coefficients of the stencil per level (instead of recomputing the node differences and their ratios in every
evaluation) is the remaining lever and would roughly halve the smoother cost.

A one-cell-wide lateral direction is hidden from the multigrid as in the flat solve; such a direction must be
periodic or Neumann.  Embedded-boundary terrain is solved with the EB multigrid solver and does not use this operator.
Unit tests (``Tests/Unit/LinearSolvers/ERF_GTestMLTerrainPoisson.cpp``) check the bitwise agreement with the GMRES
operator, the flux/divergence consistency, second-order convergence against a manufactured solution over a hill,
the coarse-level surfaces and metrics, the probed tridiagonal coefficients against the applied operator, the smoother
and the solver; the regression tests ``TerrainSolver_Hill2D`` / ``TerrainSolver_Hill3D`` compare the three solvers
on the RANS hill cases and ``TerrainMLMG_LShape`` runs a refined level with an L-shaped union of boxes.
