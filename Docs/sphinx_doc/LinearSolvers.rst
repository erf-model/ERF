
 .. role:: cpp(code)
    :language: c++

.. _subsec:LinearSolvers:

Linear Solvers
==============

Evolving the anelastic equation set requires solution of a Poisson equation in which we solve for the update to the perturbational pressure at cell centers.   AMReX provides several solver options: geometric multigrid, Fast Fourier Transforms (FFTs) and preconditioned GMRES.  For simulations with no terrain or grid stretching, one of the FFT options is generally the fastest solver, followed by multigrid.  We note that the multigrid solver has the option to ``ignore'' a coordinate direction if the domain is only one cell wide in that direction; this allows for efficient solution of effectively 2D problems.  Multigrid can also be used when the union of grids at a level is not in itself rectangular; the FFT solvers do not work in that general case.

For simulations using grid stretching in the vertical but flat terrain, we must use the hybrid FFT solver in which we perform 2D transforms only in the lateral directions and couple the solution in the vertical direction with a tridiagonal solve. In both these cases we use a 7-point stencil. To solve the Poisson equation on terrain-fitted coordinates with general terrain, we rely on the preconditioned GMRES solver since the stencil effectively has variable coefficients and requires 19 points.
