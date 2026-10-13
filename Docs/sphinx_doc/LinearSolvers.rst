
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
Multigrid can also be used when the union of grids at a level is not in itself rectangular (except on terrain-fitted meshes); the FFT solvers do not work in that general case.

For simulations using grid stretching in the vertical but a flat lower boundary, we must use the hybrid FFT solver in which
we perform 2D transforms only in the lateral directions and couple the solution in the vertical direction with a tridiagonal solve.
In both these cases we use a 7-point stencil.

To solve the Poisson equation using terrain-fitted coordinates with general terrain,
we rely on the FFT-preconditioned GMRES solver since the stencil effectively has variable coefficients and requires 19 points.

.. note::

   The FFT solver / preconditioner can only be used when the union of grids at a level is itself rectangular.

.. _subsec:TerrainNeumannFaces:

Neumann faces on terrain-fitted meshes
--------------------------------------

On a terrain-fitted mesh the operator that FFT-preconditioned GMRES applies is a sum of
face fluxes divided by the Jacobian :math:`J = h_\zeta` of the cell,

.. math::

   (A\phi)_{i,j,k} = \frac{1}{J_{i,j,k}} \left[
       \frac{(a_x F_x)_{i+\frac12} - (a_x F_x)_{i-\frac12}}{\Delta x}
     + \frac{(a_y F_y)_{j+\frac12} - (a_y F_y)_{j-\frac12}}{\Delta y}
     + \frac{(a_z F_z)_{k+\frac12} - (a_z F_z)_{k-\frac12}}{\Delta \zeta} \right],

where :math:`a_x, a_y, a_z` are the face areas and the flux through an :math:`x`-face is the
horizontal gradient at constant height,

.. math::

   F_x = \frac{\partial \phi}{\partial \xi} - \frac{h_\xi}{h_\zeta} \frac{\partial \phi}{\partial \zeta},

with :math:`h_\xi` the slope of the coordinate surface (likewise :math:`F_y`; :math:`F_z`
holds both cross terms).

Some faces of the region being solved must not let any correction through: on a refined
level the coarse/fine faces, where the fine level takes the normal momentum from the level
below, and on any level the faces of the domain that are neither periodic nor Dirichlet
(``Outflow`` or ``Open``).  Those include walls, inflow, ``InflowOutflow``, ``HO_Outflow``
and symmetry faces, and, when real-data boundary conditions are used, every lateral face.
These are the Neumann faces.  ERF copies the solution into the ghost cell across such a
face (an even mirror).  That makes the first term of :math:`F_x` vanish there, but not the
second wherever the terrain slopes and :math:`\phi` varies in the vertical, so on its own the
mirror lets flux cross the face.

This matters because the problem on such a region is singular.  ERF makes it solvable by
removing the :math:`J`-weighted mean of the right-hand side, :math:`J` being the left null
vector of the operator.  Summing :math:`J A\phi` over the region, the interior fluxes cancel
in pairs and only the fluxes through the region's faces remain, so the mean removal is
right only if those are zero for every :math:`\phi`.  With the mirror alone they are not:
the system has no solution, GMRES leaves a large divergence after the solve, and a refined
level over sloping terrain drifts.  ERF used that operator before erf-model/ERF#4279.  The
operator now sets the whole flux through every Neumann face to zero, so the sum vanishes and
the solve converges to its tolerance; the projection then leaves the normal momentum on those
faces unchanged.

With map factors other than one (real-data runs on a map projection), the divergence that
ERF corrects carries :math:`m_x m_y`, and so does the change the map-factor scaled fluxes make
to it; ERF therefore solves :math:`A\phi = \text{rhs} / (m_x m_y)`, and on a singular problem
first removes the mean of the right-hand side weighted by :math:`J / (m_x m_y)` (its
solvability condition).  With map factors of one both steps leave the right-hand side as it is.

Results differ from ERF before erf-model/ERF#4279 only where the coordinate surfaces slope
across a Neumann face, at any height (the ground is one of them; a smoothed mesh can slope
aloft over flat ground).  Periodic faces and Dirichlet (``Outflow`` or ``Open``) faces are
treated as before, and so are symmetry faces, where ERF mirrors the mesh too, so that the
surfaces are flat across the face.  Where the surfaces are flat across every Neumann face the
results are the same bit for bit.  The top of the domain is a Neumann face when it is a wall
(``SlipWall`` or ``NoSlipWall``): the usual terrain-following meshes flatten the coordinate
surfaces to the top, so nothing changes there, but a mesh whose top surface is not flat (for
example one read from a real-data file) changes there too.  Flat terrain does not use
this operator at all: ERF treats the mesh as stretched and solves with the hybrid FFT
solver, which this does not change.
