
 .. role:: cpp(code)
    :language: c++

.. _subsec:AMReX:

AMReX
==============

ERF is built on AMReX, a C++--based software framework that supports the development of structured mesh algorithms for solving systems of partial differential equations, with options for adaptive mesh refinement, on machines from laptops to exascale architectures. AMReX was developed in the U.S. Department of Energy’s Exascale Computing Project and is now a member project of the High Performance Software Foundation under the umbrella of the Linux Foundation.

AMReX uses an MPI+X model of hierarchical parallelism where blocks of data are distributed across MPI ranks (typically across multiple nodes).  Fine-grained parallelism at the node level (X) is achieved using
OpenMP for CPU-only machines, or CUDA, HIP or SYCL for NVIDIA, AMD or Intel GPUs, respectively. AMReX provides extensive support for kernel launching on GPU accelerators (using ParallelFor looping constructs and C++ lambda functions) and the effective use of various memory types, including managed, device, and pinned. Architecture-specific aspects of the software for GPUs are highly localized within the code, and essentially hidden from the application developer or user.

In addition to portability across architectures, AMReX provides data structures and iterators that define, allocate and efficiently operate on distributed multi-dimensional arrays.
Data is defined on disjoint logically rectangular regions of the domain known as patches (or grids or boxes); we note that unlike WRF, AMReX (and therefore ERF) does not require one patch per MPI rank, thus allowing much more general domain decomposition. Common operations, such as parallel communication and reduction operations, as well as interpolation and averaging operators between levels of refinement, are provided by the AMReX framework.

Finally, ERF currently leverages, or plans to leverage, AMReX capabilities for effective load balancing, adaptive mesh refinement, memory management, asynchronous I/O, Lagrangian particles with particle-mesh interactions, and linear solvers.

We note that ERF supports both a cmake and a gmake build system.


Linear Solvers
==============

Evolving the anelastic equation set requires solution of a Poisson equation in which we solve for the update to the perturbational pressure at cell centers.   AMReX provides several solver options: geometric multigrid, Fast Fourier Transforms (FFTs) and preconditioned GMRES.  For simulations with no terrain or grid stretching, one of the FFT options is generally the fastest solver, followed by multigrid.  We note that the multigrid solver has the option to ``ignore'' a coordinate direction if the domain is only one cell wide in that direction; this allows for efficient solution of effectively 2D problems.  Multigrid can also be used when the union of grids at a level is not in itself rectangular; the FFT solvers do not work in that general case.

For simulations using grid stretching in the vertical but flat terrain, we must use the hybrid FFT solver in which we perform 2D transforms only in the lateral directions and couple the solution in the vertical direction with a tridiagonal solve. In both these cases we use a 7-point stencil. To solve the Poisson equation on terrain-fitted coordinates with general terrain, we rely on the preconditioned GMRES solver since the stencil effectively has variable coefficients and requires 19 points.
