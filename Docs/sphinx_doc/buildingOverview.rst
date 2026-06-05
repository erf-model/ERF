.. _sec:build:overview:

Build Requirements and Dependencies
====================================

This page provides technical reference information about ERF's build requirements, dependencies, and platform support. For step-by-step build instructions, see :ref:`sec:build:systems`.

Dependencies
------------

ERF depends on the following external libraries:

* **AMReX** - Adaptive mesh refinement framework (required, provided as submodule)
* **RTE-RRTMGP** - Radiation model from E3SM (optional, provided as submodule)
* **SHOC** - Simplified Higher-Order Closure turbulence and cloud macrophysics from E3SM (optional, requires extra setup)
* **P3** - Microphysics scheme from E3SM (optional, requires extra setup)
* **NetCDF** - I/O library for reading WRF inputs and writing plotfiles (optional, system-provided)
* **HDF5** - Parallel I/O backend for NetCDF (optional, system-provided)

**AMReX Configuration**

ERF uses the AMReX framework. The default method uses an internal AMReX submodule located at ``Submodules/AMReX``. Alternatively, the build system can link against a pre-installed, external version of AMReX.

**Submodule Management**

AMReX, EKAT, NOAH-MP, and RTE-RRTMGP are available as submodules in the ERF repository. SHOC and P3 require additional setup steps (see :ref:`sec:build:systems`). Kokkos is accessed as a submodule within EKAT.

**Dependency Relationships**

Physics packages have specific dependency requirements:

.. graphviz::
   :align: center
   :caption: Physics Package Dependencies

   digraph physics_deps {
       rankdir=TB;
       node [shape=box, style="rounded,filled"];
       edge [fontsize=10];

       // Physics packages
       RRTMGP [label="RRTMGP\n(Radiation)", fillcolor=lightcoral];
       SHOC [label="SHOC\n(Turbulence)", fillcolor=lightcoral];
       P3 [label="P3\n(Microphysics)", fillcolor=lightcoral];

       // Dependencies
       NetCDF [fillcolor=lightblue];
       HDF5 [label="HDF5", fillcolor=lightblue];
       Kokkos [label="Kokkos\n(from EKAT)", fillcolor=wheat];
       MPI [label="MPI\n(required)", fillcolor=lightyellow];

       // Dependency arrows
       RRTMGP -> NetCDF [label="enforced"];
       RRTMGP -> Kokkos [label="auto-enabled"];
       RRTMGP -> MPI [label="required*"];

       SHOC -> Kokkos [label="auto-enabled"];
       SHOC -> MPI [label="required*"];

       P3 -> Kokkos [label="auto-enabled"];
       P3 -> MPI [label="required*"];

       NetCDF -> HDF5 [style=dashed, label="parallel I/O"];

       // Layout
       {rank=same; RRTMGP SHOC P3}
       {rank=same; NetCDF Kokkos MPI}
       {rank=min; HDF5}
   }

.. note::
   MPI is required by EKAT (which provides Kokkos). CMake enforces this; GNU Make assumes it's set.

Enabling RRTMGP, SHOC, or P3 automatically enables EKAT, which provides the Kokkos performance portability framework. The build system enforces these prerequisites:

* EKAT requires MPI to be enabled
* RRTMGP requires both NetCDF and MPI to be enabled

Compiler and Tool Requirements
-------------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 60 20

   * - Component
     - Description
     - Version
   * - **C++ Compiler**
     - Must support C++17 standard
     - GCC >= 8.0
   * - **C Compiler**
     - Must support C99 standard
     - Any modern version
   * - **Fortran Compiler (Optional)**
     - Must support Fortran 2003 standard
     - Any modern version
   * - **CMake**
     - Build configuration tool
     - >= 3.14 (>= 3.25 for Cray)
   * - **Python**
     - Build script dependency
     - >= 2.7 (including 3.x)
   * - **MPI (Optional)**
     - Distributed memory parallelism
     - MPICH, OpenMPI, Cray MPICH
   * - **CUDA (Optional)**
     - NVIDIA GPU support
     - >= 11.0
   * - **HIP (Optional)**
     - AMD GPU support
     - Latest versions only
   * - **SYCL (Optional)**
     - Intel GPU support
     - Latest versions only
   * - **NetCDF (Optional)**
     - I/O library (parallel-enabled)
     - >= 4.6
   * - **HDF5 (Optional)**
     - I/O backend (parallel-enabled)
     - >= 1.10

.. note::
   When ``ERF_ENABLE_NETCDF`` is enabled, the build system automatically enables ``ERF_ENABLE_HDF5``. Both libraries must be parallel-enabled versions when building with MPI.

Compiler-Specific Notes
-----------------------

**Intel Compilers**

.. warning::
   ERF compiles successfully with recent Intel compiler suites (e.g., icx version 2024.1.0). However, older versions may require reduced compiler optimization to avoid internal compiler errors.

   For example, with icpc version 19.1.2.254, most files compile with ``-O2`` (``CMAKE_BUILD_TYPE = RelWithDebInfo``), but ``TimeIntegration/ERF_advance_dycore.cpp`` may need manual compilation with ``-O1``.

**SYCL Support**

.. note::
   While ERF is designed to work with SYCL for Intel GPUs, we do not guarantee it will build and run on all Intel platforms. SYCL support is under active development.

Platform Support
----------------

**Fully Supported**

ERF is fully supported and actively tested on Linux systems, including major DOE HPC platforms:

* **Perlmutter (NERSC)** - NVIDIA A100 GPUs
* **Kestrel (NREL)** - NVIDIA H100 GPUs

**Partial Support**

Build setups are available for the following which use automatic configuration:

* **Frontier (OLCF)** - AMD MI250X GPUs
* **Polaris (ALCF)** - NVIDIA A100 GPUs
* **Aurora (ALCF)** - Intel GPUs

For platform-specific build instructions, see :ref:`sec:build:hpc`.

macOS is partially supported. Many users successfully build and run ERF on macOS, but the development team does not provide dedicated support.

**Not Supported**

Windows is not officially supported. Continuous integration tests are performed, but user support is not provided. Windows executables are available through GitHub Actions: https://github.com/erf-model/ERF/actions/workflows/windows-mpi.yml

Build System Options
--------------------

ERF supports two build systems. For detailed instructions, see :ref:`sec:build:systems`.

* **GNU Make** - Direct control over compiler flags and build variables
* **CMake** - Automated configuration and dependency detection

Additional Resources
--------------------

* **AMReX Documentation** - https://amrex-codes.github.io/amrex/docs_html/
* **WarpX HPC Documentation** - Inspired machine profile and some documentation structure https://warpx.readthedocs.io/en/latest/install/hpc.html
