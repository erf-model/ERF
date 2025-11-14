.. _sec:build:overview:

=================================
1.0 Introduction to the ERF Build Ecosystem
=================================

**1.1 ERF Capabilities**

The build system scripts are designed exclusively for compiling the ERF source code. They do not contain descriptions of the software's scientific capabilities or underlying architecture.

**1.2 AMReX Framework Dependency**

ERF is built upon the AMReX framework and can be configured in two ways. The default method, controlled by the ``ERF_USE_INTERNAL_AMREX`` option, utilizes an internal AMReX submodule located at ``Submodules/AMReX``. Alternatively, the build system can be configured to link against a pre-installed, external version of AMReX.

**1.3 Physics Packages and Requirements**

ERF supports several advanced physics packages, including RRTMGP, SHOC, and P3. Enabling any of these packages automatically enables the EKAT library dependency (``ERF_ENABLE_EKAT``). This is because RRTMGP, SHOC, and P3 use the Kokkos performance portability framework for heterogeneous computing; EKAT is the library that provides Kokkos, which is accessed as a submodule within EKAT. The build system enforces specific prerequisites for these components: EKAT requires an MPI implementation to be enabled, and the RRTMGP radiation package requires NetCDF libraries to be available.

**1.4 Compiler Standards**

The ERF build system requires a C++ compiler that is compliant with the C++17 standard. In addition, as ERF relies on the AMReX framework, it also inherits AMReX's requirements for a Fortran 2003 and a C99 compliant compiler.

**1.5 Platform Support**

ERF is fully supported and regularly tested on Linux-based operating systems, with a particular focus on Department of Energy (DOE) High-Performance Computing (HPC) systems such as Perlmutter (NVIDIA A100), Frontier (AMD MI250X), Polaris (NVIDIA A100), and Aurora (Intel GPUs). Support for these systems is streamlined through machine-specific configuration profiles located in the ``Build/machines/`` directory, which are detailed later in this document. While many users successfully build on macOS, it is not officially supported. Windows is not an officially supported platform. While continuous integration tests are performed, the development team does not provide user support for this platform.

=================================
2.0 Quick Start
=================================

**2.1 Prerequisites**

To clone the repository and build the software, a standard development environment is required. This includes Git for source control, Python (version 2.7 or newer), and a standard C++ toolchain including a compiler (e.g., GCC) and make.

=================================
3.0 Build System Choices
=================================

**3.1 GNU Make Description**

The GNU Make build system provides fine-grained control over the compilation process. It uses a series of ``Make.package`` files and configuration variables (e.g., ``USE_RRTMGP``, ``USE_NETCDF``) within the main ``Exec/Make.ERF`` file. These variables conditionally add source files to the build and pass preprocessor definitions (e.g., ``-DERF_USE_RRTMGP``) to the compiler.

**3.2 CMake Description**

The CMake build system is designed for cross-platform compatibility and robust dependency management. It uses a series of ``option()`` commands in its ``CMakeLists.txt`` file (e.g., ``option(ERF_ENABLE_MPI "Enable MPI" OFF)``) to generate user-configurable cache variables. These variables control which features are included in the build and how dependencies are linked.

**3.3 Machine-Specific Profiles**

For HPC environments, ERF provides machine-specific profiles located in the ``Build/machines/`` directory that configure the shell environment by managing modules and environment variables. Using these profiles is the canonical practice for ensuring reproducible builds. Examples for major DOE systems include ``perlmutter_erf.profile`` (NERSC Perlmutter), ``frontier_erf.profile`` (OLCF Frontier), ``polaris_erf.profile`` (ALCF Polaris), and ``aurora_erf.profile`` (ALCF Aurora).

=================================
4.0 System Requirements
=================================

.. list-table:: Component Requirements
   :header-rows: 1
   :widths: 20 60 15 5

   * - Component
     - Description
     - Default
     - Possible Values
   * - Compilers
     - A C++17 compliant compiler (e.g., GCC >= 8.0), a Fortran 2003 compiler, and a C99 compiler are required. Supported vendors include GNU, Cray, and Intel.
     - None
     - gnu/intel/cray
   * - Build Tools
     - CMake version 3.14 or newer is required to build ERF, with version 3.25 or newer recommended for full functionality on Cray systems. Additionally, Python version 2.7 or newer is necessary for the build process.
     - None
     - Various
   * - Parallelism (Optional)
     - For distributed memory parallelism, an MPI implementation (such as MPICH or OpenMPI) is required when the ``ERF_ENABLE_MPI`` option is turned ON.
     - OFF
     - ON/OFF
   * - GPU Support (Optional)
     - ERF supports GPU acceleration via multiple backends, including NVIDIA CUDA, AMD ROCm/HIP, and Intel oneAPI (for SYCL). When using CUDA, Toolkit version 11.0 or newer is required.
     - OFF
     - ON/OFF
   * - I/O Libraries (Optional)
     - Parallel HDF5 and NetCDF libraries are used for file I/O operations. If the ``ERF_ENABLE_NETCDF`` option is enabled, the build system will automatically enable ``ERF_ENABLE_HDF5`` by default.
     - OFF
     - ON/OFF

=================================
5.0 Support Statement
=================================

**5.1 Fully Supported**

ERF is fully supported and actively tested on Linux systems. This includes dedicated support for major DOE HPC platforms, including Perlmutter (NERSC), Frontier (OLCF), Polaris (ALCF), and Aurora (ALCF).

**5.2 Partial Support**

macOS is considered partially supported. While many users and developers successfully build and run ERF on this platform, the development team does not have the resources to provide full, dedicated support for Mac users.

**5.3 Not Supported**

Windows is not an officially supported platform. While continuous integration tests are performed, the development team does not provide user support for this platform. Windows executables using MPI are provided by the github workflow at ``.github/workflows/windows-mpi.yml`` and can be found in the actions tab for that workflow: https://github.com/erf-model/ERF/actions/workflows/windows-mpi.yml

=================================
6.0 Notes and Additional Resources
=================================

**6.1 WarpX/AMReX Documentation Acknowledgment**

ERF's documentation, particularly for machine profiles and HPC configurations, is inspired by best practices from the WarpX and AMReX projects.

**6.2 Documentation Index**

The complete build documentation consists of the following sections:

.. toctree::
   :maxdepth: 2

   buildingSystems
   buildingHPC
   buildingLibraryConfig
   buildingTroubleshooting

* :ref:`sec:build:systems` - A detailed technical reference comparing the GNU Make and CMake build systems, including key variables, workflows, and the rationale for choosing one over the other.
* :ref:`sec:build:hpc` - An architectural overview of ERF's strategy for building on HPC platforms, with a deep dive into the automated CMake detection system for Cray environments.
* :ref:`sec:build:library` - A guide to configuring external libraries such as NetCDF, HDF5, and the E3SM physics packages (SHOC, P3, RRTMGP), with a focus on parallel I/O requirements.
* :ref:`sec:build:troubleshooting` - A collection of solutions for common build failures, explanations of historical issues resolved by build system automation, and strategies for debugging complex configuration problems.

**6.3 AMReX Framework Documentation**

Official documentation for the underlying AMReX framework is available at https://amrex-codes.github.io/amrex/docs_html/, which is essential for understanding core build system variables (``AMREX_HOME``) and capabilities.
