.. _sec:build:systems:

=================================
ERF Build Systems and Options
=================================

The Energy Research and Forecasting (ERF) code, built upon the AMReX framework, possesses a complex dependency graph that includes multiple advanced physics packages (RTE-RRTMGP, SHOC, P3). To ensure performance portability across diverse computing architectures, ERF relies on the Kokkos library. This document provides a detailed, code-level explanation of ERF's two supported build systems, GNU Make and CMake, outlining how each system manages these dependencies to produce a functional executable.

The two build systems serve distinct purposes. GNU Make represents a traditional, application-centric system for compiling a single executable for a specific scientific run; it provides fine-grained debugging via utility targets like ``print-xxx``. CMake is the modern, recommended, framework-centric approach with policy-based automation, such as its Cray Detection System. It creates versioned libraries, manages dependencies through automated detection, and supports robust cross-platform builds, particularly on HPC systems.

This guide focuses exclusively on how the build systems are implemented and function, with direct references to the source code. It serves as a technical reference for developers and advanced users. Initial environment setup, dependency installation, and specific build examples are covered in other documentation sections.

1.0 Building with GNU Make
--------------------------

**1.1 System Overview and Use Case**

The GNU Make system provides a direct path to producing a single, case-specific executable. It uses an application-centric approach, where developers and scientists compile code for a particular problem setup (e.g., atmospheric boundary layer simulation) directly within the ``Exec/`` directory structure.

The primary use case for GNU Make system is scientific production runs and development scenarios where in-place compilation is preferred, avoiding overhead of creating and installing versioned libraries. The build is orchestrated by a ``GNUmakefile`` located within a specific case directory (such as ``Exec/ABL/``), which leverages build logic provided by the underlying AMReX framework.

**1.2 How it Works: The Orchestration Process**

The ERF GNU Make process uses a hierarchy of includes that separates user configuration from application and framework build logic. The user's ``GNUmakefile`` serves as the top layer, defining high-level build options. This file includes the application-level ``Make.ERF`` to gather all necessary ERF source files, which includes the core framework-level ``Make.defs`` and ``Make.rules`` from AMReX to execute compilation and linking.

The process follows these steps:

1. **GNUmakefile Location**: Process begins with user invoking ``make`` in application-specific directory, such as ``Exec/ABL/``. This directory contains the primary control file, ``GNUmakefile``.

2. **Set AMREX_HOME**: First step within ``GNUmakefile`` is defining ``AMREX_HOME`` variable. This points to location of AMReX submodule containing core build logic. Default path is ``$(ERF_HOME)/Submodules/AMReX``.

3. **Define Build Variables**: User defines boolean flags and other variables in ``GNUmakefile`` to control which features are compiled. These include options like ``USE_MPI``, ``USE_CUDA``, ``USE_RRTMGP``.

4. **Include ERF Sources**: ``GNUmakefile`` includes file ``Exec/Make.ERF``. This file systematically adds all necessary ERF source and header directories to build paths, populating ``VPATH_LOCATIONS`` and ``INCLUDE_LOCATIONS`` variables.

5. **Include AMReX Core Logic**: Finally, ``Make.ERF`` includes core AMReX makefiles, ``Make.defs`` and ``Make.rules``. These files contain generic logic for discovering dependencies, compiling source files, and linking object files into final executable.

**1.3 Important Make Variables**

The build is customized by setting variables within the application's ``GNUmakefile``. The following table details the most important variables.

.. list-table:: GNU Make Build Variables
   :header-rows: 1
   :widths: 20 50 15 15

   * - Variable Name
     - Description
     - Default
     - Possible Values
   * - ``AMREX_HOME``
     - Specifies path to AMReX source directory
     - ``$(ERF_HOME)/Submodules/AMReX``
     - Path string
   * - ``COMP``
     - Defines compiler suite to use
     - None
     - gnu/intel/cray
   * - ``DIM``
     - Sets dimensionality of problem (ERF is 3D only)
     - 3
     - 3
   * - ``USE_MPI``
     - Enables MPI for distributed-memory parallel execution
     - FALSE
     - TRUE/FALSE
   * - ``USE_OMP``
     - Enables OpenMP for shared-memory parallelism
     - FALSE
     - TRUE/FALSE
   * - ``USE_CUDA``
     - Enables NVIDIA GPU support via CUDA
     - FALSE
     - TRUE/FALSE
   * - ``USE_HIP``
     - Enables AMD GPU support via HIP
     - FALSE
     - TRUE/FALSE
   * - ``USE_SYCL``
     - Enables Intel GPU support via SYCL
     - FALSE
     - TRUE/FALSE
   * - ``USE_NETCDF``
     - Enables I/O support using NetCDF library
     - FALSE
     - TRUE/FALSE
   * - ``USE_NOAHMP``
     - Enables Noah-MP land surface model (requires ``USE_NETCDF=TRUE``)
     - FALSE
     - TRUE/FALSE
   * - ``USE_RRTMGP``
     - Enables RRTMGP radiation model (sets ``USE_KOKKOS=TRUE``, ``USE_NETCDF=TRUE``)
     - FALSE
     - TRUE/FALSE
   * - ``USE_SHOC``
     - Enables SHOC turbulence model (sets ``USE_KOKKOS=TRUE``)
     - FALSE
     - TRUE/FALSE
   * - ``USE_P3``
     - Enables P3 microphysics model (sets ``USE_KOKKOS=TRUE``)
     - FALSE
     - TRUE/FALSE
   * - ``USE_KOKKOS``
     - Enables Kokkos performance portability library
     - FALSE
     - TRUE/FALSE
   * - ``USE_MORR_FORT``
     - Enables Fortran-based Morrison microphysics scheme
     - FALSE
     - TRUE/FALSE
   * - ``USE_FFT``
     - Enables Fast Fourier Transform capabilities
     - FALSE
     - TRUE/FALSE

**1.4 Example GNUmakefile**

A typical ``GNUmakefile`` for an ERF application can be found at ``ERF/Exec/ABL/GNUmakefile``.

**1.5 Common Commands**

The AMReX make system provides standard utility targets for managing build environment and artifacts. For complete list of commands and advanced features, refer to `AMReX documentation <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html>`_. Common commands include:

* ``make`` - Compiles source code and creates executable in current directory
* ``make clean`` - Removes all build artifacts for all configurations
* ``make cleanconfig`` - Removes artifacts for current configuration only
* ``make print-xxx`` - Debugging tool that prints value of specified make variable

**1.6 Customization with Make.local**

For user- or site-specific customizations that should not be committed to main repository, AMReX GNU Make system provides ``Make.local`` file. This file, located at ``amrex/Tools/GNUMake/Make.local``, is the primary mechanism for overriding default build variables.

Settings defined in ``Make.local`` apply globally to all projects using that specific instance of AMReX submodule. Common use case is specifying non-default compiler version:

.. code:: shell

   CXX = g++-8
   CC  = gcc-8
   FC  = gfortran-8

2.0 Building with CMake
-----------------------

**2.1 System Overview and Rationale**

CMake is the recommended build system for ERF. It is a cross-platform tool with framework-centric philosophy. Instead of producing single executable for one use case, CMake produces versioned, exportable libraries. This design enables generation and installation of ``find_package``-compatible configuration file (``ERFConfig.cmake``), allowing other CMake-based projects to consume ERF as dependency with simple ``find_package(ERF)`` command.

Key advantages of using CMake for ERF:

* Promotes clean, out-of-source builds keeping compiled artifacts separate from source code
* Provides robust dependency detection through ``find_package`` mechanism
* Automatically detects and configures for Cray HPC systems, applying necessary workarounds

**2.2 Standard Workflow**

ERF offers several methods for invoking CMake build:

* **Simple Script**: Run ``Build/cmake.sh`` for basic, tested configuration
* **Advanced Script**: Use specialized scripts like ``Build/cmake_with_kokkos_many_cuda.sh`` for specific use cases
* **Manual Invocation**: Create build directory and invoke CMake: ``cmake -B build -S .``

Standard manual procedure from project root:

.. code:: shell

   mkdir build          # Create directory for build artifacts
   cd build            # Change to build directory
   cmake [options] ..  # Configure project
   make install        # Compile and install

**2.3 ERF CMake Options**

The ERF build is customized by passing options to cmake command using ``-D<VARIABLE>=<VALUE>`` syntax.

.. list-table:: ERF CMake Options
   :header-rows: 1
   :widths: 25 45 10 20

   * - Variable Name
     - Description
     - Default
     - Possible Values
   * - ``ERF_ENABLE_MPI``
     - Enables MPI support for parallel execution
     - OFF
     - ON/OFF
   * - ``ERF_ENABLE_OPENMP``
     - Enables OpenMP for shared-memory parallelism
     - OFF
     - ON/OFF
   * - ``ERF_ENABLE_CUDA``
     - Enables NVIDIA GPU support via CUDA (requires Toolkit >= 11.0)
     - OFF
     - ON/OFF
   * - ``ERF_ENABLE_HIP``
     - Enables AMD GPU support via HIP
     - OFF
     - ON/OFF
   * - ``ERF_ENABLE_SYCL``
     - Enables Intel GPU support via SYCL
     - OFF
     - ON/OFF
   * - ``ERF_ENABLE_NETCDF``
     - Enables NetCDF for I/O operations
     - OFF
     - ON/OFF
   * - ``ERF_ENABLE_NOAHMP``
     - Enables Noah-MP land surface model (requires ``ERF_ENABLE_NETCDF=ON``)
     - OFF
     - ON/OFF
   * - ``ERF_ENABLE_RRTMGP``
     - Enables RRTMGP radiation model (requires ``ERF_ENABLE_NETCDF=ON``, ``ERF_ENABLE_MPI=ON``)
     - OFF
     - ON/OFF
   * - ``ERF_ENABLE_SHOC``
     - Enables SHOC turbulence model (requires ``ERF_ENABLE_MPI=ON``)
     - OFF
     - ON/OFF
   * - ``ERF_ENABLE_P3``
     - Enables P3 microphysics model (requires ``ERF_ENABLE_MPI=ON``)
     - OFF
     - ON/OFF
   * - ``ERF_ENABLE_PARTICLES``
     - Enables support for Lagrangian particles
     - OFF
     - ON/OFF
   * - ``ERF_ENABLE_DOCUMENTATION``
     - Enables build target for generating Sphinx documentation
     - OFF
     - ON/OFF
   * - ``CMAKE_BUILD_TYPE``
     - Sets build configuration
     - Release
     - Release/Debug/RelWithDebInfo
   * - ``CMAKE_INSTALL_PREFIX``
     - Directory where ``make install`` places compiled artifacts
     - System-dependent
     - Path string

**2.4 Understanding Key Feature Flags**

The options above enable different ERF capabilities. The most commonly used flags and their interdependencies are:

**Core I/O and Analysis**

* ``ERF_ENABLE_NETCDF`` - Enables NetCDF I/O for reading WRF input files and writing plotfiles
* ``ERF_ENABLE_HDF5`` - Enables HDF5 parallel I/O backend (automatically enabled when NetCDF is enabled)
* ``ERF_ENABLE_FFT`` - Enables Fast Fourier Transform capabilities for spectral analysis

**Physics Packages**

* ``ERF_ENABLE_RRTMGP`` - Enables RRTMGP radiation model

  - **Requires:** ``ERF_ENABLE_NETCDF=ON`` and ``ERF_ENABLE_MPI=ON``
  - **Automatically enables:** ``ERF_ENABLE_EKAT=ON`` (provides Kokkos)

* ``ERF_ENABLE_SHOC`` - Enables SHOC turbulence and cloud macrophysics model

  - **Requires:** ``ERF_ENABLE_MPI=ON``
  - **Automatically enables:** ``ERF_ENABLE_EKAT=ON`` (provides Kokkos)
  - **Additional step:** Run ``source Build/GNU_Ekat/eamxx_clone.sh`` to initialize E3SM submodules

* ``ERF_ENABLE_P3`` - Enables P3 microphysics model

  - **Requires:** ``ERF_ENABLE_MPI=ON``
  - **Automatically enables:** ``ERF_ENABLE_EKAT=ON`` (provides Kokkos)

**GPU Acceleration**

For GPU builds, enable exactly one backend:

* ``ERF_ENABLE_CUDA`` - NVIDIA GPUs (requires CUDA Toolkit ≥ 11.0)
* ``ERF_ENABLE_HIP`` - AMD GPUs
* ``ERF_ENABLE_SYCL`` - Intel GPUs

.. note::
   The Kokkos-based physics packages (RRTMGP, SHOC, P3) support all three GPU backends through the EKAT library's Kokkos integration.

**2.5 Using Configuration Files**

For managing complex build configurations, CMake provides ``-C <file>`` option. This allows specifying build variables in CMake-syntax file:

.. code:: shell

   cmake -C path/to/config.cmake ..

Machine-specific profile scripts in ``Build/machines/`` directory are shell scripts sourced into terminal session before running CMake. They set up build environment by loading software modules and exporting environment variables. Example profiles:

* ``aurora_erf.profile``
* ``frontier_erf.profile``
* ``perlmutter_erf.profile``
* ``polaris_erf.profile``

**2.6 Advanced CMake Features**

**2.6.1 Hierarchical Logging**

For CMake versions 3.25 and newer, ERF's build configuration supports hierarchical logging:

* ``--log-level`` - Controls verbosity (VERBOSE shows dependency detection, DEBUG shows all diagnostics)
* ``--log-context`` - Enables message grouping by component

Example structured output:

.. code:: text

   [ERF.Cray] Detected Cray Programming Environment (CRAYPE_VERSION=2.7.30)
   [ERF.Cray] Setting Cray compiler wrappers...
   [ERF.Cray]   Set CMAKE_CXX_COMPILER = /opt/cray/pe/craype/default/bin/CC
   [ERF.AMReX] Using internal AMReX submodule
   [ERF.NetCDF] Found NetCDF: /opt/cray/pe/netcdf/4.9.0.9

3.0 Building the Documentation
------------------------------

**3.1 Overview and Procedure**

ERF uses Sphinx documentation generator to produce comprehensive project documentation from source files. Process is integrated into CMake build system.

Procedure to build documentation:

1. **Enable Documentation**: Set CMake option ``ERF_ENABLE_DOCUMENTATION=ON``
2. **Build Target**: Invoke docs target with ``make docs``

Process requires Sphinx installed on system. CMake configuration automatically searches for ``sphinx-build`` executable.
