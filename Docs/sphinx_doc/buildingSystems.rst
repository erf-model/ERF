.. _sec:build:systems:

Build Systems and Options
==========================

The Energy Research and Forecasting (ERF) code, built upon the AMReX framework, possesses a complex dependency graph that includes multiple advanced physics packages (RTE-RRTMGP, SHOC, P3). To ensure performance portability across diverse computing architectures, ERF relies on the Kokkos library. This document provides a detailed, code-level explanation of ERF's two supported build systems, GNU Make and CMake, outlining how each system manages these dependencies to produce a functional executable.

The two build systems serve distinct purposes. GNU Make provides fine-grained control over the compilation process and is an application-centric system for compiling a single executable for a specific scientific run. It provides debugging capabilities via utility targets like ``print-xxx``. CMake is designed for cross-platform compatibility and robust dependency management. It is a framework-centric approach with policy-based automation, such as its Cray Detection System. It creates versioned libraries, manages dependencies through automated detection, and supports testing through CTest.

This guide focuses on how the build systems are implemented and function, with direct references to the source code. It serves as a technical reference for developers and advanced users. For initial environment setup and quick start instructions, see :ref:`sec:build:overview` and the :ref:`GettingStarted` guide. For HPC-specific build instructions, see :ref:`sec:build:hpc`.

Directory Structure and Workflow
---------------------------------

ERF uses a paradigm where different executables are built in different subdirectories within the ``Exec`` directory. When using GNU Make, the user should build in the directory of the selected problem. When using CMake, separate executables are built for all problem directories listed in ``Exec/CMakeLists.txt``.

The problem directories within ``Exec`` are organized by purpose:

.. code-block:: text

   Exec/
   ├── ABL/                    # Atmospheric boundary layer (science runs)
   ├── DryRegTests/            # Dry atmospheric regression tests
   │   ├── IsentropicVortex/
   │   ├── TaylorGreenVortex/
   │   └── ...
   ├── MoistRegTests/          # Moist atmospheric regression tests
   │   ├── Bubble/
   │   ├── SquallLine/
   │   └── ...
   └── DevTests/               # Development and experimental features
       ├── MovingTerrain/
       └── ...

Each problem directory contains a README describing its purpose and functionality.

Building with GNU Make
----------------------

System Overview and Use Case
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The GNU Make system provides a direct path to producing a single, case-specific executable. It uses an application-centric approach, where developers and scientists compile code for a particular problem setup (e.g., atmospheric boundary layer simulation) directly within the ``Exec/`` directory structure.

The primary use case for the GNU Make system is scientific production runs and development scenarios where in-place compilation is preferred, avoiding overhead of creating and installing versioned libraries. The build is orchestrated by a ``GNUmakefile`` located within a specific case directory (such as ``Exec/ABL/``), which leverages build logic provided by the underlying AMReX framework.

How it Works: The Orchestration Process
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ERF GNU Make process uses a hierarchy of includes that separates user configuration from application and framework build logic. The user's ``GNUmakefile`` serves as the top layer, defining high-level build options. This file includes the application-level ``Make.ERF`` to gather all necessary ERF source files, which includes the core framework-level ``Make.defs`` and ``Make.rules`` from AMReX to execute compilation and linking.

The process follows these steps:

1. **GNUmakefile Location**: Process begins with user invoking ``make`` in application-specific directory, such as ``Exec/ABL/``. This directory contains the primary control file, ``GNUmakefile``.

2. **Set AMREX_HOME**: First step within ``GNUmakefile`` is defining ``AMREX_HOME`` variable. This points to location of AMReX submodule containing core build logic. Default path is ``$(ERF_HOME)/Submodules/AMReX``.

3. **Define Build Variables**: User defines boolean flags and other variables in ``GNUmakefile`` to control which features are compiled. These include options like ``USE_MPI``, ``USE_CUDA``, ``USE_RRTMGP``.

4. **Include ERF Sources**: ``GNUmakefile`` includes file ``Exec/Make.ERF``. This file systematically adds all necessary ERF source and header directories to build paths, populating ``VPATH_LOCATIONS`` and ``INCLUDE_LOCATIONS`` variables.

5. **Include AMReX Core Logic**: Finally, ``Make.ERF`` includes core AMReX makefiles, ``Make.defs`` and ``Make.rules``. These files contain generic logic for discovering dependencies, compiling source files, and linking object files into final executable.

Technical Implementation
~~~~~~~~~~~~~~~~~~~~~~~~

The GNU Make build system provides fine-grained control over the compilation process. It uses a series of ``Make.package`` files and configuration variables (e.g., ``USE_RRTMGP``, ``USE_NETCDF``) within the main ``Exec/Make.ERF`` file. These variables conditionally add source files to the build and pass preprocessor definitions (e.g., ``-DERF_USE_RRTMGP``) to the compiler.

The GNU Make system is particularly useful for:

* Scientific production runs where in-place compilation is preferred
* Development scenarios requiring direct control over compiler flags
* Debugging with ``make print-<variable>`` to inspect build configuration
* Building single executables without overhead of library versioning

Environment Variables
~~~~~~~~~~~~~~~~~~~~~

The GNU Make build system leverages setting environment variables for dependency directories. All dependencies except SHOC and P3 are provided as git submodules and can be populated using:

.. code-block:: bash

   git submodule update --init --recursive

Or before cloning:

.. code-block:: bash

   git clone --recursive https://github.com/erf-model/ERF.git

Although submodules are provided, dependencies can be placed externally as long as the ``<REPO_HOME>`` environment variables are set correctly. Example configuration in ``.bashrc``:

.. code-block:: bash

   export ERF_HOME=${HOME}/ERF
   export AMREX_HOME=${ERF_HOME}/Submodules/AMReX

The GNU Make system is set up to use the path to AMReX submodule by default, so it is not necessary to set the AMReX path explicitly. It is also possible to use an external version of AMReX, downloaded by running:

.. code-block:: bash

   git clone https://github.com/amrex-codes/amrex.git

in which case the ``AMREX_HOME`` environment variable must point to the location where AMReX has been downloaded, which will take precedence over the default path to the submodule. If using bash shell:

.. code-block:: bash

   export AMREX_HOME=/path/to/external/amrex

or if using tcsh:

.. code-block:: bash

   setenv AMREX_HOME /path/to/external/amrex

Building with SHOC or P3
~~~~~~~~~~~~~~~~~~~~~~~~

To build with SHOC or P3 using GNU Make:

.. code-block:: bash

   export ERF_DIR=/path/to/ERF
   source /path/to/ERF/Build/GNU_Ekat/eamxx_clone.sh
   source /path/to/ERF/Build/GNU_Ekat/ekat_build_commands.sh

Then follow the instructions below, ensuring that you have ``USE_SHOC=TRUE`` (when running with SHOC) and ``USE_P3=TRUE`` (when running with P3) in your GNUmakefile.

Build Steps
~~~~~~~~~~~

1. Navigate to the desired build directory:

   .. code-block:: bash

      cd ERF/Exec/DryRegTests/IsentropicVortex/

2. Edit the ``GNUmakefile`` to set build options:

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
        - None (required)
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
      * - ``USE_PARTICLES``
        - Enables support for Lagrangian particles
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
      * - ``USE_MULTIBLOCK``
        - Enables multiblock capability
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
      * - ``DEBUG``
        - Enables debug mode
        - FALSE
        - TRUE/FALSE
      * - ``PROFILE``
        - Includes profiling info
        - FALSE
        - TRUE/FALSE
      * - ``TINY_PROFILE``
        - Includes tiny profiling info
        - FALSE
        - TRUE/FALSE
      * - ``COMM_PROFILE``
        - Includes comm profiling info
        - FALSE
        - TRUE/FALSE
      * - ``TRACE_PROFILE``
        - Includes trace profiling info
        - FALSE
        - TRUE/FALSE

   .. note::
      **At most one of USE_OMP, USE_CUDA, USE_HIP, USE_SYCL should be set to TRUE.**

   Information on using other compilers can be found in the AMReX documentation at https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html

3. Build the executable:

   .. code-block:: bash

      make

   The name of the resulting executable (generated by the GNUmake system) encodes several of the build characteristics, including dimensionality of the problem, compiler name, and whether MPI and/or OpenMP were linked with the executable. Thus, several different build configurations may coexist simultaneously in a problem folder. For example, the default build in ``ERF/Exec/DryRegTests/IsentropicVortex`` will look like ``ERF3d.gnu.MPI.ex``, indicating that this is a 3-d version of the code, made with ``COMP=gnu``, and ``USE_MPI=TRUE``.

Example GNUmakefile
~~~~~~~~~~~~~~~~~~~

Typical ``GNUmakefile`` examples for ERF applications:

.. dropdown:: Exec/ABL/GNUmakefile
   :icon: code

   .. literalinclude:: ../../Exec/ABL/GNUmakefile
      :language: makefile

.. dropdown:: Exec/DryRegTests/IsentropicVortex/GNUmakefile
   :icon: code

   .. literalinclude:: ../../Exec/DryRegTests/IsentropicVortex/GNUmakefile
      :language: makefile

Common Commands
~~~~~~~~~~~~~~~

The AMReX make system provides standard utility targets for managing build environment and artifacts. For complete list of commands and advanced features, refer to `AMReX documentation <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html>`_. Common commands include:

* ``make`` - Compiles source code and creates executable in current directory
* ``make clean`` - Removes all build artifacts for all configurations
* ``make cleanconfig`` - Removes artifacts for current configuration only
* ``make print-<variable>`` - Debugging tool that prints value of specified make variable (e.g., ``make print-CXXFLAGS``)

Build Information
~~~~~~~~~~~~~~~~~

The build information can be accessed by typing:

.. code-block:: bash

   ./ERF*ex --describe

in the directory where the executable has been built.

Customization with Make.local
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For user- or site-specific customizations that should not be committed to main repository, the AMReX GNU Make system provides ``Make.local`` file. This file, located at ``amrex/Tools/GNUMake/Make.local``, is the primary mechanism for overriding default build variables.

Settings defined in ``Make.local`` apply globally to all projects using that specific instance of AMReX submodule. Common use case is specifying non-default compiler version:

.. code-block:: bash

   CXX = g++-8
   CC  = gcc-8
   FC  = gfortran-8

Building Documentation with GNU Make
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Building the ERF documentation with GNU Make can be completed by navigating to the ``/ERF/Docs/`` directory and executing the following command:

.. code-block:: bash

   source BuildDocs.sh

Note that both the Sphinx and Doxygen documentation will be built.

Building with CMake
-------------------

System Overview and Rationale
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CMake is the recommended build system for ERF, particularly on HPC systems. It is a cross-platform tool with framework-centric philosophy. Instead of producing single executable for one use case, CMake produces versioned, exportable libraries. This design enables generation and installation of ``find_package``-compatible configuration file (``ERFConfig.cmake``), allowing other CMake-based projects to consume ERF as dependency with simple ``find_package(ERF)`` command.

Key advantages of using CMake for ERF:

* Promotes clean, out-of-source builds keeping compiled artifacts separate from source code
* Provides robust dependency detection through ``find_package`` mechanism
* Automatically detects and configures for Cray HPC systems, applying necessary workarounds
* Supports CTest for testing and verification of ERF

Technical Implementation
~~~~~~~~~~~~~~~~~~~~~~~~

The CMake build system uses a series of ``option()`` commands in its ``CMakeLists.txt`` file (e.g., ``option(ERF_ENABLE_MPI "Enable MPI" OFF)``) to generate user-configurable cache variables. These variables control which features are included in the build and how dependencies are linked.

CMake is particularly useful for:

* HPC environments with complex module systems (automatic Cray detection)
* Cross-platform builds requiring robust dependency detection
* Projects that need ``find_package`` integration with ERF libraries
* Testing and verification through CTest integration

Prerequisites
~~~~~~~~~~~~~

Compiling with CMake involves an additional configure step before using the ``make`` command and it is expected that the user has cloned the ERF repo with the ``--recursive`` option or performed ``git submodule init; git submodule update`` in the ERF repo to populate its submodules.

Building with SHOC or P3
~~~~~~~~~~~~~~~~~~~~~~~~

To build with SHOC or P3 using CMake, you will need to make sure you have run ``git submodule update --init --recursive``, then:

.. code-block:: bash

   export ERF_DIR=/path/to/ERF
   source /path/to/ERF/Build/GNU_Ekat/eamxx_clone.sh

Then follow the guidance below, making sure to set ``ERF_ENABLE_SHOC`` and/or ``ERF_ENABLE_P3`` to TRUE.

Standard Workflow
~~~~~~~~~~~~~~~~~

ERF offers several methods for invoking CMake build:

* **Simple Script**: Run ``Build/cmake.sh`` for basic, tested configuration
* **Advanced Script**: Use specialized scripts like ``Build/cmake_with_kokkos_many_cuda.sh`` for specific use cases
* **Manual Invocation**: Create build directory and invoke CMake: ``cmake -B build -S .``

Standard manual procedure from project root:

.. code-block:: bash

   mkdir build          # Create directory for build artifacts
   cd build            # Change to build directory
   cmake [options] ..  # Configure project
   make                # Compile
   make install        # Install (optional)

An example CMake configuration script:

.. dropdown:: Build/cmake_cuda.sh
   :icon: code

   .. literalinclude:: ../../Build/cmake_cuda.sh
      :language: bash

For manual configuration, a typical command to build ERF with MPI:

.. code-block:: bash

   cmake -DCMAKE_BUILD_TYPE:STRING=Release \
         -DERF_ENABLE_MPI:BOOL=ON \
         -DCMAKE_CXX_COMPILER:STRING=mpicxx \
         -DCMAKE_C_COMPILER:STRING=mpicc \
         -DCMAKE_Fortran_COMPILER:STRING=mpifort \
         .. && make

Typically, a user will create a ``build`` directory in the project directory and execute the configuration from said directory (``cmake <options> ..``) before building. Note that CMake is able to generate makefiles for the Ninja build system as well which will allow for faster building of the executable(s).

Alternative workflows for using build scripts:

.. tab-set::

   .. tab-item:: From Build/ directory

      .. code-block:: bash

         mkdir build && cd build
         ../Build/cmake.sh
         make install

   .. tab-item:: From ERF root

      .. code-block:: bash

         ./Build/cmake_with_kokkos_many.sh

      This creates ``build_erf/`` and ``install_erf/`` directories.

ERF CMake Options
~~~~~~~~~~~~~~~~~

The ERF build is customized by passing options to cmake command using ``-D<VARIABLE>=<VALUE>`` syntax.

.. list-table:: ERF CMake Options
   :header-rows: 1
   :widths: 25 45 10 20

   * - Variable Name
     - Description
     - Default
     - Possible Values
   * - ``CMAKE_BUILD_TYPE``
     - Sets build configuration
     - Release
     - Release/Debug/RelWithDebInfo
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
   * - ``ERF_ENABLE_PARTICLES``
     - Enables support for Lagrangian particles
     - OFF
     - ON/OFF
   * - ``ERF_ENABLE_MULTIBLOCK``
     - Enables multiblock capability
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
   * - ``ERF_ENABLE_TESTS``
     - Enables CTest test suite
     - OFF
     - ON/OFF
   * - ``ERF_ENABLE_FCOMPARE``
     - Enables fcompare utility for regression testing
     - OFF
     - ON/OFF
   * - ``ERF_ENABLE_DOCUMENTATION``
     - Enables build target for generating Sphinx documentation
     - OFF
     - ON/OFF
   * - ``CMAKE_INSTALL_PREFIX``
     - Directory where ``make install`` places compiled artifacts
     - System-dependent
     - Path string

.. note::
   **At most one of ERF_ENABLE_OPENMP, ERF_ENABLE_CUDA, ERF_ENABLE_HIP and ERF_ENABLE_SYCL should be set to TRUE.**

Understanding Key Feature Flags
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

CMake Logging Options
~~~~~~~~~~~~~~~~~~~~~

ERF's CMake configuration supports hierarchical logging to help diagnose build issues and understand configuration decisions. These options are available with CMake 3.25+:

**Logging Levels**

The ``--log-level`` flag controls the verbosity of CMake output:

.. code-block:: bash

   # Quiet output (STATUS messages only)
   cmake ..

   # Show detection details
   cmake --log-level=VERBOSE ..

   # Show all diagnostics
   cmake --log-level=DEBUG ..

**Logging Context**

The ``--log-context`` flag shows the hierarchy of messages, making it easier to understand which component is reporting:

.. code-block:: bash

   # Show message hierarchy with component names
   cmake --log-context ..

   # Combine with verbose output for detailed diagnostics
   cmake --log-context --log-level=VERBOSE ..

Example output with ``--log-context``:

.. code-block:: text

   [ERF.Cray] Detected Cray Programming Environment (CRAYPE_VERSION=2.7.30)
   [ERF.Cray] Setting Cray compiler wrappers...
   [ERF.Cray]   Set CMAKE_CXX_COMPILER = /opt/cray/pe/craype/default/bin/CC
   [ERF.AMReX] Using internal AMReX submodule
   [ERF.NetCDF] Found NetCDF: /opt/cray/pe/netcdf/4.9.0.9

CMake Utility Targets
~~~~~~~~~~~~~~~~~~~~~

ERF provides several utility targets to help manage builds:

**Clean Build Artifacts**

To perform a complete clean (equivalent to ``make distclean`` in GNU Make):

.. code-block:: bash

   make distclean

This removes all CMake configuration and build artifacts, including:

- CMake cache and generated files (``CMakeCache.txt``, ``CMakeFiles/``, etc.)
- Build outputs (executables, libraries)
- Generated configuration files
- Test outputs

The install directory is preserved.

**Uninstall**

To uninstall files that were installed via ``make install`` or ``cmake --install``:

.. code-block:: bash

   make uninstall

**Show Detected Configuration**

On Cray systems, you can view the auto-detected configuration:

.. code-block:: bash

   make show-cray-config

This displays the configuration that was automatically detected and can be saved to use as a starting point for manual configuration.

Using Configuration Files
~~~~~~~~~~~~~~~~~~~~~~~~~~

On Cray systems, ERF automatically detects the system configuration and applies necessary workarounds. If you need to use a manual configuration, you can create a configuration file and use it with CMake's ``-C`` option:

.. code-block:: bash

   cmake -C path/to/config.cmake ..

ERF provides machine-specific profile examples in ``Build/machines/``:

- ``perlmutter_erf.profile`` - NERSC Perlmutter (NVIDIA A100)
- ``frontier_erf.profile`` - OLCF Frontier (AMD MI250X)
- ``polaris_erf.profile`` - ALCF Polaris (NVIDIA A100)
- ``aurora_erf.profile`` - ALCF Aurora (Intel GPUs)

These profiles show the recommended modules to load for each system. They are shell scripts that should be sourced into your terminal session before running CMake. For detailed information about using these profiles, see :ref:`sec:build:hpc`.

Building Documentation with CMake
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Building the ERF documentation with CMake can be completed by configuring with the flag ``-DERF_ENABLE_DOCUMENTATION:BOOL=ON`` and then compiling the following target:

.. code-block:: bash

   make docs

Note that both the Sphinx and Doxygen documentation will be built.

Developer Testing Tool
~~~~~~~~~~~~~~~~~~~~~~

For systematic testing of multiple build configurations, ERF provides a validation script:

.. dropdown:: setup_cmake_validation.sh
   :icon: code

   .. literalinclude:: ../Build/setup_cmake_validation.sh
      :language: bash

   **Usage:**

   .. code-block:: bash

      cd Build
      ./setup_cmake_validation.sh default    # or perlmutter, gnu_ekat
      cd ../build_default
      ./run.sh     # List available scripts
      ./run.sh 1   # Run first script in isolated directory

   This creates isolated build directories for each configuration script, preventing interference between builds. Each script runs in its own subdirectory (``script_<name>/``) at the ERF root.
