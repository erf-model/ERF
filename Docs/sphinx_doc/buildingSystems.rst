.. _sec:build:systems:

Build Systems and Options
==========================

This document explains ERF's two build systems, GNU Make and CMake, and how each manages dependencies to produce a functional executable.

The two build systems serve distinct purposes:

* **GNU Make** - Application-centric system for compiling a single executable for a specific scientific run. Provides fine-grained control over compilation and debugging via utility targets like ``print-xxx``.

* **CMake** - Framework-centric system designed for cross-platform compatibility and dependency management. Creates versioned libraries, provides automated detection (including Cray systems), and supports testing through CTest. Recommended for HPC environments.

This guide serves as a technical reference for developers and advanced users. For initial setup and quick start instructions, see :ref:`sec:build:quickstart` and :ref:`GettingStarted`. For HPC-specific instructions, see :ref:`sec:build:hpc`.

Directory Structure and Workflow
---------------------------------

ERF builds executables in subdirectories within ``Exec``. With GNU Make, build in the specific problem directory. With CMake, configure once to build all executables listed in ``Exec/CMakeLists.txt``.

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

System Overview
~~~~~~~~~~~~~~~

The GNU Make system provides a direct path to producing a single, case-specific executable. It uses an application-centric approach, where developers and scientists compile code for a particular problem setup (e.g., atmospheric boundary layer simulation) directly within the ``Exec/`` directory structure.

Primary use cases:

* Scientific production runs with in-place compilation
* Development scenarios requiring direct control over compiler flags
* Debugging build configuration
* Single executables without library versioning overhead

The build is orchestrated by a ``GNUmakefile`` in each case directory (such as ``Exec/ABL/``), which uses build logic from the AMReX framework.

.. dropdown:: How it Works: The Orchestration Process
   :icon: info

   The GNU Make process uses a hierarchy of includes separating user configuration from application and framework build logic:

   1. **GNUmakefile Location**: User invokes ``make`` in an application directory, such as ``Exec/ABL/``, which contains the ``GNUmakefile`` control file.

   2. **Set AMREX_HOME**: The ``GNUmakefile`` defines ``AMREX_HOME``, pointing to the AMReX submodule containing core build logic. Default path is ``$(ERF_HOME)/Submodules/AMReX``.

   3. **Define Build Variables**: User defines boolean flags and variables in ``GNUmakefile`` to control which features are compiled (e.g., ``USE_MPI``, ``USE_CUDA``, ``USE_RRTMGP``).

   4. **Include ERF Sources**: ``GNUmakefile`` includes ``Exec/Make.ERF``, which adds ERF source and header directories to build paths, populating ``VPATH_LOCATIONS`` and ``INCLUDE_LOCATIONS`` variables.

   5. **Include AMReX Core Logic**: ``Make.ERF`` includes core AMReX makefiles (``Make.defs`` and ``Make.rules``), which contain logic for discovering dependencies, compiling source files, and linking object files into the final executable.

.. dropdown:: Technical Implementation
   :icon: code

   The GNU Make build system uses ``Make.package`` files and configuration variables (e.g., ``USE_RRTMGP``, ``USE_NETCDF``) in the main ``Exec/Make.ERF`` file. These variables conditionally add source files to the build and pass preprocessor definitions (e.g., ``-DERF_USE_RRTMGP``) to the compiler.

   This provides:

   * Direct control over which source files are compiled
   * Fine-grained compiler flag customization
   * Debugging via ``make print-<variable>`` to inspect configuration
   * Multiple build configurations coexisting in the same directory

Build Steps
~~~~~~~~~~~

**1. Clone Repository and Populate Dependencies**

See :ref:`sec:submodules` for more details. All dependencies except SHOC and P3 are provided as git submodules:

.. code-block:: bash

   # Clone with submodules
   git clone --recursive https://github.com/erf-model/ERF.git

   # Or populate submodules in existing clone
   git submodule update --init --recursive

The GNU Make system uses the AMReX submodule path by default. To use an external AMReX installation, set ``AMREX_HOME``:

.. code-block:: bash

   # Download external AMReX
   git clone https://github.com/amrex-codes/amrex.git

   # Set environment variable (bash)
   export AMREX_HOME=/path/to/external/amrex

   # Or for tcsh
   setenv AMREX_HOME /path/to/external/amrex

**2. Setup for SHOC or P3 (Optional)**

If building with SHOC or P3, run the setup scripts:

.. code-block:: bash

   export ERF_DIR=/path/to/ERF
   source /path/to/ERF/Build/GNU_Ekat/eamxx_clone.sh
   source /path/to/ERF/Build/GNU_Ekat/ekat_build_commands.sh

Then set ``USE_SHOC=TRUE`` or ``USE_P3=TRUE`` in your GNUmakefile (step 4).

**3. Navigate to Problem Directory**

.. code-block:: bash

   cd ERF/Exec/DryRegTests/IsentropicVortex/

**4. Edit GNUmakefile**

Set build variables in the ``GNUmakefile``:

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
      **At most one of USE_OMP, USE_CUDA, USE_HIP, USE_SYCL should be TRUE.**

   For additional compiler options, see the `AMReX documentation <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html>`_.

.. dropdown:: Example GNUmakefile
   :icon: code

   Typical ``GNUmakefile`` examples:

   **Exec/ABL/GNUmakefile:**

   .. literalinclude:: ../../Exec/ABL/GNUmakefile
      :language: makefile

   **Exec/DryRegTests/IsentropicVortex/GNUmakefile:**

   .. literalinclude:: ../../Exec/DryRegTests/IsentropicVortex/GNUmakefile
      :language: makefile

**5. Build**

.. code-block:: bash

   make

The executable name encodes build characteristics (dimensionality, compiler, parallelization). For example, in ``Exec/DryRegTests/IsentropicVortex`` with ``COMP=gnu`` and ``USE_MPI=TRUE``, the executable is ``ERF3d.gnu.MPI.ex``. Multiple build configurations can coexist in the same directory.

**6. Verify Build (Optional)**

View build configuration:

.. code-block:: bash

   ./ERF*ex --describe

.. dropdown:: Common Commands
   :icon: terminal

   Standard utility targets:

   * ``make`` - Compile and create executable
   * ``make clean`` - Remove all build artifacts
   * ``make cleanconfig`` - Remove current configuration artifacts
   * ``make print-<variable>`` - Print value of make variable (e.g., ``make print-CXXFLAGS``)

   For complete documentation, see the `AMReX build guide <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html>`_.

.. dropdown:: Customization with Make.local
   :icon: tools

   For user- or site-specific customizations, create ``amrex/Tools/GNUMake/Make.local``. Settings apply globally to all projects using that AMReX instance. Example for specifying compiler version:

   .. code-block:: bash

      CXX = g++-8
      CC  = gcc-8
      FC  = gfortran-8

.. dropdown:: Running the Executable
   :icon: play

   Navigate to the executable location based on your chosen workflow:

   .. code-block:: bash

      # If built in Build/ directory:
      cd Build/Exec/ABL
      mpiexec -n 4 ./erf_abl ../../../Exec/ABL/inputs_most

      # If out-of-source build (without install):
      cd build/Exec/ABL
      mpiexec -n 4 ./erf_abl ../../../Exec/ABL/inputs_most

      # If installed to install/:
      cd install/bin
      mpiexec -n 4 ./erf_abl ../../Exec/ABL/inputs_most

      # If using cmake_with_kokkos_many.sh with defaults:
      cd install/bin  # or cd $ERF_INSTALL_DIR/bin if customized or cd $ERF_BUILD_DIR/Exec/ABL
      mpiexec -n 4 ./erf_abl ../../Exec/ABL/inputs_most

   For details on input files and job submission, see :ref:`sec:running`.

Building Documentation
~~~~~~~~~~~~~~~~~~~~~~

Build ERF documentation with GNU Make:

.. code-block:: bash

   cd ERF/Docs/
   source BuildDocs.sh

This builds both Sphinx and Doxygen documentation.

Building with CMake
-------------------

System Overview
~~~~~~~~~~~~~~~

CMake is a cross-platform build tool using a framework-centric approach. Instead of producing a single executable for one use case, CMake produces versioned, exportable libraries as well as multiple executables across the problem directory setups. This enables generation of ``find_package``-compatible configuration files (``ERFConfig.cmake``), allowing other CMake projects to consume ERF as a dependency.

CMake provides:

* Out-of-source builds keeping artifacts separate from source code
* Robust dependency detection through ``find_package``
* Automated detection and configuration for Cray HPC systems
* Testing and verification through CTest

CMake is recommended for HPC environments where automated dependency detection and robust cross-platform builds are important.

.. dropdown:: Technical Implementation
   :icon: code

   CMake uses ``option()`` commands in ``CMakeLists.txt`` (e.g., ``option(ERF_ENABLE_MPI "Enable MPI" OFF)``) to generate user-configurable cache variables. These variables control which features are included and how dependencies are linked.

   This provides:

   * Automated dependency detection through ``find_package``
   * Policy-based configuration (e.g., Cray Detection System)
   * Clean separation between source and build artifacts
   * CTest integration for testing
   * Exportable libraries for downstream projects

Build Steps
~~~~~~~~~~~

**1. Clone Repository and Prerequisites**

Clone with submodules:

.. code-block:: bash

   # Clone with submodules
   git clone --recursive https://github.com/erf-model/ERF.git

   # Or populate submodules in existing clone
   git submodule update --init --recursive

**2. Setup for SHOC or P3 (Optional)**

If building with SHOC or P3:

.. code-block:: bash

   export ERF_DIR=/path/to/ERF
   source /path/to/ERF/Build/GNU_Ekat/eamxx_clone.sh

Then configure with ``-DERF_ENABLE_SHOC=TRUE`` and/or ``-DERF_ENABLE_P3=TRUE`` (step 4).

**3. Choose Build Workflow**

ERF supports multiple CMake workflows. The main difference is directory structure and whether you use the install step.

.. tab-set::

   .. tab-item:: In Build/ Directory

      Build artifacts and executables in the Build/ directory tree. Executables appear in ``Exec/<problem>/`` subdirectories, similar to GNU Make workflow.

      .. code-block:: bash

         cd Build
         ./cmake.sh

      **Executable locations:** ``Build/Exec/ABL/erf_abl``, ``Build/Exec/RegTests/Bubble/erf_bubble``, etc.

      **Cleanup for rebuild:**

      .. code-block:: bash

         make distclean
         ./cmake.sh

   .. tab-item:: Out-of-Source Build

      Build directory separate from source. Executables in ``build/Exec/<problem>/`` subdirectories. Optionally use install step to copy all executables to a single ``install/bin/`` directory.

      .. code-block:: bash

         mkdir build && cd build
         ../Build/cmake.sh
         make install  # optional - copies to install/bin/ (may be needed for builds that require kokkos)

      **Executable locations:** ``build/Exec/ABL/erf_abl``, etc., and optionally ``install/bin/erf_abl`` (if installed)

      **Cleanup for rebuild:**

      .. code-block:: bash

         rm -rf build/ install/  # complete cleanup
         mkdir build && cd build

.. tab-item:: Script with Customizable Directories

      Uses environment variables for directory control. Defaults: build in current dir (`.`), source from parent (`..`), install to ``install/``.

      **Customize directories (optional):**

      Set ``ERF_BUILD_DIR``, ``ERF_SOURCE_DIR``, ``ERF_INSTALL_DIR``, or ``ERF_HOME`` environment variables.

      .. code-block:: bash

         # From ERF repository root (sets absolute paths)
         ERF_HOME=$(pwd) ./Build/cmake_with_kokkos_many.sh

         # Or customize individual directories
         ERF_BUILD_DIR=build ERF_INSTALL_DIR=install ./Build/cmake_with_kokkos_many.sh

         # Or use defaults from any directory
         cd ERF
         ./Build/cmake_with_kokkos_many.sh

      **Executable locations:** ``$ERF_BUILD_DIR/Exec/ABL/erf_abl`` and ``$ERF_INSTALL_DIR/bin/erf_abl`` (defaults: ``./Exec/ABL/erf_abl`` and ``install/bin/erf_abl``)

      **Cleanup for rebuild:**

      .. code-block:: bash

         # If using defaults
         rm -rf . install/  # from build directory

         # If using ERF_HOME
         rm -rf $ERF_HOME/build $ERF_HOME/install

         # Then rebuild
         ./Build/cmake_with_kokkos_many.sh

.. dropdown:: Workflow Comparison
   :icon: info

   **Executable Locations:**

   * **After make:** Executables always in ``<build_dir>/Exec/<problem>/erf_<problem>`` (mirrors source tree)
   * **After make install:** Executables copied to ``<install_prefix>/bin/erf_<problem>`` (all in one location)

   **Install step is optional** - it copies executables from Exec subdirectories to a single bin/ directory for convenience.

   **Choose based on:**

   * Familiarity with GNU Make workflow → in Build/ directory
   * Want separate build/source trees → out-of-source build
   * Prefer all executables in one place → use install step
   * Running multiple configurations → separate build directories
   * Minimal typing from repository root → use script with -B/-S

**4. Configure with Options**

For manual configuration, use ``-D<VARIABLE>=<VALUE>`` syntax:

.. code-block:: bash

   cmake -DCMAKE_BUILD_TYPE=Release \
         -DERF_ENABLE_MPI=ON \
         -DCMAKE_CXX_COMPILER=mpicxx \
         -DCMAKE_C_COMPILER=mpicc \
         -DCMAKE_Fortran_COMPILER=mpifort \
         .. && make

CMake can also generate makefiles for the Ninja build system for faster compilation.

.. dropdown:: Example Configuration Script
   :icon: code

   **Build/cmake_cuda.sh:**

   .. literalinclude:: ../../Build/cmake_cuda.sh
      :language: bash


**CMake Options:**

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
   **At most one of ERF_ENABLE_OPENMP, ERF_ENABLE_CUDA, ERF_ENABLE_HIP, ERF_ENABLE_SYCL should be ON.**

.. dropdown:: Feature Dependencies
   :icon: checklist

   **I/O and Analysis**

   * ``ERF_ENABLE_NETCDF`` - Enables NetCDF I/O for reading WRF input files and writing plotfiles
   * ``ERF_ENABLE_HDF5`` - Automatically enabled when NetCDF is enabled (provides parallel I/O backend)
   * ``ERF_ENABLE_FFT`` - Enables Fast Fourier Transform for spectral analysis

   **Physics Packages**

   * ``ERF_ENABLE_RRTMGP`` - RRTMGP radiation model

     - Requires ``ERF_ENABLE_NETCDF=ON`` and ``ERF_ENABLE_MPI=ON``
     - Automatically enables ``ERF_ENABLE_EKAT=ON`` (provides Kokkos)

   * ``ERF_ENABLE_SHOC`` - SHOC turbulence and cloud macrophysics

     - Requires ``ERF_ENABLE_MPI=ON``
     - Automatically enables ``ERF_ENABLE_EKAT=ON`` (provides Kokkos)
     - Additional step: Run ``source Build/GNU_Ekat/eamxx_clone.sh``

   * ``ERF_ENABLE_P3`` - P3 microphysics

     - Requires ``ERF_ENABLE_MPI=ON``
     - Automatically enables ``ERF_ENABLE_EKAT=ON`` (provides Kokkos)

   * ``ERF_ENABLE_NOAHMP`` - NOAHMP land surface model

     - Requires the NetCDF Fortran library

   **GPU Acceleration**

   Enable exactly one GPU backend:

   * ``ERF_ENABLE_CUDA`` - NVIDIA GPUs (requires CUDA Toolkit ≥ 11.0)
   * ``ERF_ENABLE_HIP`` - AMD GPUs
   * ``ERF_ENABLE_SYCL`` - Intel GPUs

   .. note::
      Kokkos-based physics packages (RRTMGP, SHOC, P3) support all three GPU backends through EKAT's Kokkos integration.

.. dropdown:: Logging Options (CMake 3.25+)
   :icon: bug

   Use hierarchical logging to diagnose build issues:

   .. code-block:: bash

      # Default output
      cmake ..

      # Show dependency detection
      cmake --log-level=VERBOSE ..

      # Show all diagnostics
      cmake --log-level=DEBUG ..

      # Show message hierarchy
      cmake --log-context ..

      # Combine for detailed output
      cmake --log-context --log-level=VERBOSE ..

   Example output with ``--log-context``:

   .. code-block:: text

      [ERF.Cray] Detected Cray Programming Environment (CRAYPE_VERSION=2.7.30)
      [ERF.Cray] Setting Cray compiler wrappers...
      [ERF.Cray]   Set CMAKE_CXX_COMPILER = /opt/cray/pe/craype/default/bin/CC
      [ERF.AMReX] Using internal AMReX submodule
      [ERF.NetCDF] Found NetCDF: /opt/cray/pe/netcdf/4.9.0.9

.. dropdown:: Utility Targets
   :icon: tools

   **Clean Build Artifacts**

   Remove all CMake configuration and build artifacts:

   .. code-block:: bash

      make distclean

   Removes:

   - CMake cache and generated files (``CMakeCache.txt``, ``CMakeFiles/``, etc.)
   - Build outputs (executables, libraries)
   - Generated configuration files
   - Test outputs

   The install directory is preserved.

   **Uninstall**

   Remove files installed via ``make install``:

   .. code-block:: bash

      make uninstall

   **Show Cray Configuration**

   On Cray systems, view auto-detected configuration:

   .. code-block:: bash

      make show-cray-config

   This displays detected configuration and can be saved for manual configuration.

.. dropdown:: Configuration Files
   :icon: file

   On Cray systems, ERF automatically detects configuration. For manual configuration, use the ``-C`` option:

   .. code-block:: bash

      cmake -C path/to/config.cmake ..

   Machine-specific profiles in ``Build/machines/``:

   - ``perlmutter_erf.profile`` - NERSC Perlmutter (NVIDIA A100)
   - ``frontier_erf.profile`` - OLCF Frontier (AMD MI250X)
   - ``polaris_erf.profile`` - ALCF Polaris (NVIDIA A100)
   - ``aurora_erf.profile`` - ALCF Aurora (Intel GPUs)

   These profiles are shell scripts that load required modules and set environment variables. Source them before running CMake. For detailed information, see :ref:`sec:build:hpc`.

.. dropdown:: Developer Testing
   :icon: beaker

   For systematic testing of multiple build configurations:

   .. literalinclude:: ../Build/setup_cmake_validation.sh
      :language: bash

   **Usage:**

   .. code-block:: bash

      cd Build
      ./setup_cmake_validation.sh default    # or perlmutter, gnu_ekat
      cd ../build_default
      ./run.sh     # List available scripts
      ./run.sh 1   # Run first script

   Creates isolated build directories for each configuration script. Each script runs in its own subdirectory (``script_<name>/``) at the ERF root.

Building Documentation
~~~~~~~~~~~~~~~~~~~~~~

Build ERF documentation with CMake:

.. code-block:: bash

   cmake -DERF_ENABLE_DOCUMENTATION=ON ..
   make docs

This builds both Sphinx and Doxygen documentation.
