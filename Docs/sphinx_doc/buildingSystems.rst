.. _sec:build:systems:

Build Systems and Options
==========================

This document explains ERF's two build systems, GNU Make and CMake, and how each manages dependencies to produce a functional executable.

The two build systems serve distinct purposes:

* **GNU Make** - Application-centric system for compiling a single executable for a specific scientific run. Provides fine-grained control over compilation and debugging via utility targets like ``print-xxx``.

* **CMake** - Framework-centric system designed for cross-platform compatibility and dependency management. Creates versioned libraries, provides automated detection (including Cray systems), and supports testing through CTest. Recommended for HPC environments.

This guide serves as a technical reference for developers and advanced users. For initial setup and quick start instructions, see :ref:`sec:build:quickstart` and :ref:`GettingStarted`. For HPC concepts, see :ref:`sec:build:hpc`; for machine-specific build/run workflows, see :ref:`sec:hpc:guides`.

Directory Structure and Workflow
---------------------------------

ERF builds executables in ``ERF/Exec`` when using GNU Make. The exception is development-test work under ``ERF/.Exec_dev``, where GNU Make builds in a problem-specific directory under ``ERF/.Exec_dev``.

With CMake, one configures once to build the core libraries and the shared test executable in ``build`` (e.g., ``erf_exec``).

Regardless of how the executable is built, it can be run on regression and canonical tests that have inputs files in subdirectories of ``Exec/RegTests`` and ``Exec/CanonicalTests``.
Ctests are run on files in subdirectories of ``ERF/Tests/test_files``.

Directories within ``Exec`` are organized into cases used primarily for regression testing and cases representing more canonical tests.

.. code-block:: text

   Exec/
   ├── RegTests/               # Regression test input decks
   │   ├── IsentropicVortex/
   │   ├── TaylorGreenVortex/
   │   ├── Bubble/
   │   └── ...
   ├── CanonicalTests/         # Canonical atmospheric flow input decks
       ├── ABL/
       ├── Bomex/
       ├── SquallLine/
       ├── SuperCell/
       └── ...

Each problem directory contains a README describing its purpose and functionality.

Building with GNU Make
----------------------

System Overview
~~~~~~~~~~~~~~~

The GNU Make system provides a direct path to producing a single, non-case-specific executable that works for all cases
represented in ``ERF/Exec/RegTests`` and ``ERF/Exec/CanonicalTests``.

Depending on the cases of interest, the user must specify in ``ERF/Exec/GNUmakefile`` whether the executable should be built
with particle capability (``USE\_PARTICLES = TRUE``), with NetCDF capability (``USE\_NETCDF = TRUE``),
with radiation (``USE\_RRGMTP = TRUE``), and several other compile-time options.

The build is orchestrated by a ``GNUmakefile`` in ``ERF/Exec/`` (and similarly in ``ERF/.Exec_dev`` for development tests), which uses build logic from the AMReX framework.

.. dropdown:: How it Works: The Orchestration Process
   :icon: info

   The GNU Make process uses a hierarchy of includes separating user configuration from application and framework build logic:

   1. **GNUmakefile Location**: User invokes ``make`` in ``ERF/Exec/`` for standard builds, or in ``ERF/.Exec_dev/<test>/`` for development-test builds.

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

See :ref:`Submodules` for more details. All dependencies except EAMxx SHOC are provided as git submodules:

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

**2. Setup for EAMxx SHOC (Optional)**

If building with EAMxx SHOC, run the setup scripts:

.. code-block:: bash

   export ERF_DIR=/path/to/ERF
   source /path/to/ERF/Build/GNU_Ekat/eamxx_clone.sh
   source /path/to/ERF/Build/GNU_Ekat/ekat_build_commands.sh

Native SHOC is always built in tree and selected at runtime with ``erf.pbl_type = NATIVE_SHOC``.
Set ``USE_EAMXX_SHOC=TRUE`` in your GNUmakefile (step 4) if you need that optional package.

**3. Navigate to GNU Make Build Directory**

.. code-block:: bash

   # Standard workflow
   cd ERF/Exec

   # Development-test workflow (exception)
   # cd ERF/.Exec_dev/<test_name>

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
      * - ``USE_FIRE``
        - Enables the fire coupling (atmosphere side; fire engines plug into it, see :ref:`sec:FireCoupling`)
        - FALSE
        - TRUE/FALSE
      * - ``FIRE_HOME``
        - Fire engine directories to compile in, space separated, each with a ``Make.package`` (needs ``USE_FIRE=TRUE``; for an engine repository that holds ERF as a submodule, see :ref:`sec:FireCoupling`)
        - (none)
        - paths
      * - ``USE_MOVING_BODIES``
        - Enables the actuator core for bodies in the flow (velocity sampling at points, force spreading onto the momentum sources; see :ref:`sec:ActuatorCore`) and the structure coupling (see :ref:`sec:StructureCoupling`)
        - FALSE
        - TRUE/FALSE
      * - ``STRUCTURE_HOME``
        - Structure model directories to compile in, space separated, each with a ``Make.package`` (needs ``USE_MOVING_BODIES=TRUE``; see :ref:`sec:StructureCoupling`)
        - (none)
        - paths
      * - ``USE_DUST``
        - Enables the dust coupling (atmosphere side; dust models plug into it, see :ref:`sec:DustCoupling`); adds the dust scalars to the state
        - FALSE
        - TRUE/FALSE
      * - ``NUM_DUST_SCALARS``
        - Number of dust scalars in the state (needs ``USE_DUST=TRUE``)
        - 1
        - positive integer
      * - ``DUST_HOME``
        - Dust model directories to compile in, space separated, each with a ``Make.package`` (needs ``USE_DUST=TRUE``; see :ref:`sec:DustCoupling`)
        - (none)
        - paths
      * - ``USE_NOAHMP``
        - Enables Noah-MP land surface model (requires ``USE_NETCDF=TRUE``)
        - FALSE
        - TRUE/FALSE
      * - ``USE_RRTMGP``
        - Enables RRTMGP radiation model (sets ``USE_KOKKOS=TRUE``, ``USE_NETCDF=TRUE``)
        - FALSE
        - TRUE/FALSE
      * - Native SHOC
        - Builds the in-tree SHOC path and selects it at runtime with ``erf.pbl_type = NATIVE_SHOC``
        - Always built
        - N/A
      * - ``USE_EAMXX_SHOC``
        - Enables optional EAMxx SHOC turbulence model (sets ``USE_KOKKOS=TRUE``)
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
      * - ``USE_WSM6_FORT``
        - Enables Fortran bridge for WSM6 microphysics (requires ``PRECISION=DOUBLE``)
        - FALSE
        - TRUE/FALSE
      * - ``USE_WDM6_FORT``
        - Enables Fortran bridge for WDM6 microphysics (requires ``PRECISION=DOUBLE``)
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
      **At most one of USE_OMP, USE_CUDA, USE_HIP, USE_SYCL may be TRUE.**
      Setting more than one is rejected by the build with an error naming the
      backends that were selected, rather than silently configuring one of them.

   For additional compiler options, see the `AMReX documentation <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html>`_.

.. dropdown:: Example GNUmakefile
   :icon: code

   Typical ``GNUmakefile`` examples:

   **Exec/GNUmakefile:**

   .. literalinclude:: ../../Exec/GNUmakefile
      :language: makefile

**5. Build**

.. code-block:: bash

   make

The executable name encodes build characteristics (dimensionality, compiler, parallelization). For example, in ``Exec`` with ``COMP=gnu`` and ``USE_MPI=TRUE``, the executable is ``ERF3d.gnu.MPI.ex``. Multiple build configurations can coexist in the same directory.

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
      cd Build/Exec
      mpiexec -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most

      # If out-of-source build (without install):
      cd build/Exec
      mpiexec -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most

      # If installed to install/:
      cd install/bin
      mpiexec -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most

      # If using cmake_with_kokkos_many.sh with defaults:
      cd install/bin  # or cd $ERF_INSTALL_DIR/bin if customized or cd $ERF_BUILD_DIR/Exec
      mpiexec -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most

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

**2. Setup for EAMxx SHOC (Optional)**

If building with EAMxx SHOC:

.. code-block:: bash

   export ERF_DIR=/path/to/ERF
   source /path/to/ERF/Build/GNU_Ekat/eamxx_clone.sh

Native SHOC is always built in tree; select it at runtime with ``erf.pbl_type = NATIVE_SHOC``.
Configure with ``-DERF_ENABLE_EAMXX_SHOC=ON`` (step 4) if you need that optional package.

**3. Choose Build Workflow**

ERF supports multiple CMake workflows. The main difference is directory structure and whether you use the install step.

.. tab-set::

   .. tab-item:: In Build/ Directory

      Build artifacts and executables in the Build/ directory tree. Executables appear in ``Exec/<problem>/`` subdirectories, similar to GNU Make workflow.

      .. code-block:: bash

         cd Build
         ./cmake.sh

      **Executable locations:** ``Build/Exec/erf_exec``

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

      **Executable locations:** ``build/Exec/erf_exec``, etc., and optionally ``install/bin/erf_exec`` (if installed)

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

      **Executable locations:** ``$ERF_BUILD_DIR/Exec/erf_exec`` and ``$ERF_INSTALL_DIR/bin/erf_exec`` (defaults: ``./Exec/erf_exec`` and ``install/bin/erf_exec``)

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
   * - ``ERF_ENABLE_FIRE``
     - Enables the fire coupling (atmosphere side; fire engines plug into it, see :ref:`sec:FireCoupling`); adds the smoke scalar to the state
     - OFF
     - ON/OFF
   * - ``ERF_FIRE_ENGINE_DIRS``
     - Fire engine directories to compile in, a CMake list, each with an ``ERFFireEngine.cmake`` (needs ``ERF_ENABLE_FIRE=ON``; set by an engine repository that holds ERF as a submodule, and left unset when that repository uses ERF as a library; see :ref:`sec:FireCoupling`)
     - (none)
     - paths
   * - ``ERF_ENABLE_MOVING_BODIES``
     - Enables the actuator core for bodies in the flow (velocity sampling at points, force spreading onto the momentum sources; see :ref:`sec:ActuatorCore`) and the structure coupling (see :ref:`sec:StructureCoupling`)
     - OFF
     - ON/OFF
   * - ``ERF_STRUCTURE_MODEL_DIRS``
     - Structure model directories to compile in, a CMake list, each with an ``ERFStructureModel.cmake`` (needs ``ERF_ENABLE_MOVING_BODIES=ON``; see :ref:`sec:StructureCoupling`)
     - (none)
     - paths
   * - ``ERF_ENABLE_DUST``
     - Enables the dust coupling (atmosphere side; dust models plug into it, see :ref:`sec:DustCoupling`); adds the dust scalars to the state
     - OFF
     - ON/OFF
   * - ``ERF_NUM_DUST_SCALARS``
     - Number of dust scalars in the state (needs ``ERF_ENABLE_DUST=ON``)
     - 1
     - positive integer
   * - ``ERF_DUST_MODEL_DIRS``
     - Dust model directories to compile in, a CMake list, each with an ``ERFDustModel.cmake`` (needs ``ERF_ENABLE_DUST=ON``; see :ref:`sec:DustCoupling`)
     - (none)
     - paths
   * - ``ERF_PRECISION``
     - Floating-point precision for mesh data
     - DOUBLE
     - SINGLE/DOUBLE
   * - ``ERF_PARTICLES_PRECISION``
     - Floating-point precision for particle data; defaults to ``ERF_PRECISION``
     - DOUBLE
     - SINGLE/DOUBLE
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
   * - Native SHOC (no option)
     - Builds the in-tree SHOC path and selects it at runtime with ``erf.pbl_type = NATIVE_SHOC``
     - Always built
     - N/A
   * - ``ERF_ENABLE_EAMXX_SHOC``
     - Enables optional EAMxx SHOC turbulence model (requires ``ERF_ENABLE_MPI=ON`` and EKAT/Kokkos/EAMxx dependencies)
     - OFF
     - ON/OFF
   * - ``ERF_ENABLE_WSM6_FORT``
     - Enables Fortran bridge for WSM6 microphysics (requires ``ERF_PRECISION=DOUBLE``)
     - OFF
     - ON/OFF
   * - ``ERF_ENABLE_WDM6_FORT``
     - Enables Fortran bridge for WDM6 microphysics (requires ``ERF_PRECISION=DOUBLE``)
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
   **At most one of ERF_ENABLE_OPENMP, ERF_ENABLE_CUDA, ERF_ENABLE_HIP, ERF_ENABLE_SYCL may be ON.**
   Setting more than one is rejected at configure time with an error naming the
   backends that were selected, rather than silently configuring one of them.

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

   * Native SHOC - Native SHOC turbulence and cloud macrophysics

     - Always built
     - Selected at runtime with ``erf.pbl_type = NATIVE_SHOC``

   * ``ERF_ENABLE_EAMXX_SHOC`` - Optional EAMxx SHOC turbulence and cloud macrophysics

     - Requires ``ERF_ENABLE_MPI=ON``
     - Automatically enables ``ERF_ENABLE_EKAT=ON`` (provides Kokkos)
     - Additional step: Run ``source Build/GNU_Ekat/eamxx_clone.sh``

   * ``ERF_ENABLE_NOAHMP`` - NOAHMP land surface model

     - Requires the NetCDF Fortran library

   **GPU Acceleration**

   Enable exactly one GPU backend:

   * ``ERF_ENABLE_CUDA`` - NVIDIA GPUs (requires CUDA Toolkit ≥ 11.0)
   * ``ERF_ENABLE_HIP`` - AMD GPUs
   * ``ERF_ENABLE_SYCL`` - Intel GPUs

   .. note::
      Kokkos-based physics packages (RRTMGP, EAMxx SHOC) support all three GPU backends through EKAT's Kokkos integration.

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
   - ``aurora_erf.profile`` - ALCF Aurora (Intel GPUs); expects ``NETCDF_DIR`` for NetCDF-enabled builds

   These profiles are shell scripts that load required modules and set environment variables. Source them before running CMake. For detailed information, see :ref:`sec:build:hpc`.

.. dropdown:: Developer Testing
   :icon: beaker

   For systematic testing of multiple build configurations:

   .. literalinclude:: ../../Build/setup_cmake_validation.sh
      :language: bash

   **Usage:**

   .. code-block:: bash

      cd Build
      ./setup_cmake_validation.sh default    # or perlmutter, gnu_ekat
      cd ../build_default
      ./run.sh     # List available scripts
      ./run.sh 1   # Run first script

   Creates isolated build directories for each configuration script. Each script runs in its own subdirectory (``script_<name>/``) at the ERF root.

.. _sec:build:parent_project:

Building ERF from a Parent CMake Project
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A repository that holds a fire engine, dust model or structure model keeps ERF as a git
submodule together with the model's folders. Such a repository can build ERF in three
ways. The third, :ref:`sec:ERFAsLibrary`, is the recommended one: the
repository builds its own executable and hands ERF no model directories.
`FLARE <https://github.com/erf-model/FLARE>`_ builds this way.

The first configures ERF itself and points the hook variables at the outer folders.
`LSFIRE <https://github.com/erf-model/LSFIRE>`_ does this from its root, with ERF in
``ERF/`` and its fire engine in ``FireEngine/``:

.. code-block:: bash

   cmake -S ERF -B build -DERF_ENABLE_FIRE=ON -DERF_FIRE_ENGINE_DIRS=$PWD/FireEngine
   cmake --build build --target erf_exec

The second makes ERF part of the repository's own CMake project, for a repository that
builds more than ``erf_exec`` or adds its own tests. The parent calls ``add_subdirectory``:

.. code-block:: cmake

   cmake_minimum_required(VERSION 3.25)   # as ERF's internal AMReX
   project(MyModel CXX C)

   set(ERF_ENABLE_MOVING_BODIES ON)
   set(ERF_STRUCTURE_MODEL_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/Model")

   enable_testing()                  # in the top-level file, so ctest runs from the top build folder
   add_subdirectory(Submodules/ERF)  # builds erf_exec with the model compiled in

- ERF options may be set as normal variables before ``add_subdirectory``, as above, or
  as cache variables (``-D`` on the command line). The options ERF sets with
  ``set(... CACHE ...)``, namely ``ERF_STRUCTURE_MODEL_DIRS``, ``ERF_FIRE_ENGINE_DIRS``,
  ``ERF_DUST_MODEL_DIRS``, ``ERF_NUM_DUST_SCALARS``, ``ERF_PRECISION``,
  ``ERF_PARTICLES_PRECISION`` and ``ERF_DIM``, keep a normal variable
  from the parent only with CMake 3.21 or newer (policy CMP0126). This matters only with
  an external AMReX (``ERF_USE_INTERNAL_AMREX=OFF``), since the internal AMReX needs
  CMake 3.25. ``option()`` keys such as ``ERF_ENABLE_MOVING_BODIES`` work either way.
- Give ``ERF_FIRE_ENGINE_DIRS``, ``ERF_STRUCTURE_MODEL_DIRS`` and ``ERF_DUST_MODEL_DIRS`` absolute paths.
- ERF's paths refer to its own source and build folders, so ``erf_exec`` is in
  ``<build>/<ERF binary dir>/Exec`` and ERF's tests (with ``ERF_ENABLE_TESTS=ON``) find it there.
  The parent can refer to the executable as ``$<TARGET_FILE:erf_exec>``.
- Under a parent project ERF adds no ``uninstall`` or ``distclean`` targets and does not
  export its CMake package; those act on the top-level build and belong to the parent.
  ``cmake --install`` on the parent still installs ``erf_exec``, ERF's library and AMReX.
- ERF sets the AMReX options (``AMReX_PRECISION``, ``AMReX_SPACEDIM`` and others) from
  its own; a parent should set ERF's options, not the ``AMReX_*`` ones.
- On Cray systems ERF's compiler detection is skipped under a parent, because the
  parent's ``project()`` has already chosen the compilers.

The parent can register its own regression cases with ERF's general CTest helpers,
which ``Tests/ERFTestHelpers.cmake`` holds: ``add_test_r`` (run a deck and compare the
plotfile with a gold plotfile by ``fcompare``), ``add_test_abort`` (pass when a run stops
with a given message), ``add_test_restart_abort``, ``add_test_plotfile_header``, and the
``setup_test`` and ``resolve_test_exe`` they use. ERF's own ``Tests/CTestList.cmake``
includes the same file.

.. code-block:: cmake

   add_subdirectory(Submodules/ERF)
   set(ERF_TEST_GOLD_FILES_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/GoldFiles")
   include(Submodules/ERF/Tests/ERFTestHelpers.cmake)
   add_test_r(MyCase "" erf_exec plt00010)   # deck test_files/MyCase/MyCase.i

- Include the file after ``add_subdirectory``, in each directory that calls the helpers.
- A test's files are read from ``test_files/<test name>`` of the directory that calls the
  helper, and ``add_test_r`` compares with ``${ERF_TEST_GOLD_FILES_DIRECTORY}/<test name>``.
  A normal variable ``ERF_TEST_GOLD_FILES_DIRECTORY``, set after ``add_subdirectory`` as
  above, gives the parent its own gold folder without changing ERF's. Set before
  ``add_subdirectory``, it would also become the gold folder of ERF's own tests.
- ``add_test_r`` needs ``FCOMPARE_EXE``, ``ERF_TEST_GOLD_FILES_DIRECTORY``,
  ``ERF_TEST_FCOMPARE_RTOL`` and ``ERF_TEST_FCOMPARE_ATOL``, and ``ERF_TEST_NRANKS`` with
  ``ERF_ENABLE_MPI``; it stops the configure, naming the variable, when one is empty.
  ``ERF_ENABLE_FCOMPARE=ON`` builds ``fcompare``, and ``ERF_ENABLE_TESTS=ON`` sets the
  tolerances, the rank count and ERF's gold folder as cache entries. Without ERF's tests,
  the parent sets them itself as normal variables.
- A caller's normal variable ``TEST_FILES_DIR`` replaces ``<test name>`` in the folder the
  inputs are read from, as ERF's own list uses it; leave it unset otherwise.
- The helpers find ``erf_exec`` and ERF's checker scripts through ``ERF_BINARY_DIR`` and
  ``ERF_SOURCE_DIR``, so they work from any directory of the parent.

``Tests/ParentProject`` holds a small parent project that does this with the example
structure model (``Tests/StructureModelExample``) copied to a folder outside the ERF
tree. It has a run checked against the closed form; three cases registered with
``add_test_abort`` and ``add_test_r`` (gold in the parent's own folder), one of which must
fail, to show that ``fcompare`` compares; and, with ``ERF_ENABLE_UNIT_TESTS=ON``, a run of
the model's GoogleTest added through ``erf_structure_model_add_unit_tests``.
The ``parent-project`` CI workflow builds it and runs the tests; its ``README.md`` gives
the commands.

.. _sec:ERFAsLibrary:

ERF as a library
^^^^^^^^^^^^^^^^

A third way uses ERF as ERF uses AMReX. ERF is a library, ``ERF::erf``, and the parent
builds its own executables from its own sources (its engines or models and its own
``ERF_Prob.cpp``) with ``erf_setup_executable``. The parent hands ERF none of the
model-directory hooks: its sources are compiled into its own executable, and their
registrations (``ERF_REGISTER_FIRE_ENGINE``, ``ERF_REGISTER_DUST_MODEL``,
``ERF_REGISTER_STRUCTURE_MODEL``) run from there.

.. code-block:: cmake

   cmake_minimum_required(VERSION 3.25)
   project(MyModel CXX C)

   set(ERF_ENABLE_MOVING_BODIES ON)   # the couplings the parent's sources need

   enable_testing()
   add_subdirectory(Submodules/ERF ERF EXCLUDE_FROM_ALL)   # build only what my_model needs

   add_executable(my_model Model/MyStructureModel.cpp ERF_Prob.cpp)
   target_include_directories(my_model PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
   erf_setup_executable(my_model)   # ERF's main.cpp, ERF::erf, ERF's compile and GPU flags

- ``erf_setup_executable(<target> [NO_ERF_MAIN])`` adds ERF's ``main.cpp`` (unless
  ``NO_ERF_MAIN``, for a GoogleTest executable, say), links ``ERF::erf`` and sets ERF's
  compile flags, and the CUDA settings with ``ERF_ENABLE_CUDA``. ERF's own executables are
  set up the same way (``build_erf_exe``, the same function under its older name).
- ``ERF::erf`` carries ERF's compile definitions (``ERF_USE_FIRE``, ``ERF_USE_DUST``,
  ``ERF_NUM_DUST_SCALARS`` and the others) and include folders to the parent's sources, as
  ``AMReX::amrex`` carries ``AMREX_SPACEDIM`` to ERF's. Set ERF's options before
  ``add_subdirectory``, so that both sides see the same state layout.
- The parent's problem file plays the part of ERF's ``Exec/ERF_Prob.cpp``; a copy of
  ERF's is a fine start.
- ERF's test helpers take the parent's executable target: ``add_test_r(MyCase "" my_model
  plt00010)``, and a normal variable ``ERF_TEST_EXE`` names it for ``add_test_abort`` and
  ``add_test_restart_abort``.
- ``ERF_BUILD_LIBRARY_ONLY`` is a different mode, for building ERF together with another
  AMReX code (REMORA) into one executable; it renames some of ERF's symbols and compiles
  ERF's own problem file into the library. A parent that builds an ERF executable does not
  set it.

``Tests/ParentExecutable`` holds a small parent project that does this: its executable is
built from a copy of the example structure model and of ERF's problem file, kept outside
the ERF tree. It checks a run against the closed form, and runs ``add_test_abort`` on the
parent's executable. The ``parent-project`` CI workflow builds and runs it too.

Building Documentation
~~~~~~~~~~~~~~~~~~~~~~

Build ERF documentation with CMake:

.. code-block:: bash

   cmake -DERF_ENABLE_DOCUMENTATION=ON ..
   make docs

This builds both Sphinx and Doxygen documentation.
