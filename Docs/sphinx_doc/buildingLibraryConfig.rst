.. _sec:build:library:

=================================
A Developer's Guide to Configuring ERF Library Dependencies
=================================

1.0 Overview of ERF's External Library Ecosystem
------------------------------------------------

**1.1 Introduction: The Importance of Library Configuration**

Configuring external libraries enables ERF's scientific capabilities, from parallel I/O on supercomputing platforms to GPU-accelerated physics packages. A correctly configured dependency stack is the foundation for ERF's advanced features. The ERF build system automates much of this complexity, but understanding the underlying components and their interactions is essential for customizing builds and diagnosing issues. This guide provides a technical reference for developers to navigate the ERF library ecosystem, ensuring a robust and reproducible path from source code to executable.

**1.2 Analysis of Key Library Dependencies**

The ERF framework integrates several external libraries to provide its core functionality and advanced features. The following table summarizes these primary dependencies, their purpose, and the primary build options used to control their integration.

.. list-table:: ERF Library Dependencies
   :header-rows: 1
   :widths: 20 40 10 15 15

   * - Library
     - Description
     - Default
     - Possible Values
   * - AMReX
     - Core framework providing block-structured adaptive mesh refinement (AMR), data structures, and solvers
     - Required
     - ``-DERF_USE_INTERNAL_AMREX=ON`` / ``AMREX_HOME``
   * - NetCDF
     - High-level I/O operations, including writing plotfiles and reading initial conditions from WRF or Metgrid files
     - Optional
     - ``-DERF_ENABLE_NETCDF=ON`` / ``USE_NETCDF=TRUE``
   * - HDF5
     - Low-level backend for high-performance parallel I/O, used by NetCDF library for distributed-memory file access
     - Optional
     - ``-DERF_ENABLE_HDF5=ON`` / Via NetCDF pkg-config
   * - SHOC
     - Simplified Higher-Order Closure turbulence and cloud macrophysics model from E3SM project
     - Optional
     - ``-DERF_ENABLE_SHOC=ON`` / ``USE_SHOC=TRUE``
   * - ZFP
     - Library for high-performance, lossy data compression, usable as optional HDF5 filter
     - Optional
     - ``-DAMReX_HDF5_ZFP=ON`` / Not directly supported

**1.3 The Requirement for Parallel I/O Consistency**

For any ERF simulation leveraging Message Passing Interface (MPI) for parallel execution, there is a requirement for consistent parallel I/O stack. When parallel I/O is enabled (e.g., via ``-DERF_ENABLE_NETCDF=ON``), both NetCDF library and its underlying HDF5 dependency must be parallel-enabled versions. These specialized builds contain necessary communication hooks to coordinate I/O operations across multiple compute nodes. Attempting to link parallel ERF build against serial versions of these libraries results in linker errors or runtime failures during file operations.

2.0 Configuring Core I/O Libraries: NetCDF and HDF5
---------------------------------------------------

**2.1 Introduction: The Foundation of ERF's Data Strategy**

NetCDF and HDF5 form the foundation of ERF's I/O strategy, particularly for large-scale, data-intensive simulations on HPC platforms. They provide the mechanism for reading complex initial conditions and writing massive datasets in a scalable manner. However, their configuration requires specific, parallel-aware builds that must be correctly matched to the MPI implementation and the rest of the software stack. This section provides the technical details to manage this process successfully.

**2.2 Delineating Library Roles and Relationships**

ERF leverages NetCDF and HDF5 for distinct but complementary roles in its I/O workflow:

* **NetCDF**: Provides the high-level Application Programming Interface (API) that ERF source files use to interact with structured, self-describing data formats. When enabled, the build system compiles source files responsible for I/O, such as ``ERF_InitFromWRFInput.cpp`` for reading initial conditions and ``ERF_NCPlotFile.cpp`` for writing standard AMReX plotfiles in NetCDF format.

* **HDF5**: Serves as the low-level, high-performance parallel storage layer that NetCDF relies on in distributed-memory environments. While ERF source code primarily calls NetCDF functions, HDF5 manages the task of coordinating data writes from many MPI ranks into a single file on a parallel filesystem.

**2.3 Evaluating the CMake Detection Hierarchy**

ERF's CMake build system uses a hierarchical search process to locate the NetCDF library. This strategy prioritizes explicit user control over environmental hints and system defaults.

.. list-table:: NetCDF Detection Priority
   :header-rows: 1
   :widths: 10 25 25 40

   * - Priority
     - Method
     - Variables Checked
     - Notes
   * - 1 (Highest)
     - User-Provided CMake Options
     - ``NETCDF_DIR``
     - Passed via ``-DNETCDF_DIR=...``, overrides all other mechanisms
   * - 2 (High)
     - Environment Variables
     - ``NETCDF_DIR``
     - Sourced from shell environment via ``export NETCDF_DIR=...``
   * - 3 (Medium)
     - pkg-config Utility
     - Searches: ``netcdf``, ``netcdf-mpi``, ``netcdf_parallel``, ``netcdf-cxx4_parallel``
     - Standard detection method on Linux systems
   * - 4 (Low)
     - System Default Paths
     - ``/usr/lib``, ``/usr/include``
     - Fallback for system-wide installations

**2.4 Cross-Referencing Build System Configuration Variables**

The following tables detail the primary options and variables used to control NetCDF and HDF5 features.

.. list-table:: CMake Options for I/O Libraries
   :header-rows: 1
   :widths: 30 50 20

   * - CMake Option
     - Description
     - Example
   * - ``ERF_ENABLE_NETCDF``
     - Master switch to enable all NetCDF-related features and I/O routines
     - ``-DERF_ENABLE_NETCDF=ON``
   * - ``ERF_ENABLE_HDF5``
     - Enables HDF5 support within AMReX backend (defaults to ``ERF_ENABLE_NETCDF`` value)
     - ``-DERF_ENABLE_HDF5=ON``
   * - ``NETCDF_DIR``
     - Specifies installation prefix for custom NetCDF library
     - ``-DNETCDF_DIR=/path/to/netcdf``
   * - ``HDF5_PREFER_PARALLEL``
     - Hint for find module to prioritize parallel HDF5 libraries
     - ``-DHDF5_PREFER_PARALLEL=ON``

.. list-table:: GNU Make Variables for I/O Libraries
   :header-rows: 1
   :widths: 30 50 20

   * - GNU Make Variable
     - Description
     - Example
   * - ``USE_NETCDF``
     - Master switch to enable NetCDF features and link necessary libraries
     - ``make USE_NETCDF=TRUE``

**2.5 Analyzing the MPI Dependency Chain and Common Pitfalls**

A correct parallel I/O setup requires a complete and consistent parallel I/O dependency chain. An MPI-enabled ERF build must link to parallel-enabled NetCDF library, which must have been compiled against parallel-enabled HDF5 library. The build system enforces this by using pkg-config to search for explicitly parallel versions of libraries before falling back to standard packages.

The most common failure mode is serial/parallel library mismatch. Linking parallel ERF build against serial versions of NetCDF or HDF5 leads to linker errors reporting undefined references to MPI symbols (e.g., ``MPI_Init``, ``MPI_Comm_rank``, ``MPI_Finalize``). These errors occur because serial libraries lack necessary MPI function calls that parallel ERF executable expects.

**2.6 Platform-Specific Configuration Strategies**

The recommended method for configuring NetCDF and HDF5 varies by platform:

* **Cray Systems (Perlmutter, Frontier, Polaris)**: Load the ``cray-netcdf-hdf5parallel`` environment module. This provides correctly configured, MPI-enabled versions of both libraries. ERF's integrated Cray detection logic automatically queries the Cray compiler wrapper to discover module-specific search paths.

* **Workstations and Generic Clusters**: Use the system's package manager or environment module system to install parallel versions of NetCDF and HDF5. Ensure MPI implementation is installed first, then install I/O libraries with MPI support enabled.

**2.7 Common Configuration Scenarios**

The following examples demonstrate common configuration commands:

.. code:: shell

   # CMake: Enabling NetCDF with custom installation path
   cmake -DERF_ENABLE_NETCDF=ON -DNETCDF_DIR=/opt/netcdf-parallel ..

   # CMake: Explicitly disabling NetCDF
   cmake -DERF_ENABLE_NETCDF=OFF ..

   # GNU Make: Enabling NetCDF
   make USE_NETCDF=TRUE

3.0 Integrating Advanced Physics and Utility Libraries
------------------------------------------------------

**3.1 Introduction: Extending ERF's Scientific Capabilities**

Beyond core I/O functions, ERF integrates specialized libraries to extend its scientific modeling and data handling capabilities. These libraries provide access to advanced physics packages and data compression schemes for simulations. Most components are included as Git submodules and enabled with simple build options.

**3.2 Library Integration Details**

* **SHOC (Simplified Higher-Order Closure)**: Turbulence and cloud macrophysics model from E3SM project. Provides schemes for modeling atmospheric boundary layer processes. Depends on Kokkos performance portability library (provided via EKAT submodule) for efficient execution on CPUs and GPUs.

  - CMake Option: ``-DERF_ENABLE_SHOC=ON``
  - GNU Make Variable: ``USE_SHOC=TRUE``

* **ZFP (Floating-Point Compression)**: Library for high-performance data compression used as optional filter within HDF5 library. Reduces output data size in large-scale simulations. Provides lossy compression - use must be evaluated against simulation's scientific goals. Requires HDF5 library compiled with ZFP plugin.

  - CMake Option: ``-DAMReX_HDF5_ZFP=ON``

* **AMReX (Adaptive Mesh Refinement Framework)**: Required dependency providing core AMR data structures, solvers, and parallelization infrastructure. Default approach uses internal Git submodule for version compatibility. External build option available for environments requiring shared AMReX installation across multiple applications.

  - CMake Option: ``-DERF_USE_INTERNAL_AMREX=ON`` (Default)
  - GNU Make Variable: ``AMREX_HOME`` points to submodule directory

4.0 Troubleshooting Common Library Configuration Issues
-------------------------------------------------------

**4.1 Introduction: A Systematic Approach to Diagnosis**

Library configuration issues are common sources of build problems, particularly on complex HPC systems with module-based environments. These errors typically arise from mismatches between user's shell environment and build system expectations. This section provides systematic guidance for diagnosing and resolving frequent errors.

**4.2 Diagnosing and Resolving Common Errors**

**Detection Failures ("Library Not Found")**

This error occurs when build system cannot locate necessary header files or library files for dependency like NetCDF or HDF5.

*Solution:*

1. Verify Environment Modules: Use ``module list`` to confirm correct modules are loaded. On Cray systems, ensure ``cray-netcdf-hdf5parallel`` is present.
2. Specify Paths Manually: For non-standard locations, provide explicit path to installation prefix (e.g., ``cmake -DNETCDF_DIR=/path/to/netcdf ..``)
3. Check pkg-config: Ensure ``pkg-config`` utility is in PATH and ``PKG_CONFIG_PATH`` environment variable includes directory containing library's ``.pc`` file.

**MPI Linker Errors ("Undefined Reference to MPI_*")**

This error occurs during final linking stage when linker cannot find implementations for core MPI functions.

*Solution:*

1. Diagnosis: Indicates mismatch where parallel (MPI-enabled) ERF build is being linked against serial version of dependency (see :ref:`sec:build:library` Section 2.5).
2. Resolution: Unload serial library modules and load correct parallel-enabled versions. Perform clean build to ensure new configuration is applied.

**Stale CMake Cache (Unexpected Failures After Environment Changes)**

This issue manifests as unexpected build failures after changing environment (loading different module or switching compiler versions).

*Solution:*

1. Diagnosis: CMake caches library locations to speed reconfiguration. If environment changes, cache becomes stale.
2. Resolution: Perform clean build using ``make distclean`` or ``cmake --build . --target distclean``. Alternatively, manually remove cache files (``CMakeCache.txt``, ``CMakeFiles/``) from build directory.

**4.3 Advanced Diagnosis Using Build Logs**

For complex issues, CMake's verbose logging capabilities provide invaluable diagnostics. These logs reveal exact search paths, environment variables, and logic used by ``find_package`` commands.

Generate detailed log showing hierarchy and context of configuration messages:

.. code:: shell

   # Generate verbose, context-aware log of CMake configuration
   cmake --log-level=VERBOSE --log-context ..

This provides structured view of configuration process, allowing precise identification of configuration failures:

.. code:: text

   [ERF.Cray] Detected Cray Programming Environment (CRAYPE_VERSION=2.7.30)
   [ERF.Cray] Setting Cray compiler wrappers...
   [ERF.Cray]   Set CMAKE_CXX_COMPILER = /opt/cray/pe/craype/default/bin/CC
   [ERF.AMReX] Using internal AMReX submodule
   [ERF.NetCDF] Found NetCDF: /opt/cray/pe/netcdf/4.9.0.9
