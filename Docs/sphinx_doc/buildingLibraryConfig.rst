A Developer's Guide to Configuring ERF Library Dependencies

1. Overview of ERF's External Library Ecosystem

1.1. Introduction: The Strategic Importance of Library Configuration

Configuring external libraries is a critical, high-stakes process for building the Energy Research and Forecasting (ERF) model. A correctly configured dependency stack is the foundation that enables ERF's advanced scientific capabilities, from high-performance parallel I/O on supercomputing platforms to sophisticated, GPU-accelerated physics packages. The ERF build system is engineered to automate much of this complexity, but a clear understanding of the underlying components and their interactions is essential for customizing builds and diagnosing issues. This guide provides a definitive technical reference for developers to navigate the ERF library ecosystem, ensuring a robust and reproducible path from source code to executable. The following sections detail the key libraries and their configuration.

1.2. Analysis of Key Library Dependencies

The ERF framework integrates several external libraries to provide its core functionality and advanced features. The following table summarizes these primary dependencies, their purpose, and the primary build options used to control their integration in both CMake and GNU Make build environments.

Library	Purpose	Required?	Primary CMake Option	Primary GNU Make Variable
AMReX	The core framework providing block-structured adaptive mesh refinement (AMR), data structures, and solvers.	Yes	-DERF_USE_INTERNAL_AMREX=ON	AMREX_HOME
NetCDF	High-level I/O operations, including writing plotfiles and reading initial conditions from WRF or Metgrid files.	Optional	-DERF_ENABLE_NETCDF=ON	USE_NETCDF=TRUE
HDF5	A low-level backend for high-performance parallel I/O, used by the NetCDF library for distributed-memory file access.	Optional	-DERF_ENABLE_HDF5=ON	Indirectly via NetCDF pkg-config
SHOC	The Simplified Higher-Order Closure turbulence and cloud macrophysics model from the E3SM project.	Optional	-DERF_ENABLE_SHOC=ON	USE_SHOC=TRUE
ZFP	A library for high-performance, lossy data compression, usable as an optional HDF5 filter to reduce data volume.	Optional	-DAMReX_HDF5_ZFP=ON	(Not directly supported)

1.3. The Critical Requirement for Parallel I/O Consistency

For any ERF simulation leveraging the Message Passing Interface (MPI) for parallel execution, there is a non-negotiable requirement for a consistent parallel I/O stack. When parallel I/O is enabled (e.g., via -DERF_ENABLE_NETCDF=ON), both the NetCDF library and its underlying HDF5 dependency must be the parallel-enabled versions. These specialized builds contain the necessary communication hooks to coordinate I/O operations across multiple compute nodes. Attempting to link a parallel ERF build against serial versions of these libraries is a common point of failure, typically resulting in cryptic linker errors or immediate runtime failures during file operations.

2. Configuring Core I/O Libraries: NetCDF and HDF5

2.1. Introduction: The Backbone of ERF's Data Strategy

NetCDF and HDF5 are the cornerstones of ERF's I/O strategy, particularly for large-scale, data-intensive simulations on High-Performance Computing (HPC) platforms. They provide the mechanism for reading complex initial conditions and writing massive datasets in a scalable, efficient manner. However, their configuration is often challenging due to the strict requirement for specific, parallel-aware builds that must be correctly matched to the MPI implementation and the rest of the software stack. This section provides the necessary technical details to manage this process successfully.

2.2. Delineating Library Roles and Relationships

ERF leverages NetCDF and HDF5 for distinct but complementary roles in its I/O workflow:

* NetCDF: This library provides the high-level Application Programming Interface (API) that ERF source files use to interact with structured, self-describing data formats. When enabled, the build system compiles source files responsible for I/O, such as ERF_InitFromWRFInput.cpp for reading initial conditions and ERF_NCPlotFile.cpp for writing standard AMReX plotfiles in NetCDF format.
* HDF5: This is the low-level, high-performance parallel storage layer that NetCDF relies on in distributed-memory environments. While ERF source code primarily calls NetCDF functions, HDF5 manages the complex task of coordinating data writes from many MPI ranks into a single, coherent file on a parallel filesystem, typically by leveraging the underlying MPI-IO library.

2.3. Evaluating the CMake Detection Hierarchy

ERF's CMake build system uses a hierarchical search process to locate the NetCDF library. This strategy is designed to prioritize explicit user control over environmental hints and system defaults, providing a robust and predictable detection mechanism.

Priority	Method	Variables Checked	Notes
1 (Highest)	User-Provided CMake Options	NETCDF_DIR	Passed via -DNETCDF_DIR=..., this setting overrides all other mechanisms, giving the user absolute control. This is the most direct way to point to a custom build.
2 (High)	Environment Variables	NETCDF_DIR	Sourced from the shell environment via export NETCDF_DIR=.... This is the standard method used by environment module systems on HPC platforms to expose specific library installations.
3 (Medium)	pkg-config Utility	Searches for packages in order: netcdf, netcdf-mpi, netcdf_parallel, netcdf-cxx4_parallel.	This is the standard detection method on many Linux systems and is the mechanism used by the automated Cray detection system.
4 (Low)	System Default Paths	Standard paths like /usr/lib and /usr/include.	This is a fallback for system-wide installations typically found on developer workstations.

2.4. Cross-Referencing Build System Configuration Variables

The following tables detail the primary options and variables used to control NetCDF and HDF5 features in both the CMake and GNU Make build systems.

CMake Options

CMake Option	Purpose	Example
ERF_ENABLE_NETCDF	Master switch to enable all NetCDF-related features and I/O routines.	-DERF_ENABLE_NETCDF=ON
ERF_ENABLE_HDF5	Enables HDF5 support within the AMReX backend. This option defaults to the value of ERF_ENABLE_NETCDF and is thus enabled automatically when NetCDF is active.	-DERF_ENABLE_HDF5=ON
NETCDF_DIR	Specifies the installation prefix for a custom NetCDF library, overriding automated detection.	-DNETCDF_DIR=/path/to/netcdf
HDF5_PREFER_PARALLEL	A hint for the find module to prioritize locating parallel HDF5 libraries. This is set automatically on Cray systems for HIP builds to ensure the correct library is chosen.	-DHDF5_PREFER_PARALLEL=ON

GNU Make Variables

GNU Make Variable	Purpose	Example
USE_NETCDF	Master switch to enable NetCDF features and link against the necessary libraries found via pkg-config.	make USE_NETCDF=TRUE

2.5. Analyzing the MPI Dependency Chain and Common Pitfalls

A correct parallel I/O setup requires a complete and consistent parallel I/O dependency chain. An MPI-enabled ERF build must link to a parallel-enabled NetCDF library, which in turn must have been compiled against a parallel-enabled HDF5 library. The build system attempts to enforce this by using pkg-config to search for explicitly parallel versions of the libraries, such as netcdf-mpi, before falling back to standard packages.

The most common failure mode related to I/O libraries is a serial/parallel library mismatch. Linking a parallel ERF build against serial versions of NetCDF or HDF5 will invariably lead to build failures. The typical symptom is a series of linker errors reporting an "undefined reference" to core MPI symbols (e.g., MPI_Init, MPI_Comm_rank, MPI_Finalize). These errors occur because the serial libraries lack the necessary MPI function calls that the parallel ERF executable expects to find.

2.6. Prescribing Platform-Specific Configuration Strategies

The recommended method for configuring NetCDF and HDF5 varies by platform.

* Cray Systems (Perlmutter, Frontier, Polaris): The officially recommended method is to load the cray-netcdf-hdf5parallel environment module. This single module provides correctly configured, MPI-enabled versions of both libraries. This approach works seamlessly because ERF's integrated Cray detection logic automatically queries the Cray compiler wrapper (CC --cray-print-opts=PKG_CONFIG_PATH) to discover the module-specific search path. It then prepends this path to the PKG_CONFIG_PATH environment variable, ensuring that the pkg-config utility finds these specific parallel libraries.
* Workstations and Generic Clusters: On other systems, developers should leverage the system's package manager (e.g., apt, yum) or environment module system to install parallel versions of NetCDF and HDF5. On macOS, for example, Homebrew is a common choice. The correct procedure is to first ensure an MPI implementation is installed (e.g., brew install open-mpi) and then install the I/O libraries (brew install netcdf hdf5). Homebrew's NetCDF formula is typically compiled with MPI support if it detects an existing MPI installation.

2.7. Illustrating Common Configuration Scenarios

The following examples demonstrate common configuration commands for both build systems.

* CMake: Enabling NetCDF with a custom installation path
* CMake: Explicitly disabling NetCDF
* GNU Make: Enabling NetCDF

With the core I/O libraries configured, we can now turn to the optional libraries that extend ERF's physics and data handling capabilities.

3. Integrating Advanced Physics and Utility Libraries

3.1. Introduction: Extending ERF's Scientific Capabilities

Beyond its core I/O functions, ERF integrates several specialized libraries to extend its scientific modeling and data handling capabilities. These libraries provide access to advanced physics packages, such as turbulence models and data compression schemes, which are critical for state-of-the-art simulations. Most of these components are included as Git submodules and are enabled with simple build options, a design choice that streamlines their integration and ensures version compatibility.

3.2. Library Integration Details

* SHOC (Simplified Higher-Order Closure) SHOC is a turbulence and cloud macrophysics model developed as part of the E3SM project. It provides sophisticated schemes for modeling atmospheric boundary layer processes. SHOC depends on the Kokkos performance portability library, which is provided via the EKAT submodule, to run efficiently on both CPUs and GPUs. It is enabled with the following CMake option and GNU Make variable:
  * CMake Option: -DERF_ENABLE_SHOC=ON
  * GNU Make Variable: USE_SHOC=TRUE
* ZFP (Floating-Point Compression) ZFP is a library for high-performance data compression that can be used as an optional filter within the HDF5 library. Its primary function is to reduce the size of output data, which can be substantial in large-scale simulations. Critically, ZFP provides lossy compression, meaning its use must be carefully evaluated against the simulation's scientific goals due to the potential for loss of numerical precision. Enabling this option requires that the HDF5 library itself was compiled with the ZFP plugin, which may not be the case for standard system installations. This feature is managed by the underlying AMReX build, and its configuration is passed through the ERF build command:
  * CMake Option: -DAMReX_HDF5_ZFP=ON
* AMReX (Adaptive Mesh Refinement Framework) AMReX is a required dependency that provides the core AMR data structures, solvers, and parallelization infrastructure for ERF. The primary configuration choice for a developer is whether to use the internal Git submodule provided with the ERF repository or link against an external, pre-built version of the library. The default and strongly recommended approach is to use the internal submodule, as this guarantees version compatibility and provides the most robust and tested configuration for most users. Opting for an external build is an advanced use case, typically reserved for environments where a single, validated AMReX installation must be shared across multiple applications to ensure consistency.
  * CMake Option: -DERF_USE_INTERNAL_AMREX=ON (Default)
  * GNU Make Variable: The AMREX_HOME variable is set to point to the submodule directory.

After configuring the required and optional libraries, developers may occasionally encounter build failures. The next section provides a guide to troubleshooting the most common issues.

4. Troubleshooting Common Library Configuration Issues

4.1. Introduction: A Systematic Approach to Diagnosis

Library configuration issues are a common source of build problems, particularly on complex HPC systems with module-based environments. These errors typically arise from a mismatch between the user's shell environment and the expectations of the build system. This section provides a systematic guide to diagnosing and resolving the most frequent errors by identifying their root causes and providing clear, actionable solutions.

4.2. Diagnosing and Resolving Common Errors

* Detection Failures ("Library Not Found") This error occurs when the build system cannot locate the necessary header files or library files for a dependency like NetCDF or HDF5.
* Solution:
  1. Verify Environment Modules: Use the module list command to confirm that the correct modules are loaded. On Cray systems, this means ensuring cray-netcdf-hdf5parallel is present in the list.
  2. Specify Paths Manually: If the library is installed in a non-standard location and not managed by environment modules, provide an explicit path to its installation prefix (e.g., cmake -DNETCDF_DIR=/path/to/netcdf ..).
  3. Check pkg-config: Ensure the pkg-config utility is in your PATH and that the PKG_CONFIG_PATH environment variable is set correctly to include the directory containing the library's .pc (package-config) file.
* MPI Linker Errors ("Undefined Reference to MPI_*") This error occurs during the final linking stage of the build and indicates that the linker cannot find implementations for core MPI functions.
* Solution:
  1. Diagnosis: This error almost always indicates a mismatch where a parallel (MPI-enabled) ERF build is being linked against a serial version of a dependency, most commonly NetCDF or HDF5. This error is a direct manifestation of a broken parallel I/O dependency chain, as described in Section 2.5.
  2. Resolution: Unload any serial library modules and load the correct parallel-enabled versions. Once the environment is corrected, perform a clean build to ensure the new configuration is applied.
* Stale CMake Cache (Unexpected Failures After Environment Changes) This issue manifests as unexpected build failures that occur after you have changed your environment, such as loading a different module or switching compiler versions.
* Solution:
  1. Diagnosis: CMake caches the locations of libraries and tools to speed up reconfiguration. If the underlying environment changes, this cache can become "stale," pointing to incorrect or non-existent paths.
  2. Resolution: The safest and most reliable solution is to perform a clean build. Use the distclean target, which removes all CMake-generated files and build artifacts: make distclean or cmake --build . --target distclean. Alternatively, you can manually remove the cache files (CMakeCache.txt, CMakeFiles/) from your build directory before re-running the CMake configuration command.

4.3. Advanced Diagnosis Using Build Logs

For complex issues where the cause is not immediately obvious, CMake's verbose logging capabilities are invaluable for diagnostics. These logs reveal the exact search paths, environment variables, and logic used by the find_package commands, making it possible to trace precisely why a particular library is not being found.

To generate a detailed log showing the hierarchy and context of configuration messages, use the following command from your build directory:

# Generate a verbose, context-aware log of the CMake configuration process
cmake --log-level=VERBOSE --log-context ..


This command provides a structured, contextual view of the entire configuration process, allowing you to trace the decision-making logic of the build system and pinpoint the exact source of a configuration failure. The output provides a clear sequence of events:

[ERF.Cray] Detected Cray Programming Environment (CRAYPE_VERSION=2.7.30)
[ERF.Cray] Setting Cray compiler wrappers...
[ERF.Cray]   Set CMAKE_CXX_COMPILER = /opt/cray/pe/craype/default/bin/CC
[ERF.AMReX] Using internal AMReX submodule
[ERF.NetCDF] Found NetCDF: /opt/cray/pe/netcdf/4.9.0.9
