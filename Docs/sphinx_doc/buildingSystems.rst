ERF Build Systems and Options

The Energy Research and Forecasting (ERF) code, built upon the AMReX framework, possesses a complex dependency graph that includes multiple advanced physics packages (RTE-RRTMGP, SHOC, P3). To ensure performance portability across diverse computing architectures, ERF relies on the Kokkos library. This document provides a detailed, code-level explanation of ERF's two supported build systems, GNU Make and CMake, outlining how each system manages these dependencies to produce a functional executable.

The two build systems serve distinct purposes. GNU Make represents a traditional, application-centric system ideal for compiling a single executable for a specific scientific run; it excels at fine-grained debugging via utility targets like print-xxx. In contrast, CMake is the modern, recommended, framework-centric approach whose strength lies in policy-based automation, such as its Cray Detection System. It is designed to create versioned libraries, manage dependencies through automated detection, and support robust cross-platform builds, particularly on High-Performance Computing (HPC) systems.

This guide focuses exclusively on how the build systems are implemented and function, with direct references to the source code. It is intended to be a technical reference for developers and advanced users. Initial environment setup, dependency installation, and specific build examples are covered in other documentation sections.

The following sections begin with a detailed examination of the GNU Make system.

1. Building with GNU Make

1.1. System Overview and Use Case

The GNU Make system is a strategically important build method for ERF, offering a fast and straightforward path to producing a single, case-specific executable. It is designed around an application-centric approach, where developers and scientists compile the code for a particular problem setup (e.g., an atmospheric boundary layer simulation) directly within the Exec/ directory structure.

The primary use case for the GNU Make system is for scientific production runs and development scenarios where in-place compilation is preferred, avoiding the overhead of creating and installing versioned libraries. The entire build is orchestrated by a GNUmakefile located within a specific case directory (such as Exec/ABL/), which in turn leverages the powerful and generic build logic provided by the underlying AMReX framework.

The following sections provide a technical breakdown of how this system orchestrates the build process.

1.2. How it Works: The Orchestration Process

The ERF GNU Make process is architected as a hierarchy of includes that separates user configuration from application and framework build logic. The user's GNUmakefile serves as the top layer, defining high-level build options. This file includes the application-level Make.ERF to gather all necessary ERF source files, which in turn includes the core framework-level Make.defs and Make.rules from AMReX to execute the compilation and linking process.

The process unfolds in the following steps:

1. GNUmakefile Location The process begins with the user invoking make in an application-specific directory, such as Exec/ABL/. This directory contains the primary control file, GNUmakefile.
2. Set AMREX_HOME The first crucial step within the GNUmakefile is defining the AMREX_HOME variable. This variable points to the location of the AMReX submodule, which contains the core build logic. The default path is $(ERF_HOME)/Submodules/AMReX.
3. Define Build Variables The user defines a series of boolean flags and other variables in the GNUmakefile to control which features are compiled into the final executable. These include options like USE_MPI, USE_CUDA, USE_RRTMGP, and others that enable or disable specific physics packages and parallel programming models.
4. Include ERF Sources The GNUmakefile includes the file Exec/Make.ERF. This file systematically adds all necessary ERF source and header directories to the build paths, populating the VPATH_LOCATIONS and INCLUDE_LOCATIONS variables used by the make system.
5. Include AMReX Core Logic Finally, Make.ERF includes the core AMReX makefiles, Make.defs and Make.rules. These files contain the generic logic for discovering dependencies, compiling C++, C, and Fortran source files, and linking the resulting object files into the final executable.

This layered approach cleanly separates user configuration, application-specific source management, and framework-level build logic.

1.3. Important Make Variables

The build is customized by setting variables within the application's GNUmakefile. The following table details the most important variables, covering both core AMReX options and ERF-specific flags that control the inclusion of physics packages and other features.

Variable	Description	Default/Example
AMREX_HOME	Specifies the path to the AMReX source directory.	$(ERF_HOME)/Submodules/AMReX
COMP	Defines the compiler suite to use (e.g., gnu, intel, cray).	gnu
DIM	Sets the dimensionality of the problem. ERF is 3D only.	3
USE_MPI	Enables MPI for distributed-memory parallel execution.	FALSE
USE_OMP	Enables OpenMP for shared-memory parallelism.	FALSE
USE_CUDA	Enables NVIDIA GPU support via CUDA. Defines the ERF_USE_CUDA preprocessor macro.	FALSE
USE_HIP	Enables AMD GPU support via HIP. Defines the ERF_USE_HIP preprocessor macro.	FALSE
USE_SYCL	Enables Intel GPU support via SYCL. Defines the ERF_USE_SYCL preprocessor macro.	FALSE
USE_NETCDF	Enables I/O support using the NetCDF library. Defines the ERF_USE_NETCDF preprocessor macro.	FALSE
USE_NOAHMP	Enables the Noah-MP land surface model. This option requires USE_NETCDF=TRUE.	FALSE
USE_RRTMGP	Enables the RRTMGP radiation model. This option automatically sets USE_KOKKOS=TRUE and USE_NETCDF=TRUE, and defines ERF_USE_RRTMGP.	FALSE
USE_SHOC	Enables the SHOC turbulence model. This option automatically sets USE_KOKKOS=TRUE and defines ERF_USE_SHOC.	FALSE
USE_P3	Enables the P3 microphysics model. This option automatically sets USE_KOKKOS=TRUE and defines ERF_USE_P3.	FALSE
USE_KOKKOS	Enables the Kokkos performance portability library, which is a dependency for RRTMGP, SHOC, and P3. Defines ERF_USE_KOKKOS and links against Kokkos libraries.	FALSE
USE_MORR_FORT	Enables the Fortran-based Morrison microphysics scheme. Defines the ERF_USE_MORR_FORT preprocessor macro.	FALSE
USE_FFT	Enables Fast Fourier Transform capabilities. Defines the ERF_USE_FFT preprocessor macro.	FALSE

The next section points to an example file where these variables are put into practice.

1.4. Example GNUmakefile

A typical GNUmakefile for an ERF application can be found at ERF/Exec/ABL/GNUmakefile.

1.5. Common Commands

The AMReX make system provides several standard utility targets for managing the build environment and its artifacts. For a complete list of commands and advanced features, refer to the official AMReX documentation. The most common commands are described below.

* make Compiles the source code and creates the executable in the current directory.
* make clean Removes all build artifacts for all configurations. This includes executables, object files (.o), and temporary build directories.
* make cleanconfig Removes only the artifacts for the current configuration, which is defined by variables like COMP, USE_MPI, etc. This allows you to clean a specific build type (e.g., the MPI debug build) while leaving other compiled configurations (e.g., a serial release build) intact.
* make print-xxx A powerful debugging tool that prints the value of a specified make variable (xxx) to the console. For example, make print-CXXFLAGS will display the C++ compiler flags being used for the current configuration.

1.6. Customization with Make.local

For user- or site-specific customizations that should not be committed to the main repository, the AMReX GNU Make system provides the Make.local file. This file, located at amrex/Tools/GNUMake/Make.local, is the primary mechanism for overriding default build variables such as compiler paths or flags.

Any settings defined in Make.local apply globally to all projects that use that specific instance of the AMReX submodule. This makes it an ideal location for site-wide policies or personal developer preferences.

A common use case is to specify a non-default compiler version. For example, to force the build system to use g++-8, you would add the following lines to your Make.local file:

CXX = g++-8
CC  = gcc-8
FC  = gfortran-8


While the GNU Make system is robust for application-focused builds, CMake is the recommended system for its flexibility and robust dependency management.

2. Building with CMake

2.1. System Overview and Rationale

CMake is the recommended build system for ERF. It is a powerful, cross-platform tool designed with a framework-centric philosophy. Instead of producing a single executable for one use case, CMake is architected to produce versioned, exportable libraries. This design is primarily enabled by its ability to generate and install a find_package-compatible configuration file (ERFConfig.cmake), which allows any other CMake-based project to consume ERF as a dependency with a simple find_package(ERF) command. This makes it ideal for complex "superbuilds," integration into larger scientific workflows, and managing dependencies on diverse HPC platforms.

Key advantages of using CMake for ERF include:

* It promotes clean, out-of-source builds, which keep all compiled artifacts separate from the source code. This practice prevents clutter in the source tree and allows multiple build configurations to coexist safely.
* It provides robust dependency detection through its find_package mechanism. This is used to locate and configure essential libraries like MPI (CMakeLists.txt:328-406) and NetCDF (CMakeLists.txt:408-426), ensuring that the correct headers and libraries are used during compilation.
* It automatically detects and configures for Cray HPC systems, applying necessary workarounds and generating a cray_detected_config.cmake file that caches the auto-detected system settings, ensuring reproducibility and providing a template for manual configuration.

The following section outlines the standard procedures for building ERF using CMake.

2.2. Standard Workflow

ERF offers several methods for invoking a CMake build, catering to different user needs from simple, pre-configured scripts to standard manual command-line invocation for advanced users.

* Simple Script The most straightforward method is to run Build/cmake.sh. This script provides a basic, tested configuration designed to work on most standard Linux systems. It is an excellent starting point for new users.
* Advanced Script For more complex build scenarios, ERF provides specialized scripts like Build/cmake_with_kokkos_many_cuda.sh. These scripts are tailored for specific use cases, such as enabling the Kokkos-based physics packages on NVIDIA GPU architectures, and they pre-configure the necessary options.
* Manual Invocation The standard practice is to create a build directory (e.g., mkdir build), change into it (cd build), and then invoke CMake pointing to the source directory (cmake ..). An equivalent single command from the project root is cmake -B build -S .. This provides maximum control over the configuration.

The out-of-source build practice is highly recommended because it keeps the source directory pristine, allowing multiple, independent build configurations to coexist. To help manage this, ERF's CMake configuration provides uninstall and distclean targets. The standard manual procedure, executed from the project's root directory, can be summarized as follows:

1. mkdir build Create a new directory to contain all build artifacts, such as object files and the final executable.
2. cd build Change into the newly created build directory. All subsequent commands are run from here.
3. cmake [options] .. Configure the project. The [options] are flags used to customize the build, specified with the syntax -D<VARIABLE>=<VALUE>. The .. argument tells CMake that the source code is located one directory level up.
4. make install Compile the source code and install the resulting executable and any libraries into the directory specified by the CMAKE_INSTALL_PREFIX variable.

The next section details the specific options available for customizing the ERF build.

2.3. ERF CMake Options

The ERF build is customized by passing options to the cmake command using the -D<VARIABLE>=<VALUE> syntax. The table below lists the primary CMake options available for configuring ERF.

Variable Name	Description	Default	Possible Values
ERF_ENABLE_MPI	Enables MPI support for parallel execution.	OFF	ON, OFF
ERF_ENABLE_OPENMP	Enables OpenMP for shared-memory parallelism.	OFF	ON, OFF
ERF_ENABLE_CUDA	Enables NVIDIA GPU support via CUDA. Requires CUDA Toolkit >= 11.0.	OFF	ON, OFF
ERF_ENABLE_HIP	Enables AMD GPU support via HIP.	OFF	ON, OFF
ERF_ENABLE_SYCL	Enables Intel GPU support via SYCL.	OFF	ON, OFF
ERF_ENABLE_NETCDF	Enables NetCDF for I/O operations.	OFF	ON, OFF
ERF_ENABLE_NOAHMP	Enables the Noah-MP land surface model. Requires ERF_ENABLE_NETCDF=ON.	OFF	ON, OFF
ERF_ENABLE_RRTMGP	Enables the RRTMGP radiation model. This option requires ERF_ENABLE_NETCDF=ON. It also enables the EKAT submodule, which provides the Kokkos library, and therefore requires ERF_ENABLE_MPI=ON.	OFF	ON, OFF
ERF_ENABLE_SHOC	Enables the SHOC turbulence model. This option enables the EKAT submodule, which provides the Kokkos library, and therefore requires ERF_ENABLE_MPI=ON.	OFF	ON, OFF
ERF_ENABLE_P3	Enables the P3 microphysics model. This option enables the EKAT submodule, which provides the Kokkos library, and therefore requires ERF_ENABLE_MPI=ON.	OFF	ON, OFF
ERF_ENABLE_PARTICLES	Enables support for Lagrangian particles.	OFF	ON, OFF
ERF_ENABLE_DOCUMENTATION	Enables the build target for generating Sphinx documentation.	OFF	ON, OFF
CMAKE_BUILD_TYPE	Sets the build configuration. Defaults to Release if not specified, following standard CMake behavior.	Release	Release, Debug, RelWithDebInfo
CMAKE_INSTALL_PREFIX	The directory where make install will place the compiled artifacts. ERF's build scripts often override this to a local ./install directory for convenience.	System-dependent (e.g., /usr/local)	/path/to/install

For complex configurations, managing these options on the command line can be cumbersome. The next section introduces a more robust file-based approach.

2.4. Using Configuration Files

For managing complex build configurations, CMake provides the -C <file> option. This allows users to specify build variables in a CMake-syntax file rather than on the command line. Using a configuration file is more reproducible, less error-prone, and easier to share than a long command-line invocation. The command syntax is as follows: cmake -C path/to/config.cmake ..

This approach is distinct from the machine-specific profile scripts ERF provides in the Build/machines/ directory. These profiles are shell scripts intended to be sourced into your terminal session (e.g., source machines/perlmutter_erf.profile) before running CMake. Their primary function is to set up the build environment by loading the correct software modules (module load ...) and exporting environment variables (NETCDF_DIR, CRAY_MPICH_DIR, etc.) that the build system depends on. This is the canonical practice on HPC systems to ensure builds are repeatable and portable across users. Example profiles include:

* aurora_erf.profile
* frontier_erf.profile
* perlmutter_erf.profile
* polaris_erf.profile

2.5. Advanced CMake Features

2.5.1. Hierarchical Logging

For CMake versions 3.25 and newer, ERF's build configuration supports hierarchical logging. This advanced feature greatly aids in diagnosing build issues and understanding the decisions made by the build system during the configuration phase. It provides a clear, structured view of which component is responsible for each configuration message.

Two primary flags control this feature:

* --log-level Controls the verbosity of the CMake output. The most common levels are VERBOSE, which shows details of dependency detection, and DEBUG, which shows all diagnostic information.
* --log-context Enables message grouping by component, annotating each line of output with its origin (e.g., [ERF.Cray], [ERF.AMReX]). This makes it easy to trace the configuration logic for a specific part of the build.

An example of the structured output generated by using --log-context is shown below:

[ERF.Cray] Detected Cray Programming Environment (CRAYPE_VERSION=2.7.30)
[ERF.Cray] Setting Cray compiler wrappers...
[ERF.Cray]   Set CMAKE_CXX_COMPILER = /opt/cray/pe/craype/default/bin/CC
[ERF.AMReX] Using internal AMReX submodule
[ERF.NetCDF] Found NetCDF: /opt/cray/pe/netcdf/4.9.0.9


This output clearly shows the sequence of events: the Cray detection module identifies the environment and sets the appropriate compiler, after which the AMReX and NetCDF configuration steps report their status.

3. Building the Documentation

3.1. Overview and Procedure

ERF uses the Sphinx documentation generator to produce comprehensive project documentation from source files. The process for building the documentation is integrated directly into the CMake build system, providing a streamlined workflow for developers.

The procedure to build the documentation is straightforward and consists of two main steps:

1. Enable Documentation The documentation build is an optional target and must first be enabled during the CMake configuration step. This is done by setting the CMake option ERF_ENABLE_DOCUMENTATION=ON.
2. Build the Target Once the project is configured with documentation enabled, the documentation can be generated by invoking the docs target. For example, with the Makefiles generator, you would run make docs.

This process requires that Sphinx is installed on the system. The CMake configuration will automatically search for the sphinx-build executable and will report an error if it cannot be found.
