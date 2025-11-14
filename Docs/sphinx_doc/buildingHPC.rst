ERF on HPC Systems: A Guide to the Cray Detection Framework

Introduction and Environment Setup

Building sophisticated scientific software on High-Performance Computing (HPC) systems is an inherently complex task. The landscape is characterized by a diversity of supercomputing architectures, proprietary compiler toolchains, and specialized environment module systems, all of which create significant challenges for software compilation. The ERF build system is engineered to automate and manage this complexity, providing a robust and reproducible path from source code to executable.

A mandatory first step before any compilation attempt on an HPC platform is the correct preparation of the shell environment. Loading the appropriate modules and setting environment variables is critical, as this process ensures that the build system can locate the correct compiler wrappers, essential parallel libraries like MPI, hardware-specific toolchains such as the CUDA or ROCm toolkits, and I/O libraries like NetCDF and HDF5. Without a properly configured environment, builds are almost certain to fail with cryptic errors related to missing dependencies.

To standardize this crucial setup process, ERF provides a set of machine profile files. These are pre-configured scripts that encapsulate the system-specific knowledge required to establish a correct build environment on a given supercomputer. By using these profiles, developers and users can reliably configure their shell sessions, ensuring a consistent and reproducible foundation for compilation.

This document explains the architecture of ERF's HPC build strategy, with a particular focus on the automated Cray detection system integrated within its CMake build framework. It is intended to provide a high-level understanding of the automation logic and design philosophy rather than serve as a step-by-step tutorial for building the code.

We will begin by examining the role and usage of the machine profile files that form the bedrock of the HPC build process.

1. Machine Profile Files

Machine profile files are the foundational layer of the ERF build process on HPC systems. They are strategically designed to abstract away the intricate details of diverse supercomputer ecosystems by managing environment modules and variables in a standardized, reproducible manner. A correctly configured environment is the single most critical prerequisite for any successful HPC build, and these profiles are the primary tool for achieving that configuration.

1.1. Purpose and Usage

Machine profile files are shell scripts located in the Build/machines/ directory. Their purpose is to prepare the shell environment for compilation by loading required software modules, exporting necessary environment variables, and defining paths to dependencies.

The primary functions of these profiles are to:

* Execute a series of module load commands. The core of the environment on Cray systems is the PrgEnv-* module (e.g., PrgEnv-gnu or PrgEnv-cray), which dictates the underlying compiler suite that the Cray wrappers (cc, CC) will use. This is then supplemented by modules for hardware acceleration (e.g., craype-accel-nvidia80), parallel libraries (e.g., cray-netcdf-hdf5parallel), and development tools (e.g., cmake).
* Run export commands to set environment variables that guide the build system, such as hints for locating libraries or specifying target hardware.

To use a profile, a user must source the appropriate script from their shell. This executes the commands within the current shell session, modifying its environment. Once the profile has been sourced, the shell is ready for the subsequent cmake configuration and build commands. For example: source machines/perlmutter_erf.profile.

1.2. Available Profiles

ERF provides pre-configured machine profiles for several major DOE HPC systems:

* perlmutter_erf.profile (NERSC Perlmutter)
* frontier_erf.profile (OLCF Frontier)
* aurora_erf.profile (ALCF Aurora)
* polaris_erf.profile (ALCF Polaris)
* kestrel_erf.profile (NERSC Kestrel)

1.3. Customizing for a New System

To build ERF on an unsupported HPC system, you can create a custom profile. The recommended approach is to copy an existing profile for a machine with a similar architecture (e.g., another system using NVIDIA GPUs or one based on AMD hardware) and modify the module load and export commands to match the target system's specific software environment and requirements.

While manual profile files are essential for environment preparation, they address only the state of the shell. They cannot manage the in-build configuration complexities that arise once CMake begins its process. Sourcing a profile sets environment variables, but CMake operates in its own configuration space and requires explicit information about compiler flags, library paths, and dependency relationships that module load commands do not directly provide to its internal logic. To bridge this gap, the integrated Cray detection system provides the next layer of automation.

2. The Cray Detection System (CMake)

The Cray Detection System is a sophisticated, two-phase automation mechanism within ERF's CMake build system. It is designed to programmatically identify Cray Programming Environments on HPC platforms like Perlmutter, Frontier, and Polaris, and then apply a series of necessary configurations and workarounds. The primary goal of this system is to shield users from the immense platform-specific complexity of these machines, delivering a consistent, reliable, and reproducible build experience with minimal manual intervention.

The system automates several critical configuration steps:

* It sets the Cray compiler wrappers (cc, CC, ftn) as the default compilers for C, C++, and Fortran.
* It inspects the $CRAY_ACCEL_TARGET environment variable to determine the specific GPU architecture (e.g., NVIDIA A100, AMD MI250X).
* It automatically sets the correct GPU architecture flags for key dependencies, including AMReX (AMReX_CUDA_ARCH, AMReX_AMD_ARCH) and Kokkos-based physics packages (Kokkos_ARCH_*).
* It configures support for GPU-aware MPI and intelligently finds the correct parallel versions of libraries like NetCDF and HDF5, which standard methods often fail to locate.

2.1. Overview of the Two-Phase Process

This two-phase design is a classic workaround for a fundamental constraint in CMake's architecture: the toolchain must be validated before the project() command is invoked, but complex, project-specific logic can only run after project(). Our system splits detection accordingly, satisfying CMake's procedural needs while enabling powerful, post-project automation.

* Phase 1: Pre-Project Detection This phase is handled by CrayCompilerDetection.cmake and executes before the project() call in the main CMakeLists.txt file. It performs a lightweight check of environment variables like $CRAYPE_VERSION and $CRAY_MPICH_DIR to determine if it is running in a Cray environment. If a Cray system is detected, it immediately sets the CMAKE_C_COMPILER and CMAKE_CXX_COMPILER cache variables to cc and CC, respectively. This ensures that CMake's initial compiler identification and test compilations succeed with the correct Cray compiler wrappers.
* Phase 2: Post-Project Configuration After the project() command has successfully run and the compilers are known, the second phase is initiated by CrayDetection.cmake. This module performs a much deeper inspection of the environment. It detects the specific GPU architecture from $CRAY_ACCEL_TARGET, configures dependencies for the unique Cray ecosystem, and applies a series of automated workarounds for known build issues related to GPU-aware MPI, parallel NetCDF, and interactions between CUDA and Kokkos-based packages.

2.2. Controlling Auto-Detection

The Cray auto-detection system is enabled by default whenever a Cray environment is identified. For advanced users, debugging, or situations requiring a fully manual configuration, this automation can be disabled by setting a single CMake option: -DERF_DISABLE_CRAY_AUTO_FIXES=ON. Setting this flag will prevent the system from applying any of the automated workarounds, giving the user complete control over the configuration.

2.3. Key Detection Steps and Automated Fixes

The Cray detection system is engineered to automatically resolve several common and complex configuration challenges encountered on Cray supercomputers. These automated workarounds are all examples of the Post-Project Configuration (Phase 2) described previously. They can only be applied after the compilers have been successfully identified and validated during Phase 1, reinforcing the architectural necessity of the two-phase design. The following table summarizes the most important automated workarounds.

Problem	Detection Clue	Automated Solution
GPU Architecture<br>Manually specifying the correct GPU architecture flags for multiple libraries (AMReX, Kokkos) is error-prone.	The $CRAY_ACCEL_TARGET environment variable is set (e.g., to nvidia80 or amd_gfx90a).	The system automatically maps $CRAY_ACCEL_TARGET to the correct architecture flags for each dependency. For example, nvidia80 sets AMReX_CUDA_ARCH to 8.0 and enables Kokkos_ARCH_AMPERE80.
CUDA+EKAT Builds<br>When using CUDA with Kokkos-based physics packages (RRTMGP, SHOC, P3), the CUDA compiler (nvcc) is not automatically aware of the MPI include paths set by the Cray environment. This leads to compilation failures where MPI headers cannot be found during the CUDA compilation phase.	ERF_ENABLE_CUDA=ON and ERF_ENABLE_RRTMGP, ERF_ENABLE_SHOC, or ERF_ENABLE_P3 is enabled.	The system executes CC --cray-print-opts=cflags to retrieve the missing compiler flags and prepends them to the CMAKE_CUDA_FLAGS variable, ensuring the CUDA compiler has access to the required MPI headers.
GPU-Aware MPI<br>Correctly linking the MPI GPU Transport Layer (GTL) libraries is complex and system-dependent.	The $MPICH_GPU_SUPPORT_ENABLED=1 environment variable is set.	The system identifies the base MPI library and automatically links the appropriate GTL library, such as -lmpi_gtl_cuda for NVIDIA GPUs or -lmpi_gtl_hsa for AMD GPUs, by adding them to CMAKE_CXX_STANDARD_LIBRARIES and CMAKE_CUDA_STANDARD_LIBRARIES.
Parallel NetCDF/HDF5<br>Standard CMake find_package calls often fail to locate the parallel versions of NetCDF and HDF5 provided through Cray's environment modules.	ERF_ENABLE_NETCDF=ON.	The system retrieves the Cray-specific search path by running CC --cray-print-opts=PKG_CONFIG_PATH and prepends this path to the $PKG_CONFIG_PATH environment variable. This guides the pkg-config tool to find the correct parallel library installations.

2.4. Manual Configuration via -C Files

For advanced users requiring ultimate control, the auto-detection system provides an "escape hatch" that leverages transparency. This workflow allows you to inspect, modify, and re-apply the system's configuration manually.

1. After a successful configuration run, the auto-detection logic saves all of its detected settings to a file named cray_detected_config.cmake in the build directory.
2. The build system provides a convenient utility target, make show-cray-config, which prints the contents of this file to the console without requiring you to locate it manually. This reveals the exact output of the auto-detection logic.
3. This generated cray_detected_config.cmake file serves as the perfect template for creating a custom manual configuration. You can copy it, modify the settings as needed, and then pass it to CMake using the -C option (e.g., cmake -C path/to/my_config.cmake ..) to apply a fully manual configuration.

While the Cray Detection System provides powerful automation within CMake, ERF maintains flexibility by offering users a choice of build tools. The following section details the provided build scripts and compares the strategic advantages of the modern, automated CMake system against the traditional, explicit control offered by GNU Make.

3. Build Scripts and System Comparisons

To further simplify the user experience, ERF provides a collection of tested build scripts that encapsulate common configurations for various systems and architectures. The choice between using these scripts with the modern, automated CMake system and the traditional, explicit GNU Make system depends on the user's specific goals, familiarity with the tools, and the level of control required.

3.1. Script Reference Table

The following table summarizes the build scripts provided with ERF and indicates which combinations have been tested and verified on each major platform. Verified builds are marked with the git commit hash at which they were last successfully tested, providing a clear record of known-good configurations.

Build Script	Perlmutter	Frontier	Aurora	Polaris	Kestrel	RegtestCPU	RegtestGPU
cmake.sh	f8665c28 (ABL)	f8665c28 (ABL)	f8665c28	f8665c28 (ABL)	Untested	f8665c28 (ABL)	f8665c28 (ABL)
cmake_with_kokkos_many.sh	Untested	Untested	Untested	Untested	Untested	f8665c28 (ABL)	Untested
cmake_with_kokkos_many_cuda.sh	f8665c28 (ABL)	—	—	Untested	Untested	—	f8665c28 (ABL)
cmake_with_kokkos_many_noradiation_hip.sh	—	f8665c28 (ABL bad fcompare result)	—	—	—	—	—
cmake_with_kokkos_many_hip.sh	—	f8665c28 + rrtmpg_workarounds_kokkos	—	—	—	—	—
cmake_with_kokkos_many_sycl.sh	—	—	f8665c28 + rrtmpg_workarounds_kokkos	—	—	—	—
Perlmutter/build_erf_with_shoc_cuda_Perlmutter.sh	Untested	Untested	Untested	Untested	Untested	—	Untested

3.2. Script Descriptions

* cmake.sh: A basic CMake build script for CPU-only configurations that works on most systems.
* cmake_with_kokkos_many.sh: A Kokkos-enabled build intended for CPU architectures.
* cmake_with_kokkos_many_cuda.sh: A Kokkos-enabled build optimized for NVIDIA GPUs using the CUDA backend.
* cmake_with_kokkos_many_noradiation_hip.sh: A Kokkos-enabled build for AMD GPUs (HIP backend) with the radiation model disabled.
* cmake_with_kokkos_many_hip.sh: A full Kokkos-enabled build for AMD GPUs using the HIP backend.
* cmake_with_kokkos_many_sycl.sh: A Kokkos-enabled build for Intel GPUs using the SYCL backend.
* Perlmutter/build_erf_with_shoc_cuda_Perlmutter.sh: A specialized script for building the SHOC turbulence model with CUDA support on Perlmutter.

3.3. GNU Make vs. CMake on HPC

Both GNU Make and CMake are fully supported build systems, and the choice between them is a strategic decision based on the desired level of abstraction and control. CMake provides a high-level, policy-based abstraction that manages complexity automatically, which is ideal for integration and reproducibility. In contrast, GNU Make offers a low-level, imperative control model that gives experts direct, transparent command over every flag and dependency, making it ideal for fine-grained debugging and performance tuning.

Build System	Strengths & Best Use Cases
CMake	Strengths: Superbuilds and handling the complex dependency graph (e.g., ERF -> EKAT -> Kokkos) are much easier. Automates configurations via the Cray Detection System. Generates and installs a find_package-compatible configuration (ERFConfig.cmake), allowing other CMake projects to easily consume ERF as a dependency.<br>Best For: Most users on supported HPC systems (Perlmutter, Frontier, Polaris, Aurora). Developers building external applications that depend on ERF.
GNU Make	Strengths: Leverages existing AMReX infrastructure and has extra utility targets (e.g., print-xxx) for debugging the build. Offers transparent and direct control over every build variable by manually setting them in a GNUmakefile. This is a strength for debugging but a weakness in managing complexity.<br>Best For: Advanced users, performance tuning experts, and those building on HPC systems not yet fully supported by the CMake auto-detection logic.

3.4. GNU Make Configuration Example

Configuration with GNU Make is handled by setting variables directly in a GNUmakefile or on the command line, such as make COMP=cray USE_MPI=TRUE USE_CUDA=TRUE. This approach provides direct and explicit control over every aspect of the build, leveraging the powerful underlying AMReX make system.

While HPC systems are the primary target for large-scale runs, ERF is also designed to be built and run on local workstations for development and testing.

4. Workstation Builds

Building ERF on a local workstation is essential for development, debugging, and running smaller test cases. The process is typically simpler than on an HPC system, as it does not require managing environment modules or navigating complex, vendor-specific toolchains.

The RegtestCPU and RegtestGPU columns in the script reference table indicate configurations that are regularly tested and verified to work in common workstation environments.

macOS-Specific Considerations

Building on macOS is supported, though it comes with specific considerations. The following practices are recommended to ensure a smooth build process:

* GNU Make is often recommended for its simplicity and the direct control it offers, which can be beneficial in the more heterogeneous macOS environment.
* A Make.local file can be used to explicitly specify a consistent compiler suite (e.g., GCC installed via Homebrew) for both C++ and Fortran.
* This approach helps avoid potential conflicts between the default Clang/Xcode C++ compiler and a separately installed Fortran compiler (e.g., gfortran), which can lead to linking issues that arise from mixing compiler toolchains.
