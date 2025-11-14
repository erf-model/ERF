
1. Resolved Issues (Historical)

Documenting issues that have already been resolved is a crucial practice for mature software projects. This section demonstrates the build system's increasing robustness by chronicling historical complexities that are now handled automatically. More importantly, it serves as a valuable diagnostic tool for advanced users. By understanding the automated fixes for known issues on platforms like Cray, users can more effectively diagnose new problems if they suspect the automation is not functioning as expected.

1.1. Cray Compiler and GPU-Aware MPI Linker Errors (Fix 4)

On Cray systems with GPU-aware MPI enabled (MPICH_GPU_SUPPORT_ENABLED=1), historical builds often failed during the final linking stage. The root cause was the Cray linker's default use of the --as-needed flag, which aggressively discards libraries it deems unused during an initial pass. This process would incorrectly remove the essential MPI GPU Transport Layer (GTL) libraries required for communication between the CPU and GPU.

The build system now automatically detects this configuration. It identifies the correct MPI base library (e.g., mpi_gnu_123) and the required GTL library for the target hardware (e.g., mpi_gtl_cuda for NVIDIA GPUs or mpi_gtl_hsa for AMD GPUs). These libraries are programmatically added to the CMAKE_CXX_STANDARD_LIBRARIES and CMAKE_CUDA_STANDARD_LIBRARIES variables, ensuring they are present during the final link and preventing these errors.

* Phase 1 (Pre-Project): The logic in CrayCompilerDetection.cmake runs before CMake's project() command. It performs initial detection of GPU-aware MPI and appends the required GTL libraries (e.g., -lmpi_gnu_123 -lmpi_gtl_cuda) to CMAKE_CXX_STANDARD_LIBRARIES and CMAKE_CUDA_STANDARD_LIBRARIES. This is a critical architectural step that ensures CMake’s own compiler test compilations succeed using the Cray wrappers.
* Phase 2 (Post-Project): After the project is configured, the more robust logic in CrayDetection.cmake (Fix 4) confirms and applies this fix, solidifying it as a core part of the Cray automation.

This two-phase strategy is a fundamental design choice that correctly addresses both CMake's procedural requirements and the complexities of the Cray environment.

1.2. NetCDF Parallel Library Detection on Cray Systems (Fix 5-6)

Standard CMake find_package calls often fail on Cray systems because critical libraries like parallel NetCDF and HDF5 are not in default system paths. Instead, their locations are managed by environment modules. Previously, this required users to manually provide hints to the build system via the NETCDF_DIR variable.

To automate this, the build system now circumvents the limitations of find_package by querying the Cray compiler wrapper directly. It executes CC --cray-print-opts=PKG_CONFIG_PATH to discover the correct search path provided by the loaded Cray environment modules. This path is then prepended to the PKG_CONFIG_PATH environment variable, which reliably guides the pkg-config utility to find the correct parallel versions of NetCDF and HDF5.

* CrayDetection.cmake: The implementation for "Fix 5-6" is verified in this file, which explicitly executes ${CMAKE_CXX_COMPILER} --cray-print-opts=PKG_CONFIG_PATH and prepends the output to the PKG_CONFIG_PATH environment variable.

1.3. Manual GPU Architecture Specification for Kokkos/AMReX

Previously, users were required to manually specify the target GPU architecture for all relevant dependencies, such as AMReX and Kokkos-based physics packages. This process was error-prone, as a mismatch between the specified architecture and the actual hardware could lead to compilation failures or non-functional executables.

This is now fully automated. The build system inspects the $CRAY_ACCEL_TARGET environment variable (e.g., nvidia80 for NVIDIA A100 or amd_gfx90a for AMD MI250X), which is set when a user loads the appropriate craype-accel-* module. This variable is then used to programmatically set the correct GPU architecture flags for all dependencies, including AMReX_CUDA_ARCH for AMReX and the appropriate Kokkos_ARCH_* variable (e.g., setting the CMake option Kokkos_ARCH_AMPERE80 to ON). This automation eliminates a significant source of user error.


* CrayDetection.cmake: The script contains distinct logic blocks that inspect $ENV{CRAY_ACCEL_TARGET} and correctly map its value (e.g., nvidia80) to the required dependency-specific flags, such as set(AMReX_CUDA_ARCH "8.0" ...) and set(Kokkos_ARCH_AMPERE80 ON ...).


--------------------------------------------------------------------------------


2. Troubleshooting Common Issues

While the build system automates many platform-specific complexities, a successful compilation still depends on a correctly configured user environment. This section addresses the most frequent user-side configuration errors that can lead to build failures and provides clear, actionable solutions to resolve them.

2.1. Missing craype-accel-* Module on GPU Builds

For any CUDA or HIP build on a Cray system, the Cray compiler wrappers require a craype-accel-* module to be loaded. This action sets the essential $CRAY_ACCEL_TARGET environment variable, which informs the compiler about the specific GPU architecture to target.

The ERF build system explicitly checks for this variable's presence before the project() command is even called. If the variable is not defined during a GPU-enabled build, CMake will issue a fatal error with clear instructions to load the correct module.

Recommended Solution Load the appropriate module for the system's hardware (e.g., craype-accel-nvidia80 for NVIDIA A100 GPUs). The most reliable way to do this is by sourcing the correct environment script from the Build/machines/ directory before running CMake.

✓ GOOD

Justification: This crucial check is implemented as described.

* CrayCompilerDetection.cmake: This script, which runs before project(), contains logic that checks if CRAYPE_VERSION is defined (indicating a Cray system) and CRAY_ACCEL_TARGET is not defined. If this condition is met for a CUDA or HIP build, it triggers a FATAL_ERROR with a detailed, helpful message. The notes_perlmutter.txt and notes_polaris.txt files confirm craype-accel-nvidia80 is loaded on those systems.

This is an excellent example of a preventative check that stops a build early with a clear, actionable error message.

2.2. NetCDF or HDF5 Libraries Not Found

Enabling NetCDF I/O with -DERF_ENABLE_NETCDF=ON causes the build system to also set -DERF_ENABLE_HDF5=ON by default. The SetAmrexOptions.cmake script then reads the value of ERF_ENABLE_HDF5 and sets AMReX_HDF5=ON, signaling to the AMReX sub-build that HDF5 support is required. The build may fail if the necessary libraries cannot be located.

Recommended Solution On Cray systems, load the cray-netcdf-hdf5parallel module. The Cray detection logic in the ERF build system is specifically designed to locate and configure the parallel libraries provided by this module.

For non-Cray or custom builds, you must point the build system to your NetCDF installation. This is best done by setting the NETCDF_DIR environment variable or CMake variable to the root of your NetCDF installation directory.

✓ GOOD

Justification: The chain of variable settings is accurate, and the recommended solution is correct.

* CMakeLists.txt: The line option(ERF_ENABLE_HDF5 "Enable HDF5 IO" ${ERF_ENABLE_NETCDF}) confirms that enabling NetCDF also enables HDF5 by default.
* CMake/SetAmrexOptions.cmake: The line set(AMReX_HDF5 "${ERF_ENABLE_HDF5}" ...) verifies that this setting is propagated to the underlying AMReX build.
* Machine Notes (notes_*.txt): The module list output for Perlmutter, Frontier, and Polaris all show the cray-netcdf-hdf5parallel module loaded, aligning with the recommended solution.


--------------------------------------------------------------------------------


3. General Debugging Strategies

This section provides guidance for empowering users to diagnose novel or complex issues not covered by the specific troubleshooting cases above. The strategies outlined here represent fundamental skills for debugging software builds on advanced computing platforms and can help illuminate the root cause of unexpected failures.

3.1. CMake Debugging

ERF's CMake configuration includes powerful logging features to trace the build process and understand how decisions are made.

* Verbose Logging: Use the --log-level=VERBOSE flag with the cmake command (e.g., cmake --log-level=VERBOSE ..). This will print detailed information, including the library search paths that CMake is checking and the specific configuration decisions it makes.
* Contextual Logging: Use the --log-context flag to display a component hierarchy for each message (e.g., [ERF.Cray]). This is extremely useful for pinpointing the origin of specific messages, warnings, or errors.
* Cache Inspection: After a configuration attempt, the CMakeCache.txt file in the build directory contains the final values of all configuration variables. This is the definitive source for checking what paths, flags, and options were ultimately used.

📝 BETTER IN DOCS

Justification: While the information is correct, the official building.rst documentation provides a clearer and more practical explanation.

* building.rst (in building.md): The "Logging Levels" and "Logging Context" sections provide explicit, copy-pasteable command examples for using --log-level=VERBOSE and --log-context. It also includes a sample of the structured output, which is highly effective for showing the user what to expect and how to interpret the results. The troubleshooting guide's text is purely descriptive, whereas the official documentation is instructional.

3.2. GNU Make Debugging

Both GNU Make and CMake are fully supported build systems, but the choice between them is a strategic decision based on the desired level of abstraction and control. While CMake provides high-level automation ideal for reproducibility, GNU Make offers a low-level, imperative control model. It is the preferred system for experts who require direct, transparent command over every compiler flag and library path, making it ideal for fine-grained debugging and performance tuning of complex, system-specific issues. For these users, AMReX provides a simple but powerful utility for inspecting the value of any makefile variable.

* To check the value of a variable, use the command make print-<variable_name>. For example, running make print-CXXFLAGS will display the exact C++ compiler flags being used in the build. This is the primary method for debugging paths, flags, and build options in the GNU Make environment.

✓ GOOD

Justification: This is an accurate and useful tip.

* build_systems_skeleton_code_filled.md: The "Common Commands" subsection confirms the existence and purpose of the make print-xxx utility target.
* BuildingAMReX.rst.txt: The official AMReX documentation also highlights make print-xxx as a primary debugging tool, stating it is "very useful for debugging and tweaking the make system."

3.3. Compiler and Linker Debugging

When a build fails during compilation or linking, it is often helpful to see the exact commands being executed.

* Verbose Build: Pass a verbose flag to the build command (e.g., make VERBOSE=1). This will print the full compiler and linker commands, including all flags and library paths, making it easier to spot issues like missing include paths or incorrect library links.
* Shared Library Dependencies: On Linux systems, the ldd command can be run on a compiled executable to check for missing shared library dependencies. This is a quick way to diagnose runtime errors caused by an incomplete environment.

✓ GOOD

Justification: This is standard, correct, and valuable advice for debugging build systems. make VERBOSE=1 is a fundamental tool for inspecting the final compilation commands generated by CMake, and ldd is an essential utility for diagnosing runtime linking failures on Linux systems.


--------------------------------------------------------------------------------


4. Verifying a Successful Build

A successful compilation signals that the code is syntactically correct, but it does not guarantee a functional executable. A final verification or "smoke test" is essential to confirm that the application initializes correctly, runs without crashing, and produces valid results. This step validates that all necessary libraries were linked correctly and are compatible with the user's runtime environment.

4.1. Procedure

1. Navigate to an Example Case. After a successful build, locate the executable. A make install command places it in the install/bin/ directory relative to your build directory. A simple make places it in an Exec/ subdirectory within the build directory (e.g., build/Exec/ABL/erf_abl). To run a test, navigate to the example case's subdirectory within the build directory (e.g., cd build/Exec/ABL/).
2. Execute a Short Run. Run one of the example cases for a small number of timesteps to confirm that the simulation can initialize and advance. The exact MPI launcher command (srun, mpiexec, mpirun) varies by system. The following is a representative command, assuming you are in the correct build subdirectory and the executable is in the current path:
3. Run the Test Suite (Optional). For more comprehensive verification, you can run the built-in regression tests. This requires configuring the build with the -DERF_ENABLE_TESTS=ON flag. The tests can then be executed from the build directory using the command ctest -L regression -VV. This will run the regression test suite with verbose output.

✓ GOOD

Justification: The described procedure is accurate and aligns with tested workflows.

* notes_*.txt: The test logs for Perlmutter, Frontier, Polaris, and other systems all follow the pattern of changing into the Exec/ABL/ subdirectory of the build directory and executing the compiled binary with a system-specific launcher (srun, mpiexec).
* CMakeLists.txt: The install(TARGETS ... RUNTIME DESTINATION bin) command confirms that make install places executables in a bin directory under the installation prefix. The check for if(ERF_ENABLE_TESTS) followed by include(CTest) and add_subdirectory(Tests) confirms the procedure for enabling and running the test suite.

The guidance is practical and directly reflects documented, verified procedures.


--------------------------------------------------------------------------------


5. Getting Help

Community support and project development depend on high-quality, reproducible bug reports. A well-formed issue enables developers to efficiently diagnose and resolve problems for the entire user community. This section provides a guide to contributing effectively to the project by submitting clear and complete issues.

5.1. Before Submitting an Issue

Please search the existing GitHub Issues on the project repository to see if your problem has already been reported. This helps avoid duplicate reports and may provide you with an immediate solution.

✓ GOOD

Justification: This is standard best practice for any open-source project and helps manage developer resources effectively.

5.2. Information to Provide in Your GitHub Issue

If your issue has not been reported, please create a new issue on the ERF GitHub repository. Include the following information to ensure we can reproduce and diagnose the problem effectively:

* The exact build command used. This includes the full cmake command with all flags, or the contents of the build script you executed.
* A description of the system. Include the operating system, compiler versions, and the complete output of the module list command.
* The full, unedited error message. Please do not summarize or screenshot the error. Instead, save the complete terminal output to a text file and attach it to the issue.
* The CMakeCache.txt file. This file is located in your build directory and contains the final configuration that CMake used. Attaching it provides a complete snapshot of your build environment.

✓ GOOD

Justification: This list represents the essential set of information required for effective remote debugging of a complex build system. Each item is critical for reproducing and understanding the user's environment and the specific failure mode.


--------------------------------------------------------------------------------


6. Contributing Fixes

Contributions from the user community are a vital part of a healthy open-source ecosystem. Users who investigate a problem and discover a solution are strongly encouraged to share their work by submitting a Pull Request. Improvements to the build system, additions of new machine-specific profiles, and enhancements to documentation are particularly valuable contributions.

6.1. Contribution Encouragement

If you discover a solution to a build problem, we strongly encourage you to share it with the community by submitting a Pull Request. Community contributions are highly valued and are critical to the sustainability of the project. Improvements to the build system, additions of new machine profile scripts, and enhancements to this troubleshooting documentation are especially welcome.

✓ GOOD

Justification: This is a positive and encouraging call to action that fosters a collaborative open-source environment.
