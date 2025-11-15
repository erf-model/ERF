.. _sec:build:troubleshooting:

=================================
ERF Build Troubleshooting Guide
=================================

1.0 Resolved Issues (Historical)
--------------------------------

Documenting resolved issues demonstrates the build system's increasing robustness by chronicling historical complexities now handled automatically. This serves as a diagnostic tool for advanced users. By understanding automated fixes for known issues on platforms like Cray, users can more effectively diagnose new problems if automation is not functioning as expected.

**1.1 Cray Compiler and GPU-Aware MPI Linker Errors (Fix 4)**

On Cray systems with GPU-aware MPI enabled (``MPICH_GPU_SUPPORT_ENABLED=1``), historical builds failed during final linking stage. Root cause was Cray linker's default use of ``--as-needed`` flag, which aggressively discards libraries deemed unused during initial pass. This process incorrectly removed essential MPI GPU Transport Layer (GTL) libraries required for CPU-GPU communication.

Build system now automatically detects this configuration. It identifies correct MPI base library (e.g., ``mpi_gnu_123``) and required GTL library for target hardware (e.g., ``mpi_gtl_cuda`` for NVIDIA GPUs or ``mpi_gtl_hsa`` for AMD GPUs). These libraries are programmatically added to ``CMAKE_CXX_STANDARD_LIBRARIES`` and ``CMAKE_CUDA_STANDARD_LIBRARIES`` variables.

* **Phase 1 (Pre-Project)**: Logic in ``CrayCompilerDetection.cmake`` runs before CMake's ``project()`` command. Performs initial detection of GPU-aware MPI and appends required GTL libraries to standard library variables.
* **Phase 2 (Post-Project)**: After project configuration, logic in ``CrayDetection.cmake`` (Fix 4) confirms and applies this fix.

This two-phase strategy addresses both CMake's procedural requirements and complexities of Cray environment.

**1.2 NetCDF Parallel Library Detection on Cray Systems (Fix 5-6)**

Standard CMake ``find_package`` calls often fail on Cray systems because libraries like parallel NetCDF and HDF5 are not in default system paths. Their locations are managed by environment modules. Previously required manual hints via ``NETCDF_DIR`` variable.

Build system now queries Cray compiler wrapper directly. Executes ``CC --cray-print-opts=PKG_CONFIG_PATH`` to discover correct search path from loaded Cray environment modules. This path is prepended to ``PKG_CONFIG_PATH`` environment variable, guiding ``pkg-config`` utility to find correct parallel library installations.

**1.3 Manual GPU Architecture Specification for Kokkos/AMReX**

Previously, users manually specified target GPU architecture for all dependencies. This error-prone process could lead to compilation failures or non-functional executables.

Now fully automated. Build system inspects ``$CRAY_ACCEL_TARGET`` environment variable (e.g., ``nvidia80`` for NVIDIA A100 or ``amd_gfx90a`` for AMD MI250X), set when loading appropriate ``craype-accel-*`` module. Variable is used to programmatically set correct GPU architecture flags for all dependencies.

2.0 Troubleshooting Common Issues
---------------------------------

While build system automates many platform-specific complexities, successful compilation depends on correctly configured user environment. This section addresses frequent user-side configuration errors and provides actionable solutions.

**2.1 Missing craype-accel-* Module on GPU Builds**

For CUDA or HIP builds on Cray systems, Cray compiler wrappers require ``craype-accel-*`` module. This sets essential ``$CRAY_ACCEL_TARGET`` environment variable.

ERF build system checks for this variable before ``project()`` command. If undefined during GPU-enabled build, CMake issues fatal error with instructions.

**Recommended Solution**: Load appropriate module for system's hardware (e.g., ``craype-accel-nvidia80`` for NVIDIA A100 GPUs). Most reliable method is sourcing correct environment script from ``Build/machines/`` directory before running CMake.

**2.2 NetCDF or HDF5 Libraries Not Found**

Enabling NetCDF I/O with ``-DERF_ENABLE_NETCDF=ON`` causes build system to set ``-DERF_ENABLE_HDF5=ON`` by default. ``SetAmrexOptions.cmake`` script reads ``ERF_ENABLE_HDF5`` value and sets ``AMReX_HDF5=ON``. Build fails if necessary libraries cannot be located.

**Recommended Solution**: On Cray systems, load ``cray-netcdf-hdf5parallel`` module. Cray detection logic in ERF build system locates and configures parallel libraries from this module.

For non-Cray or custom builds, point build system to NetCDF installation by setting ``NETCDF_DIR`` environment variable or CMake variable to NetCDF installation directory root.

**2.3 Memory Issues During Compilation**

**Symptom:** CUDA or HIP compilation fails with out-of-memory errors, particularly when building GPU-enabled physics packages (RRTMGP, SHOC, P3).

**Common Error Messages:**

.. code-block:: none

   nvcc fatal: Memory allocation failure
   c++: fatal error: Killed signal terminated program cc1plus
   internal compiler error: Killed (program cc1plus)

**Root Cause:** On HPC systems that allow partial node allocations (e.g., Kestrel, some SLURM configurations), the default memory assignment may be insufficient for GPU compilation, which is significantly more memory-intensive than CPU compilation.

**Solution:**

When requesting compute nodes for compilation, explicitly request sufficient memory:

**Option 1: Request exclusive node access** (recommended for HPC systems)

.. code-block:: bash

   # In your SLURM script or interactive session request
   #SBATCH --exclusive

This allocates the entire node's memory to your job.

**Option 2: Request specific memory amount**

.. code-block:: bash

   # Request total memory for the allocation
   #SBATCH --mem=240G

   # OR request memory per CPU
   #SBATCH --mem-per-cpu=4G

**Option 3: Compile on login nodes** (if permitted by site policy)

Some sites allow compilation on login nodes, which typically have more available memory. Check your site's policy before using this approach.

**Prevention:**

For development workflows involving frequent recompilation:

1. Use ``make -j4`` instead of ``make -j`` to limit parallel compilation jobs
2. Compile incrementally rather than from scratch when possible
3. Request an interactive session with adequate memory for the build

.. note::
   This issue is particularly common on systems like Kestrel where partial node allocations are the default. Always include ``--exclusive`` or explicit memory requests in your build job scripts.

3.0 General Debugging Strategies
--------------------------------

This section provides guidance for diagnosing novel or complex issues not covered by specific troubleshooting cases above.

**3.1 CMake Debugging**

ERF's CMake configuration includes logging features to trace build process:

* **Verbose Logging**: Use ``--log-level=VERBOSE`` flag with cmake command (e.g., ``cmake --log-level=VERBOSE ..``). Prints detailed information including library search paths.
* **Contextual Logging**: Use ``--log-context`` flag to display component hierarchy for each message (e.g., ``[ERF.Cray]``).
* **Cache Inspection**: After configuration attempt, ``CMakeCache.txt`` file in build directory contains final values of all configuration variables.

**3.2 GNU Make Debugging**

GNU Make offers low-level, imperative control model. Preferred system for experts requiring direct command over every compiler flag and library path. AMReX provides utility for inspecting makefile variables:

* To check variable value, use ``make print-<variable_name>``. For example, ``make print-CXXFLAGS`` displays exact C++ compiler flags being used.

**3.3 Compiler and Linker Debugging**

When build fails during compilation or linking:

* **Verbose Build**: Pass verbose flag to build command (e.g., ``make VERBOSE=1``). Prints full compiler and linker commands.
* **Shared Library Dependencies**: On Linux systems, ``ldd`` command on compiled executable checks for missing shared library dependencies.

4.0 Verifying a Successful Build
--------------------------------

Successful compilation signals syntactically correct code but does not guarantee functional executable. Final verification is essential to confirm application initializes correctly, runs without crashing, and produces valid results. This validates all libraries were linked correctly and are compatible with runtime environment.

**4.1 Procedure**

1. **Navigate to Example Case**: After successful build, locate executable. ``make install`` places it in ``install/bin/`` directory relative to build directory. Simple ``make`` places it in ``Exec/`` subdirectory within build directory.

2. **Execute Short Run**: Run example case for small number of timesteps to confirm simulation can initialize and advance. MPI launcher command varies by system:

.. code:: shell

   mpiexec -n 4 ./erf_abl inputs max_step=10

3. **Run Test Suite (Optional)**: For comprehensive verification, run built-in regression tests. Configure build with ``-DERF_ENABLE_TESTS=ON``, then execute from build directory:

.. code:: shell

   ctest -L regression -VV

5.0 Getting Help
----------------

Community support depends on high-quality, reproducible bug reports. Well-formed issue enables developers to efficiently diagnose and resolve problems.

**5.1 Before Submitting an Issue**

Search existing GitHub Issues on project repository to see if problem has been reported. This avoids duplicate reports and may provide immediate solution.

**5.2 Information to Provide in Your GitHub Issue**

If issue has not been reported, create new issue on ERF GitHub repository. Include following information:

* Exact build command used (full cmake command with all flags, or contents of build script)
* System description (operating system, compiler versions, complete output of ``module list``)
* Full, unedited error message (save complete terminal output to text file and attach to issue)
* ``CMakeCache.txt`` file from build directory (provides complete snapshot of build environment)

6.0 Contributing Fixes
----------------------

Contributions from user community are vital part of open-source ecosystem. Users who discover solutions are encouraged to share their work by submitting Pull Request.

**6.1 Contribution Encouragement**

If you discover solution to build problems, share it with community by submitting a Pull Request. Community contributions are valued and necessary for project sustainability. Improvements to the build system, additions of new machine profile scripts, and enhancements to troubleshooting documentation are especially welcome.
