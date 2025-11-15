.. _sec:build:troubleshooting:

Build Troubleshooting
======================

This guide helps diagnose and resolve common ERF build issues. For library-specific configuration issues, see :ref:`sec:build:library`. For HPC-specific problems, see :ref:`sec:build:hpc`.

Quick Diagnostic Guide
----------------------

**Where's the problem?**

.. dropdown:: Build won't start - CMake configuration fails
   :icon: alert

   **Common causes:**

   * Missing ``craype-accel-*`` module on Cray GPU builds → See :ref:`troubleshoot-cray-accel`
   * NetCDF/HDF5 not found → See :ref:`troubleshoot-netcdf`
   * Wrong compiler detected → Check ``module list`` and reload environment

   **Quick checks:**

   .. code-block:: bash

      # Verify modules loaded
      module list

      # Check environment variables
      echo $CRAY_ACCEL_TARGET  # Should be set for GPU builds
      echo $NETCDF_DIR         # Should point to NetCDF install

.. dropdown:: Build starts but compilation fails
   :icon: alert

   **Common causes:**

   * Out of memory during CUDA/HIP compilation → See :ref:`troubleshoot-memory`
   * Incompatible compiler flags → Check ``make print-CXXFLAGS``
   * Missing source files → Verify submodules: ``git submodule update --init --recursive``

   **Quick checks:**

   .. code-block:: bash

      # Check available memory
      free -h

      # Reduce parallel jobs
      make -j4  # instead of make -j

.. dropdown:: Linking fails
   :icon: alert

   **Common causes:**

   * MPI linker errors → Parallel/serial library mismatch, see :ref:`troubleshoot-mpi-linker`
   * GPU-aware MPI on Cray → Now auto-fixed, see :ref:`resolved-gpu-mpi`
   * Missing libraries → Check ``ldd`` on executable

   **Quick checks:**

   .. code-block:: bash

      # Check if parallel libraries loaded
      module -t list | grep -i parallel

.. dropdown:: Build succeeds but executable fails
   :icon: alert

   **Common causes:**

   * Wrong GPU architecture compiled → Check ``$CRAY_ACCEL_TARGET``
   * Missing runtime libraries → Use ``ldd`` to check dependencies
   * MPI misconfiguration → Try simple MPI hello world

   **Verification:**

   .. code-block:: bash

      # Check library dependencies
      ldd ./ERF3d.*.ex

      # Try short test run
      mpiexec -n 4 ./ERF3d.*.ex inputs max_step=10

.. _troubleshooting-common:

Common Issues
-------------

.. _troubleshoot-cray-accel:

Missing craype-accel Module (Cray GPU Builds)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Symptom:** CMake fatal error during GPU build on Cray systems.

**Error Message:**

.. code-block:: text

   CMake Error: CRAY_ACCEL_TARGET not set for GPU build

**Cause:** Cray compiler wrappers require ``craype-accel-*`` module for GPU builds. This sets ``$CRAY_ACCEL_TARGET`` environment variable.

**Solution:**

Load the appropriate module for your hardware:

.. code-block:: bash

   # NVIDIA A100 (Perlmutter, Polaris)
   module load craype-accel-nvidia80

   # AMD MI250X (Frontier)
   module load craype-accel-amd-gfx90a

   # Intel GPUs (Aurora)
   module load craype-accel-intel-gpu

**Best Practice:** Source the machine profile before building:

.. code-block:: bash

   source Build/machines/perlmutter_erf.profile
   cmake -DERF_ENABLE_CUDA=ON ..

.. _troubleshoot-netcdf:

NetCDF or HDF5 Not Found
~~~~~~~~~~~~~~~~~~~~~~~~~

**Symptom:** CMake cannot locate NetCDF or HDF5 libraries.

**Error Message:**

.. code-block:: text

   Could not find NetCDF
   Could not find HDF5

**Cause:** Libraries not in standard search paths or parallel versions not loaded.

**Solution:**

.. tab-set::

   .. tab-item:: Cray Systems

      Load the parallel NetCDF module:

      .. code-block:: bash

         module load cray-netcdf-hdf5parallel
         cmake -DERF_ENABLE_NETCDF=ON ..

      The Cray detection system automatically finds the libraries.

   .. tab-item:: Custom Installation

      Specify the NetCDF installation path:

      .. code-block:: bash

         # Via CMake option
         cmake -DERF_ENABLE_NETCDF=ON \
               -DNETCDF_DIR=/opt/netcdf-parallel \
               ..

         # Or via environment variable
         export NETCDF_DIR=/opt/netcdf-parallel
         cmake -DERF_ENABLE_NETCDF=ON ..

   .. tab-item:: Workstation

      Install parallel versions:

      .. code-block:: bash

         # Ubuntu/Debian
         sudo apt install libnetcdf-mpi-dev libhdf5-mpi-dev

         # Then build
         cmake -DERF_ENABLE_NETCDF=ON ..

For detailed NetCDF detection logic, see :ref:`sec:build:library`.

.. _troubleshoot-mpi-linker:

MPI Linker Errors
~~~~~~~~~~~~~~~~~

**Symptom:** Linker errors with undefined MPI symbols.

**Error Message:**

.. code-block:: text

   undefined reference to `MPI_Init'
   undefined reference to `MPI_Comm_rank'
   undefined reference to `MPI_Finalize'

**Cause:** Parallel ERF build linking against serial NetCDF/HDF5 libraries.

**Diagnosis:**

.. code-block:: bash

   # Check which NetCDF is loaded
   module list | grep -i netcdf

   # Should show "parallel" in the name
   # Bad:  netcdf/4.9.0
   # Good: cray-netcdf-hdf5parallel/4.9.0.3

**Solution:**

1. Unload serial libraries and load parallel versions:

   .. code-block:: bash

      module unload netcdf hdf5
      module load cray-netcdf-hdf5parallel

2. Clean and rebuild:

   .. code-block:: bash

      make distclean
      cmake ..
      make

For more details, see :ref:`sec:build:library` Section "Parallel I/O Consistency".

.. _troubleshoot-memory:

Memory Issues During Compilation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Symptom:** Compilation fails with out-of-memory errors, especially during CUDA/HIP builds.

**Error Messages:**

.. code-block:: text

   nvcc fatal: Memory allocation failure
   c++: fatal error: Killed signal terminated program cc1plus
   internal compiler error: Killed (program cc1plus)

**Cause:** GPU compilation is memory-intensive. On systems allowing partial node allocations, default memory may be insufficient.

**Solution:**

.. tab-set::

   .. tab-item:: Request Exclusive Node

      Recommended for HPC systems:

      .. code-block:: bash

         # In SLURM script
         #SBATCH --exclusive

         # Or interactive session
         salloc --exclusive -N 1

   .. tab-item:: Request Specific Memory

      .. code-block:: bash

         # Total memory
         #SBATCH --mem=240G

         # Or per CPU
         #SBATCH --mem-per-cpu=4G

   .. tab-item:: Limit Parallel Jobs

      .. code-block:: bash

         # Instead of make -j (unlimited)
         make -j4  # Limit to 4 parallel jobs

.. note::
   This is particularly common on Kestrel where partial node allocations are default. Always include ``--exclusive`` or explicit memory requests in build scripts.

**Stale CMake Cache**
~~~~~~~~~~~~~~~~~~~~~

**Symptom:** Unexpected failures after changing modules or compilers.

**Cause:** CMake caches library locations. When the environment changes, the cache becomes stale.

**Solution:**

.. code-block:: bash

   # Full clean
   make distclean

   # Or manually remove cache
   rm -rf CMakeCache.txt CMakeFiles/

   # Reconfigure
   cmake ..

.. _resolved-issues:

Resolved Issues (Automated)
----------------------------

These issues are now handled automatically by the build system. This section documents them for advanced users debugging automation failures.

.. _resolved-gpu-mpi:

Cray GPU-Aware MPI Linking
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Historical Problem:** On Cray systems with GPU-aware MPI (``MPICH_GPU_SUPPORT_ENABLED=1``), builds failed during linking. The Cray linker's ``--as-needed`` flag incorrectly removed MPI GPU Transport Layer (GTL) libraries.

**Automated Solution:** The build system now:

1. **Phase 1 (Pre-Project):** Detects GPU-aware MPI in ``CrayCompilerDetection.cmake``
2. Identifies correct MPI base library (e.g., ``mpi_gnu_123``)
3. Identifies required GTL library:
   
   - ``mpi_gtl_cuda`` for NVIDIA GPUs
   - ``mpi_gtl_hsa`` for AMD GPUs

4. **Phase 2 (Post-Project):** Confirms and applies fix in ``CrayDetection.cmake``
5. Adds libraries to ``CMAKE_CXX_STANDARD_LIBRARIES`` and ``CMAKE_CUDA_STANDARD_LIBRARIES``

**If automation fails:** Check that ``MPICH_GPU_SUPPORT_ENABLED=1`` is set and appropriate ``craype-accel-*`` module is loaded.

**NetCDF/HDF5 Detection on Cray**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Historical Problem:** Standard ``find_package`` calls failed because parallel NetCDF/HDF5 are in module-managed paths, not system paths.

**Automated Solution:** The build system queries the Cray compiler wrapper:

.. code-block:: bash

   CC --cray-print-opts=PKG_CONFIG_PATH

This discovers the correct search path from loaded modules and prepends it to ``PKG_CONFIG_PATH``, guiding ``pkg-config`` to find parallel libraries.

**If automation fails:** Manually load ``cray-netcdf-hdf5parallel`` and verify with:

.. code-block:: bash

   pkg-config --list-all | grep netcdf

**GPU Architecture Auto-Detection**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Historical Problem:** Users manually specified GPU architecture for all dependencies, leading to errors and incompatibilities.

**Automated Solution:** The build system:

1. Inspects ``$CRAY_ACCEL_TARGET`` (e.g., ``nvidia80`` for A100, ``amd_gfx90a`` for MI250X)
2. Maps to correct architecture flags for all dependencies:
   
   - AMReX: ``AMReX_CUDA_ARCH``, ``AMReX_AMD_ARCH``
   - Kokkos: ``Kokkos_ARCH_AMPERE80``, ``Kokkos_ARCH_VEGA90A``

**If automation fails:** Check that ``craype-accel-*`` module is loaded and ``$CRAY_ACCEL_TARGET`` is set:

.. code-block:: bash

   echo $CRAY_ACCEL_TARGET

Debugging Tools
---------------

CMake Debugging
~~~~~~~~~~~~~~~

Use CMake's logging features to trace the build process:

.. code-block:: bash

   # Verbose logging
   cmake --log-level=VERBOSE ..

   # With context (shows component hierarchy)
   cmake --log-context --log-level=VERBOSE ..

   # Debug level (even more detail)
   cmake --log-level=DEBUG ..

**Example output:**

.. code-block:: text

   [ERF.Cray] Detected Cray Programming Environment (CRAYPE_VERSION=2.7.30)
   [ERF.Cray] Setting Cray compiler wrappers...
   [ERF.AMReX] Using internal AMReX submodule
   [ERF.NetCDF] Found NetCDF: /opt/cray/pe/netcdf/4.9.0.9

**Inspect cache:**

.. code-block:: bash

   # View all cached variables
   cmake -LAH | less

   # Or examine CMakeCache.txt directly
   grep NETCDF CMakeCache.txt

GNU Make Debugging
~~~~~~~~~~~~~~~~~~

Use AMReX's built-in utilities:

.. code-block:: bash

   # Print specific variable value
   make print-CXXFLAGS
   make print-LIBRARIES
   make print-INCLUDE_LOCATIONS

   # Verbose build (show more)
   make VERBOSE=1

   # help build (show full variables)
   make help

Compiler and Linker Debugging
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Verbose compilation:**

.. code-block:: bash

   # CMake
   make VERBOSE=1

   # GNU Make (already verbose by default)
   make

**Check library dependencies:**

.. code-block:: bash

   # After build succeeds
   ldd ./ERF3d.*.ex

   # Check for specific library
   ldd ./ERF3d.*.ex | grep netcdf

**Check symbols in executable:**

.. code-block:: bash

   # Check for NetCDF symbols
   nm ERF3d.*.ex | grep nc_

   # Check for MPI symbols
   nm ERF3d.*.ex | grep MPI_

Verifying a Successful Build
-----------------------------

A successful compilation doesn't guarantee a functional executable. Final verification ensures correct initialization and execution.

**Quick Verification**
~~~~~~~~~~~~~~~~~~~~~~

Run a short test to confirm the executable works:

.. code-block:: bash

   # Navigate to executable location
   cd build/install/bin  # For cmake --install
   # or
   cd Exec/ABL  # For in-place builds

   # Run short simulation
   mpiexec -n 4 ./ERF3d.*.ex inputs max_step=10

**Comprehensive Verification**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run the regression test suite:

.. code-block:: bash

   # Configure with tests enabled
   cmake -DERF_ENABLE_TESTS=ON ..
   make

   # Run regression tests
   cd build
   ctest -L regression -VV

**Check Build Configuration**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

View the build configuration:

.. code-block:: bash

   ./ERF3d.*.ex --describe

This shows:

* AMReX version
* Compiler versions
* Enabled features (MPI, GPU, NetCDF, etc.)
* GPU architecture (if applicable)

Getting Help
------------

**Before Submitting an Issue**

1. Search `existing GitHub Issues <https://github.com/erf-model/ERF/issues>`_
2. Check this troubleshooting guide and :ref:`sec:build:library`
3. Try the diagnostic commands in this guide

**Creating a Good Bug Report**

If the issue hasn't been reported, create a new issue with:

**Required Information:**

.. code-block:: text

   **System Information:**
   - OS: [e.g., Ubuntu 22.04, Perlmutter/CrayOS]
   - Compiler: [output of `gcc --version` or `CC --version`]
   - MPI: [output of `mpirun --version`]
   - Modules: [output of `module list`]

   **Build Command:**
   [Paste your complete cmake command or build script]

   **Error Message:**
   [Complete, unedited terminal output]

**Attach Files:**

* ``CMakeCache.txt`` from build directory
* Full build log (redirect output: ``make 2>&1 | tee build.log``)

**Diagnostic Information:**

Run and include output from:

.. code-block:: bash

   # CMake diagnostic
   cmake --log-level=VERBOSE --log-context .. 2>&1 | tee cmake_verbose.log

   # Environment check
   echo $CRAY_ACCEL_TARGET
   echo $NETCDF_DIR
   module list

Contributing Fixes
------------------

If you discover a solution to a build problem, contribute it back to the community!

**What to Contribute:**

* New machine profiles (``Build/machines/*.profile``)
* Build system improvements
* Documentation enhancements
* Troubleshooting examples

**How to Contribute:**

1. Fork the `ERF repository <https://github.com/erf-model/ERF>`_
2. Create a feature branch
3. Make your changes
4. Submit a Pull Request

See the project's contribution guidelines for details.

.. note::
   Community contributions are essential for ERF's sustainability. Your solutions help other users and improve the project for everyone.
