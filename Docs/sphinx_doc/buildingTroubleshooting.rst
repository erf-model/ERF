.. _sec:build:troubleshooting:

Build Troubleshooting
======================

This guide helps diagnose and resolve ERF build issues. For library configuration problems, see :ref:`sec:build:library`. For HPC-specific issues, see :ref:`sec:build:hpc`.

Quick Diagnostic
----------------

**Where's the problem?**

.. dropdown:: CMake configuration fails
   :icon: alert
   :color: warning

   **Common causes:**

   * Missing ``craype-accel-*`` module on Cray GPU builds → See :ref:`troubleshoot-cray-accel`
   * NetCDF/HDF5 not found → See :ref:`sec:build:library`
   * Wrong compiler detected → Check ``module list``

   **Quick checks:**

   .. code-block:: bash

      module list
      echo $CRAY_ACCEL_TARGET
      echo $NETCDF_DIR

.. dropdown:: Compilation fails
   :icon: alert
   :color: warning

   **Common causes:**

   * Out of memory during CUDA/HIP compilation → See :ref:`troubleshoot-memory`
   * Missing source files → Check ``git submodule update --init --recursive``
   * Stale CMake cache → See :ref:`troubleshoot-cache`

   **Quick fixes:**

   .. code-block:: bash

      # Check memory
      free -h

      # Update submodules
      git submodule update --init --recursive

      # Reduce parallel jobs
      make -j4

.. dropdown:: Linking fails
   :icon: alert
   :color: warning

   **Common causes:**

   * Parallel/serial library mismatch → See :ref:`sec:build:library` for MPI linker errors
   * Missing libraries → Check ``ldd ./ERF3d.*.ex``
   * GPU-aware MPI issues (now auto-fixed on Cray)

.. dropdown:: Executable fails to run
   :icon: alert
   :color: warning

   **Common causes:**

   * Wrong GPU architecture
   * Missing runtime libraries
   * MPI misconfiguration

   **Verification:**

   .. code-block:: bash

      # Check dependencies
      ldd ./ERF3d.*.ex

      # Try short run
      mpiexec -n 4 ./ERF3d.*.ex inputs max_step=10

   → See :ref:`verification` for full verification steps

Build Process Issues
--------------------

.. _troubleshoot-cray-accel:

Missing craype-accel Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Symptom:** CMake error during GPU build on Cray systems.

**Error:**

.. code-block:: text

   CMake Error: CRAY_ACCEL_TARGET not set for GPU build

**Cause:** GPU builds on Cray require ``craype-accel-*`` module to set ``$CRAY_ACCEL_TARGET``.

**Solution:**

Load the module for your hardware:

.. code-block:: bash

   # NVIDIA A100 (Perlmutter, Polaris)
   module load craype-accel-nvidia80

   # AMD MI250X (Frontier)
   module load craype-accel-amd-gfx90a

   # Intel GPUs (Aurora)
   module load craype-accel-intel-gpu

**Best practice:** Use machine profiles:

.. code-block:: bash

   source Build/machines/perlmutter_erf.profile
   cmake -DERF_ENABLE_CUDA=ON ..

.. _troubleshoot-memory:

Out of Memory During Compilation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Symptom:** Compilation killed with memory errors.

**Error:**

.. code-block:: text

   nvcc fatal: Memory allocation failure
   c++: fatal error: Killed signal terminated program

**Cause:** GPU compilation requires more memory than default allocation on partial-node systems.

**Solution:**

.. tab-set::

   .. tab-item:: Exclusive Node

      .. code-block:: bash

         # SLURM script
         #SBATCH --exclusive

         # Interactive
         salloc --exclusive -N 1

   .. tab-item:: Specific Memory

      .. code-block:: bash

         #SBATCH --mem=240G
         # or
         #SBATCH --mem-per-cpu=4G

   .. tab-item:: Limit Parallel Jobs

      .. code-block:: bash

         make -j4  # Instead of make -j

.. note::
   Common on Kestrel where partial node allocations are default. Always use ``--exclusive`` or explicit memory requests.

.. _troubleshoot-cache:

Stale CMake Cache
~~~~~~~~~~~~~~~~~

**Symptom:** Unexpected failures after changing modules or compilers.

**Cause:** CMake caches library locations that become invalid when environment changes.

**Solution:**

.. code-block:: bash

   make distclean
   cmake ..
   make

Or manually:

.. code-block:: bash

   rm -rf CMakeCache.txt CMakeFiles/
   cmake ..

.. _debugging-tools:

Debugging Tools
---------------

CMake Debugging
~~~~~~~~~~~~~~~

.. code-block:: bash

   # Verbose output
   cmake --log-level=VERBOSE ..

   # With context (shows hierarchy)
   cmake --log-context --log-level=VERBOSE ..

**Example output:**

.. code-block:: text

   [ERF.Cray] Detected Cray Programming Environment
   [ERF.Cray] Setting Cray compiler wrappers...
   [ERF.NetCDF] Found NetCDF: /opt/cray/pe/netcdf/4.9.0.9

**Inspect cache:**

.. code-block:: bash

   cmake -LAH | less
   grep NETCDF CMakeCache.txt

GNU Make Debugging
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Print variable values
   make print-CXXFLAGS
   make print-LIBRARIES

   # Verbose build
   make VERBOSE=1

Library Dependencies
~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Check linked libraries
   ldd ./ERF3d.*.ex | grep netcdf

   # Check for symbols
   nm ERF3d.*.ex | grep nc_
   nm ERF3d.*.ex | grep MPI_

.. _verification:

Verifying Successful Builds
----------------------------

Quick Test
~~~~~~~~~~

.. code-block:: bash

   # Run short simulation
   cd build/install/bin  # or Exec/ABL for in-place builds
   mpiexec -n 4 ./ERF3d.*.ex inputs max_step=10

Regression Tests
~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Configure with tests
   cmake -DERF_ENABLE_TESTS=ON ..
   make

   # Run tests
   ctest -L regression -VV

Check Build Info
~~~~~~~~~~~~~~~~

.. code-block:: bash

   ./ERF3d.*.ex --describe

Shows compiler versions, enabled features, and GPU architecture.

.. dropdown:: Resolved Issues (Automated)
   :icon: check-circle
   :color: success

   These issues are now handled automatically by the build system.

   **Cray GPU-Aware MPI Linking**

   **Historical problem:** Linking failed with GPU-aware MPI due to Cray's ``--as-needed`` flag removing GTL libraries.

   **Automated solution:**

   1. Detects GPU-aware MPI (``MPICH_GPU_SUPPORT_ENABLED=1``)
   2. Identifies MPI base library (e.g., ``mpi_gnu_123``)
   3. Identifies required GTL library:
      - ``mpi_gtl_cuda`` for NVIDIA
      - ``mpi_gtl_hsa`` for AMD
   4. Adds to ``CMAKE_CXX_STANDARD_LIBRARIES`` and ``CMAKE_CUDA_STANDARD_LIBRARIES``

   **If automation fails:** Check ``MPICH_GPU_SUPPORT_ENABLED=1`` is set and ``craype-accel-*`` module loaded.

   **NetCDF/HDF5 Detection on Cray**

   **Historical problem:** ``find_package`` failed because parallel libraries in module-managed paths.

   **Automated solution:**

   Queries Cray compiler wrapper:

   .. code-block:: bash

      CC --cray-print-opts=PKG_CONFIG_PATH

   Prepends path to ``PKG_CONFIG_PATH``, enabling ``pkg-config`` to find parallel libraries.

   **If automation fails:** Load ``cray-netcdf-hdf5parallel`` manually.

   **GPU Architecture Auto-Detection**

   **Historical problem:** Users manually specified architecture for all dependencies.

   **Automated solution:**

   1. Reads ``$CRAY_ACCEL_TARGET`` (e.g., ``nvidia80``, ``amd_gfx90a``)
   2. Maps to architecture flags:
      - AMReX: ``AMReX_CUDA_ARCH=8.0``
      - Kokkos: ``Kokkos_ARCH_AMPERE80=ON``

   **If automation fails:** Check ``craype-accel-*`` module loaded:

   .. code-block:: bash

      echo $CRAY_ACCEL_TARGET

Getting Help
------------

**Before submitting an issue:**

1. Search `existing issues <https://github.com/erf-model/ERF/issues>`_
2. Check this guide and :ref:`sec:build:library`
3. Run diagnostic commands above

**Creating a bug report:**

Include this information in your `GitHub issue <https://github.com/erf-model/ERF/issues/new>`_:

.. code-block:: text

   **System:**
   - OS: [e.g., Perlmutter/CrayOS, Ubuntu 22.04]
   - Compiler: [gcc --version or CC --version]
   - MPI: [mpirun --version]
   - Modules: [module list]

   **Build command:**
   [Complete cmake command or script]

   **Error:**
   [Complete, unedited terminal output]

**Attach files:**

* ``CMakeCache.txt``
* Build log: ``make 2>&1 | tee build.log``

**Diagnostic output:**

.. code-block:: bash

   cmake --log-level=VERBOSE --log-context .. 2>&1 | tee cmake_verbose.log
   echo $CRAY_ACCEL_TARGET
   echo $NETCDF_DIR
   module list

Contributing Fixes
------------------

If you solve a build problem, contribute your solution!

**Contributions welcome:**

* Machine profiles (``Build/machines/*.profile``)
* Build system improvements
* Documentation enhancements
* Troubleshooting examples

**How to contribute:**

1. Fork `ERF repository <https://github.com/erf-model/ERF>`_
2. Create feature branch
3. Make changes
4. Submit Pull Request

See contribution guidelines in the repository.

.. note::
   Community contributions are essential. Your solutions help other users and improve ERF for everyone.
