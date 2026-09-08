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
   * CUDA toolkit and ``cray-mpich`` versions incompatible → See :ref:`troubleshoot-cuda-version`
   * NetCDF/HDF5 not found → See :ref:`sec:build:library`
   * Wrong compiler detected → Check ``module list``

   **Quick checks:**

   .. code-block:: bash

      module list
      echo $CRAY_ACCEL_TARGET
      echo $NETCDF_DIR

.. dropdown:: Cray GPU link error: ``cannot find -lcudart``
   :icon: info
   :color: warning

   If configure/link fails with ``cannot find -lcudart``, check what the Cray wrapper is injecting:

   .. code-block:: bash

      CC --cray-print-opts=libs | grep -E 'cuda|mpi_gtl|mpich|libsci'

   If CUDA ``-L`` paths are missing (or stale), reload your machine profile/module stack and reconfigure from a clean build directory.

.. dropdown:: Cray GPU link error: ``undefined reference to cudaMalloc@libcudart.so.NN``
   :icon: info
   :color: warning

   A versioned ``@libcudart.so.NN`` suffix, or a ``libcudart.so.NN ... not found`` warning naming ``libmpi_gtl_cuda.so``, means the loaded CUDA toolkit and ``cray-mpich`` disagree on the CUDA runtime major version. CMake reports this as a broken C++ compiler.

   .. code-block:: bash

      readelf -d $CRAY_MPICH_ROOTDIR/gtl/lib/libmpi_gtl_cuda.so | grep -o 'libcudart\.so\.[0-9]*'
      ls $CUDA_HOME/lib64/libcudart.so.*.*

   A major-version disagreement between those two is the fault → See :ref:`troubleshoot-cuda-version`

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

.. _troubleshoot-cuda-version:

CUDA Toolkit / Cray MPICH Version Mismatch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Symptom:** On a Cray GPU build, CMake reports that the C++ compiler itself is broken, or Kokkos fails to compile. Neither message mentions modules, so both are easy to misread as a broken toolchain or an upstream bug.

**Error:** the mismatch surfaces in one of two places, depending on which direction it goes.

If the CUDA toolkit is *older* than the Cray MPICH GTL library expects, configure dies at ``project()``:

.. code-block:: text

   ld: warning: libcudart.so.13, needed by .../libmpi_gtl_cuda.so, not found
   ld: .../libmpi_gtl_cuda.so: undefined reference to `cudaMalloc@libcudart.so.13'
   CMake Error: The C++ compiler ".../CC" is not able to compile a simple test program.

If the CUDA toolkit is *newer* than the vendored Kokkos supports, configure succeeds and the compile fails instead:

.. code-block:: text

   Kokkos_Cuda_Instance.hpp: error: argument of type "size_t" is incompatible
       with parameter of type "const cudaGraphEdgeData *"
   Kokkos_Cuda_Instance.hpp: error: no suitable constructor exists to convert
       from "int" to "cudaMemLocation"

**Cause:** two independent version constraints must hold at once.

* Cray MPICH's GPU Transport Layer (``libmpi_gtl_cuda.so``) records a specific CUDA runtime soname as a hard ``DT_NEEDED`` dependency. The wrapper links it unconditionally (not ``--as-needed``), so that soname must resolve even for a trivial test program. ``cray-mpich/9.1.0`` and later need ``libcudart.so.13``; ``8.1.28`` through ``9.0.1`` need ``libcudart.so.12``. A soname major-version mismatch cannot be papered over with ``-L``, ``-rpath``, or any other link flag.
* The Kokkos vendored under ``Submodules/ekat/extern/kokkos`` (4.5.1) calls ``cudaGraphAddDependencies``, ``cudaMemAdvise``, and ``cudaMemPrefetchAsync`` with their pre-CUDA-13 signatures and no ``CUDA_VERSION`` guards. CUDA 13 changed all three, adopting the ``_v2`` forms that were introduced alongside them during 12.x.

Together these mean an EKAT-enabled build needs CUDA 12, which in turn restricts Cray MPICH to a ``libcudart.so.12`` build.

**Diagnosis:** three commands establish whether this is the problem.

.. code-block:: bash

   # 1. Which CUDA runtime does the loaded MPI's GTL library demand?
   readelf -d $CRAY_MPICH_ROOTDIR/gtl/lib/libmpi_gtl_cuda.so | grep -o 'libcudart\.so\.[0-9]*'

   # 2. Which CUDA runtime is actually loaded?
   ls $CUDA_HOME/lib64/libcudart.so.*.*

   # 3. Confirm the dependency is genuinely unresolvable
   ldd -r $CRAY_MPICH_ROOTDIR/gtl/lib/libmpi_gtl_cuda.so | grep 'not found'

If (1) and (2) disagree in major version, that is the fault. These two loops then print the full compatibility matrix for the system, so a working pair can be read off directly:

.. code-block:: bash

   # CUDA runtimes installed on the system
   for d in /opt/nvidia/hpc_sdk/Linux_x86_64/*/cuda/*/lib64; do
     printf '%s: ' "$d"
     ls $d/libcudart.so.*.* 2>/dev/null | xargs -r -n1 basename | tr '\n' ' '; echo
   done

   # CUDA runtime each installed cray-mpich demands
   for g in /opt/cray/pe/mpich/*/gtl/lib/libmpi_gtl_cuda.so; do
     printf '%s -> ' "$g"; readelf -d "$g" | grep -o 'libcudart\.so\.[0-9]*'
   done

**Solution:** source the machine profile, which pins a known-good pair:

.. code-block:: bash

   source Build/machines/perlmutter_erf.profile

Or set the pair by hand and verify it before configuring:

.. code-block:: bash

   module load cudatoolkit/12.9
   module swap cray-mpich cray-mpich/9.0.1

   # expect libcudart.so.12, matching cudatoolkit/12.9
   readelf -d $CRAY_MPICH_ROOTDIR/gtl/lib/libmpi_gtl_cuda.so | grep -o 'libcudart\.so\.[0-9]*'

Then reconfigure in a **fresh** build directory, not the one that failed — ERF's Cray detection writes ``CMAKE_*_STANDARD_LIBRARIES`` into the cache with ``FORCE``, so those values do not re-derive on a reconfigure. See :ref:`troubleshoot-cache`.

.. note::
   The default ``cray-mpich`` advances with the CPE release, so a system default that worked before an upgrade can start requiring a newer CUDA runtime than the pinned toolkit provides. This is why the machine profiles pin both halves explicitly; scripts that load ``cray-mpich`` without a version are exposed to the default moving underneath them.

.. tip::
   If a build does not need EKAT/SHOC, ``-DERF_ENABLE_EKAT=OFF`` removes the Kokkos constraint entirely and allows the newer CUDA plus newer Cray MPICH pairing. AMReX guards these same three CUDA APIs on ``CUDART_VERSION >= 13000``, so the AMReX-only path builds against either toolkit.

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
   cd build/install/bin
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
