.. _sec:build:hpc:

Building on HPC Systems
========================

Building ERF on High-Performance Computing (HPC) systems requires managing diverse architectures, proprietary compiler toolchains, and specialized environment modules. The ERF build system automates this complexity through machine profile files and the Cray Detection System.

This document explains ERF's HPC build strategy. For basic build instructions, see :ref:`sec:build:systems`. For initial setup, see :ref:`sec:build:quickstart`.

.. note::
   **HPC system environments change regularly.** Module names, versions, and
   availability vary as systems are updated. Always verify current module names
   with ``module avail`` before building. The profile scripts in
   ``Build/machines/`` are maintained for major systems but may require
   adjustment after system upgrades.

.. note::
   **CMake workflow:** The examples below use out-of-source builds (``mkdir build && cd build``), which keeps build artifacts separate from source code. This is recommended for HPC systems where you may need multiple build configurations or want to preserve a clean source tree. For other workflow options, see :ref:`sec:build:systems`.

Machine Profile Files
----------------------

Machine profile files prepare the shell environment for compilation by loading required software modules and setting environment variables. They are located in ``Build/machines/`` and provide a standardized, reproducible way to configure the build environment. See https://warpx.readthedocs.io/en/latest/install/hpc.html for detailed examples that include their module load and use environment hints for setting compilers.

**Purpose**

Profile files execute module load commands and export environment variables. On Cray systems, the core is the ``PrgEnv-*`` module (e.g., ``PrgEnv-gnu`` or ``PrgEnv-cray``), which dictates the compiler suite that Cray wrappers (``cc``, ``CC``) will use. This is supplemented by modules for hardware acceleration (e.g., ``craype-accel-nvidia80``), parallel libraries (e.g., ``cray-netcdf-hdf5parallel``), and development tools (e.g., ``cmake``).

**Usage**

Source the appropriate profile from your shell:

.. code-block:: bash

   source Build/machines/perlmutter_erf.profile

This modifies the current shell session. The shell is then ready for subsequent build commands.

**Available Profiles**

ERF provides pre-configured profiles for major DOE HPC systems:

* ``perlmutter_erf.profile`` - NERSC Perlmutter (NVIDIA A100)
* ``frontier_erf.profile`` - OLCF Frontier (AMD MI250X)
* ``aurora_erf.profile`` - ALCF Aurora (Intel GPUs)
* ``polaris_erf.profile`` - ALCF Polaris (NVIDIA A100)
* ``kestrel_erf.profile`` - NREL Kestrel (NVIDIA H100)

**Customizing for New Systems**

To build on an unsupported HPC system, copy an existing profile for a similar architecture and modify the module load and export commands to match the target system's software environment.

.. dropdown:: Example Profile: Perlmutter
   :icon: code

   .. literalinclude:: ../../Build/machines/perlmutter_erf.profile
      :language: bash

The Cray Detection System (CMake)
----------------------------------

The Cray Detection System is a two-phase automation mechanism within ERF's CMake build system. It identifies Cray Programming Environments on HPC platforms and applies necessary configurations and workarounds, shielding users from platform-specific complexity.

**What It Automates**

* Sets Cray compiler wrappers (``cc``, ``CC``, ``ftn``) as default compilers
* Inspects ``$CRAY_ACCEL_TARGET`` to determine GPU architecture (e.g., NVIDIA A100, AMD MI250X)
* Automatically sets correct GPU architecture flags for AMReX (``AMReX_CUDA_ARCH``, ``AMReX_AMD_ARCH``) and Kokkos (``Kokkos_ARCH_*``)
* Configures support for GPU-aware MPI
* Finds parallel versions of NetCDF and HDF5

**Controlling Auto-Detection**

Cray auto-detection is enabled by default when a Cray environment is identified. To disable for manual configuration:

.. code-block:: bash

   cmake -DERF_ENABLE_CRAY_AUTO_FIXES=OFF ..

.. dropdown:: Technical Details: Two-Phase Process
   :icon: info

   The system works around a CMake constraint: the toolchain must be validated before ``project()`` is invoked, but complex logic can only run after ``project()``.

   **Phase 1: Pre-Project Detection**

   Handled by ``CrayCompilerDetection.cmake``. Executes before ``project()`` in the main ``CMakeLists.txt``:

   * Checks environment variables (``$CRAYPE_VERSION``, ``$CRAY_MPICH_DIR``)
   * If detected, sets ``CMAKE_C_COMPILER`` and ``CMAKE_CXX_COMPILER`` to ``cc`` and ``CC``
   * Ensures CMake's compiler identification succeeds

   **Phase 2: Post-Project Configuration**

   Handled by ``CrayDetection.cmake``. After ``project()`` has run:

   * Detects specific GPU architecture from ``$CRAY_ACCEL_TARGET``
   * Configures dependencies for the Cray ecosystem
   * Applies automated workarounds for known build issues

.. dropdown:: Automated Workarounds
   :icon: tools

   The Cray detection system automatically resolves common configuration challenges:

   .. list-table::
      :header-rows: 1
      :widths: 25 25 50

      * - Problem
        - Detection Clue
        - Automated Solution
      * - GPU Architecture flags
        - ``$CRAY_ACCEL_TARGET`` set (e.g., ``nvidia80``)
        - Maps to correct flags (e.g., ``AMReX_CUDA_ARCH=8.0``, ``Kokkos_ARCH_AMPERE80=ON``)
      * - CUDA+EKAT MPI paths
        - ``ERF_ENABLE_CUDA=ON`` + EKAT physics enabled
        - Executes ``CC --cray-print-opts=cflags`` and prepends to ``CMAKE_CUDA_FLAGS``
      * - GPU-Aware MPI linking
        - ``$MPICH_GPU_SUPPORT_ENABLED=1``
        - Links appropriate GTL library (``-lmpi_gtl_cuda`` or ``-lmpi_gtl_hsa``)
      * - Parallel NetCDF/HDF5
        - ``ERF_ENABLE_NETCDF=ON``
        - Retrieves Cray search path via ``CC --cray-print-opts=PKG_CONFIG_PATH``

.. dropdown:: Manual Configuration
   :icon: file

   For complete control, the auto-detection system provides transparency:

   1. After configuration, detected settings are saved to ``cray_detected_config.cmake`` in the build directory
   2. View contents with ``make show-cray-config``
   3. Copy and modify as a template for custom configuration
   4. Apply with ``cmake -C path/to/my_config.cmake ..``

Build Scripts Reference
-----------------------

ERF provides tested build scripts for common configurations. The following table shows which scripts have been verified on each system.

The test procedure is documented in :download:`notes_test.sh <figures/notes_test.sh>`.

.. list-table:: Verified Build Scripts by System
   :header-rows: 1
   :widths: 35 10 10 10 10 10 10 10

   * - Build Script
     - Perlmutter
     - Frontier
     - Aurora
     - Polaris
     - Kestrel
     - RegtestCPU
     - RegtestGPU
   * - ``cmake.sh``
     - `22e12035 <https://github.com/erf-model/ERF/commit/22e12035>`_ (ABL)
     - `22e12035 <https://github.com/erf-model/ERF/commit/22e12035>`_ (ABL)
     - `22e12035 <https://github.com/erf-model/ERF/commit/22e12035>`_
     - `f8665c28 <https://github.com/jmsexton03/ERF/commit/f8665c28>`_ (ABL)
     - Untested
     - `f8665c28 <https://github.com/jmsexton03/ERF/commit/f8665c28>`_ (ABL)
     - `f8665c28 <https://github.com/jmsexton03/ERF/commit/f8665c28>`_ (ABL)
   * - ``cmake_with_kokkos_many.sh``
     - Untested
     - Untested
     - Untested
     - Untested
     - Untested
     - `22e12035 <https://github.com/erf-model/ERF/commit/22e12035>`_ (ABL)
     - Untested
   * - ``cmake_with_kokkos_many_cuda.sh``
     - `f8665c28 <https://github.com/jmsexton03/ERF/commit/f8665c28>`_ (ABL)
     - —
     - —
     - Untested
     - Untested
     - —
     - `f8665c28 <https://github.com/jmsexton03/ERF/commit/f8665c28>`_ (ABL)
   * - ``cmake_with_kokkos_many_noradiation_hip.sh``
     - —
     - `22e12035 <https://github.com/erf-model/ERF/commit/22e12035>`_ (compiled)
     - —
     - —
     - —
     - —
     - —
   * - ``cmake_with_kokkos_many_hip.sh``
     - —
     - `22e12035 <https://github.com/erf-model/ERF/commit/22e12035>`_
     - —
     - —
     - —
     - —
     - —
   * - ``cmake_with_kokkos_many_sycl.sh``
     - —
     - —
     - `22e12035 <https://github.com/erf-model/ERF/commit/22e12035>`_
     - —
     - —
     - —
     - —
   * - ``build_erf_with_shoc.sh``
     - `22e12035 <https://github.com/erf-model/ERF/commit/22e12035>`_
     - `22e12035 <https://github.com/erf-model/ERF/commit/22e12035>`_
     - `22e12035 <https://github.com/erf-model/ERF/commit/22e12035>`_
     - Untested
     - Untested
     - Tested
     - Tested
   * - ``build_erf_with_shoc_cuda.sh``
     - `22e12035 <https://github.com/erf-model/ERF/commit/22e12035>`_
     - —
     - —
     - Untested
     - Untested
     - —
     - Untested
   * - ``build_erf_with_shoc_hip.sh``
     - —
     - `22e12035 <https://github.com/erf-model/ERF/commit/22e12035>`_
     - —
     - Untested
     - Untested
     - —
     - Untested
   * - ``build_erf_with_shoc_sycl.sh``
     - —
     - —
     - `22e12035 <https://github.com/erf-model/ERF/commit/22e12035>`_
     - Untested
     - Untested
     - —
     - Untested
   * - ``Perlmutter/build_erf_with_shoc_cuda_Perlmutter.sh``
     - Untested
     - Untested
     - Untested
     - Untested
     - Untested
     - —
     - Untested

.. note::
   **Reading the table:**

   * **Commit hashes** link to GitHub and indicate last successful test
   * **(ABL)** indicates tested with Atmospheric Boundary Layer case and validated with ``fcompare``
   * **(compiled)** indicates build succeeded but runtime results need investigation
   * **Tested** means successful build/run without specific commit hash
   * **Untested** means not verified but may work with modifications
   * **—** indicates configuration not applicable to that system

**Script Descriptions**

* ``cmake.sh`` - Basic CPU-only build
* ``cmake_with_kokkos_many.sh`` - Kokkos-enabled CPU build
* ``cmake_with_kokkos_many_cuda.sh`` - NVIDIA GPUs with CUDA
* ``cmake_with_kokkos_many_noradiation_hip.sh`` - AMD GPUs (HIP) without radiation
* ``cmake_with_kokkos_many_hip.sh`` - AMD GPUs (HIP) full build
* ``cmake_with_kokkos_many_sycl.sh`` - Intel GPUs with SYCL
* ``build_erf_with_shoc.sh`` - Automated SHOC workflow (CPU)
* ``build_erf_with_shoc_{cuda,hip,sycl}.sh`` - Automated SHOC workflow with GPU backend
* ``Perlmutter/build_erf_with_shoc_cuda_Perlmutter.sh`` - SHOC with CUDA on Perlmutter

.. note::
   The GPU SHOC scripts can be auto-generated. Set ``BACKEND=CUDA`` (or ``HIP``/``SYCL``), then:

   .. code-block:: bash

      sed "/ERF_ENABLE_MPI/a\      -DERF_ENABLE_${BACKEND}:BOOL=ON \\\\" Build/cmake_with_shoc.sh > Build/cmake_with_shoc_${BACKEND,,}.sh
      sed "s/cmake_with_shoc.sh/cmake_with_shoc_${BACKEND,,}.sh/" Build/build_erf_with_shoc.sh > Build/build_erf_with_shoc_${BACKEND,,}.sh
      chmod +x Build/cmake_with_shoc_${BACKEND,,}.sh Build/build_erf_with_shoc_${BACKEND,,}.sh

.. dropdown:: Build Script Examples
   :icon: code

   **Basic CPU Build:**

   .. literalinclude:: ../../Build/cmake.sh
      :language: bash

   **CUDA + Kokkos:**

   .. literalinclude:: ../../Build/cmake_with_kokkos_many_cuda.sh
      :language: bash

   **HIP + Kokkos:**

   .. literalinclude:: ../../Build/cmake_with_kokkos_many_hip.sh
      :language: bash

   **SYCL + Kokkos:**

   .. literalinclude:: ../../Build/cmake_with_kokkos_many_sycl.sh
      :language: bash

   **SHOC (CPU):**

   .. literalinclude:: ../../Build/cmake_with_shoc.sh
      :language: bash

   Requires E3SM submodules: ``source Build/GNU_Ekat/eamxx_clone.sh``

   **Automated SHOC Workflow:**

   .. literalinclude:: ../../Build/build_erf_with_shoc.sh
      :language: bash

   Usage: ``source Build/build_erf_with_shoc.sh``

For detailed CMake options and dependencies, see :ref:`sec:build:systems`.

.. dropdown:: GNU Make vs CMake on HPC
   :icon: info

   Both build systems are fully supported. The choice depends on desired abstraction level and control.

   .. list-table::
      :header-rows: 1
      :widths: 20 80

      * - Build System
        - Strengths & Best Use Cases
      * - CMake
        - **Strengths:** Automates complex dependency graphs (ERF → EKAT → Kokkos). Cray Detection System handles configurations. Generates ``find_package``-compatible configuration (``ERFConfig.cmake``). **Best For:** Most users on supported HPC systems (Perlmutter, Frontier, Polaris, Aurora). Developers building applications that depend on ERF.
      * - GNU Make
        - **Strengths:** Uses existing AMReX infrastructure. Utility targets (e.g., ``print-xxx``) for debugging. Direct control over every build variable. **Best For:** Advanced users, performance tuning experts, systems not yet supported by CMake auto-detection.

Workstation Builds
------------------

Building on a local workstation is simpler than on HPC systems, as it doesn't require managing environment modules or vendor-specific toolchains. This is essential for development, debugging, and running smaller test cases.

The RegtestCPU and RegtestGPU columns in the build scripts table indicate configurations tested on workstation environments based on the nightly regression test setups.

**macOS Builds**

Building on macOS is supported with specific considerations:

* GNU Make is often recommended for simplicity and direct control
* Use ``Make.local`` (``amrex/Tools/GNUMake/Make.local``) to specify consistent compiler suite (e.g., GCC from Homebrew) for both C++ and Fortran
* This avoids conflicts between default Clang/Xcode C++ compiler and separately installed Fortran compiler

Example ``Make.local`` for macOS:

.. code-block:: bash

   CXX = g++-13
   CC  = gcc-13
   FC  = gfortran-13

Running on HPC Systems
----------------------

After building ERF, job submission requires system-specific batch scripts.

**General HPC Resources**

For comprehensive job submission examples across DOE systems, see the `WarpX HPC Documentation <https://warpx.readthedocs.io/en/latest/install/hpc.html>`__. ERF follows similar patterns to WarpX as both are built on AMReX.

Perlmutter (NERSC)
~~~~~~~~~~~~~~~~~~

.. tab-set::

   .. tab-item:: Building with GNU Make

      Simple build using GNU compiler and CUDA:

      .. code-block:: bash

         # Load environment
         module load PrgEnv-gnu cudatoolkit cray-mpich cray-netcdf-hdf5parallel

         # Navigate to problem directory
         cd ${ERF_HOME}/Exec/ABL

         # Build
         make -j4 COMP=gnu USE_MPI=TRUE USE_CUDA=TRUE

      This produces an executable like ``ERF3d.gnu.MPI.CUDA.ex``.

   .. tab-item:: Building with CMake

      Using the provided build script:

      .. code-block:: bash

         # Load environment
         source $ERF_HOME/Build/machines/perlmutter_erf.profile

         # Configure and build (out-of-source)
         mkdir build && cd build
         ../Build/cmake_with_kokkos_many_cuda.sh

      **Executable location:** ``build/Exec/ABL/erf_abl`` (or ``install/bin/erf_abl`` if installed)

      Or manual configuration:

      .. code-block:: bash

         cmake -DCMAKE_BUILD_TYPE=Release \
               -DERF_ENABLE_MPI=ON \
               -DERF_ENABLE_CUDA=ON \
               -DERF_ENABLE_NETCDF=ON \
               -DERF_ENABLE_RRTMGP=ON \
               ..
         make -j4

   .. tab-item:: Basic GPU Job

      This example runs ERF on 4 nodes with GPU-aware MPI enabled.

      **Before submitting**, load the environment:

      .. code-block:: bash

         source $ERF_HOME/Build/machines/perlmutter_erf.profile

         # Run from scratch filesystem with executable and inputs in same directory
         cd $PSCRATCH/ERF/rundir

         # Verify paths before launching
         ls -lh ./ERF3d.*.ex inputs

      **Job submission script:**

      .. code-block:: bash

         #!/bin/bash
         #SBATCH --account=m4106_g
         #SBATCH --nodes=4
         #SBATCH --ntasks-per-node=4
         #SBATCH --gpus-per-node=4
         #SBATCH --gpu-bind=none
         #SBATCH --time=00:30:00
         #SBATCH --constraint=gpu&hbm40g
         #SBATCH --job-name=ERF
         #SBATCH --output=erf_%j.out

         # GPU-aware MPI optimizations
         export MPICH_OFI_NIC_POLICY=GPU
         export MPICH_GPU_SUPPORT_ENABLED=1
         export SLURM_CPU_BIND="cores"

         # Launch with CUDA device ordering
         srun -n 16 --cpus-per-task=4 --cpu-bind=cores bash -c "
           export CUDA_VISIBLE_DEVICES=\$((3-SLURM_LOCALID));
           ./ERF3d.gnu.MPI.CUDA.ex inputs amrex.use_gpu_aware_mpi=1"

      Submit with: ``sbatch job_script.sh``

   .. tab-item:: 80GB GPU Nodes

      For the 256 nodes with 80GB HBM per GPU, replace:

      .. code-block:: bash

         #SBATCH --constraint=gpu&hbm40g

      with:

      .. code-block:: bash

         #SBATCH --constraint=gpu&hbm80g

.. dropdown:: Advanced: AMReX Scaling Tests
   :icon: beaker

   Production script demonstrating NIC pinning optimization from the `AMReX FFT scaling repository <https://github.com/WeiqunZhang/amrex-scaling/tree/main/fft/perlmutter>`__:

   .. remote-include:: https://raw.githubusercontent.com/WeiqunZhang/amrex-scaling/refs/heads/main/fft/perlmutter/2025-02-06/run-4.sh
      :language: bash

   Key features:

   * Uses ``--gpus-per-task=1`` for fine-grained GPU binding
   * Compares default NIC policy vs ``MPICH_OFI_NIC_POLICY=GPU``
   * Demonstrates multiple runs with different configurations

.. dropdown:: WarpX Production Script
   :icon: code

   Reference implementation from the `WarpX project <https://github.com/BLAST-WarpX/warpx>`__:

   .. remote-include:: https://raw.githubusercontent.com/BLAST-WarpX/warpx/refs/heads/development/Tools/machines/perlmutter-nersc/perlmutter_gpu.sbatch
      :language: bash

Kestrel (NREL)
~~~~~~~~~~~~~~

The `Kestrel <https://nrel.github.io/HPC/Documentation/Systems/Kestrel/>`__ cluster is an HPE Cray system with Intel Xeon Sapphire Rapids CPU nodes (104 cores) and NVIDIA H100 GPUs (4 per node).

.. note::
   Kestrel has **separate login nodes for GPU work**. Access GPU login nodes via ``kestrel-gpu.hpc.nrel.gov``. GPU jobs should only be submitted from GPU login nodes.

**Building**

.. tab-set::

   .. tab-item:: GNU Make (CPU)

      .. code-block:: bash

         # Reset to default environment
         module restore

         # Build with Cray compilers
         cd ${ERF_HOME}/Exec/ABL
         make realclean
         make -j COMP=cray

   .. tab-item:: GNU Make (GPU)

      .. code-block:: bash

         # Load GPU modules
         module purge
         module load PrgEnv-gnu/8.5.0
         module load cuda/12.3
         module load craype-x86-milan

         # Build
         cd ${ERF_HOME}/Exec/ABL
         make realclean
         make -j COMP=gnu USE_CUDA=TRUE

   .. tab-item:: CMake (GPU)

      .. code-block:: bash

         # Load GPU modules
         module purge
         module load PrgEnv-gnu/8.5.0
         module load cuda/12.3
         module load craype-x86-milan

         # Configure and build
         mkdir build && cd build
         cmake -DCMAKE_BUILD_TYPE=Release \
               -DERF_ENABLE_MPI=ON \
               -DERF_ENABLE_CUDA=ON \
               ..
         make -j

.. warning::
   **System updates on Kestrel periodically change required modules.** Verify current module names with ``module avail`` before building.

**Memory Allocation**

Kestrel allows partial node allocations. For memory-intensive operations (e.g., CUDA compilation):

.. code-block:: bash

   # Option 1: Request exclusive node access
   #SBATCH --exclusive

   # Option 2: Request specific memory
   #SBATCH --mem=240G
   # or
   #SBATCH --mem-per-cpu=2G

Without these flags, CUDA compilation may fail due to insufficient memory.

**Performance and Cost Considerations**

GPU node hours on Kestrel are charged at **10× the rate** of CPU node hours. Understanding performance trade-offs is essential for efficient use of allocations.

**Typical Performance Characteristics:**

* GPU nodes (4× H100): **10-20× faster** than CPU nodes (96-104 cores)
* Best efficiency with **>1M cells per GPU**
* Smaller problems may not fully utilize GPU capability

**When to Use GPU vs CPU:**

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Problem Size
     - GPU Nodes (10× cost)
     - CPU Nodes (1× cost)
   * - < 500K cells/GPU
     - May not justify 10× cost
     - **Recommended**
   * - 500K - 1M cells/GPU
     - Marginal benefit
     - Consider for development
   * - > 1M cells/GPU
     - **Recommended** - cost effective
     - Slower time to solution
   * - > 5M cells/GPU
     - **Excellent utilization**
     - May exceed wall-time limits

**Recommendations:**

1. **Profile your specific case** - Performance varies with physics packages and I/O frequency
2. **Development on CPU** - Use CPU nodes for code development and small test cases
3. **Production on GPU** - Use GPU nodes for production runs with well-sized domains
4. **Monitor utilization** - Check GPU usage with ``nvidia-smi`` to verify saturation

.. note::
   The 10-20× performance gain typically justifies the 10× cost increase for production runs, providing faster time-to-solution and 1-2× better overall cost efficiency measured in allocation units per simulation.
