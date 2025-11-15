.. _sec:build:hpc:

=================================
ERF on HPC Systems: A Guide to the Cray Detection Framework
=================================

1.0 Introduction and Environment Setup
--------------------------------------

Building ERF on High-Performance Computing (HPC) systems requires managing diverse architectures, proprietary compiler toolchains, and specialized environment modules. The ERF build system automates this complexity, providing a reproducible path from source code to executable.

A mandatory first step before any compilation attempt on an HPC platform is the correct preparation of the shell environment. Loading the appropriate modules and setting environment variables ensures that the build system can locate the correct compiler wrappers, essential parallel libraries like MPI, hardware-specific toolchains such as the CUDA or ROCm toolkits, and I/O libraries like NetCDF and HDF5. Without a properly configured environment, builds will fail with errors related to missing dependencies.

To standardize this setup process, ERF provides a set of machine profile files. These are pre-configured scripts that encapsulate the system-specific knowledge required to establish a correct build environment on a given supercomputer. By using these profiles, developers and users can reliably configure their shell sessions, ensuring a consistent and reproducible foundation for compilation.

This document explains the architecture of ERF's HPC build strategy, with a particular focus on the automated Cray detection system integrated within its CMake build framework. It provides a high-level understanding of the automation logic and design philosophy rather than serving as a step-by-step tutorial for building the code.

2.0 Machine Profile Files
-------------------------

Machine profile files form the foundational layer of the ERF build process on HPC systems. They abstract the details of diverse supercomputer ecosystems by managing environment modules and variables in a standardized, reproducible manner. A correctly configured environment is the prerequisite for any successful HPC build, and these profiles are the primary tool for achieving that configuration.

.. note::
   **HPC system environments change regularly.** Module names, versions, and
   availability vary as systems are updated. Always verify current module names
   with ``module avail`` before building. The profile scripts in
   ``Build/machines/`` are maintained for major systems but may require
   adjustment after system upgrades.

**2.1 Purpose and Usage**

Machine profile files are shell scripts located in the ``Build/machines/`` directory. Their purpose is to prepare the shell environment for compilation by loading required software modules, exporting necessary environment variables, and defining paths to dependencies.

The primary functions of these profiles include:

* Executing module load commands. The core of the environment on Cray systems is the ``PrgEnv-*`` module (e.g., ``PrgEnv-gnu`` or ``PrgEnv-cray``), which dictates the underlying compiler suite that the Cray wrappers (``cc``, ``CC``) will use. This is supplemented by modules for hardware acceleration (e.g., ``craype-accel-nvidia80``), parallel libraries (e.g., ``cray-netcdf-hdf5parallel``), and development tools (e.g., ``cmake``).
* Running export commands to set environment variables that guide the build system, such as hints for locating libraries or specifying target hardware.

To use a profile, source the appropriate script from your shell. This executes the commands within the current shell session, modifying its environment. Once the profile has been sourced, the shell is ready for the subsequent cmake configuration and build commands. For example:

.. code:: shell

   source machines/perlmutter_erf.profile

**2.2 Available Profiles**

ERF provides pre-configured machine profiles for several major DOE HPC systems:

* ``perlmutter_erf.profile`` (NERSC Perlmutter)
* ``frontier_erf.profile`` (OLCF Frontier)
* ``aurora_erf.profile`` (ALCF Aurora)
* ``polaris_erf.profile`` (ALCF Polaris)
* ``kestrel_erf.profile`` (NERSC Kestrel)

**2.3 Customizing for a New System**

To build ERF on an unsupported HPC system, create a custom profile. The recommended approach is to copy an existing profile for a machine with a similar architecture (e.g., another system using NVIDIA GPUs or one based on AMD hardware) and modify the module load and export commands to match the target system's specific software environment and requirements.

While manual profile files are essential for environment preparation, they address only the state of the shell. They cannot manage the in-build configuration complexities that arise once CMake begins its process. Sourcing a profile sets environment variables, but CMake operates in its own configuration space and requires explicit information about compiler flags, library paths, and dependency relationships that module load commands do not directly provide. To bridge this gap, the integrated Cray detection system provides the next layer of automation.

3.0 The Cray Detection System (CMake)
-------------------------------------

The Cray Detection System is a two-phase automation mechanism within ERF's CMake build system. It programmatically identifies Cray Programming Environments on HPC platforms like Perlmutter, Frontier, and Polaris, and then applies necessary configurations and workarounds. The primary goal is to shield users from platform-specific complexity, delivering a consistent, reliable, and reproducible build experience with minimal manual intervention.

The system automates several configuration steps:

* Sets the Cray compiler wrappers (``cc``, ``CC``, ``ftn``) as the default compilers for C, C++, and Fortran
* Inspects the ``$CRAY_ACCEL_TARGET`` environment variable to determine the specific GPU architecture (e.g., NVIDIA A100, AMD MI250X)
* Automatically sets the correct GPU architecture flags for key dependencies, including AMReX (``AMReX_CUDA_ARCH``, ``AMReX_AMD_ARCH``) and Kokkos-based physics packages (``Kokkos_ARCH_*``)
* Configures support for GPU-aware MPI and finds the correct parallel versions of libraries like NetCDF and HDF5

**3.1 Overview of the Two-Phase Process**

This two-phase design works around a fundamental constraint in CMake's architecture: the toolchain must be validated before the ``project()`` command is invoked, but complex, project-specific logic can only run after ``project()``. The system splits detection accordingly:

* **Phase 1: Pre-Project Detection** - Handled by ``CrayCompilerDetection.cmake`` and executes before the ``project()`` call in the main ``CMakeLists.txt`` file. It performs a lightweight check of environment variables like ``$CRAYPE_VERSION`` and ``$CRAY_MPICH_DIR`` to determine if running in a Cray environment. If detected, it immediately sets the ``CMAKE_C_COMPILER`` and ``CMAKE_CXX_COMPILER`` cache variables to ``cc`` and ``CC``, ensuring CMake's compiler identification succeeds.

* **Phase 2: Post-Project Configuration** - After the ``project()`` command has run and compilers are known, ``CrayDetection.cmake`` performs deeper inspection. It detects the specific GPU architecture from ``$CRAY_ACCEL_TARGET``, configures dependencies for the Cray ecosystem, and applies automated workarounds for known build issues.

**3.2 Controlling Auto-Detection**

The Cray auto-detection system is enabled by default when a Cray environment is identified. For advanced users, debugging, or situations requiring manual configuration, this automation can be disabled by setting: ``-DERF_DISABLE_CRAY_AUTO_FIXES=ON``

**3.3 Key Detection Steps and Automated Fixes**

The Cray detection system automatically resolves several common configuration challenges on Cray supercomputers. These automated workarounds apply during Post-Project Configuration (Phase 2), after compilers have been identified and validated during Phase 1.

.. list-table:: Automated Cray Workarounds
   :header-rows: 1
   :widths: 25 25 50

   * - Problem
     - Detection Clue
     - Automated Solution
   * - GPU Architecture: Manually specifying correct GPU architecture flags for multiple libraries is error-prone
     - ``$CRAY_ACCEL_TARGET`` environment variable is set (e.g., to ``nvidia80`` or ``amd_gfx90a``)
     - System maps ``$CRAY_ACCEL_TARGET`` to correct architecture flags. For example, ``nvidia80`` sets ``AMReX_CUDA_ARCH`` to 8.0 and enables ``Kokkos_ARCH_AMPERE80``
   * - CUDA+EKAT Builds: CUDA compiler (nvcc) not aware of MPI include paths in Cray environment
     - ``ERF_ENABLE_CUDA=ON`` and ``ERF_ENABLE_RRTMGP``, ``ERF_ENABLE_SHOC``, or ``ERF_ENABLE_P3`` enabled
     - System executes ``CC --cray-print-opts=cflags`` to retrieve missing compiler flags and prepends them to ``CMAKE_CUDA_FLAGS``
   * - GPU-Aware MPI: Correctly linking MPI GPU Transport Layer (GTL) libraries is complex and system-dependent
     - ``$MPICH_GPU_SUPPORT_ENABLED=1`` environment variable is set
     - System identifies base MPI library and automatically links appropriate GTL library (``-lmpi_gtl_cuda`` for NVIDIA, ``-lmpi_gtl_hsa`` for AMD)
   * - Parallel NetCDF/HDF5: Standard CMake ``find_package`` calls fail to locate parallel versions on Cray systems
     - ``ERF_ENABLE_NETCDF=ON``
     - System retrieves Cray-specific search path via ``CC --cray-print-opts=PKG_CONFIG_PATH`` and prepends to ``$PKG_CONFIG_PATH`` environment variable

**3.4 Manual Configuration via -C Files**

For advanced users requiring complete control, the auto-detection system provides transparency through configuration files:

1. After successful configuration, auto-detection logic saves detected settings to ``cray_detected_config.cmake`` in the build directory
2. Use ``make show-cray-config`` to print contents to console without locating file manually
3. Copy and modify this file as template for custom manual configuration, then apply using ``cmake -C path/to/my_config.cmake ..``

4.0 Build Scripts and System Comparisons
----------------------------------------

ERF provides tested build scripts that encapsulate common configurations for various systems and architectures. The choice between the modern, automated CMake system and the traditional GNU Make system depends on specific goals, familiarity with tools, and level of control required.

**4.1 Script Reference Table**

The following table summarizes build scripts provided with ERF and indicates which combinations have been tested and verified on each major platform. Verified builds are marked with the git commit hash at which they were last successfully tested.

.. list-table:: ERF Build Scripts by System
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
     - f8665c28 (ABL)
     - f8665c28 (ABL)
     - f8665c28
     - f8665c28 (ABL)
     - Untested
     - f8665c28 (ABL)
     - f8665c28 (ABL)
   * - ``cmake_with_kokkos_many.sh``
     - Untested
     - Untested
     - Untested
     - Untested
     - Untested
     - f8665c28 (ABL)
     - Untested
   * - ``cmake_with_kokkos_many_cuda.sh``
     - f8665c28 (ABL)
     - —
     - —
     - Untested
     - Untested
     - —
     - f8665c28 (ABL)
   * - ``cmake_with_kokkos_many_noradiation_hip.sh``
     - —
     - f8665c28 (compiled)
     - —
     - —
     - —
     - —
     - —
   * - ``cmake_with_kokkos_many_hip.sh``
     - —
     - f8665c28
     - —
     - —
     - —
     - —
     - —
   * - ``cmake_with_kokkos_many_sycl.sh``
     - —
     - —
     - f8665c28
     - —
     - —
     - —
     - —
   * - ``build_erf_with_shoc.sh``
     - Tested
     - Untested
     - Untested
     - Untested
     - Untested
     - Tested
     - Tested
   * - ``Perlmutter/build_erf_with_shoc_cuda_Perlmutter.sh``
     - Untested
     - Untested
     - Untested
     - Untested
     - Untested
     - —
     - Untested

.. note::
   **Reading the verification table:**

   * **Commit hashes** (e.g., ``f8665c28``) indicate the last git commit where the script was successfully tested
   * **(ABL)** indicates the build was tested by running the Atmospheric Boundary Layer test case and validated with AMReX's ``fcompare`` utility to verify correctness
   * **(compiled)** indicates the build succeeded but runtime results were either not collected or showed unexpected behavior requiring further investigation
   * **Tested** means the script built and ran successfully but lacks a specific commit hash record
   * **Untested** means the script has not been verified on that system but may work with appropriate modifications
   * **—** (em dash) indicates the configuration is not applicable to that system

   Users who successfully test unlisted configurations are encouraged to report their results to help expand this table.

**4.2 Script Descriptions**

* ``cmake.sh``: Basic CMake build script for CPU-only configurations that works on most systems
* ``cmake_with_kokkos_many.sh``: Kokkos-enabled build intended for CPU architectures
* ``cmake_with_kokkos_many_cuda.sh``: Kokkos-enabled build optimized for NVIDIA GPUs using the CUDA backend
* ``cmake_with_kokkos_many_noradiation_hip.sh``: Kokkos-enabled build for AMD GPUs (HIP backend) with radiation model disabled
* ``cmake_with_kokkos_many_hip.sh``: Full Kokkos-enabled build for AMD GPUs using the HIP backend
* ``cmake_with_kokkos_many_sycl.sh``: Kokkos-enabled build for Intel GPUs using the SYCL backend
* ``Perlmutter/build_erf_with_shoc_cuda_Perlmutter.sh``: Specialized script for building SHOC turbulence model with CUDA support on Perlmutter

**4.3 GNU Make vs. CMake on HPC**

Both GNU Make and CMake are fully supported build systems. The choice between them is based on the desired level of abstraction and control.

.. list-table:: Build System Comparison
   :header-rows: 1
   :widths: 20 80

   * - Build System
     - Strengths & Best Use Cases
   * - CMake
     - **Strengths:** Superbuilds and handling complex dependency graphs (e.g., ERF → EKAT → Kokkos) are easier. Automates configurations via Cray Detection System. Generates and installs ``find_package``-compatible configuration (``ERFConfig.cmake``), allowing other CMake projects to consume ERF as dependency. **Best For:** Most users on supported HPC systems (Perlmutter, Frontier, Polaris, Aurora). Developers building external applications that depend on ERF.
   * - GNU Make
     - **Strengths:** Leverages existing AMReX infrastructure and has utility targets (e.g., ``print-xxx``) for debugging. Offers transparent and direct control over every build variable by manually setting them in ``GNUmakefile``. **Best For:** Advanced users, performance tuning experts, and those building on HPC systems not yet fully supported by CMake auto-detection logic.

**4.4 GNU Make Configuration Example**

Configuration with GNU Make is handled by setting variables directly in a ``GNUmakefile`` or on the command line:

.. code:: shell

   make COMP=cray USE_MPI=TRUE USE_CUDA=TRUE

This approach provides direct control over every aspect of the build, leveraging the AMReX make system.

**4.5 Build Script Examples**

The following dropdowns show complete build scripts demonstrating different feature combinations:

.. dropdown:: cmake.sh - Basic build
   :color: light
   :icon: code
   :animate: fade-in-slide-down

   Basic configuration.

   .. literalinclude:: ../../../Build/cmake.sh
      :language: bash

.. dropdown:: cmake_with_kokkos_many_cuda.sh - CUDA + Kokkos superbuild
   :color: light
   :icon: code
   :animate: fade-in-slide-down

   Full GPU build for NVIDIA hardware with FFT, NETCDF, HDF5, and RRTMGP packages.

   .. literalinclude:: ../../../Build/cmake_with_kokkos_many_cuda.sh
      :language: bash

.. dropdown:: cmake_with_kokkos_many_hip.sh - HIP + Kokkos
   :color: light
   :icon: code
   :animate: fade-in-slide-down

   Full GPU build for AMD hardware with FFT, NETCDF, HDF5, and RRTMGP packages.

   .. literalinclude:: ../../../Build/cmake_with_kokkos_many_hip.sh
      :language: bash

.. dropdown:: cmake_with_kokkos_many_sycl.sh - SYCL + Kokkos
   :color: light
   :icon: code
   :animate: fade-in-slide-down

   Full GPU build for Intel hardware with FFT, NETCDF, HDF5, and RRTMGP packages.

   .. literalinclude:: ../../../Build/cmake_with_kokkos_many_sycl.sh
      :language: bash

.. dropdown:: cmake_with_shoc.sh - CPU build with SHOC
   :color: light
   :icon: code
   :animate: fade-in-slide-down

   CPU-only build with SHOC turbulence model.

   .. note::
      Requires E3SM submodules: ``source Build/GNU_Ekat/eamxx_clone.sh``

   .. literalinclude:: ../../../Build/cmake_with_shoc.sh
      :language: bash

.. dropdown:: build_erf_with_shoc.sh - Automated SHOC workflow
   :color: light
   :icon: code
   :animate: fade-in-slide-down

   Automated script that clones E3SM dependencies and invokes CMake configuration.

   Usage: ``source Build/build_erf_with_shoc.sh``

   .. literalinclude:: ../../../Build/build_erf_with_shoc.sh
      :language: bash

See Section 2.4 in :ref:`buildingSystems` for detailed explanations of CMake flags and their interdependencies.

5.0 Workstation Builds
----------------------

Building ERF on a local workstation is essential for development, debugging, and running smaller test cases. The process is simpler than on HPC systems, as it does not require managing environment modules or navigating vendor-specific toolchains.

The RegtestCPU and RegtestGPU columns in the script reference table indicate configurations that are regularly tested and verified to work in common workstation environments.

**5.1 macOS-Specific Considerations**

Building on macOS is supported, though it comes with specific considerations. The following practices are recommended:

* GNU Make is often recommended for its simplicity and direct control
* A ``Make.local`` file can be used to specify a consistent compiler suite (e.g., GCC installed via Homebrew) for both C++ and Fortran
* This approach helps avoid potential conflicts between the default Clang/Xcode C++ compiler and a separately installed Fortran compiler

6.0 Running ERF on HPC Systems
------------------------------

After building ERF, job submission requires system-specific batch scripts. This section provides tested examples and references to production workflows.

**6.1 General HPC Resources**

For comprehensive job submission examples across multiple DOE systems, see the `WarpX HPC Documentation <https://warpx.readthedocs.io/en/latest/install/hpc.html>`__, which provides detailed scripts for Perlmutter, Frontier, Polaris, and other systems. ERF follows similar patterns to WarpX as both are built on AMReX.

**6.2 Perlmutter (NERSC) - GPU Nodes**

.. tab-set::

   .. tab-item:: Basic GPU Job

      This example runs ERF on 4 nodes with GPU-aware MPI enabled. **Before submitting**, ensure you have sourced the machine profile or loaded required modules:

      .. code-block:: bash

         # Load environment
         source $ERF_HOME/Build/machines/perlmutter_erf.profile

         # Or manually load modules:
         # module load PrgEnv-gnu cudatoolkit cray-mpich cray-netcdf-hdf5parallel

      Job submission script:

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

   .. tab-item:: Advanced: AMReX Scaling Tests

      Production script demonstrating NIC pinning optimization from the `AMReX FFT scaling repository <https://github.com/WeiqunZhang/amrex-scaling/tree/main/fft/perlmutter>`__:

      .. dropdown:: run-4.sh (fetched from amrex-scaling repo)
         :color: light
         :icon: code
         :animate: fade-in-slide-down

         .. remote-include:: https://raw.githubusercontent.com/WeiqunZhang/amrex-scaling/refs/heads/main/fft/perlmutter/2025-02-06/run-4.sh
            :language: bash

      Key features:

      * Uses ``--gpus-per-task=1`` for fine-grained GPU binding
      * Compares default NIC policy vs ``MPICH_OFI_NIC_POLICY=GPU``
      * Demonstrates multiple runs with different configurations

   .. tab-item:: WarpX Production Script

      Reference implementation from the `WarpX project <https://github.com/BLAST-WarpX/warpx>`__:

      .. dropdown:: perlmutter_gpu.sbatch (fetched from WarpX repo)
         :color: light
         :icon: code
         :animate: fade-in-slide-down

         .. remote-include:: https://raw.githubusercontent.com/BLAST-WarpX/warpx/refs/heads/development/Tools/machines/perlmutter-nersc/perlmutter_gpu.sbatch
            :language: bash

   .. tab-item:: 80GB GPU Nodes

      For the 256 nodes with 80GB HBM per GPU, replace ``--constraint=gpu&hbm40g`` with ``--constraint=gpu&hbm80g`` in the scripts above.

**6.3 Kestrel (NREL)**

The `Kestrel <https://nrel.github.io/HPC/Documentation/Systems/Kestrel/>`__ cluster is an HPE Cray system with Intel Xeon Sapphire Rapids CPU nodes (104 cores) and a GPU partition with NVIDIA H100 GPUs (4 per node).

.. note::
   Kestrel has **separate login nodes for GPU work**. Access GPU login nodes via ``kestrel-gpu.hpc.nrel.gov``. GPU jobs should only be submitted from GPU login nodes.

**CPU Builds**

For CPU-only builds, default modules are typically sufficient:

.. code-block:: bash

   # Reset to default environment if needed
   module restore

   # Build with Cray compilers
   make realclean
   make -j COMP=cray

**GPU Builds**

For GPU builds, load the following modules:

.. code-block:: bash

   module purge
   module load PrgEnv-gnu/8.5.0
   module load cuda/12.3
   module load craype-x86-milan

   # Build
   make realclean
   make -j COMP=gnu USE_CUDA=TRUE

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

GPU node hours on Kestrel are charged at **10× the rate** of CPU node hours. Understanding the performance trade-offs is essential for efficient use of your allocation.

**Typical Performance Characteristics**

For ERF simulations on Kestrel:

* GPU nodes (4× H100): **10-20× faster** than CPU nodes (96-104 cores)
* Best efficiency achieved with **>1M cells per GPU**
* Smaller problems may not fully utilize GPU capability

**When to Use GPU vs CPU Nodes**

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

**Recommendations**

1. **Profile your specific case** - Performance varies with physics packages and I/O frequency
2. **Development on CPU** - Use CPU nodes for code development and small test cases
3. **Production on GPU** - Use GPU nodes for production runs with well-sized domains
4. **Monitor utilization** - Check GPU usage with ``nvidia-smi`` to verify saturation

.. note::
   The 10-20× performance gain typically justifies the 10× cost increase for production runs, providing faster time-to-solution and 1-2× better overall cost efficiency measured in allocation units per simulation.
