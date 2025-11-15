.. _sec:build:library:

Library Configuration
=====================

Configuring external libraries enables ERF's scientific capabilities, from parallel I/O on supercomputing platforms to GPU-accelerated physics packages. This guide provides technical details for developers to navigate the ERF library ecosystem, ensuring a robust and reproducible path from source code to executable.

For basic build instructions, see :ref:`sec:build:systems`. For HPC-specific configurations, see :ref:`sec:build:hpc`.

Library Dependencies Overview
------------------------------

ERF integrates several external libraries to provide core functionality and advanced features:

.. list-table:: ERF Library Dependencies
   :header-rows: 1
   :widths: 20 40 15 25

   * - Library
     - Description
     - Requirement
     - Build Options
   * - **AMReX**
     - Block-structured adaptive mesh refinement (AMR), data structures, and solvers
     - Required
     - ``-DERF_USE_INTERNAL_AMREX=ON`` / ``AMREX_HOME``
   * - **NetCDF**
     - High-level I/O operations for plotfiles and initial conditions (WRF/Metgrid files)
     - Optional
     - ``-DERF_ENABLE_NETCDF=ON`` / ``USE_NETCDF=TRUE``
   * - **HDF5**
     - Low-level backend for parallel I/O, used by NetCDF for distributed-memory file access
     - Optional
     - ``-DERF_ENABLE_HDF5=ON`` / Via NetCDF pkg-config
   * - **SHOC**
     - Simplified Higher-Order Closure turbulence and cloud macrophysics from E3SM
     - Optional
     - ``-DERF_ENABLE_SHOC=ON`` / ``USE_SHOC=TRUE``
   * - **P3**
     - Predicted Particle Properties microphysics from E3SM
     - Optional
     - ``-DERF_ENABLE_P3=ON`` / ``USE_P3=TRUE``
   * - **RRTMGP**
     - Rapid Radiative Transfer Model for GCMs (radiation)
     - Optional
     - ``-DERF_ENABLE_RRTMGP=ON`` / ``USE_RRTMGP=TRUE``
   * - **ZFP**
     - High-performance lossy data compression as HDF5 filter
     - Optional
     - ``-DAMReX_HDF5_ZFP=ON``

**Parallel I/O Consistency Requirement**

For MPI-enabled ERF builds with parallel I/O (``ERF_ENABLE_NETCDF=ON``), both NetCDF and HDF5 must be parallel-enabled versions. These specialized builds contain communication hooks to coordinate I/O operations across multiple compute nodes. Linking a parallel ERF build against serial versions of these libraries results in linker errors or runtime failures during file operations.

NetCDF and HDF5 Configuration
------------------------------

NetCDF and HDF5 form the foundation of ERF's I/O strategy for large-scale, data-intensive simulations on HPC platforms. They provide mechanisms for reading complex initial conditions and writing massive datasets in a scalable manner.

**Library Roles**

* **NetCDF** - High-level API for structured, self-describing data formats. ERF uses NetCDF for:
  
  - Reading initial conditions: ``ERF_InitFromWRFInput.cpp``
  - Writing plotfiles: ``ERF_NCPlotFile.cpp``

* **HDF5** - Low-level parallel storage layer that NetCDF relies on in distributed-memory environments. Manages coordinating data writes from many MPI ranks into a single file on parallel filesystems.

**CMake Detection Hierarchy**

ERF's CMake build system uses a hierarchical search to locate NetCDF:

.. list-table:: NetCDF Detection Priority
   :header-rows: 1
   :widths: 10 25 25 40

   * - Priority
     - Method
     - Variables Checked
     - Notes
   * - 1 (Highest)
     - User-Provided CMake Options
     - ``NETCDF_DIR``
     - Passed via ``-DNETCDF_DIR=...``, overrides all other mechanisms
   * - 2 (High)
     - Environment Variables
     - ``NETCDF_DIR``
     - Sourced from shell environment via ``export NETCDF_DIR=...``
   * - 3 (Medium)
     - pkg-config Utility
     - Searches: ``netcdf``, ``netcdf-mpi``, ``netcdf_parallel``, ``netcdf-cxx4_parallel``
     - Standard detection method on Linux systems
   * - 4 (Low)
     - System Default Paths
     - ``/usr/lib``, ``/usr/include``
     - Fallback for system-wide installations

**Build Options**

.. list-table:: CMake Options for I/O Libraries
   :header-rows: 1
   :widths: 30 50 20

   * - CMake Option
     - Description
     - Example
   * - ``ERF_ENABLE_NETCDF``
     - Enable all NetCDF-related features and I/O routines
     - ``-DERF_ENABLE_NETCDF=ON``
   * - ``ERF_ENABLE_HDF5``
     - Enable HDF5 support in AMReX backend (defaults to ``ERF_ENABLE_NETCDF`` value)
     - ``-DERF_ENABLE_HDF5=ON``
   * - ``NETCDF_DIR``
     - Specify installation prefix for custom NetCDF library
     - ``-DNETCDF_DIR=/path/to/netcdf``
   * - ``HDF5_PREFER_PARALLEL``
     - Hint to prioritize parallel HDF5 libraries
     - ``-DHDF5_PREFER_PARALLEL=ON``

.. list-table:: GNU Make Variables for I/O Libraries
   :header-rows: 1
   :widths: 30 50 20

   * - GNU Make Variable
     - Description
     - Example
   * - ``USE_NETCDF``
     - Enable NetCDF features and link libraries
     - ``make USE_NETCDF=TRUE``

**Platform-Specific Configuration**

* **Cray Systems (Perlmutter, Frontier, Polaris)** - Load ``cray-netcdf-hdf5parallel`` module. ERF's Cray detection automatically discovers module-specific search paths.

  .. code-block:: bash

     module load cray-netcdf-hdf5parallel
     cmake -DERF_ENABLE_NETCDF=ON ..

* **Workstations and Generic Clusters** - Install parallel NetCDF and HDF5 via package manager:

  .. code-block:: bash

     # Ubuntu/Debian
     sudo apt install libnetcdf-mpi-dev libhdf5-mpi-dev

     # Fedora/RHEL
     sudo dnf install netcdf-mpich-devel hdf5-mpich-devel

  Then configure:

  .. code-block:: bash

     cmake -DERF_ENABLE_NETCDF=ON ..

**Configuration Examples**

.. tab-set::

   .. tab-item:: CMake with custom path

      .. code-block:: bash

         cmake -DERF_ENABLE_NETCDF=ON \
               -DNETCDF_DIR=/opt/netcdf-parallel \
               ..

   .. tab-item:: CMake with environment variable

      .. code-block:: bash

         export NETCDF_DIR=/opt/netcdf-parallel
         cmake -DERF_ENABLE_NETCDF=ON ..

   .. tab-item:: GNU Make

      .. code-block:: bash

         make USE_NETCDF=TRUE

   .. tab-item:: Disable NetCDF

      .. code-block:: bash

         cmake -DERF_ENABLE_NETCDF=OFF ..

**Verifying NetCDF Configuration**

After building, verify NetCDF is correctly linked:

.. code-block:: bash

   # Check if NetCDF symbols are in the executable
   nm ERF3d.*.ex | grep nc_

   # Or use ldd to see linked libraries
   ldd ERF3d.*.ex | grep netcdf

Physics and Utility Libraries
------------------------------

Beyond core I/O, ERF integrates specialized libraries for advanced physics modeling and data handling.

**SHOC and P3 (E3SM Physics)**

SHOC (Simplified Higher-Order Closure) provides turbulence and cloud macrophysics modeling. P3 (Predicted Particle Properties) provides microphysics. Both depend on Kokkos (via EKAT submodule) for CPU/GPU portability.

**Prerequisites:**

.. code-block:: bash

   # Initialize E3SM submodules
   export ERF_DIR=/path/to/ERF
   source $ERF_DIR/Build/GNU_Ekat/eamxx_clone.sh

**Configuration:**

.. tab-set::

   .. tab-item:: CMake - SHOC

      .. code-block:: bash

         cmake -DERF_ENABLE_SHOC=ON \
               -DERF_ENABLE_MPI=ON \
               ..

   .. tab-item:: CMake - P3

      .. code-block:: bash

         cmake -DERF_ENABLE_P3=ON \
               -DERF_ENABLE_MPI=ON \
               ..

   .. tab-item:: GNU Make - SHOC

      .. code-block:: bash

         source Build/GNU_Ekat/ekat_build_commands.sh
         make USE_SHOC=TRUE USE_MPI=TRUE

   .. tab-item:: GNU Make - P3

      .. code-block:: bash

         source Build/GNU_Ekat/ekat_build_commands.sh
         make USE_P3=TRUE USE_MPI=TRUE

.. note::
   SHOC and P3 both require ``ERF_ENABLE_MPI=ON`` (or ``USE_MPI=TRUE``). They automatically enable ``ERF_ENABLE_EKAT=ON`` which provides Kokkos.

**RRTMGP (Radiation)**

RRTMGP (Rapid Radiative Transfer Model) provides radiation calculations.

**Prerequisites:**

Requires NetCDF and MPI:

.. code-block:: bash

   cmake -DERF_ENABLE_RRTMGP=ON \
         -DERF_ENABLE_NETCDF=ON \
         -DERF_ENABLE_MPI=ON \
         ..

**ZFP (Compression)**

ZFP provides lossy floating-point compression as an HDF5 filter. Reduces output data size in large-scale simulations.

.. warning::
   ZFP provides **lossy compression** - evaluate impact on your scientific goals before enabling.

**Configuration:**

Requires HDF5 compiled with ZFP plugin:

.. code-block:: bash

   cmake -DERF_ENABLE_NETCDF=ON \
         -DAMReX_HDF5_ZFP=ON \
         ..

**AMReX (Core Framework)**

AMReX provides core AMR infrastructure. ERF uses an internal submodule by default for version compatibility.

**Using Internal AMReX (Default):**

.. code-block:: bash

   # Submodule automatically used
   cmake -DERF_USE_INTERNAL_AMREX=ON ..  # Default

**Using External AMReX:**

.. code-block:: bash

   # Download external AMReX
   git clone https://github.com/AMReX-Codes/amrex.git
   cd amrex && mkdir build && cd build
   cmake -DCMAKE_INSTALL_PREFIX=/path/to/install ..
   make install

   # Point ERF to external build
   cmake -DERF_USE_INTERNAL_AMREX=OFF \
         -DAMReX_ROOT=/path/to/install \
         ..

For GNU Make, set ``AMREX_HOME``:

.. code-block:: bash

   export AMREX_HOME=/path/to/external/amrex
   make

Troubleshooting
---------------

Common library configuration issues and solutions.

**Library Not Found**

**Symptom:** CMake reports "Could not find NetCDF" or similar.

**Solutions:**

1. **Verify modules loaded:**

   .. code-block:: bash

      module list  # Check for cray-netcdf-hdf5parallel

2. **Specify path manually:**

   .. code-block:: bash

      cmake -DNETCDF_DIR=/path/to/netcdf ..

3. **Check pkg-config:**

   .. code-block:: bash

      pkg-config --list-all | grep netcdf
      echo $PKG_CONFIG_PATH

**MPI Linker Errors**

**Symptom:** "Undefined reference to MPI_Init" or other MPI symbols.

**Diagnosis:** Parallel ERF build is linking against serial versions of NetCDF/HDF5.

**Solutions:**

1. **Unload serial libraries:**

   .. code-block:: bash

      module unload netcdf hdf5
      module load cray-netcdf-hdf5parallel

2. **Clean build:**

   .. code-block:: bash

      make distclean
      cmake ..
      make

**Stale CMake Cache**

**Symptom:** Unexpected failures after changing modules or compilers.

**Diagnosis:** CMake cached old library locations.

**Solutions:**

1. **Full clean:**

   .. code-block:: bash

      make distclean
      # Or manually:
      rm -rf CMakeCache.txt CMakeFiles/

2. **Reconfigure:**

   .. code-block:: bash

      cmake ..

**Advanced Diagnosis**

For complex issues, use CMake's verbose logging:

.. code-block:: bash

   cmake --log-level=VERBOSE --log-context ..

Example output:

.. code-block:: text

   [ERF.Cray] Detected Cray Programming Environment (CRAYPE_VERSION=2.7.30)
   [ERF.Cray] Setting Cray compiler wrappers...
   [ERF.AMReX] Using internal AMReX submodule
   [ERF.NetCDF] Found NetCDF: /opt/cray/pe/netcdf/4.9.0.9

This reveals exact search paths and logic used by ``find_package`` commands.

**Verification Commands**

After successful build, verify library linkage:

.. code-block:: bash

   # Check linked libraries
   ldd ERF3d.*.ex

   # Check for specific symbols
   nm ERF3d.*.ex | grep -i "netcdf\|hdf5\|mpi"

   # Run executable with --describe
   ./ERF3d.*.ex --describe
