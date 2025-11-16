.. _sec:build:library:

Library Configuration
=====================

Configuring external libraries enables ERF's scientific capabilities, from parallel I/O on supercomputing platforms to GPU-accelerated physics packages. This guide provides technical details for managing ERF's library dependencies.

For basic build instructions, see :ref:`sec:build:systems`. For HPC-specific configurations, see :ref:`sec:build:hpc`.

Library Dependencies Overview
------------------------------

ERF integrates external libraries for core functionality and advanced physics:

.. list-table:: ERF Library Dependencies
   :header-rows: 1
   :widths: 20 40 15 25

   * - Library
     - Description
     - Requirement
     - Build Options
   * - **AMReX**
     - Block-structured adaptive mesh refinement framework
     - Required
     - ``-DERF_USE_INTERNAL_AMREX=ON`` / ``AMREX_HOME``
   * - **NetCDF**
     - High-level I/O for plotfiles and initial conditions
     - Optional
     - ``-DERF_ENABLE_NETCDF=ON`` / ``USE_NETCDF=TRUE``
   * - **HDF5**
     - Parallel I/O backend for NetCDF
     - Optional
     - ``-DERF_ENABLE_HDF5=ON``
   * - **SHOC**
     - Turbulence and cloud macrophysics (:ref:`SHOC`)
     - Optional
     - ``-DERF_ENABLE_SHOC=ON`` / ``USE_SHOC=TRUE``
   * - **P3**
     - Microphysics (:ref:`Microphysics`)
     - Optional
     - ``-DERF_ENABLE_P3=ON`` / ``USE_P3=TRUE``
   * - **RRTMGP**
     - Radiation model
     - Optional
     - ``-DERF_ENABLE_RRTMGP=ON`` / ``USE_RRTMGP=TRUE``
   * - **ZFP**
     - Lossy data compression for HDF5
     - Optional
     - ``-DAMReX_HDF5_ZFP=ON``

**Parallel I/O Consistency**

For MPI-enabled builds with parallel I/O (``ERF_ENABLE_NETCDF=ON``), both NetCDF and HDF5 must be parallel-enabled. Linking a parallel ERF build against serial libraries results in linker errors or runtime failures.

NetCDF and HDF5
---------------

NetCDF and HDF5 provide ERF's I/O capabilities for large-scale simulations.

**What They Do**

* **NetCDF** - High-level API for structured data:

  - Reads initial conditions from WRF/Metgrid files
  - Writes AMReX plotfiles in NetCDF format

* **HDF5** - Low-level parallel storage:

  - Coordinates writes from multiple MPI ranks
  - Manages parallel filesystem operations

**NetCDF Detection Hierarchy**

Both build systems search for NetCDF libraries in priority order, trying multiple pkg-config variants before falling back to system paths.

.. tab-set::

   .. tab-item:: CMake

      CMake uses ``FindNetCDF.cmake`` with the following search order:

      .. list-table:: CMake NetCDF Detection Priority
         :header-rows: 1
         :widths: 10 25 25 40

         * - Priority
           - Method
           - Variables/Commands
           - Notes
         * - 1
           - CMake Option
           - ``NETCDF_DIR``
           - ``-DNETCDF_DIR=...`` (highest priority)
         * - 2
           - Environment Variable
           - ``NETCDF_DIR``
           - ``export NETCDF_DIR=...``
         * - 3
           - pkg-config
           - ``netcdf``
           - First pkg-config variant (with ``MPICH_DIR`` for Cray)
         * - 4
           - pkg-config
           - ``netcdf-mpi``
           - Parallel NetCDF variant
         * - 5
           - pkg-config
           - ``netcdf_parallel``
           - Alternative parallel naming
         * - 6
           - pkg-config
           - ``netcdf-cxx4_parallel``
           - C++ parallel variant
         * - 7
           - Manual Search
           - ``find_library()``, ``find_path()``
           - Searches ``/usr/lib``, ``/usr/include``

      **Cray Systems:** CMake automatically augments ``PKG_CONFIG_PATH`` with ``$MPICH_DIR/lib/pkgconfig`` when detecting Cray environments.

   .. tab-item:: GNU Make

      GNU Make uses ``Make.ERF`` with the following search order:

      .. list-table:: GNU Make NetCDF Detection Priority
         :header-rows: 1
         :widths: 10 25 25 40

         * - Priority
           - Method
           - Variables/Commands
           - Notes
         * - 1
           - pkg-config
           - ``netcdf``
           - With ``MPICH_DIR`` augmentation for Cray systems
         * - 2
           - pkg-config
           - ``netcdf-mpi``
           - Parallel NetCDF variant
         * - 3
           - pkg-config
           - ``netcdf-cxx4_parallel``
           - C++ parallel variant
         * - 4
           - pkg-config
           - ``netcdf_parallel``
           - Alternative parallel naming
         * - 5
           - Error
           - ``$(error NetCDF not found)``
           - No fallback to system paths

      **Cray Systems:** GNU Make manually augments ``PKG_CONFIG_PATH`` with ``$MPICH_DIR/lib/pkgconfig`` when ``MPICH_DIR`` is set.

   .. tab-item:: Key Differences

      .. list-table::
         :header-rows: 1
         :widths: 30 35 35

         * - Feature
           - CMake
           - GNU Make
         * - **User Path Override**
           - ``-DNETCDF_DIR=...`` or ``export NETCDF_DIR``
           - Not supported (pkg-config only)
         * - **Cray Integration**
           - Automatic ``MPICH_DIR`` detection
           - Manual ``MPICH_DIR`` check
         * - **Fallback Mechanism**
           - System path search via ``find_library()``
           - Build fails with error
         * - **pkg-config Variants**
           - 4 variants tried
           - 4 variants tried (different order)
         * - **Diagnostic Output**
           - Detailed logs with ``--log-level=VERBOSE``
           - Build error shows tried variants

**Build Options**

.. dropdown:: CMake Options for I/O Libraries
   :icon: code

   .. list-table::
      :header-rows: 1
      :widths: 30 50 20

      * - CMake Option
        - Description
        - Example
      * - ``ERF_ENABLE_NETCDF``
        - Master switch to enable all NetCDF-related features and I/O routines
        - ``-DERF_ENABLE_NETCDF=ON``
      * - ``ERF_ENABLE_HDF5``
        - Enables HDF5 support within AMReX backend (defaults to ``ERF_ENABLE_NETCDF`` value)
        - ``-DERF_ENABLE_HDF5=ON``
      * - ``NETCDF_DIR``
        - Specifies installation prefix for custom NetCDF library
        - ``-DNETCDF_DIR=/path/to/netcdf``
      * - ``HDF5_PREFER_PARALLEL``
        - Hint for find module to prioritize parallel HDF5 libraries
        - ``-DHDF5_PREFER_PARALLEL=ON``

.. dropdown:: GNU Make Variables for I/O Libraries
   :icon: code

   .. list-table::
      :header-rows: 1
      :widths: 30 50 20

      * - Variable
        - Description
        - Example
      * - ``USE_NETCDF``
        - Enable NetCDF
        - ``make USE_NETCDF=TRUE``

**Platform Configuration**

.. tab-set::

   .. tab-item:: Cray Systems

      Load the parallel NetCDF module:

      .. code-block:: bash

         module load cray-netcdf-hdf5parallel
         cmake -DERF_ENABLE_NETCDF=ON ..

      ERF's Cray detection automatically discovers module paths.

   .. tab-item:: Workstations

      Install parallel libraries via package manager:

      .. code-block:: bash

         # Ubuntu/Debian
         sudo apt install libnetcdf-mpi-dev libhdf5-mpi-dev

         # Fedora/RHEL
         sudo dnf install netcdf-mpich-devel hdf5-mpich-devel

      Then configure:

      .. code-block:: bash

         cmake -DERF_ENABLE_NETCDF=ON ..

   .. tab-item:: Custom Installation

      Specify installation path:

      .. code-block:: bash

         cmake -DERF_ENABLE_NETCDF=ON \
               -DNETCDF_DIR=/opt/netcdf-parallel \
               ..

      Or use environment variable:

      .. code-block:: bash

         export NETCDF_DIR=/opt/netcdf-parallel
         cmake -DERF_ENABLE_NETCDF=ON ..

**Verification**

After building, verify NetCDF linkage:

.. code-block:: bash

   # Check linked libraries
   ldd ERF3d.*.ex | grep netcdf

   # Check for NetCDF symbols
   nm ERF3d.*.ex | grep nc_

   # View build configuration
   ./ERF3d.*.ex --describe

Physics Libraries
-----------------

ERF integrates specialized libraries for advanced atmospheric physics modeling.

SHOC (Turbulence)
~~~~~~~~~~~~~~~~~

SHOC (Simplified Higher-Order Closure) provides turbulence and cloud macrophysics from the E3SM project. For theory and implementation details, see :ref:`SHOC` in the PBL schemes documentation.

**Prerequisites:**

.. code-block:: bash

   # Initialize E3SM submodules
   export ERF_DIR=/path/to/ERF
   source $ERF_DIR/Build/GNU_Ekat/eamxx_clone.sh

**Configuration:**

.. tab-set::

   .. tab-item:: CMake

      .. code-block:: bash

         cmake -DERF_ENABLE_SHOC=ON \
               -DERF_ENABLE_MPI=ON \
               ..

   .. tab-item:: GNU Make

      .. code-block:: bash

         source Build/GNU_Ekat/ekat_build_commands.sh
         make USE_SHOC=TRUE USE_MPI=TRUE

.. note::
   SHOC requires MPI and automatically enables EKAT (provides Kokkos for GPU support).

P3 (Microphysics)
~~~~~~~~~~~~~~~~~

P3 (Predicted Particle Properties) provides microphysics modeling from E3SM. For theory details, see :ref:`Microphysics` documentation.

**Prerequisites:**

.. code-block:: bash

   # Initialize E3SM submodules (same as SHOC)
   export ERF_DIR=/path/to/ERF
   source $ERF_DIR/Build/GNU_Ekat/eamxx_clone.sh

**Configuration:**

.. tab-set::

   .. tab-item:: CMake

      .. code-block:: bash

         cmake -DERF_ENABLE_P3=ON \
               -DERF_ENABLE_MPI=ON \
               ..

   .. tab-item:: GNU Make

      .. code-block:: bash

         source Build/GNU_Ekat/ekat_build_commands.sh
         make USE_P3=TRUE USE_MPI=TRUE

.. note::
   P3 requires MPI and automatically enables EKAT (provides Kokkos).

RRTMGP (Radiation)
~~~~~~~~~~~~~~~~~~

RRTMGP (Rapid Radiative Transfer Model for GCMs) provides radiation calculations.

**Prerequisites:**

Requires NetCDF and MPI:

.. code-block:: bash

   cmake -DERF_ENABLE_RRTMGP=ON \
         -DERF_ENABLE_NETCDF=ON \
         -DERF_ENABLE_MPI=ON \
         ..

Or with GNU Make:

.. code-block:: bash

   make USE_RRTMGP=TRUE USE_NETCDF=TRUE USE_MPI=TRUE

.. note::
   RRTMGP automatically enables EKAT (provides Kokkos) and requires both NetCDF and MPI.

.. dropdown:: Advanced: ZFP Compression
   :icon: tools

   ZFP provides lossy floating-point compression as an HDF5 filter, reducing output data size.

   .. warning::
      ZFP provides **lossy compression** - evaluate impact on your scientific goals before enabling.

   **Requirements:**

   * HDF5 compiled with ZFP plugin

   **Configuration:**

   .. code-block:: bash

      cmake -DERF_ENABLE_NETCDF=ON \
            -DAMReX_HDF5_ZFP=ON \
            ..

AMReX
-----

AMReX provides ERF's core AMR infrastructure, data structures, and parallelization. ERF uses an internal submodule by default for version compatibility.

**Default (Internal Submodule):**

.. code-block:: bash

   # Automatically uses internal AMReX
   cmake ..  # Default: ERF_USE_INTERNAL_AMREX=ON

**External AMReX:**

For shared AMReX installations:

.. code-block:: bash

   # Build and install external AMReX
   git clone https://github.com/AMReX-Codes/amrex.git
   cd amrex && mkdir build && cd build
   cmake -DCMAKE_INSTALL_PREFIX=/opt/amrex ..
   make install

   # Configure ERF to use external build
   cmake -DERF_USE_INTERNAL_AMREX=OFF \
         -DAMReX_ROOT=/opt/amrex \
         ..

**GNU Make:**

.. code-block:: bash

   # Default (uses submodule)
   make

   # External AMReX
   export AMREX_HOME=/path/to/external/amrex
   make

Common Library Issues
---------------------

**Library Not Found**

**Symptom:** CMake reports "Could not find NetCDF" or similar.

**Solutions:**

1. Verify modules:

   .. code-block:: bash

      module list  # Check for cray-netcdf-hdf5parallel

2. Specify path:

   .. code-block:: bash

      cmake -DNETCDF_DIR=/path/to/netcdf ..

3. Check pkg-config:

   .. code-block:: bash

      pkg-config --list-all | grep netcdf
      echo $PKG_CONFIG_PATH

**MPI Linker Errors**

**Symptom:** "Undefined reference to MPI_Init" or other MPI symbols.

**Cause:** Parallel ERF linking against serial NetCDF/HDF5.

**Solutions:**

1. Load parallel libraries:

   .. code-block:: bash

      module unload netcdf hdf5
      module load cray-netcdf-hdf5parallel

2. Clean rebuild:

   .. code-block:: bash

      make distclean
      cmake ..
      make

**Advanced Diagnosis**

Use verbose logging to see CMake's search process:

.. code-block:: bash

   cmake --log-level=VERBOSE --log-context ..

Example output:

.. code-block:: text

   [ERF.Cray] Detected Cray Programming Environment
   [ERF.Cray] Setting Cray compiler wrappers...
   [ERF.AMReX] Using internal AMReX submodule
   [ERF.NetCDF] Found NetCDF: /opt/cray/pe/netcdf/4.9.0.9

**Verifying Library Linkage**

Confirm library linkage after successful build:

.. code-block:: bash

   # Check linked libraries
   ldd ERF3d.*.ex

   # Check for symbols
   nm ERF3d.*.ex | grep -i "netcdf\|hdf5\|mpi"

   # View build info
   ./ERF3d.*.ex --describe
