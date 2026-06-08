.. _sec:aurora_build_run:

Aurora (ALCF): Build and Run with SYCL
======================================

Copy-paste instructions for building and running ERF on ALCF Aurora using Intel GPUs (SYCL backend).

For general HPC build concepts, see :ref:`sec:build:hpc`. For library configuration details, see :ref:`sec:build:library`.

Prerequisites
-------------

- Aurora compute allocation (PBSPro scheduler)
- Access to required modules (loaded via ``Build/machines/aurora_erf.profile``)
- Environment variables set by the modules after sourcing the profile:

  - ``NETCDF_C_ROOT`` — path to NetCDF-C installation (set by ``netcdf-c`` module)
  - ``NETCDF_FORTRAN_ROOT`` — path to NetCDF-Fortran installation (set by ``netcdf-fortran`` module)
  - ``HDF5_ROOT`` — path to HDF5 installation (set by ``hdf5`` module)

NetCDF Requirements
~~~~~~~~~~~~~~~~~~~

Aurora provides system NetCDF modules (``netcdf-c``, ``netcdf-fortran``, ``netcdf-cxx4``) that are
loaded automatically by ``Build/machines/aurora_erf.profile``. No user-built NetCDF installation
is required.

Quick checks for your NetCDF installation after sourcing the profile:

.. code-block:: bash

   # C headers and library (required)
   ls $NETCDF_C_ROOT/include/netcdf.h
   ls $NETCDF_C_ROOT/lib64/libnetcdf*

   # Fortran module and library (required for Noah-MP)
   ls $NETCDF_FORTRAN_ROOT/include/netcdf.mod
   ls $NETCDF_FORTRAN_ROOT/lib64/libnetcdff*

   # HDF5 library
   ls $HDF5_ROOT/lib/libhdf5*

Quick Checks
~~~~~~~~~~~~

Verify your environment before building:

.. code-block:: bash

   module list
   which cmake mpicc mpicxx mpifort icpx
   echo $NETCDF_C_ROOT
   echo $NETCDF_FORTRAN_ROOT
   echo $HDF5_ROOT

.. tab-set::

   .. tab-item:: Interactive

      .. rubric:: Interactive Build and Run

      **1) Get an interactive allocation**

      .. code-block:: bash

         qsub -I -A <PROJECT> -q <QUEUE> -l select=<NODES> -l walltime=<HH:MM:SS> -l filesystems=<FILESYSTEMS>

      Example with typical debug settings:

      .. code-block:: bash

         qsub -I -A <PROJECT> -q debug -l select=1 -l walltime=1:00:00 -l filesystems=flare

      **2) Load software environment**

      Source the Aurora machine profile (recommended):

      .. code-block:: bash

         export ERF_HOME=<PATH_TO_ERF>
         source $ERF_HOME/Build/machines/aurora_erf.profile

      .. dropdown:: Alternative: Manual module loads
         :icon: info

         If the profile is outdated or unavailable, load modules manually:

         .. code-block:: bash

            module load mpich/opt/4.2.3-intel
            module load hdf5/1.14.6
            module load netcdf-cxx4
            module load netcdf-c
            module load netcdf-fortran
            module load python/3.10.14
            module load cmake

            # Intel compilers for MPICH wrappers
            export MPICH_CC=icx
            export MPICH_CXX=icpx
            export MPICH_FC=ifx
            export MPICH_F90=ifx

            # Derive NETCDF_FORTRAN_ROOT if not set by the module
            if [[ -z "${NETCDF_FORTRAN_ROOT:-}" ]]; then
              export NETCDF_FORTRAN_ROOT=$(module show netcdf-fortran 2>&1 | \
                sed -n 's/.*setenv("NETCDF_FORTRAN_ROOT","\([^"]*\)").*/\1/p')
            fi

      **3) Set search paths**

      .. code-block:: bash

         export CPPFLAGS="-I${NETCDF_C_ROOT}/include -I${NETCDF_FORTRAN_ROOT}/include -I${HDF5_ROOT}/include"
         export LDFLAGS="-L${NETCDF_C_ROOT}/lib64 -L${NETCDF_FORTRAN_ROOT}/lib64 -L${HDF5_ROOT}/lib"
         export LD_LIBRARY_PATH="${NETCDF_C_ROOT}/lib64:${NETCDF_FORTRAN_ROOT}/lib64:${HDF5_ROOT}/lib:${LD_LIBRARY_PATH:-}"

      **4) Configure and build**

      .. code-block:: bash

         cd $ERF_HOME
         git submodule update --init --recursive

         cmake -S . -B build_aurora \
           -DCMAKE_INSTALL_PREFIX="$(pwd)/install_aurora" \
           -DCMAKE_BUILD_TYPE=Release \
           -DCMAKE_C_COMPILER=mpicc \
           -DCMAKE_CXX_COMPILER=mpicxx \
           -DCMAKE_Fortran_COMPILER=mpifort \
           -DCMAKE_PREFIX_PATH="${NETCDF_C_ROOT};${NETCDF_FORTRAN_ROOT};${HDF5_ROOT}" \
           -DCMAKE_CXX_FLAGS="-fsycl-max-parallel-link-jobs=8 --offload-compress -flink-huge-device-code" \
           -DERF_ENABLE_MPI=ON \
           -DERF_ENABLE_NETCDF=ON \
           -DERF_ENABLE_HDF5=ON \
           -DERF_ENABLE_NOAHMP=ON \
           -DERF_ENABLE_RRTMGP=ON \
           -DERF_ENABLE_SYCL=ON \
           -DAMReX_GPU_BACKEND=SYCL \
           -DAMReX_INTEL_ARCH=pvc \
           -DAMReX_SYCL_AOT=ON \
           -DAMReX_SYCL_SPLIT_KERNEL=NO \
           -DKokkos_ENABLE_SERIAL=ON \
           -DKokkos_ENABLE_SYCL=ON \
           -DKokkos_ENABLE_SYCL_RELOCATABLE_DEVICE_CODE=ON \
           -DKokkos_ARCH_INTEL_PVC=ON

         cmake --build build_aurora -j 10

      .. note::
         The first full build takes significant time — SYCL AOT compilation for Intel PVC is
         expensive. If compilation runs out of memory on the login node, submit the build
         as an interactive PBS job (``qsub -I -l select=1 -q debug ...``) or reduce parallelism:
         ``cmake --build build_aurora -j 4``.

      **CMake flag reference:**

      .. list-table::
         :header-rows: 1
         :widths: 40 60

         * - Flag
           - Purpose
         * - ``-DERF_ENABLE_MPI=ON``
           - MPI parallelism
         * - ``-DERF_ENABLE_NETCDF=ON``
           - NetCDF I/O (needed for real cases / WPS initialization)
         * - ``-DERF_ENABLE_HDF5=ON``
           - Parallel HDF5 I/O
         * - ``-DERF_ENABLE_NOAHMP=ON``
           - Noah-MP land-surface model
         * - ``-DERF_ENABLE_RRTMGP=ON``
           - RRTMGP radiation package
         * - ``-DERF_ENABLE_SYCL=ON``
           - GPU offload via SYCL
         * - ``-DAMReX_GPU_BACKEND=SYCL``
           - AMReX SYCL backend
         * - ``-DAMReX_INTEL_ARCH=pvc``
           - Target Intel Data Center GPU Max (PVC)
         * - ``-DAMReX_SYCL_AOT=ON``
           - Ahead-of-time compilation for PVC (slower build, faster runtime)
         * - ``-DAMReX_SYCL_SPLIT_KERNEL=NO``
           - Disable kernel splitting (required for large kernels on PVC)
         * - ``-DKokkos_ENABLE_SERIAL=ON``
           - Kokkos serial backend (always required alongside SYCL)
         * - ``-DKokkos_ENABLE_SYCL=ON``
           - Kokkos SYCL backend (used by RRTMGP and Noah-MP)
         * - ``-DKokkos_ENABLE_SYCL_RELOCATABLE_DEVICE_CODE=ON``
           - Required for linking Kokkos SYCL device code across translation units
         * - ``-DKokkos_ARCH_INTEL_PVC=ON``
           - Kokkos target architecture for Intel PVC

      **5) Run**

      From within a PBS allocation (``$PBS_NODEFILE`` must be set):

      .. code-block:: bash

         cd $ERF_HOME

         # Compute node count from PBS
         NNODES=$(wc -l < $PBS_NODEFILE)

         # MPI layout (Aurora: 12 tiles/GPUs per node)
         NRANKS=12         # MPI ranks per node (one per GPU tile)
         NDEPTH=8          # Hardware threads per rank
         NTHREADS=1        # OpenMP threads per rank
         NTOTRANKS=$((NNODES * NRANKS))

         echo "NNODES=$NNODES NTOTRANKS=$NTOTRANKS NRANKS=$NRANKS NTHREADS=$NTHREADS"

         mpiexec --np ${NTOTRANKS} --ppn ${NRANKS} --depth=${NDEPTH} --cpu-bind=depth \
           -env OMP_NUM_THREADS=${NTHREADS} \
           build_aurora/Exec/<case>/erf_<case> <PATH_TO_INPUTS_FILE>

      **MPI layout parameters:**

      .. list-table::
         :header-rows: 1
         :widths: 20 60 20

         * - Parameter
           - Description
           - Default
         * - ``NRANKS``
           - MPI ranks per node (Aurora has 12 GPU tiles per node)
           - 12
         * - ``NDEPTH``
           - Hardware threads per rank (controls rank spacing)
           - 8
         * - ``NTHREADS``
           - OpenMP threads per rank (``OMP_NUM_THREADS``)
           - 1

   .. tab-item:: Batch

      .. rubric:: Batch Submission

      Save the following as ``submit_erf_aurora.pbs`` (edit placeholders):

      .. code-block:: bash

         #!/bin/bash -l
         #PBS -A <PROJECT>
         #PBS -q <QUEUE>
         #PBS -l select=<NODES>
         #PBS -l walltime=<HH:MM:SS>
         #PBS -l filesystems=home:flare
         #PBS -N erf_aurora
         #PBS -j oe
         #PBS -o erf_${PBS_JOBID}.out

         set -euo pipefail

         # -------------------------------------------------------------------
         # User configuration (edit these)
         # -------------------------------------------------------------------
         export ERF_HOME=<PATH_TO_ERF>
         CASE_DIR=<PATH_TO_RUN_DIR>
         EXE=$ERF_HOME/build_aurora/Exec/<case>/erf_<case>
         INPUTS=<inputs_filename>

         # MPI layout (Aurora: 12 GPU tiles per node)
         NRANKS=12         # MPI ranks per node
         NDEPTH=8          # Hardware threads per rank
         NTHREADS=1        # OpenMP threads per rank

         # -------------------------------------------------------------------
         # Load software environment
         # -------------------------------------------------------------------
         source $ERF_HOME/Build/machines/aurora_erf.profile

         export CPPFLAGS="-I${NETCDF_C_ROOT}/include -I${NETCDF_FORTRAN_ROOT}/include -I${HDF5_ROOT}/include"
         export LDFLAGS="-L${NETCDF_C_ROOT}/lib64 -L${NETCDF_FORTRAN_ROOT}/lib64 -L${HDF5_ROOT}/lib"
         export LD_LIBRARY_PATH="${NETCDF_C_ROOT}/lib64:${NETCDF_FORTRAN_ROOT}/lib64:${HDF5_ROOT}/lib:${LD_LIBRARY_PATH:-}"
         export OMP_NUM_THREADS=${NTHREADS}

         # -------------------------------------------------------------------
         # Validate
         # -------------------------------------------------------------------
         if [[ ! -x "$EXE" ]]; then
           echo "ERROR: executable not found: $EXE"
           exit 1
         fi

         if [[ ! -f "$CASE_DIR/$INPUTS" ]]; then
           echo "ERROR: input file not found: $CASE_DIR/$INPUTS"
           exit 1
         fi

         # -------------------------------------------------------------------
         # Run
         # -------------------------------------------------------------------
         cd "$CASE_DIR"

         NNODES=$(wc -l < $PBS_NODEFILE)
         NTOTRANKS=$((NNODES * NRANKS))

         echo "NNODES=$NNODES NTOTRANKS=$NTOTRANKS NRANKS=$NRANKS NTHREADS=$NTHREADS"

         mpiexec --np ${NTOTRANKS} --ppn ${NRANKS} --depth=${NDEPTH} --cpu-bind=depth \
           "$EXE" "$INPUTS"

      **Submit the job:**

      .. code-block:: bash

         qsub submit_erf_aurora.pbs

      **Monitor:**

      .. code-block:: bash

         qstat -u $USER
         qstat -f <jobid>

      .. tip::
         **Build once, run many:** Build ERF once (interactively or in a dedicated build job),
         then reuse the executable in subsequent batch scripts. SYCL AOT compilation is expensive
         — rebuilding per job wastes allocation time.

Example Inputs File
-------------------

The ABL case (``inputs_ml_most``) demonstrates a typical configuration:

.. code-block:: text

   # Problem setup
   erf.prob_name = "ABL"
   max_step = 4000

   # Domain: 1024^3 m on 64^3 grid
   geometry.prob_extent = 1024 1024 1024
   amr.n_cell           = 64 64 64
   geometry.is_periodic = 1 1 0

   # Boundary conditions
   zlo.type = "surface_layer"
   zhi.type = "SlipWall"

   # MOST surface layer
   erf.most.z0        = 0.1
   erf.most.zref      = 8.0
   erf.most.surf_temp = 1.1

   # Time stepping
   erf.fixed_dt = 1.0

   # Refinement
   amr.max_level      = 1
   amr.ref_ratio_vect = 20 20 1
   erf.coupling_type  = "TwoWay"

   # Output
   erf.plot_file_1 = plt
   erf.plot_int_1  = 1
   erf.check_file  = chk
   erf.check_int   = 100

   # Physics
   erf.les_type = "Smagorinsky"
   erf.Cs       = 0.1

See :ref:`sec:running` for complete input file documentation.

Expected Output
---------------

After a successful run, the working directory will contain:

- ``plt*`` — Plotfiles (AMReX format) at intervals set by ``erf.plot_int_1``
- ``chk*`` — Checkpoint files at intervals set by ``erf.check_int``
- ``Backtrace.*`` — Stack traces (if errors occurred)

Troubleshooting
---------------

.. dropdown:: CMake configuration fails
   :icon: alert
   :color: warning

   **Check environment:**

   .. code-block:: bash

      module list
      which cmake mpicc mpicxx mpifort icpx
      echo $NETCDF_C_ROOT
      echo $NETCDF_FORTRAN_ROOT
      echo $HDF5_ROOT
      ls $NETCDF_C_ROOT/include/netcdf.h
      ls $NETCDF_FORTRAN_ROOT/include/netcdf.mod

   **Common causes:**

   - ``NETCDF_C_ROOT`` or ``NETCDF_FORTRAN_ROOT`` not set — re-source ``Build/machines/aurora_erf.profile``
   - Modules not loaded (run ``module list`` to verify ``netcdf-c``, ``netcdf-fortran``, ``hdf5`` are present)
   - Stale CMake cache — remove ``build_aurora/`` and reconfigure

.. dropdown:: Compilation runs out of memory
   :icon: alert
   :color: warning

   Reduce parallel compilation:

   .. code-block:: bash

      cmake --build build_aurora -j 4

   Or submit the build as an interactive PBS job to get a full compute node:

   .. code-block:: bash

      qsub -I -A <PROJECT> -q debug -l select=1 -l walltime=1:00:00 -l filesystems=flare

.. dropdown:: mpiexec fails with PBS errors
   :icon: alert
   :color: warning

   Ensure you are inside a PBS allocation:

   .. code-block:: bash

      echo $PBS_NODEFILE
      cat $PBS_NODEFILE

   If empty or unset, you are not in a PBS job context. Use ``qsub -I ...`` for interactive or submit a batch script.

.. dropdown:: Runtime errors or crashes
   :icon: alert
   :color: warning

   **Quick diagnostic run:**

   .. code-block:: bash

      mpiexec --np 12 --ppn 12 --depth=8 --cpu-bind=depth \
        -env OMP_NUM_THREADS=1 \
        build_aurora/Exec/<case>/erf_<case> <inputs_file> max_step=10

   Check ``Backtrace.*`` files for stack traces.

.. dropdown:: NETCDF_FORTRAN_ROOT not set after sourcing profile
   :icon: alert
   :color: warning

   If the ``netcdf-fortran`` module does not export ``NETCDF_FORTRAN_ROOT``, derive it manually:

   .. code-block:: bash

      export NETCDF_FORTRAN_ROOT=$(module show netcdf-fortran 2>&1 | \
        sed -n 's/.*setenv("NETCDF_FORTRAN_ROOT","\([^"]*\)").*/\1/p')
      echo $NETCDF_FORTRAN_ROOT   # should be non-empty

For additional troubleshooting, see :ref:`sec:build:troubleshooting`.
