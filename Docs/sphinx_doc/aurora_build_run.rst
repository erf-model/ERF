.. _sec:aurora_build_run:

Aurora (ALCF): Build and Run with SYCL
======================================

Copy-paste instructions for building and running ERF on ALCF Aurora using Intel GPUs (SYCL backend).

For general HPC build concepts, see :ref:`sec:build:hpc`. For library configuration details, see :ref:`sec:build:library`.

Prerequisites
-------------

- Aurora compute allocation (PBSPro scheduler)
- Access to required modules (MPICH, HDF5, CMake, Python)
- Environment variables:

  - ``ERF_HOME`` — path to ERF source directory
  - ``NETCDF_DIR`` — path to parallel NetCDF installation prefix (for NetCDF-enabled builds)

NetCDF Requirements
~~~~~~~~~~~~~~~~~~~

ERF requires a user-provided parallel NetCDF installation (set via ``NETCDF_DIR``).

- **C-only NetCDF** is sufficient for standard ERF NetCDF I/O
- **NetCDF-Fortran** is required only when enabling features that use it (e.g., Noah-MP land surface model)

For an example parallel NetCDF build with Fortran support, see:
`Lab-Notebooks/ERF-Noah-Coupling/software/netcdf <https://github.com/Lab-Notebooks/ERF-Noah-Coupling/tree/main/software/netcdf>`_

Quick checks for your NetCDF installation:

.. code-block:: bash

   # C headers and library (required)
   ls $NETCDF_DIR/include/netcdf.h
   ls $NETCDF_DIR/lib/libnetcdf*

   # Fortran module and library (only if using Noah-MP or other Fortran consumers)
   ls $NETCDF_DIR/include/netcdf.mod
   ls $NETCDF_DIR/lib/libnetcdff*

Quick Checks
~~~~~~~~~~~~

Verify your environment before building:

.. code-block:: bash

   module list
   which cmake mpicc mpicxx mpifort
   echo $ERF_HOME
   echo $NETCDF_DIR

.. tab-set::

   .. tab-item:: Interactive

      Interactive Build and Run
      -------------------------

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
            module load python/3.10.14
            module load cmake

            # Intel compilers for MPICH wrappers
            export MPICH_CC=icx
            export MPICH_CXX=icpx
            export MPICH_FC=ifx
            export MPICH_F90=ifx

      **3) Set environment variables**

      .. code-block:: bash

         # MPI and HDF5 paths (auto-detected from loaded modules)
         export MPI_HOME=$(which mpicc | sed s'/\/bin\/mpicc'//)
         export HDF5_HOME=$(which h5pfc | sed s'/\/bin\/h5pfc'//)
         export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HDF5_HOME/lib
         export C_INCLUDE_PATH=$C_INCLUDE_PATH:$HDF5_HOME/include

         # NetCDF path (edit this)
         export NETCDF_DIR=<PATH_TO_NETCDF>

      **4) Configure and build**

      .. code-block:: bash

         cd $ERF_HOME
         mkdir -p tmp_build_dir && cd tmp_build_dir

         export CRAYPE_LINK_TYPE=dynamic

         cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
            -DCMAKE_CXX_COMPILER:STRING=mpicxx \
            -DCMAKE_C_COMPILER:STRING=mpicc \
            -DCMAKE_Fortran_COMPILER:STRING=mpifort \
            -DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo \
            -DERF_DIM:STRING=3 \
            -DERF_ENABLE_PARTICLES:BOOL=ON \
            -DERF_ENABLE_MPI:BOOL=ON \
            -DERF_ENABLE_CUDA:BOOL=OFF \
            -DERF_ENABLE_SYCL:BOOL=ON \
            -DERF_ENABLE_TESTS:BOOL=OFF \
            -DERF_ENABLE_FCOMPARE:BOOL=ON \
            -DERF_ENABLE_DOCUMENTATION:BOOL=OFF \
            -DERF_ENABLE_NETCDF:BOOL=ON \
            -DNETCDF_DIR=$NETCDF_DIR \
            -DAMReX_PARALLEL_LINK_JOBS=8 \
            $ERF_HOME

         make -j

      .. note::
         If compilation runs out of memory, reduce parallelism: ``make -j4`` instead of ``make -j``.

      **5) Run**

      From within a PBS allocation (``$PBS_NODEFILE`` must be set):

      .. code-block:: bash

         cd $ERF_HOME/tmp_build_dir

         # Compute node count from PBS
         NNODES=$(wc -l < $PBS_NODEFILE)

         # MPI layout
         NRANKS=4          # MPI ranks per node
         NDEPTH=16         # Hardware threads per rank (spacing)
         NTHREADS=1        # OpenMP threads per rank
         NTOTRANKS=$((NNODES * NRANKS))

         echo "NNODES=$NNODES NTOTRANKS=$NTOTRANKS NRANKS=$NRANKS NTHREADS=$NTHREADS"

         mpiexec --np ${NTOTRANKS} -ppn ${NRANKS} -d ${NDEPTH} --cpu-bind depth \
            -env OMP_NUM_THREADS=${NTHREADS} \
            ./Exec/erf_exec <PATH_TO_INPUTS_FILE>

      **MPI layout parameters:**

      .. list-table::
         :header-rows: 1
         :widths: 20 60 20

         * - Parameter
           - Description
           - Default
         * - ``NRANKS``
           - MPI ranks per node
           - 4
         * - ``NDEPTH``
           - Hardware threads per rank (controls rank spacing)
           - 16
         * - ``NTHREADS``
           - OpenMP threads per rank (``OMP_NUM_THREADS``)
           - 1

   .. tab-item:: Batch

      Batch Submission
      ----------------

      Save the following as ``submit_erf_aurora.pbs`` (edit placeholders):

      .. code-block:: bash

         #!/bin/bash
         #PBS -A <PROJECT>
         #PBS -q <QUEUE>
         #PBS -l select=<NODES>
         #PBS -l walltime=<HH:MM:SS>
         #PBS -l filesystems=<FILESYSTEMS>
         #PBS -N erf_aurora
         #PBS -j oe
         #PBS -o erf_${PBS_JOBID}.out

         # -------------------------------------------------------------------
         # User configuration (edit these)
         # -------------------------------------------------------------------
         export ERF_HOME=<PATH_TO_ERF>
         export NETCDF_DIR=<PATH_TO_NETCDF>
         INPUTS_FILE=<PATH_TO_INPUTS_FILE>

         # MPI layout
         NRANKS=4          # MPI ranks per node
         NDEPTH=16         # Hardware threads per rank (spacing)
         NTHREADS=1        # OpenMP threads per rank

         # -------------------------------------------------------------------
         # Load software environment
         # -------------------------------------------------------------------
         source $ERF_HOME/Build/machines/aurora_erf.profile

         # MPI and HDF5 paths (auto-detected from loaded modules)
         export MPI_HOME=$(which mpicc | sed s'/\/bin\/mpicc'//)
         export HDF5_HOME=$(which h5pfc | sed s'/\/bin\/h5pfc'//)
         export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HDF5_HOME/lib
         export C_INCLUDE_PATH=$C_INCLUDE_PATH:$HDF5_HOME/include

         export CRAYPE_LINK_TYPE=dynamic

         # -------------------------------------------------------------------
         # Build (optional: skip if already built)
         # -------------------------------------------------------------------
         cd $ERF_HOME
         mkdir -p tmp_build_dir && cd tmp_build_dir

         cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
            -DCMAKE_CXX_COMPILER:STRING=mpicxx \
            -DCMAKE_C_COMPILER:STRING=mpicc \
            -DCMAKE_Fortran_COMPILER:STRING=mpifort \
            -DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo \
            -DERF_DIM:STRING=3 \
            -DERF_ENABLE_PARTICLES:BOOL=ON \
            -DERF_ENABLE_MPI:BOOL=ON \
            -DERF_ENABLE_CUDA:BOOL=OFF \
            -DERF_ENABLE_SYCL:BOOL=ON \
            -DERF_ENABLE_TESTS:BOOL=OFF \
            -DERF_ENABLE_FCOMPARE:BOOL=ON \
            -DERF_ENABLE_DOCUMENTATION:BOOL=OFF \
            -DERF_ENABLE_NETCDF:BOOL=ON \
            -DNETCDF_DIR=$NETCDF_DIR \
            -DAMReX_PARALLEL_LINK_JOBS=8 \
            $ERF_HOME

         make -j

         # -------------------------------------------------------------------
         # Run
         # -------------------------------------------------------------------
         NNODES=$(wc -l < $PBS_NODEFILE)
         NTOTRANKS=$((NNODES * NRANKS))

         echo "NNODES=$NNODES NTOTRANKS=$NTOTRANKS NRANKS=$NRANKS NTHREADS=$NTHREADS"

         mpiexec --np ${NTOTRANKS} -ppn ${NRANKS} -d ${NDEPTH} --cpu-bind depth \
            -env OMP_NUM_THREADS=${NTHREADS} \
            ./Exec/erf_exec ${INPUTS_FILE}

      **Submit the job:**

      .. code-block:: bash

         qsub submit_erf_aurora.pbs

      **Monitor:**

      .. code-block:: bash

         qstat -u $USER

      .. tip::
         **Build once, run many:** For production workflows, build ERF once (interactively
         or in a dedicated build job), then remove the build section from subsequent
         batch scripts. This saves allocation time.

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

See :ref:`InputFiles` for complete input file documentation.

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
      which cmake mpicc mpicxx mpifort
      echo $NETCDF_DIR
      ls $NETCDF_DIR/lib

   **Common causes:**

   - ``NETCDF_DIR`` not set or pointing to wrong location
   - Modules not loaded (run ``module list`` to verify)
   - Stale CMake cache — remove ``tmp_build_dir`` and reconfigure

.. dropdown:: Compilation runs out of memory
   :icon: alert
   :color: warning

   Reduce parallel compilation:

   .. code-block:: bash

      make -j4

   Or request more memory in your PBS allocation.

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

      mpiexec --np 4 -ppn 4 -d 16 --cpu-bind depth \
         -env OMP_NUM_THREADS=1 \
         ./Exec/erf_exec <inputs_file> max_step=10

   Check ``Backtrace.*`` files for stack traces.

For additional troubleshooting, see :ref:`sec:build:troubleshooting`.
