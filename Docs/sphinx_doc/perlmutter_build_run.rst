.. _sec:perlmutter_build_run:

Perlmutter (NERSC)
==================

Build and run guidance for ERF on NERSC Perlmutter.
For shared build concepts (machine profiles, Cray detection, script reference), see :ref:`sec:build:hpc`.

.. tab-set::

   .. tab-item:: Building with GNU Make

      Simple build using GNU compiler and CUDA:

      .. code-block:: bash

         # Load environment
         module load PrgEnv-gnu cudatoolkit cray-mpich cray-netcdf-hdf5parallel

         # Navigate to Exec
         cd ${ERF_HOME}/Exec

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

      **Executable location:** ``build/Exec/erf_exec`` (or ``install/bin/erf_exec`` if installed)

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
         mkdir -p $PSCRATCH/ERF/rundir
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
