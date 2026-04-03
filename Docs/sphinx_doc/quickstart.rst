.. _sec:build:quickstart:

Quickstart: Clone-Build-Run
============================

Copy-paste commands for common build scenarios.

Need HPC quickstart commands with machine profiles? Jump to :ref:`sec:build:quickstart:hpc-cmake`.

Quickstart on HPC Systems (Machine Profiles)
--------------------------------------------

Suggested path for most users (HPC or not): use the automated script workflow below.
Choose other workflows only when you need tighter control (see :ref:`sec:build:systems` for details).

.. tab-set::

   .. tab-item:: Automated Clone-Build-Run Script w/ CMake

      This is the recommended default:

      .. code-block:: bash

         export MACHINE_NAME=<my_machine>
         export SLURM_ACCOUNT=<my_account>
         export SALLOC_ACCOUNT=$SLURM_ACCOUNT
         export SBATCH_ACCOUNT=$SLURM_ACCOUNT
         
         git clone --recursive git@github.com:erf-model/ERF.git
         cd ERF
         source Build/machines/${MACHINE_NAME}_erf.profile
         ./Build/cmake_with_kokkos_many.sh

         # Run from install directory
         cd install/bin
         # Try local run
         mpiexec -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most

         # Preferred on shared HPC systems: submit batch script
         sbatch run_${MACHINE_NAME}_erf.sbatch

         # Interactive alternative:
         # salloc -A <account> -C gpu -N 1 -t 00:30:00
         # srun -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most

Use these when you are on an HPC system and want a machine-profile + build-script workflow. Each tab includes optional commented allocation/launch commands.
.. tab-set::

   .. tab-item:: Perlmutter (NERSC)

      .. code-block:: bash

         git clone --recursive https://github.com/erf-model/ERF.git
         cd ERF
         source Build/machines/perlmutter_erf.profile

         ERF_HOME=$(pwd) ./Build/cmake_with_kokkos_many_cuda.sh
         cd install/bin

         # Preferred on shared HPC systems: submit batch script
         sbatch run_perlmutter_erf.sbatch

         # Interactive alternative:
         # salloc -A <account> -C gpu -N 1 -t 00:30:00
         # srun -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most

      :download:`Download build+run snippet <scripts/quickstart/perlmutter_quickstart.sh>`
      Full guide: :doc:`perlmutter_build_run`

   .. tab-item:: Kestrel (NREL)

      .. code-block:: bash

         git clone --recursive https://github.com/erf-model/ERF.git
         cd ERF
         source Build/machines/kestrel_erf.profile

         ERF_HOME=$(pwd) ./Build/cmake_with_kokkos_many_cuda.sh
         cd install/bin

         # Preferred on shared HPC systems: submit batch script
         sbatch ../../Docs/sphinx_doc/scripts/quickstart/run.erf.aw.job_arena

         # Interactive alternative:
         # salloc -p gpu-h100s --nodes=1 --gpus-per-node=4 --time=00:30:00
         # srun -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most

      :download:`Download Kestrel custom build script <scripts/quickstart/cmake_kestrel_ERF.sh>`
      :download:`Download Kestrel sbatch example <scripts/quickstart/run.erf.aw.job_arena>`
      Full guide: :doc:`kestrel_build_run`

   .. tab-item:: Aurora (ALCF)

      .. code-block:: bash

         git clone --recursive https://github.com/erf-model/ERF.git
         cd ERF
         source Build/machines/aurora_erf.profile
         export NETCDF_DIR=<path-to-netcdf>

         ERF_HOME=$(pwd) ./Build/cmake_with_kokkos_many_sycl.sh
         cd install/bin

         # Preferred on shared HPC systems: submit batch script
         qsub submit_erf_aurora.pbs

         # Interactive alternative:
         # qsub -I -A <PROJECT> -q debug -l select=1 -l walltime=1:00:00 -l filesystems=flare
         # mpiexec -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most

      :download:`Download build+run snippet <scripts/quickstart/aurora_quickstart.sh>`
      Full guide: :doc:`aurora_build_run`

   .. tab-item:: Frontier (OLCF)

      .. code-block:: bash

         git clone --recursive https://github.com/erf-model/ERF.git
         cd ERF
         source Build/machines/frontier_erf.profile

         ERF_HOME=$(pwd) ./Build/cmake_with_kokkos_many_hip.sh
         cd install/bin

         # Preferred on shared HPC systems: submit batch script
         sbatch run_frontier_erf.sbatch

         # Interactive alternative:
         # salloc -A <project> -p batch -N 1 -t 00:30:00
         # srun -n 8 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most

      :download:`Download build+run snippet <scripts/quickstart/frontier_quickstart.sh>`
      More machine context: :ref:`sec:build:hpc`

   .. tab-item:: Polaris (ALCF)

      .. code-block:: bash

         git clone --recursive https://github.com/erf-model/ERF.git
         cd ERF
         source Build/machines/polaris_erf.profile

         ERF_HOME=$(pwd) ./Build/cmake_with_kokkos_many_cuda.sh
         cd install/bin

         # Preferred on shared HPC systems: submit batch script
         qsub submit_erf_polaris.pbs

         # Interactive alternative:
         # qsub -I -A <PROJECT> -q debug -l select=1 -l walltime=1:00:00
         # mpiexec -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most

      :download:`Download build+run snippet <scripts/quickstart/polaris_quickstart.sh>`
      More machine context: :ref:`sec:build:hpc`


Build with GNU Make
-------------------

Clone, build, and run with GNU Make:

.. tab-set::

   .. tab-item:: CPU Build

      .. code-block:: bash

         git clone --recursive git@github.com:erf-model/ERF.git
         cd ERF/Exec
         make COMP=gnu USE_MPI=TRUE
         mpiexec -n 4 ./ERF3d.gnu.TPROF.MPI.ex CanonicalTests/ABL/inputs_most

   .. tab-item:: GPU Build

      .. code-block:: bash

         git clone --recursive git@github.com:erf-model/ERF.git
         cd ERF/Exec
         make COMP=gnu USE_MPI=TRUE USE_CUDA=TRUE
         mpiexec -n 4 ./ERF3d.gnu.TPROF.MPI.CUDA.ex  CanonicalTests/ABL/inputs_most

Default GNU Make builds are done in ``ERF/Exec``. Use ``ERF/.Exec_dev/<test_name>`` only when working on development tests.

Build with CMake
----------------

Suggested path for most users (HPC or not): use the automated script workflow below.
Choose other workflows only when you need tighter control (see :ref:`sec:build:systems` for details).

.. tab-set::

   .. tab-item:: Automated Script (Suggested)

      This is the recommended default:

      .. code-block:: bash

         git clone --recursive git@github.com:erf-model/ERF.git
         cd ERF
         ./Build/cmake_with_kokkos_many.sh

         # Run from install directory
         cd install/bin  # or $ERF_INSTALL_DIR/bin if customized
         mpiexec -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most

      .. dropdown:: Optional: Customize build/source/install directories
         :icon: info

         Set ``ERF_BUILD_DIR``, ``ERF_SOURCE_DIR``, ``ERF_INSTALL_DIR``, or ``ERF_HOME`` to override defaults.

         .. code-block:: bash

            ERF_HOME=$(pwd) ./Build/cmake_with_kokkos_many.sh

            ERF_BUILD_DIR=$(pwd)/build_custom \
            ERF_INSTALL_DIR=$(pwd)/install_custom \
            ./Build/cmake_with_kokkos_many.sh

   .. tab-item:: Out-of-Source Build

      Separate build directory with optional install:

      .. code-block:: bash

         git clone --recursive git@github.com:erf-model/ERF.git
         cd ERF
         mkdir build && cd build
         ../Build/cmake.sh
         make -j
         make install  # optional

         # Run from build tree (without install)
         cd Exec
         mpiexec -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most

      .. dropdown:: Optional: Install and run from install directory
         :icon: info

         .. code-block:: bash

            make install
            cd ../../install/bin
            mpiexec -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most

   .. tab-item:: In Build/ Directory

      Build in Build/ directory (similar to GNU Make structure):

      .. code-block:: bash

         git clone --recursive git@github.com:erf-model/ERF.git
         cd ERF/Build
         ./cmake.sh
         make -j

         # Run from Exec subdirectory
         cd Exec
         mpiexec -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most

.. _sec:build:quickstart:hpc-cmake:

Cleanup
-------

Remove build artifacts:

.. tab-set::

   .. tab-item:: GNU Make

      .. code-block:: bash

         # From ``ERF/Exec`` (or ``ERF/.Exec_dev/<test_name>`` for dev-test builds)
         make clean           # Remove build artifacts
         make realclean       # Same as clean
         make cleanconfig     # Remove configuration only

   .. tab-item:: In Build/ Directory

      .. code-block:: bash

         # From ERF/Build (CMake generates Makefiles, so make targets work directly)
         make distclean       # Clean all CMake artifacts
         make uninstall       # Uninstall based on install-manifest.txt

         # If configuration failed, manually remove
         rm -rf CMakeFiles/ CMakeCache.txt

   .. tab-item:: Out-of-Source Build

      .. code-block:: bash

         # Using cmake --build (from ERF root)
         cmake --build build --target distclean  # Clean all CMake artifacts
         cmake --build build --target uninstall  # Uninstall based on install-manifest.txt

         # From inside build/ directory, Makefile targets also work:
         # cd build && make distclean

         # Or remove directories entirely
         rm -rf build/ install/

   .. tab-item:: Automated Script

      .. code-block:: bash

         # Using cmake --build (from ERF root)
         cmake --build build --target distclean  # Clean all CMake artifacts
         cmake --build build --target uninstall  # Uninstall based on install-manifest.txt

         # From inside build/ directory, Makefile targets also work:
         # cd build && make distclean

         # Or remove directories entirely
         rm -rf build/ install/

For detailed build options and troubleshooting, see :ref:`sec:build:systems`.
