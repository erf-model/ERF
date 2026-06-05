.. _sec:build:quickstart:

Quickstart: Clone-Build-Run
============================

Copy-paste commands for common build scenarios. If you are on a supported HPC system,
start with the HPC section below — it is the fastest path to a working build. If you
are on a workstation or an unlisted system, start with :ref:`sec:build:quickstart:cmake`
or :ref:`sec:build:quickstart:gnumake`.

.. _sec:build:quickstart:hpc-cmake:

Quickstart on HPC Systems
--------------------------

Each tab covers only the machine-specific setup and build command. Shared guidance for
directory layout and batch vs interactive runs is listed once below the tabs.

.. tab-set::

   .. tab-item:: Perlmutter (NERSC)

      .. code-block:: bash

         git clone --recursive https://github.com/erf-model/ERF.git
         cd ERF
         source Build/machines/perlmutter_erf.profile

         ERF_HOME=$(pwd) ./Build/cmake_with_kokkos_many_cuda.sh
         cd install/bin
         # Next: choose a run mode in "Run After Build" below

      :download:`Download build+run snippet <scripts/quickstart/perlmutter_quickstart.sh>`
      :download:`Download Perlmutter sbatch example <scripts/quickstart/run_perlmutter_erf.sbatch>`
      Full guide: :doc:`perlmutter_build_run`

   .. tab-item:: Kestrel (NREL)

      .. code-block:: bash

         git clone --recursive https://github.com/erf-model/ERF.git
         cd ERF
         source Build/machines/kestrel_erf.profile

         ERF_HOME=$(pwd) ./Build/cmake_with_kokkos_many_cuda.sh
         cd install/bin
         # Next: choose a run mode in "Run After Build" below

      :download:`Download Kestrel custom build script <scripts/quickstart/cmake_kestrel_ERF.sh>`
      :download:`Download Kestrel sbatch example <scripts/quickstart/run.erf.aw.job_arena>`
      Full guide: :doc:`kestrel_build_run`

   .. tab-item:: Frontier (OLCF)

      .. code-block:: bash

         git clone --recursive https://github.com/erf-model/ERF.git
         cd ERF
         source Build/machines/frontier_erf.profile

         ERF_HOME=$(pwd) ./Build/cmake_with_kokkos_many_hip.sh
         cd install/bin
         # Next: choose a run mode in "Run After Build" below

      :download:`Download build+run snippet <scripts/quickstart/frontier_quickstart.sh>`
      :download:`Download Frontier sbatch example <scripts/quickstart/run_frontier_erf.sbatch>`
      More machine context: :ref:`sec:build:hpc`

   .. tab-item:: Aurora (ALCF)

      .. code-block:: bash

         git clone --recursive https://github.com/erf-model/ERF.git
         cd ERF
         source Build/machines/aurora_erf.profile
         export NETCDF_DIR=<path-to-netcdf>

         ERF_HOME=$(pwd) ./Build/cmake_with_kokkos_many_sycl.sh
         cd install/bin
         # Next: choose a run mode in "Run After Build" below

      :download:`Download build+run snippet <scripts/quickstart/aurora_quickstart.sh>`
      :download:`Download Aurora PBS example <scripts/quickstart/submit_erf_aurora.pbs>`
      Full guide: :doc:`aurora_build_run`

Not on a listed machine? See :ref:`sec:build:hpc` for machine profile customization,
or continue below for generic build workflows.

.. dropdown:: Run After Build: interactive smoke test or batch job
   :icon: rocket

   Use interactive mode for a quick check (2-10 steps). Use batch mode for normal
   production runs and longer jobs.
   Use interactive allocations primarily for debugging. Use batch jobs for normal
   testing and production runs.

   .. tab-set::

      .. tab-item:: Interactive Smoke Test (Fast)

         Slurm (Perlmutter/Kestrel/Frontier):

         .. code-block:: bash

            # Account/project is usually required on shared systems
            salloc -A <account_or_project> -N 1 -t 00:30:00
            srun -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most \
              amr.max_step=10

         PBS (Aurora):

         .. code-block:: bash

            qsub -I -A <PROJECT> -q debug -l select=1 -l walltime=1:00:00
            mpiexec -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most \
              amr.max_step=10

      .. tab-item:: Batch Scripts Provided

         .. list-table::
            :header-rows: 1

            * - System
              - Scheduler
              - Provided script
              - Launch command
            * - Perlmutter
              - Slurm
              - ``run_perlmutter_erf.sbatch``
              - ``sbatch ../../Docs/sphinx_doc/scripts/quickstart/run_perlmutter_erf.sbatch``
            * - Kestrel
              - Slurm
              - ``run.erf.aw.job_arena``
              - ``sbatch ../../Docs/sphinx_doc/scripts/quickstart/run.erf.aw.job_arena``
            * - Frontier
              - Slurm
              - ``run_frontier_erf.sbatch``
              - ``sbatch ../../Docs/sphinx_doc/scripts/quickstart/run_frontier_erf.sbatch``
            * - Aurora
              - PBS
              - ``submit_erf_aurora.pbs``
              - ``qsub ../../Docs/sphinx_doc/scripts/quickstart/submit_erf_aurora.pbs``

.. _sec:build:quickstart:gnumake:

Build with GNU Make
-------------------

For workstations or systems where CMake is not preferred:

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
         mpiexec -n 4 ./ERF3d.gnu.TPROF.MPI.CUDA.ex CanonicalTests/ABL/inputs_most

Default GNU Make builds are done in ``ERF/Exec``. Use ``ERF/.Exec_dev/<test_name>``
only when working on development tests.

.. _sec:build:quickstart:cmake:

Build with CMake
----------------

Suggested path for most users not on a listed HPC system. Choose other workflows only
when you need tighter control (see :ref:`sec:build:systems` for details).

.. tab-set::

   .. tab-item:: Automated Script (Suggested)

      .. code-block:: bash

         git clone --recursive git@github.com:erf-model/ERF.git
         cd ERF
         # Customize (optional): export ERF_BUILD_DIR=... ERF_SOURCE_DIR=... ERF_INSTALL_DIR=... or ERF_HOME=...
         ERF_HOME=$(pwd) ./Build/cmake_with_kokkos_many.sh

         # Run from install directory
         cd install/bin
         mpiexec -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most

      .. dropdown:: Optional: Customize build/source/install directories
         :icon: info

         Set ``ERF_BUILD_DIR``, ``ERF_SOURCE_DIR``, ``ERF_INSTALL_DIR``, or ``ERF_HOME``
         to override defaults.

         .. code-block:: bash

            ERF_HOME=$(pwd) ./Build/cmake_with_kokkos_many.sh

            ERF_BUILD_DIR=$(pwd)/build_custom \
            ERF_INSTALL_DIR=$(pwd)/install_custom \
            ./Build/cmake_with_kokkos_many.sh

   .. tab-item:: Out-of-Source Build

      .. code-block:: bash

         git clone --recursive git@github.com:erf-model/ERF.git
         cd ERF
         mkdir build && cd build
         ../Build/cmake.sh
         make -j

         # Run from build tree
         cd Exec
         mpiexec -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most

      .. dropdown:: Optional: Install and run from install directory
         :icon: info

         .. code-block:: bash

            make install
            cd ../../install/bin
            mpiexec -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most

   .. tab-item:: In Build/ Directory

      .. code-block:: bash

         git clone --recursive git@github.com:erf-model/ERF.git
         cd ERF/Build
         ./cmake.sh
         make -j

         # Run from Exec subdirectory
         cd Exec
         mpiexec -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most

Cleanup
-------

For a quick rebuild after source changes, ``make clean`` and rebuilding in place works
for most cases. If you have changed CMake configuration options or are hitting strange
errors, removing the build directory entirely and starting fresh is more reliable —
see :ref:`sec:build:troubleshooting`.

.. tab-set::

   .. tab-item:: GNU Make

      .. code-block:: bash

         # From ERF/Exec (or ERF/.Exec_dev/<test_name> for dev-test builds)
         make clean           # Remove build artifacts
         make realclean       # Same as clean
         make cleanconfig     # Remove configuration only

   .. tab-item:: In Build/ Directory

      .. code-block:: bash

         # From ERF/Build
         make distclean       # Clean all CMake artifacts
         make uninstall       # Uninstall based on install-manifest.txt

         # If configuration failed, manually remove
         rm -rf CMakeFiles/ CMakeCache.txt

   .. tab-item:: Out-of-Source / Automated Script

      .. code-block:: bash

         # From ERF root
         cmake --build build --target distclean
         cmake --build build --target uninstall

         # Or remove directories entirely
         rm -rf build/ install/

For detailed build options and troubleshooting, see :ref:`sec:build:systems`.
