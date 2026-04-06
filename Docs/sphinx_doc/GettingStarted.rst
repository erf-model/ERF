.. _GettingStarted:

Getting Started
===============

Clone, build, and run ERF in a few steps:

.. code-block:: bash

   git clone --recursive https://github.com/erf-model/ERF.git
   cd ERF
   # Customize (optional): export ERF_BUILD_DIR=... ERF_SOURCE_DIR=... ERF_INSTALL_DIR=... or ERF_HOME=...
   ERF_HOME=$(pwd) ./Build/cmake_with_kokkos_many.sh

   # Run from install directory
   cd install/bin
   mpiexec -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most

.. dropdown:: On an HPC system? Start here instead.
   :icon: rocket

   Each tab sources a machine profile that handles environment setup, then
   runs the appropriate build script for that machine's GPU backend. For
   download links and full job submission scripts, see :ref:`sec:build:quickstart:hpc-cmake`.

   .. tab-set::

      .. tab-item:: Perlmutter (NERSC)

         .. code-block:: bash

            git clone --recursive https://github.com/erf-model/ERF.git
            cd ERF
            source Build/machines/perlmutter_erf.profile
            ERF_HOME=$(pwd) ./Build/cmake_with_kokkos_many_cuda.sh
            cd install/bin
            sbatch ../../Docs/sphinx_doc/scripts/quickstart/run_perlmutter_erf.sbatch

      .. tab-item:: Kestrel (NREL)

         .. code-block:: bash

            git clone --recursive https://github.com/erf-model/ERF.git
            cd ERF
            source Build/machines/kestrel_erf.profile
            ERF_HOME=$(pwd) ./Build/cmake_with_kokkos_many_cuda.sh
            cd install/bin
            sbatch ../../Docs/sphinx_doc/scripts/quickstart/run.erf.aw.job_arena

      .. tab-item:: Frontier (OLCF)

         .. code-block:: bash

            git clone --recursive https://github.com/erf-model/ERF.git
            cd ERF
            source Build/machines/frontier_erf.profile
            ERF_HOME=$(pwd) ./Build/cmake_with_kokkos_many_hip.sh
            cd install/bin
            sbatch ../../Docs/sphinx_doc/scripts/quickstart/run_frontier_erf.sbatch

      .. tab-item:: Aurora (ALCF)

         .. code-block:: bash

            git clone --recursive https://github.com/erf-model/ERF.git
            cd ERF
            source Build/machines/aurora_erf.profile
            export NETCDF_DIR=<path-to-netcdf>
            ERF_HOME=$(pwd) ./Build/cmake_with_kokkos_many_sycl.sh
            cd install/bin
            qsub ../../Docs/sphinx_doc/scripts/quickstart/submit_erf_aurora.pbs

   Not on a listed machine? See :ref:`sec:build:hpc` for machine profile
   customization.

----

Need more detail? Choose the path that fits your situation:

- **More copy-paste options** — GNU Make, alternative CMake workflows, and
  full HPC scripts with download links → :doc:`quickstart`
- **Step-by-step walkthrough** of what each stage is doing and why — work
  through :doc:`submodule`, :doc:`building`, :doc:`InputFiles`, and
  :doc:`testing` in order
- **Machine or library specific setup**, advanced configuration, or
  troubleshooting → :ref:`building:configuration`
  (*Beyond the Basics: Machines and Libraries*)

.. toctree::
   :maxdepth: 1
   :hidden:

   quickstart
   submodule
   building
   InputFiles
   testing
