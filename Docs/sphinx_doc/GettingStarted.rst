.. role:: cpp(code)
   :language: c++

.. _GettingStarted:

Getting Started
===============

Start with :doc:`quickstart`, then choose the path that matches your setup:

1. **Automated script (suggested):** use ``Build with CMake -> Automated Script (Suggested)``. This is the default path for most users.
   
   .. code-block:: bash

      git clone --recursive git@github.com:erf-model/ERF.git
      cd ERF
      # source Build/machines/<my_machine>_erf.profile
      ./Build/cmake_with_kokkos_many.sh
      cd install/bin
      mpiexec -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most

   .. dropdown:: HPC machine-specific quickstart variants
      :icon: rocket

      For machine-profile-driven commands, see :ref:`sec:build:quickstart:hpc-cmake` in :doc:`quickstart`.

      * :doc:`perlmutter_build_run`
      * :doc:`kestrel_build_run`
      * :doc:`aurora_build_run`

2. **Clone-build-run single page:** stay in :doc:`quickstart` and use the copy-paste commands directly.
3. **More detailed step-by-step sections:** continue through :doc:`submodule`, :doc:`building`, :doc:`InputFiles`, and :doc:`testing`.

More detailed machine- and library-specific information is in :ref:`building:configuration` (Beyond the Basics: Machines and Libraries).

.. toctree::
   :maxdepth: 1

   quickstart
   submodule
   building
   InputFiles
   testing
