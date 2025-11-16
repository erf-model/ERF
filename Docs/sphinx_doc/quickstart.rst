.. _sec:build:quickstart:

Quickstart: Clone-Build-Run
============================

Copy-paste commands for common build scenarios.

Build with GNU Make
-------------------

Clone, build, and run with GNU Make:

.. tab-set::

   .. tab-item:: CPU Build

      .. code-block:: bash

         git clone --recursive git@github.com:erf-model/ERF.git
         cd ERF/Exec/ABL
         make COMP=gnu USE_MPI=TRUE
         mpiexec -n 4 ./ERF3d.gnu.TPROF.MPI.ex inputs_most

   .. tab-item:: GPU Build

      .. code-block:: bash

         git clone --recursive git@github.com:erf-model/ERF.git
         cd ERF/Exec/ABL
         make COMP=gnu USE_MPI=TRUE USE_CUDA=TRUE
         mpiexec -n 4 ./ERF3d.gnu.TPROF.MPI.CUDA.ex inputs_most

Build with CMake
----------------

Clone, build, and run with CMake. Choose workflow based on preference (see :ref:`sec:build:systems` for details):

.. tab-set::

   .. tab-item:: In Build/ Directory

      Build in Build/ directory (similar to GNU Make structure):

      .. code-block:: bash

         git clone --recursive git@github.com:erf-model/ERF.git
         cd ERF/Build
         ./cmake.sh
         make -j

         # Run from Exec subdirectory
         cd Exec/ABL
         mpiexec -n 4 ./erf_abl ../../../Exec/ABL/inputs_most

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
         cd Exec/ABL
         mpiexec -n 4 ./erf_abl ../../../Exec/ABL/inputs_most

         # Or run from install directory
         cd ../../install/bin
         mpiexec -n 4 ./erf_abl ../../../Exec/ABL/inputs_most

.. tab-item:: Automated Script

      Auto-creates build and install directories. Customize with environment variables or edit script defaults.

      **Customize directories (optional):**

      Set ``ERF_BUILD_DIR``, ``ERF_SOURCE_DIR``, ``ERF_INSTALL_DIR`` environment variables.

      .. code-block:: bash

         git clone --recursive git@github.com:erf-model/ERF.git
         cd ERF

         # Set directory with env var
         ERF_HOME=$(pwd) ./Build/cmake_with_kokkos_many.sh

         # Use defaults (build in current dir, source from .., install to install/)
         # mkdir build; cd build; ../Build/cmake_with_kokkos_many.sh; cd ../

         # Run from install directory
         cd install/bin  # or $ERF_INSTALL_DIR/bin if customized
         mpiexec -n 4 ./erf_abl ../../Exec/ABL/inputs_most

Cleanup
-------

Remove build artifacts:

.. tab-set::

   .. tab-item:: GNU Make

      .. code-block:: bash

         # From problem directory (e.g., Exec/ABL)
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
