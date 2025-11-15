.. _sec:build:quickstart:


Quickstart clone-build-run-clean instructions
=============================================

Build with GNUMake
------------------

* Clone the repository with submodules and navigate to the ABL execution directory
* Compile using GNU compiler with MPI support (optionally with CUDA for GPU builds)
* Run the executable with 4 MPI processes using the provided input file

.. tab-set::

   .. tab-item:: CPU Build

      Clone, build with GNU/MPI, and run on CPU.

      .. code:: shell

          git clone --recursive git@github.com/erf-model/ERF # note some submodules require / expect ssh keys
          cd ERF/Exec/ABL
          make COMP=gnu USE_MPI=TRUE
          mpiexec -n 4 ./ERF3d.gnu.TPROF.MPI.ex inputs_most

   .. tab-item:: GPU Build

      Clone, build with GNU/MPI/CUDA, and run on GPU.

      .. code:: shell

          git clone --recursive git@github.com/erf-model/ERF # note some submodules require / expect ssh keys
          cd ERF/Exec/ABL
          make COMP=gnu USE_MPI=TRUE USE_CUDA=TRUE
          mpiexec -n 4 ./ERF3d.gnu.TPROF.MPI.CUDA.ex inputs_most

Build with CMake
----------------

* Clone the repository and use the provided CMake configuration scripts in the Build directory
* Build the project and optionally install to a dedicated directory for cleaner organization
* Execute the compiled binary from either the build or install location with MPI

.. tab-set::

   .. tab-item:: CPU Build

      Clone, configure with CMake script, build, and run on CPU.

      .. code:: shell

          git clone --recursive git@github.com/erf-model/ERF # note some submodules require / expect ssh keys
          cd Build
          ./cmake.sh

          cd Exec/ABL
          mpiexec -n 4 ./erf_abl ../../../Exec/ABL/inputs_most
          # Extra step for install directory:
          # make install
          # cd install/bin
          # mpiexec -n 4 ./erf_abl ../../../Exec/ABL/inputs_most

   .. tab-item:: GPU Build

      Clone, configure with CUDA CMake script, install, and run on GPU.

      .. code:: shell

          git clone --recursive git@github.com/erf-model/ERF # note some submodules require / expect ssh keys
          cd Build
          ./cmake_cuda.sh
          make install
          cd install/bin
          mpiexec -n 4 ./erf_abl inputs_most

   .. tab-item:: CPU with Kokkos

      Clone, configure with Kokkos-enabled CMake script, install, and run on CPU.

      .. code:: shell

          git clone --recursive git@github.com/erf-model/ERF # note some submodules require / expect ssh keys
          cd Build
          ./cmake_with_kokkos_many.sh

          cd install_erf/bin
          mpiexec -n 4 ./erf_abl ../../../Exec/ABL/inputs_most

Cleanup Commands
----------------

* GNUMake clean commands remove build artifacts from the current directory
* CMake clean commands operate from the Build directory and can uninstall or clean build artifacts
* Use distclean for complete removal of all CMake-generated files and build products

.. tab-set::

   .. tab-item:: GNUMake Clean

      Remove build artifacts using GNUMake clean targets in the current directory.

      .. code:: shell

          make clean # Standard
          make USE_MPI=TRUE cleanconfig # Just the GNUMakefile and command line flags cleanup
          make realclean # Same as clean

   .. tab-item:: CMake Clean from Build

      Uninstall or clean CMake build from the ERF/Build directory.

      .. code:: shell

          cd ../../
          cd ERF/Build
          make uninstall # Uninstall based on install-manifest.txt
          make distclean # Clean up anything cmake configured or used for building

   .. tab-item:: CMake Clean from build_erf

      Uninstall or clean CMake build using the build_erf subdirectory.

      .. code:: shell

          cd ../../
          cd ERF/Build
          cmake --build build_erf --target uninstall # Uninstall based on install-manifest.txt
          cmake --build build_erf --target distclean # Clean up anything cmake configured or used for building

   .. tab-item:: CMake Clean from build

      Uninstall or clean CMake build using the build subdirectory.

      .. code:: shell

          cd ../../
          cd ERF/Build
          cmake --build build --target uninstall # Uninstall based on install-manifest.txt
          cmake --build build --target distclean # Clean up anything cmake configured or used for building
