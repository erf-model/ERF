.. _sec:build:quickstart:


Quickstart clone-build-run-clean instructions
=============================================

.. tab-set::

   .. tab-item:: Quickstart GNUMake CPU

      .. code:: shell

          git clone --recursive git@github.com/erf-model/ERF # note some submodules require / expect ssh keys
          cd ERF/Exec/ABL
          make COMP=gnu USE_MPI=TRUE
          mpiexec -n 4 ./ERF3d.gnu.TPROF.MPI.ex inputs_most

   .. tab-item:: Quickstart GNUMake GPU

      .. code:: shell

          git clone --recursive git@github.com/erf-model/ERF # note some submodules require / expect ssh keys
          cd ERF/Exec/ABL
          make COMP=gnu USE_MPI=TRUE USE_CUDA=TRUE
          mpiexec -n 4 ./ERF3d.gnu.TPROF.MPI.CUDA.ex inputs_most

.. tab-set::

   .. tab-item:: Quickstart CMake CPU

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

   .. tab-item:: Quickstart CMake GPU

      .. code:: shell

          git clone --recursive git@github.com/erf-model/ERF # note some submodules require / expect ssh keys
          cd Build
          ./cmake_cuda.sh
          make install
          cd install/bin
          mpiexec -n 4 ./erf_abl inputs_most

   .. tab-item:: Quickstart CMake CPU Many

      .. code:: shell

          git clone --recursive git@github.com/erf-model/ERF # note some submodules require / expect ssh keys
          cd Build
          ./cmake_with_kokkos_many.sh

          cd install_erf/bin
          mpiexec -n 4 ./erf_abl ../../../Exec/ABL/inputs_most

.. tab-set::

   .. tab-item:: GNUMake same dir

      .. code:: shell

          make clean # Standard
          make USE_MPI=TRUE cleanconfig # Just the GNUMakefile and command line flags cleanup
          make realclean # Same as clean

   .. tab-item:: ERF/Build dir

      .. code:: shell

          cd ../../
          cd ERF/Build
          make uninstall # Uninstall based on install-manifest.txt
          make distclean # Clean up anything cmake configured or used for building

   .. tab-item:: ERF/Build/build_erf dir

      .. code:: shell

          cd ../../
          cd ERF/Build
          cmake --build build_erf --target uninstall # Uninstall based on install-manifest.txt
          cmake --build build_erf --target distclean # Clean up anything cmake configured or used for building

   .. tab-item:: ERF/build dir

      .. code:: shell

          cd ../../
          cd ERF/Build
          cmake --build build --target uninstall # Uninstall based on install-manifest.txt
          cmake --build build --target distclean # Clean up anything cmake configured or used for building
