#!/bin/bash

# Use the default cudatoolkit and cray-mpich: the default cray-mpich's GPU transport
# layer (libmpi_gtl_cuda) needs the CUDA runtime major version of the default
# cudatoolkit, and the bundled Kokkos builds against it. Pinning one without the
# other breaks the build or the run (see Docs/sphinx_doc/buildingTroubleshooting.rst).
# cray-netcdf-hdf5parallel/4.9.2.3 is the first netcdf built against the
# cray-mpich 9.1 library name (libmpi_gnu.so.12); the older default links
# libmpi_gnu_123.so.12, which loads cray-mpich 9.0.1 alongside 9.1.0.
module load cmake cray-hdf5-parallel cray-netcdf-hdf5parallel/4.9.2.3 cray-libsci

#module load gcc-native/13.2
#module load cray-mpich/8.1.30
#module load cray-hdf5-parallel/1.14.3.1
#module load cray-netcdf-hdf5parallel/4.9.0.13
#module load cmake/3.30.2
#module load cray-libsci/24.07.0
#module load cray-parallel-netcdf/1.12.3.13

# Automatically included with module load gpu
# export MPICH_GPU_SUPPORT_ENABLED=1
