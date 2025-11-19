#!/bin/bash

module load gcc-native/13.2 cmake cudatoolkit cray-hdf5-parallel cray-netcdf-hdf5parallel cray-libsci

#module load gcc-native/13.2
#module load cray-mpich/8.1.30
#module load cray-hdf5-parallel/1.14.3.1
#module load cray-netcdf-hdf5parallel/4.9.0.13
#module load cmake/3.30.2
#module load cray-libsci/24.07.0
#module load cray-parallel-netcdf/1.12.3.13

# Automatically included with module load gpu
# export MPICH_GPU_SUPPORT_ENABLED=1