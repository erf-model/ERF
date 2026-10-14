#!/bin/bash
#
# ERF machine profile for NLR Kestrel GPU environment.
# Source this from the ERF repo root:
#   source Build/machines/kestrel_erf_gpu.profile
#
# This profile is intended for the "module profile + many-script" path:
#   ERF_HOME=$(pwd) ./Build/cmake_with_kokkos_many_cuda.sh
#
# Kestrel module names/versions can change; check current availability with:
#   module avail

module reset 
module load cpe-stack
module load PrgEnv-gnu
module load craype-x86-genoa
module load cmake
module load cudatoolkit
module load craype-accel-nvidia90
export MPICH_GPU_SUPPORT_ENABLED=1 

module load cray-hdf5-parallel/1.14.3.5 cray-parallel-netcdf/1.12.3.17
module load cray-netcdf-hdf5parallel
module load cray-fftw/3.3.10.10
