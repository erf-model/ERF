#!/bin/bash
#
# ERF machine profile for NREL Kestrel GPU environment.
# Source this from the ERF repo root:
#   source Build/machines/kestrel_erf.profile
#
# This profile is intended for the "module profile + many-script" path:
#   ERF_HOME=$(pwd) ./Build/cmake_with_kokkos_many_cuda.sh
#
# Kestrel module names/versions can change; check current availability with:
#   module avail

# module purge
module load PrgEnv-gnu/8.5.0
module load craype-x86-milan
module load cuda/12.9
module load cmake
#module load craype-accel-nvidia90
export CRAY_ACCEL_TARGET=nvidia90

# Prefer NREL /nopt library stack for ERF (provides netcdf_par.h).
export HDF5_ROOT=/nopt/nrel/apps/cpu_stack/libraries-craympich/hdf5/hdf5-1.14.6
export HDF5_DIR="${HDF5_ROOT}"
export NETCDF_DIR=/nopt/nrel/apps/cpu_stack/libraries-craympich/netcdf/netcdf-4.9.3
export CMAKE_PREFIX_PATH="${HDF5_ROOT}:${NETCDF_DIR}:${CMAKE_PREFIX_PATH:-}"

# Enable GPU-aware MPI when available.
export MPICH_GPU_SUPPORT_ENABLED=1
