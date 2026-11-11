#!/bin/bash
#
# ERF machine profile for NLR Kestrel GPU environment.
# Source this from the ERF repo root:
#   source Build/machines/kestrel_erf_cpu.profile
#
# This profile is intended for the "module profile + many-script" path:
#   ERF_HOME=$(pwd) ./Build/cmake_with_kokkos_many.sh
#
# Kestrel module names/versions can change; check current availability with:
#   module avail

ml reset 
ml cpe-stack 
module load PrgEnv-gnu
module load craype-x86-spr
module load cmake

module load cray-hdf5-parallel/1.14.3.5 cray-parallel-netcdf/1.12.3.17
module load cray-netcdf-hdf5parallel
module load cray-fftw/3.3.10.10
