#!/bin/bash

module load gcc/12.2.0
module load cray-mpich
module load cray-hdf5-parallel
module load cray-netcdf-hdf5parallel
module load cmake
module load cray-libsci
module load cray-parallel-netcdf
module tablelist

# Set environment for GPU-aware MPI with GTL
export MPICH_GPU_SUPPORT_ENABLED=1

# --- GPU ARCHITECTURE - Set for NVIDIA A100 ---
KOKKOS_GPU_ARCH="AMPERE80"  # A100 = Ampere 80
CMAKE_CUDA_ARCH="80"

#CC=$(which cc) CXX=$(which CC) FC=$(which ftn)

cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
      -DCMAKE_CUDA_STANDARD_LIBRARIES="-lmpi_gtl_cuda" \
      -DCMAKE_CXX_STANDARD_LIBRARIES="-lmpi_gtl_cuda" \
      -DCMAKE_CXX_FLAGS="$(CC --cray-print-opts=cflags)" \
      -DCMAKE_C_FLAGS="$(cc --cray-print-opts=cflags)" \
      -DCMAKE_Fortran_FLAGS="$(ftn --cray-print-opts=cflags)" \
      -DCMAKE_CUDA_FLAGS="$(CC --cray-print-opts=cflags)" \
      -DCMAKE_EXE_LINKER_FLAGS="$(CC --cray-print-opts=libs) $(cc --cray-print-opts=libs) $(ftn --cray-print-opts=libs)" \
      -DMPIEXEC_PREFLAGS:STRING=--oversubscribe \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DCMAKE_CXX_COMPILER:STRING=$(which CC) \
      -DCMAKE_C_COMPILER:STRING=$(which cc) \
      -DCMAKE_Fortran_COMPILER:STRING=$(which ftn) \
      -DERF_DIM:STRING=3 \
      -DERF_ENABLE_FFT:BOOL=ON \
      -DERF_ENABLE_NETCDF:BOOL=ON \
      -DERF_ENABLE_KOKKOS:BOOL=ON \
      -DERF_ENABLE_RRTMGP:BOOL=ON \
      -DERF_ENABLE_EKAT:BOOL=ON \
      -DERF_ENABLE_MPI:BOOL=ON \
      -DERF_ENABLE_CUDA:BOOL=ON \
      -DAMReX_CUDA_ARCH=8.0 \
      -DERF_ENABLE_TESTS:BOOL=ON \
      -DERF_ENABLE_FCOMPARE:BOOL=ON \
      -DERF_ENABLE_DOCUMENTATION:BOOL=OFF \
      -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
      -B build_kokkos_ekat_gtl_sci ..

cmake --build build_kokkos_ekat_gtl_sci -j10 -v
cmake --install build_kokkos_ekat_gtl_sci
