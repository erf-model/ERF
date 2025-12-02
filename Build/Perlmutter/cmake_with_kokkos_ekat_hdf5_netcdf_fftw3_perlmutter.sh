#!/bin/bash

module load gcc-native/13.2
module load cray-mpich/8.1.30
module load cray-hdf5-parallel/1.14.3.1
module load cray-netcdf-hdf5parallel/4.9.0.13
module load cmake/3.30.2
module load cray-libsci/24.07.0
module load cray-parallel-netcdf/1.12.3.13
module tablelist

# Set environment for GPU-aware MPI with GTL
export MPICH_GPU_SUPPORT_ENABLED=1

# --- GPU ARCHITECTURE - Set for NVIDIA A100 ---
KOKKOS_GPU_ARCH="AMPERE80"  # A100 = Ampere 80
CMAKE_CUDA_ARCH="80"

CC=$(which cc) CXX=$(which CC) FC=$(which ftn)

CRAY_LIBS_CLEAN=$(CC --cray-print-opts=libs | sed 's/-Wl,--as-needed,//g; s/,--no-as-needed//g; s/,-l/ -l/g')
CRAY_LIBS_CLEAN="$CRAY_LIBS_CLEAN $(cc --cray-print-opts=libs | sed 's/-Wl,--as-needed,//g; s/,--no-as-needed//g; s/,-l/ -l/g')"
CRAY_LIBS_CLEAN="$CRAY_LIBS_CLEAN $(ftn --cray-print-opts=libs | sed 's/-Wl,--as-needed,//g; s/,--no-as-needed//g; s/,-l/ -l/g')"

echo $CRAY_LIBS_CLEAN

cmake -DCMAKE_INSTALL_PREFIX:PATH=./install_erf \
      -DCMAKE_CUDA_STANDARD_LIBRARIES="-lmpi_gnu_123 -lmpi_gtl_cuda" \
      -DCMAKE_CXX_STANDARD_LIBRARIES="-lmpi_gnu_123 -lmpi_gtl_cuda" \
      -DCMAKE_CXX_FLAGS="$(CC --cray-print-opts=cflags)" \
      -DCMAKE_C_FLAGS="$(cc --cray-print-opts=cflags)" \
      -DCMAKE_Fortran_FLAGS="$(ftn --cray-print-opts=cflags)" \
      -DCMAKE_CUDA_FLAGS="$(CC --cray-print-opts=cflags)" \
      -DCMAKE_EXE_LINKER_FLAGS="-Wl,--no-as-needed $CRAY_LIBS_CLEAN" \
      -DMPIEXEC_PREFLAGS:STRING=--oversubscribe \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DCMAKE_CXX_COMPILER:STRING=$(which CC) \
      -DCMAKE_C_COMPILER:STRING=$(which cc) \
      -DCMAKE_Fortran_COMPILER:STRING=$(which ftn) \
      -DERF_DIM:STRING=3 \
      -DERF_ENABLE_FFT:BOOL=ON \
      -DERF_ENABLE_NETCDF:BOOL=ON \
      -DERF_ENABLE_RRTMGP:BOOL=ON \
      -DERF_ENABLE_MPI:BOOL=ON \
      -DERF_ENABLE_CUDA:BOOL=ON \
      -DAMReX_CUDA_ARCH=8.0 \
      -DERF_ENABLE_TESTS:BOOL=ON \
      -DERF_ENABLE_FCOMPARE:BOOL=ON \
      -DERF_ENABLE_DOCUMENTATION:BOOL=OFF \
      -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
      -B build_erf ..

cmake --build build_erf -j10 -v
cmake --install build_erf --prefix=install_erf
