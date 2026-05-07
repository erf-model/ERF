#!/bin/bash

# --- Log file setup ---
LOG_FILE="cmake_build_$(date '+%Y%m%d_%H%M%S').log"
exec > >(tee -a "${LOG_FILE}") 2>&1
echo "Logging output to: ${LOG_FILE}"
echo "Build started at: $(date)"
echo "==========================================="

# --- Clean the build directory first ---
rm -rf CMakeCache.txt CMakeFiles/ install

# --- Define paths ---
HDF5_INSTALL_DIR="/nopt/nrel/apps/cpu_stack/libraries-craympich/hdf5/hdf5-1.14.6"
NETCDF_INSTALL_DIR="/nopt/nrel/apps/cpu_stack/libraries-craympich/netcdf/netcdf-4.9.3"
CRAY_MPI_BIN_DIR="/opt/cray/pe/mpich/8.1.28/ofi/gnu/10.3/bin"

NETCDF_DIR="${NETCDF_INSTALL_DIR}"
NETCDF_INCLUDE_PATH="${NETCDF_INSTALL_DIR}/include"

# --- Runtime library paths ---
export LD_LIBRARY_PATH="${HDF5_INSTALL_DIR}/lib:${LD_LIBRARY_PATH}"
export CMAKE_LIBRARY_PATH="${HDF5_INSTALL_DIR}/lib:${CMAKE_LIBRARY_PATH}"

# --- GPU architecture for NVIDIA H100 ---
KOKKOS_GPU_ARCH="HOPPER90"
CMAKE_CUDA_ARCH="90"

# --- CUDA Toolkit ---
CUDA_ROOT="/nopt/cuda/12.9"
export PATH="${CUDA_ROOT}/bin:${PATH}"
export LD_LIBRARY_PATH="${CUDA_ROOT}/lib64:${LD_LIBRARY_PATH}"

# --- GPU-aware MPI with GTL ---
export MPICH_GPU_SUPPORT_ENABLED=1

GTL_LIBS="-L/opt/cray/pe/mpich/8.1.28/gtl/lib -lmpi_gtl_cuda"
GTL_LIBS="${GTL_LIBS} -L/opt/cray/pe/mpich/8.1.28/ofi/gnu/10.3/lib -lmpi_gnu_103"
GTL_LIBS="${GTL_LIBS} -L/nopt/cuda/12.9/lib64/stubs -lcuda"

CMAKE_ARGS=(
    "-DCMAKE_INSTALL_PREFIX:PATH=./install"
    "-DCMAKE_CUDA_STANDARD_LIBRARIES=${GTL_LIBS}"
    "-DCMAKE_CXX_STANDARD_LIBRARIES=${GTL_LIBS}"
    "-DCMAKE_C_COMPILER=${CRAY_MPI_BIN_DIR}/mpicc"
    "-DCMAKE_CXX_COMPILER=${CRAY_MPI_BIN_DIR}/mpicxx"
    "-DCMAKE_Fortran_COMPILER=${CRAY_MPI_BIN_DIR}/mpifort"
    "-DCMAKE_CUDA_HOST_COMPILER=${CRAY_MPI_BIN_DIR}/mpicxx"
    "-DMPIEXEC_PREFLAGS:STRING=--oversubscribe"
    "-DCMAKE_BUILD_TYPE:STRING=Release"
    "-DERF_DIM:STRING=3"
    "-DERF_ENABLE_MPI:BOOL=ON"
    "-DERF_ENABLE_CUDA:BOOL=ON"
    "-DERF_ENABLE_HDF5:BOOL=ON"
    "-DERF_ENABLE_NETCDF:BOOL=ON"
    "-DERF_ENABLE_FFT:BOOL=ON"
    "-DERF_ENABLE_KOKKOS:BOOL=ON"
    "-DERF_ENABLE_RRTMGP:BOOL=ON"
    "-DERF_ENABLE_NOAHMP:BOOL=ON"
    "-DERF_ENABLE_P3=ON"
    "-DERF_ENABLE_SHOC=ON"
    "-DERF_ENABLE_EKAT:BOOL=ON"
    "-DERF_ENABLE_TESTS:BOOL=ON"
    "-DERF_ENABLE_FCOMPARE:BOOL=ON"
    "-DERF_ENABLE_DOCUMENTATION:BOOL=OFF"
    "-DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON"
    "-DHAVE_MPICH=TRUE"
    "-DERF_CHECK_MODULES=OFF"
    "-DMPICH_SKIP_MPICXX=TRUE"
    "-DCMAKE_CUDA_ARCHITECTURES=${CMAKE_CUDA_ARCH}"
    "-DCUDAToolkit_ROOT=${CUDA_ROOT}"
    "-DKokkos_ARCH_${KOKKOS_GPU_ARCH}=ON"
    "-DCMAKE_CXX_FLAGS=-I${NETCDF_INCLUDE_PATH}"
    "-DHDF5_ROOT=${HDF5_INSTALL_DIR}"
    "-DHDF5_DIR=${HDF5_INSTALL_DIR}"
    "-DHDF5_INCLUDE_DIR=${HDF5_INSTALL_DIR}/include"
    "-DHDF5_LIBRARY_DIR=${HDF5_INSTALL_DIR}/lib"
    "-DHDF5_HL_LIBRARY=${HDF5_INSTALL_DIR}/lib/libhdf5_hl.so"
    "-DHDF5_C_LIBRARY=${HDF5_INSTALL_DIR}/lib/libhdf5.so"
    "-DCMAKE_INSTALL_RPATH=${HDF5_INSTALL_DIR}/lib"
    "-DCMAKE_BUILD_WITH_INSTALL_RPATH=ON"
    "-DNETCDF_DIR=${NETCDF_DIR}"
    "-DCMAKE_PREFIX_PATH=${HDF5_INSTALL_DIR};${NETCDF_INSTALL_DIR}"
    ".."
)

echo "Running CMake with the following arguments:"
printf "  %s\n" "${CMAKE_ARGS[@]}"
echo "-------------------------------------------"

cmake "${CMAKE_ARGS[@]}" && make -j8
BUILD_STATUS=$?

echo "==========================================="
if [ ${BUILD_STATUS} -eq 0 ]; then
    echo "Build SUCCEEDED at: $(date)"
else
    echo "Build FAILED at: $(date) (exit code: ${BUILD_STATUS})"
fi
echo "Full log saved to: ${LOG_FILE}"

exit ${BUILD_STATUS}
