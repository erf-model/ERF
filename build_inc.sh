#!/bin/bash
set -euo pipefail

# Keep the compiler, MPI, math, I/O, and CMake stack identical across
# incremental builds. Swap site defaults to the pinned Perlmutter stack.
export LMOD_PAGER=cat
module swap gcc-native/14 gcc-native/13.2 2>/dev/null || module load gcc-native/13.2
module swap cray-mpich/9.0.1 cray-mpich/8.1.30 2>/dev/null || module load cray-mpich/8.1.30
module swap cray-libsci/25.09.0 cray-libsci/24.07.0 2>/dev/null || module load cray-libsci/24.07.0
module load cray-hdf5-parallel/1.14.3.1
module load cray-netcdf-hdf5parallel/4.9.0.13
module load cmake/3.30.2
module load cray-parallel-netcdf/1.12.3.13

for required_module in gcc-native/13.2 cray-mpich/8.1.30 cray-libsci/24.07.0 cray-hdf5-parallel/1.14.3.1 cray-netcdf-hdf5parallel/4.9.0.13 cmake/3.30.2 cray-parallel-netcdf/1.12.3.13; do
    module is-loaded "$required_module" || { echo "Required module is not loaded: $required_module" >&2; exit 1; }
done
module list

# Set environment for GPU-aware MPI with GTL
export MPICH_GPU_SUPPORT_ENABLED=1

# --- GPU ARCHITECTURE - Set for NVIDIA A100 ---
KOKKOS_GPU_ARCH="AMPERE80"  # A100 = Ampere 80
CMAKE_CUDA_ARCH="80"

CC=$(which cc) CXX=$(which CC) FC=$(which ftn)

CRAY_LIBS_CLEAN=$(CC --cray-print-opts=libs | sed 's/-Wl,--as-needed,//g; s/,--no-as-needed//g; s/,-l/ -l/g')
CRAY_LIBS_CLEAN="$CRAY_LIBS_CLEAN $(cc --cray-print-opts=libs | sed 's/-Wl,--as-needed,//g; s/,--no-as-needed//g; s/,-l/ -l/g')"
CRAY_LIBS_CLEAN="$CRAY_LIBS_CLEAN $(ftn --cray-print-opts=libs | sed 's/-Wl,--as-needed,//g; s/,--no-as-needed//g; s/,-l/ -l/g')"

echo "$CRAY_LIBS_CLEAN"

BUILD_DIR="${ERF_BUILD_DIR:-build_real_gpu}"
INSTALL_DIR="${ERF_INSTALL_DIR:-install_real_gpu}"
PRECISION="${ERF_PRECISION:-DOUBLE}"

# --- Only run cmake configure if the build dir doesn't exist yet, or if
#     the caller explicitly forces a reconfigure with --reconfigure ---
FORCE_RECONFIGURE=0

for arg do
    case "$arg" in
        --single|--precision=SINGLE)
            BUILD_DIR="build_single_gpu"
            INSTALL_DIR="install_single_gpu"
            PRECISION="SINGLE"
            ;;
        --double|--precision=DOUBLE)
            BUILD_DIR="build_real_gpu"
            INSTALL_DIR="install_real_gpu"
            PRECISION="DOUBLE"
            ;;
        --reconfigure)
            FORCE_RECONFIGURE=1
            ;;
        *)
            echo "Unknown option: $arg" >&2
            echo "Usage: $0 [--single|--double] [--reconfigure]" >&2
            exit 2
            ;;
    esac
done

echo ">>> Using pinned modules for $PRECISION precision in $BUILD_DIR <<<"

if [[ ! -d "$BUILD_DIR" || ! -f "$BUILD_DIR/CMakeCache.txt" || $FORCE_RECONFIGURE -eq 1 ]]; then
    echo ">>> Running CMake configure step (fresh or forced reconfigure) <<<"
    cmake -U "HDF5_*" -U "MPI_*" -U "NETCDF_*" -U "NetCDF_*" \
          -DCMAKE_INSTALL_PREFIX:PATH=./${INSTALL_DIR} \
          -DCMAKE_CUDA_STANDARD_LIBRARIES="-lmpi_gnu_123 -lmpi_gtl_cuda" \
          -DCMAKE_CXX_STANDARD_LIBRARIES="-lmpi_gnu_123 -lmpi_gtl_cuda" \
          -DCMAKE_CXX_FLAGS="$(CC --cray-print-opts=cflags)" \
          -DCMAKE_C_FLAGS="$(cc --cray-print-opts=cflags)" \
          -DCMAKE_Fortran_FLAGS="$(ftn --cray-print-opts=cflags)" \
          -DCMAKE_CUDA_FLAGS="$(CC --cray-print-opts=cflags)" \
          -DCMAKE_EXE_LINKER_FLAGS="-Wl,--no-as-needed $CRAY_LIBS_CLEAN" \
          -DMPIEXEC_PREFLAGS:STRING=--oversubscribe \
          -DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo \
          -DCMAKE_CXX_COMPILER:STRING=$(which CC) \
          -DCMAKE_C_COMPILER:STRING=$(which cc) \
          -DCMAKE_Fortran_COMPILER:STRING=$(which ftn) \
          -DERF_DIM:STRING=3 \
          -DERF_ENABLE_FFT:BOOL=ON \
          -DERF_ENABLE_NETCDF:BOOL=ON \
          -DERF_ENABLE_RRTMGP:BOOL=ON \
          -DERF_ENABLE_NOAHMP:BOOL=ON \
          -DERF_ENABLE_SHOC:BOOL=OFF \
          -DERF_ENABLE_P3:BOOL=OFF \
          -DERF_ENABLE_MPI:BOOL=ON \
          -DERF_ENABLE_CUDA:BOOL=ON \
          -DAMReX_CUDA_ARCH=8.0 \
          -DERF_PRECISION:STRING="$PRECISION" \
          -DERF_ENABLE_TESTS:BOOL=ON \
          -DERF_ENABLE_FCOMPARE:BOOL=ON \
          -DERF_ENABLE_DOCUMENTATION:BOOL=OFF \
          -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
          -B "$BUILD_DIR" ..
else
    echo ">>> Existing build dir with CMakeCache.txt found — skipping configure, doing incremental build <<<"
fi

cmake --build "$BUILD_DIR" -j10 -v
cmake --install "$BUILD_DIR" --prefix="$INSTALL_DIR"
