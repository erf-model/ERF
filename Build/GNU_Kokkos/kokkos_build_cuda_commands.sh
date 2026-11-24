#!/bin/bash

# Build and install the Kokkos submodule for a CUDA build of ERF with RRTMGP.
#
# This is optional: the GNU Make build builds Kokkos itself under
# Submodules/kokkos/build/<config> when ERF_KOKKOS_HOME is not set, and the CMake
# build adds the submodule directly. Use this script to make one install that
# several builds share, then point them at it with ERF_KOKKOS_HOME (GNU Make) or
# Kokkos_ROOT (CMake).
#
# Set KOKKOS_ARCH to the Kokkos name of the GPU (e.g. AMPERE80 for A100,
# HOPPER90 for H100); without it Kokkos autodetects the GPU of the build node.

# Set path to ERF
: "${ERF_DIR:=/home/ERF}"

# The Kokkos source, and the build and install directories (not in the source)
SRC_DIR=${ERF_DIR}/Submodules/kokkos
WORK_DIR=${SRC_DIR}/build/standalone-cuda
INSTALL_DIR=${WORK_DIR}/install

# nvcc_wrapper's host compiler: the MPI C++ wrapper, never nvcc itself
export NVCC_WRAPPER_DEFAULT_COMPILER=${NVCC_WRAPPER_DEFAULT_COMPILER:-$(command -v mpicxx)}

ARCH_FLAG=""
if [ -n "${KOKKOS_ARCH}" ]; then
  ARCH_FLAG="-D Kokkos_ARCH_${KOKKOS_ARCH}:BOOL=ON"
fi

rm -rf ${WORK_DIR}/CMakeFiles ${WORK_DIR}/CMakeCache.txt

# These match the options Exec/Make.ERF builds Kokkos with
cmake -S ${SRC_DIR} -B ${WORK_DIR}                                 \
  -D CMAKE_BUILD_TYPE:STRING=Release                              \
  -D CMAKE_CXX_COMPILER:STRING=${SRC_DIR}/bin/nvcc_wrapper        \
  -D CMAKE_CXX_STANDARD:STRING=20                                 \
  -D CMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR}                     \
  -D CMAKE_INSTALL_LIBDIR:PATH=lib                                \
  -D Kokkos_ENABLE_SERIAL:BOOL=ON                                 \
  -D Kokkos_ENABLE_CUDA:BOOL=ON                                   \
  -D Kokkos_ENABLE_CUDA_LAMBDA:BOOL=ON                            \
  -D Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE:BOOL=ON           \
  -D Kokkos_ENABLE_IMPL_VIEW_LEGACY:BOOL=ON                       \
  -D Kokkos_ENABLE_DEPRECATED_CODE:BOOL=OFF                       \
  -D Kokkos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=OFF                \
  ${ARCH_FLAG}

cmake --build ${WORK_DIR} --target install --parallel 4

export ERF_KOKKOS_HOME=${INSTALL_DIR}
