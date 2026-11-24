#!/bin/bash

# Build and install the Kokkos submodule for a CUDA build of ERF with RRTMGP
# on Perlmutter (A100, Cray PE). Source this script so that ERF_KOKKOS_HOME is set.
#
# This is optional: the GNU Make build builds Kokkos itself under
# Submodules/kokkos/build/<config> when ERF_KOKKOS_HOME is not set.

# Set ERF_DIR relative to the script location
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ERF_DIR=$(cd "${SCRIPT_DIR}/../.." && pwd)
echo "ERF_DIR set to: $ERF_DIR"

# The Kokkos source, and the build and install directories (not in the source)
SRC_DIR=${ERF_DIR}/Submodules/kokkos
WORK_DIR=${SRC_DIR}/build/standalone-perlmutter
INSTALL_DIR=${WORK_DIR}/install

# nvcc_wrapper's host compiler: Cray's C++ wrapper, never nvcc itself
export NVCC_WRAPPER_DEFAULT_COMPILER=$(command -v CC)

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
  -D Kokkos_ARCH_AMPERE80:BOOL=ON

# Build and install
cmake --build ${WORK_DIR} --target install --parallel 8

export ERF_KOKKOS_HOME=${INSTALL_DIR}
echo "Kokkos installed in ${ERF_KOKKOS_HOME}"
