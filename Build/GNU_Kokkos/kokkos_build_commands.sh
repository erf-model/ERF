#!/bin/bash

# Build and install the Kokkos submodule for a CPU (Serial) build of ERF with RRTMGP.
#
# This is optional: the GNU Make build builds Kokkos itself under
# Submodules/kokkos/build/<config> when ERF_KOKKOS_HOME is not set, and the CMake
# build adds the submodule directly. Use this script to make one install that
# several builds share, then point them at it with ERF_KOKKOS_HOME (GNU Make) or
# Kokkos_ROOT (CMake). Set KOKKOS_ENABLE_OPENMP=ON for an OpenMP build.

# Set path to ERF
: "${ERF_DIR:=/home/ERF}"

# The Kokkos source, and the build and install directories (not in the source)
SRC_DIR=${ERF_DIR}/Submodules/kokkos
WORK_DIR=${SRC_DIR}/build/standalone
INSTALL_DIR=${WORK_DIR}/install

rm -rf ${WORK_DIR}/CMakeFiles ${WORK_DIR}/CMakeCache.txt

# These match the options Exec/Make.ERF builds Kokkos with
cmake -S ${SRC_DIR} -B ${WORK_DIR}                             \
  -D CMAKE_BUILD_TYPE:STRING=Release                          \
  -D CMAKE_CXX_COMPILER:STRING=mpicxx                         \
  -D CMAKE_CXX_STANDARD:STRING=20                             \
  -D CMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR}                 \
  -D CMAKE_INSTALL_LIBDIR:PATH=lib                            \
  -D Kokkos_ENABLE_SERIAL:BOOL=ON                             \
  -D Kokkos_ENABLE_OPENMP:BOOL=${KOKKOS_ENABLE_OPENMP:-OFF}   \
  -D Kokkos_ENABLE_IMPL_VIEW_LEGACY:BOOL=ON                   \
  -D Kokkos_ENABLE_DEPRECATED_CODE:BOOL=OFF                   \
  -D Kokkos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=OFF

cmake --build ${WORK_DIR} --target install --parallel 4

export ERF_KOKKOS_HOME=${INSTALL_DIR}
