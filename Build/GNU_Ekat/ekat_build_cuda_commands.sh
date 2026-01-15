#!/bin/bash

# Set path to ERF
${ERF_DIR:=/home/ERF}

# Set the path to ekat source
SRC_DIR=${ERF_DIR}/Submodules/ekat
cd ${SRC_DIR}

# Cannot build/install in the src dir
WORK_DIR=${SRC_DIR}/build
INSTALL_DIR=${WORK_DIR}/install
mkdir -p ${INSTALL_DIR}
cd ${WORK_DIR}

# Locations of our machine files and the generic ones
MY_MACH_FILE_DIR=${ERF_DIR}/Build/GNU_Ekat
MACH_FILE_DIR=${ERF_DIR}/Submodules/ekat/cmake/machine-files/kokkos

rm -rf CMakeFiles
rm -f  CMakeCache.txt

cmake \
  -D MACH_FILE_DIR="${MACH_FILE_DIR}"               \
  -C ${MY_MACH_FILE_DIR}/my_ekat_cuda.cmake         \
  -D CMAKE_BUILD_TYPE:STRING=Release                \
  -D EKAT_ENABLE_MPI:BOOL=ON                        \
  -D CMAKE_CXX_COMPILER:STRING=mpicxx               \
  -D CMAKE_Fortran_COMPILER:STRING=mpif90           \
  -D CMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR}       \
  -D EKAT_ENABLE_TESTS:BOOL=ON                      \
  -D EKAT_ENABLE_KOKKOS:BOOL=ON                     \
  -D EKAT_ENABLE_LOGGING:BOOL=ON                    \
  -D EKAT_TEST_DOUBLE_PRECISION:BOOL=ON             \
  -D EKAT_TEST_SINGLE_PRECISION:BOOL=ON             \
  ${SRC_DIR}

make install

export EKAT_HOME=${INSTALL_DIR}
