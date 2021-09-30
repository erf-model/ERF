############################ BASE ######################################

cmake_minimum_required (VERSION 3.14 FATAL_ERROR)
project(ERF CXX C Fortran)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")
include(CMakePackageConfigHelpers)

########################## OPTIONS #####################################

#General options for all executables in the project
set(ERF_DIM "3" CACHE STRING "Number of physical dimensions")
option(ERF_ENABLE_DOCUMENTATION "Build documentation" OFF)
option(ERF_ENABLE_AMREX_EB "Enable EB in AMReX library" OFF)
option(ERF_ENABLE_ALL_WARNINGS "Enable all compiler warnings" OFF)
option(ERF_ENABLE_TESTS "Enable regression and unit tests" OFF)
option(ERF_ENABLE_MASA "Enable tests that require MASA" OFF)
option(ERF_ENABLE_NETCDF "Enable NetCDF IO" OFF)
option(ERF_ENABLE_FCOMPARE "Enable building fcompare when not testing" OFF)
set(ERF_PRECISION "DOUBLE" CACHE STRING "Floating point precision SINGLE or DOUBLE")

#Options for performance
option(ERF_ENABLE_MPI "Enable MPI" OFF)
option(ERF_ENABLE_OPENMP "Enable OpenMP" OFF)
option(ERF_ENABLE_CUDA "Enable CUDA" OFF)

#Options for C++
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(ERF_ENABLE_CUDA)
  enable_language(CUDA)
  if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "9.0")
    message(FATAL_ERROR "Your nvcc version is ${CMAKE_CUDA_COMPILER_VERSION} which is unsupported."
      "Please use CUDA toolkit version 9.0 or newer.")
  endif()
endif()

if(NOT ERF_DIM EQUAL 3)
  message(FATAL_ERROR "ERF is only supported in 3D.")
endif()

# Configure measuring code coverage in tests
option(CODECOV "Enable code coverage profiling" OFF)
if(CODECOV)
  # Only supports GNU
  if(NOT CMAKE_CXX_COMPILER_ID MATCHES GNU)
    message(WARNING "CODECOV is only support with GNU Compilers. The current C++ compiler is ${CMAKE_CXX_COMPILER_ID}")
  endif()
  if(NOT CMAKE_C_COMPILER_ID MATCHES GNU)
    message(WARNING "CODECOV is only support with GNU Compilers. The current C compiler is ${CMAKE_C_COMPILER_ID}")
  endif()

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
endif()

########################### AMReX #####################################

set(AMREX_SUBMOD_LOCATION "${CMAKE_SOURCE_DIR}/Submodules/AMReX")
include(${CMAKE_SOURCE_DIR}/CMake/SetAmrexOptions.cmake)
list(APPEND CMAKE_MODULE_PATH "${AMREX_SUBMOD_LOCATION}/Tools/CMake")
add_subdirectory(${AMREX_SUBMOD_LOCATION})
include(${CMAKE_SOURCE_DIR}/CMake/SetAmrexCompileFlags.cmake)

########################### MASA #####################################

if(ERF_ENABLE_MASA)
  set(CMAKE_PREFIX_PATH ${MASA_DIR} ${CMAKE_PREFIX_PATH})
  find_package(MASA QUIET REQUIRED)
  if(MASA_FOUND)
    message(STATUS "Found MASA = ${MASA_DIR}")
  endif()
endif()

########################## NETCDF ##################################

if(ERF_ENABLE_NETCDF)
  set(CMAKE_PREFIX_PATH ${NETCDF_DIR} ${CMAKE_PREFIX_PATH})
##  set(NETCDF_CXX "YES")
  find_package (NetCDF REQUIRED)
  if(NETCDF_FOUND)
     message(STATUS "Found NetCDF, NETCDF_DIR = ${NETCDF_DIR}")
  endif() 
endif()

########################### ERF #####################################

if(ERF_ENABLE_MPI)
  find_package(MPI REQUIRED)
endif()

# General information about machine, compiler, and build type
message(STATUS "ERF Information:")
message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

# Turn on rpath stuff
include(${CMAKE_SOURCE_DIR}/CMake/SetRpath.cmake)

#Build erf executables and link to amrex library
add_subdirectory(Exec)

if(ERF_ENABLE_TESTS)
  include(CTest)
  add_subdirectory("Submodules/GoogleTest")
  add_subdirectory(Tests)
endif()

if(ERF_ENABLE_DOCUMENTATION)
   add_subdirectory(Docs)
endif()
