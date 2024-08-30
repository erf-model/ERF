############################ BASE ######################################

cmake_minimum_required (VERSION 3.14 FATAL_ERROR)
project(ERF CXX C)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")
include(CMakePackageConfigHelpers)

########################## OPTIONS #####################################

#General options for all executables in the project
set(ERF_DIM "3" CACHE STRING "Number of physical dimensions")
option(ERF_ENABLE_DOCUMENTATION "Build documentation" OFF)
option(ERF_ENABLE_ALL_WARNINGS "Enable all compiler warnings" OFF)
option(ERF_ENABLE_TESTS "Enable regression and unit tests" OFF)
option(ERF_ENABLE_REGRESSION_TESTS_ONLY "Enable only regression tests" OFF)
option(ERF_USE_INTERNAL_AMREX "Add AMReX as subproject" ON)
option(ERF_ENABLE_NETCDF "Enable NetCDF IO" OFF)
option(ERF_ENABLE_HDF5 "Enable HDF5 IO" ${ERF_ENABLE_NETCDF})
option(ERF_ENABLE_PARTICLES "Enable Lagrangian particles" OFF)
option(ERF_ENABLE_FCOMPARE "Enable building fcompare when not testing" OFF)
set(ERF_PRECISION "DOUBLE" CACHE STRING "Floating point precision SINGLE or DOUBLE")

option(ERF_ENABLE_MOISTURE "Enable Full Moisture" ON)
option(ERF_ENABLE_WARM_NO_PRECIP "Enable Warm Moisture" OFF)
option(ERF_ENABLE_RRTMGP "Enable RTE-RRTMGP Radiation" OFF)

option(ERF_ENABLE_POISSON_SOLVE "Enable Poisson solve for incompressible flow" OFF)

#Options for performance
option(ERF_ENABLE_MPI "Enable MPI" OFF)
option(ERF_ENABLE_OPENMP "Enable OpenMP" OFF)
option(ERF_ENABLE_CUDA "Enable CUDA" OFF)
option(ERF_ENABLE_HIP  "Enable HIP" OFF)
option(ERF_ENABLE_SYCL "Enable SYCL" OFF)

#Options for C++
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(ERF_ENABLE_CUDA)
  enable_language(CUDA)
  if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "11.0")
    message(FATAL_ERROR "Your nvcc version is ${CMAKE_CUDA_COMPILER_VERSION} which is unsupported."
      "Please use CUDA toolkit version 11.0 or newer.")
  endif()
endif()

if(NOT ERF_DIM EQUAL 3)
  message(FATAL_ERROR "ERF is only supported in 3D.")
endif()

# Configure measuring code coverage in tests
option(CODECOVERAGE "Enable code coverage profiling" OFF)
if(CODECOVERAGE)
  # Only supports GNU
  if(NOT CMAKE_CXX_COMPILER_ID MATCHES GNU)
    message(WARNING "CODECOVERAGE is only support with GNU Compilers. The current C++ compiler is ${CMAKE_CXX_COMPILER_ID}")
  endif()
  if(NOT CMAKE_C_COMPILER_ID MATCHES GNU)
    message(WARNING "CODECOVERAGE is only support with GNU Compilers. The current C compiler is ${CMAKE_C_COMPILER_ID}")
  endif()

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
endif()

########################### AMReX #####################################

if (${ERF_USE_INTERNAL_AMREX})
    set(AMREX_SUBMOD_LOCATION "${CMAKE_SOURCE_DIR}/Submodules/AMReX")
    include(${CMAKE_SOURCE_DIR}/CMake/SetAmrexOptions.cmake)
    list(APPEND CMAKE_MODULE_PATH "${AMREX_SUBMOD_LOCATION}/Tools/CMake")
#    if (ERF_ENABLE_CUDA AND (CMAKE_VERSION VERSION_LESS 3.20))
#      include(AMReX_SetupCUDA)
#    endif()

########################### AMReX #####################################

    add_subdirectory(${AMREX_SUBMOD_LOCATION})
    set(FCOMPARE_EXE ${CMAKE_BINARY_DIR}/Submodules/AMReX/Tools/Plotfile/amrex_fcompare
      CACHE INTERNAL "Path to fcompare executable for regression tests")
else()
    set(CMAKE_PREFIX_PATH ${AMREX_DIR} ${CMAKE_PREFIX_PATH})
    list(APPEND AMREX_COMPONENTS
      "3D" "PIC" "PARTICLES" "PDOUBLE" "DOUBLE" "LSOLVERS")
    if (ERF_ENABLE_MPI)
      list(APPEND AMREX_COMPONENTS "MPI")
    endif()
    if (ERF_ENABLE_OPENMP)
      list(APPEND AMREX_COMPONENTS "OMP")
    endif()
    if (ERF_ENABLE_CUDA)
      list(APPEND AMREX_COMPONENTS "CUDA")
    endif()
    if (ERF_ENABLE_SYCL)
      list(APPEND AMREX_COMPONENTS "SYCL")
    endif()
    if (ERF_ENABLE_ROCM)
      list(APPEND AMREX_COMPONENTS "HIP")
    endif()
    if (ERF_ENABLE_HYPRE)
      list(APPEND AMREX_COMPONENTS "HYPRE")
    endif()
    if (ERF_ENABLE_TINY_PROFILE)
      list(APPEND AMREX_COMPONENTS "TINY_PROFILE")
    endif()
    separate_arguments(AMREX_COMPONENTS)
    find_package(AMReX CONFIG REQUIRED
      COMPONENTS ${AMREX_COMPONENTS})
    message(STATUS "Found AMReX = ${AMReX_DIR}")
    set(FCOMPARE_EXE ${AMReX_DIR}/../../../bin/amrex_fcompare
      CACHE INTERNAL "Path to fcompare executable for regression tests")
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

########################### RRTMGP #################################

if(ERF_ENABLE_RRTMGP)
   message(STATUS "Building RRTMGP/YAKL...")

   # YAKL_ARCH can be CUDA, HIP, SYCL, OPENMP45, or empty
   if(ERF_ENABLE_CUDA)
      set(YAKL_ARCH "CUDA")
      # CUDA_FLAGS is set the same as ERF_CUDA_FLAGS
      string(APPEND YAKL_CUDA_FLAGS " -arch sm_70")
      if(ENABLE_CUDA_FASTMATH)
        string(APPEND YAKL_CUDA_FLAGS " --use_fast_math")
      endif()
      set_cuda_architectures(AMReX_CUDA_ARCH)
   elseif(ERF_ENABLE_SYCL)
      set(YAKL_ARCH "SYCL")
      # SYCL_FLAGS is set through Macros.cmake / config_compilers.xml
      string(APPEND YAKL_SYCL_FLAGS " -fsycl")
   elseif(ERF_ENABLE_HIP)
      set(YAKL_ARCH "HIP")
      # SYCL_FLAGS is set through Macros.cmake / config_compilers.xml
      string(APPEND YAKL_HIP_FLAGS " -D__HIP_ROCclr__ -D__HIP_ARCH_GFX90A__=1 --rocm-path=${ROCM_PATH} --offload-arch=gfx90a -x hip")
   else()
      # For CPU C++ compilers duplicate flags are fine, the last ones win typically
      set(YAKL_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
      set(YAKL_ARCH "")
   endif()

   # Build YAKL as a static library
   # YAKL_HOME is YAKL's source directlry
   set(YAKL_HOME ${CMAKE_SOURCE_DIR}/Submodules/YAKL)
   # YAKL_BIN is where we're placing the YAKL library
   set(YAKL_BIN  ${CMAKE_BINARY_DIR}/yakl)
   # Build the YAKL static library
   add_subdirectory(${YAKL_HOME} ${YAKL_BIN})

   # Build the static rrtmgp library
   set(RRTMGP_BIN ${CMAKE_BINARY_DIR}/rrtmgp)
   add_subdirectory(${CMAKE_SOURCE_DIR}/Submodules/RRTMGP/cpp ${RRTMGP_BIN})

   # Find the radiation data files
   if (DEFINED ENV{ERF_RADIATION_DATA_DIR})
       set(ERF_RADIATION_DATA_DIR "$ENV{ERF_RADIATION_DATA_DIR}")
   else()
       message(FATAL_ERROR "Environment variable ERF_RADIATION_DATA_DIR not set!")
   endif()
endif()

########################### ERF #####################################

if(ERF_ENABLE_MPI)
  find_package(MPI REQUIRED)
endif()

configure_file(
  ${CMAKE_SOURCE_DIR}/Source/ERF_Config.H.in
  ${CMAKE_BINARY_DIR}/ERF_Config.H
  @ONLY
)

# General information about machine, compiler, and build type
message(STATUS "ERF Information:")
message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

# Turn on rpath stuff
include(${CMAKE_SOURCE_DIR}/CMake/SetRpath.cmake)

#Build erf executables and intermediate object library and link to amrex library
add_subdirectory(Exec)

if(ERF_ENABLE_TESTS)
  include(CTest)
  add_subdirectory(Tests)
endif()

if(ERF_ENABLE_DOCUMENTATION)
   add_subdirectory(Docs)
endif()
