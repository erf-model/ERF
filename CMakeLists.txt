############################ BASE ######################################

cmake_minimum_required (VERSION 3.14 FATAL_ERROR)

# Set main ERF context
list(APPEND CMAKE_MESSAGE_CONTEXT "ERF")

# Include Cray compiler detection BEFORE project() to set compilers
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")
include(CrayCompilerDetection)

if(ERF_ENABLE_MORR_FORT OR ERF_ENABLE_NOAHMP)
  project(ERF CXX C Fortran)
else()
  project(ERF CXX C)
endif()

message(STATUS "Configuring ERF")
message(VERBOSE "Source directory: ${CMAKE_SOURCE_DIR}")
message(VERBOSE "Build directory: ${CMAKE_BINARY_DIR}")

function(print_target_properties target_name)
  if(NOT TARGET ${target_name})
    message(STATUS "=== Target ${target_name} does not exist ===")
    return()
  endif()

  message(STATUS "=== Target: ${target_name} ===")

  # Define properties to check
  set(prop_names
    INTERFACE_INCLUDE_DIRECTORIES
    INTERFACE_LINK_OPTIONS
    INTERFACE_LINK_LIBRARIES
  )

  set(display_names
    "Include Directories"
    "Link Options"
    "Link Libraries"
  )

  list(LENGTH prop_names num_props)
  math(EXPR last_index "${num_props} - 1")

  foreach(idx RANGE ${last_index})
    list(GET prop_names ${idx} prop_name)
    list(GET display_names ${idx} display_name)

    get_target_property(values ${target_name} ${prop_name})
    if(values)
      message(VERBOSE "  ${display_name}:")
      foreach(value IN LISTS values)
        message(VERBOSE "    ${value}")
      endforeach()
    endif()
  endforeach()

  message(STATUS "")
endfunction()

# Find NVHPC package and create aliases if needed
if(ERF_ENABLE_CUDA AND ERF_ENABLE_NVHPC)
  list(APPEND CMAKE_MESSAGE_CONTEXT "NVHPC")
  message(STATUS "Configuring NVHPC CUDA support")

  find_package(NVHPC QUIET COMPONENTS MATH CUDA)

  if(NOT NVHPC_FOUND)
    message(STATUS "")
    message(STATUS "====================================================================")
    message(STATUS "NVHPC Detection Failed")
    message(STATUS "====================================================================")
    message(STATUS "")
    message(STATUS "To resolve:")
    message(STATUS "  Load NVHPC module or set NVHPC_ROOT")
    message(STATUS "")
    message(STATUS "  Example:")
    message(STATUS "    module load nvhpc")
    message(STATUS "    cmake -DNVHPC_ROOT=/path/to/nvhpc ..")
    message(STATUS "")
    message(STATUS "====================================================================")
    message(STATUS "")
    message(FATAL_ERROR "NVHPC required but not found")
  endif()

  message(STATUS "Found NVHPC")

  function(create_cuda_alias nvhpc_target cuda_name)
    if(TARGET NVHPC::${nvhpc_target} AND NOT TARGET CUDA::${cuda_name})
      add_library(CUDA::${cuda_name} ALIAS NVHPC::${nvhpc_target})
      message(STATUS "  Created alias: CUDA::${cuda_name} -> NVHPC::${nvhpc_target}")
    elseif(NOT TARGET NVHPC::${nvhpc_target})
      message(WARNING "  Cannot create alias CUDA::${cuda_name}: NVHPC::${nvhpc_target} not found")
    endif()
  endfunction()

  create_cuda_alias(CUBLAS cublas)
  create_cuda_alias(CUBLAS_STATIC cublas_static)
  create_cuda_alias(CURAND curand)
  create_cuda_alias(CURAND_STATIC curand_static)
  create_cuda_alias(CUSPARSE cusparse)
  create_cuda_alias(CUSPARSE_STATIC cusparse_static)

  list(POP_BACK CMAKE_MESSAGE_CONTEXT)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")
include(CMakePackageConfigHelpers)

# Include Cray/Perlmutter auto-detection and fixes
include(CrayDetection)
include(UtilityTargets)

########################## OPTIONS #####################################

message(VERBOSE "Configuring build options")

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

option(ERF_ENABLE_RRTMGP "Enable RTE-RRTMGP Radiation" OFF)
option(ERF_ENABLE_SHOC "Enable SHOC" OFF)

#Options for performance
option(ERF_ENABLE_MPI "Enable MPI" OFF)
option(ERF_ENABLE_OPENMP "Enable OpenMP" OFF)
option(ERF_ENABLE_CUDA "Enable CUDA" OFF)
option(ERF_ENABLE_HIP  "Enable HIP" OFF)
option(ERF_ENABLE_SYCL "Enable SYCL" OFF)
option(ERF_ENABLE_NVHPC "Enable NVHPC" OFF)

#Options for NOAH-MP
option(ERF_ENABLE_NOAHMP "Enable Noah-MP" OFF)

#Option to build tools
option(ERF_ENABLE_TOOLS "Enable building of additional tools" OFF)

#Options for C++
set(CMAKE_CXX_STANDARD 17)
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
  message(VERBOSE "Enabling code coverage profiling")
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

list(APPEND CMAKE_MESSAGE_CONTEXT "AMReX")

if (${ERF_USE_INTERNAL_AMREX})
    message(STATUS "Using internal AMReX submodule")
    set(AMREX_SUBMOD_LOCATION "${CMAKE_SOURCE_DIR}/Submodules/AMReX")
    include(${CMAKE_SOURCE_DIR}/CMake/SetAmrexOptions.cmake)
    list(APPEND CMAKE_MODULE_PATH "${AMREX_SUBMOD_LOCATION}/Tools/CMake")

    add_subdirectory(${AMREX_SUBMOD_LOCATION})

    # Print AMReX properties
    print_target_properties(amrex)

    if(WIN32)
      set(FCOMPARE_EXE ${CMAKE_BINARY_DIR}/Submodules/AMReX/Tools/Plotfile/*/amrex_fcompare.exe
        CACHE STRING "Path to fcompare executable for regression tests")
    else()
      set(FCOMPARE_EXE ${CMAKE_BINARY_DIR}/Submodules/AMReX/Tools/Plotfile/amrex_fcompare
        CACHE STRING "Path to fcompare executable for regression tests")
    endif()
    message(VERBOSE "fcompare executable: ${FCOMPARE_EXE}")
else()
    message(STATUS "Using external AMReX")

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
    message(VERBOSE "Required AMReX components: ${AMREX_COMPONENTS}")

    find_package(AMReX CONFIG QUIET COMPONENTS ${AMREX_COMPONENTS})

    if(NOT AMReX_FOUND)
      message(STATUS "")
      message(STATUS "====================================================================")
      message(STATUS "AMReX Detection Failed")
      message(STATUS "====================================================================")
      message(STATUS "")
      message(STATUS "Required components: ${AMREX_COMPONENTS}")
      message(STATUS "")
      message(STATUS "To resolve:")
      message(STATUS "")
      message(STATUS "  Option 1: Use internal AMReX (recommended):")
      message(STATUS "    cmake -DERF_USE_INTERNAL_AMREX=ON ..")
      message(STATUS "")
      message(STATUS "  Option 2: Build and install AMReX separately, then:")
      message(STATUS "    cmake -DAMReX_DIR=/path/to/amrex/lib/cmake/AMReX ..")
      message(STATUS "")
      message(STATUS "====================================================================")
      message(STATUS "")
      message(FATAL_ERROR "AMReX required but not found")
    endif()

    message(STATUS "Found AMReX: ${AMReX_DIR}")
    message(VERBOSE "  AMReX_VERSION: ${AMReX_VERSION}")

    # Print AMReX properties
    print_target_properties(amrex)

    if(WIN32)
      set(FCOMPARE_EXE ${AMReX_DIR}/../../../*/amrex_fcompare.exe
        CACHE STRING "Path to fcompare executable for regression tests")
    else()
      set(FCOMPARE_EXE ${AMReX_DIR}/../../../bin/amrex_fcompare
        CACHE STRING "Path to fcompare executable for regression tests")
    endif()
    message(VERBOSE "fcompare executable: ${FCOMPARE_EXE}")

endif()

list(POP_BACK CMAKE_MESSAGE_CONTEXT)

########################## EKAT ##################################

if(ERF_ENABLE_RRTMGP OR ERF_ENABLE_SHOC OR ERF_ENABLE_P3)
  set(ERF_ENABLE_EKAT ON)
else()
  set(ERF_ENABLE_EKAT OFF)
endif()

if(ERF_ENABLE_EKAT)
  list(APPEND CMAKE_MESSAGE_CONTEXT "EKAT")

  message(STATUS "Configuring EKAT")

  if(NOT ERF_ENABLE_MPI)
    message(FATAL_ERROR "CMake Error: MPI must be enabled if EKAT is enabled.")
  endif()

  # NOTE: EKAT provides KOKKOS, so set relevant flags for KOKKOS
  if(ERF_ENABLE_CUDA)
    message(VERBOSE "Enabling Kokkos CUDA support")
    set(Kokkos_ENABLE_CUDA ON CACHE BOOL "kokkos enable cuda")
    set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "kokkos enable cuda lambda")
    set(Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE ON CACHE BOOL "kokkos enable cuda RDC")
  endif()

  if(ERF_ENABLE_HIP)
    message(VERBOSE "Enabling Kokkos HIP support")
    set(Kokkos_ENABLE_HIP ON CACHE BOOL "kokkos enable hip")
    set(Kokkos_ENABLE_HIP_LAMBDA ON CACHE BOOL "kokkos enable hip lambda")
    set(Kokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE ON CACHE BOOL "kokkos enable hip RDC")
  endif()

  if(ERF_ENABLE_SYCL)
    message(VERBOSE "Enabling Kokkos SYCL support")
    set(Kokkos_ENABLE_SYCL ON CACHE BOOL "kokkos enable sycl")
    set(Kokkos_ENABLE_SYCL_LAMBDA ON CACHE BOOL "kokkos enable sycl lambda")
  endif()

  set(CMAKE_CXX_STANDARD 17)
  set(EKAT_ENABLE_MPI ON CACHE BOOL "ekat enable mpi")
  set(EKAT_ENABLE_TESTS OFF CACHE BOOL "ekat enable tests")
  set(EKAT_ENABLE_KOKKOS ON CACHE BOOL "ekat enable kokkos")
  set(EKAT_ENABLE_LOGGING ON CACHE BOOL "ekat enable spdlog")
  set(EKAT_BIN ${CMAKE_BINARY_DIR}/Submodules/ekat)

  # Force use of internal spdlog from EKAT
  message(STATUS "Forcing Ekat to build their own spdlog")
  set(CMAKE_DISABLE_FIND_PACKAGE_spdlog TRUE)

  message(VERBOSE "EKAT binary directory: ${EKAT_BIN}")
  add_subdirectory(${CMAKE_SOURCE_DIR}/Submodules/ekat ${EKAT_BIN})

  # Print Kokkos properties
  if(TARGET Kokkos::kokkos)
    print_target_properties(Kokkos::kokkos)
  elseif(TARGET kokkoscore)
    print_target_properties(kokkoscore)
  endif()

# Strip GPU flags from Kokkos (AMReX provides them)
if(ERF_ENABLE_HIP AND ERF_ENABLE_EKAT)
  # Clear ALL link options from Kokkos targets - AMReX provides everything we need
  foreach(kokkos_target IN ITEMS kokkoscore kokkoscontainers kokkosalgorithms Kokkos::kokkos)
    if(TARGET ${kokkos_target})
      get_target_property(kokkos_link_opts ${kokkos_target} INTERFACE_LINK_OPTIONS)
      if(kokkos_link_opts AND NOT kokkos_link_opts STREQUAL "kokkos_link_opts-NOTFOUND")
        # Clear ALL link options - AMReX provides GPU flags
        set_target_properties(${kokkos_target} PROPERTIES INTERFACE_LINK_OPTIONS "")
        message(STATUS "Cleared link options from ${kokkos_target} (AMReX provides GPU flags)")
      endif()
    endif()
  endforeach()

  # Print kokkoscore after clearing
  print_target_properties(kokkoscore)
endif()
# Strip GPU flags from Kokkos (AMReX provides them)
#if(ERF_ENABLE_HIP AND ERF_ENABLE_EKAT)
#  # Strip from ALL Kokkos targets
#  foreach(kokkos_target IN ITEMS kokkoscore kokkoscontainers kokkosalgorithms Kokkos::kokkos)
#    if(TARGET ${kokkos_target})
#      get_target_property(kokkos_link_opts ${kokkos_target} INTERFACE_LINK_OPTIONS)
#      if(kokkos_link_opts AND NOT kokkos_link_opts STREQUAL "kokkos_link_opts-NOTFOUND")
#        # Filter out GPU-related generator expressions
#        list(FILTER kokkos_link_opts EXCLUDE REGEX ".*--hip-link.*")
#        list(FILTER kokkos_link_opts EXCLUDE REGEX ".*--offload-arch=.*")
#        list(FILTER kokkos_link_opts EXCLUDE REGEX ".*-fgpu-rdc.*")
#        set_target_properties(${kokkos_target} PROPERTIES INTERFACE_LINK_OPTIONS "${kokkos_link_opts}")
#        message(STATUS "Stripped GPU flags from ${kokkos_target}")
#      endif()
#    endif()
#  endforeach()

#  # Print kokkoscore after stripping
#  print_target_properties(kokkoscore)
#endif()

list(POP_BACK CMAKE_MESSAGE_CONTEXT)
endif()

########################### MPI #####################################

if(ERF_ENABLE_MPI)
  list(APPEND CMAKE_MESSAGE_CONTEXT "MPI")
  message(STATUS "Configuring MPI")

  # Check if we're on Cray with bare MPI wrappers (which will hang)
  set(SKIP_MPI_DETECTION FALSE)

  if(DEFINED ENV{CRAYPE_VERSION} OR DEFINED ENV{CRAY_MPICH_DIR})
    # On Cray system - check if using problematic bare MPI wrappers
    if(CMAKE_CXX_COMPILER MATCHES "mpicxx" OR
       CMAKE_C_COMPILER MATCHES "mpicc" OR
       CMAKE_Fortran_COMPILER MATCHES "mpifort")
      message(STATUS "Detected bare MPI wrappers on Cray - skipping detection (would hang)")
      set(SKIP_MPI_DETECTION TRUE)
    endif()
  endif()

  if(SKIP_MPI_DETECTION)
    # Workaround: Manual MPI setup (avoids hang)
    message(VERBOSE "Manually configuring MPI (bypassing find_package)")

    # Get Cray MPICH version for informational purposes
    set(MPICH_VERSION "UNKNOWN")
    if(DEFINED ENV{CRAY_MPICH_VERSION})
      set(MPICH_VERSION "$ENV{CRAY_MPICH_VERSION}")
    elseif(DEFINED ENV{CRAY_MPICH_VER})
      set(MPICH_VERSION "$ENV{CRAY_MPICH_VER}")
    endif()

    # Cray MPICH 8.x supports MPI 3.1 standard
    set(MPI_VERSION "3.1")

    # Create MPI targets
    if(NOT TARGET MPI::MPI_CXX)
      add_library(MPI::MPI_CXX INTERFACE IMPORTED)
      message(DEBUG "Created MPI::MPI_CXX target")
    endif()

    if(NOT TARGET MPI::MPI_C)
      add_library(MPI::MPI_C INTERFACE IMPORTED)
      message(DEBUG "Created MPI::MPI_C target")
    endif()

    if(NOT TARGET MPI::MPI_Fortran)
      add_library(MPI::MPI_Fortran INTERFACE IMPORTED)
      message(DEBUG "Created MPI::MPI_Fortran target")
    endif()

    # Set MPI variables
    set(MPI_FOUND TRUE)
    set(MPI_CXX_FOUND TRUE)
    set(MPI_C_FOUND TRUE)
    set(MPI_Fortran_FOUND TRUE)
    set(MPI_C_VERSION "${MPI_VERSION}")
    set(MPI_CXX_VERSION "${MPI_VERSION}")
    set(MPI_Fortran_VERSION "${MPI_VERSION}")

    message(STATUS "Cray MPICH: ${MPICH_VERSION}")
    message(STATUS "MPI standard: ${MPI_VERSION}")
    message(VERBOSE "Created MPI::MPI_CXX and MPI::MPI_C targets")
  else()
    # Normal path: Use find_package with QUIET
    message(VERBOSE "Using find_package(MPI) for detection")
    set(_mpi_comps C CXX)
    if(ERF_ENABLE_MORR_FORT OR ERF_ENABLE_NOAHMP)
      list(APPEND _mpi_comps Fortran)
    endif()
    message(DEBUG "MPI components: ${_mpi_comps}")

    find_package(MPI QUIET COMPONENTS ${_mpi_comps})

    if(NOT MPI_FOUND)
      message(STATUS "")
      message(STATUS "====================================================================")
      message(STATUS "MPI Detection Failed")
      message(STATUS "====================================================================")
      message(STATUS "")
      message(STATUS "To resolve, load modules from your machine profile:")
      message(STATUS "")
      erf_suggest_machine_profile()
      message(STATUS "")
      message(STATUS "  Or on non-Cray systems, install MPI:")
      message(STATUS "    OpenMPI, MPICH, Intel MPI, etc.")
      message(STATUS "")
      message(STATUS "  Then configure with:")
      message(STATUS "    cmake -DMPI_C_COMPILER=mpicc -DMPI_CXX_COMPILER=mpicxx ..")
      message(STATUS "")
      message(STATUS "====================================================================")
      message(STATUS "")
      message(FATAL_ERROR "MPI required but not found")
    endif()

    message(STATUS "Found MPI")
    message(VERBOSE "  MPI_C_VERSION: ${MPI_C_VERSION}")
    message(VERBOSE "  MPI_CXX_VERSION: ${MPI_CXX_VERSION}")
  endif()

  message(STATUS "MPI configuration complete")
  list(POP_BACK CMAKE_MESSAGE_CONTEXT)
else()
  message(DEBUG "MPI not enabled")
endif()

########################## NETCDF ##################################

if(ERF_ENABLE_NETCDF OR ERF_ENABLE_TOOLS)
  list(APPEND CMAKE_MESSAGE_CONTEXT "NetCDF")

  message(STATUS "Configuring NetCDF")
  message(DEBUG "ERF_ENABLE_NETCDF: ${ERF_ENABLE_NETCDF}")
  message(DEBUG "ERF_ENABLE_TOOLS: ${ERF_ENABLE_TOOLS}")

  set(CMAKE_PREFIX_PATH ${NETCDF_DIR} ${CMAKE_PREFIX_PATH})
  find_package (NetCDF REQUIRED)
  if(NETCDF_FOUND)
     message(STATUS "Found NetCDF: ${NETCDF_DIR}")
     message(VERBOSE "NetCDF includes: ${NETCDF_INCLUDES}")
     message(VERBOSE "NetCDF libraries: ${NETCDF_LIBRARIES}")
  endif()

  message(STATUS "NetCDF configuration complete")
  list(POP_BACK CMAKE_MESSAGE_CONTEXT)
endif()

########################## NOAH-MP ##################################

if(ERF_ENABLE_NOAHMP)
  list(APPEND CMAKE_MESSAGE_CONTEXT "Noah-MP")

  message(STATUS "Configuring Noah-MP")

  if(ERF_ENABLE_NETCDF)
     set(NOAHMP_HOME ${CMAKE_SOURCE_DIR}/Submodules/Noah-MP/drivers/erf)
     set(NOAHMP_BIN  ${CMAKE_BINARY_DIR}/Submodules/Noah-MP/drivers/erf)
     message(VERBOSE "Noah-MP source: ${NOAHMP_HOME}")
     message(VERBOSE "Noah-MP binary: ${NOAHMP_BIN}")
     add_subdirectory(${NOAHMP_HOME} ${NOAHMP_BIN})
  else()
     message(FATAL_ERROR "Noah-MP requires NetCDF be enabled")
  endif()

  message(STATUS "Noah-MP configuration complete")
  list(POP_BACK CMAKE_MESSAGE_CONTEXT)
endif()

########################### RRTMGP #################################

if(ERF_ENABLE_RRTMGP)
  list(APPEND CMAKE_MESSAGE_CONTEXT "RRTMGP")

  message(STATUS "Configuring RRTMGP")

  if(NOT ERF_ENABLE_EKAT)
      message(FATAL_ERROR "EKAT must be enabled if RRTMGP is enabled")
  endif()

  if(NOT ERF_ENABLE_NETCDF)
      message(FATAL_ERROR "NetCDF must be enabled if RRTMGP is enabled")
  endif()

  message(VERBOSE "Building RRTMGP with Kokkos support")
  add_compile_options(-Wno-c++11-narrowing)
  message(DEBUG "Adding -Wno-c++11-narrowing for mo_gas_optics_rrtmgp issue on sycl/hip RRTMGP target")

  # Build the static rrtmgp library
  set(RRTMGP_BIN ${CMAKE_BINARY_DIR}/Submodules/rrtmgp)
  message(VERBOSE "RRTMGP binary directory: ${RRTMGP_BIN}")
  add_subdirectory(${CMAKE_SOURCE_DIR}/Submodules/RRTMGP/cpp ${RRTMGP_BIN})

  # Set up kokkos library and definitions
  set(RRTMGP_ENABLE_KOKKOS TRUE)
  target_compile_definitions(rrtmgp PUBLIC RRTMGP_ENABLE_KOKKOS)
  target_link_libraries(rrtmgp kokkos)
  message(DEBUG "Added Kokkos to RRTMGP target")

  # Print RRTMGP properties
  print_target_properties(rrtmgp)

  message(STATUS "RRTMGP configuration complete")
  list(POP_BACK CMAKE_MESSAGE_CONTEXT)
endif()

########################### SHOC #################################

if(ERF_ENABLE_SHOC)
  list(APPEND CMAKE_MESSAGE_CONTEXT "SHOC")

  message(STATUS "Configuring SHOC")

  if(NOT ERF_ENABLE_EKAT)
      message(FATAL_ERROR "EKAT must be enabled if SHOC is enabled")
  endif()

  # NOTE: We compile shoc src files directly
  message(VERBOSE "SHOC source files will be compiled directly")

  # E3SM Cloned Check
  if(ERF_ENABLE_SHOC OR ERF_ENABLE_P3)
    message(VERBOSE "Checking EAMxx files from E3SM")

    set(E3SM_EXPECTED_PATH "${CMAKE_SOURCE_DIR}/external/E3SM")

    if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics")
        message(FATAL_ERROR
            "E3SM provides eamxx required for SHOC/P3 but not found.\n"
            "\n"
            "This requires specific EAMxx physics source files from E3SM at commit 2eb52d9.\n"
            "\n"
            "Option 1 - Use provided script (recommended):\n"
            "  ${CMAKE_SOURCE_DIR}/Build/eamxx_clone.sh\n"
            "\n"
            "Option 2 - Symlink existing E3SM (not tested, must be at commit 2eb52d9b3d):\n"
            "  ln -s /path/to/your/E3SM ${CMAKE_SOURCE_DIR}/external/E3SM\n"
            "\n"
            "Or disable features:\n"
            "  cmake -DERF_ENABLE_SHOC=OFF -DERF_ENABLE_P3=OFF ..\n")
    endif()
  else()
    message(DEBUG "EKAT physics not enabled, skipping E3SM check")
  endif()


  message(STATUS "SHOC configuration complete")
  list(POP_BACK CMAKE_MESSAGE_CONTEXT)
endif()

########################### ERF #####################################

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
  message(STATUS "Configuring tests")
  include(CTest)
  add_subdirectory(Tests)
endif()

if(ERF_ENABLE_DOCUMENTATION)
  message(STATUS "Configuring documentation")
  add_subdirectory(Docs)
endif()

# Installation rules
list(APPEND CMAKE_MESSAGE_CONTEXT "Install")

message(STATUS "Configuring installation")

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

# Create non-object library for use as external target
add_library(erf_api)
if (BUILD_SHARED_LIBS)
  set_target_properties(erf_api PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()
if (ERF_ENABLE_CUDA)
   setup_target_for_cuda_compilation(erf_api)
endif ()
target_link_libraries(erf_api PUBLIC erf_srclib)
add_library(${PROJECT_NAME}::erf_api ALIAS erf_srclib)

# Collect all headers and make them installable with the target
file(GLOB_RECURSE ERF_INCLUDES "${CMAKE_SOURCE_DIR}/Source/*.H")

set_target_properties(
  erf_srclib PROPERTIES PUBLIC_HEADER "${ERF_INCLUDES}")

# Install ERF
install(
  TARGETS erf_api erf_srclib
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin
  ARCHIVE       DESTINATION lib
  LIBRARY       DESTINATION lib
  INCLUDES      DESTINATION include
  PUBLIC_HEADER DESTINATION include
  )

# Install all headers in include directories
#install(
#  DIRECTORY ${ERF_INCLUDE_DIRECTORIES}
#  DESTINATION include
#  FILES_MATCHING PATTERN "*.H")

# Make ERF discoverable using `find_package`
install(
  EXPORT ${PROJECT_NAME}Targets
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
  )
configure_package_config_file(
  CMake/${PROJECT_NAME}Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
  )
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
  )

message(VERBOSE "Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(VERBOSE "Install libdir: ${CMAKE_INSTALL_LIBDIR}")

list(POP_BACK CMAKE_MESSAGE_CONTEXT)

# Pop main ERF context
list(POP_BACK CMAKE_MESSAGE_CONTEXT)

message(STATUS "ERF configuration complete")
