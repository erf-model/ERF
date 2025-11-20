# ==============================================================================
# Cray System Auto-Detection and Workarounds
# ==============================================================================
# This module detects Cray systems and automatically applies workarounds for
# common build issues. Each fix corresponds to a documented checklist item.
#
# CMake 3.25+ version using standard message log levels:
#   cmake ..                      # Quiet (STATUS messages only)
#   cmake --log-level=VERBOSE ..  # Show detection details
#   cmake --log-level=DEBUG ..    # Show all diagnostics
#   cmake --log-context ..        # Show message hierarchy
#
# Options:
#   -DERF_ENABLE_CRAY_AUTO_FIXES=OFF : Disable automatic Cray system fixes
# ==============================================================================

option(ERF_ENABLE_CRAY_AUTO_FIXES "Enable automatic Cray system fixes" ON)

# Set Cray context for hierarchical logging
list(APPEND CMAKE_MESSAGE_CONTEXT "Cray")

if(NOT ERF_ENABLE_CRAY_AUTO_FIXES)
    message(STATUS "Auto-fixes disabled by user")
    list(POP_BACK CMAKE_MESSAGE_CONTEXT)
    return()
endif()

message(DEBUG "Starting Cray detection and workaround application")

# ==============================================================================
# Detect Cray Environment
# ==============================================================================

set(ERF_ON_CRAY FALSE)

message(DEBUG "Checking for Cray environment")
message(TRACE "  CMAKE_C_COMPILER: ${CMAKE_C_COMPILER}")
message(TRACE "  CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")
message(TRACE "  CRAY_MPICH_DIR: $ENV{CRAY_MPICH_DIR}")

# Check for Cray compiler wrappers
if(CMAKE_C_COMPILER MATCHES ".*cc$" AND 
   CMAKE_CXX_COMPILER MATCHES ".*CC$" AND
   DEFINED ENV{CRAY_MPICH_DIR})
    set(ERF_ON_CRAY TRUE)
    message(STATUS "Detected Cray system via compiler wrappers")
    message(VERBOSE "  C compiler: ${CMAKE_C_COMPILER}")
    message(VERBOSE "  C++ compiler: ${CMAKE_CXX_COMPILER}")
    message(VERBOSE "  CRAY_MPICH_DIR: $ENV{CRAY_MPICH_DIR}")
endif()

# Additional check for Cray environment variables
if(DEFINED ENV{CRAYPE_VERSION})
    set(ERF_ON_CRAY TRUE)
    message(STATUS "Detected Cray Programming Environment")
    message(VERBOSE "  CRAYPE_VERSION: $ENV{CRAYPE_VERSION}")
endif()

if(NOT ERF_ON_CRAY)
    message(STATUS "Not on a Cray system, skipping Cray-specific fixes")
    message(DEBUG "Detection criteria not met:")
    message(DEBUG "  Compiler wrappers cc/CC: NO")
    message(DEBUG "  CRAY_MPICH_DIR set: NO")
    message(DEBUG "  CRAYPE_VERSION set: NO")
    list(POP_BACK CMAKE_MESSAGE_CONTEXT)
    return()
endif()

# ==============================================================================
# Optional: Check for Stale Configuration
# ==============================================================================

option(ERF_CHECK_MODULES "Check for stale configuration from module changes" ON)

if(ERF_CHECK_MODULES)
    list(APPEND CMAKE_MESSAGE_CONTEXT "CrayConfigCheck")

    message(DEBUG "Starting configuration verification")

    # Detection log for issues
    set(STALE_CONFIG_LOG "")

    # Determine if this is first configure (before we cache anything)
    set(IS_FIRST_CONFIGURE FALSE)
    if(NOT DEFINED CACHED_LOADED_MODULES)
        set(IS_FIRST_CONFIGURE TRUE)
        message(DEBUG "First configure detected")
    endif()

    # Check 1: Module environment changed
    if(DEFINED ENV{LOADEDMODULES})
        set(CURRENT_MODULES "$ENV{LOADEDMODULES}")
        if(DEFINED CACHED_LOADED_MODULES)
            if(NOT "${CURRENT_MODULES}" STREQUAL "${CACHED_LOADED_MODULES}")
                message(VERBOSE "Module environment changed since last configure")
                list(APPEND STALE_CONFIG_LOG "LOADEDMODULES changed")
                list(APPEND STALE_CONFIG_LOG "  Previous: ${CACHED_LOADED_MODULES}")
                list(APPEND STALE_CONFIG_LOG "  Current:  ${CURRENT_MODULES}")
            else()
                message(DEBUG "Module environment unchanged")
            endif()
        endif()
        set(CACHED_LOADED_MODULES "${CURRENT_MODULES}" CACHE INTERNAL "Modules at configure time")
    endif()

    # Check 2: PE_ENV changed
    if(DEFINED ENV{PE_ENV})
        set(CURRENT_PE_ENV "$ENV{PE_ENV}")
        if(DEFINED CACHED_PE_ENV AND NOT "${CURRENT_PE_ENV}" STREQUAL "${CACHED_PE_ENV}")
            message(VERBOSE "PE_ENV changed: ${CACHED_PE_ENV} -> ${CURRENT_PE_ENV}")
            list(APPEND STALE_CONFIG_LOG "PE_ENV changed from ${CACHED_PE_ENV} to ${CURRENT_PE_ENV}")
        endif()
        set(CACHED_PE_ENV "${CURRENT_PE_ENV}" CACHE INTERNAL "")
    endif()

    # Check 3: Compiler version changed
    if(DEFINED CMAKE_CXX_COMPILER_VERSION)
        if(DEFINED CACHED_CXX_COMPILER_VERSION AND NOT "${CMAKE_CXX_COMPILER_VERSION}" STREQUAL "${CACHED_CXX_COMPILER_VERSION}")
            message(VERBOSE "Compiler version changed")
            list(APPEND STALE_CONFIG_LOG "Compiler version changed from ${CACHED_CXX_COMPILER_VERSION} to ${CMAKE_CXX_COMPILER_VERSION}")
        endif()
        set(CACHED_CXX_COMPILER_VERSION "${CMAKE_CXX_COMPILER_VERSION}" CACHE INTERNAL "")
    endif()

    # Check 4: CMAKE_*_STANDARD_LIBRARIES already contains MPI (from previous run)
    if(NOT IS_FIRST_CONFIGURE)  # Only check on reconfigure
        if(DEFINED CMAKE_CXX_STANDARD_LIBRARIES AND CMAKE_CXX_STANDARD_LIBRARIES)
            if(CMAKE_CXX_STANDARD_LIBRARIES MATCHES "mpi_")
                message(VERBOSE "CMAKE_CXX_STANDARD_LIBRARIES already contains MPI libraries")
                list(APPEND STALE_CONFIG_LOG "CMAKE_CXX_STANDARD_LIBRARIES pre-populated with MPI libs")
                list(APPEND STALE_CONFIG_LOG "  Found: ${CMAKE_CXX_STANDARD_LIBRARIES}")
            endif()
        endif()

        if(DEFINED CMAKE_CUDA_STANDARD_LIBRARIES AND CMAKE_CUDA_STANDARD_LIBRARIES)
            if(CMAKE_CUDA_STANDARD_LIBRARIES MATCHES "mpi_")
                message(VERBOSE "CMAKE_CUDA_STANDARD_LIBRARIES already contains MPI libraries")
                list(APPEND STALE_CONFIG_LOG "CMAKE_CUDA_STANDARD_LIBRARIES pre-populated with MPI libs")
                list(APPEND STALE_CONFIG_LOG "  Found: ${CMAKE_CUDA_STANDARD_LIBRARIES}")
            endif()
        endif()
    else()
        message(DEBUG "First configure - skipping pre-populated library check")
    endif()
    
    list(POP_BACK CMAKE_MESSAGE_CONTEXT)
endif()

# ==============================================================================
# Compiler Version Checks
# ==============================================================================

message(VERBOSE "Checking compiler versions")

# GCC Version Check (for std::filesystem support)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    message(VERBOSE "Detected GNU compiler: ${CMAKE_CXX_COMPILER_VERSION}")
    
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "8.0")
        message(FATAL_ERROR 
        "\n"
        "════════════════════════════════════════════════════════════════\n"
        "ERF requires GCC 8.0+ for C++17 <filesystem> support\n"
        "Found: GCC ${CMAKE_CXX_COMPILER_VERSION}\n"
        "════════════════════════════════════════════════════════════════\n"
        "\n"
        "On Cray systems:\n"
        "  1. Load newer compiler: module load PrgEnv-gnu gcc\n"
        "  2. Verify version: CC --version\n"
        "")
    endif()
    
    message(DEBUG "GCC ${CMAKE_CXX_COMPILER_VERSION} >= 8.0 (C++17 filesystem supported)")
    
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Cray")
    message(VERBOSE "Detected Cray compiler: ${CMAKE_CXX_COMPILER_VERSION}")
    if(DEFINED ENV{PE_ENV})
        message(DEBUG "Programming Environment: $ENV{PE_ENV}")
    endif()
else()
    message(VERBOSE "Detected ${CMAKE_CXX_COMPILER_ID} compiler: ${CMAKE_CXX_COMPILER_VERSION}")
endif()

# GPU Compiler Checks
if(ERF_ENABLE_CUDA)
    message(VERBOSE "Checking CUDA compiler configuration")
    
    if(CMAKE_CUDA_COMPILER)
        message(DEBUG "CMAKE_CUDA_COMPILER: ${CMAKE_CUDA_COMPILER}")
    else()
        message(DEBUG "CMAKE_CUDA_COMPILER not set (will auto-detect)")
    endif()
    
    if(CMAKE_CUDA_FLAGS)
        message(DEBUG "CMAKE_CUDA_FLAGS: ${CMAKE_CUDA_FLAGS}")
    endif()
    
    # Detect AMReX CUDA architecture
    message(DEBUG "Detecting CUDA architecture")
    
    if(AMReX_CUDA_ARCH)
        message(VERBOSE "AMReX_CUDA_ARCH: ${AMReX_CUDA_ARCH} (user specified)")
        
    elseif(DEFINED ENV{AMREX_CUDA_ARCH})
        set(AMReX_CUDA_ARCH "$ENV{AMREX_CUDA_ARCH}" CACHE STRING "CUDA arch from AMREX_CUDA_ARCH")
        message(VERBOSE "AMReX_CUDA_ARCH: $ENV{AMREX_CUDA_ARCH} (from environment)")
        
    elseif(DEFINED ENV{CMAKE_CUDA_ARCH})
        set(ENV_CUDA_ARCH "$ENV{CMAKE_CUDA_ARCH}")
        message(DEBUG "Found CMAKE_CUDA_ARCH: ${ENV_CUDA_ARCH}")
        
        # Convert to AMReX format
        if(ENV_CUDA_ARCH MATCHES "^[0-9][0-9]$")
            string(SUBSTRING "${ENV_CUDA_ARCH}" 0 1 MAJOR)
            string(SUBSTRING "${ENV_CUDA_ARCH}" 1 1 MINOR)
            set(DETECTED_CUDA_ARCH "${MAJOR}.${MINOR}")
            message(TRACE "Converted ${ENV_CUDA_ARCH} -> ${DETECTED_CUDA_ARCH}")
        else()
            set(DETECTED_CUDA_ARCH "${ENV_CUDA_ARCH}")
        endif()
        
        set(AMReX_CUDA_ARCH "${DETECTED_CUDA_ARCH}" CACHE STRING "CUDA arch from CMAKE_CUDA_ARCH")
        message(VERBOSE "AMReX_CUDA_ARCH: ${DETECTED_CUDA_ARCH} (from CMAKE_CUDA_ARCH)")
        
    elseif(DEFINED ENV{CRAY_ACCEL_TARGET})
        set(CRAY_ACCEL_TARGET "$ENV{CRAY_ACCEL_TARGET}")
        message(VERBOSE "CRAY_ACCEL_TARGET: ${CRAY_ACCEL_TARGET}")
        
        if(CRAY_ACCEL_TARGET STREQUAL "nvidia70")
            set(AMReX_CUDA_ARCH "7.0" CACHE STRING "CUDA arch from CRAY_ACCEL_TARGET")
            message(VERBOSE "AMReX_CUDA_ARCH: 7.0 (Tesla V100)")
        elseif(CRAY_ACCEL_TARGET STREQUAL "nvidia80")
            set(AMReX_CUDA_ARCH "8.0" CACHE STRING "CUDA arch from CRAY_ACCEL_TARGET")
            message(VERBOSE "AMReX_CUDA_ARCH: 8.0 (A100)")
        elseif(CRAY_ACCEL_TARGET STREQUAL "nvidia90")
            set(AMReX_CUDA_ARCH "9.0" CACHE STRING "CUDA arch from CRAY_ACCEL_TARGET")
            message(VERBOSE "AMReX_CUDA_ARCH: 9.0 (H100)")
        else()
            message(WARNING "Unknown CRAY_ACCEL_TARGET: ${CRAY_ACCEL_TARGET}")
        endif()
    else()
        message(WARNING "AMReX_CUDA_ARCH not detected")
        message(STATUS "  For Perlmutter: module load gpu")
        message(STATUS "  Or: export CMAKE_CUDA_ARCH=80")
        message(STATUS "  Or: -DAMReX_CUDA_ARCH=8.0")
    endif()
endif()

# AMD GPU Architecture Detection (HIP)
if(AMReX_GPU_BACKEND MATCHES "HIP" OR ERF_ENABLE_HIP)
    message(VERBOSE "Checking HIP/ROCm configuration")
    
    if(AMReX_AMD_ARCH)
        message(VERBOSE "AMReX_AMD_ARCH: ${AMReX_AMD_ARCH} (user specified)")
    elseif(DEFINED ENV{AMREX_AMD_ARCH})
        set(AMReX_AMD_ARCH "$ENV{AMREX_AMD_ARCH}" CACHE STRING "AMD arch from AMREX_AMD_ARCH")
        message(VERBOSE "AMReX_AMD_ARCH: $ENV{AMREX_AMD_ARCH} (from environment)")
    elseif(DEFINED ENV{CMAKE_AMD_ARCH})
        set(AMReX_AMD_ARCH "$ENV{CMAKE_AMD_ARCH}" CACHE STRING "AMD arch from CMAKE_AMD_ARCH")
        message(VERBOSE "AMReX_AMD_ARCH: $ENV{CMAKE_AMD_ARCH} (from CMAKE_AMD_ARCH)")
    elseif(DEFINED ENV{CRAY_ACCEL_TARGET})
        set(CRAY_ACCEL_TARGET "$ENV{CRAY_ACCEL_TARGET}")
        message(VERBOSE "CRAY_ACCEL_TARGET: ${CRAY_ACCEL_TARGET}")
        
        if(CRAY_ACCEL_TARGET STREQUAL "amd_gfx90a")
            set(AMReX_AMD_ARCH "gfx90a" CACHE STRING "AMD arch from CRAY_ACCEL_TARGET")
            message(VERBOSE "AMReX_AMD_ARCH: gfx90a (MI200)")
        elseif(CRAY_ACCEL_TARGET STREQUAL "amd_gfx908")
            set(AMReX_AMD_ARCH "gfx908" CACHE STRING "AMD arch from CRAY_ACCEL_TARGET")
            message(VERBOSE "AMReX_AMD_ARCH: gfx908 (MI100)")
        elseif(CRAY_ACCEL_TARGET STREQUAL "amd_gfx942")
            set(AMReX_AMD_ARCH "gfx942" CACHE STRING "AMD arch from CRAY_ACCEL_TARGET")
            message(VERBOSE "AMReX_AMD_ARCH: gfx942 (MI300)")
        else()
            message(WARNING "Unknown CRAY_ACCEL_TARGET: ${CRAY_ACCEL_TARGET}")
        endif()
    else()
        message(WARNING "AMReX_AMD_ARCH not detected")
        message(STATUS "  For Frontier: module load craype-accel-amd-gfx90a")
        message(STATUS "  Or: export CMAKE_AMD_ARCH=gfx90a")
    endif()
endif()

# Kokkos Architecture Detection (for EKAT physics)
if(ERF_ENABLE_RRTMGP OR ERF_ENABLE_SHOC OR ERF_ENABLE_P3)
    message(VERBOSE "EKAT-based physics enabled, checking Kokkos architecture")
    
    set(KOKKOS_ARCH_SET FALSE)
    
    # Check if user already set via CMake
    if(Kokkos_ARCH_VOLTA70 OR Kokkos_ARCH_AMPERE80 OR Kokkos_ARCH_HOPPER90 OR
       Kokkos_ARCH_VEGA90A OR Kokkos_ARCH_VEGA908 OR Kokkos_ARCH_MI300A)
        set(KOKKOS_ARCH_SET TRUE)
        message(VERBOSE "Kokkos architecture already set by user")
        
    elseif(DEFINED ENV{KOKKOS_GPU_ARCH})
        set(KOKKOS_GPU_ARCH_ENV "$ENV{KOKKOS_GPU_ARCH}")
        message(VERBOSE "KOKKOS_GPU_ARCH: ${KOKKOS_GPU_ARCH_ENV}")
        
        # Map to Kokkos arch variables
        if(KOKKOS_GPU_ARCH_ENV STREQUAL "VOLTA70")
            set(Kokkos_ARCH_VOLTA70 ON CACHE BOOL "Kokkos arch from KOKKOS_GPU_ARCH")
            set(KOKKOS_ARCH_SET TRUE)
            message(DEBUG "Mapped VOLTA70 -> Kokkos_ARCH_VOLTA70")
        elseif(KOKKOS_GPU_ARCH_ENV STREQUAL "AMPERE80")
            set(Kokkos_ARCH_AMPERE80 ON CACHE BOOL "Kokkos arch from KOKKOS_GPU_ARCH")
            set(KOKKOS_ARCH_SET TRUE)
            message(DEBUG "Mapped AMPERE80 -> Kokkos_ARCH_AMPERE80")
        elseif(KOKKOS_GPU_ARCH_ENV STREQUAL "HOPPER90")
            set(Kokkos_ARCH_HOPPER90 ON CACHE BOOL "Kokkos arch from KOKKOS_GPU_ARCH")
            set(KOKKOS_ARCH_SET TRUE)
            message(DEBUG "Mapped HOPPER90 -> Kokkos_ARCH_HOPPER90")
        elseif(KOKKOS_GPU_ARCH_ENV STREQUAL "VEGA90A")
            set(Kokkos_ARCH_VEGA90A ON CACHE BOOL "Kokkos arch from KOKKOS_GPU_ARCH")
            set(KOKKOS_ARCH_SET TRUE)
            message(DEBUG "Mapped VEGA90A -> Kokkos_ARCH_VEGA90A")
        elseif(KOKKOS_GPU_ARCH_ENV STREQUAL "VEGA908")
            set(Kokkos_ARCH_VEGA908 ON CACHE BOOL "Kokkos arch from KOKKOS_GPU_ARCH")
            set(KOKKOS_ARCH_SET TRUE)
            message(DEBUG "Mapped VEGA908 -> Kokkos_ARCH_VEGA908")
        else()
            message(WARNING "Unknown KOKKOS_GPU_ARCH: ${KOKKOS_GPU_ARCH_ENV}")
        endif()
        
    elseif(DEFINED ENV{CRAY_ACCEL_TARGET})
        set(CRAY_ACCEL_TARGET "$ENV{CRAY_ACCEL_TARGET}")
        message(DEBUG "Using CRAY_ACCEL_TARGET for Kokkos: ${CRAY_ACCEL_TARGET}")
        
        # Map NVIDIA targets
        if(CRAY_ACCEL_TARGET STREQUAL "nvidia70")
            set(Kokkos_ARCH_VOLTA70 ON CACHE BOOL "Kokkos arch from CRAY_ACCEL_TARGET")
            set(KOKKOS_ARCH_SET TRUE)
            message(VERBOSE "Set Kokkos_ARCH_VOLTA70 from CRAY_ACCEL_TARGET")
        elseif(CRAY_ACCEL_TARGET STREQUAL "nvidia80")
            set(Kokkos_ARCH_AMPERE80 ON CACHE BOOL "Kokkos arch from CRAY_ACCEL_TARGET")
            set(KOKKOS_ARCH_SET TRUE)
            message(VERBOSE "Set Kokkos_ARCH_AMPERE80 from CRAY_ACCEL_TARGET")
        elseif(CRAY_ACCEL_TARGET STREQUAL "nvidia90")
            set(Kokkos_ARCH_HOPPER90 ON CACHE BOOL "Kokkos arch from CRAY_ACCEL_TARGET")
            set(KOKKOS_ARCH_SET TRUE)
            message(VERBOSE "Set Kokkos_ARCH_HOPPER90 from CRAY_ACCEL_TARGET")
        # Map AMD targets
        elseif(CRAY_ACCEL_TARGET STREQUAL "amd_gfx90a")
            set(Kokkos_ARCH_VEGA90A ON CACHE BOOL "Kokkos arch from CRAY_ACCEL_TARGET")
            set(KOKKOS_ARCH_SET TRUE)
            message(VERBOSE "Set Kokkos_ARCH_VEGA90A from CRAY_ACCEL_TARGET")
        elseif(CRAY_ACCEL_TARGET STREQUAL "amd_gfx908")
            set(Kokkos_ARCH_VEGA908 ON CACHE BOOL "Kokkos arch from CRAY_ACCEL_TARGET")
            set(KOKKOS_ARCH_SET TRUE)
            message(VERBOSE "Set Kokkos_ARCH_VEGA908 from CRAY_ACCEL_TARGET")
        elseif(CRAY_ACCEL_TARGET STREQUAL "amd_gfx942")
            set(Kokkos_ARCH_MI300A ON CACHE BOOL "Kokkos arch from CRAY_ACCEL_TARGET")
            set(KOKKOS_ARCH_SET TRUE)
            message(VERBOSE "Set Kokkos_ARCH_MI300A from CRAY_ACCEL_TARGET")
        endif()
    endif()
    
    if(NOT KOKKOS_ARCH_SET)
        message(WARNING "Kokkos architecture not detected")
        message(STATUS "  For Perlmutter: module load gpu")
        message(STATUS "  For Frontier: module load craype-accel-amd-gfx90a")
        message(STATUS "  Or: export KOKKOS_GPU_ARCH=AMPERE80")
    else()
        message(DEBUG "Note: Kokkos will set CMAKE_CUDA_ARCHITECTURES when CUDA language is enabled")
    endif()
endif()

# ==============================================================================
# Prerequisite Checks
# ==============================================================================

message(VERBOSE "Checking prerequisites")

# CMake Version Check
set(ERF_RECOMMENDED_CMAKE_VERSION "3.24.0")
if(CMAKE_VERSION VERSION_LESS ${ERF_RECOMMENDED_CMAKE_VERSION})
    message(WARNING "CMake ${CMAKE_VERSION} < recommended ${ERF_RECOMMENDED_CMAKE_VERSION}")
    message(STATUS "  Fix: module load cmake")
    message(DEBUG "Older CMake may have issues with Cray wrappers and CUDA")
else()
    message(DEBUG "CMake ${CMAKE_VERSION} >= ${ERF_RECOMMENDED_CMAKE_VERSION}")
endif()

# CUDA Toolkit Check
if(ERF_ENABLE_CUDA)
    message(VERBOSE "Checking CUDA toolkit")
    
    set(CUDA_TOOLKIT_LOADED FALSE)
    
    if(DEFINED ENV{CUDA_HOME})
        set(CUDA_TOOLKIT_LOADED TRUE)
        message(VERBOSE "  CUDA_HOME: $ENV{CUDA_HOME}")
    endif()
    
    if(NOT CUDA_TOOLKIT_LOADED AND DEFINED ENV{CUDATOOLKIT_HOME})
        set(CUDA_TOOLKIT_LOADED TRUE)
        message(VERBOSE "  CUDATOOLKIT_HOME: $ENV{CUDATOOLKIT_HOME}")
    endif()
    
    find_program(NVCC_EXECUTABLE nvcc)
    if(NVCC_EXECUTABLE)
        set(CUDA_TOOLKIT_LOADED TRUE)
        message(VERBOSE "  nvcc: ${NVCC_EXECUTABLE}")
    endif()
    
    if(NOT CUDA_TOOLKIT_LOADED)
        message(WARNING "CUDA toolkit not detected")
        message(STATUS "  Fix: module load cuda")
    endif()
else()
    message(DEBUG "CUDA not enabled, skipping toolkit check")
endif()

# NetCDF Module Check
if(ERF_ENABLE_NETCDF)
    message(VERBOSE "Checking NetCDF")
    
    if(DEFINED ENV{NETCDF_DIR})
        message(VERBOSE "  NETCDF_DIR: $ENV{NETCDF_DIR}")
    else()
        message(STATUS "  NetCDF module not detected")
        message(STATUS "  Recommended: module load cray-netcdf-hdf5parallel")
    endif()
else()
    message(DEBUG "NetCDF not enabled")
endif()

# HDF5 Module Check
if(AMReX_HDF5)
    message(VERBOSE "Checking HDF5")
    
    if(DEFINED ENV{HDF5_DIR})
        message(VERBOSE "  HDF5_DIR: $ENV{HDF5_DIR}")
    elseif(DEFINED ENV{HDF5_ROOT})
        message(VERBOSE "  HDF5_ROOT: $ENV{HDF5_ROOT}")
    else()
        message(STATUS "  HDF5 module not detected")
        message(STATUS "  Recommended: module load cray-hdf5-parallel")
    endif()
else()
    message(DEBUG "HDF5 not enabled")
endif()

# FFTW Module Check
if(ERF_ENABLE_FFT)
    message(VERBOSE "Checking FFTW")
    
    if(DEFINED ENV{FFTW_DIR})
        message(VERBOSE "  FFTW_DIR: $ENV{FFTW_DIR}")
    elseif(DEFINED ENV{CRAY_FFTW_DIR})
        message(VERBOSE "  CRAY_FFTW_DIR: $ENV{CRAY_FFTW_DIR}")
    else()
        message(STATUS "  FFTW module not detected")
        message(STATUS "  Recommended: module load cray-fftw")
    endif()
else()
    message(DEBUG "FFTW not enabled")
endif()

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

# Environment summary (DEBUG level)
message(DEBUG "Key environment variables:")
message(DEBUG "  CRAYPE_VERSION: $ENV{CRAYPE_VERSION}")
message(DEBUG "  CRAY_MPICH_DIR: $ENV{CRAY_MPICH_DIR}")
message(DEBUG "  CUDA_HOME: $ENV{CUDA_HOME}")
message(DEBUG "  NETCDF_DIR: $ENV{NETCDF_DIR}")
message(DEBUG "  HDF5_DIR: $ENV{HDF5_DIR}")
message(DEBUG "  FFTW_DIR: $ENV{FFTW_DIR}")
message(DEBUG "  MPICH_GPU_SUPPORT_ENABLED: $ENV{MPICH_GPU_SUPPORT_ENABLED}")

# ==============================================================================
# Fix 1: CUDA + EKAT -> nvcc_wrapper complications
# ==============================================================================

if(ERF_ENABLE_CUDA AND (ERF_ENABLE_RRTMGP OR ERF_ENABLE_SHOC OR ERF_ENABLE_P3))
    message(STATUS "Applying Fix 1: CUDA+EKAT nvcc_wrapper")

    message(DEBUG "Problem: nvcc_wrapper doesn't inherit Cray include paths")
    message(DEBUG "Solution: Add flags from CC --cray-print-opts=cflags")

    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} --cray-print-opts=cflags
        OUTPUT_VARIABLE CRAY_CUDA_FLAGS
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE CRAY_CUDA_FLAGS_RESULT
    )

    if(CRAY_CUDA_FLAGS_RESULT EQUAL 0 AND CRAY_CUDA_FLAGS)
        message(VERBOSE "Adding Cray flags to CMAKE_CUDA_FLAGS")
        message(DEBUG "Flags: ${CRAY_CUDA_FLAGS}")

        if(CMAKE_CUDA_FLAGS)
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} ${CRAY_CUDA_FLAGS}" CACHE STRING "" FORCE)
        else()
            set(CMAKE_CUDA_FLAGS "${CRAY_CUDA_FLAGS}" CACHE STRING "" FORCE)
        endif()
    else()
        message(WARNING "Could not retrieve Cray CUDA flags")
        message(STATUS "  Try: -DCMAKE_CUDA_FLAGS=\"\$(CC --cray-print-opts=cflags)\"")
    endif()
else()
    message(DEBUG "Fix 1 not needed (CUDA+EKAT not both enabled)")
endif()

# ==============================================================================
# Fix 2: FCOMPARE + Cray + EKAT/Kokkos -> mpi_gnu_123 not found
# ==============================================================================

if(ERF_ENABLE_FCOMPARE AND (ERF_ENABLE_RRTMGP OR ERF_ENABLE_SHOC OR ERF_ENABLE_P3))
    message(STATUS "Applying Fix 2: fcompare linker with EKAT")

    message(DEBUG "Problem: --as-needed drops required MPI libs when EKAT is enabled")
    message(DEBUG "Solution: Clean Cray libs and add --no-as-needed")

    # Query CXX compiler for libs
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} --cray-print-opts=libs
        OUTPUT_VARIABLE CRAY_LIBS_RAW
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE CRAY_LIBS_RESULT
    )

    if(CRAY_LIBS_RESULT EQUAL 0 AND CRAY_LIBS_RAW)
        message(TRACE "Raw libs from ${CMAKE_CXX_COMPILER}: ${CRAY_LIBS_RAW}")

        # Verify C and Fortran compilers return same libs
        set(OTHER_COMPILERS "${CMAKE_C_COMPILER};${CMAKE_Fortran_COMPILER}")
        set(COMPILER_NAMES "C;Fortran")
        foreach(COMPILER COMPILER_NAME IN ZIP_LISTS OTHER_COMPILERS COMPILER_NAMES)
            if(EXISTS ${COMPILER})
                execute_process(
                    COMMAND ${COMPILER} --cray-print-opts=libs
                    OUTPUT_VARIABLE OTHER_LIBS
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_QUIET
                    RESULT_VARIABLE OTHER_RESULT
                )

                if(OTHER_RESULT EQUAL 0 AND NOT "${OTHER_LIBS}" STREQUAL "${CRAY_LIBS_RAW}")
                    message(STATUS "${COMPILER_NAME} compiler libs differ from CXX (using CXX libs only)")
            	    message(VERBOSE "${COMPILER_NAME} compiler returns different libs than CXX:\n  CXX: ${CRAY_LIBS_RAW}\n  ${COMPILER_NAME}: ${OTHER_LIBS}\n  Using only CXX libs - build may fail if language-specific libs needed")
                endif()
            endif()
        endforeach()

        # Remove problematic flags
        string(REGEX REPLACE "-Wl,--as-needed," "" CRAY_LIBS_CLEAN "${CRAY_LIBS_RAW}")
        string(REGEX REPLACE ",--no-as-needed" "" CRAY_LIBS_CLEAN "${CRAY_LIBS_CLEAN}")
        string(REGEX REPLACE ",-l" " -l" CRAY_LIBS_CLEAN "${CRAY_LIBS_CLEAN}")

	# Strip GPU-related flags that Kokkos/AMReX already handle
	string(REGEX REPLACE "--offload-arch=[^ ]* *" "" CRAY_LIBS_CLEAN "${CRAY_LIBS_CLEAN}")
	string(REGEX REPLACE "--hip-link *" "" CRAY_LIBS_CLEAN "${CRAY_LIBS_CLEAN}")
	string(REGEX REPLACE "-fgpu-rdc *" "" CRAY_LIBS_CLEAN "${CRAY_LIBS_CLEAN}")
	string(REGEX REPLACE "-xhip *" "" CRAY_LIBS_CLEAN "${CRAY_LIBS_CLEAN}")

	message(VERBOSE "Adding: -Wl,--no-as-needed + cleaned libs")
        message(DEBUG "Cleaned libs: ${CRAY_LIBS_CLEAN}")

        # Check if Fix 2 already applied
        string(FIND "${CMAKE_EXE_LINKER_FLAGS}" "-Wl,--no-as-needed" already_present)
        if(already_present EQUAL -1)
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed ${CRAY_LIBS_CLEAN}"
                CACHE STRING "" FORCE)
        else()
            message(DEBUG "Cray libs already in linker flags, skipping")
        endif()
    else()
        message(WARNING "Could not retrieve Cray library paths")
        message(STATUS "  Fcompare may fail to link")
    endif()
else()
    message(DEBUG "Fix 2 not needed (fcompare disabled or EKAT not enabled)")
endif()

# ==============================================================================
# Fix 3: CUDA math libs not found
# ==============================================================================

if(ERF_ENABLE_CUDA AND DEFINED ENV{CUDA_HOME})
    set(CUDA_MATH_PATH "$ENV{CUDA_HOME}/../../math_libs/lib64")

    message(DEBUG "Checking CUDA math libs: ${CUDA_MATH_PATH}")

    if(EXISTS ${CUDA_MATH_PATH})
        message(STATUS "Applying Fix 3: CUDA math libraries")
        message(VERBOSE "Adding: ${CUDA_MATH_PATH}")

        # Only append if not already present
        list(FIND CMAKE_PREFIX_PATH "${CUDA_MATH_PATH}" already_present)
        if(already_present EQUAL -1)
            list(APPEND CMAKE_PREFIX_PATH ${CUDA_MATH_PATH})
            set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} CACHE STRING "" FORCE)
        else()
            message(DEBUG "CUDA math path already in CMAKE_PREFIX_PATH, skipping")
        endif()
    else()
        message(WARNING "CUDA math libs not found at ${CUDA_MATH_PATH}")
        message(STATUS "  Fix: module load cuda")
    endif()
else()
    message(DEBUG "Fix 3 not needed (CUDA disabled or CUDA_HOME not set)")
endif()

# ==============================================================================
# Fix 4: GPU-aware MPI with Cray GTL
# ==============================================================================

if(ERF_ENABLE_MPI AND "$ENV{MPICH_GPU_SUPPORT_ENABLED}" STREQUAL "1")
    set(APPLY_FIX4 FALSE)
    set(GPU_TYPE "")
    set(GTL_LIB "")
    set(MPI_BASE_LIB "")

    message(VERBOSE "Detecting MPI library for GPU-aware support")

    # Try pkg-config first
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        execute_process(
            COMMAND CC --cray-print-opts=pkg_config_path
            OUTPUT_VARIABLE CRAY_PKG_CONFIG_PATH
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE CC_RESULT
        )

        if(CC_RESULT EQUAL 0 AND CRAY_PKG_CONFIG_PATH)
            set(ENV{PKG_CONFIG_PATH} "${CRAY_PKG_CONFIG_PATH}:$ENV{PKG_CONFIG_PATH}")
            message(TRACE "PKG_CONFIG_PATH: ${CRAY_PKG_CONFIG_PATH}")
        endif()

        pkg_check_modules(CRAY_MPI QUIET mpich)
        if(CRAY_MPI_FOUND)
            message(DEBUG "Found mpich via pkg-config")
            foreach(lib IN LISTS CRAY_MPI_LIBRARIES CRAY_MPI_LINK_LIBRARIES)
                if(lib MATCHES "^mpi_" AND NOT lib MATCHES "mpi_gtl")
                    set(MPI_BASE_LIB "${lib}")
                    message(DEBUG "Detected MPI base: ${MPI_BASE_LIB}")
                    break()
                endif()
            endforeach()
        endif()
    endif()

    # Fallback: Search filesystem
    if(NOT MPI_BASE_LIB)
        message(DEBUG "Falling back to filesystem search")
        set(MPI_LIB_SEARCH_PATHS "")
        if(DEFINED ENV{MPICH_DIR})
            list(APPEND MPI_LIB_SEARCH_PATHS "$ENV{MPICH_DIR}/lib")
        endif()
        if(DEFINED ENV{CRAY_MPICH_DIR})
            list(APPEND MPI_LIB_SEARCH_PATHS "$ENV{CRAY_MPICH_DIR}/lib")
        endif()

        foreach(path IN LISTS MPI_LIB_SEARCH_PATHS)
            file(GLOB mpi_libs "${path}/libmpi_*.so" "${path}/libmpi_*.a")
            foreach(lib IN LISTS mpi_libs)
                get_filename_component(libname "${lib}" NAME_WE)
                string(REGEX REPLACE "^lib" "" libname "${libname}")
                if(libname MATCHES "^mpi_(gnu|cray|intel)" AND NOT MPI_BASE_LIB)
                    set(MPI_BASE_LIB "${libname}")
                    message(DEBUG "Found MPI lib: ${MPI_BASE_LIB}")
                    break()
                endif()
            endforeach()
            if(MPI_BASE_LIB)
                break()
            endif()
        endforeach()
    endif()

    # Last resort: Heuristic
    if(NOT MPI_BASE_LIB)
        if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
            set(MPI_BASE_LIB "mpi_gnu_123")
            message(WARNING "Using heuristic MPI library: ${MPI_BASE_LIB}")
        elseif(CMAKE_CXX_COMPILER_ID MATCHES "Cray")
            set(MPI_BASE_LIB "mpi_cray")
            message(WARNING "Using heuristic MPI library: ${MPI_BASE_LIB}")
        else()
            set(MPI_BASE_LIB "mpi")
        endif()
    endif()

    # Verify MPI library exists (if checking enabled)
    if(ERF_CHECK_MODULES AND MPI_BASE_LIB)
        message(DEBUG "Verifying MPI library: ${MPI_BASE_LIB}")

        find_library(MPI_BASE_VERIFY
            NAMES ${MPI_BASE_LIB}
            PATHS
                $ENV{MPICH_DIR}/lib
                $ENV{CRAY_MPICH_DIR}/lib
            NO_DEFAULT_PATH
        )

        if(MPI_BASE_VERIFY)
            message(DEBUG "Verified MPI library exists: ${MPI_BASE_VERIFY}")
        else()
            message(VERBOSE "MPI library ${MPI_BASE_LIB} not found")
            list(APPEND STALE_CONFIG_LOG "MPI library lib${MPI_BASE_LIB}.so not found")
            list(APPEND STALE_CONFIG_LOG "  Searched in: \$MPICH_DIR/lib, \$CRAY_MPICH_DIR/lib")

            # Try to suggest correct version
            if(CMAKE_CXX_COMPILER_VERSION MATCHES "^([0-9]+)\\.([0-9]+)")
                set(EXPECTED_VER "${CMAKE_MATCH_1}${CMAKE_MATCH_2}")
                list(APPEND STALE_CONFIG_LOG "  Expected based on GCC ${CMAKE_CXX_COMPILER_VERSION}: mpi_gnu_${EXPECTED_VER}")
            endif()
        endif()
        unset(MPI_BASE_VERIFY CACHE)
    endif()

    # Determine GPU type and GTL library
    if(ERF_ENABLE_CUDA)
        set(APPLY_FIX4 TRUE)
        set(GPU_TYPE "CUDA")
        set(GTL_LIB "mpi_gtl_cuda")
    elseif(AMReX_GPU_BACKEND MATCHES "HIP")
        set(APPLY_FIX4 TRUE)
        set(GPU_TYPE "HIP")
        set(GTL_LIB "mpi_gtl_hsa")
    endif()

    if(APPLY_FIX4)
        message(STATUS "Applying Fix 4: GPU-aware MPI (${GPU_TYPE})")
        message(VERBOSE "MPI base library: ${MPI_BASE_LIB}")
        message(VERBOSE "GTL library: ${GTL_LIB}")

        set(CRAY_MPI_LIBS "-l${MPI_BASE_LIB} -l${GTL_LIB}")

        # Only append if not already present (check for MPI_BASE_LIB as marker)
        if(ERF_ENABLE_CUDA)
            string(FIND "${CMAKE_CUDA_STANDARD_LIBRARIES}" "-l${MPI_BASE_LIB}" already_present)
            if(already_present EQUAL -1)
                message(DEBUG "Adding to CMAKE_CUDA_STANDARD_LIBRARIES: ${CRAY_MPI_LIBS}")
                set(CMAKE_CUDA_STANDARD_LIBRARIES "${CMAKE_CUDA_STANDARD_LIBRARIES} ${CRAY_MPI_LIBS}"
                    CACHE STRING "" FORCE)
            else()
                message(DEBUG "CUDA libraries already contain MPI libs, skipping")
            endif()
        elseif(ERF_ENABLE_HIP)
            string(FIND "${CMAKE_HIP_STANDARD_LIBRARIES}" "-l${MPI_BASE_LIB}" already_present)
            if(already_present EQUAL -1)
                message(DEBUG "Adding to CMAKE_HIP_STANDARD_LIBRARIES: ${CRAY_MPI_LIBS}")
                set(CMAKE_HIP_STANDARD_LIBRARIES "${CMAKE_HIP_STANDARD_LIBRARIES} ${CRAY_MPI_LIBS}"
                    CACHE STRING "" FORCE)
            else()
                message(DEBUG "HIP libraries already contain MPI libs, skipping")
            endif()
        endif()

        string(FIND "${CMAKE_CXX_STANDARD_LIBRARIES}" "-l${MPI_BASE_LIB}" already_present)
        if(already_present EQUAL -1)
            message(DEBUG "Adding to CMAKE_CXX_STANDARD_LIBRARIES: ${CRAY_MPI_LIBS}")
            set(CMAKE_CXX_STANDARD_LIBRARIES "${CMAKE_CXX_STANDARD_LIBRARIES} ${CRAY_MPI_LIBS}"
                 CACHE STRING "" FORCE)
        else()
            message(DEBUG "CXX libraries already contain MPI libs, skipping")
        endif()
    endif()
else()
    message(DEBUG "Fix 4 not needed (GPU+MPI not enabled or GPU support not enabled)")
endif()

# ==============================================================================
# Fix 5-6: NetCDF with cray-netcdf-hdf5parallel
# ==============================================================================

if(ERF_ENABLE_NETCDF)
    message(STATUS "Applying Fix 5-6: NetCDF configuration")

    message(DEBUG "Setting up pkg-config path for NetCDF/MPI")

    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} --cray-print-opts=PKG_CONFIG_PATH
        OUTPUT_VARIABLE CRAY_PKG_CONFIG_PATH
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE PKG_RESULT
    )

    if(PKG_RESULT EQUAL 0 AND CRAY_PKG_CONFIG_PATH)
        message(VERBOSE "PKG_CONFIG_PATH from Cray wrapper")
        message(DEBUG "  ${CRAY_PKG_CONFIG_PATH}")

        if(DEFINED ENV{PKG_CONFIG_PATH})
            set(ENV{PKG_CONFIG_PATH} "${CRAY_PKG_CONFIG_PATH}:$ENV{PKG_CONFIG_PATH}")
        else()
            set(ENV{PKG_CONFIG_PATH} "${CRAY_PKG_CONFIG_PATH}")
        endif()
    endif()

    # Add NetCDF/HDF5 to search paths
    if(DEFINED ENV{NETCDF_DIR})
        list(APPEND CMAKE_PREFIX_PATH $ENV{NETCDF_DIR})
        message(VERBOSE "Added NETCDF_DIR to search: $ENV{NETCDF_DIR}")
    endif()

    if(DEFINED ENV{HDF5_DIR})
        list(APPEND CMAKE_PREFIX_PATH $ENV{HDF5_DIR})
        message(VERBOSE "Added HDF5_DIR to search: $ENV{HDF5_DIR}")
    endif()
else()
    message(DEBUG "Fix 5-6 not needed (NetCDF disabled)")
endif()

# ==============================================================================
# Fix 7: HDF5 parallel detection for HIP builds
# ==============================================================================

if(AMReX_GPU_BACKEND MATCHES "HIP" AND AMReX_HDF5)
    message(STATUS "Applying Fix 7: HDF5 for HIP")

    message(DEBUG "Configuring HDF5 hints for HIP build")

    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        execute_process(
            COMMAND CC --cray-print-opts=pkg_config_path
            OUTPUT_VARIABLE CRAY_PKG_CONFIG_PATH
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE CC_RESULT
        )

        if(CC_RESULT EQUAL 0 AND CRAY_PKG_CONFIG_PATH)
            set(ENV{PKG_CONFIG_PATH} "${CRAY_PKG_CONFIG_PATH}:$ENV{PKG_CONFIG_PATH}")
        endif()

        pkg_check_modules(PC_HDF5 QUIET hdf5)
        if(PC_HDF5_FOUND)
            message(VERBOSE "Found HDF5 via pkg-config: ${PC_HDF5_PREFIX}")

            set(HDF5_ROOT "${PC_HDF5_PREFIX}" CACHE PATH "HDF5 root from pkg-config")
            set(HDF5_PREFER_PARALLEL ON CACHE BOOL "Prefer parallel HDF5")
            set(HDF5_IS_PARALLEL TRUE CACHE BOOL "HDF5 is parallel")

            list(APPEND CMAKE_PREFIX_PATH "${PC_HDF5_PREFIX}")

            message(DEBUG "Set HDF5_ROOT: ${PC_HDF5_PREFIX}")
            message(DEBUG "Set HDF5_PREFER_PARALLEL: ON")
        else()
            message(WARNING "pkg-config could not find HDF5")
        endif()
    else()
        message(WARNING "PkgConfig not found, cannot auto-configure HDF5")
    endif()
else()
    message(DEBUG "Fix 7 not needed (not HIP+HDF5)")
endif()

# ==============================================================================
# Summary
# ==============================================================================

message(STATUS "Cray configuration complete")

# Track which fixes were applied
set(FIX1_ACTIVE OFF)
set(FIX2_ACTIVE OFF)
set(FIX3_ACTIVE OFF)
set(FIX4_ACTIVE OFF)
set(FIX56_ACTIVE OFF)
set(FIX7_ACTIVE OFF)

# Fix 1: CUDA + EKAT
if(ERF_ENABLE_CUDA AND (ERF_ENABLE_RRTMGP OR ERF_ENABLE_SHOC OR ERF_ENABLE_P3) AND CRAY_CUDA_FLAGS)
    set(FIX1_ACTIVE ON)
endif()

# Fix 2: fcompare + EKAT
if(ERF_ENABLE_FCOMPARE AND (ERF_ENABLE_RRTMGP OR ERF_ENABLE_SHOC OR ERF_ENABLE_P3) AND CRAY_LIBS_CLEAN)
    set(FIX2_ACTIVE ON)
endif()

# Fix 3: CUDA math libs
if(ERF_ENABLE_CUDA AND DEFINED ENV{CUDA_HOME})
    set(CUDA_MATH_CHECK "$ENV{CUDA_HOME}/../../math_libs/lib64")
    if(EXISTS ${CUDA_MATH_CHECK})
        set(FIX3_ACTIVE ON)
    endif()
endif()

# Fix 4: GPU-aware MPI
if(APPLY_FIX4)
    set(FIX4_ACTIVE ON)
endif()

# Fix 5-6: NetCDF
if(ERF_ENABLE_NETCDF)
    set(FIX56_ACTIVE ON)
endif()

# Fix 7: HDF5 for HIP
if(AMReX_GPU_BACKEND MATCHES "HIP" AND AMReX_HDF5)
    set(FIX7_ACTIVE ON)
endif()

# Show summary at VERBOSE level
message(VERBOSE "Applied fixes:")
if(FIX1_ACTIVE)
    message(VERBOSE "  Fix 1 (CUDA+EKAT): ACTIVE")
endif()
if(FIX2_ACTIVE)
    message(VERBOSE "  Fix 2 (fcompare+EKAT): ACTIVE")
endif()
if(FIX3_ACTIVE)
    message(VERBOSE "  Fix 3 (CUDA math): ACTIVE")
endif()
if(FIX4_ACTIVE)
    message(VERBOSE "  Fix 4 (GPU-aware MPI): ACTIVE")
endif()
if(FIX56_ACTIVE)
    message(VERBOSE "  Fix 5-6 (NetCDF): ACTIVE")
endif()
if(FIX7_ACTIVE)
    message(VERBOSE "  Fix 7 (HDF5+HIP): ACTIVE")
endif()

# Command-line equivalents (DEBUG level)
message(DEBUG "Command-line equivalents for active fixes:")
message(DEBUG "=====================================================================")

if(FIX1_ACTIVE)
    message(DEBUG "Fix 1 (CUDA+EKAT):")
    message(DEBUG "  -DCMAKE_CUDA_FLAGS=\"\$(CC --cray-print-opts=cflags)\"")
    message(DEBUG "")
endif()

if(FIX2_ACTIVE)
    message(DEBUG "Fix 2 (fcompare+EKAT):")
    message(DEBUG "  CRAY_LIBS=\"\$(CC --cray-print-opts=libs | sed 's/-Wl,--as-needed,//g; s/,--no-as-needed//g; s/,-l/ -l/g')\"")
    message(DEBUG "  -DCMAKE_EXE_LINKER_FLAGS=\"-Wl,--no-as-needed \$CRAY_LIBS\"")
    message(DEBUG "")
endif()

if(FIX3_ACTIVE)
    message(DEBUG "Fix 3 (CUDA math):")
    message(DEBUG "  -DCMAKE_PREFIX_PATH=\"\$CUDA_HOME/../../math_libs/lib64\"")
    message(DEBUG "")
endif()

if(FIX4_ACTIVE)
    message(DEBUG "Fix 4 (GPU-aware MPI):")
    message(DEBUG "  export MPICH_GPU_SUPPORT_ENABLED=1")
    if(ERF_ENABLE_CUDA)
        message(DEBUG "  -DCMAKE_CUDA_STANDARD_LIBRARIES=\"-l${MPI_BASE_LIB} -l${GTL_LIB}\"")
    else()
        message(DEBUG "  -DCMAKE_HIP_STANDARD_LIBRARIES=\"-l${MPI_BASE_LIB} -l${GTL_LIB}\"")
    endif()
    message(DEBUG "  -DCMAKE_CXX_STANDARD_LIBRARIES=\"-l${MPI_BASE_LIB} -l${GTL_LIB}\"")
    message(DEBUG "")
endif()

if(FIX56_ACTIVE)
    message(DEBUG "Fix 5-6 (NetCDF):")
    message(DEBUG "  export PKG_CONFIG_PATH=\"\$(CC --cray-print-opts=PKG_CONFIG_PATH):\$PKG_CONFIG_PATH\"")
    if(DEFINED ENV{NETCDF_DIR})
        message(DEBUG "  -DCMAKE_PREFIX_PATH=\"\$NETCDF_DIR\"")
    endif()
    if(DEFINED ENV{HDF5_DIR})
        message(DEBUG "  -DCMAKE_PREFIX_PATH=\"\$CMAKE_PREFIX_PATH:\$HDF5_DIR\"")
    endif()
    message(DEBUG "")
endif()

if(FIX7_ACTIVE)
    message(DEBUG "Fix 7 (HDF5+HIP):")
    message(DEBUG "  -DHDF5_ROOT=\$(pkg-config --variable=prefix hdf5)")
    message(DEBUG "  -DHDF5_PREFER_PARALLEL=ON")
    message(DEBUG "  -DHDF5_IS_PARALLEL=TRUE")
    message(DEBUG "")
endif()

if(FIX1_ACTIVE OR FIX2_ACTIVE OR FIX3_ACTIVE OR FIX4_ACTIVE OR FIX56_ACTIVE OR FIX7_ACTIVE)
    message(DEBUG "Complete manual equivalent (all active fixes):")
    message(DEBUG "=====================================================================")
    if(FIX1_ACTIVE)
        message(DEBUG "  -DCMAKE_CUDA_FLAGS=\"\$(CC --cray-print-opts=cflags)\" \\")
    endif()
    if(FIX2_ACTIVE)
        message(DEBUG "  -DCMAKE_EXE_LINKER_FLAGS=\"-Wl,--no-as-needed ${CRAY_LIBS_CLEAN}\" \\")
    endif()
    if(FIX3_ACTIVE)
        message(DEBUG "  -DCMAKE_PREFIX_PATH=\"\$CUDA_HOME/../../math_libs/lib64\" \\")
    endif()
    if(FIX4_ACTIVE)
        if(ERF_ENABLE_CUDA)
            message(DEBUG "  -DCMAKE_CUDA_STANDARD_LIBRARIES=\"-l${MPI_BASE_LIB} -l${GTL_LIB}\" \\")
        else()
            message(DEBUG "  -DCMAKE_HIP_STANDARD_LIBRARIES=\"-l${MPI_BASE_LIB} -l${GTL_LIB}\" \\")
        endif()
        message(DEBUG "  -DCMAKE_CXX_STANDARD_LIBRARIES=\"-l${MPI_BASE_LIB} -l${GTL_LIB}\" \\")
    endif()
    message(DEBUG "")
endif()

# ==============================================================================
# Generate Concise Config File
# ==============================================================================

set(CRAY_CONFIG_FILE "${CMAKE_BINARY_DIR}/cray_detected_config.cmake")

file(WRITE ${CRAY_CONFIG_FILE}
"# ==============================================================================
# Auto-detected Cray Configuration
# Generated: ${CMAKE_CURRENT_LIST_FILE}
# Date: ${CMAKE_TIMESTAMP}
# ==============================================================================
# This file shows the settings auto-detected by CrayDetection.cmake
# You can use this as a starting point for a manual config file.
#
# To use this config manually:
#
#   From build directory:
#     cmake -C ${CRAY_CONFIG_FILE} ${CMAKE_SOURCE_DIR}
#
#   From source directory:
#     cmake -C ${CRAY_CONFIG_FILE} -B ${CMAKE_BINARY_DIR}
#
# ==============================================================================

")

# System info
file(APPEND ${CRAY_CONFIG_FILE} "
# System Detection
set(ERF_ON_CRAY TRUE CACHE BOOL \"Detected Cray system\")
set(CRAYPE_VERSION \"$ENV{CRAYPE_VERSION}\" CACHE STRING \"Cray PE version\")
")

# Compiler info
file(APPEND ${CRAY_CONFIG_FILE} "
# Compiler Configuration
set(CMAKE_C_COMPILER \"${CMAKE_C_COMPILER}\" CACHE FILEPATH \"\")
set(CMAKE_CXX_COMPILER \"${CMAKE_CXX_COMPILER}\" CACHE FILEPATH \"\")
set(CMAKE_CXX_COMPILER_ID \"${CMAKE_CXX_COMPILER_ID}\" CACHE STRING \"\")
set(CMAKE_CXX_COMPILER_VERSION \"${CMAKE_CXX_COMPILER_VERSION}\" CACHE STRING \"\")
")

# GPU architectures
if(ERF_ENABLE_CUDA AND AMReX_CUDA_ARCH)
    file(APPEND ${CRAY_CONFIG_FILE} "
# CUDA Configuration
set(AMReX_CUDA_ARCH \"${AMReX_CUDA_ARCH}\" CACHE STRING \"Auto-detected\")
")
endif()

if(AMReX_AMD_ARCH)
    file(APPEND ${CRAY_CONFIG_FILE} "
# HIP Configuration
set(AMReX_AMD_ARCH \"${AMReX_AMD_ARCH}\" CACHE STRING \"Auto-detected\")
")
endif()

if(KOKKOS_ARCH_SET)
    file(APPEND ${CRAY_CONFIG_FILE} "
# Kokkos Architecture
")
    foreach(arch IN ITEMS VOLTA70 AMPERE80 HOPPER90 VEGA90A VEGA908 MI300A)
        if(Kokkos_ARCH_${arch})
            file(APPEND ${CRAY_CONFIG_FILE} "set(Kokkos_ARCH_${arch} ON CACHE BOOL \"Auto-detected\")\n")
        endif()
    endforeach()
endif()

# Applied fixes
file(APPEND ${CRAY_CONFIG_FILE} "
# Applied Fixes
")

if(FIX1_ACTIVE)
    file(APPEND ${CRAY_CONFIG_FILE} "
# Fix 1: CUDA+EKAT nvcc_wrapper flags
set(CMAKE_CUDA_FLAGS \"${CMAKE_CUDA_FLAGS}\" CACHE STRING \"\")
")
endif()

if(FIX2_ACTIVE)
    file(APPEND ${CRAY_CONFIG_FILE} "
# Fix 2: fcompare+EKAT linker flags
set(CMAKE_EXE_LINKER_FLAGS \"${CMAKE_EXE_LINKER_FLAGS}\" CACHE STRING \"\")
")
endif()

if(FIX3_ACTIVE)
    file(APPEND ${CRAY_CONFIG_FILE} "
# Fix 3: CUDA math libraries path
list(APPEND CMAKE_PREFIX_PATH \"${CUDA_MATH_PATH}\")
")
endif()

if(FIX4_ACTIVE)
    file(APPEND ${CRAY_CONFIG_FILE} "
# Fix 4: GPU-aware MPI (${GPU_TYPE})
set(CMAKE_CXX_STANDARD_LIBRARIES \"${CMAKE_CXX_STANDARD_LIBRARIES}\" CACHE STRING \"\")
")
    if(ERF_ENABLE_CUDA)
        file(APPEND ${CRAY_CONFIG_FILE} "set(CMAKE_CUDA_STANDARD_LIBRARIES \"${CMAKE_CUDA_STANDARD_LIBRARIES}\" CACHE STRING \"\")\n")
    endif()
endif()

if(FIX56_ACTIVE)
    file(APPEND ${CRAY_CONFIG_FILE} "
# Fix 5-6: NetCDF/HDF5 paths
set(ENV{PKG_CONFIG_PATH} \"$ENV{PKG_CONFIG_PATH}\")
")
    if(DEFINED ENV{NETCDF_DIR})
        file(APPEND ${CRAY_CONFIG_FILE} "list(APPEND CMAKE_PREFIX_PATH \"$ENV{NETCDF_DIR}\")\n")
    endif()
endif()

if(FIX7_ACTIVE)
    file(APPEND ${CRAY_CONFIG_FILE} "
# Fix 7: HDF5 parallel for HIP
set(HDF5_ROOT \"${HDF5_ROOT}\" CACHE PATH \"\")
set(HDF5_PREFER_PARALLEL ON CACHE BOOL \"\")
set(HDF5_IS_PARALLEL TRUE CACHE BOOL \"\")
")
endif()

message(STATUS "Generated config: ${CRAY_CONFIG_FILE}")

# Add a target to display it
add_custom_target(show-cray-config
    COMMAND ${CMAKE_COMMAND} -E echo "Displaying auto-detected Cray configuration"
    COMMAND ${CMAKE_COMMAND} -E echo "==================================================================="
    COMMAND ${CMAKE_COMMAND} -E echo "Auto-detected Cray Configuration:"
    COMMAND ${CMAKE_COMMAND} -E echo "==================================================================="
    COMMAND ${CMAKE_COMMAND} -E cat ${CRAY_CONFIG_FILE}
    COMMAND ${CMAKE_COMMAND} -E echo ""
    COMMAND ${CMAKE_COMMAND} -E echo "To use this config manually:"
    COMMAND ${CMAKE_COMMAND} -E echo "  cmake -C ${CRAY_CONFIG_FILE} ${CMAKE_SOURCE_DIR}"
    COMMAND ${CMAKE_COMMAND} -E echo "==================================================================="
    DEPENDS ${CRAY_CONFIG_FILE}
    COMMENT "Showing Cray configuration from ${CRAY_CONFIG_FILE}"
)

# ==============================================================================
# Display Configuration Verification Results
# ==============================================================================

if(ERF_CHECK_MODULES AND STALE_CONFIG_LOG)
    message(STATUS "")
    message(STATUS "====================================================================")
    message(STATUS "STALE CONFIGURATION DETECTED")
    message(STATUS "====================================================================")
    message(STATUS "")
    message(STATUS "Configuration issues found:")
    foreach(issue ${STALE_CONFIG_LOG})
        message(STATUS "  ${issue}")
    endforeach()
    message(STATUS "")
    message(STATUS "This usually happens when:")
    message(STATUS "  - You changed which modules are loaded")
    message(STATUS "  - You switched compiler versions")
    message(STATUS "  - CMake cache contains old settings")
    message(STATUS "")
    message(STATUS "To resolve, clean your build:")
    message(STATUS "")
    message(STATUS "  Recommended:")
    message(STATUS "    cmake --build . --target distclean")
    message(STATUS "")
    message(STATUS "  Or manually:")
    message(STATUS "    rm -rf CMakeCache.txt CMakeFiles/ cray_detected_config.cmake")
    message(STATUS "")
    message(STATUS "  Then reconfigure:")
    message(STATUS "    cmake ..")
    message(STATUS "")
    message(STATUS "To disable this check:")
    message(STATUS "    cmake -DERF_CHECK_MODULES=OFF ..")
    message(STATUS "")
    message(STATUS "====================================================================")
    message(STATUS "")

    # Make it a hard error instead of warning
    message(FATAL_ERROR "Stale configuration detected - clean build required")
endif()

message(DEBUG "=====================================================================")
message(DEBUG "To disable auto-fixes: -DERF_ENABLE_CRAY_AUTO_FIXES=OFF")
message(DEBUG "For verbose output: cmake --log-level=VERBOSE ..")
message(DEBUG "For debug output: cmake --log-level=DEBUG ..")

# Pop Cray context
list(POP_BACK CMAKE_MESSAGE_CONTEXT)
