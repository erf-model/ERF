# ==============================================================================
# Cray System Auto-Detection and Workarounds
# ==============================================================================
# This module detects Cray systems and automatically applies workarounds for
# common build issues. Each fix corresponds to a checklist item.
#
# Options:
#   -DERF_DISABLE_CRAY_AUTO_FIXES=ON   : Disable automatic Cray system fixes
#   -DERF_VERBOSE_CRAY_FIXES=ON        : Show detailed info for each fix
# ==============================================================================

option(ERF_DISABLE_CRAY_AUTO_FIXES "Disable automatic Cray system fixes" OFF)
option(ERF_VERBOSE_CRAY_FIXES "Show verbose output for Cray fixes" OFF)

# Helper macro for verbose messages
macro(erf_cray_verbose)
    if(ERF_VERBOSE_CRAY_FIXES)
        message(STATUS "  [VERBOSE] ${ARGN}")
    endif()
endmacro()

if(ERF_DISABLE_CRAY_AUTO_FIXES)
    message(STATUS "ERF: Cray auto-fixes disabled by user")
    return()
endif()

# ==============================================================================
# Detect Cray Environment
# ==============================================================================

set(ERF_ON_CRAY FALSE)

erf_cray_verbose("Checking for Cray environment...")

# Check for Cray compiler wrappers
if(CMAKE_C_COMPILER MATCHES ".*cc$" AND 
   CMAKE_CXX_COMPILER MATCHES ".*CC$" AND
   DEFINED ENV{CRAY_MPICH_DIR})
    set(ERF_ON_CRAY TRUE)
    message(STATUS "ERF: Detected Cray system")
    message(STATUS "  CMAKE_C_COMPILER   = ${CMAKE_C_COMPILER}")
    message(STATUS "  CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER}")
    message(STATUS "  CRAY_MPICH_DIR     = $ENV{CRAY_MPICH_DIR}")
    erf_cray_verbose("Detection method: Cray compiler wrappers (cc, CC) + CRAY_MPICH_DIR")
endif()

# Additional check for Cray environment variables
if(DEFINED ENV{CRAYPE_VERSION})
    set(ERF_ON_CRAY TRUE)
    message(STATUS "ERF: Detected Cray Programming Environment")
    message(STATUS "  CRAYPE_VERSION = $ENV{CRAYPE_VERSION}")
    erf_cray_verbose("Detection method: CRAYPE_VERSION environment variable")
endif()

if(NOT ERF_ON_CRAY)
    message(STATUS "ERF: Not on a Cray system, skipping Cray-specific fixes")
    erf_cray_verbose("CMAKE_C_COMPILER = ${CMAKE_C_COMPILER}")
    erf_cray_verbose("CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER}")
    erf_cray_verbose("CRAY_MPICH_DIR = $ENV{CRAY_MPICH_DIR}")
    erf_cray_verbose("CRAYPE_VERSION = $ENV{CRAYPE_VERSION}")
    return()
endif()

# ==============================================================================
# Compiler Version Checks
# ==============================================================================

message(STATUS "ERF: Checking compiler versions...")

# -----------------------------------------------------------------------------
# GCC Version Check (for std::filesystem support)
# -----------------------------------------------------------------------------
# ERF uses C++17 <filesystem> which requires GCC 8.0+
# Older GCC versions will fail with "fatal error: filesystem: No such file"

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    message(STATUS "  Detected GNU C++ compiler version: ${CMAKE_CXX_COMPILER_VERSION}")
    
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "8.0")
        message(FATAL_ERROR 
        "\n"
        "════════════════════════════════════════════════════════════════\n"
        "ERF requires GCC 8.0+ for C++17 <filesystem> support\n"
        "Found: GCC ${CMAKE_CXX_COMPILER_VERSION}\n"
        "════════════════════════════════════════════════════════════════\n"
        "\n"
        "On Cray systems, fix by using the Cray wrapper with a modern compiler:\n"
        "  1. Load a newer compiler module:\n"
        "       module load PrgEnv-gnu\n"
        "       module load gcc\n"
        "\n"
        "  2. Set compiler explicitly:\n"
        "       -DCMAKE_CXX_COMPILER=\$(which CC)\n"
        "     Or set environment variable:\n"
        "       export CXX=\$(which CC)\n"
        "\n"
        "  3. Verify compiler version:\n"
        "       CC --version\n"
        "")
    else()
        message(STATUS "  GCC version ${CMAKE_CXX_COMPILER_VERSION} >= 8.0")
        erf_cray_verbose("GCC version sufficient for C++17 <filesystem>")
    endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Cray")
    message(STATUS "  Detected Cray C++ compiler version: ${CMAKE_CXX_COMPILER_VERSION}")
    erf_cray_verbose("Cray compiler wrappers detected")
    
    # Cray wrappers forward to underlying compiler - check what's loaded
    if(DEFINED ENV{PE_ENV})
        message(STATUS "  Programming Environment: $ENV{PE_ENV}")
        erf_cray_verbose("PE_ENV = $ENV{PE_ENV}")
    endif()
else()
    message(STATUS "  Detected C++ compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
endif()

# -----------------------------------------------------------------------------
# GPU Compiler Checks (for CUDA builds)
# -----------------------------------------------------------------------------
# Kokkos and EKAT read CMAKE_CUDA_COMPILER and CMAKE_CUDA_FLAGS
# We need to ensure these are set correctly for Cray systems

if(ERF_ENABLE_CUDA)
    message(STATUS "")
    message(STATUS "ERF: Checking GPU compiler configuration...")
    
    # Check if CMAKE_CUDA_COMPILER is set
    if(CMAKE_CUDA_COMPILER)
        message(STATUS "  CMAKE_CUDA_COMPILER = ${CMAKE_CUDA_COMPILER}")
        erf_cray_verbose("CUDA compiler explicitly set by user or CMake")
    else()
        message(STATUS "  CMAKE_CUDA_COMPILER not set (CMake will auto-detect)")
        erf_cray_verbose("CMake will search for nvcc in PATH")
    endif()
    
    # Check if CMAKE_CUDA_FLAGS has been set
    if(CMAKE_CUDA_FLAGS)
        message(STATUS "  CMAKE_CUDA_FLAGS = ${CMAKE_CUDA_FLAGS}")
        erf_cray_verbose("CUDA flags explicitly set by user")
    else()
        message(STATUS "  CMAKE_CUDA_FLAGS not set (will be auto-configured)")
        erf_cray_verbose("Cray-specific CUDA flags will be added by Fix 1 if needed")
    endif()
    
    # -------------------------------------------------------------------------
    # Detect AMReX CUDA architecture
    # Priority: CMake var > AMREX_CUDA_ARCH env > CMAKE_CUDA_ARCH env > CRAY_ACCEL_TARGET
    # -------------------------------------------------------------------------
    
    if(AMReX_CUDA_ARCH)
        message(STATUS "  AMReX_CUDA_ARCH = ${AMReX_CUDA_ARCH} (user specified)")
        erf_cray_verbose("AMReX CUDA arch set via CMake variable")
        
    elseif(DEFINED ENV{AMREX_CUDA_ARCH})
        set(AMReX_CUDA_ARCH "$ENV{AMREX_CUDA_ARCH}" CACHE STRING "CUDA arch from AMREX_CUDA_ARCH")
        message(STATUS "  AMReX_CUDA_ARCH = $ENV{AMREX_CUDA_ARCH} (from AMREX_CUDA_ARCH)")
        erf_cray_verbose("AMReX CUDA arch from AMREX_CUDA_ARCH environment variable")
        
    elseif(DEFINED ENV{CMAKE_CUDA_ARCH})
        # Common in build scripts: CMAKE_CUDA_ARCH="80"
        set(ENV_CUDA_ARCH "$ENV{CMAKE_CUDA_ARCH}")
        
        # Convert to AMReX format (add decimal point if needed)
        if(ENV_CUDA_ARCH MATCHES "^[0-9][0-9]$")
            # Two-digit format: 70, 80, 90 -> 7.0, 8.0, 9.0
            string(SUBSTRING "${ENV_CUDA_ARCH}" 0 1 MAJOR)
            string(SUBSTRING "${ENV_CUDA_ARCH}" 1 1 MINOR)
            set(DETECTED_CUDA_ARCH "${MAJOR}.${MINOR}")
        else()
            # Already in decimal format or other format
            set(DETECTED_CUDA_ARCH "${ENV_CUDA_ARCH}")
        endif()
        
        set(AMReX_CUDA_ARCH "${DETECTED_CUDA_ARCH}" CACHE STRING "CUDA arch from CMAKE_CUDA_ARCH")
        message(STATUS "  AMReX_CUDA_ARCH = ${DETECTED_CUDA_ARCH} (from CMAKE_CUDA_ARCH=${ENV_CUDA_ARCH})")
        erf_cray_verbose("Converted CMAKE_CUDA_ARCH=${ENV_CUDA_ARCH} -> AMReX_CUDA_ARCH=${DETECTED_CUDA_ARCH}")
        
    elseif(DEFINED ENV{CRAY_ACCEL_TARGET})
        # Auto-detect from Cray accelerator module (set by 'module load gpu')
        set(CRAY_ACCEL_TARGET "$ENV{CRAY_ACCEL_TARGET}")
        message(STATUS "  Detected CRAY_ACCEL_TARGET = ${CRAY_ACCEL_TARGET}")
        
        if(CRAY_ACCEL_TARGET STREQUAL "nvidia70")
            set(AMReX_CUDA_ARCH "7.0" CACHE STRING "CUDA arch from CRAY_ACCEL_TARGET")
            message(STATUS "  AMReX_CUDA_ARCH = 7.0 (Tesla V100 from CRAY_ACCEL_TARGET)")
        elseif(CRAY_ACCEL_TARGET STREQUAL "nvidia80")
            set(AMReX_CUDA_ARCH "8.0" CACHE STRING "CUDA arch from CRAY_ACCEL_TARGET")
            message(STATUS "  AMReX_CUDA_ARCH = 8.0 (A100 from CRAY_ACCEL_TARGET)")
        elseif(CRAY_ACCEL_TARGET STREQUAL "nvidia90")
            set(AMReX_CUDA_ARCH "9.0" CACHE STRING "CUDA arch from CRAY_ACCEL_TARGET")
            message(STATUS "  AMReX_CUDA_ARCH = 9.0 (H100 from CRAY_ACCEL_TARGET)")
        else()
            message(WARNING "ERF: Unknown CRAY_ACCEL_TARGET = ${CRAY_ACCEL_TARGET}")
        endif()
        erf_cray_verbose("AMReX CUDA arch from CRAY_ACCEL_TARGET module variable")
    else()
        message(WARNING "")
        message(WARNING "ERF: AMReX_CUDA_ARCH not detected")
        message(WARNING "  For Perlmutter: module load gpu")
        message(WARNING "  Or set: export CMAKE_CUDA_ARCH=80")
        message(WARNING "  Or set: -DAMReX_CUDA_ARCH=8.0")
        message(WARNING "")
    endif()
    
    # -------------------------------------------------------------------------
    # Detect Kokkos architecture (for EKAT builds)
    # Priority: CMake var > KOKKOS_GPU_ARCH env > CRAY_ACCEL_TARGET
    # -------------------------------------------------------------------------
    
    if(ERF_ENABLE_RRTMGP OR ERF_ENABLE_SHOC OR ERF_ENABLE_P3)
        message(STATUS "")
        message(STATUS "  EKAT-based physics enabled, checking Kokkos architecture...")
        
        # Check if user already set Kokkos_ARCH_* via CMake
        set(KOKKOS_ARCH_SET FALSE)
        if(Kokkos_ARCH_VOLTA70 OR Kokkos_ARCH_AMPERE80 OR Kokkos_ARCH_HOPPER90)
            set(KOKKOS_ARCH_SET TRUE)
            message(STATUS "    Kokkos_ARCH_* already set by user")
            erf_cray_verbose("User specified Kokkos architecture via CMake variable")
            
        elseif(DEFINED ENV{KOKKOS_GPU_ARCH})
            # Detect from KOKKOS_GPU_ARCH environment variable (build scripts)
            set(KOKKOS_GPU_ARCH_ENV "$ENV{KOKKOS_GPU_ARCH}")
            message(STATUS "    Detected KOKKOS_GPU_ARCH = ${KOKKOS_GPU_ARCH_ENV}")
            
            # Map to Kokkos_ARCH_* CMake variable
            if(KOKKOS_GPU_ARCH_ENV STREQUAL "VOLTA70")
                set(Kokkos_ARCH_VOLTA70 ON CACHE BOOL "Kokkos arch from KOKKOS_GPU_ARCH")
                message(STATUS "    Set Kokkos_ARCH_VOLTA70 = ON")
                set(KOKKOS_ARCH_SET TRUE)
                erf_cray_verbose("Mapped KOKKOS_GPU_ARCH=VOLTA70 -> Kokkos_ARCH_VOLTA70=ON")
                
            elseif(KOKKOS_GPU_ARCH_ENV STREQUAL "AMPERE80")
                set(Kokkos_ARCH_AMPERE80 ON CACHE BOOL "Kokkos arch from KOKKOS_GPU_ARCH")
                message(STATUS "    Set Kokkos_ARCH_AMPERE80 = ON")
                set(KOKKOS_ARCH_SET TRUE)
                erf_cray_verbose("Mapped KOKKOS_GPU_ARCH=AMPERE80 -> Kokkos_ARCH_AMPERE80=ON")
                
            elseif(KOKKOS_GPU_ARCH_ENV STREQUAL "HOPPER90")
                set(Kokkos_ARCH_HOPPER90 ON CACHE BOOL "Kokkos arch from KOKKOS_GPU_ARCH")
                message(STATUS "    Set Kokkos_ARCH_HOPPER90 = ON")
                set(KOKKOS_ARCH_SET TRUE)
                erf_cray_verbose("Mapped KOKKOS_GPU_ARCH=HOPPER90 -> Kokkos_ARCH_HOPPER90=ON")
                
            else()
                message(WARNING "ERF: Unknown KOKKOS_GPU_ARCH = ${KOKKOS_GPU_ARCH_ENV}")
                message(WARNING "  Expected: VOLTA70, AMPERE80, or HOPPER90")
            endif()
            
        elseif(DEFINED ENV{CRAY_ACCEL_TARGET})
            # Fall back to CRAY_ACCEL_TARGET (set by 'module load gpu')
            set(CRAY_ACCEL_TARGET "$ENV{CRAY_ACCEL_TARGET}")
            
            if(CRAY_ACCEL_TARGET STREQUAL "nvidia70")
                set(Kokkos_ARCH_VOLTA70 ON CACHE BOOL "Kokkos arch from CRAY_ACCEL_TARGET")
                message(STATUS "    Set Kokkos_ARCH_VOLTA70 = ON (from CRAY_ACCEL_TARGET)")
                set(KOKKOS_ARCH_SET TRUE)
                erf_cray_verbose("Mapped CRAY_ACCEL_TARGET=nvidia70 -> Kokkos_ARCH_VOLTA70=ON")
                
            elseif(CRAY_ACCEL_TARGET STREQUAL "nvidia80")
                set(Kokkos_ARCH_AMPERE80 ON CACHE BOOL "Kokkos arch from CRAY_ACCEL_TARGET")
                message(STATUS "    Set Kokkos_ARCH_AMPERE80 = ON (from CRAY_ACCEL_TARGET)")
                set(KOKKOS_ARCH_SET TRUE)
                erf_cray_verbose("Mapped CRAY_ACCEL_TARGET=nvidia80 -> Kokkos_ARCH_AMPERE80=ON")
                
            elseif(CRAY_ACCEL_TARGET STREQUAL "nvidia90")
                set(Kokkos_ARCH_HOPPER90 ON CACHE BOOL "Kokkos arch from CRAY_ACCEL_TARGET")
                message(STATUS "    Set Kokkos_ARCH_HOPPER90 = ON (from CRAY_ACCEL_TARGET)")
                set(KOKKOS_ARCH_SET TRUE)
                erf_cray_verbose("Mapped CRAY_ACCEL_TARGET=nvidia90 -> Kokkos_ARCH_HOPPER90=ON")
            endif()
        endif()
        
        if(NOT KOKKOS_ARCH_SET)
            message(WARNING "")
            message(WARNING "ERF: Kokkos architecture not detected")
            message(WARNING "  For Perlmutter: module load gpu")
            message(WARNING "  Or set: export KOKKOS_GPU_ARCH=AMPERE80")
            message(WARNING "  Or set: -DKokkos_ARCH_AMPERE80=ON")
            message(WARNING "")
        else()
            message(STATUS "")
            message(STATUS "    Note: After Kokkos configures, CMAKE_CUDA_ARCHITECTURES")
            message(STATUS "          will be set from Kokkos_CUDA_ARCHITECTURES")
            erf_cray_verbose("Kokkos will set CMAKE_CUDA_ARCHITECTURES when CUDA language is enabled")
        endif()
    endif()
    
    message(STATUS "")
endif()

# ==============================================================================
# Prerequisite Checks
# ==============================================================================

message(STATUS "ERF: Checking Cray prerequisites...")

# -----------------------------------------------------------------------------
# CMake Version Check
# -----------------------------------------------------------------------------
# Cray systems work best with CMake 3.24.0+
# Earlier versions may have issues with Cray wrappers and CUDA when NVHPC is splayed

set(ERF_RECOMMENDED_CMAKE_VERSION "3.24.0")

if(CMAKE_VERSION VERSION_LESS ${ERF_RECOMMENDED_CMAKE_VERSION})
    message(WARNING 
        "\n"
        "ERF: CMake version ${CMAKE_VERSION} detected\n"
        "  Recommended minimum for Cray systems: ${ERF_RECOMMENDED_CMAKE_VERSION}\n"
        "  You may experience issues with Cray compiler wrappers and CUDA\n"
        "\n"
        "  To fix:\n"
        "    module load cmake\n"
        "")
    
    erf_cray_verbose("Current CMake: ${CMAKE_VERSION}")
    erf_cray_verbose("Recommended: ${ERF_RECOMMENDED_CMAKE_VERSION}+")
    erf_cray_verbose("Known issues with older CMake on Cray:")
    erf_cray_verbose("  - CUDA language detection failures")
    erf_cray_verbose("  - Incorrect compiler wrapper handling")
    erf_cray_verbose("  - Missing Cray-specific find modules")
else()
    message(STATUS "  CMake version ${CMAKE_VERSION} >= ${ERF_RECOMMENDED_CMAKE_VERSION}")
    erf_cray_verbose("CMake version check passed")
endif()

# -----------------------------------------------------------------------------
# CUDA Toolkit Check
# -----------------------------------------------------------------------------
# When building with CUDA, the cudatoolkit module should be loaded
# This sets CUDA_HOME and other necessary environment variables

if(ERF_ENABLE_CUDA)
    message(STATUS "  Checking for CUDA toolkit...")
    
    set(CUDA_TOOLKIT_LOADED FALSE)
    
    # Check for CUDA_HOME (set by cudatoolkit module)
    if(DEFINED ENV{CUDA_HOME})
        message(STATUS "    CUDA_HOME = $ENV{CUDA_HOME}")
        set(CUDA_TOOLKIT_LOADED TRUE)
        erf_cray_verbose("CUDA toolkit appears to be loaded (CUDA_HOME set)")
    endif()
    
    # Additional check for CUDATOOLKIT_HOME (alternative Cray variable)
    if(DEFINED ENV{CUDATOOLKIT_HOME})
        message(STATUS "    CUDATOOLKIT_HOME = $ENV{CUDATOOLKIT_HOME}")
        set(CUDA_TOOLKIT_LOADED TRUE)
        erf_cray_verbose("CUDA toolkit appears to be loaded (CUDATOOLKIT_HOME set)")
    endif()
    
    # Check for nvcc in PATH
    find_program(NVCC_EXECUTABLE nvcc)
    if(NVCC_EXECUTABLE)
        message(STATUS "    Found nvcc: ${NVCC_EXECUTABLE}")
        set(CUDA_TOOLKIT_LOADED TRUE)
        erf_cray_verbose("nvcc found in PATH")
    endif()
    
    # Warn if CUDA toolkit doesn't appear to be loaded
    if(NOT CUDA_TOOLKIT_LOADED)
        message(WARNING
            "\n"
            "ERF: CUDA enabled but CUDA toolkit not detected\n"
            "  Expected environment variables not found:\n"
            "    - CUDA_HOME\n"
            "    - CUDATOOLKIT_HOME\n"
            "    - nvcc in PATH\n"
            "\n"
            "  To fix:\n"
            "    module load cudatoolkit\n"
            "  Or on newer systems:\n"
            "    module load cuda\n"
            "\n"
            "  Build may fail with CUDA-related errors\n"
            "")
        
        erf_cray_verbose("CUDA_HOME = $ENV{CUDA_HOME}")
        erf_cray_verbose("CUDATOOLKIT_HOME = $ENV{CUDATOOLKIT_HOME}")
        erf_cray_verbose("nvcc search result: ${NVCC_EXECUTABLE}")
    endif()
    
    # Check CUDA architecture is set for GPU builds
    if(NOT AMReX_CUDA_ARCH AND NOT DEFINED ENV{AMREX_CUDA_ARCH})
        message(WARNING "")
        message(WARNING "ERF: CUDA enabled but GPU architecture not specified")
        message(WARNING "  Set AMReX_CUDA_ARCH for optimal performance")
        message(WARNING "  For Perlmutter A100 GPUs:")
        message(WARNING "    -DAMReX_CUDA_ARCH=8.0")
        message(WARNING "  Or set in environment:")
        message(WARNING "    export AMREX_CUDA_ARCH=8.0")
        message(WARNING "")
        
        erf_cray_verbose("AMReX_CUDA_ARCH not set (will use CMake default)")
    else()
        if(AMReX_CUDA_ARCH)
            message(STATUS "    AMReX_CUDA_ARCH = ${AMReX_CUDA_ARCH}")
        else()
            message(STATUS "    AMREX_CUDA_ARCH = $ENV{AMREX_CUDA_ARCH}")
        endif()
    endif()
else()
    erf_cray_verbose("CUDA not enabled, skipping CUDA toolkit checks")
endif()

# -----------------------------------------------------------------------------
# NetCDF Module Check
# -----------------------------------------------------------------------------

if(ERF_ENABLE_NETCDF)
    message(STATUS "  Checking for NetCDF...")
    
    set(NETCDF_LOADED FALSE)
    
    if(DEFINED ENV{NETCDF_DIR})
        message(STATUS "    NETCDF_DIR = $ENV{NETCDF_DIR}")
        set(NETCDF_LOADED TRUE)
    endif()
    
    if(NOT NETCDF_LOADED)
        message(WARNING "")
        message(WARNING "ERF: NetCDF enabled but NETCDF_DIR not set")
        message(WARNING "  To fix:")
        message(WARNING "    module load cray-netcdf-hdf5parallel")
        message(WARNING "  Or:")
        message(WARNING "    module load cray-netcdf")
        message(WARNING "")
        
        erf_cray_verbose("NETCDF_DIR not found in environment")
    endif()
else()
    erf_cray_verbose("NetCDF not enabled, skipping NetCDF checks")
endif()

# -----------------------------------------------------------------------------
# Module Environment Summary
# -----------------------------------------------------------------------------

if(ERF_VERBOSE_CRAY_FIXES)
    message(STATUS "")
    message(STATUS "[VERBOSE] Key environment variables:")
    message(STATUS "[VERBOSE]   CRAYPE_VERSION      = $ENV{CRAYPE_VERSION}")
    message(STATUS "[VERBOSE]   CRAY_MPICH_DIR      = $ENV{CRAY_MPICH_DIR}")
    message(STATUS "[VERBOSE]   MPICH_DIR           = $ENV{MPICH_DIR}")
    message(STATUS "[VERBOSE]   CUDA_HOME           = $ENV{CUDA_HOME}")
    message(STATUS "[VERBOSE]   CUDATOOLKIT_HOME    = $ENV{CUDATOOLKIT_HOME}")
    message(STATUS "[VERBOSE]   NETCDF_DIR          = $ENV{NETCDF_DIR}")
    message(STATUS "[VERBOSE]   HDF5_DIR            = $ENV{HDF5_DIR}")
    message(STATUS "[VERBOSE]   MPICH_GPU_SUPPORT   = $ENV{MPICH_GPU_SUPPORT_ENABLED}")
    message(STATUS "")
endif()

message(STATUS "")

# ==============================================================================
# Fix 1: CUDA + EKAT -> nvcc_wrapper complications (Checklist Item 1)
# ==============================================================================
# PROBLEM: When building with EKAT, we get nvcc_wrapper which can cause 
#          "mpi.h not found" errors because nvcc_wrapper doesn't know about
#          Cray's include paths
# SOLUTION: Add Cray compiler flags to CUDA compilation via --cray-print-opts

if(ERF_ENABLE_CUDA AND (ERF_ENABLE_RRTMGP OR ERF_ENABLE_SHOC OR ERF_ENABLE_P3))
    message(STATUS "ERF: [Fix 1] Applying CUDA+EKAT nvcc_wrapper fix")
    
    erf_cray_verbose("Problem: EKAT uses nvcc_wrapper which doesn't inherit Cray paths")
    erf_cray_verbose("Condition: ERF_ENABLE_CUDA=ON and (RRTMGP or SHOC or P3 enabled)")
    erf_cray_verbose("Solution: Add Cray-specific flags from 'CC --cray-print-opts=cflags'")
    
    # Get Cray-specific flags
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} --cray-print-opts=cflags
        OUTPUT_VARIABLE CRAY_CUDA_FLAGS
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE CRAY_CUDA_FLAGS_RESULT
    )
    
    if(CRAY_CUDA_FLAGS_RESULT EQUAL 0 AND CRAY_CUDA_FLAGS)
        message(STATUS "  Adding Cray flags to CUDA compilation")
        erf_cray_verbose("Retrieved flags: ${CRAY_CUDA_FLAGS}")
        erf_cray_verbose("Command used: ${CMAKE_CXX_COMPILER} --cray-print-opts=cflags")
        
        if(CMAKE_CUDA_FLAGS)
            erf_cray_verbose("Appending to existing CMAKE_CUDA_FLAGS: ${CMAKE_CUDA_FLAGS}")
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} ${CRAY_CUDA_FLAGS}" CACHE STRING "" FORCE)
        else()
            erf_cray_verbose("Setting new CMAKE_CUDA_FLAGS")
            set(CMAKE_CUDA_FLAGS "${CRAY_CUDA_FLAGS}" CACHE STRING "" FORCE)
        endif()
        
        erf_cray_verbose("Final CMAKE_CUDA_FLAGS: ${CMAKE_CUDA_FLAGS}")
    else()
        message(WARNING "ERF: Could not retrieve Cray CUDA flags")
        message(WARNING "  Command attempted: ${CMAKE_CXX_COMPILER} --cray-print-opts=cflags")
        message(WARNING "  Return code: ${CRAY_CUDA_FLAGS_RESULT}")
        message(WARNING "  You may need to set CMAKE_CUDA_FLAGS manually")
        message(WARNING "  Example: -DCMAKE_CUDA_FLAGS=\"\$(CC --cray-print-opts=cflags)\"")
    endif()
else()
    erf_cray_verbose("Fix 1 not needed (CUDA+EKAT not both enabled)")
endif()

# ==============================================================================
# Fix 2: FCOMPARE + Cray -> mpi_gnu_123 not found (Checklist Item 2)
# ==============================================================================
# PROBLEM: When building with fcompare, Cray's --as-needed linker flag causes
#          the linker to drop MPI libraries it thinks aren't needed, leading to
#          "cannot find -lmpi_gnu_123" errors
# SOLUTION: Remove --as-needed from Cray library flags and add --no-as-needed

if(ERF_ENABLE_FCOMPARE)
    message(STATUS "ERF: [Fix 2] Applying fcompare linker fix")
    
    erf_cray_verbose("Problem: Cray uses --as-needed which drops required MPI libs")
    erf_cray_verbose("Condition: ERF_ENABLE_FCOMPARE=ON")
    erf_cray_verbose("Solution: Clean Cray lib flags and add --no-as-needed")
    
    # Get Cray library paths and clean them
    set(CRAY_LIBS_CLEAN "")
    set(COMPILERS_CHECKED "")
    
    foreach(COMPILER IN ITEMS ${CMAKE_CXX_COMPILER} ${CMAKE_C_COMPILER} ${CMAKE_Fortran_COMPILER})
        if(EXISTS ${COMPILER})
            erf_cray_verbose("Checking compiler: ${COMPILER}")
            
            execute_process(
                COMMAND ${COMPILER} --cray-print-opts=libs
                OUTPUT_VARIABLE COMPILER_LIBS
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
                RESULT_VARIABLE COMPILER_LIBS_RESULT
            )
            
            if(COMPILER_LIBS_RESULT EQUAL 0)
                erf_cray_verbose("  Original libs: ${COMPILER_LIBS}")
                
                # Remove problematic --as-needed flags
                string(REGEX REPLACE "-Wl,--as-needed," "" COMPILER_LIBS "${COMPILER_LIBS}")
                string(REGEX REPLACE ",--no-as-needed" "" COMPILER_LIBS "${COMPILER_LIBS}")
                string(REGEX REPLACE ",-l" " -l" COMPILER_LIBS "${COMPILER_LIBS}")
                
                erf_cray_verbose("  Cleaned libs:  ${COMPILER_LIBS}")
                
                set(CRAY_LIBS_CLEAN "${CRAY_LIBS_CLEAN} ${COMPILER_LIBS}")
                list(APPEND COMPILERS_CHECKED ${COMPILER})
            else()
                erf_cray_verbose("  Failed to get libs from ${COMPILER}")
            endif()
        endif()
    endforeach()
    
    if(CRAY_LIBS_CLEAN)
        message(STATUS "  Adding Cray linker flags: -Wl,--no-as-needed + libs")
        erf_cray_verbose("Compilers checked: ${COMPILERS_CHECKED}")
        erf_cray_verbose("Combined cleaned libs: ${CRAY_LIBS_CLEAN}")
        erf_cray_verbose("Final linker flags: -Wl,--no-as-needed ${CRAY_LIBS_CLEAN}")
        
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed ${CRAY_LIBS_CLEAN}" 
            CACHE STRING "" FORCE)
            
        erf_cray_verbose("CMAKE_EXE_LINKER_FLAGS updated")
    else()
        message(WARNING "ERF: Could not retrieve Cray library paths")
        message(WARNING "  Fcompare may fail to link with: cannot find -lmpi_gnu_123")
        message(WARNING "  Workaround: Set CMAKE_EXE_LINKER_FLAGS manually")
        message(WARNING "  Example: -DCMAKE_EXE_LINKER_FLAGS=\"-Wl,--no-as-needed \$CRAY_LIBS_CLEAN\"")
        erf_cray_verbose("No compilers returned valid library flags")
    endif()
else()
    erf_cray_verbose("Fix 2 not needed (ERF_ENABLE_FCOMPARE=OFF)")
endif()

# ==============================================================================
# Fix 3: CUDA without cmake module -> math libs not found (Checklist Item 3)
# ==============================================================================
# PROBLEM: If 'module load cmake' isn't run, CMAKE_PREFIX_PATH may not include
#          CUDA math libraries path, causing link errors for cuBLAS, cuRAND, etc.
# SOLUTION: Add $CUDA_HOME/../../math_libs/lib64 to CMAKE_PREFIX_PATH

if(ERF_ENABLE_CUDA AND DEFINED ENV{CUDA_HOME})
    set(CUDA_MATH_PATH "$ENV{CUDA_HOME}/../../math_libs/lib64")
    
    erf_cray_verbose("Checking for CUDA math libraries...")
    erf_cray_verbose("CUDA_HOME = $ENV{CUDA_HOME}")
    erf_cray_verbose("Expected math libs path: ${CUDA_MATH_PATH}")
    
    if(EXISTS ${CUDA_MATH_PATH})
        message(STATUS "ERF: [Fix 3] Adding CUDA math libraries path")
        message(STATUS "  ${CUDA_MATH_PATH}")
        
        erf_cray_verbose("Problem: CUDA math libs may not be in default search path")
        erf_cray_verbose("Condition: ERF_ENABLE_CUDA=ON and CUDA_HOME set")
        erf_cray_verbose("Solution: Add CUDA_HOME/../../math_libs/lib64 to CMAKE_PREFIX_PATH")
        erf_cray_verbose("Path exists: YES")
        
        list(APPEND CMAKE_PREFIX_PATH ${CUDA_MATH_PATH})
        set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} CACHE STRING "" FORCE)
        
        erf_cray_verbose("CMAKE_PREFIX_PATH updated: ${CMAKE_PREFIX_PATH}")
    else()
        message(WARNING "ERF: CUDA math libs path not found at ${CUDA_MATH_PATH}")
        message(WARNING "  You may need to 'module load cuda' or set CMAKE_PREFIX_PATH manually")
        message(WARNING "  Expected libraries: cuBLAS, cuRAND, cuSPARSE, etc.")
        erf_cray_verbose("Path exists: NO")
        erf_cray_verbose("This may cause link errors for CUDA math libraries")
    endif()
else()
    if(ERF_ENABLE_CUDA AND NOT DEFINED ENV{CUDA_HOME})
        message(WARNING "ERF: CUDA enabled but CUDA_HOME not set")
        message(WARNING "  Math libraries may not be found")
        message(WARNING "  Solution: Load CUDA module or set CUDA_HOME")
        erf_cray_verbose("CUDA_HOME not set in environment")
    else()
        erf_cray_verbose("Fix 3 not needed (ERF_ENABLE_CUDA=OFF)")
    endif()
endif()

# ==============================================================================
# Fix 4: GPU-aware MPI with Cray GTL (Checklist Item 4)
# ==============================================================================
# PROBLEM: GPU-aware MPI on Cray requires linking against mpi_gtl_cuda library
#          which enables GPU Transfer Library for direct GPU-GPU communication
# SOLUTION: Detect GPU-aware MPI and add GTL libraries to link flags

# ==============================================================================
# Fix 4: GPU-aware MPI with Cray GTL (Checklist Item 4)
# ==============================================================================
# PROBLEM: GPU-aware MPI on Cray requires linking against mpi_gtl_cuda library
#          which enables GPU Transfer Library for direct GPU-GPU communication
# SOLUTION: Detect GPU-aware MPI and add GTL libraries to link flags

if(ERF_ENABLE_CUDA AND ERF_ENABLE_MPI AND "$ENV{MPICH_GPU_SUPPORT_ENABLED}" STREQUAL "1")
    message(STATUS "ERF: [Fix 4] Applying GPU-aware MPI fix (Cray GTL)")
    
    erf_cray_verbose("Problem: GPU-aware MPI needs Cray GTL libraries")
    erf_cray_verbose("Condition: CUDA + MPI + MPICH_GPU_SUPPORT_ENABLED=1")
    erf_cray_verbose("Solution: Add -lmpi_gnu_123 -lmpi_gtl_cuda to link flags")
    erf_cray_verbose("MPICH_GPU_SUPPORT_ENABLED = $ENV{MPICH_GPU_SUPPORT_ENABLED}")
    
    # Set the MPI+GTL libraries
    # Note: We use -lmpi_gnu_123 explicitly because Cray's --as-needed can drop it
    set(CRAY_MPI_LIBS "-lmpi_gnu_123 -lmpi_gtl_cuda")
    
    # Try to verify the library exists (for diagnostics)
    set(MPI_LIB_SEARCH_PATHS "")
    if(DEFINED ENV{MPICH_DIR})
        list(APPEND MPI_LIB_SEARCH_PATHS "$ENV{MPICH_DIR}/lib")
    endif()
    if(DEFINED ENV{CRAY_MPICH_DIR})
        list(APPEND MPI_LIB_SEARCH_PATHS "$ENV{CRAY_MPICH_DIR}/lib")
    endif()
    
    erf_cray_verbose("Searching for mpi_gtl_cuda library in:")
    foreach(path IN LISTS MPI_LIB_SEARCH_PATHS)
        erf_cray_verbose("  ${path}")
    endforeach()
    
    find_library(CRAY_MPI_GTL_CUDA 
        NAMES mpi_gtl_cuda
        HINTS ${MPI_LIB_SEARCH_PATHS}
        NO_DEFAULT_PATH
    )
    
    if(CRAY_MPI_GTL_CUDA)
        message(STATUS "  Found GTL library: ${CRAY_MPI_GTL_CUDA}")
        erf_cray_verbose("Library verification successful")
    else()
        message(STATUS "  GTL library not found via find_library (will rely on linker search)")
        erf_cray_verbose("Library not found in search paths, but linker may still find it")
        erf_cray_verbose("This is normal if libraries are in non-standard Cray locations")
    endif()
    
    # Apply the fix regardless of whether find_library succeeded
    # The Cray linker knows where to find these libraries
    message(STATUS "  Adding MPI+GTL libraries: ${CRAY_MPI_LIBS}")
    erf_cray_verbose("Adding to CMAKE_CUDA_STANDARD_LIBRARIES")
    erf_cray_verbose("Adding to CMAKE_CXX_STANDARD_LIBRARIES")
    
    set(CMAKE_CUDA_STANDARD_LIBRARIES "${CMAKE_CUDA_STANDARD_LIBRARIES} ${CRAY_MPI_LIBS}" 
        CACHE STRING "" FORCE)
    set(CMAKE_CXX_STANDARD_LIBRARIES "${CMAKE_CXX_STANDARD_LIBRARIES} ${CRAY_MPI_LIBS}" 
        CACHE STRING "" FORCE)
        
    erf_cray_verbose("CMAKE_CUDA_STANDARD_LIBRARIES: ${CMAKE_CUDA_STANDARD_LIBRARIES}")
    erf_cray_verbose("CMAKE_CXX_STANDARD_LIBRARIES: ${CMAKE_CXX_STANDARD_LIBRARIES}")
    
else()
    if(ERF_ENABLE_CUDA AND ERF_ENABLE_MPI)
        if(NOT DEFINED ENV{MPICH_GPU_SUPPORT_ENABLED})
            message(STATUS "")
            message(STATUS "  Note: MPICH_GPU_SUPPORT_ENABLED not set")
            message(STATUS "    For GPU-aware MPI, add to your script:")
            message(STATUS "      export MPICH_GPU_SUPPORT_ENABLED=1")
            message(STATUS "")
            erf_cray_verbose("Fix 4 not applied: MPICH_GPU_SUPPORT_ENABLED not set")
        elseif(NOT "$ENV{MPICH_GPU_SUPPORT_ENABLED}" STREQUAL "1")
            erf_cray_verbose("Fix 4 not applied: MPICH_GPU_SUPPORT_ENABLED=$ENV{MPICH_GPU_SUPPORT_ENABLED} (not '1')")
        endif()
    else()
        erf_cray_verbose("Fix 4 not needed (CUDA+MPI not both enabled)")
    endif()
endif()

# ==============================================================================
# Fix 5-6: NetCDF with cray-netcdf-hdf5parallel (Checklist Items 5-6)
# ==============================================================================
# PROBLEM 5: Cray NetCDF may use different C++ library names or structures
# PROBLEM 6: pkg-config may not find MPI/NetCDF without correct PKG_CONFIG_PATH
# SOLUTION: Set up pkg-config path and add NetCDF/HDF5 directories to search

if(ERF_ENABLE_NETCDF)
    message(STATUS "ERF: [Fix 5-6] Configuring NetCDF with Cray paths")
    
    # Get PKG_CONFIG_PATH directly from Cray wrapper
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} --cray-print-opts=PKG_CONFIG_PATH
        OUTPUT_VARIABLE CRAY_PKG_CONFIG_PATH
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE PKG_RESULT
    )
    
    if(PKG_RESULT EQUAL 0 AND CRAY_PKG_CONFIG_PATH)
        message(STATUS "  Setting PKG_CONFIG_PATH from Cray wrapper")
        
        # Append to existing PKG_CONFIG_PATH
        if(DEFINED ENV{PKG_CONFIG_PATH})
            set(ENV{PKG_CONFIG_PATH} "${CRAY_PKG_CONFIG_PATH}:$ENV{PKG_CONFIG_PATH}")
        else()
            set(ENV{PKG_CONFIG_PATH} "${CRAY_PKG_CONFIG_PATH}")
        endif()
        
        message(STATUS "  PKG_CONFIG_PATH = $ENV{PKG_CONFIG_PATH}")
    endif()
    
    message(STATUS "  PKG_CONFIG_PATH = ${PKG_CONFIG_PATH}")
    erf_cray_verbose("This allows cmake/gnumake to find MPI and NetCDF via pkg-config")
    
    # Help find NetCDF (may be named differently on Cray)
    if(DEFINED ENV{NETCDF_DIR})
        list(APPEND CMAKE_PREFIX_PATH $ENV{NETCDF_DIR})
        message(STATUS "  Added NETCDF_DIR to search path: $ENV{NETCDF_DIR}")
        erf_cray_verbose("NetCDF headers/libs will be searched in NETCDF_DIR")
    else()
        erf_cray_verbose("NETCDF_DIR not set (may still work via module)")
    endif()
    
    if(DEFINED ENV{HDF5_DIR})
        list(APPEND CMAKE_PREFIX_PATH $ENV{HDF5_DIR})
        message(STATUS "  Added HDF5_DIR to search path: $ENV{HDF5_DIR}")
        erf_cray_verbose("HDF5 headers/libs will be searched in HDF5_DIR")
    else()
        erf_cray_verbose("HDF5_DIR not set (may still work via module)")
    endif()
    
    erf_cray_verbose("CMAKE_PREFIX_PATH now includes: ${CMAKE_PREFIX_PATH}")
else()
    if(ERF_ENABLE_NETCDF)
        erf_cray_verbose("Fix 5-6 not fully applied: MPICH_DIR not set")
        message(WARNING "ERF: NetCDF enabled but MPICH_DIR not set")
        message(WARNING "  pkg-config may not find MPI libraries")
        message(WARNING "  Load MPI module or set MPICH_DIR")
    else()
        erf_cray_verbose("Fix 5-6 not needed (ERF_ENABLE_NETCDF=OFF)")
    endif()
endif()

# ==============================================================================
# Summary
# ==============================================================================

message(STATUS "")
message(STATUS "ERF: Cray system fixes summary")
message(STATUS "══════════════════════════════════════════════════════════════")

# Fix 1: CUDA + EKAT
set(FIX1_ACTIVE OFF)
if(ERF_ENABLE_CUDA AND (ERF_ENABLE_RRTMGP OR ERF_ENABLE_SHOC OR ERF_ENABLE_P3))
    set(FIX1_ACTIVE ON)
endif()
message(STATUS "  Fix 1 (CUDA+EKAT):     ${FIX1_ACTIVE}")
if(FIX1_ACTIVE AND CRAY_CUDA_FLAGS)
    message(STATUS "    Applied Cray CUDA flags:")
    message(STATUS "      ${CRAY_CUDA_FLAGS}")
    message(STATUS "")
    message(STATUS "    Command line equivalent:")
    message(STATUS "      -DCMAKE_CUDA_FLAGS=\"\$(CC --cray-print-opts=cflags)\"")
endif()

# Fix 2: fcompare
message(STATUS "")
message(STATUS "  Fix 2 (fcompare):      ${ERF_ENABLE_FCOMPARE}")
if(ERF_ENABLE_FCOMPARE AND CRAY_LIBS_CLEAN)
    message(STATUS "    Applied Cray library cleanup:")
    message(STATUS "      ${CRAY_LIBS_CLEAN}")
    message(STATUS "")
    message(STATUS "    Command line equivalent:")
    message(STATUS "      CRAY_LIBS=\"\$(CC --cray-print-opts=libs | sed 's/-Wl,--as-needed,//g; s/,--no-as-needed//g; s/,-l/ -l/g')\"")
    message(STATUS "      CRAY_LIBS=\"\$CRAY_LIBS \$(cc --cray-print-opts=libs | sed ...)\"")
    message(STATUS "      CRAY_LIBS=\"\$CRAY_LIBS \$(ftn --cray-print-opts=libs | sed ...)\"")
    message(STATUS "      -DCMAKE_EXE_LINKER_FLAGS=\"-Wl,--no-as-needed \$CRAY_LIBS\"")
    message(STATUS "")
    message(STATUS "    What was actually set:")
    message(STATUS "      CMAKE_EXE_LINKER_FLAGS=\"-Wl,--no-as-needed ${CRAY_LIBS_CLEAN}\"")
endif()

# Fix 3: CUDA math libs
set(FIX3_ACTIVE OFF)
if(ERF_ENABLE_CUDA AND DEFINED ENV{CUDA_HOME})
    set(CUDA_MATH_PATH_CHECK "$ENV{CUDA_HOME}/../../math_libs/lib64")
    if(EXISTS ${CUDA_MATH_PATH_CHECK})
        set(FIX3_ACTIVE ON)
    endif()
endif()
message(STATUS "")
message(STATUS "  Fix 3 (CUDA math):     ${FIX3_ACTIVE}")
if(FIX3_ACTIVE)
    message(STATUS "    Command line equivalent:")
    message(STATUS "      -DCMAKE_PREFIX_PATH=\"\$CUDA_HOME/../../math_libs/lib64\"")
endif()

# Fix 4: GPU-aware MPI
set(FIX4_ACTIVE OFF)
if(ERF_ENABLE_CUDA AND ERF_ENABLE_MPI AND DEFINED ENV{MPICH_GPU_SUPPORT_ENABLED})
    if(CRAY_MPI_GTL_CUDA)
        set(FIX4_ACTIVE ON)
    endif()
endif()
message(STATUS "")
message(STATUS "  Fix 4 (GPU-aware MPI): ${FIX4_ACTIVE}")
if(FIX4_ACTIVE)
    message(STATUS "    Command line equivalent:")
    message(STATUS "      export MPICH_GPU_SUPPORT_ENABLED=1")
    message(STATUS "      -DCMAKE_CUDA_STANDARD_LIBRARIES=\"-lmpi_gnu_123 -lmpi_gtl_cuda\"")
    message(STATUS "      -DCMAKE_CXX_STANDARD_LIBRARIES=\"-lmpi_gnu_123 -lmpi_gtl_cuda\"")
endif()

# Fix 5-6: NetCDF
set(FIX56_ACTIVE OFF)
if(ERF_ENABLE_NETCDF AND DEFINED ENV{MPICH_DIR})
    set(FIX56_ACTIVE ON)
endif()
message(STATUS "")
message(STATUS "  Fix 5-6 (NetCDF):      ${FIX56_ACTIVE}")
if(FIX56_ACTIVE)
    message(STATUS "    Command line equivalent:")
    message(STATUS "      export PKG_CONFIG_PATH=\"\$MPICH_DIR/lib/pkgconfig:\$PKG_CONFIG_PATH\"")
    if(DEFINED ENV{NETCDF_DIR})
        message(STATUS "      -DCMAKE_PREFIX_PATH=\"\$NETCDF_DIR\"")
    endif()
    if(DEFINED ENV{HDF5_DIR})
        message(STATUS "      -DCMAKE_PREFIX_PATH=\"\$CMAKE_PREFIX_PATH:\$HDF5_DIR\"")
    endif()
endif()

message(STATUS "")
message(STATUS "══════════════════════════════════════════════════════════════")
message(STATUS "  To disable auto-fixes: -DERF_DISABLE_CRAY_AUTO_FIXES=ON")
message(STATUS "  To see verbose output: -DERF_VERBOSE_CRAY_FIXES=ON")
message(STATUS "  To override any fix:   Set the corresponding CMAKE_* variable explicitly")
message(STATUS "")
message(STATUS "  Complete manual equivalent (all active fixes):")
message(STATUS "  ------------------------------------------------")
if(FIX1_ACTIVE)
message(STATUS "  -DCMAKE_CUDA_FLAGS=\"\$(CC --cray-print-opts=cflags)\" \\")
endif()
if(ERF_ENABLE_FCOMPARE AND CRAY_LIBS_CLEAN)
message(STATUS "  -DCMAKE_EXE_LINKER_FLAGS=\"-Wl,--no-as-needed ${CRAY_LIBS_CLEAN}\" \\")
endif()
if(FIX3_ACTIVE)
message(STATUS "  -DCMAKE_PREFIX_PATH=\"\$CUDA_HOME/../../math_libs/lib64\" \\")
endif()
if(FIX4_ACTIVE)
message(STATUS "  -DCMAKE_CUDA_STANDARD_LIBRARIES=\"-lmpi_gnu_123 -lmpi_gtl_cuda\" \\")
message(STATUS "  -DCMAKE_CXX_STANDARD_LIBRARIES=\"-lmpi_gnu_123 -lmpi_gtl_cuda\" \\")
endif()
message(STATUS "")

if(ERF_VERBOSE_CRAY_FIXES)
    message(STATUS "[VERBOSE] All Cray fixes processing complete")
    message(STATUS "[VERBOSE] Review messages above for detailed information")
    message(STATUS "[VERBOSE] The command-line equivalents above show what this module does automatically")
endif()