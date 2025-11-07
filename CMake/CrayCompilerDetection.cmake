# ==============================================================================
# Cray Compiler Detection (Pre-Project Stage)
# ==============================================================================
# This file runs BEFORE project() to detect Cray systems and set compilers
# The main CrayDetection.cmake runs AFTER project() to apply build fixes
# ==============================================================================

# Skip if user already set compilers explicitly
if(CMAKE_C_COMPILER OR CMAKE_CXX_COMPILER OR CMAKE_Fortran_COMPILER)
    message(STATUS "ERF: Compilers already specified by user, skipping Cray auto-detection")
    return()
endif()

# -----------------------------------------------------------------------------
# Detect Cray Environment (using only environment variables)
# -----------------------------------------------------------------------------

set(ERF_ON_CRAY_PREPROJECT FALSE)

# Check for Cray Programming Environment
if(DEFINED ENV{CRAYPE_VERSION})
    set(ERF_ON_CRAY_PREPROJECT TRUE)
    message(STATUS "ERF: Detected Cray Programming Environment (CRAYPE_VERSION=$ENV{CRAYPE_VERSION})")
endif()

# Additional check for Cray MPI
if(DEFINED ENV{CRAY_MPICH_DIR})
    set(ERF_ON_CRAY_PREPROJECT TRUE)
    message(STATUS "ERF: Detected Cray MPI (CRAY_MPICH_DIR=$ENV{CRAY_MPICH_DIR})")
endif()

# Check for Cray compiler module
if(DEFINED ENV{PE_ENV})
    set(ERF_ON_CRAY_PREPROJECT TRUE)
    message(STATUS "ERF: Detected Cray Programming Environment: $ENV{PE_ENV}")
endif()

if(NOT ERF_ON_CRAY_PREPROJECT)
    # Not on Cray, skip compiler setup
    return()
endif()

# -----------------------------------------------------------------------------
# Set Cray Compiler Wrappers as Defaults
# -----------------------------------------------------------------------------

message(STATUS "ERF: Setting Cray compiler wrappers...")

# Find Cray C compiler wrapper
find_program(ERF_CRAY_CC cc)
if(ERF_CRAY_CC)
    set(CMAKE_C_COMPILER "${ERF_CRAY_CC}" CACHE FILEPATH "C compiler")
    message(STATUS "  Set CMAKE_C_COMPILER = ${ERF_CRAY_CC}")
else()
    message(WARNING "ERF: On Cray system but 'cc' wrapper not found in PATH")
endif()

# Find Cray C++ compiler wrapper
find_program(ERF_CRAY_CXX CC)
if(ERF_CRAY_CXX)
    set(CMAKE_CXX_COMPILER "${ERF_CRAY_CXX}" CACHE FILEPATH "C++ compiler")
    message(STATUS "  Set CMAKE_CXX_COMPILER = ${ERF_CRAY_CXX}")
else()
    message(WARNING "ERF: On Cray system but 'CC' wrapper not found in PATH")
endif()

# Find Cray Fortran compiler wrapper (if needed)
if(ERF_ENABLE_MORR_FORT OR ERF_ENABLE_NOAHMP)
    find_program(ERF_CRAY_FC ftn)
    if(ERF_CRAY_FC)
        set(CMAKE_Fortran_COMPILER "${ERF_CRAY_FC}" CACHE FILEPATH "Fortran compiler")
        message(STATUS "  Set CMAKE_Fortran_COMPILER = ${ERF_CRAY_FC}")
    else()
        message(WARNING "ERF: On Cray system but 'ftn' wrapper not found in PATH")
    endif()
endif()

message(STATUS "")

# -----------------------------------------------------------------------------
# GPU Host Compilers (for CUDA, HIP, SYCL)
# -----------------------------------------------------------------------------

# CUDA Host Compiler - detect via environment
if(DEFINED ENV{CUDA_HOME} OR DEFINED ENV{CUDATOOLKIT_HOME} OR DEFINED ENV{CRAY_ACCEL_TARGET})
    message(STATUS "  Detected CUDA environment, configuring CUDA host compiler...")
    
    # Only set if not already specified by user
    if(NOT CMAKE_CUDA_HOST_COMPILER AND NOT DEFINED ENV{CUDAHOSTCXX})
        if(ERF_CRAY_CXX)
            set(CMAKE_CUDA_HOST_COMPILER "${ERF_CRAY_CXX}" CACHE FILEPATH "CUDA host compiler" FORCE)
            message(STATUS "    Set CMAKE_CUDA_HOST_COMPILER = ${ERF_CRAY_CXX}")
        endif()
    else()
        message(STATUS "    CUDA host compiler already set by user")
    endif()
endif()

# HIP Host Compiler - detect via ROCM environment
if(DEFINED ENV{ROCM_PATH} OR DEFINED ENV{HIP_PATH})
    message(STATUS "  Detected ROCm/HIP environment, configuring HIP host compiler...")
    
    # Only set if not already specified by user
    if(NOT CMAKE_HIP_HOST_COMPILER AND NOT DEFINED ENV{HIPHOSTCXX})
        if(ERF_CRAY_CXX)
            set(CMAKE_HIP_HOST_COMPILER "${ERF_CRAY_CXX}" CACHE FILEPATH "HIP host compiler" FORCE)
            message(STATUS "    Set CMAKE_HIP_HOST_COMPILER = ${ERF_CRAY_CXX}")
        endif()
    else()
        message(STATUS "    HIP host compiler already set by user")
    endif()
endif()

# SYCL - detect via Intel oneAPI
if(DEFINED ENV{ONEAPI_ROOT} OR DEFINED ENV{I_MPI_ROOT})
    message(STATUS "  Detected Intel oneAPI environment")
    message(STATUS "    SYCL will use CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER}")
endif()