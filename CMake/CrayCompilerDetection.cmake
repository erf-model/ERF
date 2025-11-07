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