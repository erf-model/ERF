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

# -----------------------------------------------------------------------------
# Detect Cray MPI and GTL Library Names (with smart fallbacks)
# -----------------------------------------------------------------------------

set(CRAY_MPI_LIB "")
set(CRAY_GTL_LIB "")

# Method 1: Parse from CC --cray-print-opts=libs (BEST)
if(ERF_CRAY_CXX)
    execute_process(
        COMMAND ${ERF_CRAY_CXX} --cray-print-opts=libs
        OUTPUT_VARIABLE CRAY_LIBS_OUTPUT
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE CRAY_LIBS_RESULT
    )
    
    if(CRAY_LIBS_RESULT EQUAL 0 AND CRAY_LIBS_OUTPUT)
        string(REGEX MATCH "-lmpi_gnu_[0-9]+" CRAY_MPI_LIB "${CRAY_LIBS_OUTPUT}")
        string(REGEX MATCH "-lmpi_gtl_[a-z]+" CRAY_GTL_LIB "${CRAY_LIBS_OUTPUT}")
    endif()
endif()

# Method 2: Fallback from environment variables
if(NOT CRAY_MPI_LIB AND DEFINED ENV{CRAY_MPICH_DIR})
    string(REGEX MATCH "/gnu/([0-9]+)\\.([0-9]+)" MATCH_RESULT "$ENV{CRAY_MPICH_DIR}")
    if(CMAKE_MATCH_1 AND CMAKE_MATCH_2)
        set(CRAY_MPI_LIB "-lmpi_gnu_${CMAKE_MATCH_1}${CMAKE_MATCH_2}")
    endif()
endif()

if(NOT CRAY_GTL_LIB AND DEFINED ENV{CRAY_ACCEL_TARGET})
    set(GTL_VAR "PE_MPICH_GTL_LIBS_$ENV{CRAY_ACCEL_TARGET}")
    if(DEFINED ENV{${GTL_VAR}})
        set(CRAY_GTL_LIB "$ENV{${GTL_VAR}}")
    elseif("$ENV{CRAY_ACCEL_TARGET}" MATCHES "nvidia")
        set(CRAY_GTL_LIB "-lmpi_gtl_cuda")
    elseif("$ENV{CRAY_ACCEL_TARGET}" MATCHES "amd")
        set(CRAY_GTL_LIB "-lmpi_gtl_hsa")
    endif()
endif()

# Method 3: Ultimate fallback
if(NOT CRAY_MPI_LIB)
    set(CRAY_MPI_LIB "-lmpi")
endif()

# Combine
if(CRAY_MPI_LIB AND CRAY_GTL_LIB)
    set(GTL_LIBS "${CRAY_MPI_LIB} ${CRAY_GTL_LIB}")
elseif(CRAY_MPI_LIB)
    set(GTL_LIBS "${CRAY_MPI_LIB}")
endif()

# -----------------------------------------------------------------------------
# Set Minimal Flags for Compiler Tests (using detected libraries)
# -----------------------------------------------------------------------------
message(STATUS "ERF: Setting minimal flags for compiler tests...")

if(DEFINED ENV{MPICH_GPU_SUPPORT_ENABLED} AND "$ENV{MPICH_GPU_SUPPORT_ENABLED}" STREQUAL "1")
    message(STATUS "  GPU-aware MPI detected, adding CUDA runtime for tests")
    message(STATUS "  Detected libraries: ${GTL_LIBS}")

    # APPEND to linker flags
    if(CMAKE_EXE_LINKER_FLAGS)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lcudart -lcuda" CACHE STRING "" FORCE)
    else()
        set(CMAKE_EXE_LINKER_FLAGS "-lcudart -lcuda" CACHE STRING "" FORCE)
    endif()

    # APPEND to standard libraries (use DETECTED GTL_LIBS, not hardcoded!)
    if(CMAKE_CXX_STANDARD_LIBRARIES)
        set(CMAKE_CXX_STANDARD_LIBRARIES "${CMAKE_CXX_STANDARD_LIBRARIES} ${GTL_LIBS}" CACHE STRING "" FORCE)
    else()
        set(CMAKE_CXX_STANDARD_LIBRARIES "${GTL_LIBS}" CACHE STRING "" FORCE)
    endif()

    if(CMAKE_CUDA_STANDARD_LIBRARIES)
        set(CMAKE_CUDA_STANDARD_LIBRARIES "${CMAKE_CUDA_STANDARD_LIBRARIES} ${GTL_LIBS}" CACHE STRING "" FORCE)
    else()
        set(CMAKE_CUDA_STANDARD_LIBRARIES "${GTL_LIBS}" CACHE STRING "" FORCE)
    endif()
endif()