# - Find NetCDF
# Find the native NetCDF includes and library
#
#  NETCDF_INCLUDE_DIRS - where to find netcdf.h, etc
#  NETCDF_FOUND        - True if NetCDF found
#
# The following are not for general use and are included in
# NETCDF_LIBRARIES if the corresponding option above is set.
#
#  NETCDF_LIBRARIES      - only the libraries (without the '-l')
#  NETCDF_LINK_LIBRARIES - the libraries and their absolute paths
#  NETCDF_LDFLAGS        - all required linker flags
#
# Normal usage would be:
#  find_package (NetCDF REQUIRED)
#  target_link_libraries (target_name PUBLIC ${NETCDF_LINK_LIBRARIES})

# Set FindNetCDF context
list(APPEND CMAKE_MESSAGE_CONTEXT "FindNetCDF")

message(DEBUG "Starting NetCDF detection")

# Detection log for failures
set(NETCDF_DETECTION_LOG "")

# Check cache
if(NETCDF_INCLUDES AND NETCDF_LIBRARIES)
    set(NETCDF_FIND_QUIETLY TRUE)
    message(VERBOSE "NetCDF already in cache")
    message(DEBUG "  NETCDF_INCLUDES: ${NETCDF_INCLUDES}")
    message(DEBUG "  NETCDF_LIBRARIES: ${NETCDF_LIBRARIES}")
endif()

# Build hints
set(NETCDF_INCLUDE_HINTS)
set(NETCDF_LIBRARY_HINTS)

message(DEBUG "Building search hints")

if(NETCDF_DIR)
    list(APPEND NETCDF_INCLUDE_HINTS ${NETCDF_DIR}/include)
    list(APPEND NETCDF_LIBRARY_HINTS ${NETCDF_DIR}/lib)
    message(VERBOSE "Using NETCDF_DIR: ${NETCDF_DIR}")
    list(APPEND NETCDF_DETECTION_LOG "NETCDF_DIR=${NETCDF_DIR}")
else()
    message(DEBUG "NETCDF_DIR not set")
    list(APPEND NETCDF_DETECTION_LOG "NETCDF_DIR not set")
endif()

if(DEFINED ENV{NETCDF_DIR})
    list(APPEND NETCDF_INCLUDE_HINTS $ENV{NETCDF_DIR}/include)
    list(APPEND NETCDF_LIBRARY_HINTS $ENV{NETCDF_DIR}/lib)
    message(VERBOSE "Using ENV NETCDF_DIR: $ENV{NETCDF_DIR}")
    list(APPEND NETCDF_DETECTION_LOG "ENV NETCDF_DIR=$ENV{NETCDF_DIR}")
else()
    list(APPEND NETCDF_DETECTION_LOG "ENV NETCDF_DIR not set")
endif()

if(NETCDF_INCLUDE_DIR)
    list(APPEND NETCDF_INCLUDE_HINTS ${NETCDF_INCLUDE_DIR})
endif()

if(NETCDF_LIBRARY_DIR)
    list(APPEND NETCDF_LIBRARY_HINTS ${NETCDF_LIBRARY_DIR})
endif()

# Pkg-config
message(VERBOSE "Attempting pkg-config detection")
set(ENV{PKG_CONFIG_PATH} "$ENV{MPICH_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
message(DEBUG "PKG_CONFIG_PATH: $ENV{PKG_CONFIG_PATH}")

find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
    message(DEBUG "pkg-config available")
    
    set(PKG_VARIANTS netcdf netcdf-mpi netcdf_parallel netcdf-cxx4_parallel)
    foreach(variant ${PKG_VARIANTS})
        if(NOT NETCDF_FOUND)
            message(DEBUG "  Trying: ${variant}")
            pkg_check_modules(NETCDF QUIET IMPORTED_TARGET ${variant})
            
            if(NETCDF_FOUND)
                message(VERBOSE "Found via pkg-config: ${variant}")
                message(DEBUG "  Version: ${NETCDF_VERSION}")
                list(APPEND NETCDF_DETECTION_LOG "pkg-config ${variant}: found")
                list(APPEND NETCDF_INCLUDE_HINTS ${NETCDF_INCLUDE_DIRS})
                list(APPEND NETCDF_LIBRARY_HINTS ${NETCDF_LIBRARY_DIRS})
                break()
            else()
                list(APPEND NETCDF_DETECTION_LOG "pkg-config ${variant}: not found")
            endif()
        endif()
    endforeach()
else()
    message(DEBUG "pkg-config not available")
    list(APPEND NETCDF_DETECTION_LOG "pkg-config: not available")
endif()

# Manual search
message(VERBOSE "Searching for netcdf.h and libnetcdf")
message(DEBUG "  Include hints: ${NETCDF_INCLUDE_HINTS}")
message(DEBUG "  Library hints: ${NETCDF_LIBRARY_HINTS}")

find_path(NETCDF_INCLUDES netcdf.h
    HINTS ${NETCDF_INCLUDE_HINTS}
          $ENV{NETCDF_DIR}/include)

if(NETCDF_INCLUDES)
    message(VERBOSE "Found netcdf.h: ${NETCDF_INCLUDES}")
    list(APPEND NETCDF_DETECTION_LOG "find_path: ${NETCDF_INCLUDES}")
else()
    message(DEBUG "netcdf.h not found")
    list(APPEND NETCDF_DETECTION_LOG "find_path: failed")
endif()

find_library(NETCDF_LIBRARIES_C NAMES netcdf
    HINTS ${NETCDF_LIBRARY_HINTS}
          $ENV{NETCDF_DIR}/lib)
mark_as_advanced(NETCDF_LIBRARIES_C)

if(NETCDF_LIBRARIES_C)
    message(VERBOSE "Found libnetcdf: ${NETCDF_LIBRARIES_C}")
    list(APPEND NETCDF_DETECTION_LOG "find_library: ${NETCDF_LIBRARIES_C}")
    
    # Only add HDF5 if pkg-config told us NetCDF needs it
    if(NETCDF_LINK_LIBRARIES)
        # Check if pkg-config's library list includes hdf5
        string(FIND "${NETCDF_LINK_LIBRARIES}" "hdf5" HDF5_IN_NETCDF)
        if(HDF5_IN_NETCDF GREATER -1)
            message(STATUS "NetCDF was built with HDF5 support")
            # Check if HDF5 was already found (e.g., by AMReX)
            if(TARGET hdf5::hdf5 OR HDF5_FOUND)
                list(APPEND NETCDF_LIBRARIES_C ${HDF5_LIBRARIES})
                message(STATUS "  Using HDF5 libraries (already found): ${HDF5_LIBRARIES}")
            else()
                # Fallback: use pkg-config's complete library list which includes HDF5
                set(NETCDF_LIBRARIES_C ${NETCDF_LINK_LIBRARIES})
                message(STATUS "  HDF5 not already a target, using pkg-config's complete library list:")
                message(STATUS "  NETCDF_LIBRARIES_C = ${NETCDF_LINK_LIBRARIES}")
            endif()
        else()
            message(STATUS "NetCDF has no HDF5 dependency in pkg-config")
        endif()
    endif()  # <-- THIS WAS MISSING!
    
# FALLBACK: If find_library failed but pkg-config succeeded, use pkg-config's library list
elseif(NETCDF_FOUND AND NETCDF_LINK_LIBRARIES)
    set(NETCDF_LIBRARIES_C ${NETCDF_LINK_LIBRARIES})
    message(STATUS "Using NetCDF libraries from pkg-config: ${NETCDF_LINK_LIBRARIES}")
    list(APPEND NETCDF_DETECTION_LOG "pkg-config fallback: ${NETCDF_LINK_LIBRARIES}")
else()
    message(DEBUG "libnetcdf not found")
    list(APPEND NETCDF_DETECTION_LOG "find_library: failed")
endif()

# HDF5 dependency
message(DEBUG "Checking HDF5 dependency")
if(NETCDF_LIBRARIES_C AND NETCDF_LINK_LIBRARIES)
    string(FIND "${NETCDF_LINK_LIBRARIES}" "hdf5" HDF5_IN_NETCDF)
    if(HDF5_IN_NETCDF GREATER -1)
        message(VERBOSE "NetCDF requires HDF5")
        if(TARGET hdf5::hdf5 OR HDF5_FOUND)
            list(APPEND NETCDF_LIBRARIES_C ${HDF5_LIBRARIES})
            message(DEBUG "Using HDF5: ${HDF5_LIBRARIES}")
        else()
            set(NETCDF_LIBRARIES_C ${NETCDF_LINK_LIBRARIES})
            message(DEBUG "Using pkg-config libraries with HDF5")
        endif()
    endif()
endif()

set(NetCDF_has_interfaces "YES")
set(NetCDF_libs "${NETCDF_LIBRARIES_C}")

get_filename_component(NetCDF_lib_dirs "${NETCDF_LIBRARIES_C}" PATH)

macro(NetCDF_check_interface lang header libs)
    if(NETCDF_${lang})
        find_path(NETCDF_INCLUDES_${lang} NAMES ${header}
            HINTS ${NETCDF_INCLUDE_HINTS}
                  "${NETCDF_INCLUDES}")
        find_library(NETCDF_LIBRARIES_${lang} NAMES ${libs}
            HINTS ${NETCDF_LIBRARY_HINTS}
                  "${NetCDF_lib_dirs}")
        mark_as_advanced(NETCDF_INCLUDES_${lang} NETCDF_LIBRARIES_${lang})

        if(NETCDF_INCLUDES_${lang} AND NETCDF_LIBRARIES_${lang})
            list(INSERT NetCDF_libs 0 ${NETCDF_LIBRARIES_${lang}}) # prepend so that -lnetcdf is last
        else(NETCDF_INCLUDES_${lang} AND NETCDF_LIBRARIES_${lang})
            set(NetCDF_has_interfaces "NO")
            message(STATUS "Failed to find NetCDF interface for ${lang}")
        endif(NETCDF_INCLUDES_${lang} AND NETCDF_LIBRARIES_${lang})
    endif(NETCDF_${lang})
endmacro(NetCDF_check_interface)

NetCDF_check_interface(CXX netcdfcpp.h netcdf_c++)
NetCDF_check_interface(F77 netcdf.inc netcdff)
NetCDF_check_interface(F90 netcdf.mod netcdff)

set(NETCDF_LIBRARIES "${NetCDF_libs}" CACHE STRING "All NetCDF libraries required for interface level")
set(NETCDF_LINK_LIBRARIES ${NetCDF_libs})
set(NETCDF_INCLUDE_DIRS ${NETCDF_INCLUDES})

# Check if detection failed - show helpful error BEFORE standard handler
if(NOT NETCDF_LIBRARIES_C OR NOT NETCDF_INCLUDES)
    message(STATUS "")
    message(STATUS "====================================================================")
    message(STATUS "NetCDF Detection Failed")
    message(STATUS "====================================================================")
    message(STATUS "")
    message(STATUS "Detection attempts:")
    foreach(attempt ${NETCDF_DETECTION_LOG})
        message(STATUS "  ${attempt}")
    endforeach()
    message(STATUS "")
    message(STATUS "Missing components:")
    message(STATUS "  netcdf.h:    ${NETCDF_INCLUDES}")
    message(STATUS "  libnetcdf:   ${NETCDF_LIBRARIES_C}")
    message(STATUS "")
    message(STATUS "To resolve:")
    message(STATUS "")
    message(STATUS "  On Perlmutter/NERSC:")
    message(STATUS "    module load cray-netcdf-hdf5parallel")
    message(STATUS "")
    message(STATUS "  On other Cray systems:")
    message(STATUS "    module load cray-netcdf")
    message(STATUS "")
    message(STATUS "  Or specify manually:")
    message(STATUS "    cmake -DNETCDF_DIR=/path/to/netcdf ..")
    message(STATUS "")
    message(STATUS "  Or via environment:")
    message(STATUS "    export NETCDF_DIR=/path/to/netcdf")
    message(STATUS "")
    message(STATUS "====================================================================")
    message(STATUS "")
    message(FATAL_ERROR "NetCDF not found")
endif()

# Standard find package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args (NetCDF DEFAULT_MSG NETCDF_LIBRARIES NETCDF_LINK_LIBRARIES NETCDF_INCLUDE_DIRS NETCDF_INCLUDES NetCDF_has_interfaces)

# Show diagnostics on failure
if(NOT NETCDF_FOUND)
    message(STATUS "Detection attempts:")
    foreach(attempt ${NETCDF_DETECTION_LOG})
        message(STATUS "  ${attempt}")
    endforeach()
    message(STATUS "")
    message(STATUS "To resolve:")
    message(STATUS "  Cray:      module load cray-netcdf-hdf5parallel")
    message(STATUS "  Manual:    -DNETCDF_DIR=/path/to/netcdf")
    message(STATUS "  Env var:   export NETCDF_DIR=/path/to/netcdf")
endif()

mark_as_advanced(NETCDF_LIBRARIES NETCDF_INCLUDES)

# Pop FindNetCDF context
list(POP_BACK CMAKE_MESSAGE_CONTEXT)
