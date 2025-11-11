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

if (NETCDF_INCLUDES AND NETCDF_LIBRARIES)
  # Already in cache, be silent
  set (NETCDF_FIND_QUIETLY TRUE)
endif (NETCDF_INCLUDES AND NETCDF_LIBRARIES)

# Build hints from user variables first
set(NETCDF_INCLUDE_HINTS)
set(NETCDF_LIBRARY_HINTS)

if(NETCDF_DIR)
    list(APPEND NETCDF_INCLUDE_HINTS ${NETCDF_DIR}/include)
    list(APPEND NETCDF_LIBRARY_HINTS ${NETCDF_DIR}/lib)
endif()

if(NETCDF_INCLUDE_DIR)
    list(APPEND NETCDF_INCLUDE_HINTS ${NETCDF_INCLUDE_DIR})
endif()

if(NETCDF_LIBRARY_DIR)
    list(APPEND NETCDF_LIBRARY_HINTS ${NETCDF_LIBRARY_DIR})
endif()

# Use pkg-config to get hints
set(ENV{PKG_CONFIG_PATH} "$ENV{MPICH_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
message(STATUS "PKG_CONFIG_PATH = $ENV{PKG_CONFIG_PATH}")

find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
    # Try multiple NetCDF variants in order of preference
    pkg_check_modules(NETCDF QUIET IMPORTED_TARGET netcdf)
    if(NOT NETCDF_FOUND)
        pkg_check_modules(NETCDF QUIET IMPORTED_TARGET netcdf-mpi)
    endif()
    if(NOT NETCDF_FOUND)
        pkg_check_modules(NETCDF QUIET IMPORTED_TARGET netcdf_parallel)
    endif()
    if(NOT NETCDF_FOUND)
        pkg_check_modules(NETCDF QUIET IMPORTED_TARGET netcdf-cxx4_parallel)
    endif()

    if(NETCDF_FOUND)
        message(STATUS "Found NetCDF via pkg-config: ${NETCDF_MODULE_NAME}")
        # Add pkg-config results to hints
        list(APPEND NETCDF_INCLUDE_HINTS ${NETCDF_INCLUDE_DIRS})
        list(APPEND NETCDF_LIBRARY_HINTS ${NETCDF_LIBRARY_DIRS})
    endif()
endif()

# Try CMake's find_library using hints
find_path(NETCDF_INCLUDES netcdf.h
    HINTS ${NETCDF_INCLUDE_HINTS}
          $ENV{NETCDF_DIR}/include)

find_library(NETCDF_LIBRARIES_C NAMES netcdf
    HINTS ${NETCDF_LIBRARY_HINTS}
          $ENV{NETCDF_DIR}/lib)
mark_as_advanced(NETCDF_LIBRARIES_C)

# If find_library succeeded, check if we need HDF5
if(NETCDF_LIBRARIES_C)
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
            message(STATUS "NetCDF was built without HDF5 support")
        endif()
    else()
        message(STATUS "No pkg-config information available; assuming NetCDF doesn't need HDF5")
    endif()
# FALLBACK: If find_library failed but pkg-config succeeded, use pkg-config's library list
elseif(NETCDF_FOUND AND NETCDF_LINK_LIBRARIES)
    set(NETCDF_LIBRARIES_C ${NETCDF_LINK_LIBRARIES})
    message(STATUS "Using NetCDF libraries from pkg-config: ${NETCDF_LINK_LIBRARIES}")
endif()

set(NetCDF_has_interfaces "YES") # will be set to NO if we're missing any interfaces
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

# handle the QUIETLY and REQUIRED arguments and set NETCDF_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (NetCDF DEFAULT_MSG NETCDF_LIBRARIES NETCDF_LINK_LIBRARIES NETCDF_INCLUDE_DIRS NETCDF_INCLUDES NetCDF_has_interfaces)

message(STATUS "  NETCDF_LIBRARIES = ${NETCDF_LIBRARIES}")
#message(STATUS "  NETCDF_LINK_LIBRARIES = ${NETCDF_LINK_LIBRARIES}")
#message(STATUS "  NETCDF_INCLUDE_DIRS = ${NETCDF_INCLUDE_DIRS}")
message(STATUS "  NETCDF_INCLUDES = ${NETCDF_INCLUDES}")
mark_as_advanced (NETCDF_LIBRARIES NETCDF_INCLUDES)

