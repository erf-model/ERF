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

# Build hints from user variables first, then pkg-config
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

set(ENV{PKG_CONFIG_PATH} "$ENV{MPICH_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
message(STATUS "PKG_CONFIG_PATH = $ENV{PKG_CONFIG_PATH}")

find_package(PkgConfig REQUIRED QUIET)
pkg_check_modules(NETCDF QUIET IMPORTED_TARGET netcdf)
if(NOT NETCDF_FOUND)
    pkg_check_modules(NETCDF REQUIRED IMPORTED_TARGET netcdf-cxx4_parallel)
endif()

# Add pkg-config results to hints
list(APPEND NETCDF_INCLUDE_HINTS ${NETCDF_INCLUDE_DIRS})
list(APPEND NETCDF_LIBRARY_HINTS ${NETCDF_LIBRARY_DIRS})

find_path(NETCDF_INCLUDES netcdf.h
    HINTS ${NETCDF_INCLUDE_HINTS}
          $ENV{NETCDF_DIR}/include)

find_library(NETCDF_LIBRARIES_C NAMES netcdf
    HINTS ${NETCDF_LIBRARY_HINTS}
          $ENV{NETCDF_DIR}/lib)
mark_as_advanced(NETCDF_LIBRARIES_C)

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

