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

find_package(PkgConfig REQUIRED QUIET)
pkg_check_modules(NETCDF REQUIRED IMPORTED_TARGET netcdf)

# handle the QUIETLY and REQUIRED arguments and set NETCDF_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (NetCDF DEFAULT_MSG NETCDF_LIBRARIES NETCDF_LINK_LIBRARIES NETCDF_INCLUDE_DIRS)

mark_as_advanced (NETCDF_LIBRARIES NETCDF_INCLUDES)

