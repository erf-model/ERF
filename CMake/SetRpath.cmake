# Use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# Add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Apply the RPATH to binaries in the build tree (don't wait for install)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

# The RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

if(UNIX AND NOT APPLE)
    # On Cray PE without CRAY_ADD_RPATH=yes, the wrappers won't inject RPATH.
    # Explicitly bake all module library paths into CMAKE_INSTALL_RPATH so
    # CMake cannot strip them as "implicit system paths".
    if(DEFINED ENV{CRAYPE_VERSION} AND NOT "$ENV{CRAY_ADD_RPATH}" MATCHES "^(yes|YES|1)$")
        # (label env_var subdir) triples — subdir is relative to the library root
        set(_erf_rpath_specs
            NETCDF           NETCDF_DIR          lib
            HDF5             HDF5_DIR            lib
            HDF5             HDF5_DIR            lib64
            HDF5             HDF5_ROOT           lib
            HDF5             HDF5_ROOT           lib64
            MPI              CRAY_MPICH_DIR      lib
            CUDA             CUDA_HOME           lib64
            CUDA             CUDATOOLKIT_HOME    lib64
            ROCm             ROCM_PATH           lib
            ROCm             ROCM_PATH           hip/lib
            FFTW             FFTW_DIR            lib
            FFTW             CRAY_FFTW_DIR       lib
        )
        list(LENGTH _erf_rpath_specs _n)
        math(EXPR _last "${_n} - 1")
        foreach(_i RANGE 0 ${_last} 3)
            math(EXPR _j "${_i} + 1")
            math(EXPR _k "${_i} + 2")
            list(GET _erf_rpath_specs ${_j} _env)
            list(GET _erf_rpath_specs ${_k} _subdir)
            if(DEFINED ENV{${_env}})
                list(APPEND CMAKE_INSTALL_RPATH "$ENV{${_env}}/${_subdir}")
            endif()
        endforeach()
        unset(_erf_rpath_specs)
        unset(_i)
        unset(_j)
        unset(_k)
        unset(_n)
        unset(_last)
        unset(_env)
        unset(_subdir)
    endif()

    # Force DT_RPATH (not DT_RUNPATH) so transitive dependencies (e.g., HDF5
    # via NetCDF) are resolved on non-Cray Linux and on Cray without wrappers.
    if(NOT DEFINED ENV{CRAYPE_VERSION} OR NOT "$ENV{CRAY_ADD_RPATH}" MATCHES "^(yes|YES|1)$")
        add_link_options("-Wl,--disable-new-dtags")
    endif()
endif()
