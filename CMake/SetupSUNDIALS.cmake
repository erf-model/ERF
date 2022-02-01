# Note this setup file expects to be used after AMReX options are set and before any targets are linked
set(SUNDIALS_MINIMUM_VERSION 6.0 CACHE INTERNAL "Minimum required SUNDIALS version")

# We first check if we can find an AMReX installation.
# If so, we proceed with STANDALONE mode
# If not, we proceed with SUPERBUILD MODE
find_package( SUNDIALS ${SUNDIALS_MINIMUM_VERSION} CONFIG QUIET )

#
# ~~~~~~~~~~~~~~~~~~~~~~~ STANDALONE MODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
if (SUNDIALS_FOUND)

   # Unfortunately SUNDIALS config file doesn't seem to allow to check for
   # components, so we assume that the SUNDIALS installation, if found,
   # has all the necessary components
   message(STATUS "SUNDIALS found: configuration file located at ${SUNDIALS_DIR}")

else ()
#
# ~~~~~~~~~~~~~~~~~~~~~~~ SUPERBUILD MODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

   message(STATUS "No SUNDIALS installation found: cloning SUNDIALS repo")

   if (NOT EXISTS  "${PROJECT_SOURCE_DIR}/.git")
      message(FATAL_ERROR
         "${PROJECT_SOURCE_DIR} is not a Git repo: missing .git")
   endif ()

   set(SUNDIALS_SRC_DIR "${PROJECT_SOURCE_DIR}/Submodules/SUNDIALS"
     CACHE INTERNAL "Path to SUNDIALS source (submodule)")

   if (NOT EXISTS "${SUNDIALS_SRC_DIR}/.git")
      message(STATUS "Initializing git submodule for SUNDIALS")

      find_package(Git REQUIRED)

      execute_process(
         COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive Submodules/SUNDIALS
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
         RESULT_VARIABLE GIT_SUBMOD_RESULT
         )

      if ( NOT GIT_SUBMOD_RESULT EQUAL "0")
         message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
      endif()

      unset(GIT_SUBMOD_RESULT)

   endif ()

   # Set build options for subproject
   set(EXAMPLES_ENABLE_C            OFF                        CACHE INTERNAL "" )
   set(EXAMPLES_ENABLE_CXX          OFF                        CACHE INTERNAL "" )
   set(EXAMPLES_ENABLE_CUDA         OFF                        CACHE INTERNAL "" )
   set(EXAMPLES_ENABLE_F77          OFF                        CACHE INTERNAL "" )
   set(EXAMPLES_ENABLE_F90          OFF                        CACHE INTERNAL "" )
   set(EXAMPLES_ENABLE_F2003        OFF                        CACHE INTERNAL "" )
   set(EXAMPLES_INSTALL             OFF                        CACHE INTERNAL "" )
   set(ENABLE_MPI                   OFF                        CACHE INTERNAL "" )
   set(ENABLE_OPENMP                ${ERF_ENABLE_OMP}                 CACHE INTERNAL "" )
   if (AMReX_GPU_BACKEND STREQUAL CUDA)
     set(ENABLE_CUDA                  ON                      CACHE INTERNAL "" )
     set(SUNDIALS_INDEX_SIZE          32                      CACHE INTERNAL "" )
     set(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS     ON          CACHE INTERNAL "" )
    else ()
      set(ENABLE_CUDA                  OFF                     CACHE INTERNAL "" )
      if (AMReX_GPU_BACKEND STREQUAL HIP)
	set(ENABLE_HIP                   ON                      CACHE INTERNAL "" )
	set(AMDGPU_TARGETS               ${AMReX_AMD_ARCH}       CACHE INTERNAL "" )
#        set(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS     ON          CACHE INTERNAL "" )
      endif ()
      if (AMReX_GPU_BACKEND STREQUAL SYCL)
	set(ENABLE_SYCL                  ON                      CACHE INTERNAL "" )
	set(CMAKE_CXX_STANDARD           17                      CACHE INTERNAL "" )
      endif ()
   endif ()
   set(BUILD_ARKODE                 ON                         CACHE INTERNAL "" )
   set(BUILD_KINSOL                 OFF                        CACHE INTERNAL "" )
   set(BUILD_IDA                    OFF                        CACHE INTERNAL "" )
   set(BUILD_IDAS                   OFF                        CACHE INTERNAL "" )
   set(BUILD_CVODES                 OFF                        CACHE INTERNAL "" )

   #  Add  SUNDIALS sources to the build
   add_subdirectory(${SUNDIALS_SRC_DIR})

   # This is to use the same target name uses by the sundials exported targets
   add_library(SUNDIALS::cvode      ALIAS sundials_cvode_static)
   add_library(SUNDIALS::arkode     ALIAS sundials_arkode_static)
   add_library(SUNDIALS::nvecmanyvector ALIAS sundials_nvecmanyvector_static)

   if (AMReX_GPU_BACKEND STREQUAL CUDA)
      add_library(SUNDIALS::nveccuda ALIAS sundials_nveccuda_static)
      install(TARGETS sundials_nveccuda_static
      EXPORT sundials)
   elseif (AMReX_GPU_BACKEND STREQUAL HIP)
      add_library(SUNDIALS::nvechip ALIAS sundials_nvechip_static)
      install(TARGETS sundials_nvechip_static
      EXPORT sundials)
   elseif (AMReX_GPU_BACKEND STREQUAL SYCL)
      add_library(SUNDIALS::nvecsycl ALIAS sundials_nvecsycl_static)
      install(TARGETS sundials_nvecsycl_static
      EXPORT sundials)
   else ()
     add_library(SUNDIALS::nvecserial ALIAS sundials_nvecserial_static)
     install(TARGETS sundials_nvecserial_static
     EXPORT sundials)
   endif ()

   #export(EXPORT sundials)
   set(SUNDIALS_FOUND TRUE)

endif ()
