set(erf_lib_name erf_srclib)

add_library(${erf_lib_name} OBJECT)
if (BUILD_SHARED_LIBS)
  set_target_properties(${erf_lib_name} PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

include(${PROJECT_SOURCE_DIR}/CMake/BuildERFExe.cmake)

build_erf_lib(${erf_lib_name})
# The ERF library under the name a parent project links, as ERF links AMReX::amrex
add_library(ERF::erf ALIAS ${erf_lib_name})

# This compile-time capability exists only in a separately configured
# qualification build.  The normal executable retains the max_level == 0
# guard and does not recognize the internal qualification input switch.
if(ERF_ENABLE_SBM_AMR_QUALIFICATION)
  target_compile_definitions(${erf_lib_name} PUBLIC ERF_SBM_AMR_QUALIFICATION=1)
  target_sources(${erf_lib_name} PRIVATE
    ${PROJECT_SOURCE_DIR}/Source/Microphysics/SBM/ERF_SBMAMRQualification.cpp)
endif()

if(NOT ERF_BUILD_LIBRARY_ONLY)
  # Build the shared regtests executable from Exec so it can be used by
  # both RegTests and CanonicalTests input decks.
  set(erf_exe_name erf_exec)
  add_executable(${erf_exe_name} "")
  target_sources(${erf_exe_name}
     PRIVATE
       ${PROJECT_SOURCE_DIR}/Exec/ERF_Prob.cpp
  )
  target_include_directories(${erf_exe_name} PRIVATE ${PROJECT_SOURCE_DIR}/Exec)
  build_erf_exe(${erf_exe_name})

  add_subdirectory(${PROJECT_SOURCE_DIR}/.Exec_dev ${PROJECT_BINARY_DIR}/Exec_dev)
if(ERF_ENABLE_TOOLS)
  message(STATUS "Building additional ERF tools")
  add_subdirectory(Tools)
endif()

endif() # NOT ERF_BUILD_LIBRARY_ONLY
