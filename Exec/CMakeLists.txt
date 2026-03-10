set(erf_lib_name erf_srclib)

add_library(${erf_lib_name} OBJECT)
if (BUILD_SHARED_LIBS)
  set_target_properties(${erf_lib_name} PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

include(${CMAKE_SOURCE_DIR}/CMake/BuildERFExe.cmake)

build_erf_lib(${erf_lib_name})

if(NOT ERF_BUILD_LIBRARY_ONLY)
  add_subdirectory(RegTests)
  add_subdirectory(DevTests/RICO)
  add_subdirectory(DevTests/Shoc)
  add_subdirectory(DevTests/TemperatureSourceSpatial)
  add_subdirectory(DevTests/TropicalCyclone)
  add_subdirectory(DevTests/Tornado)
  add_subdirectory(DevTests/sinusoidal_mass_flux)
  add_subdirectory(DevTests/WindFarmTests)

if(ERF_ENABLE_TOOLS)
  message(STATUS "Building additional ERF tools")
  add_subdirectory(Tools)
endif()

endif() # NOT ERF_BUILD_LIBRARY_ONLY
