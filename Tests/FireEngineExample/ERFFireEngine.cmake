# CMake: the fire engine directory's sources. ERF includes this file for each
# directory listed in ERF_FIRE_ENGINE_DIRS and calls erf_fire_engine_add with
# ERF's OBJECT library; ERF_FIRE_ENGINE_DIR holds this directory during the
# call. An engine that links a prebuilt library adds it here too, for example
#   target_link_libraries(${target} PUBLIC ${ELMFIRE_LIB})
function(erf_fire_engine_add target)
  target_sources(${target} PRIVATE ${ERF_FIRE_ENGINE_DIR}/ERF_FireEngineExample.cpp)
  target_include_directories(${target} PRIVATE ${ERF_FIRE_ENGINE_DIR})
endfunction()
