# CMake: the flood engine directory's sources. ERF includes this file for each
# directory listed in ERF_FLOOD_ENGINE_DIRS and calls erf_flood_engine_add with
# ERF's OBJECT library; ERF_FLOOD_ENGINE_DIR holds this directory during the call.
# An engine that links a prebuilt library adds it here too, for example
#   target_link_libraries(${target} PUBLIC ${MY_FLOOD_LIB})
function(erf_flood_engine_add target)
  target_sources(${target} PRIVATE ${ERF_FLOOD_ENGINE_DIR}/ERF_FloodEngineExample.cpp)
  target_include_directories(${target} PRIVATE ${ERF_FLOOD_ENGINE_DIR})
endfunction()

# Optional: the directory's GoogleTest files. With ERF_ENABLE_UNIT_TESTS, Tests/Unit
# includes this file again and calls erf_flood_engine_add_unit_tests with ERF's
# unit-test executable (erf_unit_tests), which has ERF's test main and links GoogleTest
# and the library above, so the tests can create the engine by its registered name.
function(erf_flood_engine_add_unit_tests target)
  target_sources(${target} PRIVATE ${ERF_FLOOD_ENGINE_DIR}/ERF_GTestFloodEngineExample.cpp)
endfunction()
