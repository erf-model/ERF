# CMake: the structure model directory's sources. ERF includes this file for each
# directory listed in ERF_STRUCTURE_MODEL_DIRS and calls erf_structure_model_add with
# ERF's OBJECT library; ERF_STRUCTURE_MODEL_DIR holds this directory during the call.
# A model that links a prebuilt library adds it here too, for example
#   target_link_libraries(${target} PUBLIC ${MOORDYN_LIB})
function(erf_structure_model_add target)
  target_sources(${target} PRIVATE ${ERF_STRUCTURE_MODEL_DIR}/ERF_StructureModelExample.cpp)
  target_include_directories(${target} PRIVATE ${ERF_STRUCTURE_MODEL_DIR})
endfunction()

# Optional: the directory's GoogleTest files. With ERF_ENABLE_UNIT_TESTS, Tests/Unit
# includes this file again and calls erf_structure_model_add_unit_tests with ERF's
# unit-test executable (erf_unit_tests), which has ERF's test main and links GoogleTest
# and the library above, so the tests can create the model by its registered name.
function(erf_structure_model_add_unit_tests target)
  target_sources(${target} PRIVATE ${ERF_STRUCTURE_MODEL_DIR}/ERF_GTestStructureModelExample.cpp)
endfunction()
