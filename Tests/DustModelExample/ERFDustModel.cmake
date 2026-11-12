# CMake: the dust model directory's sources. ERF includes this file for each
# directory listed in ERF_DUST_MODEL_DIRS and calls erf_dust_model_add with
# ERF's OBJECT library; ERF_DUST_MODEL_DIR holds this directory during the
# call. A model that links a prebuilt library adds it here too.
function(erf_dust_model_add target)
  target_sources(${target} PRIVATE ${ERF_DUST_MODEL_DIR}/ERF_DustModelExample.cpp)
  target_include_directories(${target} PRIVATE ${ERF_DUST_MODEL_DIR})
endfunction()
