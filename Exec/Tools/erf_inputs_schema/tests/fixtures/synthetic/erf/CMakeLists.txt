option(ERF_ENABLE_WIDGET "Build the widget" OFF)
if(ERF_ENABLE_WIDGET)
  target_compile_definitions(erf_srclib PUBLIC ERF_USE_WIDGET)
endif()
