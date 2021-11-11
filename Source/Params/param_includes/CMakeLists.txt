target_sources(${erf_lib_name}
   PRIVATE
     erf_defaults.H
     erf_params.H
     erf_queries.H
   )

target_include_directories(${erf_lib_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
