target_sources(${erf_exe_name}
   PRIVATE
     erf_defaults.H
     erf_params.H
     erf_queries.H 
   )

target_include_directories(${erf_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
