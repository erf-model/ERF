set(erf_exe_name erf_super_cell)

add_executable(${erf_exe_name} "")
target_sources(${erf_exe_name}
   PRIVATE
     ERF_prob.cpp
)

target_include_directories(${erf_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

include(${CMAKE_SOURCE_DIR}/CMake/BuildERFExe.cmake)
build_erf_exe(${erf_exe_name})
