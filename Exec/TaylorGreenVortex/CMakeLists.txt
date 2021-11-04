set(erf_exe_name taylor_green)

add_library(${erf_exe_name} "")
add_executable(${erf_exe_name}_exe "")
target_sources(${erf_exe_name}_exe
   PRIVATE
     prob.H
     prob.cpp
)

target_include_directories(${erf_exe_name}_exe PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(${erf_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

include(${CMAKE_SOURCE_DIR}/CMake/BuildERFExe.cmake)
build_erf_exe(${erf_exe_name})
