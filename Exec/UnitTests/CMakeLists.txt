set(erf_exe_name erf_unit_tests)

add_executable(${erf_exe_name} "")
target_sources(${erf_exe_name}
  PUBLIC
  unit-tests-main.cpp
  test-config.cpp
  prob.cpp
  prob.H
  prob_parm.H
  )

target_include_directories(${erf_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

include(${CMAKE_SOURCE_DIR}/CMake/BuildERFExe.cmake)
build_erf_exe(${erf_exe_name})

target_include_directories(${erf_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/GoogleTest/googletest/include)
target_link_libraries(${erf_exe_name} PRIVATE gtest)
