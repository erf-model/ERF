set(erf_exe_name erf_noahmp)
add_executable(${erf_exe_name} "")
include(${CMAKE_SOURCE_DIR}/CMake/BuildERFExe.cmake)
build_erf_exe(${erf_exe_name})
