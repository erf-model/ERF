set(erf_exe_name erf_witch_of_agnesi)

add_executable(${erf_exe_name} "")
target_sources(${erf_exe_name}
   PRIVATE
     prob.cpp
)

target_include_directories(${erf_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

include(${CMAKE_SOURCE_DIR}/CMake/BuildERFExe.cmake)
build_erf_exe(${erf_exe_name})

#find_package(AMReX)
#target_link_libraries( ${_target}  AMReX::amrex)
