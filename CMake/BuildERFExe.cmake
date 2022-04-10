function(target_link_libraries_system target visibility)
  set(libs ${ARGN})
  foreach(lib ${libs})
    get_target_property(lib_include_dirs ${lib} INTERFACE_INCLUDE_DIRECTORIES)
    target_include_directories(${target} SYSTEM ${visibility} ${lib_include_dirs})
    target_link_libraries(${target} ${visibility} ${lib})
  endforeach(lib)
endfunction(target_link_libraries_system)

function(build_erf_lib erf_lib_name)

  set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)
  set(BIN_DIR ${CMAKE_BINARY_DIR}/Source/${erf_lib_name})

  include(${CMAKE_SOURCE_DIR}/CMake/SetERFCompileFlags.cmake)
  set_erf_compile_flags(${erf_lib_name})

  set(ERF_EOS_DIR "${CMAKE_SOURCE_DIR}/Source")
  target_sources(${erf_lib_name} PRIVATE
                 ${ERF_EOS_DIR}/EOS.H)
  target_include_directories(${erf_lib_name} SYSTEM PUBLIC ${ERF_EOS_DIR})

  if(ERF_ENABLE_NETCDF)
    target_sources(${erf_lib_name} PRIVATE
                   ${SRC_DIR}/IO/NCInterface.H
                   ${SRC_DIR}/IO/NCWpsFile.H
                   ${SRC_DIR}/IO/NCPlotFile.H
                   ${SRC_DIR}/IO/NCBuildFABs.cpp
                   ${SRC_DIR}/IO/NCInterface.cpp
                   ${SRC_DIR}/IO/NCPlotFile.cpp
                   ${SRC_DIR}/IO/NCCheckpoint.cpp
                   ${SRC_DIR}/IO/NCMultiFabFile.cpp
                   ${SRC_DIR}/IO/NCColumnFile.cpp)
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_NETCDF)
  endif()
 
  target_sources(${erf_lib_name}
     PRIVATE
       ${SRC_DIR}/DataStruct.H
       ${SRC_DIR}/ERF_Constants.H
       ${SRC_DIR}/Derive.H
       ${SRC_DIR}/Derive.cpp
       ${SRC_DIR}/IndexDefines.H
       ${SRC_DIR}/Metrics.cpp
       ${SRC_DIR}/prob_common.H
       ${SRC_DIR}/ERF.H
       ${SRC_DIR}/ERF.cpp
       ${SRC_DIR}/ERF_init1d.cpp
       ${SRC_DIR}/ERF_SumIQ.cpp
       ${SRC_DIR}/ERF_Tagging.cpp
       ${SRC_DIR}/ERF_ComputeTimestep.cpp
       ${SRC_DIR}/ERF_FillPatch.cpp
       ${SRC_DIR}/ERF_TimeStepping.cpp
       ${SRC_DIR}/IO/Checkpoint.cpp
       ${SRC_DIR}/IO/ERF_ReadBndryPlanes.H
       ${SRC_DIR}/IO/ERF_ReadBndryPlanes.cpp
       ${SRC_DIR}/IO/ERF_WriteBndryPlanes.H
       ${SRC_DIR}/IO/ERF_WriteBndryPlanes.cpp
       ${SRC_DIR}/IO/Plotfile.cpp
       ${SRC_DIR}/IO/writeJobInfo.cpp
       ${SRC_DIR}/SpatialStencils/Advection.cpp
       ${SRC_DIR}/SpatialStencils/Diffusion.cpp
       ${SRC_DIR}/SpatialStencils/EddyViscosity.H
       ${SRC_DIR}/SpatialStencils/ExpansionRate.H
       ${SRC_DIR}/SpatialStencils/StrainRate.H
       ${SRC_DIR}/SpatialStencils/StressTerm.H
       ${SRC_DIR}/SpatialStencils/Interpolation.cpp
       ${SRC_DIR}/SpatialStencils/MomentumToVelocity.cpp
       ${SRC_DIR}/SpatialStencils/VelocityToMomentum.cpp
       ${SRC_DIR}/TimeIntegration/ERF_MRI.H
       ${SRC_DIR}/TimeIntegration/TimeIntegration.H
       ${SRC_DIR}/TimeIntegration/TimeIntegration_driver.cpp
       ${SRC_DIR}/TimeIntegration/TimeIntegration_rhs.cpp
       ${SRC_DIR}/TimeIntegration/TimeIntegration_fast_implicit.cpp
       ${SRC_DIR}/WPS_Interface.cpp
  )

  if(NOT "${erf_exe_name}" STREQUAL "erf_unit_tests")
    target_sources(${erf_lib_name}
       PRIVATE
         ${SRC_DIR}/main.cpp
    )
  endif()

  include(AMReXBuildInfo)
  generate_buildinfo(${erf_lib_name} ${CMAKE_SOURCE_DIR})
  target_include_directories(${erf_lib_name} PUBLIC ${AMREX_SUBMOD_LOCATION}/Tools/C_scripts)

  if(ERF_ENABLE_NETCDF)
    if(NETCDF_FOUND)
      #Link our executable to the NETCDF libraries, etc
      target_link_libraries(${erf_lib_name} PUBLIC ${NETCDF_LIBRARIES_C})
      target_include_directories(${erf_lib_name} PUBLIC ${NETCDF_INCLUDES})
    endif()
  endif()

  if(ERF_ENABLE_MPI)
    target_link_libraries(${erf_lib_name} PUBLIC $<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_CXX>)
  endif()

  #ERF include directories
  target_include_directories(${erf_lib_name} PUBLIC ${SRC_DIR})
  target_include_directories(${erf_lib_name} PUBLIC ${SRC_DIR}/SpatialStencils)
  target_include_directories(${erf_lib_name} PUBLIC ${SRC_DIR}/TimeIntegration)
  target_include_directories(${erf_lib_name} PUBLIC ${SRC_DIR}/IO)
  target_include_directories(${erf_lib_name} PUBLIC ${CMAKE_BINARY_DIR})

  if (AMReX_SUNDIALS)
    target_link_libraries(${erf_lib_name} PUBLIC SUNDIALS::cvode)
    target_link_libraries(${erf_lib_name} PUBLIC SUNDIALS::arkode)
    target_link_libraries(${erf_lib_name} PUBLIC SUNDIALS::nvecmanyvector)
  endif ()

  #Link to amrex library
  target_link_libraries_system(${erf_lib_name} PUBLIC amrex)
  if(ERF_ENABLE_CUDA)
    set(pctargets "${erf_lib_name}")
    foreach(tgt IN LISTS pctargets)
      get_target_property(ERF_SOURCES ${tgt} SOURCES)
      list(FILTER ERF_SOURCES INCLUDE REGEX "\\.cpp")
      set_source_files_properties(${ERF_SOURCES} PROPERTIES LANGUAGE CUDA)
      message(STATUS "setting cuda for ${ERF_SOURCES}")
    endforeach()
    set_target_properties(
    ${erf_lib_name} PROPERTIES
    LANGUAGE CUDA
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  endif()

  #Define what we want to be installed during a make install 
  install(TARGETS ${erf_lib_name}
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)

endfunction(build_erf_lib)

function(build_erf_exe erf_exe_name)

  set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)

  target_link_libraries(${erf_exe_name}  PUBLIC ${erf_lib_name})
  include(${CMAKE_SOURCE_DIR}/CMake/SetERFCompileFlags.cmake)
  set_erf_compile_flags(${erf_exe_name})

  target_sources(${erf_exe_name}
     PRIVATE
       ${SRC_DIR}/ERF_init_bcs.cpp

  )
  if(ERF_ENABLE_CUDA)
    set(pctargets "${erf_exe_name}")
    foreach(tgt IN LISTS pctargets)
      get_target_property(ERF_SOURCES ${tgt} SOURCES)
      list(FILTER ERF_SOURCES INCLUDE REGEX "\\.cpp")
      set_source_files_properties(${ERF_SOURCES} PROPERTIES LANGUAGE CUDA)
      message(STATUS "setting cuda for ${ERF_SOURCES}")
    endforeach()
    set_target_properties(
    ${erf_exe_name} PROPERTIES
    LANGUAGE CUDA
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  endif()

  install(TARGETS ${erf_exe_name} 
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)

endfunction()
