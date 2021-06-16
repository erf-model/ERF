function(build_erf_exe erf_exe_name)

  set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)
  set(BIN_DIR ${CMAKE_BINARY_DIR}/Source/${erf_exe_name})

  include(${CMAKE_SOURCE_DIR}/CMake/SetERFCompileFlags.cmake)

  add_subdirectory(${SRC_DIR}/Params ${BIN_DIR}/Params/${erf_exe_name})

  set(ERF_EOS_DIR "${CMAKE_SOURCE_DIR}/Source")
  target_sources(${erf_exe_name} PRIVATE
                 ${ERF_EOS_DIR}/EOS.H)
  target_include_directories(${erf_exe_name} SYSTEM PRIVATE ${ERF_EOS_DIR})

  if(ERF_ENABLE_MASA)
    target_sources(${erf_exe_name} PRIVATE
                   ${SRC_DIR}/MMS.cpp)
    target_compile_definitions(${erf_exe_name} PRIVATE ERF_USE_MASA)
  endif()
  
  target_sources(${erf_exe_name}
     PRIVATE
       ${SRC_DIR}/Advance.cpp
       ${SRC_DIR}/BCfill.cpp
       ${SRC_DIR}/Bld.cpp
       ${SRC_DIR}/Constants.H
       ${SRC_DIR}/Derive.H
       ${SRC_DIR}/Derive.cpp
       ${SRC_DIR}/External.cpp
       ${SRC_DIR}/Forcing.H
       ${SRC_DIR}/Forcing.cpp
#      ${SRC_DIR}/GradUtil.H
       ${SRC_DIR}/IndexDefines.H
#      ${SRC_DIR}/IndexDefines.cpp
       ${SRC_DIR}/IO.H
       ${SRC_DIR}/IO.cpp
       ${SRC_DIR}/Transport.H
       ${SRC_DIR}/TransportParams.H
       ${SRC_DIR}/ERF.H
       ${SRC_DIR}/ERF.cpp
       ${SRC_DIR}/Problem.H
       ${SRC_DIR}/ProblemDerive.H
       ${SRC_DIR}/Setup.cpp
       ${SRC_DIR}/Sources.cpp
       ${SRC_DIR}/SumIQ.cpp
       ${SRC_DIR}/SumUtils.cpp
       ${SRC_DIR}/Tagging.H
       ${SRC_DIR}/Tagging.cpp
       ${SRC_DIR}/Utilities.H
       ${SRC_DIR}/Transport.cpp
       ${SRC_DIR}/TransportParams.cpp
       ${SRC_DIR}/RK3/RK3.H
       ${SRC_DIR}/RK3/CalcAdvFlux.cpp
       ${SRC_DIR}/RK3/CalcDiffFlux.cpp
       ${SRC_DIR}/RK3/MomentumToVelocity.cpp
       ${SRC_DIR}/RK3/VelocityToMomentum.cpp
       ${SRC_DIR}/RK3/RK3_driver.cpp
       ${SRC_DIR}/RK3/RK3_stage.cpp
       ${SRC_DIR}/RK3/RK3_utils.cpp
       ${SRC_DIR}/RK3/Interpolation.cpp
       ${SRC_DIR}/RK3/Advection.cpp
       ${SRC_DIR}/RK3/Diffusion.cpp
       ${SRC_DIR}/RK3/StrainRate.cpp
       ${SRC_DIR}/RK3/EddyViscosity.cpp
  )

  if(NOT "${erf_exe_name}" STREQUAL "erf_unit_tests")
    target_sources(${erf_exe_name}
       PRIVATE
         ${SRC_DIR}/main.cpp
    )
  endif()
  
  include(AMReXBuildInfo)
  generate_buildinfo(${erf_exe_name} ${CMAKE_SOURCE_DIR})
  target_include_directories(${erf_exe_name} PUBLIC ${AMREX_SUBMOD_LOCATION}/Tools/C_scripts)

  if(ERF_ENABLE_MASA)
    if(MASA_FOUND)
      #Link our executable to the MASA libraries, etc
      target_link_libraries(${erf_exe_name} PRIVATE ${MASA_LIBRARY})
      target_compile_definitions(${erf_exe_name} PRIVATE USE_MASA DO_PROBLEM_POST_TIMESTEP DO_PROBLEM_POST_INIT)
      target_include_directories(${erf_exe_name} SYSTEM PRIVATE ${MASA_INCLUDE_DIRS})
      target_include_directories(${erf_exe_name} SYSTEM PRIVATE ${MASA_MOD_DIRS})
    endif()
  endif()

  if(ERF_ENABLE_MPI)
    target_link_libraries(${erf_exe_name} PUBLIC $<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_CXX>)
  endif()

  #ERF include directories
  target_include_directories(${erf_exe_name} PRIVATE ${SRC_DIR})
  target_include_directories(${erf_exe_name} PRIVATE ${SRC_DIR}/RK3)
  target_include_directories(${erf_exe_name} PRIVATE ${CMAKE_BINARY_DIR})
  
  #Link to amrex library
  target_link_libraries(${erf_exe_name} PRIVATE amrex)

  if(ERF_ENABLE_CUDA)
    set(pctargets "${erf_exe_name}")
    foreach(tgt IN LISTS pctargets)
      get_target_property(ERF_SOURCES ${tgt} SOURCES)
      list(FILTER ERF_SOURCES INCLUDE REGEX "\\.cpp")
      set_source_files_properties(${ERF_SOURCES} PROPERTIES LANGUAGE CUDA)
    endforeach()
  endif()
 
  #Define what we want to be installed during a make install 
  install(TARGETS ${erf_exe_name}
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)

endfunction()
