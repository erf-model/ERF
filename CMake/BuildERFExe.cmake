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

  target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_MOISTURE)

  if(ERF_ENABLE_MULTIBLOCK)
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_MULTIBLOCK)
  endif()

  if(ERF_ENABLE_WARM_NO_PRECIP)
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_WARM_NO_PRECIP)
  endif()

  if(ERF_ENABLE_PARTICLES)
    target_sources(${erf_lib_name} PRIVATE
                   ${SRC_DIR}/Particles/ERFPCEvolve.cpp
                   ${SRC_DIR}/Particles/ERFPCInitializations.cpp
                   ${SRC_DIR}/Particles/ERFPCUtils.cpp
                   ${SRC_DIR}/Particles/ERFTracers.cpp)
    target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/Particles>)
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_PARTICLES)
  endif()

  if(ERF_ENABLE_FFT)
    target_sources(${erf_lib_name} PRIVATE
                   ${SRC_DIR}/LinearSolvers/ERF_SolveWithFFT.cpp)
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_FFT)
  endif()

  if(ERF_ENABLE_NETCDF)
    target_sources(${erf_lib_name} PRIVATE
                   ${SRC_DIR}/IO/ERF_NCInterface.cpp
                   ${SRC_DIR}/IO/ERF_NCPlotFile.cpp
                   ${SRC_DIR}/IO/ERF_NCCheckpoint.cpp
                   ${SRC_DIR}/IO/ERF_NCMultiFabFile.cpp
                   ${SRC_DIR}/IO/ERF_ReadFromMetgrid.cpp
                   ${SRC_DIR}/IO/ERF_ReadFromWRFBdy.cpp
                   ${SRC_DIR}/IO/ERF_ReadFromWRFInput.cpp
                   ${SRC_DIR}/IO/ERF_NCColumnFile.cpp)
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_NETCDF)
  endif()

  if(ERF_ENABLE_NOAH)
    target_include_directories(${erf_lib_name} PUBLIC
                               $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/LandSurfaceModel/NOAH>
                               $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Submodules/NOAH-MP/drivers/hrldas>)
    target_sources(${erf_lib_name} PRIVATE
                   ${SRC_DIR}/LandSurfaceModel/NOAH/ERF_NOAH.cpp)
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_NOAH)
    target_link_libraries_system(${erf_lib_name} PUBLIC NoahMP::noahmp)
  endif()

  if(ERF_ENABLE_RRTMGP)
    target_sources(${erf_lib_name} PRIVATE
                   ${SRC_DIR}/Utils/ERF_Orbit.cpp
                   ${SRC_DIR}/Radiation/ERF_InitRRTMGP.cpp
                   ${SRC_DIR}/Radiation/ERF_FinalizeRRTMGP.cpp
                   ${SRC_DIR}/Radiation/ERF_RunLongWaveRRTMGP.cpp
                   ${SRC_DIR}/Radiation/ERF_RunShortWaveRRTMGP.cpp
                   ${SRC_DIR}/Radiation/ERF_CloudRadProps.cpp
                   ${SRC_DIR}/Radiation/ERF_AeroRadProps.cpp
                   ${SRC_DIR}/Radiation/ERF_Optics.cpp
                   ${SRC_DIR}/Radiation/ERF_Radiation.cpp
                   ${SRC_DIR}/Radiation/ERF_Albedo.cpp
                   ${CMAKE_SOURCE_DIR}/Submodules/RRTMGP/cpp/examples/mo_load_coefficients.cpp
                   ${CMAKE_SOURCE_DIR}/Submodules/RRTMGP/cpp/extensions/fluxes_byband/mo_fluxes_byband_kernels.cpp
                  )

    # The interface code needs to know about the RRTMGP includes
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_RRTMGP)

    target_include_directories(${erf_lib_name} SYSTEM PUBLIC
                               ${CMAKE_SOURCE_DIR}/Submodules/RRTMGP/cpp/extensions/fluxes_byband
                               ${CMAKE_SOURCE_DIR}/Submodules/RRTMGP/cpp/extensions/cloud_optics
                               ${CMAKE_SOURCE_DIR}/Submodules/RRTMGP/cpp/examples
                              )
  endif()

  target_sources(${erf_lib_name}
     PRIVATE
       ${SRC_DIR}/ERF_Derive.cpp
       ${SRC_DIR}/ERF.cpp
       ${SRC_DIR}/ERF_MakeNewArrays.cpp
       ${SRC_DIR}/ERF_MakeNewLevel.cpp
       ${SRC_DIR}/ERF_ReadWaves.cpp
       ${SRC_DIR}/ERF_Tagging.cpp
       ${SRC_DIR}/Advection/ERF_AdvectionSrcForMom.cpp
       ${SRC_DIR}/Advection/ERF_AdvectionSrcForState.cpp
       ${SRC_DIR}/Advection/ERF_AdvectionSrcForOpenBC.cpp
       ${SRC_DIR}/BoundaryConditions/ERF_ABLMost.cpp
       ${SRC_DIR}/BoundaryConditions/ERF_MOSTAverage.cpp
       ${SRC_DIR}/BoundaryConditions/ERF_BoundaryConditionsCons.cpp
       ${SRC_DIR}/BoundaryConditions/ERF_BoundaryConditionsXvel.cpp
       ${SRC_DIR}/BoundaryConditions/ERF_BoundaryConditionsYvel.cpp
       ${SRC_DIR}/BoundaryConditions/ERF_BoundaryConditionsZvel.cpp
       ${SRC_DIR}/BoundaryConditions/ERF_BoundaryConditionsBaseState.cpp
       ${SRC_DIR}/BoundaryConditions/ERF_BoundaryConditionsBndryReg.cpp
       ${SRC_DIR}/BoundaryConditions/ERF_BoundaryConditionsRealbdy.cpp
       ${SRC_DIR}/BoundaryConditions/ERF_FillPatch.cpp
       ${SRC_DIR}/BoundaryConditions/ERF_FillCoarsePatch.cpp
       ${SRC_DIR}/BoundaryConditions/ERF_FillIntermediatePatch.cpp
       ${SRC_DIR}/BoundaryConditions/ERF_FillBdyCCVels.cpp
       ${SRC_DIR}/BoundaryConditions/ERF_FillPatcher.cpp
       ${SRC_DIR}/BoundaryConditions/ERF_PhysBCFunct.cpp
       ${SRC_DIR}/Diffusion/ERF_DiffusionSrcForMom_N.cpp
       ${SRC_DIR}/Diffusion/ERF_DiffusionSrcForMom_T.cpp
       ${SRC_DIR}/Diffusion/ERF_DiffusionSrcForState_N.cpp
       ${SRC_DIR}/Diffusion/ERF_DiffusionSrcForState_T.cpp
       ${SRC_DIR}/Diffusion/ERF_ComputeStress_N.cpp
       ${SRC_DIR}/Diffusion/ERF_ComputeStress_T.cpp
       ${SRC_DIR}/Diffusion/ERF_ComputeStrain_N.cpp
       ${SRC_DIR}/Diffusion/ERF_ComputeStrain_T.cpp
       ${SRC_DIR}/Diffusion/ERF_ComputeTurbulentViscosity.cpp
       ${SRC_DIR}/EB/ERF_InitEB.cpp
       ${SRC_DIR}/EB/ERF_EBBox.cpp
       ${SRC_DIR}/EB/ERF_EBRegular.cpp
       ${SRC_DIR}/EB/ERF_InitEB.cpp
       ${SRC_DIR}/EB/ERF_WriteEBSurface.cpp 
       ${SRC_DIR}/Initialization/ERF_InitBCs.cpp
       ${SRC_DIR}/Initialization/ERF_InitCustom.cpp
       ${SRC_DIR}/Initialization/ERF_InitFromHSE.cpp
       ${SRC_DIR}/Initialization/ERF_InitFromInputSounding.cpp
       ${SRC_DIR}/Initialization/ERF_InitFromWRFInput.cpp
       ${SRC_DIR}/Initialization/ERF_InitFromMetgrid.cpp
       ${SRC_DIR}/Initialization/ERF_InitGeowind.cpp
       ${SRC_DIR}/Initialization/ERF_InitRayleigh.cpp
       ${SRC_DIR}/Initialization/ERF_InitSponge.cpp
       ${SRC_DIR}/Initialization/ERF_InitUniform.cpp
       ${SRC_DIR}/Initialization/ERF_Init1D.cpp
       ${SRC_DIR}/Initialization/ERF_InitTurbPert.cpp
       ${SRC_DIR}/IO/ERF_Checkpoint.cpp
       ${SRC_DIR}/IO/ERF_ReadBndryPlanes.cpp
       ${SRC_DIR}/IO/ERF_WriteBndryPlanes.cpp
       ${SRC_DIR}/IO/ERF_Write1DProfiles.cpp
       ${SRC_DIR}/IO/ERF_Write1DProfiles_stag.cpp
       ${SRC_DIR}/IO/ERF_WriteScalarProfiles.cpp
       ${SRC_DIR}/IO/ERF_Plotfile.cpp
       ${SRC_DIR}/IO/ERF_WriteJobInfo.cpp
       ${SRC_DIR}/IO/ERF_ConsoleIO.cpp
       ${SRC_DIR}/PBL/ERF_ComputeDiffusivityMYNN25.cpp
       ${SRC_DIR}/PBL/ERF_ComputeDiffusivityYSU.cpp
       ${SRC_DIR}/SourceTerms/ERF_ApplySpongeZoneBCs.cpp
       ${SRC_DIR}/SourceTerms/ERF_ApplySpongeZoneBCs_ReadFromFile.cpp
       ${SRC_DIR}/SourceTerms/ERF_MakeBuoyancy.cpp
       ${SRC_DIR}/SourceTerms/ERF_AddThinBodySources.cpp
       ${SRC_DIR}/SourceTerms/ERF_MakeMomSources.cpp
       ${SRC_DIR}/SourceTerms/ERF_MakeSources.cpp
       ${SRC_DIR}/SourceTerms/ERF_MoistSetRhs.cpp
       ${SRC_DIR}/SourceTerms/ERF_NumericalDiffusion.cpp
       ${SRC_DIR}/SourceTerms/ERF_ForestDrag.cpp
       ${SRC_DIR}/SourceTerms/ERF_TerrainDrag.cpp       
       ${SRC_DIR}/TimeIntegration/ERF_ComputeTimestep.cpp
       ${SRC_DIR}/TimeIntegration/ERF_Advance.cpp
       ${SRC_DIR}/TimeIntegration/ERF_TimeStep.cpp
       ${SRC_DIR}/TimeIntegration/ERF_AdvanceDycore.cpp
       ${SRC_DIR}/TimeIntegration/ERF_AdvanceMicrophysics.cpp
       ${SRC_DIR}/TimeIntegration/ERF_AdvanceLSM.cpp
       ${SRC_DIR}/TimeIntegration/ERF_AdvanceRadiation.cpp
       ${SRC_DIR}/TimeIntegration/ERF_MakeFastCoeffs.cpp
       ${SRC_DIR}/TimeIntegration/ERF_MakeTauTerms.cpp
       ${SRC_DIR}/TimeIntegration/ERF_SlowRhsPre.cpp
       ${SRC_DIR}/TimeIntegration/ERF_SlowRhsPost.cpp
       ${SRC_DIR}/TimeIntegration/ERF_FastRhs_N.cpp
       ${SRC_DIR}/TimeIntegration/ERF_FastRhs_T.cpp
       ${SRC_DIR}/TimeIntegration/ERF_FastRhs_MT.cpp
       ${SRC_DIR}/Utils/ERF_AverageDown.cpp
       ${SRC_DIR}/Utils/ERF_ChopGrids.cpp
       ${SRC_DIR}/Utils/ERF_MomentumToVelocity.cpp
       ${SRC_DIR}/LinearSolvers/ERF_PoissonSolve.cpp
       ${SRC_DIR}/LinearSolvers/ERF_PoissonSolve_tb.cpp
       ${SRC_DIR}/LinearSolvers/ERF_PoissonWallDist.cpp
       ${SRC_DIR}/LinearSolvers/ERF_ComputeDivergence.cpp 
       ${SRC_DIR}/LinearSolvers/ERF_SolveWithEBMLMG.cpp
       ${SRC_DIR}/LinearSolvers/ERF_SolveWithGMRES.cpp
       ${SRC_DIR}/LinearSolvers/ERF_SolveWithMLMG.cpp
       ${SRC_DIR}/LinearSolvers/ERF_TerrainPoisson.cpp
       ${SRC_DIR}/Utils/ERF_TerrainMetrics.cpp
       ${SRC_DIR}/Utils/ERF_VelocityToMomentum.cpp
       ${SRC_DIR}/Utils/ERF_InteriorGhostCells.cpp
       ${SRC_DIR}/Utils/ERF_TimeAvgVel.cpp
       ${SRC_DIR}/Microphysics/SAM/ERF_InitSAM.cpp
       ${SRC_DIR}/Microphysics/SAM/ERF_CloudSAM.cpp
       ${SRC_DIR}/Microphysics/SAM/ERF_IceFall.cpp
       ${SRC_DIR}/Microphysics/SAM/ERF_Precip.cpp
       ${SRC_DIR}/Microphysics/SAM/ERF_PrecipFall.cpp
       ${SRC_DIR}/Microphysics/SAM/ERF_UpdateSAM.cpp
       ${SRC_DIR}/Microphysics/Kessler/ERF_InitKessler.cpp
       ${SRC_DIR}/Microphysics/Kessler/ERF_Kessler.cpp
       ${SRC_DIR}/Microphysics/Kessler/ERF_UpdateKessler.cpp
       ${SRC_DIR}/Microphysics/SatAdj/ERF_InitSatAdj.cpp
       ${SRC_DIR}/Microphysics/SatAdj/ERF_SatAdj.cpp
       ${SRC_DIR}/Microphysics/SatAdj/ERF_UpdateSatAdj.cpp
       ${SRC_DIR}/WindFarmParametrization/Fitch/ERF_AdvanceFitch.cpp
       ${SRC_DIR}/WindFarmParametrization/EWP/ERF_AdvanceEWP.cpp
       ${SRC_DIR}/WindFarmParametrization/SimpleActuatorDisk/ERF_AdvanceSimpleAD.cpp
       ${SRC_DIR}/WindFarmParametrization/GeneralActuatorDisk/ERF_AdvanceGeneralAD.cpp
       ${SRC_DIR}/LandSurfaceModel/SLM/ERF_SLM.cpp
       ${SRC_DIR}/LandSurfaceModel/MM5/ERF_MM5.cpp
  )

  include(AMReXBuildInfo)
  generate_buildinfo(${erf_lib_name} ${CMAKE_SOURCE_DIR})
if (${ERF_USE_INTERNAL_AMREX})
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${AMREX_SUBMOD_LOCATION}/Tools/C_scripts>)
endif()

  if(ERF_ENABLE_NETCDF)
    if(NETCDF_FOUND)
      #Link our executable to the NETCDF libraries, etc
      target_link_libraries(${erf_lib_name} PUBLIC ${NETCDF_LINK_LIBRARIES})
      target_include_directories(${erf_lib_name} PUBLIC ${NETCDF_INCLUDE_DIRS})
    endif()
  endif()

  if(ERF_ENABLE_RRTMGP)
    target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/Radiation>)
  endif()

  if(ERF_ENABLE_MPI)
    target_link_libraries(${erf_lib_name} PUBLIC $<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_CXX>)
  endif()

  #ERF include directories
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/Advection>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/BoundaryConditions>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/DataStructs>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/Diffusion>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/EB>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/Initialization>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/IO>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/LinearSolvers>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/PBL>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/SourceTerms>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/TimeIntegration>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/Utils>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/Microphysics>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/Microphysics/Null>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/Microphysics/SAM>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/Microphysics/Kessler>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/Microphysics/SatAdj>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/WindFarmParametrization>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/WindFarmParametrization/Null>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/WindFarmParametrization/Fitch>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/WindFarmParametrization/EWP>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/WindFarmParametrization/SimpleActuatorDisk>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/WindFarmParametrization/GeneralActuatorDisk>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/LandSurfaceModel>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/LandSurfaceModel/Null>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/LandSurfaceModel/SLM>)
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/LandSurfaceModel/MM5>)

  if(ERF_ENABLE_RRTMGP)
     target_link_libraries(${erf_lib_name} PUBLIC yakl)
     target_link_libraries(${erf_lib_name} PUBLIC rrtmgp)
  endif()

  #Link to amrex library
  target_link_libraries_system(${erf_lib_name} PUBLIC AMReX::amrex)
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

  if(NOT "${erf_exe_name}" STREQUAL "erf_unit_tests")
  target_sources(${erf_exe_name}
     PRIVATE
       ${SRC_DIR}/main.cpp
  )
  endif()

  target_link_libraries(${erf_exe_name}  PUBLIC ${erf_lib_name})
  include(${CMAKE_SOURCE_DIR}/CMake/SetERFCompileFlags.cmake)
  set_erf_compile_flags(${erf_exe_name})

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
