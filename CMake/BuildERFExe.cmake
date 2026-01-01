function(target_link_libraries_system target visibility)
  set(libs ${ARGN})
  foreach(lib ${libs})
    get_target_property(lib_include_dirs ${lib} INTERFACE_INCLUDE_DIRECTORIES)
    target_include_directories(${target} SYSTEM ${visibility} ${lib_include_dirs})
    target_link_libraries(${target} ${visibility} ${lib})
  endforeach(lib)
endfunction(target_link_libraries_system)

# Link library but only propagate include directories, not link options
function(target_link_libraries_includes_only target visibility lib)
  # Link the library (this target will use it)
  target_link_libraries(${target} PRIVATE ${lib})
  
  # But propagate includes with specified visibility
  get_target_property(lib_include_dirs ${lib} INTERFACE_INCLUDE_DIRECTORIES)
  if(lib_include_dirs)
    target_include_directories(${target} ${visibility} ${lib_include_dirs})
  endif()
endfunction()

function(build_erf_lib erf_lib_name)

  set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)
  set(BIN_DIR ${CMAKE_BINARY_DIR}/Source/${erf_lib_name})

  include(${CMAKE_SOURCE_DIR}/CMake/SetERFCompileFlags.cmake)
  set_erf_compile_flags(${erf_lib_name})

  target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_MOISTURE)

  # NOTE: EKAT provides KOKKOS
  if(ERF_ENABLE_EKAT)
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_KOKKOS)
    target_include_directories(${erf_lib_name} PUBLIC
                               $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/PhysicsInterfaces>
                               $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Submodules/ekat/src/pack>
                               $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Submodules/ekat/src/algorithm>
                              )
  endif()

  if(ERF_ENABLE_IMPLICIT_W)
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_IMPLICIT_W)
  endif()

  if(ERF_ENABLE_MULTIBLOCK)
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_MULTIBLOCK)
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

  ########################## NETCDF ##################################
  if(ERF_ENABLE_NETCDF)
    target_sources(${erf_lib_name} PRIVATE
                   ${SRC_DIR}/Initialization/ERF_InitFromWRFInput.cpp
                   ${SRC_DIR}/Initialization/ERF_InitFromMetgrid.cpp
                   ${SRC_DIR}/Initialization/ERF_InitFromNCFile.cpp
                   ${SRC_DIR}/IO/ERF_NCInterface.cpp
                   ${SRC_DIR}/IO/ERF_NCPlotFile.cpp
                   ${SRC_DIR}/IO/ERF_ReadFromMetgrid.cpp
                   ${SRC_DIR}/IO/ERF_ReadFromWRFBdy.cpp
                   ${SRC_DIR}/IO/ERF_ReadFromWRFInput.cpp
                   ${SRC_DIR}/IO/ERF_ReadFromWRFLow.cpp
                   ${SRC_DIR}/IO/ERF_NCColumnFile.cpp)
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_NETCDF)
  endif()

  ########################## NOAHMP ##################################
  if(ERF_ENABLE_NOAHMP)
    target_include_directories(${erf_lib_name} PUBLIC
                               $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/LandSurfaceModel/Noah-MP>
                               $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Submodules/Noah-MP/drivers/erf>)
    target_sources(${erf_lib_name} PRIVATE
                   ${SRC_DIR}/LandSurfaceModel/Noah-MP/ERF_NOAHMP.cpp)
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_NOAHMP)
    target_link_libraries_system(${erf_lib_name} PUBLIC NoahMP::noahmp)
  endif()

  ########################### GPU defs for KOKKOS #################################
  if(ERF_ENABLE_CUDA)
      target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_CUDA)
  endif()
  if(ERF_ENABLE_HIP)
      target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_HIP)
  endif()
  if(ERF_ENABLE_SYCL)
      target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_SYCL)
  endif()

  ########################### RRTMGP #################################
  if(ERF_ENABLE_RRTMGP)
    target_include_directories(${erf_lib_name} PUBLIC
                               $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/PhysicsInterfaces/Radiation>
                               $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Submodules/RRTMGP/cpp>
                               $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Submodules/RRTMGP/cpp/rrtmgp>
                               $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Submodules/RRTMGP/cpp/rrtmgp/kernels>
                               $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Submodules/RRTMGP/cpp/rte>
                               $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Submodules/RRTMGP/cpp/rte/kernels>
                               $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Submodules/RRTMGP/cpp/examples>
                               $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Submodules/RRTMGP/cpp/examples/all-sky>
                               $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Submodules/RRTMGP/cpp/extensions/cloud_optics>
                               $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Submodules/RRTMGP/cpp/extensions/fluxes_byband>
                              )
    target_sources(${erf_lib_name} PRIVATE
                   ${SRC_DIR}/PhysicsInterfaces/Radiation/ERF_RRTMGP_Interface.cpp
                   ${SRC_DIR}/PhysicsInterfaces/Radiation/ERF_Radiation.cpp
                   ${CMAKE_SOURCE_DIR}/Submodules/RRTMGP/cpp/rrtmgp/mo_rrtmgp_util_reorder.cpp
                   ${CMAKE_SOURCE_DIR}/Submodules/RRTMGP/cpp/rrtmgp/kernels/mo_gas_optics_kernels.cpp
                   ${CMAKE_SOURCE_DIR}/Submodules/RRTMGP/cpp/rte/expand_and_transpose.cpp
                   ${CMAKE_SOURCE_DIR}/Submodules/RRTMGP/cpp/rte/kernels/mo_fluxes_broadband_kernels.cpp
                   ${CMAKE_SOURCE_DIR}/Submodules/RRTMGP/cpp/rte/kernels/mo_optical_props_kernels.cpp
                   ${CMAKE_SOURCE_DIR}/Submodules/RRTMGP/cpp/rte/kernels/mo_rte_solver_kernels.cpp
                   ${CMAKE_SOURCE_DIR}/Submodules/RRTMGP/cpp/examples/mo_load_coefficients.cpp
                   ${CMAKE_SOURCE_DIR}/Submodules/RRTMGP/cpp/examples/all-sky/mo_garand_atmos_io.cpp
                   ${CMAKE_SOURCE_DIR}/Submodules/RRTMGP/cpp/examples/all-sky/mo_load_cloud_coefficients.cpp
                   ${CMAKE_SOURCE_DIR}/Submodules/RRTMGP/cpp/extensions/fluxes_byband/mo_fluxes_byband_kernels.cpp
                  )
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_RRTMGP)
    target_compile_definitions(${erf_lib_name} PUBLIC RRTMGP_ENABLE_KOKKOS)
  endif()

  ########################### SHOC #################################
  if(ERF_ENABLE_SHOC)
    target_include_directories(${erf_lib_name} PUBLIC
                               $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/PhysicsInterfaces/Shoc>
                               $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src>
                               $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics>
                               $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/share>
                               $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc>
                               $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti>
                               $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/impl>
                              )
    target_sources(${erf_lib_name} PRIVATE
                   ${SRC_DIR}/PhysicsInterfaces/Shoc/ERF_ShocInterface.cpp
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/share/physics_saturation.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/disp/shoc_assumed_pdf_disp.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/disp/shoc_check_tke_disp.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/disp/shoc_compute_shoc_temperature_disp.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/disp/shoc_compute_shoc_vapor_disp.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/disp/shoc_diag_obklen_disp.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/disp/shoc_diag_second_shoc_moments_disp.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/disp/shoc_diag_third_shoc_moments_disp.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/disp/shoc_energy_fixer_disp.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/disp/shoc_energy_integrals_disp.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/disp/shoc_grid_disp.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/disp/shoc_length_disp.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/disp/shoc_pblintd_disp.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/disp/shoc_tke_disp.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/disp/shoc_update_host_dse_disp.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/disp/shoc_update_prognostics_implicit_disp.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_adv_sgs_tke.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_assumed_pdf.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_calc_shoc_varorcovar.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_calc_shoc_vertflux.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_check_length_scale_shoc_length.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_check_tke.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_clipping_diag_third_shoc_moments.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_compute_brunt_shoc_length.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_compute_diag_third_shoc_moment.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_compute_l_inf_shoc_length.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_compute_shoc_mix_shoc_length.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_compute_shoc_temperature.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_compute_shoc_vapor.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_compute_shr_prod.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_compute_tmpi.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_diag_obklen.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_diag_second_moments.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_diag_second_moments_lbycond.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_diag_second_moments_srf.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_diag_second_moments_ubycond.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_diag_second_shoc_moments.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_diag_third_shoc_moments.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_dp_inverse.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_eddy_diffusivities.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_energy_fixer.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_energy_integrals.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_grid.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_integ_column_stability.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_isotropic_ts.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_linear_interp.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_length.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_main.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_pblintd.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_pblintd_check_pblh.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_pblintd_cldcheck.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_pblintd_height.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_pblintd_init_pot.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_pblintd_surf_temp.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_tridiag_solver.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_tke.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_update_host_dse.cpp>
                   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/E3SM/components/eamxx/src/physics/shoc/eti/shoc_update_prognostics_implicit.cpp>
                  )
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_SHOC)
    target_compile_definitions(${erf_lib_name} PUBLIC SCREAM_SHOC_SMALL_KERNELS)              
  endif()

  if(ERF_ENABLE_MORR_FORT)
  target_sources(${erf_lib_name}
     PRIVATE
       ${SRC_DIR}/Microphysics/Morrison/ERF_module_mp_morr_two_moment.F90
       ${SRC_DIR}/Microphysics/Morrison/ERF_module_mp_morr_two_moment_isohelper.F90
       ${SRC_DIR}/Microphysics/Morrison/ERF_module_model_constants.F90
       )
  target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_MORR_FORT)
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
       ${SRC_DIR}/Advection/ERF_AdvectionSrcForMom_ConstantDz.cpp
       ${SRC_DIR}/Advection/ERF_AdvectionSrcForMom_StretchedDz.cpp
       ${SRC_DIR}/Advection/ERF_AdvectionSrcForMom_EB.cpp
       ${SRC_DIR}/Advection/ERF_AdvectionSrcForMom_TF.cpp
       ${SRC_DIR}/Advection/ERF_AdvectionSrcForState.cpp
       ${SRC_DIR}/Advection/ERF_AdvectionSrcForOpenBC.cpp
       ${SRC_DIR}/BoundaryConditions/ERF_SurfaceLayer.cpp
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
       ${SRC_DIR}/Diffusion/ERF_DiffusionSrcForMom.cpp
       ${SRC_DIR}/Diffusion/ERF_DiffusionSrcForMom_EB.cpp
       ${SRC_DIR}/Diffusion/ERF_DiffusionSrcForState_N.cpp
       ${SRC_DIR}/Diffusion/ERF_DiffusionSrcForState_S.cpp
       ${SRC_DIR}/Diffusion/ERF_DiffusionSrcForState_T.cpp
       ${SRC_DIR}/Diffusion/ERF_ImplicitDiff_N.cpp
       ${SRC_DIR}/Diffusion/ERF_ImplicitDiff_S.cpp
       ${SRC_DIR}/Diffusion/ERF_ImplicitDiff_T.cpp
       ${SRC_DIR}/Diffusion/ERF_ComputeStress_EB.cpp
       ${SRC_DIR}/Diffusion/ERF_ComputeStress_N.cpp
       ${SRC_DIR}/Diffusion/ERF_ComputeStress_S.cpp
       ${SRC_DIR}/Diffusion/ERF_ComputeStress_T.cpp
       ${SRC_DIR}/Diffusion/ERF_ComputeStrain_EB.cpp
       ${SRC_DIR}/Diffusion/ERF_ComputeStrain_N.cpp
       ${SRC_DIR}/Diffusion/ERF_ComputeStrain_S.cpp
       ${SRC_DIR}/Diffusion/ERF_ComputeStrain_T.cpp
       ${SRC_DIR}/Diffusion/ERF_ComputeTurbulentViscosity.cpp
       ${SRC_DIR}/EB/ERF_EBAdvectionSrcForState.cpp
       ${SRC_DIR}/EB/ERF_EBAux.cpp
       ${SRC_DIR}/EB/ERF_EB.cpp
       ${SRC_DIR}/EB/ERF_EBCutCell.cpp
       ${SRC_DIR}/EB/ERF_EBRedistribute.cpp
       ${SRC_DIR}/Initialization/ERF_InitBCs.cpp
       ${SRC_DIR}/Initialization/ERF_InitCustom.cpp
       ${SRC_DIR}/Initialization/ERF_InitFromHSE.cpp
       ${SRC_DIR}/Initialization/ERF_InitFromInputSounding.cpp
       ${SRC_DIR}/Initialization/ERF_InitGeowind.cpp
       ${SRC_DIR}/Initialization/ERF_InitRayleigh.cpp
       ${SRC_DIR}/Initialization/ERF_InitSponge.cpp
       ${SRC_DIR}/Initialization/ERF_InitUniform.cpp
       ${SRC_DIR}/Initialization/ERF_Init1D.cpp
       ${SRC_DIR}/Initialization/ERF_InitTurbPert.cpp
       ${SRC_DIR}/Initialization/ERF_InitImmersedForcing.cpp
       ${SRC_DIR}/IO/ERF_Checkpoint.cpp
       ${SRC_DIR}/IO/ERF_ReadBndryPlanes.cpp
       ${SRC_DIR}/IO/ERF_WriteBndryPlanes.cpp
       ${SRC_DIR}/IO/ERF_TrackerOutput.cpp
       ${SRC_DIR}/IO/ERF_Write1DProfiles.cpp
       ${SRC_DIR}/IO/ERF_Write1DProfiles_stag.cpp
       ${SRC_DIR}/IO/ERF_WriteScalarProfiles.cpp
       ${SRC_DIR}/IO/ERF_Plotfile.cpp
       ${SRC_DIR}/IO/ERF_WriteSubvolume.cpp
       ${SRC_DIR}/IO/ERF_WriteJobInfo.cpp
       ${SRC_DIR}/IO/ERF_ConsoleIO.cpp
       ${SRC_DIR}/LinearSolvers/ERF_PoissonSolve.cpp
       ${SRC_DIR}/LinearSolvers/ERF_PoissonSolve_tb.cpp
       ${SRC_DIR}/LinearSolvers/ERF_PoissonWallDist.cpp
       ${SRC_DIR}/LinearSolvers/ERF_ComputeDivergence.cpp
       ${SRC_DIR}/LinearSolvers/ERF_FillZeroAreaFaceFluxes.cpp
       ${SRC_DIR}/LinearSolvers/ERF_ImposeBCsOnPhi.cpp
       ${SRC_DIR}/LinearSolvers/ERF_SolveWithEBMLMG.cpp
       ${SRC_DIR}/LinearSolvers/ERF_SolveWithGMRES.cpp
       ${SRC_DIR}/LinearSolvers/ERF_SolveWithMLMG.cpp
       ${SRC_DIR}/LinearSolvers/ERF_TerrainPoisson.cpp
       ${SRC_DIR}/Microphysics/Morrison/ERF_InitMorrison.cpp
       ${SRC_DIR}/Microphysics/Morrison/ERF_AdvanceMorrison.cpp
       ${SRC_DIR}/Microphysics/Morrison/ERF_UpdateMorrison.cpp
       ${SRC_DIR}/Microphysics/Morrison/ERF_Morrison_Plot.cpp
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
       ${SRC_DIR}/PBL/ERF_ComputeDiffusivityMYJ.cpp
       ${SRC_DIR}/PBL/ERF_ComputeDiffusivityMYNN25.cpp
       ${SRC_DIR}/PBL/ERF_ComputeDiffusivityMYNNEDMF.cpp
       ${SRC_DIR}/PBL/ERF_ComputeDiffusivityYSU.cpp
       ${SRC_DIR}/PBL/ERF_ComputeDiffusivityMRF.cpp
       ${SRC_DIR}/SourceTerms/ERF_ApplySpongeZoneBCs.cpp
       ${SRC_DIR}/SourceTerms/ERF_ApplySpongeZoneBCs_ReadFromFile.cpp
       ${SRC_DIR}/SourceTerms/ERF_ApplyBndryForcing_Forecast.cpp
       ${SRC_DIR}/SourceTerms/ERF_MakeBuoyancy.cpp
       ${SRC_DIR}/SourceTerms/ERF_AddThinBodySources.cpp
       ${SRC_DIR}/SourceTerms/ERF_MakeGradP.cpp
       ${SRC_DIR}/SourceTerms/ERF_MakeMomSources.cpp
       ${SRC_DIR}/SourceTerms/ERF_MakeSources.cpp
       ${SRC_DIR}/SourceTerms/ERF_MoistSetRhs.cpp
       ${SRC_DIR}/SourceTerms/ERF_NumericalDiffusion.cpp
       ${SRC_DIR}/SourceTerms/ERF_ForestDrag.cpp
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
       ${SRC_DIR}/TimeIntegration/ERF_Substep_NS.cpp
       ${SRC_DIR}/TimeIntegration/ERF_Substep_T.cpp
       ${SRC_DIR}/TimeIntegration/ERF_Substep_MT.cpp
       ${SRC_DIR}/Utils/ERF_AverageDown.cpp
       ${SRC_DIR}/Utils/ERF_ChopGrids.cpp
       ${SRC_DIR}/Utils/ERF_ConvertForProjection.cpp
       ${SRC_DIR}/Utils/ERF_InitZLevels.cpp
       ${SRC_DIR}/Utils/ERF_MakeSubdomains.cpp
       ${SRC_DIR}/Utils/ERF_MomentumToVelocity.cpp
       ${SRC_DIR}/Utils/ERF_TerrainMetrics.cpp
       ${SRC_DIR}/Utils/ERF_VelocityToMomentum.cpp
       ${SRC_DIR}/Utils/ERF_InteriorGhostCells.cpp
       ${SRC_DIR}/Utils/ERF_ThinBodyWallDist.cpp
       ${SRC_DIR}/Utils/ERF_TimeAvgVel.cpp
       ${SRC_DIR}/Utils/ERF_VolWgtSum.cpp
       ${SRC_DIR}/Utils/ERF_WeatherDataInterpolation.cpp
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
      target_link_libraries(${erf_lib_name} PUBLIC ${NETCDF_LINK_LIBRARIES})
      target_include_directories(${erf_lib_name} PUBLIC ${NETCDF_INCLUDE_DIRS})
    endif()
  endif()

  if(ERF_ENABLE_EKAT)
    target_link_libraries(${erf_lib_name} PUBLIC ekat::AllLibs spdlog::spdlog Kokkos::kokkos)
  endif()
#  if(ERF_ENABLE_EKAT)
#    # Link privately to avoid duplicate link flags, but propagate includes
#    target_link_libraries(${erf_lib_name} PRIVATE ekat::AllLibs spdlog::spdlog Kokkos::kokkos)

#    # Manually propagate include directories (but not link options)
#    foreach(lib IN ITEMS ekat::AllLibs spdlog::spdlog Kokkos::kokkos)
#        if(TARGET ${lib})
#            get_target_property(inc_dirs ${lib} INTERFACE_INCLUDE_DIRECTORIES)
#            if(inc_dirs)
#                target_include_directories(${erf_lib_name} PUBLIC ${inc_dirs})
#            endif()
#        endif()
#    endforeach()
#  endif()

  if(ERF_ENABLE_MPI)
    target_link_libraries(${erf_lib_name} PUBLIC $<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_CXX>)
    if(ERF_ENABLE_MORR_FORT OR ERF_ENABLE_NOAHMP)
      target_link_libraries(${erf_lib_name} PUBLIC $<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_Fortran>)
    endif()
  endif()

  # Workaround for gcc-8 where std::filesystem is in libstdc++fs. Starting with
  # gcc-9 std::filesystem is part of libstdc++.
  target_link_libraries(${erf_lib_name} PUBLIC $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:stdc++fs>)

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
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/Microphysics/Morrison>)
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
  target_include_directories(${erf_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/PhysicsInterfaces/Radiation/>)

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

target_link_libraries(${erf_exe_name} PUBLIC ${erf_lib_name})
include(${CMAKE_SOURCE_DIR}/CMake/SetERFCompileFlags.cmake)
set_erf_compile_flags(${erf_exe_name})

  if(ERF_ENABLE_EKAT)
    # Link privately to avoid duplicate link flags, but propagate includes
    target_link_libraries(${erf_exe_name} PRIVATE ekat::AllLibs spdlog::spdlog Kokkos::kokkos)

    # Manually propagate include directories (but not link options)
    foreach(lib IN ITEMS ekat::AllLibs spdlog::spdlog Kokkos::kokkos)
        if(TARGET ${lib})
            get_target_property(inc_dirs ${lib} INTERFACE_INCLUDE_DIRECTORIES)
            if(inc_dirs)
                target_include_directories(${erf_exe_name} PUBLIC ${inc_dirs})
            endif()
        endif()
    endforeach()
  endif()

  # Remove duplicate HIP link flags from executable
  if(ERF_ENABLE_HIP AND ERF_ENABLE_EKAT)
      get_target_property(exe_link_opts ${erf_exe_name} LINK_OPTIONS)
      if(exe_link_opts)
          # Deduplicate by converting to list and back
          list(REMOVE_DUPLICATES exe_link_opts)

          # Keep only one --hip-link and one --offload-arch
          set(filtered "")
          set(seen_hip FALSE)
          set(seen_offload FALSE)
          foreach(opt IN LISTS exe_link_opts)
              if(opt STREQUAL "--hip-link")
                  if(NOT seen_hip)
                      list(APPEND filtered ${opt})
                      set(seen_hip TRUE)
                  endif()
              elseif(opt MATCHES "^--offload-arch=")
                  if(NOT seen_offload)
                      list(APPEND filtered ${opt})
                      set(seen_offload TRUE)
                  endif()
              else()
                  list(APPEND filtered ${opt})
              endif()
          endforeach()
          set_target_properties(${erf_exe_name} PROPERTIES LINK_OPTIONS "${filtered}")
       endif()
  endif()

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
