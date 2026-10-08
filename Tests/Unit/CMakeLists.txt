set(BUILD_GMOCK OFF CACHE BOOL
    "Build Google Mock targets for the vendored GoogleTest dependency")
set(INSTALL_GTEST OFF CACHE BOOL
    "Install targets from the vendored GoogleTest dependency")
if(MSVC)
  set(gtest_force_shared_crt ON CACHE BOOL
      "Use the shared CRT when building the vendored GoogleTest dependency")
endif()

add_subdirectory(${PROJECT_SOURCE_DIR}/Submodules/googletest
                 ${PROJECT_BINARY_DIR}/Submodules/googletest
                 EXCLUDE_FROM_ALL)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|IntelLLVM")
  foreach(_erf_gtest_target gtest gtest_main)
    if(TARGET ${_erf_gtest_target})
      target_compile_options(${_erf_gtest_target}
        PRIVATE
          -Wno-unknown-warning-option
      )
    endif()
  endforeach()
endif()

include(GoogleTest)

set(ERF_GTEST_REPORT_ROOT
    "${PROJECT_BINARY_DIR}/test-results/gtest")

function(erf_register_gtest_suite target)
  set(options)
  set(oneValueArgs
      DISCOVERY_TIMEOUT
      EXECUTION_TIMEOUT
      SKIP_RETURN_CODE
      TEST_PREFIX
      TEST_SUFFIX
      WORKING_DIRECTORY
      XML_OUTPUT_DIR)
  set(multiValueArgs LABELS)
  cmake_parse_arguments(ERF_GTEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(NOT TARGET ${target})
    message(FATAL_ERROR "Cannot register unknown gtest suite target: ${target}")
  endif()

  if(NOT ERF_GTEST_WORKING_DIRECTORY)
    set(ERF_GTEST_WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  endif()

  if(NOT ERF_GTEST_LABELS)
    set(ERF_GTEST_LABELS unit)
  endif()

  if(NOT ERF_GTEST_DISCOVERY_TIMEOUT)
    set(ERF_GTEST_DISCOVERY_TIMEOUT 60)
  endif()

  if(NOT DEFINED ERF_GTEST_EXECUTION_TIMEOUT)
    set(ERF_GTEST_EXECUTION_TIMEOUT 600)
  endif()

  if(NOT DEFINED ERF_GTEST_SKIP_RETURN_CODE)
    set(ERF_GTEST_SKIP_RETURN_CODE 77)
  endif()

  if(NOT ERF_GTEST_XML_OUTPUT_DIR)
    set(ERF_GTEST_XML_OUTPUT_DIR
        "${ERF_GTEST_REPORT_ROOT}/${target}")
  endif()

  file(MAKE_DIRECTORY "${ERF_GTEST_XML_OUTPUT_DIR}")

  list(JOIN ERF_GTEST_LABELS ";" _erf_gtest_labels)
  # GoogleTest's POST_BUILD discovery forwards TEST_PROPERTIES through a
  # Windows command line. Keep a multi-label value as one CMake argument so
  # the semicolon is not mistaken for an argument separator.
  string(REPLACE ";" "\\;" _erf_gtest_labels_property
         "${_erf_gtest_labels}")

  set(_erf_gtest_discover_extra_args)
  if(DEFINED ERF_GTEST_TEST_PREFIX AND NOT "${ERF_GTEST_TEST_PREFIX}" STREQUAL "")
    list(APPEND _erf_gtest_discover_extra_args
      TEST_PREFIX "${ERF_GTEST_TEST_PREFIX}")
  endif()
  if(DEFINED ERF_GTEST_TEST_SUFFIX AND NOT "${ERF_GTEST_TEST_SUFFIX}" STREQUAL "")
    list(APPEND _erf_gtest_discover_extra_args
      TEST_SUFFIX "${ERF_GTEST_TEST_SUFFIX}")
  endif()

  gtest_discover_tests(${target}
    DISCOVERY_MODE POST_BUILD
    WORKING_DIRECTORY ${ERF_GTEST_WORKING_DIRECTORY}
    DISCOVERY_TIMEOUT ${ERF_GTEST_DISCOVERY_TIMEOUT}
    XML_OUTPUT_DIR "${ERF_GTEST_XML_OUTPUT_DIR}"
    ${_erf_gtest_discover_extra_args}
    PROPERTIES
      LABELS "${_erf_gtest_labels_property}"
      TIMEOUT ${ERF_GTEST_EXECUTION_TIMEOUT}
      SKIP_RETURN_CODE ${ERF_GTEST_SKIP_RETURN_CODE}
  )
endfunction()
set(erf_lib_name erf_srclib)
set(erf_exe_name erf_unit_tests)
set(erf_parallel_exe_name erf_parallel_tests)

add_executable(${erf_exe_name})
target_sources(${erf_exe_name}
  PRIVATE
    ${PROJECT_SOURCE_DIR}/Exec/ERF_Prob.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ERF_GTestMain.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ERF_GTestHelloWorld.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ERF_GTestRuntimeStubs.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestSurfaceDiagnosticSource.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestSurfaceLayerStress.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestSurfaceLayer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestMOSTAverage.cpp
    # Motivation: the MOST reference height must use the true stretched cell
    # centers, and a query on a mesh face must be found rather than abort.
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestMOSTAverageZref.cpp
    # Motivation: a BTF fine level's mesh must not depend on how its grids are split in z.
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ERF_GTestTerrainFineColumns.cpp
    # Motivation: the input sponge's per-height reduction read every k from every box, so a
    # refined region that does not reach the top was read outside its data.
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ERF_GTestReduceToMaxPerHeight.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestWrfBoundaryDensity.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestWrfBoundaryMoisture.cpp
    # Motivation: face-local wall metadata must not inherit a later wall's
    # density representation or treat explicit zero as absent.
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestWallBCMetadata.cpp
    # Motivation: surface-model weights, output mapping, and model-field scaling
    # must remain correct for host and device execution paths.
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestSurfaceModel.cpp
    # Motivation: where two open lateral boundaries meet, each boundary cell must be
    # written by exactly one open-BC kernel, and by the same one however the grid is split.
    ${CMAKE_CURRENT_SOURCE_DIR}/Advection/ERF_GTestOpenBCCorners.cpp
    # Motivation: reusable scalar advection must honor explicit source, flux, and
    # tendency components instead of inferring primitive storage from native cons layout.
    ${CMAKE_CURRENT_SOURCE_DIR}/Advection/ERF_GTestScalarAdvection.cpp
    # Motivation: the z-momentum box must be shrunk only at the bottom and top of the
    # domain, so that a grid split in z leaves no w-face without a source term.
    ${CMAKE_CURRENT_SOURCE_DIR}/Advection/ERF_GTestZmomBoxDecomposition.cpp
    # Motivation: the chamber profile and perturbation must remain analytic,
    # bounded, and independent of MPI decomposition and tile shape.
    ${CMAKE_CURRENT_SOURCE_DIR}/Prob/ERF_GTestCloudChamber.cpp
    # Motivation: production scalar diffusion must use the exact low/high
    # physical-face stencil in every rotated Cartesian direction.
    ${CMAKE_CURRENT_SOURCE_DIR}/Diffusion/ERF_GTestStateDiffusionDirichlet.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Diffusion/ERF_GTestTKEBuoyancySource.cpp
    # Motivation: reusable scalar diffusion must honor explicit scalar, coefficient,
    # flux, and tendency components while preserving ERF's N/S/T mapped discretizations.
    ${CMAKE_CURRENT_SOURCE_DIR}/Diffusion/ERF_GTestScalarDiffusion.cpp
    # Motivation: the k-eqn RANS closure (AL01) must keep its stability
    # functions finite and its length-scale limits exact for any Rt.
    ${CMAKE_CURRENT_SOURCE_DIR}/Diffusion/ERF_GTestRANSClosure.cpp
    # Motivation: the opt-in WRF Smagorinsky-2D limits must reproduce WRF's smag2d_km numbers
    # (slope factor, deformation threshold, cap and their order), limit K_h before the scalar
    # diffusivities are derived from it, and the diffusive time-step rates must follow their
    # documented formulas, with the vertical part dropped only where every implicit stage is
    # shown to stay bounded.
    ${CMAKE_CURRENT_SOURCE_DIR}/Diffusion/ERF_GTestSmag2DLimiters.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Diagnostics/ERF_GTestSurfaceFluxDiagnostics.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Diagnostics/ERF_GTestSeaLevelPressure.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/ERF_GTestCanopyBiophysics.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/ERF_GTestForestDrag.cpp
    # Motivation: a calm cell above an immersed wall gave the wall law u* = 0, an Obukhov length
    # of zero and z / 0 (erf-model/ERF#4016); the stability estimate is now bounded as on flat ground.
    ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/ERF_GTestImmersedWallStability.cpp
    # Motivation: a refined level's buildings must keep level 0's numbers for the
    # materials by building, even when its tallest column is a spire level 0
    # leaves fluid, or the level holds only some of the buildings.
    ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/ERF_GTestIBSEBLevelNumbering.cpp
    # Motivation: observation nudging must read a station file exactly, relax
    # to the closed-form target with the right weight, sigma band and rate cap,
    # and put its sources on the right faces over flat, fitted and immersed
    # terrain.
    ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/ERF_GTestObsNudging.cpp
    # Motivation: the two-stream radiation sweep must keep k = 0 as the
    # surface, cool to space from the top layer, use Exner temperature,
    # and reflect (1 - eps) of the downwelling LW at a gray surface.
    ${CMAKE_CURRENT_SOURCE_DIR}/Radiation/ERF_GTestTwoStreamKernels.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Radiation/ERF_GTestTwoStreamColumn.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Radiation/ERF_GTestTwoStreamRadiationDriver.cpp
    # Motivation: with erf.ibseb.sun_mode = two_stream the building faces must
    # see the two-stream columns' sun: the same zenith, a morning sun in the
    # east, and the start date read the same way.
    ${CMAKE_CURRENT_SOURCE_DIR}/Radiation/ERF_GTestIBSEBTwoStreamSun.cpp
    # Motivation: the surface energy balance must take H and LE from a land
    # model, else from the surface layer's applied flux in W/m^2, else from
    # the scalar defaults; LE stays at its default without a moisture model.
    ${CMAKE_CURRENT_SOURCE_DIR}/Radiation/ERF_GTestSEBTurbulentFlux.cpp
    # Motivation: with the moisture coupling, the force-restore surface evaporates at
    # beta times the potential rate, beta linear between Noah-MP's wilting point and
    # field capacity for the soil category, so it dries as its water content falls.
    ${CMAKE_CURRENT_SOURCE_DIR}/Radiation/ERF_GTestSEBSoilMoisture.cpp
    # Motivation: the forcing the two-stream model hands Noah-MP must be the
    # downwelling SW and LW at the surface interface and a zenith cosine floored
    # at zero, copied into the land model's k = 0 plane and nowhere else.
    ${CMAKE_CURRENT_SOURCE_DIR}/Radiation/ERF_GTestTwoStreamLandForcing.cpp
    # Motivation: RRTMGP consumes absolute surface temperature; a
    # SurfaceLayer fallback must convert theta at the physical surface
    # pressure diagnosed from the lowest atmospheric cell,
    # respect land/water precedence, and preserve valid land LSM temperatures.
    ${CMAKE_CURRENT_SOURCE_DIR}/Radiation/ERF_GTestRRTMGP_SurfaceTemperature.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Radiation/ERF_GTestOrbitalGeometry.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestPlotfile2D.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestPlotfileSelection.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestStationSampler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestProvenance.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestCheckpoint.cpp
    # Motivation: a deck may include several files with FILE = <name>; the duplicate-input check
    # must not refuse it.
    ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestInputsName.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ERF_GTestEpochTime.cpp
    # Motivation: a station placed by latitude/longitude must land where the
    # lat/lon arrays put it on a rotated grid and on a grid across 180 degrees.
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ERF_GTestLatLonMap.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ERF_GTestGridUtils.cpp
    # Motivation: a BTF fine level's mesh must not depend on how its grids are split in z.
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ERF_GTestTerrainFineColumns.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/HashRng/ERF_GTestHashRngScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/HashRng/ERF_GTestHashRngKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/EOS/ERF_GTestEOSScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/EOS/ERF_GTestEOSKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/Interpolation/ERF_GTestInterpolationScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/Interpolation/ERF_GTestInterpolationKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/Interpolation/ERF_GTestInterpolationPolicyKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/Interpolation/ERF_GTestInterpolationWENOScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/Interpolation/ERF_GTestInterpolationWENOKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/Microphysics/ERF_GTestMicrophysicsUtilsScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/Microphysics/ERF_GTestMicrophysicsUtilsKernel.cpp
    # Motivation: supported Eulerian microphysics must consume anelastic p0
    # explicitly while preserving the compressible EOS and SAM plane contract.
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/Microphysics/ERF_GTestAnelasticMicrophysicsWiring.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SatAdj/ERF_GTestSatAdjThermo.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SatAdj/ERF_GTestSatAdjCell.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SatAdj/ERF_GTestSatAdjMultiFab.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SatAdj/ERF_GTestSatAdjBoxCoverage.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SatAdj/ERF_GTestSatAdjTemperatureDiagnostics.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SAM/ERF_GTestSAMScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SAM/ERF_GTestSAMKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SAM/ERF_GTestSAMPhysicalProperties.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/Kessler/ERF_GTestKesslerScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/Kessler/ERF_GTestKesslerKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/Kessler/ERF_GTestKesslerPhysicalProperties.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/Morrison/ERF_GTestMorrisonPhysicalProperties.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SBM/ERF_GTestSBMFoundation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SBM/ERF_GTestSBMRemapping.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/LandSurfaceModel/SLM/ERF_GTestSLMInterface.cpp
)
if(ERF_DIM STREQUAL "3")
  # Motivation: the auxiliary mapped state must reproduce ERF's computational
  # divergence, exact host stage recipes, canonical metric transfers, and ledger.
  target_sources(${erf_exe_name} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/AuxiliaryState/ERF_GTestAuxiliaryMappedTransport.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SBM/ERF_GTestSBMTransport.cpp)
endif()
if(ERF_ENABLE_NETCDF)
  target_sources(${erf_exe_name} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestTerrainNetCDF.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/ERF_GTestForestDragNetCDF.cpp
  )
  set(ERF_FOREST_TEST_FIXTURE_DIR
      "${CMAKE_CURRENT_BINARY_DIR}/forest_fixtures/BellForest")
  file(MAKE_DIRECTORY "${ERF_FOREST_TEST_FIXTURE_DIR}")
  file(COPY
    ${PROJECT_SOURCE_DIR}/Tests/test_files/BellForest/forest_lai_bell.nc
    ${PROJECT_SOURCE_DIR}/Tests/test_files/BellForest/forest_height_bell.nc
    DESTINATION "${ERF_FOREST_TEST_FIXTURE_DIR}")
  target_compile_definitions(${erf_exe_name} PRIVATE
    ERF_FOREST_TEST_FIXTURE_DIR="${ERF_FOREST_TEST_FIXTURE_DIR}")
endif()
target_sources(${erf_exe_name} PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestIntervalMeansRestart.cpp
)
# These exercise the super-droplet phase-change kernels, whose headers are
# only available when particles are enabled.
if(ERF_ENABLE_PARTICLES)
  target_sources(${erf_exe_name} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SDM/TerminalVelocity/ERF_GTestTerminalVelocityScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SDM/TerminalVelocity/ERF_GTestTerminalVelocityKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SDM/IceDeposition/ERF_GTestIceDepositionScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SDM/IceDeposition/ERF_GTestIceDepositionKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SDM/Condensation/ERF_GTestCondensationScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SDM/Condensation/ERF_GTestCondensationKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SDM/Melting/ERF_GTestMeltingScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SDM/Melting/ERF_GTestMeltingKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SDM/Freezing/ERF_GTestFreezingScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SDM/Freezing/ERF_GTestFreezingKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SDM/Coalescence/ERF_GTestCoalescenceScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SDM/Coalescence/ERF_GTestCoalescenceKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SDM/Riming/ERF_GTestRimingScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SDM/Riming/ERF_GTestRimingKernel.cpp
  )
endif()

if(ERF_ENABLE_NOAHMP)
  target_sources(${erf_exe_name} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/LandSurfaceModel/Noah-MP/ERF_GTestNoahMPResultPolicy.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/LandSurfaceModel/Noah-MP/ERF_GTestNoahMPLevelDriver.cpp
  )
endif()
build_erf_exe(${erf_exe_name} NO_ERF_MAIN)
target_link_libraries(${erf_exe_name} PRIVATE GTest::gtest)

erf_register_gtest_suite(${erf_exe_name}
  LABELS unit
  DISCOVERY_TIMEOUT 60
  EXECUTION_TIMEOUT 600
)

add_test(NAME erf_unit_tests_shuffle_repeat
  COMMAND ${erf_exe_name}
    --gtest_shuffle
    --gtest_repeat=3
    --gtest_output=xml:${ERF_GTEST_REPORT_ROOT}/${erf_exe_name}/shuffle-repeat.xml
)

# Exercise M3 through ERF's native level,
# advance, checkpoint, and exact-schema restart paths.  Header min/max values
# make the manufactured uniform spectrum and projected qc/qr independently
# observable without adding test-only diagnostics to the application.
if(TARGET erf_exec AND ERF_DIM STREQUAL "3" AND
   ERF_PRECISION STREQUAL "DOUBLE" AND NOT AMReX_HDF5)
  add_dependencies(${erf_exe_name} erf_exec)
  add_executable(erf_sbm_checkpoint_corruptor
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SBM/ERF_SBMCorrupter.cpp
    ${PROJECT_SOURCE_DIR}/Exec/ERF_Prob.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ERF_GTestRuntimeStubs.cpp)
  build_erf_exe(erf_sbm_checkpoint_corruptor NO_ERF_MAIN)

  add_test(NAME SBM_M3_ZeroCarrier_ERF_Smoke
    COMMAND ${CMAKE_COMMAND}
      -DERF_EXECUTABLE=$<TARGET_FILE:erf_exec>
      -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_zero_carrier
      -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}
      -P ${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/RunM3Smoke.cmake)
  set_tests_properties(SBM_M3_ZeroCarrier_ERF_Smoke PROPERTIES
    LABELS "integration;sbm"
    RUN_SERIAL TRUE
    TIMEOUT 120)

  add_test(NAME SBM_M3_ZeroCarrier_2M_ERF_Smoke
    COMMAND ${CMAKE_COMMAND}
      -DERF_EXECUTABLE=$<TARGET_FILE:erf_exec>
      -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_zero_carrier_2m
      -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}
      -DEXPECTED_SPECTRUM_VALUES=1.50000000000000004e-06,6.00000000000000015e-06,1.80000000000000005e-05,4.80000000000000012e-05,1.00000000000000000e+12,2.00000000000000000e+12,3.00000000000000000e+12,4.00000000000000000e+12,
      -DEXPECTED_CORE_VALUES=1.00000000000000000e+00,2.99999999999999886e+02,0.00000000000000000e+00,0.00000000000000000e+00,0.00000000000000000e+00,7.50000000000000019e-06,6.60000000000000051e-05,
      -P ${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/RunM3Smoke.cmake)
  set_tests_properties(SBM_M3_ZeroCarrier_2M_ERF_Smoke PROPERTIES
    LABELS "integration;sbm"
    RUN_SERIAL TRUE
    TIMEOUT 120)

  function(add_sbm_expected_failure_test test_name input_file diagnostic)
    add_test(NAME ${test_name}
      COMMAND ${CMAKE_COMMAND}
        -DERF_EXECUTABLE=$<TARGET_FILE:erf_exec>
        -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/${input_file}
        -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}
        "-DEXPECTED_DIAGNOSTIC=${diagnostic}"
        ${ARGN}
        -P ${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/RunExpectedFailure.cmake)
    set_tests_properties(${test_name} PROPERTIES
      LABELS "integration;sbm;negative"
      RUN_SERIAL TRUE
      TIMEOUT 120)
  endfunction()

  add_sbm_expected_failure_test(SBM_Unrealizable_2M_Init inputs_unrealizable_2m
    "SBM fixture initial state is not canonical for selected spectral representation: population=0 representation=two-moment bin=0")
  add_sbm_expected_failure_test(SBM_M3_Donor_CFL inputs_zero_carrier
    "SBM M3 donor outgoing-demand condition failed"
    "-DOVERRIDE=erf.sbm_test_donor_cfl_fault=true")
  add_test(NAME SBM_M3_Adaptive_Anelastic_Donor_CFL
    COMMAND ${CMAKE_COMMAND}
      -DERF_EXECUTABLE=$<TARGET_FILE:erf_exec>
      -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_adaptive_anelastic_donor
      -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}
      -DEXPECT_SECOND_CHECKPOINT=true
      -DSKIP_RESTART=true
      -P ${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/RunM3Smoke.cmake)
  set_tests_properties(SBM_M3_Adaptive_Anelastic_Donor_CFL PROPERTIES
    LABELS "integration;sbm"
    RUN_SERIAL TRUE
    TIMEOUT 120)
  add_sbm_expected_failure_test(SBM_M3_Fixed_Unsafe_Donor_DT
    inputs_adaptive_anelastic_donor
    "SBM M3 fixed timestep exceeds donor positivity limit"
    "-DOVERRIDE=erf.fixed_dt=0.2")

  add_test(NAME SBM_M3_Adaptive_Anelastic_Stretched_Donor_CFL
    COMMAND ${CMAKE_COMMAND}
      -DERF_EXECUTABLE=$<TARGET_FILE:erf_exec>
      -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_adaptive_anelastic_stretched_donor
      -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}
      -DEXPECT_SECOND_CHECKPOINT=true
      -DSKIP_RESTART=true
      -P ${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/RunM3Smoke.cmake)
  set_tests_properties(SBM_M3_Adaptive_Anelastic_Stretched_Donor_CFL PROPERTIES
    LABELS "integration;sbm"
    RUN_SERIAL TRUE
    TIMEOUT 120)
  add_sbm_expected_failure_test(SBM_M3_Fixed_Unsafe_Anelastic_Stretched_Donor_DT
    inputs_adaptive_anelastic_stretched_donor
    "SBM M3 fixed timestep exceeds donor positivity limit"
    "-DOVERRIDE=erf.fixed_dt=0.05")
  add_sbm_expected_failure_test(SBM_M3_Reject_Legacy_Zero_Transport inputs_zero_carrier
    "erf.sbm_zero_transport_fixture=true is obsolete"
    "-DOVERRIDE=erf.sbm_zero_transport_fixture=true")
  add_sbm_expected_failure_test(SBM_M3_Reject_Anelastic_StaticTerrain
    inputs_periodic_static_terrain
    "SBM M3 does not currently support anelastic StaticFittedMesh"
    "-DOVERRIDE=erf.anelastic=1")
  add_sbm_expected_failure_test(SBM_Native_QC_Write_Fault inputs_zero_carrier
    "SBM owns compact qc/qr; native advection write denied for component range"
    "-DOVERRIDE=erf.sbm_test_native_qc_write_fault=true")
  add_sbm_expected_failure_test(SBM_Invalid_Nonmonotonic_Edges inputs_invalid_edges_nonmonotonic
    "Invalid SBM configuration: edges must be strictly increasing")
  add_sbm_expected_failure_test(SBM_Invalid_Zero_Edge_AutoPivot inputs_invalid_zero_edge_autopivot
    "Invalid SBM configuration: pivots must be finite and positive")
  add_sbm_expected_failure_test(SBM_Implicit_Substepping inputs_zero_carrier
    "SBM M3 requires acoustic substepping_type=None"
    "-DOVERRIDE=erf.substepping_type=Implicit")
  add_sbm_expected_failure_test(SBM_Nonperiodic_Geometry inputs_nonperiodic
    "SBM M3 currently requires triply periodic geometry")
  add_sbm_expected_failure_test(SBM_Immersed_Buildings inputs_zero_carrier
    "SBM M3 requires buildings_type=None"
    "-DOVERRIDE=erf.buildings_type=ImmersedForcing")
  add_sbm_expected_failure_test(SBM_M3_Moving_Terrain inputs_zero_carrier
    "SBM M3 requires static non-EB geometry"
    "-DOVERRIDE=erf.terrain_type=MovingFittedMesh")

  add_test(NAME SBM_M3_Stretched_Static_Mesh
    COMMAND ${CMAKE_COMMAND}
      -DERF_EXECUTABLE=$<TARGET_FILE:erf_exec>
      -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_stretched
      -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}
      -P ${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/RunM3Smoke.cmake)
  set_tests_properties(SBM_M3_Stretched_Static_Mesh PROPERTIES
    LABELS "integration;sbm"
    RUN_SERIAL TRUE
    TIMEOUT 120)

  add_test(NAME SBM_M3_Periodic_StaticTerrain_Mesh
    COMMAND ${CMAKE_COMMAND}
      -DERF_EXECUTABLE=$<TARGET_FILE:erf_exec>
      -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_periodic_static_terrain
      -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}
      -DSKIP_RESTART=true
      -DNUMERIC_HEADER_COMPARE=ON
      -DHEADER_COMPARE_SIGDIGITS=12
      -DHEADER_COMPARE_ULPS=8
      -DEXPECTED_SPECTRUM_VALUES=9.99967817714348408e-07,1.99993563542869682e-06,2.99990345314304565e-06,3.99987127085739363e-06,
      -DEXPECTED_SPECTRUM_MAX_VALUES=1.00000128039821617e-06,2.00000256079643235e-06,3.00000384119464895e-06,4.00000512159286470e-06,
      -DEXPECTED_CORE_VALUES=9.99967817714348506e-01,2.99990345314304420e+02,0.00000000000000000e+00,0.00000000000000000e+00,0.00000000000000000e+00,2.99990345314304522e-06,6.99977472400043928e-06,
      -DEXPECTED_CORE_MAX_VALUES=1.00000128039821612e+00,3.00000384119464741e+02,0.00000000000000000e+00,0.00000000000000000e+00,0.00000000000000000e+00,3.00000384119464852e-06,7.00000896278751364e-06,
      -P ${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/RunM3Smoke.cmake)
  set_tests_properties(SBM_M3_Periodic_StaticTerrain_Mesh PROPERTIES
    LABELS "integration;sbm"
    RUN_SERIAL TRUE
    TIMEOUT 120)

  if(AMReX_PRECISION STREQUAL "DOUBLE" AND
     NOT ERF_ENABLE_CUDA AND NOT ERF_ENABLE_HIP AND NOT ERF_ENABLE_SYCL)
    add_executable(erf_sbm_checkpoint_compare
      ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SBM/ERF_SBMCheckpointCompare.cpp)
    target_link_libraries(erf_sbm_checkpoint_compare PRIVATE AMReX::amrex)
    target_include_directories(erf_sbm_checkpoint_compare PRIVATE ${PROJECT_SOURCE_DIR}/Source)
    target_compile_features(erf_sbm_checkpoint_compare PRIVATE cxx_std_17)

    add_test(NAME SBM_M3_StaticTerrain_Restart_Equivalence
      COMMAND ${CMAKE_COMMAND}
        -DERF_EXECUTABLE=$<TARGET_FILE:erf_exec>
        -DCOMPARE_EXECUTABLE=$<TARGET_FILE:erf_sbm_checkpoint_compare>
        -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_periodic_static_terrain
        -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}
        -P ${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/RunM3RestartEquivalence.cmake)
    set_tests_properties(SBM_M3_StaticTerrain_Restart_Equivalence PROPERTIES
      LABELS "integration;sbm"
      RUN_SERIAL TRUE
      TIMEOUT 180)
  endif()

  foreach(_mode IN ITEMS 1m 2m)
    set(_advection_expected_values)
    if(_mode STREQUAL "2m")
      list(APPEND _advection_expected_values
        -DEXPECTED_SPECTRUM_VALUES=1.50000000000000004e-06,6.00000000000000015e-06,1.80000000000000005e-05,4.80000000000000012e-05,1.00000000000000000e+12,2.00000000000000000e+12,3.00000000000000000e+12,4.00000000000000000e+12,
        -DEXPECTED_CORE_VALUES=1.00000000000000000e+00,2.99999999999999886e+02,0.00000000000000000e+00,0.00000000000000000e+00,0.00000000000000000e+00,7.50000000000000019e-06,6.60000000000000051e-05,)
    endif()
    add_test(NAME SBM_M3_PeriodicAdvection_${_mode}
      COMMAND ${CMAKE_COMMAND}
        -DERF_EXECUTABLE=$<TARGET_FILE:erf_exec>
        -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_periodic_advection_${_mode}
        -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}
        ${_advection_expected_values}
        -P ${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/RunM3Smoke.cmake)
    set_tests_properties(SBM_M3_PeriodicAdvection_${_mode} PROPERTIES
      LABELS "integration;sbm"
      RUN_SERIAL TRUE
      TIMEOUT 120)
  endforeach()

  add_test(NAME SBM_Inadmissible_2M_Restart
    COMMAND ${CMAKE_COMMAND}
      -DERF_EXECUTABLE=$<TARGET_FILE:erf_exec>
      -DTEST_HELPER=$<TARGET_FILE:erf_sbm_checkpoint_corruptor>
      -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_zero_carrier_2m
      -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}
      "-DEXPECTED_DIAGNOSTIC=SBM authoritative restart state is inadmissible.*level=0.*population=0.*bin=0.*constraint=canonical-two-moment-bin-state"
      -P ${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/RunAuthoritativeRestartFailure.cmake)
  set_tests_properties(SBM_Inadmissible_2M_Restart PROPERTIES
    LABELS "integration;sbm;negative"
    RUN_SERIAL TRUE
    TIMEOUT 180)

  if(NOT ERF_ENABLE_CUDA AND NOT ERF_ENABLE_HIP AND NOT ERF_ENABLE_SYCL)
    add_test(NAME Auxiliary_InertTracer_ERF_Smoke
      COMMAND ${CMAKE_COMMAND}
        -DERF_EXECUTABLE=$<TARGET_FILE:erf_exec>
        -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/AuxiliaryState/inputs_inert_tracer
        -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}
        -P ${PROJECT_SOURCE_DIR}/Tests/test_files/AuxiliaryState/RunInertTracerSmoke.cmake)
    set_tests_properties(Auxiliary_InertTracer_ERF_Smoke PROPERTIES
      LABELS "integration;auxiliary"
      RUN_SERIAL TRUE
      TIMEOUT 120)
  endif()
endif()
set_tests_properties(erf_unit_tests_shuffle_repeat
  PROPERTIES
    LABELS "stress"
    TIMEOUT 600
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_subdirectory(Shoc)

if(ERF_ENABLE_MPI)
  if(NOT TARGET MPI::MPI_CXX)
    find_package(MPI REQUIRED COMPONENTS CXX)
  endif()

  if(NOT MPIEXEC_EXECUTABLE)
    message(FATAL_ERROR
      "MPIEXEC_EXECUTABLE is required to register parallel unit tests")
  endif()

  if(NOT DEFINED ERF_PARALLEL_TEST_NRANKS OR "${ERF_PARALLEL_TEST_NRANKS}" STREQUAL "")
    set(ERF_PARALLEL_TEST_NRANKS "2")
  endif()

  # Convert variables that may contain space-separated arguments used to build
  # MPI_COMMANDS in CTestList.txt to lists compatible with COMMAND options below
  # to handle resource managers with subcommands e.g., flux run
  separate_arguments(erf_mpiexec_executable UNIX_COMMAND "${MPIEXEC_EXECUTABLE}")
  separate_arguments(erf_mpiexec_numproc_flag UNIX_COMMAND "${MPIEXEC_NUMPROC_FLAG}")
  separate_arguments(erf_mpiexec_preflags UNIX_COMMAND "${MPIEXEC_PREFLAGS}")
  separate_arguments(erf_mpiexec_postflags UNIX_COMMAND "${MPIEXEC_POSTFLAGS}")


  set(erf_parallel_test_mpiexec_preflags)

  if(DEFINED ERF_PARALLEL_TEST_MPIEXEC_PREFLAGS AND
     NOT "${ERF_PARALLEL_TEST_MPIEXEC_PREFLAGS}" STREQUAL "")
    separate_arguments(erf_parallel_test_user_mpiexec_preflags
      UNIX_COMMAND "${ERF_PARALLEL_TEST_MPIEXEC_PREFLAGS}")
    list(APPEND erf_parallel_test_mpiexec_preflags
      ${erf_parallel_test_user_mpiexec_preflags})
  endif()

  if(ERF_PARALLEL_TEST_ALLOW_OVERSUBSCRIBE)
    execute_process(
      COMMAND ${erf_mpiexec_executable} --version
      OUTPUT_VARIABLE erf_mpiexec_version_out
      ERROR_VARIABLE erf_mpiexec_version_err
      RESULT_VARIABLE erf_mpiexec_version_result
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_STRIP_TRAILING_WHITESPACE
    )

    string(TOLOWER
      "${erf_mpiexec_version_out} ${erf_mpiexec_version_err}"
      erf_mpiexec_version_text)

    if(erf_mpiexec_version_result EQUAL 0 AND
       erf_mpiexec_version_text MATCHES "open mpi|openrte|prrte")
      list(APPEND erf_parallel_test_mpiexec_preflags
        "--map-by" ":OVERSUBSCRIBE")
      message(STATUS
        "Parallel GoogleTest MPI oversubscription enabled for Open MPI/PRRTE")
    else()
      message(STATUS
        "ERF_PARALLEL_TEST_ALLOW_OVERSUBSCRIBE is ON, but mpiexec does not look like Open MPI/PRRTE; no oversubscribe flag was added.")
    endif()
  endif()

  add_executable(${erf_parallel_exe_name})
  target_sources(${erf_parallel_exe_name}
    PRIVATE
      ${PROJECT_SOURCE_DIR}/Exec/ERF_Prob.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/ERF_GTestParallelMain.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/ERF_GTestParallelHelloWorld.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/ERF_GTestRuntimeStubs.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestMOSTAverageParallel.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestSurfaceLayerParallel.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestProvenanceParallel.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SatAdj/ERF_GTestSatAdjParallel.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SAM/ERF_GTestSAMParallel.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SBM/ERF_GTestSBMTransportParallel.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/Utils/HashRng/ERF_GTestHashRngParallel.cpp
      # Motivation: planar surface-layer fields on a BoxArray split in z hold
      # duplicate 2D boxes, and every copy must be filled from the computed
      # surface copy, across ranks and periodic boundaries.
      ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ERF_GTestPlanarBoundaryParallel.cpp
      # Motivation: rebalance_columns integrated each box of a BoxArray split in
      # z from its own lowest cell, so the hydrostatic base state depended on
      # where the columns were split.
      ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ERF_GTestRebalanceColumnsParallel.cpp
      # Motivation: the hydrostatic base state over terrain integrated each box of
      # a level split in z from the ghost cell below it, so a fine level's base
      # state depended on where its boxes were split.
      ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ERF_GTestColumnBandsParallel.cpp
      # The column boxes are distributed afresh, so check the copies across ranks too
      ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ERF_GTestTerrainFineColumns.cpp
      # The boxes are spread over the ranks, so a k may be held by another rank only
      ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ERF_GTestReduceToMaxPerHeight.cpp
  )
  if(ERF_ENABLE_NETCDF)
    # The terrain reader opens the fixture only on the I/O rank and then
    # broadcasts the result.  Keep this source in the MPI executable so the
    # WPS edge/corner regression exercises every rank.
    target_sources(${erf_parallel_exe_name} PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestTerrainNetCDF.cpp
    )
  endif()
  build_erf_exe(${erf_parallel_exe_name} NO_ERF_MAIN)
  target_link_libraries(${erf_parallel_exe_name}
    PRIVATE
      GTest::gtest
      MPI::MPI_CXX
  )

  set(erf_parallel_ctest_file
      "${CMAKE_CURRENT_BINARY_DIR}/${erf_parallel_exe_name}_tests_$<CONFIG>.cmake")
  set(erf_parallel_ctest_wrapper
      "${CMAKE_CURRENT_BINARY_DIR}/${erf_parallel_exe_name}_tests.cmake")

  # CTest includes this stable wrapper. Multi-config generators build one
  # generated test file per configuration; CTest chooses the requested file
  # through CTEST_CONFIGURATION_TYPE. A single-config build falls back to
  # its only generated file when ctest is invoked without -C.
  file(WRITE "${erf_parallel_ctest_wrapper}" [=[
set(_erf_parallel_test_files)
if(DEFINED CTEST_CONFIGURATION_TYPE AND
   NOT "${CTEST_CONFIGURATION_TYPE}" STREQUAL "")
  set(_erf_parallel_test_file
      "${CMAKE_CURRENT_LIST_DIR}/erf_parallel_tests_tests_${CTEST_CONFIGURATION_TYPE}.cmake")
  if(NOT EXISTS "${_erf_parallel_test_file}")
    message(FATAL_ERROR
      "No generated MPI CTest file exists for configuration '${CTEST_CONFIGURATION_TYPE}'")
  endif()
  include("${_erf_parallel_test_file}")
else()
  file(GLOB _erf_parallel_test_files
       "${CMAKE_CURRENT_LIST_DIR}/erf_parallel_tests_tests_*.cmake")
  list(LENGTH _erf_parallel_test_files _erf_parallel_test_file_count)
  if(_erf_parallel_test_file_count EQUAL 1)
    include("${_erf_parallel_test_files}")
  elseif(_erf_parallel_test_file_count EQUAL 0)
    message(FATAL_ERROR "No generated MPI CTest file exists")
  else()
    message(FATAL_ERROR
      "Multiple generated MPI CTest files exist; use ctest -C <configuration>")
  endif()
endif()
]=])

  # Encode CMake lists as comma-separated values for the -D arguments below.
  # COMMAND_EXPAND_LISTS would otherwise split a single -D value into several
  # unrelated command-line arguments.
  string(REPLACE ";" "," erf_parallel_test_mpiexec_preflags_arg
         "${erf_parallel_test_mpiexec_preflags}")
  string(REPLACE ";" "," erf_parallel_test_ranks_arg
         "${ERF_PARALLEL_TEST_NRANKS}")
  string(REPLACE ";" "," erf_mpiexec_preflags_arg "${MPIEXEC_PREFLAGS}")
  string(REPLACE ";" "," erf_mpiexec_postflags_arg "${MPIEXEC_POSTFLAGS}")

  get_target_property(erf_parallel_test_emulator
    ${erf_parallel_exe_name} CROSSCOMPILING_EMULATOR)
  if(NOT erf_parallel_test_emulator OR
     erf_parallel_test_emulator MATCHES "-NOTFOUND$")
    set(erf_parallel_test_emulator)
  endif()
  # Encode target-property lists as comma-separated -D values. This preserves
  # semicolon-delimited list boundaries and spaces inside emulator arguments,
  # matching the existing MPI flag encoding used by this custom generator.
  string(REPLACE ";" "," erf_parallel_test_emulator_arg
         "${erf_parallel_test_emulator}")

  add_custom_command(TARGET ${erf_parallel_exe_name} POST_BUILD
    COMMAND ${CMAKE_COMMAND}
      -DTEST_EXECUTABLE=$<TARGET_FILE:${erf_parallel_exe_name}>
      -DOUTPUT_FILE=${erf_parallel_ctest_file}
      -DREPORT_DIRECTORY=${ERF_GTEST_REPORT_ROOT}/${erf_parallel_exe_name}
      -DMPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}
      -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
      -DMPIEXEC_PREFLAGS=${erf_mpiexec_preflags_arg}
      -DMPIEXEC_POSTFLAGS=${erf_mpiexec_postflags_arg}
      -DERF_MPI_PREFLAGS=${erf_parallel_test_mpiexec_preflags_arg}
      -DERF_MPI_RANKS=${erf_parallel_test_ranks_arg}
      -DERF_CROSSCOMPILING_EMULATOR=${erf_parallel_test_emulator_arg}
      -DTEST_TIMEOUT=600
      -P ${CMAKE_CURRENT_LIST_DIR}/GenerateERFParallelCTest.cmake
    COMMAND_EXPAND_LISTS
    VERBATIM
  )

  set_property(DIRECTORY APPEND PROPERTY TEST_INCLUDE_FILES
    "${erf_parallel_ctest_wrapper}")

  add_executable(erf_mpi_gtest_listener_self_test
    ERF_GTestMPIListenerSelfTest.cpp
  )
  target_link_libraries(erf_mpi_gtest_listener_self_test
    PRIVATE
      GTest::gtest
      MPI::MPI_CXX
  )

  get_target_property(erf_listener_self_test_emulator
    erf_mpi_gtest_listener_self_test CROSSCOMPILING_EMULATOR)
  if(NOT erf_listener_self_test_emulator OR
     erf_listener_self_test_emulator MATCHES "-NOTFOUND$")
    set(erf_listener_self_test_emulator)
  endif()
  string(REPLACE ";" "," erf_listener_self_test_emulator_arg
         "${erf_listener_self_test_emulator}")

  add_test(NAME erf_mpi_gtest_listener_self_test
    COMMAND ${CMAKE_COMMAND}
      -DMPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}
      -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
      -DMPIEXEC_PREFLAGS=${erf_mpiexec_preflags_arg}
      -DMPIEXEC_POSTFLAGS=${erf_mpiexec_postflags_arg}
      -DERF_MPI_PREFLAGS=${erf_parallel_test_mpiexec_preflags_arg}
      -DERF_CROSSCOMPILING_EMULATOR=${erf_listener_self_test_emulator_arg}
      -DTEST_EXECUTABLE=$<TARGET_FILE:erf_mpi_gtest_listener_self_test>
      -P ${CMAKE_CURRENT_LIST_DIR}/RunERFMPIListenerSelfTest.cmake
  )
  set_tests_properties(erf_mpi_gtest_listener_self_test
    PROPERTIES
      LABELS "infrastructure"
      PROCESSORS 2
      TIMEOUT 120
  )
endif()
