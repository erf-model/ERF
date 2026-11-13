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
    # Motivation: the ghost cells beside an ingested boundary plane along a periodic
    # direction must take the plane at the periodic image, not the near side of the seam.
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestBndryPlaneIndex.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestMultiScalarSlots.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestInflowProfile.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestSurfaceLayerStress.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestSurfaceLayer.cpp
    # Motivation: the roughness inverted from surface fluxes must stay inside its bounds.
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestInvertedRoughness.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestMOSTAverage.cpp
    # Motivation: the MOST flux kernels must stay bounded and correctly signed at the
    # stable fold of a prescribed flux and when the log-law bracket would change sign.
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestMOSTIteration.cpp
    # Motivation: the QNSE stable-side functions must match their closed forms, keep the
    # default and the unstable side untouched, and move the stable fold the way the analysis says.
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestMOSTStableFunctions.cpp
    # Motivation: the MOST reference height must use the true stretched cell
    # centers, and a query on a mesh face must be found rather than abort.
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestMOSTAverageZref.cpp
    # Motivation: a BTF fine level's mesh must not depend on how its grids are split in z.
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ERF_GTestTerrainFineColumns.cpp
    # Motivation: the input sponge's per-height reduction read every k from every box, so a
    # refined region that does not reach the top was read outside its data.
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ERF_GTestReduceToMaxPerHeight.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestWrfBoundaryDensity.cpp
    # Motivation: the wrfbdy time-interpolation weight must follow the elapsed time to roundoff,
    # not to the ~2.4e-7 s spacing of doubles at epoch times.
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestBdyTimeWeights.cpp
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
    # Motivation: with an EB that cuts no cell, EB momentum advection must equal the non-EB
    # operator, map factors included: x-face mass fluxes divide by m_y, y-face ones by m_x.
    ${CMAKE_CURRENT_SOURCE_DIR}/Advection/ERF_GTestEBMomMapFactors.cpp
    # Motivation: the z-momentum box must be shrunk only at the bottom and top of the
    # domain, so that a grid split in z leaves no w-face without a source term.
    ${CMAKE_CURRENT_SOURCE_DIR}/Advection/ERF_GTestZmomBoxDecomposition.cpp
    # Motivation: the chamber profile and perturbation must remain analytic,
    # bounded, and independent of MPI decomposition and tile shape.
    ${CMAKE_CURRENT_SOURCE_DIR}/Prob/ERF_GTestCloudChamber.cpp
    # Motivation: the ABL warm bubble (prob.bubble_*) must follow WRF's cos^2 form on each
    # cell's own height, and a compressible bubble must be made at constant pressure.
    ${CMAKE_CURRENT_SOURCE_DIR}/Prob/ERF_GTestABLWarmBubble.cpp
    # Motivation: RICO's moisture source must be continuous at prob.moisture_source_cutoff and
    # default to RICO's 2980 m (erf-model/ERF#4626).
    ${CMAKE_CURRENT_SOURCE_DIR}/Prob/ERF_GTestRICOMoistureSource.cpp
    # Motivation: production scalar diffusion must use the exact low/high
    # physical-face stencil in every rotated Cartesian direction.
    ${CMAKE_CURRENT_SOURCE_DIR}/Diffusion/ERF_GTestStateDiffusionDirichlet.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Diffusion/ERF_GTestTKEBuoyancySource.cpp
    # Motivation: the YSUNew PBL height must interpolate a bulk Richardson crossing at the top
    # cell, give a column with no crossing the top of the scan, and give a surface cell already
    # above Ribcr the floor without dividing 0 by 0 (erf-model/ERF#4458).
    ${CMAKE_CURRENT_SOURCE_DIR}/Diffusion/ERF_GTestPBLRibCrossing.cpp
    # Motivation: reusable scalar diffusion must honor explicit scalar, coefficient,
    # flux, and tendency components while preserving ERF's N/S/T mapped discretizations.
    ${CMAKE_CURRENT_SOURCE_DIR}/Diffusion/ERF_GTestScalarDiffusion.cpp
    # Motivation: the k-eqn RANS closure (AL01) must keep its stability
    # functions finite and its length-scale limits exact for any Rt.
    ${CMAKE_CURRENT_SOURCE_DIR}/Diffusion/ERF_GTestRANSClosure.cpp
    # Motivation: terrain-fitted momentum diffusion must apply K_h to the projected horizontal
    # stresses of the zeta-face flux and take the cell-centred du/dz from the cell's own faces
    # (erf-model/ERF#4214), so the operator is exact for bilinear fields and dissipates energy.
    ${CMAKE_CURRENT_SOURCE_DIR}/Diffusion/ERF_GTestTerrainStress.cpp
    # Motivation: erf.implicit_terrain_metric must put K_h times the terrain-metric factor into
    # the implicit vertical solves on the interior faces and take the same compact term out of
    # the explicit fluxes, and nothing on the bottom and top faces.
    ${CMAKE_CURRENT_SOURCE_DIR}/Diffusion/ERF_GTestTerrainImplicitMetric.cpp
    # Motivation: the implicit vertical momentum solve must close a no-slip z wall with the
    # second-order one-sided gradient the explicit strain uses, or the steady state of a
    # wall-bounded flow (PoiseuilleFlow) drifts at the wall cells.
    ${CMAKE_CURRENT_SOURCE_DIR}/Diffusion/ERF_GTestImplicitWallStencil.cpp
    # Motivation: the implicit vertical momentum solves must carry the YSU countergradient flux
    # through the lower face of the top row and both faces of a coarse/fine row, or the column
    # loses momentum through the face below the top cell (erf-model/ERF#4822).
    ${CMAKE_CURRENT_SOURCE_DIR}/Diffusion/ERF_GTestImplicitMomCountergradient.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/LinearSolvers/ERF_GTestTerrainPoisson.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/LinearSolvers/ERF_GTestTerrainMLMG.cpp
    # Motivation: the opt-in WRF Smagorinsky-2D limits must reproduce WRF's smag2d_km numbers
    # (slope factor, deformation threshold, cap and their order), limit K_h before the scalar
    # diffusivities are derived from it, and the diffusive time-step rates must follow their
    # documented formulas, with the vertical part dropped only where every implicit stage is
    # shown to stay bounded.
    ${CMAKE_CURRENT_SOURCE_DIR}/Diffusion/ERF_GTestSmag2DLimiters.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Diagnostics/ERF_GTestSurfaceFluxDiagnostics.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Diagnostics/ERF_GTestSeaLevelPressure.cpp
    # Motivation: a Noah-MP 2 m tile with zero weight must not reject the bundle, and a
    # station must leave -999 neighbours out of its interpolation rather than average them.
    ${CMAKE_CURRENT_SOURCE_DIR}/Diagnostics/ERF_GTestSurfaceMissingValues.cpp
    # Motivation: the stable-side z/L in the 2 m profile must stay bounded for a very
    # small positive Obukhov length.
    ${CMAKE_CURRENT_SOURCE_DIR}/Diagnostics/ERF_GTestNearSurfaceStability.cpp
    # Motivation: the 10 m winds must take the surface-layer momentum profile's speed
    # in the first-cell wind's direction, turn earth-relative with the grid, and be
    # missing where the profile cannot be evaluated.
    ${CMAKE_CURRENT_SOURCE_DIR}/Diagnostics/ERF_GTestNearSurfaceWind.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/ERF_GTestCanopyBiophysics.cpp
    # Motivation: the k-eqn RANS canopy terms must keep the kynema-sgf
    # formulas, cap L only inside the forest and never remove more than k.
    ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/ERF_GTestCanopyTurbulence.cpp
    # Motivation: buoyancy dispatch must preserve selector semantics, active
    # anelastic formulations, and neutral dry type-4 primitive-theta behavior.
    ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/ERF_GTestBuoyancy.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/ERF_GTestForestDrag.cpp
    # Motivation: a calm cell above an immersed wall gave the wall law u* = 0, an Obukhov length
    # of zero and z / 0 (erf-model/ERF#4016); the stability estimate is now bounded as on flat ground.
    ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/ERF_GTestImmersedWallStability.cpp
    # Motivation: the building wall function's Louis (1979) factors must be 1
    # at neutral, follow the published closed forms on both sides and weaken
    # mixing monotonically with the bulk Richardson number.
    ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/ERF_GTestIBSEBLouis.cpp
    # Motivation: the immersed-forcing temperature kernels applied the MOST similarity functions whatever
    # erf.if_stability_correction said, while the momentum kernels gate on it (erf-model/ERF#4931); off,
    # every kernel is the log law.
    ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/ERF_GTestImmersedScalarStability.cpp
    # Motivation: the fraction-stress immersed wall law (erf.if_wall_form = fraction_stress) must
    # force the wall cell (the fluid cell on a face-aligned surface too, which the legacy law leaves
    # unforced), close its bottom face to diffusion and measure the RANS wall distance from the surface.
    ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/ERF_GTestImmersedFractionStress.cpp
    # Motivation: the implicit solves must carry no flux through the wall face of the fraction-stress
    # wall law, and a partial wall cell must dissipate like flat ground's first cell.
    ${CMAKE_CURRENT_SOURCE_DIR}/Diffusion/ERF_GTestImmersedWallFace.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/ERF_GTestImmersedBuildingCutoff.cpp
    # Motivation (erf-model/ERF#4887): the building kernels and the implicit projection read each
    # face and its masks through one header; with erf.if_implicit_projection the kernels keep the
    # wall law only, and the drag they drop must equal, face by face, the weight the projection
    # applies (checked with and without the snap and the face-centred blanking).
    ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/ERF_GTestImmersedBuildingMasks.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/ERF_GTestImmersedBuildingImplicit.cpp
    # Motivation: a refined level's buildings must keep level 0's numbers for the
    # materials by building, even when its tallest column is a spire level 0
    # leaves fluid, or the level holds only some of the buildings.
    ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/ERF_GTestIBSEBLevelNumbering.cpp
    # Motivation: the ray cast of the building faces reads a per-rank patch of the
    # column map with a halo sized for the shallowest ray, instead of a map of the whole
    # domain on every rank (erf-model/ERF#4020): a ray that leaves the patch is in open
    # country, a patch spanning a periodic domain wraps round, the buildings are
    # labelled across boxes and periodic seams in one numbering whatever the
    # decomposition, and a refined level sees the buildings of the level below.
    ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/ERF_GTestIBSEBRayMap.cpp
    # Motivation: with erf.ibseb.radiation = two_stream every face must read its own column
    # at its own height (roof interface, wall mean of its cell's two), on a level of several
    # boxes whose columns all differ.
    ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/ERF_GTestIBSEBTwoStreamFaces.cpp
    # Motivation: the building faces took the air temperature and the Exner function of
    # their fluid cell without its water vapour (erf-model/ERF#4481); they now match the
    # two-stream column's equation of state for the same moist cell.
    ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/ERF_GTestIBSEBMoistAir.cpp
    # Motivation: with erf.ibseb.two_stream_buildings the ground balance takes the beam
    # off the ground the buildings shade; the shadow must be the buildings' geometry,
    # across a periodic edge too, and under an overhang, in every box of a level.
    ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/ERF_GTestIBSEBGroundShadow.cpp
    # Motivation: with erf.ibseb.ground_view_shadow a face takes the beam the ground
    # reflects through the sunlit part of its ground view; that part must be the
    # buildings' geometry (a closed form behind a long wall, a count from the geometry
    # for a building set), and the answer must not depend on where a rank's map ends.
    ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/ERF_GTestIBSEBGroundViewShadow.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/ERF_GTestIBSEBGroundWalls.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/ERF_GTestIBSEBRoofs.cpp
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
    ${CMAKE_CURRENT_SOURCE_DIR}/Radiation/ERF_GTestTwoStreamRoof.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Radiation/ERF_GTestTwoStreamRadiationDriver.cpp
    # Motivation: the Simple longwave scheme must stay inside the column when the
    # moist layer reaches the domain top (erf-model/ERF#4484), and give the top
    # cell its closed-form heating there.
    ${CMAKE_CURRENT_SOURCE_DIR}/Radiation/ERF_GTestRadiationSimple.cpp
    # Motivation: with erf.ibseb.two_stream_buildings the sweep writes no heating in a
    # building's solid cells, and the ground balance holds a footprint and takes the
    # absorbed beam (not the diffuse light) off the shaded ground, once per sweep.
    ${CMAKE_CURRENT_SOURCE_DIR}/Radiation/ERF_GTestTwoStreamBuildingCover.cpp
    # Motivation: with erf.ibseb.sun_mode = two_stream the building faces must
    # see the two-stream columns' sun: the same zenith, a morning sun in the
    # east, and the start date read the same way.
    ${CMAKE_CURRENT_SOURCE_DIR}/Radiation/ERF_GTestIBSEBTwoStreamSun.cpp
    # Motivation: the building faces' solar declination and equation of time
    # must be Spencer's (1971) published series on his 365-day angle, which
    # puts the June-solstice sun near 23.44 deg, and the date must advance
    # at 00:00 UTC in a run of several days.
    ${CMAKE_CURRENT_SOURCE_DIR}/Radiation/ERF_GTestIBSEBSolarDeclination.cpp
    # Motivation: the surface energy balance must take H and LE from a land
    # model, else from the surface layer's applied flux in W/m^2, else from
    # the scalar defaults; LE stays at its default without a moisture model.
    ${CMAKE_CURRENT_SOURCE_DIR}/Radiation/ERF_GTestSEBTurbulentFlux.cpp
    # Motivation: with the moisture coupling, the force-restore surface evaporates at
    # beta times the potential rate, beta linear between Noah-MP's wilting point and
    # field capacity for the soil category, so it dries as its water content falls.
    ${CMAKE_CURRENT_SOURCE_DIR}/Radiation/ERF_GTestSEBSoilMoisture.cpp
    # Motivation: the force-restore hydrology must add rain to the surface layer, keep the
    # column's water budget closed with runoff, drainage and any deficit, and relax the
    # drainage and the deep temperature exactly over the step (NP89, MN96), and carry the
    # rain of a step the balance skips into the next.
    ${CMAKE_CURRENT_SOURCE_DIR}/Radiation/ERF_GTestForceRestoreHydrology.cpp
    # Motivation: the forcing the two-stream model hands Noah-MP must be the
    # downwelling SW and LW at the surface interface and a zenith cosine floored
    # at zero, copied into the land model's k = 0 plane and nowhere else.
    ${CMAKE_CURRENT_SOURCE_DIR}/Radiation/ERF_GTestTwoStreamLandForcing.cpp
    # Motivation: the building faces (erf.ibseb.radiation = two_stream) must get the
    # Beer-Lambert beam at every interface up to the highest one they read, blended
    # clear/cloudy, no beam with the shortwave off or the sun down, only for the step
    # of the sweep, and each face its own column's fluxes at its own height.
    ${CMAKE_CURRENT_SOURCE_DIR}/Radiation/ERF_GTestTwoStreamCanopyForcing.cpp
    # Motivation: RRTMGP consumes absolute surface temperature; a
    # SurfaceLayer fallback must convert theta at the physical surface
    # pressure diagnosed from the lowest atmospheric cell,
    # respect land/water precedence, and preserve valid land LSM temperatures.
    ${CMAKE_CURRENT_SOURCE_DIR}/Radiation/ERF_GTestRRTMGP_SurfaceTemperature.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Radiation/ERF_GTestOrbitalGeometry.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Radiation/ERF_GTestTopographicShading.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestPlotfile2D.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestCloudOptics.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestPlotfileSelection.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestStationSampler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestProvenance.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestCheckpoint.cpp
    # Motivation: a deck may include several files with FILE = <name>; the duplicate-input check
    # must not refuse it.
    ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestInputsName.cpp
    # Motivation: a restart used to write its start-up boundary plane as step 0, replacing the
    # first plane and leaving a time.dat that erf.input_bndry_planes cannot read.
    ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestWriteBndryPlanes.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ERF_GTestEpochTime.cpp
    # Motivation: SLM's initial soil state from a met_em file must interpolate the source
    # dataset's soil layers in depth to the Noah layers, and the monthly vegetation to the
    # start date the way real.exe does, across the year's end and in leap years.
    ${CMAKE_CURRENT_SOURCE_DIR}/Initialization/ERF_GTestMetgridLandUtils.cpp
    # Motivation: erf.metgrid_force_sfc_k must drop the same origin levels whether or not the
    # column has origin levels below the surface (erf-model/ERF#4607).
    ${CMAKE_CURRENT_SOURCE_DIR}/Initialization/ERF_GTestMetgridForceSfc.cpp
    # Parent-driven inflow masks must follow height and wind reversals, including
    # calm periods, without applying corner increments twice.
    ${CMAKE_CURRENT_SOURCE_DIR}/Initialization/ERF_GTestAutomaticPerturbation.cpp
    # Motivation: a station placed by latitude/longitude must land where the
    # lat/lon arrays put it on a rotated grid and on a grid across 180 degrees.
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ERF_GTestLatLonMap.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ERF_GTestGridUtils.cpp
    # Motivation: erf.initial_dz without a stretching ratio must solve for the ratio that
    # closes the grid on prob_hi, and an explicit ratio must not be solved for.
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ERF_GTestInitZLevels.cpp
    # Motivation: a malformed, mis-ordered or too small terrain file must stop the run with
    # the problem named, and steep terrain next to an inflow or outflow face must be reported.
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ERF_GTestTerrainFileChecks.cpp
    # Motivation: a BTF fine level's mesh must not depend on how its grids are split in z.
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ERF_GTestTerrainFineColumns.cpp
    # Motivation: the average down over immersed-forcing levels must weight the fine cells by
    # their fluid mass and leave coarse cells that are solid on their own level alone.
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ERF_GTestImmersedAverageDown.cpp
    # Motivation: over immersed-forcing terrain with the fraction-stress wall law the RANS wall
    # distance must be measured from the surface the wall law uses, column by column.
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ERF_GTestImmersedWallDist.cpp
    # Motivation: the coarse/fine surface-mismatch warning must count only the coarse cells along
    # the coarse/fine faces whose fine solid fractions spread by at least the threshold.
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ERF_GTestImmersedCFFaces.cpp
    # Motivation: the log-law coarse/fine fill must keep the flux of every coarse face and share
    # the wall face's flux by ln(1 + z/z0), closing the solid sub-faces of immersed terrain.
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestCFLogLawFill.cpp
    # Motivation: the relax band of the coarse/fine fill patcher must follow the top face of a
    # level that stops below the domain top, and not appear along the physical top or bottom.
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestFillPatcherRelaxBand.cpp
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
    # Motivation: a finer level's land state must be interpolated land from land and water
    # from water (erf-model/ERF#4390), keeping each coarse mean and cell_cons_interp elsewhere.
    ${CMAKE_CURRENT_SOURCE_DIR}/LandSurfaceModel/ERF_GTestTypedInterp.cpp
)
if(ERF_DIM STREQUAL "3")
  # Motivation: the auxiliary mapped state must reproduce ERF's computational
  # divergence, exact host stage recipes, canonical metric transfers, and ledger.
  target_sources(${erf_exe_name} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/AuxiliaryState/ERF_GTestAuxiliaryMappedTransport.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SBM/ERF_GTestSBMTransport.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SBM/ERF_GTestSBMBoundary.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SBM/ERF_GTestSBMParticleMixing.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SBM/ERF_GTestSBMMixingTransport.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SBM/ERF_GTestSBMAMRTransfer.cpp)
endif()
if(ERF_ENABLE_NETCDF)
  target_sources(${erf_exe_name} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestNCInterface.cpp
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
# These exercise super-droplet projection and process kernels, whose headers are
# only available when particles are enabled.
if(ERF_ENABLE_PARTICLES)
  target_sources(${erf_exe_name} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SDM/ERF_GTestAnelasticThermodynamics.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SDM/ERF_GTestCellVolume.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SDM/ERF_GTestSpeciesProjection.cpp
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
    ${CMAKE_CURRENT_SOURCE_DIR}/LandSurfaceModel/Noah-MP/ERF_GTestNoahMPSnowFields.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/LandSurfaceModel/Noah-MP/ERF_GTestNoahMPCheckpointLayout.cpp
  )
endif()
if(ERF_ENABLE_MOVING_BODIES)
  # Motivation: the actuator velocity sampler must reproduce a field linear in x, y and the
  # physical height exactly, from each component's own staggered grid, on uniform, stretched
  # and terrain-following meshes, and find every point in exactly one box however the domain
  # is split.
  target_sources(${erf_exe_name} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/MovingBodies/ERF_GTestActuatorSampling.cpp
  )
  # Motivation: a point force spread onto the momentum sources must integrate back to the
  # force exactly on every component, with the terrain cell volumes, under the ground
  # cut-off, for several points, and however the domain is split.
  target_sources(${erf_exe_name} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/MovingBodies/ERF_GTestActuatorSpreading.cpp
  )
  # Motivation: the wake sampling lines must sit at the requested distances behind a body
  # along its axis, their running average must be the mean of the samples (an analytic wake
  # comes back exactly), and that average must survive a checkpoint round trip.
  target_sources(${erf_exe_name} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/MovingBodies/ERF_GTestWakeLines.cpp
  )
  # Motivation: a body's running statistics must be the exact mean, rms, minimum and maximum
  # of the samples since the averaging start and must survive a checkpoint round trip.
  target_sources(${erf_exe_name} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/MovingBodies/ERF_GTestRunningStats.cpp
  )
  # Motivation: a structure model that registers itself from its own file must be listed and
  # created by the name erf.structures.model gives; an unknown name must stop the run with
  # the names it accepts, and a name registered twice or an empty name must stop it too.
  target_sources(${erf_exe_name} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/MovingBodies/ERF_GTestStructureRegistry.cpp
  )
  # Motivation: a structure model's source restricted onto the coarser levels must put the
  # same force on the air as on the model's level, over a hill (with the mesh Jacobian),
  # one and two levels down, with anisotropic ratios and a fine level split into boxes and
  # across a periodic boundary; without the Jacobian the hill must change the force.
  target_sources(${erf_exe_name} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/MovingBodies/ERF_GTestStructureCoarseSource.cpp
  )
  # A structure model directory may define erf_structure_model_add_unit_tests(target) in its
  # ERFStructureModel.cmake to build its own GoogleTest files into this executable, with
  # ERF's test main; ERF_STRUCTURE_MODEL_DIR holds the directory during the call. The file
  # is included again here because CMake keeps only the last definition of a function, and
  # the empty definition before each include keeps a directory without the hook from
  # calling the one of the directory before it. Called before build_erf_exe, so that a
  # CUDA build compiles the added files as CUDA too. A relative directory is resolved
  # against Exec/, as build_erf_lib resolves it, so both name the same folder. The file
  # must therefore be safe to include twice: definitions only, no targets or options.
  foreach(_structure_dir IN LISTS ERF_STRUCTURE_MODEL_DIRS)
    get_filename_component(_structure_dir "${_structure_dir}" ABSOLUTE BASE_DIR "${ERF_SOURCE_DIR}/Exec")
    function(erf_structure_model_add_unit_tests target)
    endfunction()
    set(ERF_STRUCTURE_MODEL_DIR "${_structure_dir}")
    include("${_structure_dir}/ERFStructureModel.cmake")
    erf_structure_model_add_unit_tests(${erf_exe_name})
  endforeach()
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

if(TARGET erf_exec AND ERF_DIM STREQUAL "3")
  foreach(_factory_choice IN ITEMS default native)
    if(_factory_choice STREQUAL "native")
      set(_expected_factory "native face-centered")
      set(_factory_override erf.buildings_fc_factories=true)
    else()
      set(_expected_factory "cell-centered")
      set(_factory_override)
    endif()
    add_test(NAME ImmersedBuildingFactories_${_factory_choice}
      COMMAND ${CMAKE_COMMAND}
        -DERF_EXECUTABLE=$<TARGET_FILE:erf_exec>
        -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/ImmersedForcing/inputs_building_fc_factories
        "-DEXPECTED_FACTORY=${_expected_factory}"
        -DFACTORY_OVERRIDES=${_factory_override}
        -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/RunImmersedBuildingFactories.cmake)
    set_tests_properties(ImmersedBuildingFactories_${_factory_choice} PROPERTIES
      LABELS "integration;immersed-forcing"
      RUN_SERIAL TRUE
      TIMEOUT 120)
  endforeach()

  if(NOT AMReX_HDF5)
    # At 50 m the masks are active; at 60 m the masks are zero and winds stay unchanged.
    add_executable(erf_immersed_cutoff_check
      ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/ERF_CheckImmersedBuildingCutoff.cpp)
    # AMReX's headers need nvcc in a CUDA build
    if(ERF_ENABLE_CUDA)
      setup_target_for_cuda_compilation(erf_immersed_cutoff_check)
    endif()
    target_link_libraries(erf_immersed_cutoff_check PRIVATE AMReX::amrex)
    target_compile_features(erf_immersed_cutoff_check PRIVATE cxx_std_17)
    add_test(NAME ImmersedBuildingResolutionCutoff
      COMMAND ${CMAKE_COMMAND}
        -DERF_EXECUTABLE=$<TARGET_FILE:erf_exec>
        -DCHECK_EXECUTABLE=$<TARGET_FILE:erf_immersed_cutoff_check>
        -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/ImmersedForcing/inputs_building_fc_factories
        -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/SourceTerms/RunImmersedBuildingCutoff.cmake)
    set_tests_properties(ImmersedBuildingResolutionCutoff PROPERTIES
      LABELS "integration;immersed-forcing"
      RUN_SERIAL TRUE
      TIMEOUT 120)
  endif()
endif()

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

  # The dust coupling's dust scalars and the fire coupling's smoke scalar add
  # components ahead of the moisture ones in the checkpoint (NSCALARS in
  # ERF_IndexDefines.H); every RunM3Smoke test passes the count so the runner
  # can insert them in the expected values.
  set(_sbm_extra_scalars 0)
  if(ERF_ENABLE_FIRE)
    math(EXPR _sbm_extra_scalars "${_sbm_extra_scalars} + 1")
  endif()
  if(ERF_ENABLE_DUST)
    math(EXPR _sbm_extra_scalars "${_sbm_extra_scalars} + ${ERF_NUM_DUST_SCALARS}")
  endif()

  add_test(NAME SBM_M3_ZeroCarrier_ERF_Smoke
    COMMAND ${CMAKE_COMMAND}
      -DERF_EXECUTABLE=$<TARGET_FILE:erf_exec>
      -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_zero_carrier
      -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}
      -DEXTRA_SCALAR_COMPONENTS=${_sbm_extra_scalars}
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
      -DEXTRA_SCALAR_COMPONENTS=${_sbm_extra_scalars}
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
  # Multilevel SBM outside the qualified two-level envelope is rejected at
  # construction, in production and qualification builds alike.
  add_sbm_expected_failure_test(SBM_AMR_Reject_ThreeLevels
    inputs_amr_production_2level
    "SBM supports at most two AMR levels"
    "-DOVERRIDE=amr.max_level=2")
  add_sbm_expected_failure_test(SBM_AMR_Reject_RefRatio4
    inputs_amr_production_2level
    "SBM AMR supports amr.ref_ratio = 2 in every direction only"
    "-DOVERRIDE=amr.ref_ratio=4")
  add_sbm_expected_failure_test(SBM_AMR_Reject_OneWayCoupling
    inputs_amr_production_2level
    "SBM AMR requires erf.coupling_type = TwoWay"
    "-DOVERRIDE=erf.coupling_type=OneWay")
  add_sbm_expected_failure_test(SBM_AMR_Reject_StretchedDz
    inputs_amr_production_2level
    "SBM AMR requires a ConstantDz mesh without terrain"
    "-DOVERRIDE=erf.grid_stretching_ratio=1.1 erf.initial_dz=0.025")
  # (ERF itself makes levels of different compressibility one-way coupled.)
  add_sbm_expected_failure_test(SBM_AMR_Reject_MixedIntegrators
    inputs_amr_production_2level
    "SBM AMR requires the same compressibility and erf.timestepping_type on every level"
    "-DOVERRIDE=erf.compressibility_type=Anelastic erf.timestepping_type=RK2 MidPoint")
  add_sbm_expected_failure_test(SBM_AMR_Reject_AnelasticStratified
    inputs_amr_production_2level
    "SBM AMR with anelastic levels requires a constant base density"
    "-DOVERRIDE=erf.init_type=Isentropic erf.compressibility_type=Anelastic erf.timestepping_type=MidPoint")
  add_test(NAME SBM_M3_Adaptive_Anelastic_Donor_CFL
    COMMAND ${CMAKE_COMMAND}
      -DERF_EXECUTABLE=$<TARGET_FILE:erf_exec>
      -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_adaptive_anelastic_donor
      -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}
      "-DEXPECTED_STDOUT=SBM transport host integrator: AnelasticHeun"
      -DEXPECT_SECOND_CHECKPOINT=true
      -DSKIP_RESTART=true
      -DEXTRA_SCALAR_COMPONENTS=${_sbm_extra_scalars}
      -P ${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/RunM3Smoke.cmake)
  set_tests_properties(SBM_M3_Adaptive_Anelastic_Donor_CFL PROPERTIES
    LABELS "integration;sbm"
    RUN_SERIAL TRUE
    TIMEOUT 120)
  # The same deck with ERF's default anelastic integrator, MidPoint
  add_test(NAME SBM_M3_Adaptive_Anelastic_MidPoint_Donor_CFL
    COMMAND ${CMAKE_COMMAND}
      -DERF_EXECUTABLE=$<TARGET_FILE:erf_exec>
      -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_adaptive_anelastic_donor
      -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}
      -DRUN_ARGS=erf.timestepping_type=MidPoint
      "-DEXPECTED_STDOUT=SBM transport host integrator: AnelasticMidPoint"
      -DEXPECT_SECOND_CHECKPOINT=true
      -DSKIP_RESTART=true
      -DEXTRA_SCALAR_COMPONENTS=${_sbm_extra_scalars}
      -P ${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/RunM3Smoke.cmake)
  set_tests_properties(SBM_M3_Adaptive_Anelastic_MidPoint_Donor_CFL PROPERTIES
    LABELS "integration;sbm"
    RUN_SERIAL TRUE
    TIMEOUT 120)
  add_sbm_expected_failure_test(SBM_M3_Fixed_Unsafe_Donor_DT
    inputs_adaptive_anelastic_donor
    "SBM M3 fixed timestep exceeds donor positivity limit"
    "-DOVERRIDE=erf.fixed_dt=0.2")
  add_sbm_expected_failure_test(SBM_M3_Reject_PseudoIncompressible
    inputs_adaptive_anelastic_donor
    "SBM M3 does not yet support erf.compressibility_type = PseudoIncompressible"
    "-DOVERRIDE=erf.compressibility_type=PseudoIncompressible")

  add_test(NAME SBM_M3_Adaptive_Anelastic_Stretched_Donor_CFL
    COMMAND ${CMAKE_COMMAND}
      -DERF_EXECUTABLE=$<TARGET_FILE:erf_exec>
      -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_adaptive_anelastic_stretched_donor
      -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}
      "-DEXPECTED_STDOUT=SBM transport host integrator: AnelasticHeun"
      -DEXPECT_SECOND_CHECKPOINT=true
      -DSKIP_RESTART=true
      -DEXTRA_SCALAR_COMPONENTS=${_sbm_extra_scalars}
      -P ${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/RunM3Smoke.cmake)
  set_tests_properties(SBM_M3_Adaptive_Anelastic_Stretched_Donor_CFL PROPERTIES
    LABELS "integration;sbm"
    RUN_SERIAL TRUE
    TIMEOUT 120)
  # The stretched deck with MidPoint: non-uniform cell measures and the dt/2 first stage together
  add_test(NAME SBM_M3_Adaptive_Anelastic_Stretched_MidPoint_Donor_CFL
    COMMAND ${CMAKE_COMMAND}
      -DERF_EXECUTABLE=$<TARGET_FILE:erf_exec>
      -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_adaptive_anelastic_stretched_donor
      -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}
      -DRUN_ARGS=erf.timestepping_type=MidPoint
      "-DEXPECTED_STDOUT=SBM transport host integrator: AnelasticMidPoint"
      -DEXPECT_SECOND_CHECKPOINT=true
      -DSKIP_RESTART=true
      -DEXTRA_SCALAR_COMPONENTS=${_sbm_extra_scalars}
      -P ${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/RunM3Smoke.cmake)
  set_tests_properties(SBM_M3_Adaptive_Anelastic_Stretched_MidPoint_Donor_CFL PROPERTIES
    LABELS "integration;sbm"
    RUN_SERIAL TRUE
    TIMEOUT 120)
  add_sbm_expected_failure_test(SBM_M3_Fixed_Unsafe_Anelastic_Stretched_Donor_DT
    inputs_adaptive_anelastic_stretched_donor
    "SBM M3 fixed timestep exceeds donor positivity limit"
    "-DOVERRIDE=erf.fixed_dt=0.08")
  # The limit on this deck is 1/14.4 s = 0.0694 s, from the outflow rate U/dx + V/dy + W/dz of the
  # 0.1 thick bottom cell: the estimate must see the same mapped faces as the transport
  add_test(NAME SBM_M3_Fixed_Safe_Anelastic_Stretched_Donor_DT
    COMMAND $<TARGET_FILE:erf_exec>
      ${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_adaptive_anelastic_stretched_donor
      erf.fixed_dt=0.068 erf.check_int=-1 amrex.call_addr2line=0
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(SBM_M3_Fixed_Safe_Anelastic_Stretched_Donor_DT PROPERTIES
    LABELS "integration;sbm"
    RUN_SERIAL TRUE
    TIMEOUT 120)
  add_sbm_expected_failure_test(SBM_M3_Reject_Legacy_Zero_Transport inputs_zero_carrier
    "erf.sbm_zero_transport_fixture=true is obsolete"
    "-DOVERRIDE=erf.sbm_zero_transport_fixture=true")
  add_sbm_expected_failure_test(SBM_M3_Reject_Anelastic_StaticTerrain
    inputs_periodic_static_terrain
    "SBM M3 does not currently support anelastic StaticFittedMesh"
    "-DOVERRIDE=erf.compressibility_type=Anelastic erf.timestepping_type=RK2")
  add_sbm_expected_failure_test(SBM_Native_QC_Write_Fault inputs_zero_carrier
    "SBM owns compact qc/qr; native advection write denied for component range"
    "-DOVERRIDE=erf.sbm_test_native_qc_write_fault=true")
  add_sbm_expected_failure_test(SBM_Invalid_Nonmonotonic_Edges inputs_invalid_edges_nonmonotonic
    "Invalid SBM configuration: edges must be strictly increasing")
  add_sbm_expected_failure_test(SBM_Invalid_Zero_Edge_AutoPivot inputs_invalid_zero_edge_autopivot
    "Invalid SBM configuration: pivots must be finite and positive")
  add_sbm_expected_failure_test(SBM_Implicit_Substepping inputs_zero_carrier
    "SBM M3 does not support acoustic substepping"
    "-DOVERRIDE=erf.timestepping_type=RK3_ImplicitSubstep")
  add_sbm_expected_failure_test(SBM_Immersed_Buildings inputs_zero_carrier
    "SBM M3 requires buildings_type=None"
    "-DOVERRIDE=erf.buildings_type=ImmersedForcing")
  add_sbm_expected_failure_test(SBM_M3_Moving_Terrain inputs_zero_carrier
    "SBM M3 requires static non-EB geometry"
    "-DOVERRIDE=erf.terrain_type=MovingFittedMesh")

  # Physical faces: walls carry nothing, outflow loses only what its ledger
  # records, inflow brings a complete spectrum. Every run checks its leaf
  # inventory against the boundary ledger after each coarse step.
  function(add_sbm_boundary_budget_test test_name input_file expect_transfer)
    add_test(NAME ${test_name}
      COMMAND ${CMAKE_COMMAND}
        -DERF_EXECUTABLE=$<TARGET_FILE:erf_exec>
        -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/${input_file}
        -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}
        -DEXPECT_TRANSFER=${expect_transfer}
        ${ARGN}
        -P ${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/RunBoundaryBudget.cmake)
    set_tests_properties(${test_name} PROPERTIES
      LABELS "integration;sbm"
      RUN_SERIAL TRUE
      TIMEOUT 300)
  endfunction()
  set(_sbm_boundary_2m
    "erf.sbm_moment_mode=2|erf.sbm_edges=1.0e-18 2.0e-18 4.0e-18 8.0e-18 1.6e-17|erf.sbm_fixture_initial_state=1.5e-6 6.0e-6 1.8e-5 4.8e-5 1.0e12 2.0e12 3.0e12 4.0e12")
  add_sbm_boundary_budget_test(SBM_Boundary_Walls_SlipCompressible inputs_boundary_walls zero)
  add_sbm_boundary_budget_test(SBM_Boundary_Walls_NoSlipCompressible inputs_boundary_walls zero
    "-DRUN_ARGS=zlo.type=NoSlipWall|zhi.type=NoSlipWall")
  add_sbm_boundary_budget_test(SBM_Boundary_Walls_AnelasticMidPoint inputs_boundary_walls zero
    "-DRUN_ARGS=erf.compressibility_type=Anelastic|erf.timestepping_type=MidPoint")
  add_sbm_boundary_budget_test(SBM_Boundary_InflowOutflow_1M_Compressible
    inputs_boundary_inflow_outflow nonzero)
  add_sbm_boundary_budget_test(SBM_Boundary_InflowOutflow_2M_Compressible
    inputs_boundary_inflow_outflow nonzero
    "-DRUN_ARGS=${_sbm_boundary_2m}|xlo.sbm_inflow_z=1.2e-6 3.0e-6 3.6e-5 2.4e-5 0.8e12 1.0e12 6.0e12 2.0e12")
  add_sbm_boundary_budget_test(SBM_Boundary_InflowOutflow_AnelasticHeun
    inputs_boundary_inflow_outflow nonzero
    "-DRUN_ARGS=erf.compressibility_type=Anelastic|erf.timestepping_type=RK2")
  add_sbm_boundary_budget_test(SBM_Boundary_InflowOutflow_AnelasticMidPoint
    inputs_boundary_inflow_outflow nonzero
    "-DRUN_ARGS=erf.compressibility_type=Anelastic|erf.timestepping_type=MidPoint")
  add_sbm_boundary_budget_test(SBM_Boundary_InflowOutflow_HostTypes
    inputs_boundary_inflow_outflow nonzero
    "-DRUN_ARGS=xlo.type=inflow_outflow|xhi.type=ho_outflow")
  add_sbm_boundary_budget_test(SBM_Boundary_OutflowOnly inputs_boundary_inflow_outflow nonzero
    "-DRUN_ARGS=xlo.sbm_inflow_z=0 0 0 0")
  add_sbm_boundary_budget_test(SBM_Boundary_InflowOutflow_Stretched
    inputs_boundary_inflow_outflow nonzero
    "-DRUN_ARGS=erf.grid_stretching_ratio=1.1|erf.initial_dz=0.25")
  add_sbm_boundary_budget_test(SBM_Boundary_InflowOutflow_StaticTerrain
    inputs_boundary_static_terrain nonzero)
  # A front entering an empty domain binds the limiter, which drains cells to
  # zero; without its roundoff margin a drained cell rounds below zero.
  add_sbm_boundary_budget_test(SBM_Boundary_InflowFront_LimiterBinds
    inputs_boundary_inflow_outflow nonzero
    "-DRUN_ARGS=erf.sbm_fixture_initial_state=0 0 0 0|max_step=8")
  # The qualification build's executable accepts only its own two-level decks.
  if(NOT ERF_ENABLE_SBM_AMR_QUALIFICATION)
    add_sbm_boundary_budget_test(SBM_Boundary_AMR_Outflow_Compressible
      inputs_boundary_amr_2level nonzero -DEXPECTED_LEVELS=2)
    add_sbm_boundary_budget_test(SBM_Boundary_AMR_Outflow_AnelasticMidPoint
      inputs_boundary_amr_2level nonzero -DEXPECTED_LEVELS=2
      "-DRUN_ARGS=erf.compressibility_type=Anelastic|erf.timestepping_type=MidPoint")
    add_sbm_boundary_budget_test(SBM_Boundary_AMR_Inflow_Compressible
      inputs_boundary_amr_2level nonzero -DEXPECTED_LEVELS=2
      "-DRUN_ARGS=erf.box1.in_box_lo=0.0 0.25 0.0|erf.box1.in_box_hi=0.375 0.75 0.375|erf.box1.move_speed_x=0")
  endif()

  add_sbm_expected_failure_test(SBM_Boundary_Reject_ReverseOutflow
    inputs_boundary_inflow_outflow
    "SBM advective outflow xhi received an inward dry-air carrier"
    "-DOVERRIDE=prob.U_0=-0.8 xlo.velocity=-0.8 0.6 0.0")
  add_sbm_expected_failure_test(SBM_Boundary_Reject_MissingInflow
    inputs_boundary_inflow_outflow
    "SBM prescribed inflow at xhi requires xhi.sbm_inflow_z"
    "-DOVERRIDE=xhi.type=Inflow xhi.velocity=0.8 0.6 0.0 xhi.density=1.0 xhi.theta=300.0")
  add_sbm_expected_failure_test(SBM_Boundary_Reject_CompactInflow
    inputs_boundary_inflow_outflow
    "SBM owns compact qc/qr: xlo.qc cannot define or override a spectral inflow"
    "-DOVERRIDE=xlo.qc=1.e-6")
  add_sbm_expected_failure_test(SBM_Boundary_Reject_IncompleteInflow
    inputs_boundary_inflow_outflow
    "SBM inflow spectrum xlo.sbm_inflow_z is not admissible: it has 2 values"
    "-DOVERRIDE=xlo.sbm_inflow_z=4.e-6 3.e-6")
  add_sbm_expected_failure_test(SBM_Boundary_Reject_IncompleteTwoMomentGroup
    inputs_boundary_inflow_outflow
    "SBM inflow spectrum xlo.sbm_inflow_z is not admissible: population 0 bin 0 is not canonical"
    "-DOVERRIDE=erf.sbm_moment_mode=2 erf.sbm_edges=1.0e-18 2.0e-18 4.0e-18 8.0e-18 1.6e-17 erf.sbm_fixture_initial_state=1.5e-6 6.0e-6 1.8e-5 4.8e-5 1.0e12 2.0e12 3.0e12 4.0e12 xlo.sbm_inflow_z=1.2e-6 3.0e-6 3.6e-5 2.4e-5 0 1.0e12 6.0e12 2.0e12")
  add_sbm_expected_failure_test(SBM_Boundary_Reject_InflowOnWall
    inputs_boundary_inflow_outflow
    "zlo.sbm_inflow_z is set but zlo is not an inflow boundary"
    "-DOVERRIDE=zlo.sbm_inflow_z=4.e-6 3.e-6 2.e-6 1.e-6")
  add_sbm_expected_failure_test(SBM_Boundary_Reject_Open
    inputs_boundary_inflow_outflow
    "SBM physical boundary xhi is not supported: open boundaries"
    "-DOVERRIDE=xhi.type=Open")
  add_sbm_expected_failure_test(SBM_Boundary_Reject_Symmetry
    inputs_boundary_inflow_outflow
    "SBM physical boundary zlo is not supported: ERF does not hold the normal dry-air carrier"
    "-DOVERRIDE=zlo.type=Symmetry zhi.type=Symmetry")
  add_sbm_expected_failure_test(SBM_Boundary_Reject_BoundaryPlanes
    inputs_boundary_inflow_outflow
    "SBM does not support erf.input_bndry_planes"
    "-DOVERRIDE=erf.input_bndry_planes=1")
  add_sbm_expected_failure_test(SBM_Budget_Reject_WithoutSBM
    inputs_boundary_walls
    "erf.sbm_budget_check requires erf.moisture_model = SBM"
    "-DOVERRIDE=erf.moisture_model=None")

  # Particle mixing: ERF's LES eddy diffusivity, or a constant reference one,
  # mixes the spectrum through the combined advection-diffusion FCT on a single
  # level with periodic, wall and surface-layer faces. Walls carry nothing and
  # every coarse step closes the inventory budget.
  set(_sbm_mixing_host "-DEXPECTED_STDOUT=SBM particle mixing: particle-mixing=host-passive-scalar-v1")
  set(_sbm_mixing_constant
    "erf.les_type=None|erf.sbm_particle_mixing=ConstantReference|erf.sbm_particle_mixing_kh=2.e-2|erf.sbm_particle_mixing_kv=1.e-2")
  add_sbm_boundary_budget_test(SBM_Mixing_Walls_Smagorinsky inputs_mixing_walls zero
    ${_sbm_mixing_host})
  add_sbm_boundary_budget_test(SBM_Mixing_Walls_Deardorff inputs_mixing_walls zero
    ${_sbm_mixing_host} "-DRUN_ARGS=erf.les_type=Deardorff|erf.tke_min=1.e-2")
  add_sbm_boundary_budget_test(SBM_Mixing_Walls_NoSlipAnelasticMidPoint inputs_mixing_walls zero
    ${_sbm_mixing_host}
    "-DRUN_ARGS=zlo.type=NoSlipWall|zhi.type=NoSlipWall|erf.compressibility_type=Anelastic|erf.timestepping_type=MidPoint")
  add_sbm_boundary_budget_test(SBM_Mixing_SurfaceLayer_Smagorinsky inputs_mixing_walls zero
    ${_sbm_mixing_host} "-DRUN_ARGS=zlo.type=surface_layer|erf.most.z0=0.01")
  add_sbm_boundary_budget_test(SBM_Mixing_Stretched_ConstantReference inputs_mixing_walls zero
    "-DEXPECTED_STDOUT=SBM particle mixing: particle-mixing=constant-reference-v1"
    "-DRUN_ARGS=${_sbm_mixing_constant}|erf.grid_stretching_ratio=1.1|erf.initial_dz=0.08")
  add_sbm_boundary_budget_test(SBM_Mixing_Walls_2M inputs_mixing_walls zero
    ${_sbm_mixing_host}
    "-DRUN_ARGS=${_sbm_boundary_2m}|erf.sbm_fixture_initial_state_xhi=3.0e-6 1.2e-5 9.0e-6 2.4e-5 2.0e12 4.0e12 1.5e12 2.0e12")
  # An empty half carried by the flow binds the limiter with mixing active.
  add_sbm_boundary_budget_test(SBM_Mixing_EmptyHalf_LimiterBinds inputs_mixing_walls zero
    ${_sbm_mixing_host} "-DRUN_ARGS=erf.sbm_fixture_initial_state_xhi=0 0 0 0|prob.U_0=0.5|max_step=8")
  # Adaptive steps: the lagged eddy estimate, and a reference diffusivity whose
  # combined limit (6.9e-5 s) is below the acoustic one (1.2e-4 s).
  add_sbm_boundary_budget_test(SBM_Mixing_Adaptive_Smagorinsky inputs_mixing_walls zero
    ${_sbm_mixing_host} "-DRUN_ARGS=erf.fixed_dt=-1|erf.cfl=0.8|max_step=6")
  set(_sbm_mixing_constant_adaptive
    "erf.fixed_dt=-1|erf.cfl=0.8|erf.les_type=None|erf.sbm_particle_mixing=ConstantReference|erf.sbm_particle_mixing_kh=40|erf.sbm_particle_mixing_kv=40")
  add_sbm_boundary_budget_test(SBM_Mixing_Adaptive_ConstantReference_DiffusionBound
    inputs_mixing_walls zero "-DRUN_ARGS=${_sbm_mixing_constant_adaptive}|max_step=6")

  add_sbm_expected_failure_test(SBM_Mixing_Reject_PBL inputs_mixing_walls
    "SBM does not support a PBL closure .erf.pbl_type = MYNN25 .level 0.."
    "-DOVERRIDE=erf.les_type=None erf.sbm_particle_mixing=None erf.pbl_type=MYNN25 zlo.type=surface_layer")
  add_sbm_expected_failure_test(SBM_Mixing_Reject_Smagorinsky2D inputs_mixing_walls
    "SBM particle mixing HostPassiveScalar does not support erf.les_type = Smagorinsky2D"
    "-DOVERRIDE=erf.les_type=Smagorinsky2D")
  add_sbm_expected_failure_test(SBM_Mixing_Reject_LESWithoutMixing inputs_mixing_walls
    "SBM particle mixing None does not mix the spectrum with ERF's LES closure"
    "-DOVERRIDE=erf.sbm_particle_mixing=None")
  add_sbm_expected_failure_test(SBM_Mixing_Reject_LESWithConstantReference inputs_mixing_walls
    "SBM particle mixing ConstantReference does not mix the spectrum with ERF's LES closure"
    "-DOVERRIDE=erf.sbm_particle_mixing=ConstantReference erf.sbm_particle_mixing_kh=1.e-2")
  add_sbm_expected_failure_test(SBM_Mixing_Reject_UnknownPolicy inputs_mixing_walls
    "erf.sbm_particle_mixing must be None, ConstantReference or HostPassiveScalar, not 'Eddy'"
    "-DOVERRIDE=erf.sbm_particle_mixing=Eddy")
  add_sbm_expected_failure_test(SBM_Mixing_Reject_HostWithReference inputs_mixing_walls
    "reference diffusivities apply only to ConstantReference, not HostPassiveScalar"
    "-DOVERRIDE=erf.sbm_particle_mixing_kv=1.e-2")
  add_sbm_expected_failure_test(SBM_Mixing_Reject_NegativeReference inputs_mixing_walls
    "ConstantReference diffusivities must be finite and nonnegative"
    "-DOVERRIDE=erf.les_type=None erf.sbm_particle_mixing=ConstantReference erf.sbm_particle_mixing_kh=-1.e-2")
  add_sbm_expected_failure_test(SBM_Mixing_Reject_StaticTerrain inputs_periodic_static_terrain
    "SBM particle mixing has no qualified monotone low-order operator on a terrain-fitted"
    "-DOVERRIDE=erf.les_type=Smagorinsky erf.Cs=0.2 erf.sbm_particle_mixing=HostPassiveScalar")
  # The qualification build's executable accepts only its own two-level decks.
  if(NOT ERF_ENABLE_SBM_AMR_QUALIFICATION)
    add_sbm_expected_failure_test(SBM_Mixing_Reject_AMR inputs_amr_production_2level
      "SBM particle mixing HostPassiveScalar is not yet qualified with AMR"
      "-DOVERRIDE=erf.les_type=Smagorinsky erf.Cs=0.2 erf.sbm_particle_mixing=HostPassiveScalar")
  endif()
  add_sbm_expected_failure_test(SBM_Mixing_Reject_InflowOutflow inputs_boundary_inflow_outflow
    "SBM particle mixing HostPassiveScalar is qualified with periodic and wall faces only, not the"
    "-DOVERRIDE=erf.les_type=Smagorinsky erf.Cs=0.2 erf.sbm_particle_mixing=HostPassiveScalar")
  add_sbm_expected_failure_test(SBM_Mixing_Stage_Combined_CFL inputs_mixing_walls
    "SBM combined advection-diffusion outgoing-demand condition failed"
    "-DOVERRIDE=erf.sbm_test_donor_cfl_fault=true")
  # dx = 1/8 and K = 100 m^2/s: the diffusive outflow alone allows about 2.6e-5 s
  add_sbm_expected_failure_test(SBM_Mixing_Fixed_Unsafe_DT inputs_mixing_walls
    "SBM fixed timestep exceeds the combined advection-diffusion positivity limit"
    "-DOVERRIDE=erf.les_type=None erf.sbm_particle_mixing=ConstantReference erf.sbm_particle_mixing_kh=100 erf.sbm_particle_mixing_kv=100")

  add_test(NAME SBM_M3_Stretched_Static_Mesh
    COMMAND ${CMAKE_COMMAND}
      -DERF_EXECUTABLE=$<TARGET_FILE:erf_exec>
      -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_stretched
      -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}
      -DEXTRA_SCALAR_COMPONENTS=${_sbm_extra_scalars}
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
      -DEXTRA_SCALAR_COMPONENTS=${_sbm_extra_scalars}
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

    # The same across a restart that re-makes the level-0 grids (issue 4225): the spectrum
    # is read on the checkpoint's 2x2 grids and redistributed onto the restart's 4x4, and
    # the run must continue as the straight run did.
    add_test(NAME SBM_M3_StaticTerrain_RegridOnRestart
      COMMAND ${CMAKE_COMMAND}
        -DERF_EXECUTABLE=$<TARGET_FILE:erf_exec>
        -DCOMPARE_EXECUTABLE=$<TARGET_FILE:erf_sbm_checkpoint_compare>
        -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_periodic_static_terrain
        -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}
        "-DRESTART_ARGS=erf.regrid_level_0_on_restart=1|amr.max_grid_size_x=2|amr.max_grid_size_y=2"
        -DALLOW_DIFF_GRIDS=ON
        -P ${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/RunM3RestartEquivalence.cmake)
    set_tests_properties(SBM_M3_StaticTerrain_RegridOnRestart PROPERTIES
      LABELS "integration;sbm"
      RUN_SERIAL TRUE
      TIMEOUT 180)

    # Production two-level SBM: a refined box that moves one coarse cell per
    # step is regridded every step, and a restart from the middle of the run
    # reproduces the continuous run bitwise on both levels.  The qualification
    # build's executable accepts only its own decks.
    if(NOT ERF_ENABLE_SBM_AMR_QUALIFICATION)
      foreach(_sbm_host IN ITEMS Compressible AnelasticMidPoint)
        set(_sbm_host_args)
        if(_sbm_host STREQUAL "AnelasticMidPoint")
          set(_sbm_host_args -DANELASTIC=MidPoint)
        endif()
        add_test(NAME SBM_AMR_Production_Regrid_Restart_${_sbm_host}
          COMMAND ${CMAKE_COMMAND}
            -DERF_EXECUTABLE=$<TARGET_FILE:erf_exec>
            -DCOMPARE_EXECUTABLE=$<TARGET_FILE:erf_sbm_checkpoint_compare>
            -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_amr_production_2level
            -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}/sbm_amr_production_${_sbm_host}
            -DPRODUCTION=ON
            -DMODE=equivalence
            -DRESTART_STEP=2
            -DFINAL_STEP=4
            -DEXPECTED_FINAL_LEVELS=2
            ${_sbm_host_args}
            -P ${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/RunAMRRestart.cmake)
        set_tests_properties(SBM_AMR_Production_Regrid_Restart_${_sbm_host} PROPERTIES
          LABELS "integration;sbm"
          RUN_SERIAL TRUE
          TIMEOUT 300)
      endforeach()
    endif()

    # A restart carries the boundary ledger and re-reads the inflow it was
    # fingerprinted with; a changed inflow is a different schema.
    function(add_sbm_boundary_restart_test test_name input_file mode)
      add_test(NAME ${test_name}
        COMMAND ${CMAKE_COMMAND}
          -DERF_EXECUTABLE=$<TARGET_FILE:erf_exec>
          -DCOMPARE_EXECUTABLE=$<TARGET_FILE:erf_sbm_checkpoint_compare>
          -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/${input_file}
          -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}/sbm_boundary_restart_${test_name}
          -DPRODUCTION=ON
          -DMODE=${mode}
          -DRESTART_STEP=2
          -DFINAL_STEP=4
          -DEXPECT_BOUNDARY_LEDGER=ON
          ${ARGN}
          -P ${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/RunAMRRestart.cmake)
      set_tests_properties(${test_name} PROPERTIES
        LABELS "integration;sbm"
        RUN_SERIAL TRUE
        TIMEOUT 300)
    endfunction()
    add_sbm_boundary_restart_test(SBM_Boundary_Restart_Equivalence
      inputs_boundary_inflow_outflow equivalence -DEXPECTED_FINAL_LEVELS=1)
    add_sbm_boundary_restart_test(SBM_Boundary_Restart_ChangedInflow
      inputs_boundary_inflow_outflow changed_restart_inputs
      "-DRESTART_ARGS=xlo.sbm_inflow_z=4.e-6 3.e-6 2.e-6 2.e-6"
      "-DEXPECTED_DIAGNOSTIC=SBM restart schema does not exactly match")
    if(NOT ERF_ENABLE_SBM_AMR_QUALIFICATION)
      add_sbm_boundary_restart_test(SBM_Boundary_AMR_Restart_Equivalence_Compressible
        inputs_boundary_amr_2level equivalence -DEXPECTED_FINAL_LEVELS=2)
      add_sbm_boundary_restart_test(SBM_Boundary_AMR_Restart_Equivalence_AnelasticMidPoint
        inputs_boundary_amr_2level equivalence -DEXPECTED_FINAL_LEVELS=2 -DANELASTIC=MidPoint)
    endif()

    # A restart reproduces the uninterrupted run bitwise, and a restart under
    # another mixing policy is a different schema. Only ConstantReference has a
    # clean reference without mixing (the same flow), so only it checks that
    # mixing changes the spectrum; turning LES off would also change the flow.
    function(add_sbm_particle_mixing_test test_name)
      add_test(NAME ${test_name}
        COMMAND ${CMAKE_COMMAND}
          -DERF_EXECUTABLE=$<TARGET_FILE:erf_exec>
          -DCOMPARE_EXECUTABLE=$<TARGET_FILE:erf_sbm_checkpoint_compare>
          -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_mixing_walls
          -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}
          -DFINAL_STEP=4
          ${ARGN}
          -P ${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/RunParticleMixing.cmake)
      set_tests_properties(${test_name} PROPERTIES
        LABELS "integration;sbm"
        RUN_SERIAL TRUE
        TIMEOUT 300)
    endfunction()
    add_sbm_particle_mixing_test(SBM_Mixing_Smagorinsky_Restart
      -DRESTART_STEP=2
      "-DMISMATCH_ARGS=erf.les_type=None|erf.sbm_particle_mixing=None")
    # At rest Deardorff's eddy diffusivity comes from tke_min alone and the flow
    # stays at rest, so the run without LES differs from it only by the mixing.
    add_sbm_particle_mixing_test(SBM_Mixing_Deardorff_AtRest_Effect_Restart
      "-DRUN_ARGS=erf.les_type=Deardorff|erf.tke_min=1.e-2|prob.test_flow_amp=0"
      "-DREFERENCE_ARGS=erf.les_type=None|erf.sbm_particle_mixing=None"
      -DRESTART_STEP=2)
    add_sbm_particle_mixing_test(SBM_Mixing_ConstantReference_Effect_Restart
      "-DRUN_ARGS=${_sbm_mixing_constant}"
      "-DREFERENCE_ARGS=erf.sbm_particle_mixing=None|erf.sbm_particle_mixing_kh=0|erf.sbm_particle_mixing_kv=0"
      -DRESTART_STEP=2
      "-DMISMATCH_ARGS=erf.sbm_particle_mixing_kh=3.e-2")
    # Reference coefficients need no restart cap, so an adaptive restart is bitwise.
    add_sbm_particle_mixing_test(SBM_Mixing_ConstantReference_Adaptive_Restart
      "-DRUN_ARGS=${_sbm_mixing_constant_adaptive}"
      -DRESTART_STEP=2)
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
        -DEXTRA_SCALAR_COMPONENTS=${_sbm_extra_scalars}
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

if(TARGET erf_exec AND ERF_DIM STREQUAL "3" AND
   ERF_PRECISION STREQUAL "SINGLE" AND NOT AMReX_HDF5 AND
   NOT ERF_ENABLE_CUDA AND NOT ERF_ENABLE_HIP AND NOT ERF_ENABLE_SYCL)
  add_executable(erf_sbm_checkpoint_compare
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SBM/ERF_SBMCheckpointCompare.cpp)
  target_link_libraries(erf_sbm_checkpoint_compare PRIVATE AMReX::amrex)
  target_include_directories(erf_sbm_checkpoint_compare PRIVATE ${PROJECT_SOURCE_DIR}/Source)
  target_compile_features(erf_sbm_checkpoint_compare PRIVATE cxx_std_17)

  add_test(NAME SBM_FLOAT_FractionalTimeRestart
    COMMAND ${CMAKE_COMMAND}
      -DERF_EXECUTABLE=$<TARGET_FILE:erf_exec>
      -DCOMPARE_EXECUTABLE=$<TARGET_FILE:erf_sbm_checkpoint_compare>
      -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_periodic_advection_1m
      -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}
      -P ${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/RunFloatRestartTime.cmake)
  set_tests_properties(SBM_FLOAT_FractionalTimeRestart PROPERTIES
    LABELS "integration;sbm;single"
    RUN_SERIAL TRUE
    TIMEOUT 180)
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
      ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SBM/ERF_GTestSBMAMRTransferParallel.cpp
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
      ${CMAKE_CURRENT_SOURCE_DIR}/Radiation/ERF_GTestTopographicShading.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestCloudOptics.cpp
      # The immersed average down copies the fine sums to the coarse layout: check it across ranks
      ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ERF_GTestImmersedAverageDown.cpp
      # The mismatch count reduces over the ranks: check it across ranks
      ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ERF_GTestImmersedCFFaces.cpp
  )
  if(ERF_ENABLE_NETCDF)
    # The terrain reader opens the fixture only on the I/O rank and then
    # broadcasts the result.  Keep this source in the MPI executable so the
    # WPS edge/corner regression exercises every rank.
    target_sources(${erf_parallel_exe_name} PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestTerrainNetCDF.cpp
    )
  endif()
  if(ERF_ENABLE_MOVING_BODIES)
    # Motivation: with the domain distributed over ranks, each actuator point is sampled by
    # the one rank owning its cell and summed across ranks; the result must not depend on
    # the rank count or on which rank holds the point.
    target_sources(${erf_parallel_exe_name} PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/MovingBodies/ERF_GTestActuatorSampling.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/MovingBodies/ERF_GTestActuatorSpreading.cpp
    )
  endif()
  if(ERF_DIM STREQUAL "3")
    # The coefficient check and the rate maximum reduce over the ranks: check them across ranks
    target_sources(${erf_parallel_exe_name} PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SBM/ERF_GTestSBMParticleMixing.cpp
      # Mixing reads neighbours across box and rank boundaries in every stage
      ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SBM/ERF_GTestSBMMixingTransport.cpp
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

if(TARGET erf_exec AND ERF_DIM STREQUAL "3" AND NOT AMReX_HDF5)
  add_executable(erf_cloud_optics_plotfile_check
    ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_CloudOpticsPlotfileCheck.cpp)
  target_link_libraries(erf_cloud_optics_plotfile_check PRIVATE AMReX::amrex)
  if(ERF_ENABLE_CUDA)
    setup_target_for_cuda_compilation(erf_cloud_optics_plotfile_check)
  endif()
  add_test(NAME CloudOptics_RadiationOff_Plot2D
    COMMAND ${CMAKE_COMMAND}
      -DERF_EXECUTABLE=$<TARGET_FILE:erf_exec>
      -DCHECK_EXECUTABLE=$<TARGET_FILE:erf_cloud_optics_plotfile_check>
      -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Exec/RegTests/Bubble/inputs_BF02_moist_bubble
      -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}/cloud_optics_writer
      -P ${CMAKE_CURRENT_SOURCE_DIR}/IO/RunCloudOptics.cmake)
  set_tests_properties(CloudOptics_RadiationOff_Plot2D PROPERTIES
    LABELS "integration;cloud_optics" TIMEOUT 60 RUN_SERIAL TRUE)
endif()
