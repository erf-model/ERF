if(ERF_ENABLE_TESTS OR ERF_ENABLE_UNIT_TESTS)
  if(NOT DEFINED ERF_TEST_NRANKS OR "${ERF_TEST_NRANKS}" STREQUAL "")
    set(ERF_TEST_NRANKS 2 CACHE STRING "Number of MPI ranks to use for each regression test" FORCE)
  else()
    set(ERF_TEST_NRANKS "${ERF_TEST_NRANKS}" CACHE STRING "Number of MPI ranks to use for each regression test")
  endif()

  if(NOT DEFINED ERF_PARALLEL_TEST_NRANKS OR "${ERF_PARALLEL_TEST_NRANKS}" STREQUAL "")
    set(ERF_PARALLEL_TEST_NRANKS "2" CACHE STRING
        "Semicolon-separated MPI rank counts for parallel GoogleTest unit tests" FORCE)
  else()
    set(ERF_PARALLEL_TEST_NRANKS "${ERF_PARALLEL_TEST_NRANKS}" CACHE STRING
        "Semicolon-separated MPI rank counts for parallel GoogleTest unit tests")
  endif()

  option(ERF_PARALLEL_TEST_ALLOW_OVERSUBSCRIBE
    "Allow oversubscribed MPI launches for parallel GoogleTest unit tests"
    OFF)

  set(ERF_PARALLEL_TEST_MPIEXEC_PREFLAGS "" CACHE STRING
      "Extra mpiexec preflags for parallel GoogleTest unit tests")
endif()

if(ERF_ENABLE_TESTS)
  # Additional testing options
  set(ERF_TEST_FCOMPARE_RTOL "2.0e-10" CACHE STRING "fcompare relative tolerance")
  set(ERF_TEST_FCOMPARE_ATOL "2.0e-10" CACHE STRING "fcompare absolute tolerance")
  set(ERF_TEST_GOLD_FILES_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/ERFGoldFiles" CACHE PATH "Path to ERF gold files")
  set(ERF_TEST_EXTRA_FILES_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/test_files/" CACHE PATH "Path to additional input files")
  option(ERF_TEST_ENABLE_EXTRA_SDM_TESTS "Enable additional SDM tests that require machine-specific gold files" OFF)
  option(ERF_TEST_ENABLE_EXTRA_LSM_TESTS "Enable long-running LSM regression tests" OFF)
  option(ERF_TEST_SDM_SMOKE_GATED "Run the machine-specific SDM tests to completion, without gold comparison, under the smoke label" OFF)

  # The SEB multi-level parity check is a python script rather than another compiled
  # checker: it asserts an exact relation between two fextract slices, so it needs no
  # AMReX linkage of its own.
  find_package(Python3 COMPONENTS Interpreter QUIET)
  set(TWO_STREAM_SEB_PARITY_CHECKER
      "${CMAKE_CURRENT_SOURCE_DIR}/check_two_stream_seb_parity.py")
  set(TWO_STREAM_SEB_FLUX_SOURCE_CHECKER
      "${CMAKE_CURRENT_SOURCE_DIR}/check_two_stream_seb_flux_source.py")
  set(TWO_STREAM_SEB_MOISTURE_CHECKER
      "${CMAKE_CURRENT_SOURCE_DIR}/check_two_stream_seb_moisture.py")
  set(TWO_STREAM_NOAHMP_LEVELS_CHECKER
      "${CMAKE_CURRENT_SOURCE_DIR}/check_two_stream_noahmp_levels.py")
  if(Python3_Interpreter_FOUND)
    set(ERF_TEST_PYTHON "${Python3_EXECUTABLE}")
  else()
    set(ERF_TEST_PYTHON "")
  endif()

  message(STATUS "ERF testing configuration summary:")
  message(STATUS "   Number of ranks               = ${ERF_TEST_NRANKS}")
  message(STATUS "   fcompare executable           = ${FCOMPARE_EXE}")
  message(STATUS "   fextrema executable           = ${FEXTREMA_EXE}")
  message(STATUS "   fextract executable           = ${FEXTRACT_EXE}")
  message(STATUS "   test python interpreter       = ${ERF_TEST_PYTHON}")
  message(STATUS "   comparison relative tolerance = ${ERF_TEST_FCOMPARE_RTOL}")
  message(STATUS "   comparison absolute tolerance = ${ERF_TEST_FCOMPARE_ATOL}")
  message(STATUS "   Gold file directory           = ${ERF_TEST_GOLD_FILES_DIRECTORY}")
  message(STATUS "   Extra SDM tests enabled       = ${ERF_TEST_ENABLE_EXTRA_SDM_TESTS}")
  message(STATUS "   Extra LSM tests enabled       = ${ERF_TEST_ENABLE_EXTRA_LSM_TESTS}")
  message(STATUS "   Gated SDM smoke tests         = ${ERF_TEST_SDM_SMOKE_GATED}")

  add_executable(erf_shoc_plotfile_check
    ${CMAKE_CURRENT_SOURCE_DIR}/ShocPlotfileCheck.cpp)
  target_link_libraries(erf_shoc_plotfile_check PRIVATE AMReX::amrex)
  target_compile_features(erf_shoc_plotfile_check PRIVATE cxx_std_17)
  set(SHOC_PLOTFILE_CHECKER
      "${CMAKE_CURRENT_BINARY_DIR}/erf_shoc_plotfile_check${CMAKE_EXECUTABLE_SUFFIX}")

  add_executable(erf_shoc_microphysics_differential
    ${CMAKE_CURRENT_SOURCE_DIR}/ShocMicrophysicsDifferential.cpp)
  target_link_libraries(erf_shoc_microphysics_differential PRIVATE AMReX::amrex)
  target_compile_features(erf_shoc_microphysics_differential PRIVATE cxx_std_17)
  set(SHOC_MICROPHYSICS_DIFFERENTIAL
      "${CMAKE_CURRENT_BINARY_DIR}/erf_shoc_microphysics_differential${CMAKE_EXECUTABLE_SUFFIX}")

  add_executable(erf_shoc_mutation_differential
    ${CMAKE_CURRENT_SOURCE_DIR}/ShocMutationDifferential.cpp)
  target_link_libraries(erf_shoc_mutation_differential PRIVATE AMReX::amrex)
  target_compile_features(erf_shoc_mutation_differential PRIVATE cxx_std_17)
  set(SHOC_MUTATION_DIFFERENTIAL
      "${CMAKE_CURRENT_BINARY_DIR}/erf_shoc_mutation_differential${CMAKE_EXECUTABLE_SUFFIX}")

  add_executable(erf_shoc_gold_differential
    ${CMAKE_CURRENT_SOURCE_DIR}/ShocGoldDifferential.cpp)
  target_link_libraries(erf_shoc_gold_differential PRIVATE AMReX::amrex)
  target_compile_features(erf_shoc_gold_differential PRIVATE cxx_std_17)
  set(SHOC_GOLD_DIFFERENTIAL
      "${CMAKE_CURRENT_BINARY_DIR}/erf_shoc_gold_differential${CMAKE_EXECUTABLE_SUFFIX}")

  # Regression motivation: these three rotated manufactured chambers exercise
  # the x/y/z wall-diffusion paths, including the historical x-high branch,
  # with one independent stationary-state oracle.
  add_executable(erf_anelastic_wall_diffusion_check
    ${CMAKE_CURRENT_SOURCE_DIR}/AnelasticWallDiffusionCheck.cpp)
  target_link_libraries(erf_anelastic_wall_diffusion_check PRIVATE AMReX::amrex)
  target_compile_features(erf_anelastic_wall_diffusion_check PRIVATE cxx_std_17)
  set(ANELASTIC_WALL_DIFFUSION_CHECKER
      "${CMAKE_CURRENT_BINARY_DIR}/erf_anelastic_wall_diffusion_check${CMAKE_EXECUTABLE_SUFFIX}")

  # Motivation: the average down over immersed-forcing levels must give every covered coarse
  # cell the fluid-mass-weighted average of the fine cells over it; checked on the plotfile
  # instead of a gold file, with a leg that proves the plain average fails the check.
  add_executable(erf_immersed_average_down_check
    ${CMAKE_CURRENT_SOURCE_DIR}/ImmersedAverageDownCheck.cpp)
  target_link_libraries(erf_immersed_average_down_check PRIVATE AMReX::amrex)
  target_compile_features(erf_immersed_average_down_check PRIVATE cxx_std_17)

  # Motivation: a refined patch whose coarse/fine faces cross an immersed surface must leave no net
  # heat source; checked on the plotfile instead of a gold file, with a leg that proves the explicit
  # solid freeze fails the check.
  add_executable(erf_immersed_cf_heat_balance_check
    ${CMAKE_CURRENT_SOURCE_DIR}/ImmersedCFHeatBalanceCheck.cpp)
  target_link_libraries(erf_immersed_cf_heat_balance_check PRIVATE AMReX::amrex)
  target_compile_features(erf_immersed_cf_heat_balance_check PRIVATE cxx_std_17)

  # Motivation: the fraction-stress immersed wall law on a face-aligned surface must reproduce flat
  # ground; checked by comparing two runs' profiles, with a leg that proves the legacy law fails.
  add_executable(erf_immersed_flat_equivalence_check
    ${CMAKE_CURRENT_SOURCE_DIR}/ImmersedFlatEquivalenceCheck.cpp)
  target_link_libraries(erf_immersed_flat_equivalence_check PRIVATE AMReX::amrex)
  target_compile_features(erf_immersed_flat_equivalence_check PRIVATE cxx_std_17)

  # Motivation: the chamber regressions must prove the initializer is wired
  # into a real 3-D run and that the SatAdj deck forms finite cloud water.
  add_executable(erf_cloud_chamber_check
    ${CMAKE_CURRENT_SOURCE_DIR}/CloudChamberCheck.cpp)
  target_link_libraries(erf_cloud_chamber_check PRIVATE AMReX::amrex)
  target_include_directories(erf_cloud_chamber_check PRIVATE
    $<TARGET_PROPERTY:erf_srclib,INTERFACE_INCLUDE_DIRECTORIES>)
  target_compile_features(erf_cloud_chamber_check PRIVATE cxx_std_17)
  set(CLOUD_CHAMBER_CHECKER "$<TARGET_FILE:erf_cloud_chamber_check>")

  # Run the typed Cloud Chamber parser directly in a child process so its
  # allowlist abort can be tested without the earlier global solver guards.
  add_executable(erf_cloud_chamber_parser_check
    ${CMAKE_CURRENT_SOURCE_DIR}/CloudChamberParserCheck.cpp)
  target_link_libraries(erf_cloud_chamber_parser_check PRIVATE AMReX::amrex)
  target_include_directories(erf_cloud_chamber_parser_check PRIVATE
    $<TARGET_PROPERTY:erf_srclib,INTERFACE_INCLUDE_DIRECTORIES>)
  target_compile_features(erf_cloud_chamber_parser_check PRIVATE cxx_std_17)
  set(CLOUD_CHAMBER_PARSER_CHECKER "$<TARGET_FILE:erf_cloud_chamber_parser_check>")

  # Motivation: the TwoStream column sweep must keep the surface at k = 0 and
  # cool to space from the top layer; this checker reads qsrc_sw / qsrc_lw
  # from a plotfile and verifies that vertical structure without gold files.
  add_executable(erf_two_stream_radiation_check
    ${CMAKE_CURRENT_SOURCE_DIR}/TwoStreamRadiationCheck.cpp)
  target_link_libraries(erf_two_stream_radiation_check PRIVATE AMReX::amrex)
  target_compile_features(erf_two_stream_radiation_check PRIVATE cxx_std_17)
  set(TWO_STREAM_RADIATION_CHECKER
      "${CMAKE_CURRENT_BINARY_DIR}/erf_two_stream_radiation_check${CMAKE_EXECUTABLE_SUFFIX}")
  # Motivation: the bulk wall fixed-dt safety contract must fail loudly when
  # a requested fixed timestep exceeds the analytically known wall limit.
  add_executable(erf_cloud_chamber_wall_dt_guard_check
    ${CMAKE_CURRENT_SOURCE_DIR}/CloudChamberWallDtGuardCheck.cpp)
  target_link_libraries(erf_cloud_chamber_wall_dt_guard_check PRIVATE AMReX::amrex)
  target_include_directories(erf_cloud_chamber_wall_dt_guard_check PRIVATE
    $<TARGET_PROPERTY:erf_srclib,INTERFACE_INCLUDE_DIRECTORIES>)
  target_compile_features(erf_cloud_chamber_wall_dt_guard_check PRIVATE cxx_std_17)

  # Motivation: station output is compared against itself (the same point asked
  # for two ways) and, for observation nudging, against closed forms; this
  # checker reads the station time series ERF writes.  Plain C++, no AMReX; the
  # tests name it by $<TARGET_FILE:erf_station_series_check>.
  add_executable(erf_station_series_check
    ${CMAKE_CURRENT_SOURCE_DIR}/StationSeriesCheck.cpp)
  target_compile_features(erf_station_series_check PRIVATE cxx_std_17)

  if(ERF_ENABLE_CUDA)
    foreach(_erf_test_tool IN ITEMS
        erf_shoc_plotfile_check
        erf_shoc_microphysics_differential
        erf_shoc_mutation_differential
        erf_shoc_gold_differential
        erf_anelastic_wall_diffusion_check
        erf_cloud_chamber_check
        erf_cloud_chamber_parser_check
        erf_two_stream_radiation_check
        erf_cloud_chamber_wall_dt_guard_check
        erf_immersed_average_down_check
        erf_immersed_cf_heat_balance_check
        erf_immersed_flat_equivalence_check)
      setup_target_for_cuda_compilation(${_erf_test_tool})
    endforeach()
  endif()

  include(${CMAKE_CURRENT_SOURCE_DIR}/CTestList.cmake)

  # Keep the end-to-end two-level SBM trajectory reachable only from the
  # separately compiled qualification executable. The production executable
  # does not register or enable this internal AMR route.
  if(ERF_ENABLE_SBM_AMR_QUALIFICATION AND ERF_DIM STREQUAL "3" AND
     ERF_PRECISION STREQUAL "DOUBLE" AND ERF_ENABLE_MPI)
    if(NOT ERF_TEST_PYTHON)
      message(FATAL_ERROR
        "Split SBM qualification requires Python 3 for operation-aware decomposition parity")
    endif()
    resolve_test_exe("" "erf_exec" SBM_AMR_QUALIFICATION_EXECUTABLE)
    set(_sbm_amr_qualification_script
        "${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/RunAMRQualification.cmake")
    set(_sbm_amr_qualification_root
        "${CMAKE_CURRENT_BINARY_DIR}/test_files/SBM_AMR_Qualification")

    add_test(NAME SBM_AMR_Qualification_1M
      COMMAND ${CMAKE_COMMAND}
        -DERF_EXECUTABLE=${SBM_AMR_QUALIFICATION_EXECUTABLE}
        -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_amr_qualification_1m
        -DTEST_ROOT=${_sbm_amr_qualification_root}/1m
        -DMPIEXEC=${MPIEXEC_EXECUTABLE}
        -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
        -DMPI_RANKS=1
        -DEXPECTED_COMPONENTS=4
        -DEXPECTED_FINE_SUBSTEPS=3
        -P ${_sbm_amr_qualification_script})
    set_tests_properties(SBM_AMR_Qualification_1M PROPERTIES
      LABELS "integration;sbm;sbm_amr_qualification"
      RUN_SERIAL TRUE
      TIMEOUT 300)

    add_test(NAME SBM_AMR_Qualification_ZeroFlow_1M
      COMMAND ${CMAKE_COMMAND}
        -DERF_EXECUTABLE=${SBM_AMR_QUALIFICATION_EXECUTABLE}
        -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_amr_qualification_1m
        -DTEST_ROOT=${_sbm_amr_qualification_root}/zero_flow_1m
        -DMPIEXEC=${MPIEXEC_EXECUTABLE}
        -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
        -DMPI_RANKS=1
        -DEXPECTED_COMPONENTS=4
        -DEXPECTED_STEPS=1
        -DEXPECTED_FINE_SUBSTEPS=3
        -DMAX_STEPS=1
        -DZERO_FLOW=ON
        -P ${_sbm_amr_qualification_script})
    set_tests_properties(SBM_AMR_Qualification_ZeroFlow_1M PROPERTIES
      LABELS "integration;sbm;sbm_amr_qualification"
      RUN_SERIAL TRUE
      TIMEOUT 300)

    add_test(NAME SBM_AMR_Qualification_2M_Attached_MPI
      COMMAND ${CMAKE_COMMAND}
        -DERF_EXECUTABLE=${SBM_AMR_QUALIFICATION_EXECUTABLE}
        -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_amr_qualification_2m_attached_empty
        -DTEST_ROOT=${_sbm_amr_qualification_root}/2m_attached
        -DMPIEXEC=${MPIEXEC_EXECUTABLE}
        -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
        -DMPI_RANKS=2
        -DEXPECTED_COMPONENTS=12
        -DEXPECTED_FINE_SUBSTEPS=3
        -DEXPECT_EMPTY_FINE_RANK=ON
        -P ${_sbm_amr_qualification_script})
    set_tests_properties(SBM_AMR_Qualification_2M_Attached_MPI PROPERTIES
      LABELS "integration;sbm;sbm_amr_qualification;parallel"
      RUN_SERIAL TRUE
      PROCESSORS 2
      TIMEOUT 300)

    add_test(NAME SBM_AMR_Qualification_1M_Split_1Rank
      COMMAND ${CMAKE_COMMAND}
        -DERF_EXECUTABLE=${SBM_AMR_QUALIFICATION_EXECUTABLE}
        -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_amr_qualification_1m_split
        -DTEST_ROOT=${_sbm_amr_qualification_root}/1m_split_1rank
        -DMPIEXEC=${MPIEXEC_EXECUTABLE}
        -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
        -DMPI_RANKS=1
        -DEXPECTED_COMPONENTS=4
        -DEXPECTED_FINE_SUBSTEPS=3
        -DSPATIAL_SPLIT_FIXTURE=ON
        -P ${_sbm_amr_qualification_script})
    set_tests_properties(SBM_AMR_Qualification_1M_Split_1Rank PROPERTIES
      LABELS "integration;sbm;sbm_amr_qualification"
      RUN_SERIAL TRUE
      TIMEOUT 300)

    add_test(NAME SBM_AMR_Qualification_1M_Split_2Rank
      COMMAND ${CMAKE_COMMAND}
        -DERF_EXECUTABLE=${SBM_AMR_QUALIFICATION_EXECUTABLE}
        -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_amr_qualification_1m_split
        -DTEST_ROOT=${_sbm_amr_qualification_root}/1m_split_2rank
        -DMPIEXEC=${MPIEXEC_EXECUTABLE}
        -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
        -DMPI_RANKS=2
        -DEXPECTED_COMPONENTS=4
        -DEXPECTED_FINE_SUBSTEPS=3
        -DSPATIAL_SPLIT_FIXTURE=ON
        -DEXPECT_SPLIT=ON
        -DPYTHON_EXECUTABLE=${ERF_TEST_PYTHON}
        -DREFERENCE_LOG=${_sbm_amr_qualification_root}/1m_split_1rank/qualification.log
        -P ${_sbm_amr_qualification_script})
    set_tests_properties(SBM_AMR_Qualification_1M_Split_2Rank PROPERTIES
      LABELS "integration;parallel;sbm;sbm_amr_qualification"
      RUN_SERIAL TRUE
      PROCESSORS 2
      DEPENDS SBM_AMR_Qualification_1M_Split_1Rank
      TIMEOUT 300)

    add_test(NAME SBM_AMR_Qualification_RankLocalGhostStoreFailure_2Rank
      COMMAND ${CMAKE_COMMAND}
        -DERF_EXECUTABLE=${SBM_AMR_QUALIFICATION_EXECUTABLE}
        -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_amr_qualification_1m_split
        -DTEST_ROOT=${_sbm_amr_qualification_root}/rank_local_ghost_store_failure_2rank
        -DMPIEXEC=${MPIEXEC_EXECUTABLE}
        -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
        -DMPI_RANKS=2
        -DEXPECTED_COMPONENTS=4
        -DSPATIAL_SPLIT_FIXTURE=ON
        -DFAILURE_TEST=ON
        -P ${_sbm_amr_qualification_script})
    set_tests_properties(SBM_AMR_Qualification_RankLocalGhostStoreFailure_2Rank PROPERTIES
      LABELS "integration;parallel;sbm;sbm_amr_qualification;negative"
      RUN_SERIAL TRUE
      PROCESSORS 2
      TIMEOUT 60)

    add_test(NAME SBM_AMR_Qualification_2M_Attached_Split_1Rank
      COMMAND ${CMAKE_COMMAND}
        -DERF_EXECUTABLE=${SBM_AMR_QUALIFICATION_EXECUTABLE}
        -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_amr_qualification_2m_attached_split
        -DTEST_ROOT=${_sbm_amr_qualification_root}/2m_attached_split_1rank
        -DMPIEXEC=${MPIEXEC_EXECUTABLE}
        -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
        -DMPI_RANKS=1
        -DEXPECTED_COMPONENTS=12
        -DEXPECTED_FINE_SUBSTEPS=3
        -DSPATIAL_SPLIT_FIXTURE=ON
        -P ${_sbm_amr_qualification_script})
    set_tests_properties(SBM_AMR_Qualification_2M_Attached_Split_1Rank PROPERTIES
      LABELS "integration;sbm;sbm_amr_qualification"
      RUN_SERIAL TRUE
      TIMEOUT 300)

    add_test(NAME SBM_AMR_Qualification_2M_Attached_Split_2Rank
      COMMAND ${CMAKE_COMMAND}
        -DERF_EXECUTABLE=${SBM_AMR_QUALIFICATION_EXECUTABLE}
        -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_amr_qualification_2m_attached_split
        -DTEST_ROOT=${_sbm_amr_qualification_root}/2m_attached_split_2rank
        -DMPIEXEC=${MPIEXEC_EXECUTABLE}
        -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
        -DMPI_RANKS=2
        -DEXPECTED_COMPONENTS=12
        -DEXPECTED_FINE_SUBSTEPS=3
        -DSPATIAL_SPLIT_FIXTURE=ON
        -DEXPECT_SPLIT=ON
        -DPYTHON_EXECUTABLE=${ERF_TEST_PYTHON}
        -DREFERENCE_LOG=${_sbm_amr_qualification_root}/2m_attached_split_1rank/qualification.log
        -P ${_sbm_amr_qualification_script})
    set_tests_properties(SBM_AMR_Qualification_2M_Attached_Split_2Rank PROPERTIES
      LABELS "integration;parallel;sbm;sbm_amr_qualification"
      RUN_SERIAL TRUE
      PROCESSORS 2
      DEPENDS SBM_AMR_Qualification_2M_Attached_Split_1Rank
      TIMEOUT 300)

    # Measure the actual coarse-parent halo path at the requested spectral
    # widths and process counts with the nonzero carrier from the profile
    # deck.  This keeps complete fine-step timings representative of transport;
    # the split runs also assert remote parent and fine-fine donor ownership.
    foreach(_sbm_profile_ncomp IN ITEMS 32 128 256)
      foreach(_sbm_profile_ranks IN ITEMS 1 2 4)
        set(_sbm_profile_test
            "SBM_AMR_Profile_N${_sbm_profile_ncomp}_R${_sbm_profile_ranks}")
        set(_sbm_profile_root
            "${_sbm_amr_qualification_root}/profile_n${_sbm_profile_ncomp}_r${_sbm_profile_ranks}")
        set(_sbm_profile_args
            -DERF_EXECUTABLE=${SBM_AMR_QUALIFICATION_EXECUTABLE}
            -DINPUT_FILE=${PROJECT_SOURCE_DIR}/Tests/test_files/SBM/inputs_amr_qualification_profile
            -DTEST_ROOT=${_sbm_profile_root}
            -DMPIEXEC=${MPIEXEC_EXECUTABLE}
            -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
            -DMPI_RANKS=${_sbm_profile_ranks}
            -DEXPECTED_COMPONENTS=${_sbm_profile_ncomp}
            -DEXPECTED_STEPS=1
            -DEXPECTED_FINE_SUBSTEPS=3
            -DSBM_NBINS=${_sbm_profile_ncomp}
            -DMAX_STEPS=1
            -DSPATIAL_SPLIT_FIXTURE=ON)
        if(_sbm_profile_ranks GREATER 1)
          list(APPEND _sbm_profile_args
               -DEXPECT_SPLIT=ON
               -DPYTHON_EXECUTABLE=${ERF_TEST_PYTHON}
               -DREFERENCE_LOG=${_sbm_amr_qualification_root}/profile_n${_sbm_profile_ncomp}_r1/qualification.log)
        endif()
        list(APPEND _sbm_profile_args -P ${_sbm_amr_qualification_script})
        add_test(NAME ${_sbm_profile_test}
          COMMAND ${CMAKE_COMMAND} ${_sbm_profile_args})
        set_tests_properties(${_sbm_profile_test} PROPERTIES
          LABELS "integration;parallel;sbm;sbm_amr_qualification;sbm_amr_profile"
          RUN_SERIAL TRUE
          PROCESSORS ${_sbm_profile_ranks}
          TIMEOUT 600)
        if(_sbm_profile_ranks GREATER 1)
          set_tests_properties(${_sbm_profile_test} PROPERTIES
            DEPENDS SBM_AMR_Profile_N${_sbm_profile_ncomp}_R1)
        endif()
      endforeach()
    endforeach()
  endif()

  if(ERF_ENABLE_MPI)
    # Open MPI 4 keeps one session directory per user under TMPDIR; concurrent mpiexec
    # launches race on it. Open MPI 5 makes the directory unique per launch itself.
    get_property(_erf_tests DIRECTORY PROPERTY TESTS)
    foreach(_erf_test IN LISTS _erf_tests)
      set_property(TEST ${_erf_test} APPEND PROPERTY ENVIRONMENT
        "TMPDIR=${CMAKE_CURRENT_BINARY_DIR}/test_files/${_erf_test}")
    endforeach()
  endif()
endif()

if(ERF_ENABLE_UNIT_TESTS)
  add_subdirectory(Unit)
endif()
