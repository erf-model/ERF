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
  option(ERF_TEST_ENABLE_EXTRA_SDM_TESTS "Enable additional SDM tests that require machine-specific gold files" OFF)
  option(ERF_TEST_SDM_SMOKE_GATED "Run the machine-specific SDM tests to completion, without gold comparison, under the smoke label" OFF)

  message(STATUS "ERF testing configuration summary:")
  message(STATUS "   Number of ranks               = ${ERF_TEST_NRANKS}")
  message(STATUS "   fcompare executable           = ${FCOMPARE_EXE}")
  message(STATUS "   comparison relative tolerance = ${ERF_TEST_FCOMPARE_RTOL}")
  message(STATUS "   comparison absolute tolerance = ${ERF_TEST_FCOMPARE_ATOL}")
  message(STATUS "   Gold file directory           = ${ERF_TEST_GOLD_FILES_DIRECTORY}")
  message(STATUS "   Extra SDM tests enabled       = ${ERF_TEST_ENABLE_EXTRA_SDM_TESTS}")
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

  # Motivation: the chamber regressions must prove the initializer is wired
  # into a real 3-D run and that the SatAdj deck forms finite cloud water.
  add_executable(erf_cloud_chamber_check
    ${CMAKE_CURRENT_SOURCE_DIR}/CloudChamberCheck.cpp)
  target_link_libraries(erf_cloud_chamber_check PRIVATE AMReX::amrex)
  target_include_directories(erf_cloud_chamber_check PRIVATE
    ${PROJECT_SOURCE_DIR}/Source)
  target_compile_features(erf_cloud_chamber_check PRIVATE cxx_std_17)
  set(CLOUD_CHAMBER_CHECKER
      "${CMAKE_CURRENT_BINARY_DIR}/erf_cloud_chamber_check${CMAKE_EXECUTABLE_SUFFIX}")

  if(ERF_ENABLE_CUDA)
    foreach(_erf_test_tool IN ITEMS
        erf_shoc_plotfile_check
        erf_shoc_microphysics_differential
        erf_shoc_mutation_differential
        erf_shoc_gold_differential
        erf_anelastic_wall_diffusion_check
        erf_cloud_chamber_check)
      setup_target_for_cuda_compilation(${_erf_test_tool})
    endforeach()
  endif()

  include(${CMAKE_CURRENT_SOURCE_DIR}/CTestList.cmake)
endif()

if(ERF_ENABLE_UNIT_TESTS)
  add_subdirectory(Unit)
endif()
