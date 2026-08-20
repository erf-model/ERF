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

  if(NOT ERF_GTEST_EXECUTION_TIMEOUT)
    set(ERF_GTEST_EXECUTION_TIMEOUT 900)
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
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestWrfBoundaryDensity.cpp
    # Motivation: face-local wall metadata must not inherit a later wall's
    # density representation or treat explicit zero as absent.
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions/ERF_GTestWallBCMetadata.cpp
    # Motivation: the chamber profile and perturbation must remain analytic,
    # bounded, and independent of MPI decomposition and tile shape.
    ${CMAKE_CURRENT_SOURCE_DIR}/Prob/ERF_GTestCloudChamber.cpp
    # Motivation: production scalar diffusion must use the exact low/high
    # physical-face stencil in every rotated Cartesian direction.
    ${CMAKE_CURRENT_SOURCE_DIR}/Diffusion/ERF_GTestStateDiffusionDirichlet.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Diagnostics/ERF_GTestSurfaceFluxDiagnostics.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Diagnostics/ERF_GTestSeaLevelPressure.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestPlotfile2D.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestPlotfileSelection.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestProvenance.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ERF_GTestEpochTime.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/EOS/ERF_GTestEOSScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/EOS/ERF_GTestEOSKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/Interpolation/ERF_GTestInterpolationScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/Interpolation/ERF_GTestInterpolationKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/Interpolation/ERF_GTestInterpolationPolicyKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/Interpolation/ERF_GTestInterpolationWENOScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/Interpolation/ERF_GTestInterpolationWENOKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/Microphysics/ERF_GTestMicrophysicsUtilsScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/Microphysics/ERF_GTestMicrophysicsUtilsKernel.cpp
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
)
# These exercise the super-droplet phase-change kernels, whose headers are
# only available when particles are enabled.
if(ERF_ENABLE_PARTICLES)
  target_sources(${erf_exe_name} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/TerminalVelocity/ERF_GTestTerminalVelocityScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/TerminalVelocity/ERF_GTestTerminalVelocityKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/IceDeposition/ERF_GTestIceDepositionScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/IceDeposition/ERF_GTestIceDepositionKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/Condensation/ERF_GTestCondensationScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/Condensation/ERF_GTestCondensationKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/Melting/ERF_GTestMeltingScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/Melting/ERF_GTestMeltingKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/Freezing/ERF_GTestFreezingScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/Freezing/ERF_GTestFreezingKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/Coalescence/ERF_GTestCoalescenceScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/Coalescence/ERF_GTestCoalescenceKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/Riming/ERF_GTestRimingScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/Riming/ERF_GTestRimingKernel.cpp
  )
endif()

if(ERF_ENABLE_NOAHMP)
  target_sources(${erf_exe_name} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/LandSurfaceModel/Noah-MP/ERF_GTestNoahMPResultPolicy.cpp
  )
endif()
build_erf_exe(${erf_exe_name} NO_ERF_MAIN)
target_link_libraries(${erf_exe_name} PRIVATE GTest::gtest)

erf_register_gtest_suite(${erf_exe_name}
  LABELS unit
  DISCOVERY_TIMEOUT 60
  EXECUTION_TIMEOUT 900
)

add_test(NAME erf_unit_tests_shuffle_repeat
  COMMAND ${erf_exe_name}
    --gtest_shuffle
    --gtest_repeat=3
    --gtest_output=xml:${ERF_GTEST_REPORT_ROOT}/${erf_exe_name}/shuffle-repeat.xml
)
set_tests_properties(erf_unit_tests_shuffle_repeat
  PROPERTIES
    LABELS "stress"
    TIMEOUT 2700
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
      ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestProvenanceParallel.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SatAdj/ERF_GTestSatAdjParallel.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SAM/ERF_GTestSAMParallel.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestProvenanceParallel.cpp
  )
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
      -DTEST_TIMEOUT=1800
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
