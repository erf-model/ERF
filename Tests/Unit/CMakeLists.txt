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
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/Kessler/ERF_GTestKesslerScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/Kessler/ERF_GTestKesslerKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/Kessler/ERF_GTestKesslerPhysicalProperties.cpp
)
build_erf_exe(${erf_exe_name} NO_ERF_MAIN)
target_link_libraries(${erf_exe_name} PRIVATE GTest::gtest)

include(GoogleTest)
gtest_discover_tests(${erf_exe_name}
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS "unit"
)

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

  set(erf_parallel_test_mpiexec_preflags)

  if(DEFINED ERF_PARALLEL_TEST_MPIEXEC_PREFLAGS AND
     NOT "${ERF_PARALLEL_TEST_MPIEXEC_PREFLAGS}" STREQUAL "")
    list(APPEND erf_parallel_test_mpiexec_preflags
      ${ERF_PARALLEL_TEST_MPIEXEC_PREFLAGS})
  endif()

  if(ERF_PARALLEL_TEST_ALLOW_OVERSUBSCRIBE)
    execute_process(
      COMMAND ${MPIEXEC_EXECUTABLE} --version
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
  )
  build_erf_exe(${erf_parallel_exe_name} NO_ERF_MAIN)
  target_link_libraries(${erf_parallel_exe_name}
    PRIVATE
      GTest::gtest
      MPI::MPI_CXX
  )

  foreach(erf_parallel_nprocs IN LISTS ERF_PARALLEL_TEST_NRANKS)
    if("${erf_parallel_nprocs}" STREQUAL "")
      message(FATAL_ERROR
        "ERF_PARALLEL_TEST_NRANKS contains an empty rank count")
    endif()

    if(NOT "${erf_parallel_nprocs}" MATCHES "^[1-9][0-9]*$")
      message(FATAL_ERROR
        "ERF_PARALLEL_TEST_NRANKS entries must be positive integers; got '${erf_parallel_nprocs}'")
    endif()

    set(erf_parallel_test_name
        "${erf_parallel_exe_name}_np${erf_parallel_nprocs}")

    add_test(NAME ${erf_parallel_test_name}
      COMMAND
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} ${erf_parallel_nprocs}
        ${MPIEXEC_PREFLAGS}
        ${erf_parallel_test_mpiexec_preflags}
        $<TARGET_FILE:${erf_parallel_exe_name}>
        ${MPIEXEC_POSTFLAGS}
      COMMAND_EXPAND_LISTS
    )

    set_tests_properties(${erf_parallel_test_name}
      PROPERTIES
        LABELS "parallel"
        PROCESSORS ${erf_parallel_nprocs}
        ENVIRONMENT "ERF_PARALLEL_TEST_NPROCS=${erf_parallel_nprocs}"
    )
  endforeach()
endif()
