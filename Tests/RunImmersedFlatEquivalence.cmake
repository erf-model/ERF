# Run the flat-ground deck and the immersed-plane deck of one test directory (COMMON_OPTIONS on
# both runs, RUNTIME_OPTIONS on the plane run only), then compare their
# plotfiles with Tests/ImmersedFlatEquivalenceCheck.cpp (the plane run's profiles from cell
# K_OFFSET up against the flat run's). With EXPECT_CHECKER_FAIL the check must fail (the leg that
# runs the legacy wall law), so the test proves the checker sees the difference.
# -DX= defines X as empty, so test for a value, not for DEFINED
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

foreach(arg NRANKS TEST_EXE INPUT_FLAT INPUT_PLANE WORKING_DIRECTORY SIMULATION_LOG CHECKER_LOG CHECKER
            PLOTFILE_FLAT PLOTFILE_PLANE K_OFFSET TOL K_COUNT EXPECT_CHECKER_FAIL)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunImmersedFlatEquivalence.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunImmersedFlatEquivalence.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunImmersedFlatEquivalence.cmake: ERF executable")
erf_resolve_executable(CHECKER "${CHECKER}" CONFIG "${CONFIG}"
    CONTEXT "RunImmersedFlatEquivalence.cmake: checker")

erf_mpi_launcher_command(mpi_launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunImmersedFlatEquivalence.cmake")
erf_mpi_launcher_command(mpi_launch_one
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunImmersedFlatEquivalence.cmake")

separate_arguments(runtime_options UNIX_COMMAND "${RUNTIME_OPTIONS}")
separate_arguments(common_options UNIX_COMMAND "${COMMON_OPTIONS}")
execute_process(
    COMMAND ${mpi_launch} ${TEST_EXE} ${INPUT_FLAT} ${common_options}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${SIMULATION_LOG}.flat"
    ERROR_FILE "${SIMULATION_LOG}.flat"
    RESULT_VARIABLE flat_result)
if(NOT flat_result EQUAL 0)
    message(FATAL_ERROR "flat-ground simulation failed: ${flat_result} (see ${SIMULATION_LOG}.flat)")
endif()
execute_process(
    COMMAND ${mpi_launch} ${TEST_EXE} ${INPUT_PLANE} ${common_options} ${runtime_options}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${SIMULATION_LOG}"
    ERROR_FILE "${SIMULATION_LOG}"
    RESULT_VARIABLE plane_result)
if(NOT plane_result EQUAL 0)
    message(FATAL_ERROR "immersed-plane simulation failed: ${plane_result} (see ${SIMULATION_LOG})")
endif()

execute_process(
    COMMAND ${mpi_launch_one} ${CHECKER} ${PLOTFILE_FLAT} ${PLOTFILE_PLANE} ${K_OFFSET} ${TOL} ${K_COUNT}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${CHECKER_LOG}"
    ERROR_FILE "${CHECKER_LOG}"
    RESULT_VARIABLE checker_result)
file(READ "${CHECKER_LOG}" checker_output)
if(EXPECT_CHECKER_FAIL)
    if(checker_result EQUAL 0)
        message(FATAL_ERROR "the legacy wall law passed the flat-ground comparison; the check cannot see the difference\n${checker_output}")
    elseif(NOT checker_result EQUAL 1 OR NOT checker_output MATCHES "differs from flat ground")
        message(FATAL_ERROR "the checker did not reach the comparison (exit ${checker_result})\n${checker_output}")
    endif()
    message(STATUS "the legacy wall law fails the comparison, as expected:\n${checker_output}")
elseif(NOT checker_result EQUAL 0)
    message(FATAL_ERROR "immersed flat-equivalence check failed: ${checker_result}\n${checker_output}")
else()
    message(STATUS "${checker_output}")
endif()
