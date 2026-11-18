# Run a two-level immersed-forcing case, then check its plotfile with
# Tests/ImmersedAverageDownCheck.cpp. With EXPECT_CHECKER_FAIL the check must fail (the leg that
# runs the plain average), so the test proves the checker sees the difference.
# -DX= defines X as empty, so test for a value, not for DEFINED
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

foreach(arg NRANKS TEST_EXE INPUT WORKING_DIRECTORY SIMULATION_LOG CHECKER_LOG CHECKER
            PLOTFILE TOL EXPECT_CHECKER_FAIL)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunImmersedAverageDown.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunImmersedAverageDown.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunImmersedAverageDown.cmake: ERF executable")
erf_resolve_executable(CHECKER "${CHECKER}" CONFIG "${CONFIG}"
    CONTEXT "RunImmersedAverageDown.cmake: checker")

erf_mpi_launcher_command(mpi_launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunImmersedAverageDown.cmake")
erf_mpi_launcher_command(mpi_launch_one
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunImmersedAverageDown.cmake")

separate_arguments(runtime_options UNIX_COMMAND "${RUNTIME_OPTIONS}")
execute_process(
    COMMAND ${mpi_launch} ${TEST_EXE} ${INPUT} ${runtime_options}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${SIMULATION_LOG}"
    ERROR_FILE "${SIMULATION_LOG}"
    RESULT_VARIABLE simulation_result)
if(NOT simulation_result EQUAL 0)
    message(FATAL_ERROR "simulation failed: ${simulation_result} (see ${SIMULATION_LOG})")
endif()

execute_process(
    COMMAND ${mpi_launch_one} ${CHECKER} ${PLOTFILE} ${TOL}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${CHECKER_LOG}"
    ERROR_FILE "${CHECKER_LOG}"
    RESULT_VARIABLE checker_result)
file(READ "${CHECKER_LOG}" checker_output)
if(EXPECT_CHECKER_FAIL)
    if(checker_result EQUAL 0)
        message(FATAL_ERROR "the plain average passed the fluid-mass check; the check cannot see the difference\n${checker_output}")
    elseif(NOT checker_result EQUAL 1 OR NOT checker_output MATCHES "not the fluid-mass-weighted fine average")
        message(FATAL_ERROR "the checker did not reach the comparison (exit ${checker_result})\n${checker_output}")
    endif()
    message(STATUS "the plain average fails the check, as expected:\n${checker_output}")
elseif(NOT checker_result EQUAL 0)
    message(FATAL_ERROR "immersed average-down check failed: ${checker_result}\n${checker_output}")
else()
    message(STATUS "${checker_output}")
endif()
