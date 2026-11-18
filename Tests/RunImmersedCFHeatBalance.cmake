# Run a one- or two-level immersed-forcing case, then check its plotfile with
# Tests/ImmersedCFHeatBalanceCheck.cpp. With EXPECT_CHECKER_FAIL the check must fail (the leg that
# runs the explicit solid freeze), so the test proves the checker sees the difference. With
# EXPECT_WARNING the run must have printed the coarse/fine surface-mismatch warning (and without it,
# must not have).
# -DX= defines X as empty, so test for a value, not for DEFINED
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

foreach(arg NRANKS TEST_EXE INPUT WORKING_DIRECTORY SIMULATION_LOG CHECKER_LOG CHECKER
            PLOTFILE THETA0 TOL EXPECT_CHECKER_FAIL EXPECT_WARNING)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunImmersedCFHeatBalance.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunImmersedCFHeatBalance.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunImmersedCFHeatBalance.cmake: ERF executable")
erf_resolve_executable(CHECKER "${CHECKER}" CONFIG "${CONFIG}"
    CONTEXT "RunImmersedCFHeatBalance.cmake: checker")

erf_mpi_launcher_command(mpi_launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunImmersedCFHeatBalance.cmake")
erf_mpi_launcher_command(mpi_launch_one
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunImmersedCFHeatBalance.cmake")

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
# With EXPECT_WARNING the start-up warning of erf.if_cf_mismatch_warning must have been printed
file(READ "${SIMULATION_LOG}" simulation_output)
if(EXPECT_WARNING)
    if(NOT simulation_output MATCHES "coarse cells along the lateral coarse/fine faces")
        message(FATAL_ERROR "the coarse/fine surface-mismatch warning was not printed (see ${SIMULATION_LOG})")
    endif()
elseif(simulation_output MATCHES "coarse cells along the lateral coarse/fine faces")
    message(FATAL_ERROR "the coarse/fine surface-mismatch warning was printed at the default threshold (see ${SIMULATION_LOG})")
endif()

execute_process(
    COMMAND ${mpi_launch_one} ${CHECKER} ${PLOTFILE} ${THETA0} ${TOL}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${CHECKER_LOG}"
    ERROR_FILE "${CHECKER_LOG}"
    RESULT_VARIABLE checker_result)
file(READ "${CHECKER_LOG}" checker_output)
if(EXPECT_CHECKER_FAIL)
    if(checker_result EQUAL 0)
        message(FATAL_ERROR "the explicit freeze passed the heat-balance check; the check cannot see the difference\n${checker_output}")
    elseif(NOT checker_result EQUAL 1 OR NOT checker_output MATCHES "drifted by more than tol")
        message(FATAL_ERROR "the checker did not reach the comparison (exit ${checker_result})\n${checker_output}")
    endif()
    message(STATUS "the explicit freeze fails the check, as expected:\n${checker_output}")
elseif(NOT checker_result EQUAL 0)
    message(FATAL_ERROR "immersed heat-balance check failed: ${checker_result}\n${checker_output}")
else()
    message(STATUS "${checker_output}")
endif()
