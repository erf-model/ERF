include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")

foreach(required MPIEXEC MPIEXEC_NUMPROC_FLAG NRANKS TEST_EXE INPUT
                 WORKING_DIRECTORY SIMULATION_LOG CHECKER_LOG CHECKER)
    if(NOT DEFINED ${required})
        message(FATAL_ERROR "RunCloudChamberActiveCondensation.cmake missing ${required}")
    endif()
endforeach()

erf_mpi_launcher_command(mpi_launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunCloudChamberActiveCondensation.cmake")
erf_mpi_launcher_command(check_launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunCloudChamberActiveCondensation.cmake checker")

file(READ "${INPUT}" template)
file(REMOVE_RECURSE "${WORKING_DIRECTORY}/active" "${WORKING_DIRECTORY}/control")
file(MAKE_DIRECTORY "${WORKING_DIRECTORY}/active" "${WORKING_DIRECTORY}/control")

function(run_case model case_dir log_path)
    set(input_path "${case_dir}/CloudChamber_ActiveCondensation.i")
    string(REPLACE "@MODEL@" "${model}" case_input "${template}")
    file(WRITE "${input_path}" "${case_input}")
    execute_process(
        COMMAND ${mpi_launch} "${TEST_EXE}" "${input_path}" amrex.call_addr2line=0
        WORKING_DIRECTORY "${case_dir}"
        OUTPUT_FILE "${log_path}"
        ERROR_FILE "${log_path}"
        TIMEOUT 600
        RESULT_VARIABLE result)
    if(NOT "${result}" STREQUAL "0")
        message(FATAL_ERROR "${model} chamber run failed or timed out (${result}); see ${log_path}")
    endif()
endfunction()

run_case(Kessler_NoRain "${WORKING_DIRECTORY}/active" "${SIMULATION_LOG}")
run_case(MoistNoCondensation "${WORKING_DIRECTORY}/control" "${WORKING_DIRECTORY}/control.simulation.log")

execute_process(
    COMMAND ${check_launch} "${CHECKER}" active_condensation Kessler_NoRain
            "${WORKING_DIRECTORY}/active/plt00000"
            "${WORKING_DIRECTORY}/active/plt01000"
            "${WORKING_DIRECTORY}/control/plt01000"
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${CHECKER_LOG}"
    ERROR_FILE "${CHECKER_LOG}"
    TIMEOUT 120
    RESULT_VARIABLE checker_result)
if(NOT "${checker_result}" STREQUAL "0")
    message(FATAL_ERROR "Active-condensation oracle failed or timed out (${checker_result}); see ${CHECKER_LOG}")
endif()
