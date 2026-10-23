include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")

foreach(required MPIEXEC MPIEXEC_NUMPROC_FLAG TEST_EXE INPUT WORKING_DIRECTORY
                 EXPECTED_MODEL EXPECTED_TEXT LOG_FILE)
    if(NOT DEFINED ${required})
        message(FATAL_ERROR "RunCloudChamberParserReject.cmake missing ${required}")
    endif()
endforeach()

erf_mpi_launcher_command(parser_launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunCloudChamberParserReject.cmake")

file(READ "${INPUT}" template)
file(REMOVE_RECURSE "${WORKING_DIRECTORY}/positive" "${WORKING_DIRECTORY}/rejected")
file(MAKE_DIRECTORY "${WORKING_DIRECTORY}/positive" "${WORKING_DIRECTORY}/rejected")

function(run_parser model case_dir log_path)
    string(REPLACE "@MODEL@" "${model}" case_input "${template}")
    set(input_path "${case_dir}/CloudChamber_ActiveCondensation.i")
    file(WRITE "${input_path}" "${case_input}")
    execute_process(
        COMMAND ${parser_launch} "${TEST_EXE}" "${input_path}" amrex.call_addr2line=0
        WORKING_DIRECTORY "${case_dir}"
        OUTPUT_FILE "${log_path}"
        ERROR_FILE "${log_path}"
        TIMEOUT 30
        RESULT_VARIABLE result)
    set(${model}_result "${result}" PARENT_SCOPE)
endfunction()

run_parser(Kessler_NoRain "${WORKING_DIRECTORY}/positive"
           "${WORKING_DIRECTORY}/positive.log")
if(NOT "${Kessler_NoRain_result}" STREQUAL "0")
    message(FATAL_ERROR "Supported parser positive control failed (${Kessler_NoRain_result})")
endif()

run_parser("${EXPECTED_MODEL}" "${WORKING_DIRECTORY}/rejected" "${LOG_FILE}")
if(NOT "${${EXPECTED_MODEL}_result}" MATCHES "^-?[0-9]+$")
    message(FATAL_ERROR
        "Parser rejection timed out or launcher failed (${${EXPECTED_MODEL}_result}); see ${LOG_FILE}")
endif()
if("${${EXPECTED_MODEL}_result}" STREQUAL "0")
    message(FATAL_ERROR "Cloud Chamber parser accepted unsupported ${EXPECTED_MODEL}")
endif()

file(READ "${LOG_FILE}" failure_output)
string(FIND "${failure_output}" "${EXPECTED_TEXT}" expected_position)
if(expected_position EQUAL -1)
    message(FATAL_ERROR
        "Cloud Chamber parser failed without its allowlist diagnostic '${EXPECTED_TEXT}'; see ${LOG_FILE}")
endif()
