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
    if(model STREQUAL "dry")
        string(REPLACE "erf.moisture_model = dry\n" "" case_input "${case_input}")
        string(REPLACE "prob.initial_relative_humidity = 1.0\n" "" case_input "${case_input}")
    endif()
    set(input_path "${case_dir}/CloudChamber_ActiveCondensation.i")
    file(WRITE "${input_path}" "${case_input}")
    file(REMOVE "${log_path}")
    execute_process(
        COMMAND ${parser_launch} "${TEST_EXE}" "${input_path}" amrex.call_addr2line=0
        WORKING_DIRECTORY "${case_dir}"
        OUTPUT_FILE "${log_path}"
        ERROR_FILE "${log_path}"
        TIMEOUT 30
        RESULT_VARIABLE result)
    set(${model}_result "${result}" PARENT_SCOPE)
endfunction()

foreach(model IN ITEMS Kessler_NoRain kEsSlEr_NoRaIn MoIsTnOcOnDeNsAtIoN dry)
    run_parser("${model}" "${WORKING_DIRECTORY}/positive"
               "${WORKING_DIRECTORY}/positive.log")
    if(NOT "${${model}_result}" STREQUAL "0")
        message(FATAL_ERROR "Supported parser positive control ${model} failed (${${model}_result}); see ${WORKING_DIRECTORY}/positive.log")
    endif()
endforeach()

run_parser("${EXPECTED_MODEL}" "${WORKING_DIRECTORY}/rejected" "${LOG_FILE}")
set(parser_result "${${EXPECTED_MODEL}_result}")
string(TOLOWER "${parser_result}" parser_result_lower)
if("${parser_result}" STREQUAL "" OR "${parser_result_lower}" MATCHES "timeout")
    message(FATAL_ERROR
        "Parser rejection timed out or did not terminate (${parser_result}); see ${LOG_FILE}")
endif()
if("${parser_result}" STREQUAL "0")
    message(FATAL_ERROR "Cloud Chamber parser accepted rejected ${EXPECTED_MODEL}")
endif()
if(NOT "${parser_result}" MATCHES "^-?[0-9]+$" AND
   NOT "${parser_result_lower}" STREQUAL "subprocess aborted")
    message(FATAL_ERROR
        "Parser rejection had an unrecognized result (${parser_result}); see ${LOG_FILE}")
endif()

file(READ "${LOG_FILE}" failure_output)
string(FIND "${failure_output}" "${EXPECTED_TEXT}" expected_position)
if(expected_position EQUAL -1)
    message(FATAL_ERROR
    "Cloud Chamber parser failed without its expected diagnostic '${EXPECTED_TEXT}'; see ${LOG_FILE}")
endif()
