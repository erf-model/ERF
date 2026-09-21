# Run a deck that is expected to fail, and require the failure message to name the
# problem.  Input validation is only useful if it tells the user which key is wrong,
# so the message itself is what is tested, not just the non-zero exit status.
foreach(arg TEST_EXE INPUT WORKING_DIRECTORY LOG EXPECTED_MESSAGE)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunExpectAbort.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()

separate_arguments(runtime_options  UNIX_COMMAND "${RUNTIME_OPTIONS}")
separate_arguments(mpiexec_preflags UNIX_COMMAND "${MPIEXEC_PREFLAGS}")

if(NOT "${MPIEXEC}" STREQUAL "")
    if("${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
        message(FATAL_ERROR "RunExpectAbort.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
    endif()
    set(launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${mpiexec_preflags})
else()
    set(launch)
endif()

execute_process(
    COMMAND ${launch} ${TEST_EXE} ${INPUT} ${runtime_options}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${LOG}"
    ERROR_FILE "${LOG}"
    RESULT_VARIABLE run_result)

if(run_result EQUAL 0)
    message(FATAL_ERROR "RunExpectAbort.cmake: the run succeeded, but it was expected to fail with "
                        "'${EXPECTED_MESSAGE}' (see ${LOG})")
endif()

file(READ "${LOG}" log_text)
string(FIND "${log_text}" "${EXPECTED_MESSAGE}" found_at)
if(found_at EQUAL -1)
    message(FATAL_ERROR "RunExpectAbort.cmake: the run failed, but without saying "
                        "'${EXPECTED_MESSAGE}' (see ${LOG})")
endif()
message(STATUS "RunExpectAbort: failed as expected with '${EXPECTED_MESSAGE}'")
