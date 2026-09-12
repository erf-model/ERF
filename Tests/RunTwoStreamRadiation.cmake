if(NOT DEFINED TEST_EXE OR NOT DEFINED INPUT OR
   NOT DEFINED WORKING_DIRECTORY OR NOT DEFINED SIMULATION_LOG OR
   NOT DEFINED CHECKER_LOG OR NOT DEFINED CHECKER OR NOT DEFINED PLOTFILE)
    message(FATAL_ERROR "RunTwoStreamRadiation.cmake missing required argument")
endif()

# Build the launcher prefix. A build without MPI passes an empty MPIEXEC, in
# which case the executables are run directly instead of through a launcher.
function(two_stream_launcher nranks out_var)
    if(DEFINED MPIEXEC AND NOT "${MPIEXEC}" STREQUAL "")
        set(launcher ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${nranks} ${MPIEXEC_PREFLAGS})
    else()
        set(launcher "")
    endif()
    set(${out_var} "${launcher}" PARENT_SCOPE)
endfunction()

# Echo a log to the test output so a CI failure is diagnosable without having
# to fetch the attached files.
function(two_stream_report_log label path)
    if(EXISTS "${path}")
        file(READ "${path}" contents)
        message(STATUS "---- ${label} (${path}) ----\n${contents}\n---- end ${label} ----")
    else()
        message(STATUS "---- ${label}: ${path} was never written ----")
    endif()
endfunction()

if(NOT DEFINED NRANKS OR "${NRANKS}" STREQUAL "")
    set(NRANKS 1)
endif()

two_stream_launcher(${NRANKS} simulation_launcher)
execute_process(
    COMMAND ${simulation_launcher} ${TEST_EXE} ${INPUT}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${SIMULATION_LOG}"
    ERROR_FILE "${SIMULATION_LOG}"
    RESULT_VARIABLE simulation_result)
if(NOT simulation_result EQUAL 0)
    two_stream_report_log("simulation log" "${SIMULATION_LOG}")
    message(FATAL_ERROR "TwoStream radiation simulation failed: ${simulation_result}")
endif()

two_stream_launcher(1 checker_launcher)
execute_process(
    COMMAND ${checker_launcher} ${CHECKER} ${PLOTFILE}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${CHECKER_LOG}"
    ERROR_FILE "${CHECKER_LOG}"
    RESULT_VARIABLE checker_result)
if(NOT checker_result EQUAL 0)
    two_stream_report_log("checker log" "${CHECKER_LOG}")
    message(FATAL_ERROR "TwoStream radiation column check failed: ${checker_result}")
endif()
