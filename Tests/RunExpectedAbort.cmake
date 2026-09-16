# Run a deck that must stop at start-up and check the abort message.
# The run must exit with a non-zero status, and its output must contain every string in
# EXPECTED (a ;-separated list), so a crash for an unrelated reason does not pass.
# -DX= defines X as empty, so test for a value, not for DEFINED
foreach(arg TEST_EXE INPUT WORKING_DIRECTORY LOG EXPECTED)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunExpectedAbort.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()

separate_arguments(runtime_options  UNIX_COMMAND "${RUNTIME_OPTIONS}")
separate_arguments(mpiexec_preflags UNIX_COMMAND "${MPIEXEC_PREFLAGS}")

set(launch_command "${TEST_EXE}")
if(NOT "${MPIEXEC}" STREQUAL "")
    set(launch_command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${mpiexec_preflags} "${TEST_EXE}")
endif()

execute_process(
    COMMAND ${launch_command} ${INPUT} ${runtime_options}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    RESULT_VARIABLE result
    OUTPUT_VARIABLE stdout
    ERROR_VARIABLE stderr)

set(output "${stdout}\n${stderr}")
file(WRITE "${LOG}" "${output}")

if("${result}" STREQUAL "0")
    message(FATAL_ERROR "RunExpectedAbort.cmake: the run returned success; it must abort at start-up")
endif()

foreach(expected IN LISTS EXPECTED)
    string(FIND "${output}" "${expected}" position)
    if(position EQUAL -1)
        message(FATAL_ERROR "RunExpectedAbort.cmake: missing expected abort text '${expected}'\n${output}")
    endif()
endforeach()
message(STATUS "RunExpectedAbort: aborted as expected (exit ${result})")
