if(NOT DEFINED ERF_EXECUTABLE OR NOT EXISTS "${ERF_EXECUTABLE}")
  message(FATAL_ERROR "ERF_EXECUTABLE must name a built erf_exec")
endif()
if(NOT DEFINED INPUT_FILE OR NOT EXISTS "${INPUT_FILE}")
  message(FATAL_ERROR "INPUT_FILE must name an ERF input deck")
endif()
if(NOT DEFINED TEST_ROOT OR NOT DEFINED EXPECTED_DIAGNOSTIC)
  message(FATAL_ERROR "TEST_ROOT and EXPECTED_DIAGNOSTIC are required")
endif()

string(RANDOM LENGTH 12 ALPHABET 0123456789abcdef _run_id)
set(_run_dir "${TEST_ROOT}/sbm_expected_failure_${_run_id}")
file(MAKE_DIRECTORY "${_run_dir}")
set(_command "${ERF_EXECUTABLE}" "${INPUT_FILE}")
if(DEFINED OVERRIDE AND NOT "${OVERRIDE}" STREQUAL "")
  list(APPEND _command "${OVERRIDE}")
endif()
execute_process(
  COMMAND ${_command}
  WORKING_DIRECTORY "${_run_dir}"
  RESULT_VARIABLE _result
  OUTPUT_VARIABLE _stdout
  ERROR_VARIABLE _stderr)
if("${_result}" STREQUAL "0")
  message(FATAL_ERROR "ERF unexpectedly succeeded; expected '${EXPECTED_DIAGNOSTIC}'")
endif()
set(_output "${_stdout}\n${_stderr}")
if(NOT _output MATCHES "${EXPECTED_DIAGNOSTIC}")
  message(FATAL_ERROR
    "ERF failed without the expected diagnostic '${EXPECTED_DIAGNOSTIC}' (exit ${_result}).\nstdout:\n${_stdout}\nstderr:\n${_stderr}")
endif()
message(STATUS "Observed expected ERF failure: ${EXPECTED_DIAGNOSTIC}")
