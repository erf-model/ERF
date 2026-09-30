if(NOT DEFINED ERF_EXECUTABLE OR NOT EXISTS "${ERF_EXECUTABLE}")
  message(FATAL_ERROR "ERF_EXECUTABLE must name a built erf_exec")
endif()
if(NOT DEFINED TEST_HELPER OR NOT EXISTS "${TEST_HELPER}")
  message(FATAL_ERROR "TEST_HELPER must name the SBM checkpoint corruption helper")
endif()
if(NOT DEFINED INPUT_FILE OR NOT EXISTS "${INPUT_FILE}")
  message(FATAL_ERROR "INPUT_FILE must name the 2M SBM zero-transport input deck")
endif()
if(NOT DEFINED TEST_ROOT OR NOT DEFINED EXPECTED_DIAGNOSTIC)
  message(FATAL_ERROR "TEST_ROOT and EXPECTED_DIAGNOSTIC are required")
endif()

string(RANDOM LENGTH 12 ALPHABET 0123456789abcdef _run_id)
set(_run_dir "${TEST_ROOT}/sbm_restart_realizability_${_run_id}")
file(MAKE_DIRECTORY "${_run_dir}")
set(_checkpoint "${_run_dir}/sbm_source_chk00001")

execute_process(
  COMMAND "${ERF_EXECUTABLE}" "${INPUT_FILE}"
    "erf.check_file=sbm_source_chk" "erf.plot_int_1=-1"
  WORKING_DIRECTORY "${_run_dir}"
  RESULT_VARIABLE _initial_result
  OUTPUT_VARIABLE _initial_stdout
  ERROR_VARIABLE _initial_stderr)
if(NOT "${_initial_result}" STREQUAL "0")
  message(FATAL_ERROR
    "Valid 2M checkpoint setup failed (exit ${_initial_result}).\nstdout:\n${_initial_stdout}\nstderr:\n${_initial_stderr}")
endif()
if(NOT EXISTS "${_checkpoint}/SBM_Schema" OR
   NOT EXISTS "${_checkpoint}/Level_0/SBMSpectrum_H" OR
   NOT EXISTS "${_checkpoint}/Level_0/Cell_H")
  message(FATAL_ERROR "Valid 2M checkpoint setup did not produce the expected state files")
endif()

execute_process(
  COMMAND "${CMAKE_COMMAND}" -E env "SBM_CHECKPOINT_DIR=${_checkpoint}" "${TEST_HELPER}"
  WORKING_DIRECTORY "${_run_dir}"
  RESULT_VARIABLE _corrupt_result
  OUTPUT_VARIABLE _corrupt_stdout
  ERROR_VARIABLE _corrupt_stderr)
if(NOT "${_corrupt_result}" STREQUAL "0")
  message(FATAL_ERROR
    "Finite number-moment corruption helper failed (exit ${_corrupt_result}).\nstdout:\n${_corrupt_stdout}\nstderr:\n${_corrupt_stderr}")
endif()

execute_process(
  COMMAND "${ERF_EXECUTABLE}" "${INPUT_FILE}"
    "erf.restart=${_checkpoint}"
    "erf.check_file=sbm_rejected_restart_chk"
    "erf.plot_int_1=-1" "max_step=2" "stop_time=2.e-4"
  WORKING_DIRECTORY "${_run_dir}"
  RESULT_VARIABLE _restart_result
  OUTPUT_VARIABLE _restart_stdout
  ERROR_VARIABLE _restart_stderr)
if("${_restart_result}" STREQUAL "0")
  message(FATAL_ERROR "ERF accepted an inadmissible authoritative 2M restart")
endif()
set(_restart_output "${_restart_stdout}\n${_restart_stderr}")
if(NOT _restart_output MATCHES "${EXPECTED_DIAGNOSTIC}")
  message(FATAL_ERROR
    "ERF rejected the corrupted checkpoint without the expected authoritative-state diagnostic '${EXPECTED_DIAGNOSTIC}' (exit ${_restart_result}).\nstdout:\n${_restart_stdout}\nstderr:\n${_restart_stderr}")
endif()
message(STATUS "Observed expected authoritative-state restart rejection: ${EXPECTED_DIAGNOSTIC}")
