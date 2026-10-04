if(NOT DEFINED ERF_EXECUTABLE OR NOT EXISTS "${ERF_EXECUTABLE}")
  message(FATAL_ERROR "ERF_EXECUTABLE must name a built erf_exec")
endif()
if(NOT DEFINED COMPARE_EXECUTABLE OR NOT EXISTS "${COMPARE_EXECUTABLE}")
  message(FATAL_ERROR "COMPARE_EXECUTABLE must name the SBM checkpoint comparator")
endif()
if(NOT DEFINED INPUT_FILE OR NOT EXISTS "${INPUT_FILE}")
  message(FATAL_ERROR "INPUT_FILE must name the static-terrain SBM fixture")
endif()
if(NOT DEFINED TEST_ROOT)
  message(FATAL_ERROR "TEST_ROOT must name an isolated test output directory")
endif()

string(RANDOM LENGTH 12 ALPHABET 0123456789abcdef _run_id)
set(_run_dir "${TEST_ROOT}/sbm_restart_equivalence_${_run_id}")
file(MAKE_DIRECTORY "${_run_dir}")

function(run_erf description)
  execute_process(
    COMMAND "${ERF_EXECUTABLE}" "${INPUT_FILE}" ${ARGN}
    WORKING_DIRECTORY "${_run_dir}"
    RESULT_VARIABLE _result
    OUTPUT_VARIABLE _stdout
    ERROR_VARIABLE _stderr)
  if(NOT "${_result}" STREQUAL "0")
    message(FATAL_ERROR
      "${description} failed (exit ${_result}).\nstdout:\n${_stdout}\nstderr:\n${_stderr}")
  endif()
endfunction()

run_erf("Continuous two-step SBM run"
  "max_step=2" "stop_time=2.e-5"
  "erf.check_file=sbm_continuous_chk" "erf.check_int=1"
  "erf.plot_int_1=-1")

set(_continuous_step_one "${_run_dir}/sbm_continuous_chk00001")
set(_continuous_step_two "${_run_dir}/sbm_continuous_chk00002")
foreach(_checkpoint IN ITEMS "${_continuous_step_one}" "${_continuous_step_two}")
  if(NOT EXISTS "${_checkpoint}/Level_0/SBMSpectrum_H" OR
     NOT EXISTS "${_checkpoint}/Level_0/Cell_H")
    message(FATAL_ERROR "Continuous run did not write the expected SBM checkpoint: ${_checkpoint}")
  endif()
endforeach()

run_erf("One-step restart from continuous checkpoint 1"
  "erf.restart=${_continuous_step_one}"
  "max_step=2" "stop_time=2.e-5"
  "erf.check_file=sbm_restarted_chk" "erf.check_int=1"
  "erf.plot_int_1=-1")

set(_restarted_step_two "${_run_dir}/sbm_restarted_chk00002")
if(NOT EXISTS "${_restarted_step_two}/Level_0/SBMSpectrum_H" OR
   NOT EXISTS "${_restarted_step_two}/Level_0/Cell_H")
  message(FATAL_ERROR "Restarted run did not write step-2 checkpoint: ${_restarted_step_two}")
endif()

execute_process(
  COMMAND "${COMPARE_EXECUTABLE}" "${_continuous_step_two}" "${_restarted_step_two}"
  WORKING_DIRECTORY "${_run_dir}"
  RESULT_VARIABLE _compare_result
  OUTPUT_VARIABLE _compare_stdout
  ERROR_VARIABLE _compare_stderr)
if(NOT "${_compare_result}" STREQUAL "0")
  message(FATAL_ERROR
    "Full checkpoint state comparison failed (exit ${_compare_result}).\n"
    "stdout:\n${_compare_stdout}\nstderr:\n${_compare_stderr}")
endif()
message(STATUS "Continuous and restarted step-2 SBM spectrum and qc/qr agree.\n${_compare_stdout}")
