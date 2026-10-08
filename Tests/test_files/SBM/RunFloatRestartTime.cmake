if(NOT DEFINED ERF_EXECUTABLE OR NOT EXISTS "${ERF_EXECUTABLE}")
  message(FATAL_ERROR "ERF_EXECUTABLE must name a built SINGLE erf_exec")
endif()
if(NOT DEFINED COMPARE_EXECUTABLE OR NOT EXISTS "${COMPARE_EXECUTABLE}")
  message(FATAL_ERROR "COMPARE_EXECUTABLE must name the SINGLE SBM checkpoint comparator")
endif()
if(NOT DEFINED INPUT_FILE OR NOT EXISTS "${INPUT_FILE}")
  message(FATAL_ERROR "INPUT_FILE must name the periodic SBM input deck")
endif()
if(NOT DEFINED TEST_ROOT)
  message(FATAL_ERROR "TEST_ROOT must name an isolated test output directory")
endif()

string(RANDOM LENGTH 12 ALPHABET 0123456789abcdef _run_id)
set(_run_dir "${TEST_ROOT}/sbm_float_restart_${_run_id}")
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
  if("${_stdout}\n${_stderr}" MATCHES "SBM begin-step time mismatch")
    message(FATAL_ERROR "${description} reported an SBM begin-step time mismatch")
  endif()
endfunction()

# Three 0.1 s SINGLE steps accumulate to a double time that is not representable
# in amrex::Real. Zero velocity keeps this advection-only fixture inside its
# donor-CFL bound while exercising the production checkpoint/restart path.
run_erf("Continuous four-step FLOAT SBM run"
  "max_step=4" "stop_time=0.5" "erf.fixed_dt=0.1"
  "prob.U_0=0.0" "prob.V_0=0.0" "prob.W_0=0.0"
  "erf.check_file=sbm_float_continuous_chk" "erf.check_int=1"
  "erf.plot_int_1=-1")

set(_continuous_step_three "${_run_dir}/sbm_float_continuous_chk00003")
set(_continuous_step_four "${_run_dir}/sbm_float_continuous_chk00004")
foreach(_checkpoint IN ITEMS "${_continuous_step_three}" "${_continuous_step_four}")
  if(NOT EXISTS "${_checkpoint}/Level_0/SBMSpectrum_H" OR
     NOT EXISTS "${_checkpoint}/Level_0/Cell_H")
    message(FATAL_ERROR "FLOAT run did not write the expected SBM checkpoint: ${_checkpoint}")
  endif()
endforeach()

file(READ "${_continuous_step_three}/Header" _step_three_header)
if(NOT _step_three_header MATCHES "0\\.3000000000000000[0-9]+")
  message(FATAL_ERROR
    "The step-three checkpoint does not preserve the expected non-float-representable double time; header:\n${_step_three_header}")
endif()

run_erf("FLOAT restart continuation from fractional-time checkpoint"
  "erf.restart=${_continuous_step_three}"
  "max_step=4" "stop_time=0.5" "erf.fixed_dt=0.1"
  "prob.U_0=0.0" "prob.V_0=0.0" "prob.W_0=0.0"
  "erf.check_file=sbm_float_restarted_chk" "erf.check_int=1"
  "erf.plot_int_1=-1")

set(_restarted_step_four "${_run_dir}/sbm_float_restarted_chk00004")
if(NOT EXISTS "${_restarted_step_four}/Level_0/SBMSpectrum_H" OR
   NOT EXISTS "${_restarted_step_four}/Level_0/Cell_H")
  message(FATAL_ERROR "FLOAT restart did not advance and write step 4: ${_restarted_step_four}")
endif()

execute_process(
  COMMAND "${COMPARE_EXECUTABLE}" "${_continuous_step_four}" "${_restarted_step_four}"
  WORKING_DIRECTORY "${_run_dir}"
  RESULT_VARIABLE _compare_result
  OUTPUT_VARIABLE _compare_stdout
  ERROR_VARIABLE _compare_stderr)
if(NOT "${_compare_result}" STREQUAL "0")
  message(FATAL_ERROR
    "FLOAT restart state comparison failed (exit ${_compare_result}).\n"
    "stdout:\n${_compare_stdout}\nstderr:\n${_compare_stderr}")
endif()
message(STATUS
  "FLOAT restart advanced beyond elapsed time 0.30000000000000004 s; spectrum and projected qc/qr match the continuous run.\n${_compare_stdout}")
