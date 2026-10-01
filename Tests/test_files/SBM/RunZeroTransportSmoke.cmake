if(NOT DEFINED ERF_EXECUTABLE OR NOT EXISTS "${ERF_EXECUTABLE}")
  message(FATAL_ERROR "ERF_EXECUTABLE must name a built erf_exec")
endif()
if(NOT DEFINED INPUT_FILE OR NOT EXISTS "${INPUT_FILE}")
  message(FATAL_ERROR "INPUT_FILE must name the SBM zero-transport input deck")
endif()
if(NOT DEFINED TEST_ROOT)
  message(FATAL_ERROR "TEST_ROOT must name an isolated test output directory")
endif()
if(NOT DEFINED EXPECTED_SPECTRUM_VALUES)
  set(EXPECTED_SPECTRUM_VALUES
    "9.99999999999999955e-07,1.99999999999999991e-06,3.00000000000000008e-06,3.99999999999999982e-06,")
endif()
if(NOT DEFINED EXPECTED_CORE_VALUES)
  set(EXPECTED_CORE_VALUES
    "1.00000000000000000e+00,2.99999999999999886e+02,0.00000000000000000e+00,0.00000000000000000e+00,0.00000000000000000e+00,3.00000000000000008e-06,6.99999999999999989e-06,")
endif()

string(RANDOM LENGTH 12 ALPHABET 0123456789abcdef _run_id)
set(_run_dir "${TEST_ROOT}/sbm_zero_transport_${_run_id}")
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

function(assert_checkpoint checkpoint)
  set(_checkpoint_dir "${_run_dir}/${checkpoint}")
  foreach(_file IN ITEMS
      "${_checkpoint_dir}/SBM_Schema"
    "${_checkpoint_dir}/Level_0/SBMSpectrum_H"
      "${_checkpoint_dir}/Level_0/Cell_H")
    if(NOT EXISTS "${_file}")
      message(FATAL_ERROR "SBM smoke checkpoint is missing ${_file}")
    endif()
  endforeach()

  file(READ "${_checkpoint_dir}/SBM_Schema" _schema)
  if(NOT _schema MATCHES "transport=zero-transport-fixture-v1")
    message(FATAL_ERROR "Checkpoint has the wrong SBM fixture schema: ${_schema}")
  endif()

  file(READ "${_checkpoint_dir}/Level_0/SBMSpectrum_H" _spectrum_header)
  string(FIND "${_spectrum_header}" "${EXPECTED_SPECTRUM_VALUES}" _spectrum_found)
  if(_spectrum_found EQUAL -1)
    message(FATAL_ERROR "SBM checkpoint spectrum differs from the manufactured state: ${_checkpoint_dir}")
  endif()

  file(READ "${_checkpoint_dir}/Level_0/Cell_H" _cell_header)
  string(FIND "${_cell_header}" "${EXPECTED_CORE_VALUES}" _projection_found)
  if(_projection_found EQUAL -1)
    message(FATAL_ERROR "SBM compact qc/qr checkpoint values do not match the fixed projection: ${_checkpoint_dir}")
  endif()
endfunction()

run_erf("Initial zero-transport run")
assert_checkpoint("sbm_smoke_chk00001")

run_erf("Exact-schema SBM restart"
  "erf.restart=${_run_dir}/sbm_smoke_chk00001"
  "erf.check_file=sbm_restart_chk"
  "erf.plot_int_1=-1"
  "max_step=2"
  "stop_time=2.e-4")
assert_checkpoint("sbm_restart_chk00002")
