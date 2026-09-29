if(NOT DEFINED ERF_EXECUTABLE OR NOT EXISTS "${ERF_EXECUTABLE}")
  message(FATAL_ERROR "ERF_EXECUTABLE must name a built erf_exec")
endif()
if(NOT DEFINED INPUT_FILE OR NOT EXISTS "${INPUT_FILE}")
  message(FATAL_ERROR "INPUT_FILE must name the auxiliary inert-tracer input deck")
endif()
if(NOT DEFINED TEST_ROOT)
  message(FATAL_ERROR "TEST_ROOT must name an isolated test output directory")
endif()

set(_evidence_file "${TEST_ROOT}/auxiliary_m2_evidence.txt")
file(REMOVE "${_evidence_file}")
string(RANDOM LENGTH 12 ALPHABET 0123456789abcdef _run_id)
set(_run_dir "${TEST_ROOT}/auxiliary_inert_tracer_${_run_id}")
file(MAKE_DIRECTORY "${_run_dir}")

execute_process(
  COMMAND "${ERF_EXECUTABLE}" "${INPUT_FILE}"
  WORKING_DIRECTORY "${_run_dir}"
  RESULT_VARIABLE _result
  OUTPUT_VARIABLE _stdout
  ERROR_VARIABLE _stderr)
set(_run_output "${_stdout}\n${_stderr}")
if(NOT "${_result}" STREQUAL "0")
  message(FATAL_ERROR
    "Auxiliary inert-tracer ERF smoke failed (exit ${_result}).\nstdout:\n${_stdout}\nstderr:\n${_stderr}")
endif()

if(NOT _run_output MATCHES "AUX_M2_INITIALIZED schema=erf-auxiliary-inert-tracer-m2-v1 components=1")
  message(FATAL_ERROR "ERF run did not initialize the non-SBM M2 tracer.\n${_run_output}")
endif()
string(REGEX MATCHALL "AUX_M2_STAGE [^\r\n]*" _stage_records "${_run_output}")
list(LENGTH _stage_records _stage_count)
if(NOT _stage_count EQUAL 3)
  message(FATAL_ERROR "Expected exactly 3 compressible M2 stages; found ${_stage_count}.\n${_run_output}")
endif()

set(_max_stage_rate_delta 0.0)
foreach(_expected_stage RANGE 0 2)
  list(GET _stage_records ${_expected_stage} _stage_record)
  if(NOT _stage_record MATCHES "method=CompressibleRK3 stage=${_expected_stage}([^0-9]|$)")
    message(FATAL_ERROR "Unexpected M2 stage order: ${_stage_record}")
  endif()
  if(NOT _stage_record MATCHES "finite_state=true")
    message(FATAL_ERROR "M2 stage reported a nonfinite state: ${_stage_record}")
  endif()
  if(NOT _stage_record MATCHES "carrier_max=([0-9.eE+-]+)")
    message(FATAL_ERROR "M2 stage did not report a host carrier magnitude: ${_stage_record}")
  endif()
  set(_carrier_max "${CMAKE_MATCH_1}")
  if(NOT _carrier_max GREATER 0.0)
    message(FATAL_ERROR "M2 stage used a zero host carrier: ${_stage_record}")
  endif()
  if(NOT _stage_record MATCHES "rate_stage_delta=([0-9.eE+-]+)")
    message(FATAL_ERROR "M2 stage did not report its mapped rate change: ${_stage_record}")
  endif()
  set(_stage_rate_delta "${CMAKE_MATCH_1}")
  if(_stage_rate_delta GREATER _max_stage_rate_delta)
    set(_max_stage_rate_delta "${_stage_rate_delta}")
  endif()
endforeach()
if(NOT _max_stage_rate_delta GREATER 0.0)
  message(FATAL_ERROR "M2 mapped face rates did not vary over the host stages.\n${_run_output}")
endif()

string(REGEX MATCH "AUX_M2_LEDGER [^\r\n]*" _ledger_record "${_run_output}")
if(_ledger_record STREQUAL "")
  message(FATAL_ERROR "ERF run did not report a completed-step M2 ledger.\n${_run_output}")
endif()
if(NOT _ledger_record MATCHES "max_abs_residual=([0-9.eE+-]+)")
  message(FATAL_ERROR "M2 ledger omitted its absolute residual: ${_ledger_record}")
endif()
set(_max_abs_residual "${CMAKE_MATCH_1}")
if(NOT _ledger_record MATCHES "max_scaled_residual=([0-9.eE+-]+)")
  message(FATAL_ERROR "M2 ledger omitted its scaled residual: ${_ledger_record}")
endif()
set(_max_scaled_residual "${CMAKE_MATCH_1}")
if(NOT _ledger_record MATCHES "roundoff_tolerance=([0-9.eE+-]+)")
  message(FATAL_ERROR "M2 ledger omitted its roundoff tolerance: ${_ledger_record}")
endif()
set(_roundoff_tolerance "${CMAKE_MATCH_1}")
if(NOT _max_scaled_residual LESS_EQUAL _roundoff_tolerance)
  message(FATAL_ERROR "M2 ledger residual exceeds roundoff tolerance: ${_ledger_record}")
endif()

file(WRITE "${_evidence_file}"
  "stage_count=${_stage_count}\n"
  "max_stage_rate_delta=${_max_stage_rate_delta}\n"
  "max_abs_residual=${_max_abs_residual}\n"
  "max_scaled_residual=${_max_scaled_residual}\n"
  "roundoff_tolerance=${_roundoff_tolerance}\n"
  "${_ledger_record}\n")
message(STATUS "Auxiliary M2 smoke evidence: ${_evidence_file}")
message(STATUS "${_ledger_record}")
