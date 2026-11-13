# Run an SBM particle-mixing deck and check what the mixing changes.
# The run must succeed and its final SBM_Schema carry the mixing lines.
# Optional checks, each list of inputs separated by "|":
#   REFERENCE_ARGS  a run with these overrides (no mixing) must differ;
#   RESTART_STEP    a restart from that checkpoint must reproduce the final one;
#   MISMATCH_ARGS   a restart with these overrides must fail the schema check.
foreach(_required IN ITEMS ERF_EXECUTABLE COMPARE_EXECUTABLE INPUT_FILE TEST_ROOT FINAL_STEP)
  if(NOT DEFINED ${_required})
    message(FATAL_ERROR "${_required} is required for the SBM particle-mixing test")
  endif()
endforeach()

string(RANDOM LENGTH 12 ALPHABET 0123456789abcdef _run_id)
set(_root "${TEST_ROOT}/sbm_particle_mixing_${_run_id}")

function(_split out text)
  set(_list)
  if(NOT "${text}" STREQUAL "")
    string(REPLACE "|" ";" _list "${text}")
  endif()
  set(${out} "${_list}" PARENT_SCOPE)
endfunction()

function(_run dir result_var stdout_var)
  file(MAKE_DIRECTORY "${dir}")
  execute_process(
    COMMAND "${ERF_EXECUTABLE}" "${INPUT_FILE}" ${ARGN} "amrex.call_addr2line=0"
    WORKING_DIRECTORY "${dir}"
    TIMEOUT 240
    RESULT_VARIABLE _result
    OUTPUT_VARIABLE _stdout
    ERROR_VARIABLE _stderr)
  set(${result_var} "${_result}" PARENT_SCOPE)
  set(${stdout_var} "${_stdout}${_stderr}" PARENT_SCOPE)
endfunction()

function(_compare first second result_var output_var)
  execute_process(
    COMMAND "${COMPARE_EXECUTABLE}" "${first}" "${second}"
    TIMEOUT 120
    RESULT_VARIABLE _result
    OUTPUT_VARIABLE _stdout
    ERROR_VARIABLE _stderr)
  set(${result_var} "${_result}" PARENT_SCOPE)
  set(${output_var} "${_stdout}${_stderr}" PARENT_SCOPE)
endfunction()

_split(_run_args "${RUN_ARGS}")
string(LENGTH "${FINAL_STEP}" _digits)
math(EXPR _pad "5 - ${_digits}")
string(REPEAT "0" ${_pad} _zeros)
set(_final "chk${_zeros}${FINAL_STEP}")

_run("${_root}/mixing" _result _stdout ${_run_args})
if(NOT "${_result}" STREQUAL "0")
  message(FATAL_ERROR "The particle-mixing run failed (${_result}):\n${_stdout}")
endif()
if(NOT _stdout MATCHES "SBM particle mixing: particle-mixing=")
  message(FATAL_ERROR "The run did not report its particle mixing policy:\n${_stdout}")
endif()
string(REGEX MATCHALL "SBM_BUDGET step=[0-9]+ [^\n]*" _budget "${_stdout}")
if(NOT _budget)
  message(FATAL_ERROR "The particle-mixing run printed no SBM_BUDGET lines:\n${_stdout}")
endif()
foreach(_line IN LISTS _budget)
  string(REGEX MATCH "outward_transfer=([^ ]+)" _match "${_line}")
  if(NOT CMAKE_MATCH_1 STREQUAL "0")
    message(FATAL_ERROR "Particle mixing moved material through a wall:\n${_line}")
  endif()
endforeach()
set(_mixing_checkpoint "${_root}/mixing/${_final}")
file(READ "${_mixing_checkpoint}/SBM_Schema" _schema)
if(NOT _schema MATCHES "\nparticle-mixing=[^\n]+\ndiffusion=two-point-low-host-high-group-fct-v1\n")
  message(FATAL_ERROR "The SBM restart schema lacks the particle-mixing lines:\n${_schema}")
endif()

if(DEFINED REFERENCE_ARGS)
  _split(_reference_args "${REFERENCE_ARGS}")
  _run("${_root}/reference" _result _stdout ${_run_args} ${_reference_args})
  if(NOT "${_result}" STREQUAL "0")
    message(FATAL_ERROR "The reference run failed (${_result}):\n${_stdout}")
  endif()
  _compare("${_mixing_checkpoint}" "${_root}/reference/${_final}" _result _output)
  if("${_result}" STREQUAL "0")
    message(FATAL_ERROR "Particle mixing left the spectrum unchanged:\n${_output}")
  endif()
endif()

if(DEFINED RESTART_STEP)
  string(LENGTH "${RESTART_STEP}" _digits)
  math(EXPR _pad "5 - ${_digits}")
  string(REPEAT "0" ${_pad} _zeros)
  _run("${_root}/restart" _result _stdout ${_run_args}
       "amr.restart=${_root}/mixing/chk${_zeros}${RESTART_STEP}")
  if(NOT "${_result}" STREQUAL "0")
    message(FATAL_ERROR "The particle-mixing restart failed (${_result}):\n${_stdout}")
  endif()
  _compare("${_mixing_checkpoint}" "${_root}/restart/${_final}" _result _output)
  if(NOT "${_result}" STREQUAL "0")
    message(FATAL_ERROR "The restart did not reproduce the uninterrupted run:\n${_output}")
  endif()
  string(REGEX MATCHALL "max_abs_difference=[^ ]+" _differences "${_output}")
  list(REMOVE_DUPLICATES _differences)
  if(NOT _differences STREQUAL "max_abs_difference=0")
    message(FATAL_ERROR "The restart is not bitwise equal to the uninterrupted run:\n${_output}")
  endif()
endif()

if(DEFINED MISMATCH_ARGS)
  _split(_mismatch_args "${MISMATCH_ARGS}")
  _run("${_root}/mismatch" _result _stdout ${_run_args} ${_mismatch_args}
       "amr.restart=${_mixing_checkpoint}")
  if("${_result}" STREQUAL "0" OR
     NOT _stdout MATCHES "SBM restart schema does not exactly match")
    message(FATAL_ERROR "A restart with another particle mixing was not rejected:\n${_stdout}")
  endif()
endif()
message(STATUS "SBM particle mixing checks passed")
