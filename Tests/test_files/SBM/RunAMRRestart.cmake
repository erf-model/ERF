# PRODUCTION runs the normal executable, which prints no qualification
# inventories, without a launcher unless MPIEXEC is given.
set(_required_inputs ERF_EXECUTABLE INPUT_FILE TEST_ROOT RESTART_STEP FINAL_STEP MODE)
if(NOT PRODUCTION)
  list(APPEND _required_inputs MPIEXEC MPIEXEC_NUMPROC_FLAG MPI_RANKS EXPECTED_COMPONENTS)
endif()
foreach(_required IN LISTS _required_inputs)
  if(NOT DEFINED ${_required})
    message(FATAL_ERROR "${_required} is required for the SBM AMR restart qualification")
  endif()
endforeach()
if(NOT DEFINED MPI_RANKS)
  set(MPI_RANKS 1)
endif()

string(RANDOM LENGTH 12 ALPHABET 0123456789abcdef _run_id)
set(_run_dir "${TEST_ROOT}/sbm_amr_restart_${_run_id}")
set(_continuous_dir "${_run_dir}/continuous")
set(_restarted_dir "${_run_dir}/restarted")
file(MAKE_DIRECTORY "${_continuous_dir}" "${_restarted_dir}")

set(_host_args)
if(DEFINED ANELASTIC)
  set(_host_args "erf.compressibility_type=Anelastic" "erf.timestepping_type=${ANELASTIC}")
endif()
# RUN_ARGS reaches both runs and RESTART_ARGS the restart alone, each "|"-separated.
if(DEFINED RUN_ARGS)
  string(REPLACE "|" ";" _run_args "${RUN_ARGS}")
  list(APPEND _host_args ${_run_args})
endif()
set(_restart_args)
if(DEFINED RESTART_ARGS)
  string(REPLACE "|" ";" _restart_args "${RESTART_ARGS}")
endif()
set(_launcher)
if(DEFINED MPIEXEC AND NOT "${MPIEXEC}" STREQUAL "")
  set(_launcher "${MPIEXEC}" "${MPIEXEC_NUMPROC_FLAG}" "${MPI_RANKS}")
endif()

function(run_erf description working_dir log_name result_var)
  execute_process(
    COMMAND ${_launcher}
            "${ERF_EXECUTABLE}" "${INPUT_FILE}" "max_step=${FINAL_STEP}"
            "erf.check_int=1" "erf.check_file=chk" ${_host_args} ${ARGN}
    WORKING_DIRECTORY "${working_dir}"
    TIMEOUT 240
    RESULT_VARIABLE _result
    OUTPUT_VARIABLE _stdout
    ERROR_VARIABLE _stderr)
  set(_log "${_stdout}\n${_stderr}")
  file(WRITE "${working_dir}/${log_name}" "${_log}")
  set(${result_var} "${_result}" PARENT_SCOPE)
  set(_run_log "${_log}" PARENT_SCOPE)
endfunction()

run_erf("Continuous run" "${_continuous_dir}" continuous.log _continuous_result)
if(NOT "${_continuous_result}" STREQUAL "0")
  message(FATAL_ERROR "Continuous SBM AMR run failed (${_continuous_result}); see ${_continuous_dir}/continuous.log")
endif()
set(_continuous_log "${_run_log}")
string(LENGTH "${RESTART_STEP}" _restart_digits)
string(REPEAT "0" 5 _zeros)
math(EXPR _pad "5 - ${_restart_digits}")
string(SUBSTRING "${_zeros}" 0 ${_pad} _restart_pad)
set(_checkpoint "${_continuous_dir}/chk${_restart_pad}${RESTART_STEP}")
string(LENGTH "${FINAL_STEP}" _final_digits)
math(EXPR _pad "5 - ${_final_digits}")
string(SUBSTRING "${_zeros}" 0 ${_pad} _final_pad)
set(_final_name "chk${_final_pad}${FINAL_STEP}")
if(NOT EXISTS "${_checkpoint}/SBM_Schema")
  message(FATAL_ERROR "Continuous run did not write ${_checkpoint}")
endif()

if(MODE STREQUAL "missing_fine_spectrum")
  file(GLOB _fine_spectrum "${_checkpoint}/Level_1/SBMSpectrum*")
  if(NOT _fine_spectrum)
    message(FATAL_ERROR "The checkpoint has no level-1 spectrum to remove: ${_checkpoint}")
  endif()
  file(REMOVE ${_fine_spectrum})
elseif(MODE STREQUAL "schema_mismatch")
  # Same component count, different regrid policy identity.
  file(READ "${_checkpoint}/SBM_Schema" _schema)
  string(REPLACE "regrid=retained-copy-plus-carrier-relative-prolongation-v1"
                 "regrid=altered-policy-v9" _altered "${_schema}")
  if(_altered STREQUAL _schema)
    message(FATAL_ERROR "The checkpoint schema has no regrid policy to alter:\n${_schema}")
  endif()
  file(WRITE "${_checkpoint}/SBM_Schema" "${_altered}")
elseif(MODE STREQUAL "corruption")
  foreach(_required IN ITEMS CORRUPTER CORRUPTION CORRUPT_LEVEL)
    if(NOT DEFINED ${_required})
      message(FATAL_ERROR "${_required} is required for a corruption restart test")
    endif()
  endforeach()
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env "SBM_CHECKPOINT_DIR=${_checkpoint}"
            "SBM_CHECKPOINT_LEVEL=${CORRUPT_LEVEL}" "SBM_CORRUPTION=${CORRUPTION}"
            "${CORRUPTER}"
    WORKING_DIRECTORY "${_run_dir}"
    TIMEOUT 30
    RESULT_VARIABLE _corrupt_result
    OUTPUT_VARIABLE _corrupt_stdout
    ERROR_VARIABLE _corrupt_stderr)
  if(NOT "${_corrupt_result}" STREQUAL "0")
    message(FATAL_ERROR
      "Checkpoint corruption failed (${_corrupt_result}).\n${_corrupt_stdout}\n${_corrupt_stderr}")
  endif()
elseif(NOT MODE STREQUAL "equivalence" AND NOT MODE STREQUAL "changed_restart_inputs")
  message(FATAL_ERROR "Unknown SBM AMR restart mode '${MODE}'")
endif()

run_erf("Restarted run" "${_restarted_dir}" restarted.log _restart_result
        "erf.restart=${_checkpoint}" "amrex.call_addr2line=0" ${_restart_args})
set(_restarted_log "${_run_log}")

if(NOT MODE STREQUAL "equivalence")
  if(_restart_result MATCHES "[Tt]imeout" OR "${_restart_result}" STREQUAL "0")
    message(FATAL_ERROR
      "The ${MODE} restart was not rejected (${_restart_result}); see ${_restarted_dir}/restarted.log")
  endif()
  if(NOT _restarted_log MATCHES "${EXPECTED_DIAGNOSTIC}")
    message(FATAL_ERROR
      "The ${MODE} restart failed without '${EXPECTED_DIAGNOSTIC}'; see ${_restarted_dir}/restarted.log")
  endif()
  if(_restarted_log MATCHES "SBM_AMR_QUAL step=")
    message(FATAL_ERROR "A step was advanced from the rejected ${MODE} checkpoint")
  endif()
  message(STATUS "SBM AMR ${MODE} restart rejected: ${EXPECTED_DIAGNOSTIC}")
  return()
endif()

if(NOT "${_restart_result}" STREQUAL "0")
  message(FATAL_ERROR "Restarted SBM AMR run failed (${_restart_result}); see ${_restarted_dir}/restarted.log")
endif()

# Composite continuity: the restored inventory is the continuous run's
# inventory at the checkpoint step, and both runs end with the same one.
if(PRODUCTION)
  if("${_continuous_log}${_restarted_log}" MATCHES "SBM_AMR_QUAL")
    message(FATAL_ERROR "The production restart test ran a qualification executable")
  endif()
  if(NOT DEFINED COMPARE_EXECUTABLE)
    message(FATAL_ERROR "A production restart test requires COMPARE_EXECUTABLE")
  endif()
else()
string(REGEX MATCHALL "SBM_AMR_QUAL_INVENTORY component=[0-9]+ initial=[^ ]+ final=[^ ]+"
       _continuous_inventory "${_continuous_log}")
string(REGEX MATCHALL "SBM_AMR_QUAL_INVENTORY component=[0-9]+ initial=[^ ]+ final=[^ ]+"
       _restarted_inventory "${_restarted_log}")
string(REGEX MATCHALL "SBM_AMR_QUAL_INITIAL_INVENTORY component=[0-9]+ value=[^ \n]+"
       _restored_inventory "${_restarted_log}")
list(LENGTH _restored_inventory _restored_count)
if(NOT _restored_count EQUAL EXPECTED_COMPONENTS)
  message(FATAL_ERROR "The restarted run recorded ${_restored_count} restored inventories")
endif()
math(EXPR _restart_first "(${RESTART_STEP} - 1) * ${EXPECTED_COMPONENTS}")
math(EXPR _continuous_last_first "(${FINAL_STEP} - 1) * ${EXPECTED_COMPONENTS}")
math(EXPR _restarted_last_first "(${FINAL_STEP} - ${RESTART_STEP} - 1) * ${EXPECTED_COMPONENTS}")
math(EXPR _last_component "${EXPECTED_COMPONENTS} - 1")
foreach(_component RANGE 0 ${_last_component})
  math(EXPR _at_restart "${_restart_first} + ${_component}")
  math(EXPR _continuous_last "${_continuous_last_first} + ${_component}")
  math(EXPR _restarted_last "${_restarted_last_first} + ${_component}")
  list(GET _continuous_inventory ${_at_restart} _line)
  string(REGEX REPLACE ".* final=" "" _checkpoint_value "${_line}")
  list(GET _restored_inventory ${_component} _line)
  string(REGEX REPLACE ".* value=" "" _restored_value "${_line}")
  if(NOT _restored_value STREQUAL _checkpoint_value)
    message(FATAL_ERROR
      "Component ${_component} restored inventory ${_restored_value} differs from the checkpoint-step inventory ${_checkpoint_value}")
  endif()
  list(GET _continuous_inventory ${_continuous_last} _line)
  string(REGEX REPLACE ".* final=" "" _continuous_value "${_line}")
  list(GET _restarted_inventory ${_restarted_last} _line)
  string(REGEX REPLACE ".* final=" "" _restarted_value "${_line}")
  if(NOT _restarted_value STREQUAL _continuous_value)
    message(FATAL_ERROR
      "Component ${_component} final inventory ${_restarted_value} differs from the continuous run's ${_continuous_value}")
  endif()
endforeach()
endif()

if(DEFINED COMPARE_EXECUTABLE)
  execute_process(
    COMMAND "${COMPARE_EXECUTABLE}" "${_continuous_dir}/${_final_name}"
            "${_restarted_dir}/${_final_name}"
    WORKING_DIRECTORY "${_run_dir}"
    RESULT_VARIABLE _compare_result
    OUTPUT_VARIABLE _compare_stdout
    ERROR_VARIABLE _compare_stderr)
  if(NOT "${_compare_result}" STREQUAL "0")
    message(FATAL_ERROR
      "Continuous and restarted checkpoints differ (${_compare_result}).\n${_compare_stdout}\n${_compare_stderr}")
  endif()
  # A same-decomposition restart is bitwise, on every level the final
  # checkpoint holds.
  string(REGEX MATCHALL "max_abs_difference=[^ \n]+" _differences "${_compare_stdout}")
  list(REMOVE_DUPLICATES _differences)
  if(NOT _differences STREQUAL "max_abs_difference=0")
    message(FATAL_ERROR "The restarted checkpoint is not bitwise identical:\n${_compare_stdout}")
  endif()
  file(GLOB _final_levels LIST_DIRECTORIES true "${_continuous_dir}/${_final_name}/Level_*")
  list(LENGTH _final_levels _final_level_count)
  math(EXPR _finest "${_final_level_count} - 1")
  if(NOT _compare_stdout MATCHES "level ${_finest} SBMSpectrum")
    message(FATAL_ERROR "The comparison did not cover level ${_finest}:\n${_compare_stdout}")
  endif()
  if(DEFINED EXPECTED_FINAL_LEVELS AND NOT _final_level_count EQUAL EXPECTED_FINAL_LEVELS)
    message(FATAL_ERROR
      "The final checkpoint has ${_final_level_count} levels, not ${EXPECTED_FINAL_LEVELS}")
  endif()
endif()
# A run with physical faces carries its boundary ledger across the restart.
set(_continuous_ledger "${_continuous_dir}/${_final_name}/SBM_BoundaryLedger")
set(_restarted_ledger "${_restarted_dir}/${_final_name}/SBM_BoundaryLedger")
if(EXPECT_BOUNDARY_LEDGER)
  if(NOT EXISTS "${_continuous_ledger}" OR NOT EXISTS "${_restarted_ledger}")
    message(FATAL_ERROR "A final checkpoint has no SBM boundary ledger")
  endif()
  file(READ "${_continuous_ledger}" _continuous_ledger_text)
  file(READ "${_restarted_ledger}" _restarted_ledger_text)
  if(NOT _continuous_ledger_text STREQUAL _restarted_ledger_text)
    message(FATAL_ERROR "The restarted boundary ledger differs from the continuous run's:\n"
                        "${_continuous_ledger_text}\n---\n${_restarted_ledger_text}")
  endif()
elseif(EXISTS "${_continuous_ledger}")
  message(FATAL_ERROR "A triply periodic run wrote an SBM boundary ledger")
endif()
message(STATUS
  "SBM AMR restart from step ${RESTART_STEP} reproduces step ${FINAL_STEP}: ranks=${MPI_RANKS}; ${_run_dir}")
