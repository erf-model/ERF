# Run an SBM deck with erf.sbm_budget_check = true and check its boundary
# budget lines. EXPECT_TRANSFER=zero requires every face transfer to vanish;
# EXPECT_TRANSFER=nonzero requires some face transfer. RUN_ARGS holds extra
# inputs overrides separated by "|". EXPECTED_LEVELS checks the finest
# level count the budget reported.
foreach(_required IN ITEMS ERF_EXECUTABLE INPUT_FILE TEST_ROOT EXPECT_TRANSFER)
  if(NOT DEFINED ${_required})
    message(FATAL_ERROR "${_required} is required for the SBM boundary budget test")
  endif()
endforeach()

string(RANDOM LENGTH 12 ALPHABET 0123456789abcdef _run_id)
set(_run_dir "${TEST_ROOT}/sbm_boundary_budget_${_run_id}")
file(MAKE_DIRECTORY "${_run_dir}")
set(_run_args)
if(DEFINED RUN_ARGS)
  string(REPLACE "|" ";" _run_args "${RUN_ARGS}")
endif()
execute_process(
  COMMAND "${ERF_EXECUTABLE}" "${INPUT_FILE}" ${_run_args} "amrex.call_addr2line=0"
  WORKING_DIRECTORY "${_run_dir}"
  TIMEOUT 240
  RESULT_VARIABLE _result
  OUTPUT_VARIABLE _stdout
  ERROR_VARIABLE _stderr)
if(NOT "${_result}" STREQUAL "0")
  message(FATAL_ERROR
    "SBM boundary budget run failed (${_result}).\nstdout:\n${_stdout}\nstderr:\n${_stderr}")
endif()

string(REGEX MATCHALL "SBM_BUDGET step=[0-9]+ [^\n]*" _lines "${_stdout}")
list(LENGTH _lines _count)
if(_count EQUAL 0)
  message(FATAL_ERROR "The run printed no SBM_BUDGET lines.\nstdout:\n${_stdout}")
endif()
if(DEFINED EXPECTED_LEVELS AND
   NOT _stdout MATCHES "SBM_BUDGET_REFERENCE step=0 levels=${EXPECTED_LEVELS} ")
  message(FATAL_ERROR "The budget reference did not see ${EXPECTED_LEVELS} levels.\nstdout:\n${_stdout}")
endif()
set(_nonzero FALSE)
foreach(_line IN LISTS _lines)
  string(REGEX MATCH "outward_transfer=([^ ]+)" _match "${_line}")
  if(NOT CMAKE_MATCH_1 STREQUAL "0")
    set(_nonzero TRUE)
  endif()
endforeach()
if(EXPECT_TRANSFER STREQUAL "zero" AND _nonzero)
  message(FATAL_ERROR "A wall-bounded run recorded a boundary transfer:\n${_lines}")
endif()
if(EXPECT_TRANSFER STREQUAL "nonzero" AND NOT _nonzero)
  message(FATAL_ERROR "The run recorded no boundary transfer:\n${_lines}")
endif()
if(DEFINED EXPECTED_STDOUT AND NOT _stdout MATCHES "${EXPECTED_STDOUT}")
  message(FATAL_ERROR "The run did not print '${EXPECTED_STDOUT}'.\nstdout:\n${_stdout}")
endif()
list(GET _lines -1 _last)
message(STATUS "SBM boundary budget closed over ${_count} component-steps; last: ${_last}")
