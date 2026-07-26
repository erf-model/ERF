cmake_minimum_required(VERSION 3.24)

if(NOT DEFINED MPIEXEC OR NOT DEFINED TEST_EXE OR NOT DEFINED INPUT)
  message(FATAL_ERROR "RunShocRegression.cmake requires MPIEXEC, TEST_EXE, and INPUT")
endif()

set(_run_command "${MPIEXEC}" "${MPIEXEC_NUMPROC_FLAG}" "${NRANKS}")
if(DEFINED MPIEXEC_PREFLAGS AND NOT "${MPIEXEC_PREFLAGS}" STREQUAL "")
  separate_arguments(_mpi_preflags UNIX_COMMAND "${MPIEXEC_PREFLAGS}")
  list(APPEND _run_command ${_mpi_preflags})
endif()
list(APPEND _run_command "${TEST_EXE}" "${INPUT}")
if(DEFINED RUNTIME_OPTIONS AND NOT "${RUNTIME_OPTIONS}" STREQUAL "")
  separate_arguments(_runtime_options UNIX_COMMAND "${RUNTIME_OPTIONS}")
  list(APPEND _run_command ${_runtime_options})
endif()

execute_process(
  COMMAND ${_run_command}
  WORKING_DIRECTORY "${WORKING_DIRECTORY}"
  OUTPUT_FILE "${LOG}"
  ERROR_FILE "${LOG}"
  RESULT_VARIABLE _run_result)
if(NOT _run_result EQUAL 0)
  message(FATAL_ERROR "SHOC simulation failed with exit code ${_run_result}; see ${LOG}")
endif()

set(_checker_command "${MPIEXEC}" "${MPIEXEC_NUMPROC_FLAG}" "1")
if(DEFINED MPIEXEC_PREFLAGS AND NOT "${MPIEXEC_PREFLAGS}" STREQUAL "")
  list(APPEND _checker_command ${_mpi_preflags})
endif()
list(APPEND _checker_command
  "${CHECKER}"
  --mode "${CHECK_MODE}"
  --initial "${INITIAL}"
  --midpoint "${MIDPOINT}"
  --final "${FINAL}")
execute_process(COMMAND ${_checker_command} RESULT_VARIABLE _checker_result)
if(NOT _checker_result EQUAL 0)
  message(FATAL_ERROR "SHOC property checker failed with exit code ${_checker_result}")
endif()

if(DEFINED SKIP_GOLD AND SKIP_GOLD)
  message(STATUS "SHOC property stage passed; gold comparison skipped")
  return()
endif()

if(NOT DEFINED GOLD_COMPARISON OR "${GOLD_COMPARISON}" STREQUAL "")
  set(GOLD_COMPARISON "fcompare")
endif()

set(_compare_command "${MPIEXEC}" "${MPIEXEC_NUMPROC_FLAG}" "1")
if(DEFINED MPIEXEC_PREFLAGS AND NOT "${MPIEXEC_PREFLAGS}" STREQUAL "")
  list(APPEND _compare_command ${_mpi_preflags})
endif()
if(GOLD_COMPARISON STREQUAL "fcompare")
  message(STATUS "SHOC property stage passed; using strict amrex_fcompare gold comparison")
  list(APPEND _compare_command
    "${FCOMPARE}"
    --abort_if_not_all_found
    -a -r "${RTOL}" --abs_tol "${ATOL}"
    "${GOLD}" "${FINAL}")
elseif(GOLD_COMPARISON STREQUAL "field_aware")
  if(NOT DEFINED GOLD_DIFFERENTIAL OR "${GOLD_DIFFERENTIAL}" STREQUAL "" OR
     NOT DEFINED GOLD_MODE OR "${GOLD_MODE}" STREQUAL "")
    message(FATAL_ERROR "field_aware gold comparison requires GOLD_DIFFERENTIAL and GOLD_MODE")
  endif()
  message(STATUS "SHOC property stage passed; using field-aware gold comparison mode ${GOLD_MODE}")
  list(APPEND _compare_command
    "${GOLD_DIFFERENTIAL}"
    --mode "${GOLD_MODE}"
    --gold "${GOLD}"
    --candidate "${FINAL}")
else()
  message(FATAL_ERROR "unknown SHOC gold comparison path '${GOLD_COMPARISON}'")
endif()
execute_process(COMMAND ${_compare_command} RESULT_VARIABLE _compare_result)
if(NOT _compare_result EQUAL 0)
  message(FATAL_ERROR "SHOC gold comparison failed with exit code ${_compare_result}")
endif()
