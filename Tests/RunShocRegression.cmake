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
if(NOT _run_result EQUAL 0 AND
   NOT (_run_result EQUAL 143 AND EXISTS "${FINAL}"))
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
if(DEFINED EXPECT_CHECKER_FAILURE AND EXPECT_CHECKER_FAILURE)
  if(_checker_result EQUAL 0)
    message(FATAL_ERROR "SHOC negative control unexpectedly passed the property checker")
  endif()
  message(STATUS "SHOC negative control correctly failed the property checker")
  return()
endif()
if(NOT _checker_result EQUAL 0)
  message(FATAL_ERROR "SHOC property checker failed with exit code ${_checker_result}")
endif()

set(_compare_command "${MPIEXEC}" "${MPIEXEC_NUMPROC_FLAG}" "1")
if(DEFINED MPIEXEC_PREFLAGS AND NOT "${MPIEXEC_PREFLAGS}" STREQUAL "")
  list(APPEND _compare_command ${_mpi_preflags})
endif()
list(APPEND _compare_command
  "${FCOMPARE}"
  --abort_if_not_all_found
  -a -r "${RTOL}" --abs_tol "${ATOL}"
  "${GOLD}" "${FINAL}")
execute_process(COMMAND ${_compare_command} RESULT_VARIABLE _compare_result)
if(NOT _compare_result EQUAL 0)
  message(FATAL_ERROR "SHOC gold comparison failed with exit code ${_compare_result}")
endif()
