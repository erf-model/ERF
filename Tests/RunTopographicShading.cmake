foreach(required TEST_EXE PYTHON_EXE CASE_DIR WORKING_DIRECTORY)
  if(NOT DEFINED ${required} OR "${${required}}" STREQUAL "")
    message(FATAL_ERROR "RunTopographicShading needs ${required}")
  endif()
endforeach()

# This is the test's own scratch directory, beneath the CMake build tree.
set(run "${WORKING_DIRECTORY}/run")
file(REMOVE_RECURSE "${run}")
execute_process(COMMAND "${PYTHON_EXE}" "${CASE_DIR}/prepare.py" "${run}"
  RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "Hill preparation failed:\n${output}\n${error}")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
erf_mpi_launcher_command(launcher
  LAUNCHER "${MPIEXEC}" NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
  NRANKS "${NRANKS}" PREFLAGS "${MPIEXEC_PREFLAGS}"
  CONTEXT "RunTopographicShading.cmake")
execute_process(COMMAND ${launcher} "${TEST_EXE}" inputs
  WORKING_DIRECTORY "${run}" RESULT_VARIABLE result
  OUTPUT_FILE "${run}/simulation.log" ERROR_FILE "${run}/simulation.log"
  TIMEOUT 600)
if(NOT result EQUAL 0 OR NOT EXISTS "${run}/surface00003_d01.nc")
  message(FATAL_ERROR "Hill run failed (${result}); see ${run}/simulation.log")
endif()

execute_process(COMMAND "${PYTHON_EXE}" "${CASE_DIR}/check.py"
  "${run}/surface00001_d01.nc" --held "${run}/surface00002_d01.nc"
  RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
message(STATUS "${output}")
if(NOT result EQUAL 0)
  message(FATAL_ERROR "Hill check failed:\n${output}\n${error}")
endif()
execute_process(COMMAND "${PYTHON_EXE}" "${CASE_DIR}/check.py"
  "${run}/surface00003_d01.nc"
  RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
message(STATUS "${output}")
if(NOT result EQUAL 0)
  message(FATAL_ERROR "Hill check after radiation update failed:\n${output}\n${error}")
endif()
