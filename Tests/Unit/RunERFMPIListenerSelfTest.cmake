# Run the intentionally failing MPI listener executable and validate that its
# failure is collected on the root rank. This is a CTest wrapper, not a death
# test, so MPI launchers remain in control of process setup and teardown.

foreach(required_var
    MPIEXEC_EXECUTABLE MPIEXEC_NUMPROC_FLAG MPIEXEC_PREFLAGS MPIEXEC_POSTFLAGS
    ERF_MPI_PREFLAGS TEST_EXECUTABLE)
  if(NOT DEFINED ${required_var})
    message(FATAL_ERROR "${required_var} is required")
  endif()
endforeach()

string(REPLACE "," " " _mpi_preflags_string "${ERF_MPI_PREFLAGS}")

separate_arguments(_mpiexec UNIX_COMMAND "${MPIEXEC_EXECUTABLE}")
separate_arguments(_numproc_flag UNIX_COMMAND "${MPIEXEC_NUMPROC_FLAG}")
separate_arguments(_preflags UNIX_COMMAND "${MPIEXEC_PREFLAGS}")
separate_arguments(_postflags UNIX_COMMAND "${MPIEXEC_POSTFLAGS}")
separate_arguments(_erf_mpi_preflags UNIX_COMMAND "${_mpi_preflags_string}")

execute_process(
  COMMAND ${_mpiexec} ${_numproc_flag} 2
          ${_preflags} ${_erf_mpi_preflags}
          "${TEST_EXECUTABLE}" ${_postflags}
  RESULT_VARIABLE _result
  OUTPUT_VARIABLE _stdout
  ERROR_VARIABLE _stderr
)
set(_output "${_stdout}\n${_stderr}")

if(_result EQUAL 0)
  message(FATAL_ERROR
    "MPI listener self-test unexpectedly passed; output:\n${_output}")
endif()
if(NOT _output MATCHES "NonRootFailureIsReported")
  message(FATAL_ERROR "Self-test output omitted the test name:\n${_output}")
endif()
if(NOT _output MATCHES "rank=1")
  message(FATAL_ERROR "Self-test output omitted the failing rank:\n${_output}")
endif()
if(NOT _output MATCHES "intentional non-root failure")
  message(FATAL_ERROR "Self-test output omitted the distinctive message:\n${_output}")
endif()

message(STATUS "MPI listener self-test validated non-root failure collection")
