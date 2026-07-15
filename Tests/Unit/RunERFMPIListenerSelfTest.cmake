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

string(REPLACE "," ";" _erf_mpi_preflags "${ERF_MPI_PREFLAGS}")

if(NOT DEFINED ERF_CROSSCOMPILING_EMULATOR OR
   "${ERF_CROSSCOMPILING_EMULATOR}" STREQUAL "" OR
   ERF_CROSSCOMPILING_EMULATOR MATCHES "-NOTFOUND$")
  set(_crosscompiling_emulator)
else()
  string(REPLACE "," ";" _crosscompiling_emulator
         "${ERF_CROSSCOMPILING_EMULATOR}")
endif()

separate_arguments(_mpiexec UNIX_COMMAND "${MPIEXEC_EXECUTABLE}")
separate_arguments(_numproc_flag UNIX_COMMAND "${MPIEXEC_NUMPROC_FLAG}")
separate_arguments(_preflags UNIX_COMMAND "${MPIEXEC_PREFLAGS}")
separate_arguments(_postflags UNIX_COMMAND "${MPIEXEC_POSTFLAGS}")

execute_process(
  COMMAND ${_mpiexec} ${_numproc_flag} 2
          ${_preflags} ${_erf_mpi_preflags}
          ${_crosscompiling_emulator}
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
