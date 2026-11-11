cmake_minimum_required(VERSION 3.24)

include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")

# Write boundary planes with WRITER, read them with READER (erf.bndry_file = BndryFiles in the
# working directory) and compare the reader's plotfile PLTFILE with the gold directory GOLD
# with fcompare at the given tolerances. The writer and reader run on NRANKS ranks; fcompare on
# one.

foreach(_required MPIEXEC MPIEXEC_NUMPROC_FLAG NRANKS TEST_EXE FCOMPARE WRITER READER
                  WORKING_DIRECTORY PLTFILE GOLD RTOL ATOL)
  if(NOT DEFINED ${_required} OR "${${_required}}" STREQUAL "")
    message(FATAL_ERROR "RunBndryPlanesGold.cmake requires ${_required}")
  endif()
endforeach()

erf_mpi_launcher_command(_mpi_run
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunBndryPlanesGold.cmake")
erf_mpi_launcher_command(_mpi_one
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunBndryPlanesGold.cmake")

file(REMOVE_RECURSE "${WORKING_DIRECTORY}/BndryFiles" "${WORKING_DIRECTORY}/${PLTFILE}")

execute_process(
  COMMAND ${_mpi_run} "${TEST_EXE}" "${WRITER}" amrex.call_addr2line=0
  WORKING_DIRECTORY "${WORKING_DIRECTORY}"
  OUTPUT_FILE "${WORKING_DIRECTORY}/writer.log"
  ERROR_FILE "${WORKING_DIRECTORY}/writer.log"
  RESULT_VARIABLE _result)
if(NOT _result EQUAL 0)
  message(FATAL_ERROR "the plane writer failed (${_result}); see ${WORKING_DIRECTORY}/writer.log")
endif()

execute_process(
  COMMAND ${_mpi_run} "${TEST_EXE}" "${READER}" amrex.call_addr2line=0
  WORKING_DIRECTORY "${WORKING_DIRECTORY}"
  OUTPUT_FILE "${WORKING_DIRECTORY}/reader.log"
  ERROR_FILE "${WORKING_DIRECTORY}/reader.log"
  RESULT_VARIABLE _result)
if(NOT _result EQUAL 0)
  message(FATAL_ERROR "the plane reader failed (${_result}); see ${WORKING_DIRECTORY}/reader.log")
endif()
if(NOT EXISTS "${WORKING_DIRECTORY}/${PLTFILE}")
  message(FATAL_ERROR "the plane reader wrote no ${PLTFILE}; see ${WORKING_DIRECTORY}/reader.log")
endif()

execute_process(
  COMMAND ${_mpi_one} "${FCOMPARE}" --abort_if_not_all_found -a -r ${RTOL} --abs_tol ${ATOL}
          "${GOLD}" "${WORKING_DIRECTORY}/${PLTFILE}"
  WORKING_DIRECTORY "${WORKING_DIRECTORY}"
  OUTPUT_FILE "${WORKING_DIRECTORY}/fcompare.log"
  ERROR_FILE "${WORKING_DIRECTORY}/fcompare.log"
  RESULT_VARIABLE _result)
if(NOT _result EQUAL 0)
  message(FATAL_ERROR "the reader's ${PLTFILE} differs from the gold ${GOLD}; see ${WORKING_DIRECTORY}/fcompare.log")
endif()
message(STATUS "the reader's ${PLTFILE} matches the gold")
