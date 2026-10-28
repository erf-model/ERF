cmake_minimum_required(VERSION 3.24)

include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")

# Write boundary planes from a moist, periodic run, read them on an inflow face
# of a dry run that starts at rest, and require one field of the reader's last
# plotfile to have taken the planes' values: its largest value must lie in
# [MIN, MAX]. The writer's planes carry u = 10 m/s and qv = 5 g/kg; the reader
# starts with u = 0 and qv = 0, so a field the reader never fills from the
# planes stays near zero and fails MIN. The bounds leave room for the start-up
# transient: the 10 m/s inflow into air at rest sends a pressure wave in, and
# the largest u behind it is about 11 m/s.
#
# READ_VARS is the reader's erf.bndry_input_var_names and FIELD the plotfile
# variable checked.

foreach(_required MPIEXEC MPIEXEC_NUMPROC_FLAG NRANKS TEST_EXE FEXTREMA WRITER READER
                  WORKING_DIRECTORY READ_VARS FIELD MIN MAX)
  if(NOT DEFINED ${_required} OR "${${_required}}" STREQUAL "")
    message(FATAL_ERROR "RunBndryPlanesIngest.cmake requires ${_required}")
  endif()
endforeach()
if(NOT EXISTS "${FEXTREMA}")
  message(FATAL_ERROR "RunBndryPlanesIngest.cmake: FEXTREMA=${FEXTREMA} does not exist; "
                      "set the FEXTREMA_EXE cache entry to the amrex_fextrema built "
                      "alongside amrex_fcompare")
endif()

erf_mpi_launcher_command(_mpi_run
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunBndryPlanesIngest.cmake")
erf_mpi_launcher_command(_mpi_one
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunBndryPlanesIngest.cmake")

file(REMOVE_RECURSE "${WORKING_DIRECTORY}/BndryFiles")
file(GLOB _stale "${WORKING_DIRECTORY}/reader_plt*")
if(_stale)
  file(REMOVE_RECURSE ${_stale})
endif()

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
  COMMAND ${_mpi_run} "${TEST_EXE}" "${READER}" "erf.bndry_input_var_names=${READ_VARS}"
          erf.plot_file_1=reader_plt amrex.call_addr2line=0
  WORKING_DIRECTORY "${WORKING_DIRECTORY}"
  OUTPUT_FILE "${WORKING_DIRECTORY}/reader.log"
  ERROR_FILE "${WORKING_DIRECTORY}/reader.log"
  RESULT_VARIABLE _result)
if(NOT _result EQUAL 0)
  message(FATAL_ERROR "the plane reader failed (${_result}); see ${WORKING_DIRECTORY}/reader.log")
endif()

file(GLOB _plts LIST_DIRECTORIES true "${WORKING_DIRECTORY}/reader_plt*")
list(SORT _plts)
list(LENGTH _plts _nplt)
if(_nplt EQUAL 0)
  message(FATAL_ERROR "the plane reader wrote no plotfile; see ${WORKING_DIRECTORY}/reader.log")
endif()
list(GET _plts -1 _plt)

execute_process(
  COMMAND ${_mpi_one} "${FEXTREMA}" -v "${FIELD}" "${_plt}"
  WORKING_DIRECTORY "${WORKING_DIRECTORY}"
  OUTPUT_VARIABLE _extrema
  ERROR_VARIABLE _extrema_err
  RESULT_VARIABLE _result)
file(WRITE "${WORKING_DIRECTORY}/fextrema.log" "${_extrema}${_extrema_err}")
# Some launchers label task output with a per-line "<rank>: " prefix; strip it.
string(REGEX REPLACE "\n[0-9]+:[ \t]*" "\n " _extrema "\n${_extrema}")
# fextrema prints one line per variable: " name   min   max"
if(NOT _result EQUAL 0 OR
   NOT "${_extrema}" MATCHES "\n ${FIELD}[ \t]+([^ \t\n]+)[ \t]+([^ \t\n]+)[ \t]*\n")
  message(FATAL_ERROR "fextrema could not read ${FIELD} from ${_plt}; "
                      "see ${WORKING_DIRECTORY}/fextrema.log")
endif()
set(_max "${CMAKE_MATCH_2}")
# Written so that a value that is not a number (nan, inf, text) fails: CMake's
# LESS and GREATER are both false for it.
if(NOT (_max GREATER_EQUAL MIN AND _max LESS_EQUAL MAX))
  message(FATAL_ERROR "reading the planes of '${READ_VARS}': the largest ${FIELD} is ${_max}, "
                      "outside [${MIN}, ${MAX}]; the reader did not take the planes' values")
endif()
message(STATUS "reading the planes of '${READ_VARS}': the largest ${FIELD} is ${_max}, "
               "inside [${MIN}, ${MAX}]")
