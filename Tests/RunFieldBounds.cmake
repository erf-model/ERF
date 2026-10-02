# Run one ERF deck and require the extrema of one plotfile variable to lie within [LO, HI].
# Arguments (-D): MPIEXEC, MPIEXEC_NUMPROC_FLAG, MPIEXEC_PREFLAGS, NRANKS, TEST_EXE, INPUT,
# WORKING_DIRECTORY, FEXTREMA, PLTFILE, VARIABLE, LO, HI.
cmake_minimum_required(VERSION 3.20)
include(${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake)

foreach(arg TEST_EXE INPUT WORKING_DIRECTORY FEXTREMA PLTFILE VARIABLE LO HI)
    if(NOT DEFINED ${arg})
        message(FATAL_ERROR "RunFieldBounds.cmake: missing required argument ${arg}")
    endif()
endforeach()

erf_mpi_launcher_command(_launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunFieldBounds.cmake")
erf_mpi_launcher_command(_launch_one
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunFieldBounds.cmake")

# a plotfile left by an earlier run must not satisfy the check
file(REMOVE_RECURSE "${WORKING_DIRECTORY}/${PLTFILE}")

execute_process(
    COMMAND ${_launch} ${TEST_EXE} ${INPUT}
    WORKING_DIRECTORY ${WORKING_DIRECTORY}
    OUTPUT_FILE ${WORKING_DIRECTORY}/simulation.log
    ERROR_FILE  ${WORKING_DIRECTORY}/simulation.log
    RESULT_VARIABLE run_result)
if(NOT run_result EQUAL 0)
    message(FATAL_ERROR "RunFieldBounds.cmake: the run failed (${run_result}); see simulation.log")
endif()

set(PLT "${WORKING_DIRECTORY}/${PLTFILE}")
if(NOT EXISTS "${PLT}")
    message(FATAL_ERROR "RunFieldBounds.cmake: the run wrote no plotfile ${PLT}")
endif()
execute_process(
    COMMAND ${_launch_one} ${FEXTREMA} -v "${VARIABLE}" "${PLT}"
    OUTPUT_VARIABLE extrema_out
    ERROR_VARIABLE extrema_err
    RESULT_VARIABLE extrema_result)
if(NOT extrema_result EQUAL 0)
    message(FATAL_ERROR "RunFieldBounds.cmake: fextrema failed on ${PLT}: ${extrema_result}\n${extrema_err}")
endif()
if(NOT extrema_out MATCHES "${VARIABLE}[ \t]+([-+0-9.eE]+)[ \t]+([-+0-9.eE]+)")
    message(FATAL_ERROR "RunFieldBounds.cmake: cannot read the ${VARIABLE} extrema from:\n${extrema_out}")
endif()
set(V_MIN "${CMAKE_MATCH_1}")
set(V_MAX "${CMAKE_MATCH_2}")
message(STATUS "RunFieldBounds.cmake: ${VARIABLE} in [${V_MIN}, ${V_MAX}], required [${LO}, ${HI}]")
if(V_MIN LESS LO OR V_MAX GREATER HI)
    message(FATAL_ERROR "RunFieldBounds.cmake: ${VARIABLE} leaves [${LO}, ${HI}]: min ${V_MIN}, max ${V_MAX}")
endif()
