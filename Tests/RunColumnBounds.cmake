# Run one ERF deck and require, in the vertical column through the center of the domain, every
# cell at or above Z_MIN to lie within [LO, HI] for each listed plotfile variable.  Cells below
# Z_MIN are skipped, so a deck can keep a wall at the bottom (the RANS wall distance needs one)
# while the check covers only the cells its influence has not reached.
# Arguments (-D): MPIEXEC, MPIEXEC_NUMPROC_FLAG, MPIEXEC_PREFLAGS, NRANKS, TEST_EXE, INPUT,
# WORKING_DIRECTORY, FEXTRACT, PLTFILE, Z_MIN, VARIABLES, LOS, HIS (the last three are
# comma-separated lists of the same length; a ';' would split the add_test argument).
# Optional: Z_MAX skips the cells above it as well, so that a check can single out one cell
# or one layer of the column; N_CELLS is the number of cells the check must cover, so that a
# change of the grid cannot widen or narrow it without notice; OPTIONS is a space-separated
# list of inputs appended to the command line, so that several tests can run one deck with
# different settings.
cmake_minimum_required(VERSION 3.20)
include(${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake)

foreach(arg TEST_EXE INPUT WORKING_DIRECTORY FEXTRACT PLTFILE Z_MIN VARIABLES LOS HIS)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunColumnBounds.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
string(REPLACE "," ";" VARIABLES "${VARIABLES}")
string(REPLACE "," ";" LOS "${LOS}")
string(REPLACE "," ";" HIS "${HIS}")
list(LENGTH VARIABLES n_vars)
list(LENGTH LOS n_los)
list(LENGTH HIS n_his)
if(NOT n_vars EQUAL n_los OR NOT n_vars EQUAL n_his)
    message(FATAL_ERROR "RunColumnBounds.cmake: VARIABLES, LOS and HIS must have the same length")
endif()

erf_mpi_launcher_command(_launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunColumnBounds.cmake")
erf_mpi_launcher_command(_launch_one
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunColumnBounds.cmake")

separate_arguments(run_options UNIX_COMMAND "${OPTIONS}")
set(range_text "at or above z = ${Z_MIN}")
if(NOT "${Z_MAX}" STREQUAL "")
    if(Z_MAX LESS Z_MIN)
        message(FATAL_ERROR "RunColumnBounds.cmake: Z_MAX = ${Z_MAX} is below Z_MIN = ${Z_MIN}")
    endif()
    set(range_text "with z in [${Z_MIN}, ${Z_MAX}]")
endif()

# a plotfile left by an earlier run must not satisfy the check
file(REMOVE_RECURSE "${WORKING_DIRECTORY}/${PLTFILE}")

execute_process(
    COMMAND ${_launch} ${TEST_EXE} ${INPUT} ${run_options}
    WORKING_DIRECTORY ${WORKING_DIRECTORY}
    OUTPUT_FILE ${WORKING_DIRECTORY}/simulation.log
    ERROR_FILE  ${WORKING_DIRECTORY}/simulation.log
    RESULT_VARIABLE run_result)
if(NOT run_result EQUAL 0)
    message(FATAL_ERROR "RunColumnBounds.cmake: the run failed (${run_result}); see simulation.log")
endif()

set(PLT "${WORKING_DIRECTORY}/${PLTFILE}")
if(NOT EXISTS "${PLT}")
    message(FATAL_ERROR "RunColumnBounds.cmake: the run wrote no plotfile ${PLT}")
endif()

string(REPLACE ";" " " var_string "${VARIABLES}")
set(column_file "${WORKING_DIRECTORY}/column.txt")
file(REMOVE "${column_file}")
execute_process(
    COMMAND ${_launch_one} ${FEXTRACT} -d 2 -p 17 -v "${var_string}" -s "${column_file}" "${PLT}"
    OUTPUT_VARIABLE extract_out
    ERROR_VARIABLE extract_err
    RESULT_VARIABLE extract_result)
if(NOT extract_result EQUAL 0 OR NOT EXISTS "${column_file}")
    message(FATAL_ERROR "RunColumnBounds.cmake: fextract failed on ${PLT}: ${extract_result}\n${extract_err}")
endif()

file(STRINGS "${column_file}" lines)
set(n_checked 0)
foreach(line IN LISTS lines)
    if(line MATCHES "^[ \t]*#" OR line MATCHES "^[ \t]*$")
        continue()
    endif()
    string(STRIP "${line}" line)
    string(REGEX REPLACE "[ \t]+" ";" fields "${line}")
    list(GET fields 0 z)
    if(z LESS Z_MIN)
        continue()
    endif()
    if(NOT "${Z_MAX}" STREQUAL "" AND z GREATER Z_MAX)
        continue()
    endif()
    math(EXPR last "${n_vars} - 1")
    foreach(iv RANGE ${last})
        math(EXPR icol "${iv} + 1")
        list(GET fields ${icol} value)
        list(GET VARIABLES ${iv} name)
        list(GET LOS ${iv} lo)
        list(GET HIS ${iv} hi)
        # written so that a value that is not a number (nan, inf, empty) fails too
        if(NOT (value GREATER_EQUAL lo AND value LESS_EQUAL hi))
            message(FATAL_ERROR "RunColumnBounds.cmake: ${name} = ${value} at z = ${z} is outside [${lo}, ${hi}]")
        endif()
    endforeach()
    math(EXPR n_checked "${n_checked} + 1")
endforeach()
if(n_checked EQUAL 0)
    message(FATAL_ERROR "RunColumnBounds.cmake: no cell ${range_text} in ${column_file}")
endif()
if(NOT "${N_CELLS}" STREQUAL "" AND NOT n_checked EQUAL N_CELLS)
    message(FATAL_ERROR "RunColumnBounds.cmake: ${n_checked} cells ${range_text} in ${column_file}, expected ${N_CELLS}")
endif()
message(STATUS "RunColumnBounds.cmake: ${n_checked} cells ${range_text} within bounds for ${VARIABLES}")
