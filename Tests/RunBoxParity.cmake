# Run a deck twice, once on a single box on one rank and once on a BoxArray split into
# several boxes, and compare the final plotfiles with fcompare. The results must not depend
# on how the domain is decomposed. The reference runs on one rank because on more ranks
# amr.refine_grid_layout splits even a single box.
# -DX= defines X as empty, so test for a value, not for DEFINED
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")

foreach(arg NRANKS TEST_EXE INPUT WORKING_DIRECTORY FCOMPARE PLTFILE RTOL ATOL)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunBoxParity.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunBoxParity.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

separate_arguments(common_options    UNIX_COMMAND "${COMMON_OPTIONS}")
separate_arguments(reference_options UNIX_COMMAND "${REFERENCE_OPTIONS}")
separate_arguments(split_options     UNIX_COMMAND "${SPLIT_OPTIONS}")

set(REF_DIR   "${WORKING_DIRECTORY}/one_box")
set(SPLIT_DIR "${WORKING_DIRECTORY}/split")
file(REMOVE_RECURSE "${REF_DIR}" "${SPLIT_DIR}")
file(MAKE_DIRECTORY "${REF_DIR}" "${SPLIT_DIR}")

# MPIEXEC may be a multi-word command such as "flux run"; the helper splits
# it, validates the program and applies MPIEXEC_PREFLAGS. An empty MPIEXEC
# yields an empty prefix, so the runs stay serial.
erf_mpi_launcher_command(launch_one
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunBoxParity.cmake")
erf_mpi_launcher_command(launch_split
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunBoxParity.cmake")

execute_process(
    COMMAND ${launch_one} ${TEST_EXE} ${INPUT} ${common_options} ${reference_options}
    WORKING_DIRECTORY "${REF_DIR}"
    OUTPUT_FILE "${REF_DIR}/simulation.log"
    ERROR_FILE "${REF_DIR}/simulation.log"
    RESULT_VARIABLE ref_result)
if(NOT ref_result EQUAL 0)
    message(FATAL_ERROR "Single-box simulation failed: ${ref_result}")
endif()

execute_process(
    COMMAND ${launch_split} ${TEST_EXE} ${INPUT} ${common_options} ${split_options}
    WORKING_DIRECTORY "${SPLIT_DIR}"
    OUTPUT_FILE "${SPLIT_DIR}/simulation.log"
    ERROR_FILE "${SPLIT_DIR}/simulation.log"
    RESULT_VARIABLE split_result)
if(NOT split_result EQUAL 0)
    message(FATAL_ERROR "Split-box simulation failed: ${split_result}")
endif()

# The comparison proves nothing unless the split run is in fact split: the number of boxes
# on level 0 is the first number on the fifth line of Level_0/Cell_H
function(count_level0_boxes plotfile out_var)
    set(_cell_h "${plotfile}/Level_0/Cell_H")
    if(NOT EXISTS "${_cell_h}")
        message(FATAL_ERROR "RunBoxParity.cmake: no ${_cell_h}")
    endif()
    file(STRINGS "${_cell_h}" _lines LIMIT_COUNT 5)
    list(GET _lines 4 _boxes_line)
    if(NOT _boxes_line MATCHES "^\\(([0-9]+) ")
        message(FATAL_ERROR "RunBoxParity.cmake: cannot read the box count from ${_cell_h}: '${_boxes_line}'")
    endif()
    set(${out_var} ${CMAKE_MATCH_1} PARENT_SCOPE)
endfunction()
count_level0_boxes("${REF_DIR}/${PLTFILE}"   ref_boxes)
count_level0_boxes("${SPLIT_DIR}/${PLTFILE}" split_boxes)
if(NOT split_boxes GREATER ref_boxes)
    message(FATAL_ERROR "RunBoxParity.cmake: the split run has ${split_boxes} boxes on level 0 and the reference ${ref_boxes}; the split run must have more, or the comparison is trivial")
endif()
message(STATUS "RunBoxParity: reference ${ref_boxes} boxes, split ${split_boxes} boxes on level 0")

execute_process(
    COMMAND ${launch_one} ${FCOMPARE} --abort_if_not_all_found --allow_diff_grids
            --rel_tol ${RTOL} --abs_tol ${ATOL}
            ${REF_DIR}/${PLTFILE} ${SPLIT_DIR}/${PLTFILE}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/parity.log"
    ERROR_FILE "${WORKING_DIRECTORY}/parity.log"
    RESULT_VARIABLE parity_result)
if(NOT parity_result EQUAL 0)
    message(FATAL_ERROR "Single-box and split-box plotfiles differ: ${parity_result}")
endif()

# Optional: the data logs (erf.data_log) of the two runs must agree.  fcompare never reads them,
# and planar diagnostics can depend on the decomposition while the plotfiles do not.  The logs are
# compared numerically, not byte for byte: they print six significant digits, and a planar sum
# reduces in a different order for each decomposition, so a value near a rounding boundary can
# print a different last digit without anything being wrong.
if(NOT "${DATALOG}" STREQUAL "")
    foreach(dir "${REF_DIR}" "${SPLIT_DIR}")
        if(NOT EXISTS "${dir}/${DATALOG}")
            message(FATAL_ERROR "RunBoxParity.cmake: no data log ${dir}/${DATALOG}")
        endif()
    endforeach()
    file(STRINGS "${REF_DIR}/${DATALOG}" ref_log)
    list(LENGTH ref_log ref_lines)
    if(ref_lines LESS 2)
        message(FATAL_ERROR "RunBoxParity.cmake: data log ${DATALOG} has ${ref_lines} lines; the comparison would be trivial")
    endif()
    include("${CMAKE_CURRENT_LIST_DIR}/CompareDataLogs.cmake")
    # datprecision in Source/ERF.H, and a couple of units of the last digit of tolerance
    erf_compare_data_logs("${REF_DIR}/${DATALOG}" "${SPLIT_DIR}/${DATALOG}" 6 2 logs_agree log_message)
    if(NOT logs_agree)
        message(FATAL_ERROR "RunBoxParity.cmake: data log ${DATALOG} differs between the single-box "
                            "and split runs: ${log_message}")
    endif()
    message(STATUS "RunBoxParity: data log ${DATALOG} agrees (${ref_lines} lines)")
endif()
