# Run a deck twice, once on a single box on one rank and once on a BoxArray split into
# several boxes, and compare the final plotfiles with fcompare. The results must not depend
# on how the domain is decomposed. The reference runs on one rank because on more ranks
# amr.refine_grid_layout splits even a single box.
# -DX= defines X as empty, so test for a value, not for DEFINED
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
separate_arguments(mpiexec_preflags  UNIX_COMMAND "${MPIEXEC_PREFLAGS}")

set(REF_DIR   "${WORKING_DIRECTORY}/one_box")
set(SPLIT_DIR "${WORKING_DIRECTORY}/split")
file(REMOVE_RECURSE "${REF_DIR}" "${SPLIT_DIR}")
file(MAKE_DIRECTORY "${REF_DIR}" "${SPLIT_DIR}")

if(NOT "${MPIEXEC}" STREQUAL "")
    set(launch_one   ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1         ${mpiexec_preflags})
    set(launch_split ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NRANKS} ${mpiexec_preflags})
else()
    set(launch_one)
    set(launch_split)
endif()

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
