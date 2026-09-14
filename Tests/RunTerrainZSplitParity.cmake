# Run a two-level terrain deck with whole-height fine grids and with the fine grids split
# in z, and require the two plotfiles to agree bit for bit.
if(NOT DEFINED MPIEXEC OR NOT DEFINED MPIEXEC_NUMPROC_FLAG OR
   NOT DEFINED NRANKS OR NOT DEFINED TEST_EXE OR NOT DEFINED INPUT OR
   NOT DEFINED WORKING_DIRECTORY OR NOT DEFINED FCOMPARE OR NOT DEFINED PLTFILE)
    message(FATAL_ERROR "RunTerrainZSplitParity.cmake missing required argument")
endif()

set(FULL_DIR "${WORKING_DIRECTORY}/full_columns")
set(SPLIT_DIR "${WORKING_DIRECTORY}/split_in_z")
file(REMOVE_RECURSE "${FULL_DIR}" "${SPLIT_DIR}")
file(MAKE_DIRECTORY "${FULL_DIR}" "${SPLIT_DIR}")

foreach(RUN IN ITEMS "full_columns|1024 1024" "split_in_z|1024 16")
    string(REPLACE "|" ";" RUN "${RUN}")
    list(GET RUN 0 RUN_NAME)
    list(GET RUN 1 MAX_GRID_SIZE_Z)
    execute_process(
        COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NRANKS} ${MPIEXEC_PREFLAGS}
                ${TEST_EXE} ${INPUT} "amr.max_grid_size_z=${MAX_GRID_SIZE_Z}"
        WORKING_DIRECTORY "${WORKING_DIRECTORY}/${RUN_NAME}"
        OUTPUT_FILE "${WORKING_DIRECTORY}/${RUN_NAME}/simulation.log"
        ERROR_FILE "${WORKING_DIRECTORY}/${RUN_NAME}/simulation.log"
        RESULT_VARIABLE run_result)
    if(NOT run_result EQUAL 0)
        message(FATAL_ERROR "${RUN_NAME} simulation failed: ${run_result}")
    endif()
endforeach()

# The comparison only means something if the split run really has more level-1 boxes: an
# override that stopped taking effect would give two identical plotfiles, which agree.
foreach(RUN_NAME IN ITEMS full_columns split_in_z)
    set(CELL_H "${WORKING_DIRECTORY}/${RUN_NAME}/${PLTFILE}/Level_1/Cell_H")
    if(NOT EXISTS "${CELL_H}")
        message(FATAL_ERROR "${RUN_NAME} wrote no level-1 data (${CELL_H})")
    endif()
    file(READ "${CELL_H}" CELL_H_TEXT)
    # The BoxArray is written as "(<number of boxes> 0" followed by one line per box
    if(NOT CELL_H_TEXT MATCHES "\n\\(([0-9]+) 0\n")
        message(FATAL_ERROR "Cannot read the level-1 BoxArray from ${CELL_H}")
    endif()
    set(NBOXES_${RUN_NAME} "${CMAKE_MATCH_1}")
endforeach()
if(NOT NBOXES_split_in_z GREATER NBOXES_full_columns)
    message(FATAL_ERROR "The split run has ${NBOXES_split_in_z} level-1 boxes and the "
                        "whole-column run ${NBOXES_full_columns}: the fine grids were not split in z")
endif()

# The BoxArrays differ, so allow different grids; any difference at all fails
execute_process(
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
            ${FCOMPARE} --abort_if_not_all_found -a -r 0 --abs_tol 0
            ${FULL_DIR}/${PLTFILE} ${SPLIT_DIR}/${PLTFILE}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/parity.log"
    ERROR_FILE "${WORKING_DIRECTORY}/parity.log"
    RESULT_VARIABLE parity_result)
if(NOT parity_result EQUAL 0)
    message(FATAL_ERROR "Fine-level plotfiles depend on the z split: ${parity_result}")
endif()
