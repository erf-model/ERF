if(NOT DEFINED MPIEXEC OR NOT DEFINED MPIEXEC_NUMPROC_FLAG OR
   NOT DEFINED NRANKS OR NOT DEFINED TEST_EXE OR NOT DEFINED INPUT OR
   NOT DEFINED WORKING_DIRECTORY OR NOT DEFINED CHECKER)
    message(FATAL_ERROR "RunCloudChamberParity.cmake missing required argument")
endif()

set(OFF_DIR "${WORKING_DIRECTORY}/budget_off")
set(ON_DIR "${WORKING_DIRECTORY}/budget_on")
file(REMOVE_RECURSE "${OFF_DIR}" "${ON_DIR}")
file(MAKE_DIRECTORY "${OFF_DIR}" "${ON_DIR}")

execute_process(
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NRANKS} ${MPIEXEC_PREFLAGS}
            ${TEST_EXE} ${INPUT}
            amr.n_cell=16 16 16 amr.max_grid_size=8 max_step=4
            erf.plot_int_1=4 erf.cloud_chamber_budget_interval=-1
    WORKING_DIRECTORY "${OFF_DIR}"
    OUTPUT_FILE "${OFF_DIR}/simulation.log"
    ERROR_FILE "${OFF_DIR}/simulation.log"
    RESULT_VARIABLE off_result)
if(NOT off_result EQUAL 0)
    message(FATAL_ERROR "Cloud Chamber budget-off simulation failed: ${off_result}")
endif()

execute_process(
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NRANKS} ${MPIEXEC_PREFLAGS}
            ${TEST_EXE} ${INPUT}
            amr.n_cell=16 16 16 amr.max_grid_size=8 max_step=4
            erf.plot_int_1=4 erf.cloud_chamber_budget_interval=1
    WORKING_DIRECTORY "${ON_DIR}"
    OUTPUT_FILE "${ON_DIR}/simulation.log"
    ERROR_FILE "${ON_DIR}/simulation.log"
    RESULT_VARIABLE on_result)
if(NOT on_result EQUAL 0)
    message(FATAL_ERROR "Cloud Chamber budget-on simulation failed: ${on_result}")
endif()

execute_process(
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
            ${CHECKER} parity ${OFF_DIR}/plt00004 ${ON_DIR}/plt00004
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/parity.log"
    ERROR_FILE "${WORKING_DIRECTORY}/parity.log"
    RESULT_VARIABLE parity_result)
if(NOT parity_result EQUAL 0)
    message(FATAL_ERROR "Cloud Chamber budget parity failed: ${parity_result}")
endif()
