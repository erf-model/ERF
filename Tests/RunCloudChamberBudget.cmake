if(NOT DEFINED MPIEXEC OR NOT DEFINED MPIEXEC_NUMPROC_FLAG OR
   NOT DEFINED NRANKS OR NOT DEFINED TEST_EXE OR NOT DEFINED INPUT OR
   NOT DEFINED WORKING_DIRECTORY OR NOT DEFINED CHECKER OR NOT DEFINED MODE)
    message(FATAL_ERROR "RunCloudChamberBudget.cmake missing required argument")
endif()

file(REMOVE_RECURSE "${WORKING_DIRECTORY}/plt00000"
                    "${WORKING_DIRECTORY}/plt00002"
                    "${WORKING_DIRECTORY}/plt00004"
                    "${WORKING_DIRECTORY}/plt00006"
                    "${WORKING_DIRECTORY}/cloud_chamber_budget.dat")

execute_process(
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NRANKS} ${MPIEXEC_PREFLAGS}
            ${TEST_EXE} ${INPUT}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/simulation.log"
    ERROR_FILE "${WORKING_DIRECTORY}/simulation.log"
    RESULT_VARIABLE simulation_result)
if(NOT simulation_result EQUAL 0)
    message(FATAL_ERROR "Cloud Chamber budget simulation failed: ${simulation_result}")
endif()

execute_process(
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
            ${CHECKER} ${MODE} ${WORKING_DIRECTORY}/plt00000
            ${WORKING_DIRECTORY}/plt00006 ${WORKING_DIRECTORY}/cloud_chamber_budget.dat
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/checker.log"
    ERROR_FILE "${WORKING_DIRECTORY}/checker.log"
    RESULT_VARIABLE checker_result)
if(NOT checker_result EQUAL 0)
    message(FATAL_ERROR "Cloud Chamber budget check failed: ${checker_result}")
endif()
