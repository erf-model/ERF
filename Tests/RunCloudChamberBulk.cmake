include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")

foreach(required MPIEXEC MPIEXEC_NUMPROC_FLAG NRANKS TEST_EXE INPUT
                 WORKING_DIRECTORY SIMULATION_LOG CHECKER_LOG CHECKER
                 BULK_MODEL BULK_PLOT_VARS BULK_BUDGET)
    if(NOT DEFINED ${required})
        message(FATAL_ERROR "RunCloudChamberBulk.cmake missing ${required}")
    endif()
endforeach()

if(NOT DEFINED CHECKER_MODE)
    set(CHECKER_MODE "bulk")
endif()

erf_mpi_launcher_command(mpi_launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunCloudChamberBulk.cmake")
erf_mpi_launcher_command(mpi_launch_one
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunCloudChamberBulk.cmake checker")

file(READ "${INPUT}" input_contents)
string(REPLACE "@BULK_MODEL@" "${BULK_MODEL}" input_contents "${input_contents}")
string(REPLACE "@BULK_PLOT_VARS@" "${BULK_PLOT_VARS}" input_contents "${input_contents}")
string(REPLACE "@BULK_BUDGET@" "${BULK_BUDGET}" input_contents "${input_contents}")
set(generated_input "${WORKING_DIRECTORY}/CloudChamber_Bulk_Run.i")
file(WRITE "${generated_input}" "${input_contents}")

execute_process(
    COMMAND ${mpi_launch} ${TEST_EXE} ${generated_input}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${SIMULATION_LOG}"
    ERROR_FILE "${SIMULATION_LOG}"
    RESULT_VARIABLE simulation_result)

if(DEFINED EXPECT_FAILURE AND EXPECT_FAILURE)
    if(simulation_result EQUAL 0)
        message(FATAL_ERROR "Cloud Chamber ${BULK_MODEL} was expected to fail but succeeded")
    endif()
    if(NOT DEFINED EXPECT_FAILURE_TEXT)
        message(FATAL_ERROR "Expected failure text is required for a negative regression")
    endif()
    file(READ "${SIMULATION_LOG}" simulation_output)
    string(FIND "${simulation_output}" "${EXPECT_FAILURE_TEXT}" failure_text_position)
    if(failure_text_position EQUAL -1)
        message(FATAL_ERROR
            "Cloud Chamber failed without the expected diagnostic '${EXPECT_FAILURE_TEXT}'")
    endif()
    return()
endif()

if(NOT simulation_result EQUAL 0)
    message(FATAL_ERROR "Cloud Chamber ${BULK_MODEL} simulation failed: ${simulation_result}")
endif()

execute_process(
    COMMAND ${mpi_launch_one} ${CHECKER} ${CHECKER_MODE} ${BULK_MODEL}
            "${WORKING_DIRECTORY}/plt00000" "${WORKING_DIRECTORY}/plt00002"
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${CHECKER_LOG}"
    ERROR_FILE "${CHECKER_LOG}"
    RESULT_VARIABLE checker_result)
if(NOT checker_result EQUAL 0)
    message(FATAL_ERROR "Cloud Chamber ${BULK_MODEL} property check failed: ${checker_result}")
endif()
