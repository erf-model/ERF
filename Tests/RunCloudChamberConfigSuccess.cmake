if(NOT DEFINED MPIEXEC OR NOT DEFINED MPIEXEC_NUMPROC_FLAG OR
   NOT DEFINED TEST_EXE OR NOT DEFINED INPUT OR
   NOT DEFINED WORKING_DIRECTORY OR NOT DEFINED LOG OR
   NOT DEFINED OUTPUT_DIRECTORY OR NOT DEFINED OUTPUT_ARTIFACT)
    message(FATAL_ERROR "RunCloudChamberConfigSuccess.cmake missing required argument")
endif()

if(NOT EXISTS "${MPIEXEC}")
    message(FATAL_ERROR "Cloud Chamber startup test MPI launcher is missing: ${MPIEXEC}")
endif()
if(NOT EXISTS "${TEST_EXE}")
    message(FATAL_ERROR "Cloud Chamber startup test executable is missing: ${TEST_EXE}")
endif()
if(NOT EXISTS "${INPUT}")
    message(FATAL_ERROR "Cloud Chamber startup test input is missing: ${INPUT}")
endif()

file(REMOVE_RECURSE "${OUTPUT_DIRECTORY}")
file(REMOVE "${LOG}")

execute_process(
    COMMAND "${MPIEXEC}" "${MPIEXEC_NUMPROC_FLAG}" 1 ${MPIEXEC_PREFLAGS}
            "${TEST_EXE}" "${INPUT}"
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_VARIABLE simulation_stdout
    ERROR_VARIABLE simulation_stderr
    RESULT_VARIABLE simulation_result)

set(combined_output "${simulation_stdout}\n${simulation_stderr}")
file(WRITE "${LOG}"
    "RESULT=${simulation_result}\nINPUT=${INPUT}\n\n${combined_output}")

if(NOT "${simulation_result}" STREQUAL "0")
    message(FATAL_ERROR
        "Cloud Chamber legacy startup failed with result ${simulation_result}; see ${LOG}")
endif()

if(NOT EXISTS "${OUTPUT_ARTIFACT}")
    message(FATAL_ERROR
        "Cloud Chamber legacy startup exited successfully without expected output artifact ${OUTPUT_ARTIFACT}; see ${LOG}")
endif()
