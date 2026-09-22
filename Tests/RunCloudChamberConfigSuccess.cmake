cmake_minimum_required(VERSION 3.24)

include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")

if(NOT DEFINED MPIEXEC OR NOT DEFINED MPIEXEC_NUMPROC_FLAG OR
   NOT DEFINED TEST_EXE OR NOT DEFINED INPUT OR
   NOT DEFINED WORKING_DIRECTORY OR NOT DEFINED LOG OR
   NOT DEFINED OUTPUT_DIRECTORY OR NOT DEFINED OUTPUT_ARTIFACT)
    message(FATAL_ERROR "RunCloudChamberConfigSuccess.cmake missing required argument")
endif()

if(NOT EXISTS "${TEST_EXE}")
    message(FATAL_ERROR "Cloud Chamber startup test executable is missing: ${TEST_EXE}")
endif()
if(NOT EXISTS "${INPUT}")
    message(FATAL_ERROR "Cloud Chamber startup test input is missing: ${INPUT}")
endif()

# MPIEXEC may be a multi-word command such as "flux run", so it is split and
# validated by the shared helper rather than treated as one program path.
# The helper also drops the rank flag when MPIEXEC_NUMPROC_FLAG is empty and
# splits MPIEXEC_PREFLAGS into separate arguments.
erf_mpi_launcher_command(run_command
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "Cloud Chamber startup test")

list(APPEND run_command "${TEST_EXE}" "${INPUT}")

file(REMOVE_RECURSE "${OUTPUT_DIRECTORY}")
file(REMOVE "${LOG}")

execute_process(
    COMMAND ${run_command}
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
