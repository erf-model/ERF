include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")

if(NOT DEFINED MPIEXEC OR NOT DEFINED MPIEXEC_NUMPROC_FLAG OR
   NOT DEFINED NRANKS OR NOT DEFINED TEST_EXE OR NOT DEFINED INPUT OR
   NOT DEFINED WORKING_DIRECTORY OR NOT DEFINED SIMULATION_LOG OR
   NOT DEFINED CHECKER_LOG OR NOT DEFINED CHECKER OR NOT DEFINED MODE)
    message(FATAL_ERROR "RunCloudChamber.cmake missing required argument")
endif()

# MPIEXEC may be a multi-word command such as "flux run", so the launcher
# prefix is built once by the shared helper instead of being pasted into
# each COMMAND as a single token.
erf_mpi_launcher_command(mpi_launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunCloudChamber.cmake")
erf_mpi_launcher_command(mpi_launch_one
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunCloudChamber.cmake")

execute_process(
    COMMAND ${mpi_launch}
            ${TEST_EXE} ${INPUT}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${SIMULATION_LOG}"
    ERROR_FILE "${SIMULATION_LOG}"
    RESULT_VARIABLE simulation_result)
if(NOT simulation_result EQUAL 0)
    message(FATAL_ERROR "Cloud Chamber simulation failed: ${simulation_result}")
endif()

execute_process(
    COMMAND ${mpi_launch_one}
            ${CHECKER} ${MODE} ${WORKING_DIRECTORY}/plt00000 ${WORKING_DIRECTORY}/plt00002
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${CHECKER_LOG}"
    ERROR_FILE "${CHECKER_LOG}"
    RESULT_VARIABLE checker_result)
if(NOT checker_result EQUAL 0)
    message(FATAL_ERROR "Cloud Chamber property check failed: ${checker_result}")
endif()
