include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")

if(NOT DEFINED MPIEXEC OR NOT DEFINED MPIEXEC_NUMPROC_FLAG OR
   NOT DEFINED NRANKS OR NOT DEFINED TEST_EXE OR NOT DEFINED INPUT OR
   NOT DEFINED WORKING_DIRECTORY OR NOT DEFINED CHECKER OR NOT DEFINED MODE)
    message(FATAL_ERROR "RunCloudChamberBudget.cmake missing required argument")
endif()

# MPIEXEC may be a multi-word command such as "flux run", so the launcher
# prefix is built once by the shared helper instead of being pasted into
# each COMMAND as a single token.
erf_mpi_launcher_command(mpi_launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunCloudChamberBudget.cmake")
erf_mpi_launcher_command(mpi_launch_one
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunCloudChamberBudget.cmake")

file(REMOVE_RECURSE "${WORKING_DIRECTORY}/plt00000"
                    "${WORKING_DIRECTORY}/plt00002"
                    "${WORKING_DIRECTORY}/plt00004"
                    "${WORKING_DIRECTORY}/cloud_chamber_budget.dat")

execute_process(
    COMMAND ${mpi_launch}
            ${TEST_EXE} ${INPUT}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/simulation.log"
    ERROR_FILE "${WORKING_DIRECTORY}/simulation.log"
    RESULT_VARIABLE simulation_result)
if(NOT simulation_result EQUAL 0)
    message(FATAL_ERROR "Cloud Chamber budget simulation failed: ${simulation_result}")
endif()

execute_process(
    COMMAND ${mpi_launch_one}
            ${CHECKER} ${MODE} ${WORKING_DIRECTORY}/plt00000
            ${WORKING_DIRECTORY}/plt00004 ${WORKING_DIRECTORY}/cloud_chamber_budget.dat
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/checker.log"
    ERROR_FILE "${WORKING_DIRECTORY}/checker.log"
    RESULT_VARIABLE checker_result)
if(NOT checker_result EQUAL 0)
    message(FATAL_ERROR "Cloud Chamber budget check failed: ${checker_result}")
endif()
