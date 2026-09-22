include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")

if(NOT DEFINED MPIEXEC OR NOT DEFINED MPIEXEC_NUMPROC_FLAG OR
   NOT DEFINED NRANKS OR NOT DEFINED TEST_EXE OR NOT DEFINED INPUT OR
   NOT DEFINED WORKING_DIRECTORY OR NOT DEFINED CHECKER)
    message(FATAL_ERROR "RunCloudChamberParity.cmake missing required argument")
endif()

# MPIEXEC may be a multi-word command such as "flux run", so the launcher
# prefix is built once by the shared helper instead of being pasted into
# each COMMAND as a single token.
erf_mpi_launcher_command(mpi_launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunCloudChamberParity.cmake")
erf_mpi_launcher_command(mpi_launch_one
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunCloudChamberParity.cmake")

set(OFF_DIR "${WORKING_DIRECTORY}/budget_off")
set(ON_DIR "${WORKING_DIRECTORY}/budget_on")
file(REMOVE_RECURSE "${OFF_DIR}" "${ON_DIR}")
file(MAKE_DIRECTORY "${OFF_DIR}" "${ON_DIR}")

execute_process(
    COMMAND ${mpi_launch}
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
    COMMAND ${mpi_launch}
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
    COMMAND ${mpi_launch_one}
            ${CHECKER} parity ${OFF_DIR}/plt00004 ${ON_DIR}/plt00004
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/parity.log"
    ERROR_FILE "${WORKING_DIRECTORY}/parity.log"
    RESULT_VARIABLE parity_result)
if(NOT parity_result EQUAL 0)
    message(FATAL_ERROR "Cloud Chamber budget parity failed: ${parity_result}")
endif()
