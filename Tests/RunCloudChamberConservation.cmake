if(NOT DEFINED MPIEXEC OR NOT DEFINED MPIEXEC_NUMPROC_FLAG OR
   NOT DEFINED TEST_EXE OR NOT DEFINED INPUT OR NOT DEFINED WORKING_DIRECTORY OR
   NOT DEFINED CHECKER)
    message(FATAL_ERROR "RunCloudChamberConservation.cmake missing required argument")
endif()
if(NOT DEFINED MPIEXEC_PREFLAGS)
    set(MPIEXEC_PREFLAGS "")
endif()

file(REMOVE_RECURSE "${WORKING_DIRECTORY}/cloud_chamber_water_ledger.dat")

execute_process(
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
            ${TEST_EXE} ${INPUT}
            amr.n_cell=16 16 16 amr.max_grid_size=8
            fabarray.mfiter_tile_size=8 8 8
            erf.fixed_dt=0.005 stop_time=0.005 max_step=1
            erf.plot_int_1=-1 erf.cloud_chamber_budget_interval=-1
            erf.cloud_chamber_water_ledger=true erf.v=0 amr.v=0
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/simulation.log"
    ERROR_FILE "${WORKING_DIRECTORY}/simulation.log"
    RESULT_VARIABLE simulation_result)
if(NOT simulation_result EQUAL 0)
    message(FATAL_ERROR "Cloud Chamber untiled-containment simulation failed: ${simulation_result}")
endif()

execute_process(
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
            ${CHECKER} ledger containment ${WORKING_DIRECTORY}/cloud_chamber_water_ledger.dat
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/checker.log"
    ERROR_FILE "${WORKING_DIRECTORY}/checker.log"
    RESULT_VARIABLE checker_result)
if(NOT checker_result EQUAL 0)
    message(FATAL_ERROR "Cloud Chamber untiled-containment check failed: ${checker_result}")
endif()
