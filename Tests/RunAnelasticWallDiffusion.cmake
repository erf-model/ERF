if(NOT DEFINED MPIEXEC OR NOT DEFINED MPIEXEC_NUMPROC_FLAG OR
   NOT DEFINED NRANKS OR NOT DEFINED TEST_EXE OR NOT DEFINED INPUT OR
   NOT DEFINED WORKING_DIRECTORY OR NOT DEFINED SIMULATION_LOG OR
   NOT DEFINED CHECKER_LOG OR NOT DEFINED CHECKER OR
   NOT DEFINED PLOTFILE OR NOT DEFINED AXIS OR NOT DEFINED THETA_LO OR
   NOT DEFINED THETA_HI)
    message(FATAL_ERROR "RunAnelasticWallDiffusion.cmake missing required argument")
endif()

execute_process(
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NRANKS} ${MPIEXEC_PREFLAGS}
            ${TEST_EXE} ${INPUT}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${SIMULATION_LOG}"
    ERROR_FILE "${SIMULATION_LOG}"
    RESULT_VARIABLE simulation_result)
if(NOT simulation_result EQUAL 0)
    message(FATAL_ERROR "anelastic wall diffusion simulation failed: ${simulation_result}")
endif()

execute_process(
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
            ${CHECKER} ${PLOTFILE} ${AXIS} ${THETA_LO} ${THETA_HI}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${CHECKER_LOG}"
    ERROR_FILE "${CHECKER_LOG}"
    RESULT_VARIABLE checker_result)
if(NOT checker_result EQUAL 0)
    message(FATAL_ERROR "anelastic wall diffusion property check failed: ${checker_result}")
endif()
