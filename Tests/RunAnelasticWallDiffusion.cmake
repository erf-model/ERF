if(NOT DEFINED MPIEXEC OR NOT DEFINED MPIEXEC_NUMPROC_FLAG OR
   NOT DEFINED NRANKS OR NOT DEFINED TEST_EXE OR NOT DEFINED INPUT OR
   NOT DEFINED WORKING_DIRECTORY OR NOT DEFINED SIMULATION_LOG OR
   NOT DEFINED CHECKER_LOG OR NOT DEFINED CHECKER OR
   NOT DEFINED PLOTFILE OR NOT DEFINED AXIS OR NOT DEFINED THETA_LO OR
   NOT DEFINED THETA_HI)
    message(FATAL_ERROR "RunAnelasticWallDiffusion.cmake missing required argument")
endif()

set(_have_slurm_allocation FALSE)
if(DEFINED ENV{SLURM_JOB_ID} OR DEFINED ENV{SLURM_JOBID})
    set(_have_slurm_allocation TRUE)
endif()

set(_simulation_command)
if(_have_slurm_allocation)
    list(APPEND _simulation_command
        ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NRANKS})
    if(DEFINED MPIEXEC_PREFLAGS AND NOT "${MPIEXEC_PREFLAGS}" STREQUAL "")
        separate_arguments(_mpi_preflags UNIX_COMMAND "${MPIEXEC_PREFLAGS}")
        list(APPEND _simulation_command ${_mpi_preflags})
    endif()
else()
    message(STATUS
        "No Slurm allocation detected; running anelastic wall diffusion directly with one rank")
endif()
list(APPEND _simulation_command ${TEST_EXE} ${INPUT})

execute_process(
    COMMAND ${_simulation_command}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${SIMULATION_LOG}"
    ERROR_FILE "${SIMULATION_LOG}"
    RESULT_VARIABLE simulation_result)
if(NOT simulation_result EQUAL 0)
    message(FATAL_ERROR "anelastic wall diffusion simulation failed: ${simulation_result}")
endif()

set(_checker_command)
if(_have_slurm_allocation)
    list(APPEND _checker_command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1)
    if(DEFINED MPIEXEC_PREFLAGS AND NOT "${MPIEXEC_PREFLAGS}" STREQUAL "")
        list(APPEND _checker_command ${_mpi_preflags})
    endif()
else()
    message(STATUS
        "No Slurm allocation detected; running anelastic wall diffusion checker directly")
endif()
list(APPEND _checker_command ${CHECKER} ${PLOTFILE} ${AXIS} ${THETA_LO} ${THETA_HI})

execute_process(
    COMMAND ${_checker_command}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${CHECKER_LOG}"
    ERROR_FILE "${CHECKER_LOG}"
    RESULT_VARIABLE checker_result)
if(NOT checker_result EQUAL 0)
    message(FATAL_ERROR "anelastic wall diffusion property check failed: ${checker_result}")
endif()
