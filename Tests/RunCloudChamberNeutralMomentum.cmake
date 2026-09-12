if(NOT DEFINED MPIEXEC OR NOT DEFINED MPIEXEC_NUMPROC_FLAG OR
   NOT DEFINED NRANKS OR NOT DEFINED TEST_EXE OR NOT DEFINED INPUT OR
   NOT DEFINED WORKING_DIRECTORY OR NOT DEFINED CHECKER)
    message(FATAL_ERROR "RunCloudChamberNeutralMomentum.cmake missing required argument")
endif()

set(BASELINE_DIR "${WORKING_DIRECTORY}/z0_baseline")
set(CHANGED_DIR "${WORKING_DIRECTORY}/z0_changed")
file(REMOVE_RECURSE "${BASELINE_DIR}" "${CHANGED_DIR}")
file(MAKE_DIRECTORY "${BASELINE_DIR}" "${CHANGED_DIR}")

execute_process(
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NRANKS} ${MPIEXEC_PREFLAGS}
            ${TEST_EXE} ${INPUT} xlo.z0_m=0.001
    WORKING_DIRECTORY "${BASELINE_DIR}"
    OUTPUT_FILE "${BASELINE_DIR}/simulation.log"
    ERROR_FILE "${BASELINE_DIR}/simulation.log"
    RESULT_VARIABLE baseline_result)
if(NOT baseline_result EQUAL 0)
    message(FATAL_ERROR "neutral baseline simulation failed: ${baseline_result}")
endif()

execute_process(
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NRANKS} ${MPIEXEC_PREFLAGS}
            ${TEST_EXE} ${INPUT} xlo.z0_m=0.05
    WORKING_DIRECTORY "${CHANGED_DIR}"
    OUTPUT_FILE "${CHANGED_DIR}/simulation.log"
    ERROR_FILE "${CHANGED_DIR}/simulation.log"
    RESULT_VARIABLE changed_result)
if(NOT changed_result EQUAL 0)
    message(FATAL_ERROR "neutral changed-roughness simulation failed: ${changed_result}")
endif()

execute_process(
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
            ${CHECKER} neutral_momentum
            ${BASELINE_DIR}/plt00000 ${BASELINE_DIR}/plt00060
            ${CHANGED_DIR}/plt00060
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/neutral_momentum_checker.log"
    ERROR_FILE "${WORKING_DIRECTORY}/neutral_momentum_checker.log"
    RESULT_VARIABLE checker_result)
if(NOT checker_result EQUAL 0)
    message(FATAL_ERROR "neutral momentum activation check failed: ${checker_result}")
endif()
