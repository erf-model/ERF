if(NOT DEFINED MPIEXEC OR NOT DEFINED MPIEXEC_NUMPROC_FLAG OR
   NOT DEFINED NRANKS OR NOT DEFINED TEST_EXE OR NOT DEFINED INPUT OR
   NOT DEFINED WORKING_DIRECTORY OR NOT DEFINED CHECKER OR NOT DEFINED CMAKE_COMMAND)
    message(FATAL_ERROR "RunCloudChamberOpenMP.cmake missing required argument")
endif()

set(ONE_DIR "${WORKING_DIRECTORY}/omp_1_thread")
set(TWO_DIR "${WORKING_DIRECTORY}/omp_2_threads")
file(REMOVE_RECURSE "${ONE_DIR}" "${TWO_DIR}")
file(MAKE_DIRECTORY "${ONE_DIR}" "${TWO_DIR}")

set(THREAD_COUNTS 1 2)
set(THREAD_DIRS "${ONE_DIR}" "${TWO_DIR}")
foreach(THREADS DIR IN ZIP_LISTS THREAD_COUNTS THREAD_DIRS)
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E env OMP_NUM_THREADS=${THREADS}
                ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NRANKS} ${MPIEXEC_PREFLAGS}
                ${TEST_EXE} ${INPUT}
                amr.n_cell=16 16 16 amr.max_grid_size=8 max_step=6
                erf.plot_int_1=6 erf.cloud_chamber_budget_interval=1
        WORKING_DIRECTORY "${DIR}"
        OUTPUT_FILE "${DIR}/simulation.log"
        ERROR_FILE "${DIR}/simulation.log"
        RESULT_VARIABLE simulation_result)
    if(NOT simulation_result EQUAL 0)
        message(FATAL_ERROR "Cloud Chamber OpenMP ${THREADS}-thread simulation failed: ${simulation_result}")
    endif()
endforeach()

execute_process(
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
            ${CHECKER} parity ${ONE_DIR}/plt00006 ${TWO_DIR}/plt00006
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/openmp_parity.log"
    ERROR_FILE "${WORKING_DIRECTORY}/openmp_parity.log"
    RESULT_VARIABLE parity_result)
if(NOT parity_result EQUAL 0)
    message(FATAL_ERROR "Cloud Chamber OpenMP parity failed: ${parity_result}")
endif()
