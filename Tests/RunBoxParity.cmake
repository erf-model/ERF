# Run a deck twice, once on a single box on one rank and once on a BoxArray split into
# several boxes, and compare the final plotfiles with fcompare. The results must not depend
# on how the domain is decomposed. The reference runs on one rank because on more ranks
# amr.refine_grid_layout splits even a single box.
foreach(arg MPIEXEC_NUMPROC_FLAG NRANKS TEST_EXE INPUT WORKING_DIRECTORY FCOMPARE PLTFILE RTOL ATOL)
    if(NOT DEFINED ${arg})
        message(FATAL_ERROR "RunBoxParity.cmake missing required argument ${arg}")
    endif()
endforeach()

separate_arguments(common_options    UNIX_COMMAND "${COMMON_OPTIONS}")
separate_arguments(reference_options UNIX_COMMAND "${REFERENCE_OPTIONS}")
separate_arguments(split_options     UNIX_COMMAND "${SPLIT_OPTIONS}")
separate_arguments(mpiexec_preflags  UNIX_COMMAND "${MPIEXEC_PREFLAGS}")

set(REF_DIR   "${WORKING_DIRECTORY}/one_box")
set(SPLIT_DIR "${WORKING_DIRECTORY}/split")
file(REMOVE_RECURSE "${REF_DIR}" "${SPLIT_DIR}")
file(MAKE_DIRECTORY "${REF_DIR}" "${SPLIT_DIR}")

if(NOT "${MPIEXEC}" STREQUAL "")
    set(launch_one   ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1         ${mpiexec_preflags})
    set(launch_split ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NRANKS} ${mpiexec_preflags})
else()
    set(launch_one)
    set(launch_split)
endif()

execute_process(
    COMMAND ${launch_one} ${TEST_EXE} ${INPUT} ${common_options} ${reference_options}
    WORKING_DIRECTORY "${REF_DIR}"
    OUTPUT_FILE "${REF_DIR}/simulation.log"
    ERROR_FILE "${REF_DIR}/simulation.log"
    RESULT_VARIABLE ref_result)
if(NOT ref_result EQUAL 0)
    message(FATAL_ERROR "Single-box simulation failed: ${ref_result}")
endif()

execute_process(
    COMMAND ${launch_split} ${TEST_EXE} ${INPUT} ${common_options} ${split_options}
    WORKING_DIRECTORY "${SPLIT_DIR}"
    OUTPUT_FILE "${SPLIT_DIR}/simulation.log"
    ERROR_FILE "${SPLIT_DIR}/simulation.log"
    RESULT_VARIABLE split_result)
if(NOT split_result EQUAL 0)
    message(FATAL_ERROR "Split-box simulation failed: ${split_result}")
endif()

execute_process(
    COMMAND ${launch_one} ${FCOMPARE} --abort_if_not_all_found --allow_diff_grids
            --rel_tol ${RTOL} --abs_tol ${ATOL}
            ${REF_DIR}/${PLTFILE} ${SPLIT_DIR}/${PLTFILE}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/parity.log"
    ERROR_FILE "${WORKING_DIRECTORY}/parity.log"
    RESULT_VARIABLE parity_result)
if(NOT parity_result EQUAL 0)
    message(FATAL_ERROR "Single-box and split-box plotfiles differ: ${parity_result}")
endif()
