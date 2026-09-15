if(NOT DEFINED MPIEXEC OR NOT DEFINED MPIEXEC_NUMPROC_FLAG OR
   NOT DEFINED NRANKS OR NOT DEFINED TEST_EXE OR NOT DEFINED INPUT OR
   NOT DEFINED WORKING_DIRECTORY OR NOT DEFINED LOG)
    message(FATAL_ERROR "RunSBMP2AMR.cmake missing required argument")
endif()

set(_command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NRANKS})
if(DEFINED MPIEXEC_PREFLAGS AND NOT "${MPIEXEC_PREFLAGS}" STREQUAL "")
    separate_arguments(_mpi_preflags UNIX_COMMAND "${MPIEXEC_PREFLAGS}")
    list(APPEND _command ${_mpi_preflags})
endif()
list(APPEND _command ${TEST_EXE} ${INPUT})

execute_process(
    COMMAND ${_command}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${LOG}"
    ERROR_FILE "${LOG}"
    RESULT_VARIABLE simulation_result)
if(NOT simulation_result EQUAL 0)
    message(FATAL_ERROR "SBM P2 AMR simulation failed: ${simulation_result}")
endif()

file(READ "${LOG}" simulation_log)
foreach(expected_text
        "SBM layout identity"
        "ncomp=8"
        "moment=1"
        "Coarse STEP 2 ends")
    string(FIND "${simulation_log}" "${expected_text}" expected_offset)
    if(expected_offset EQUAL -1)
        message(FATAL_ERROR "SBM P2 AMR log is missing expected text: ${expected_text}")
    endif()
endforeach()
