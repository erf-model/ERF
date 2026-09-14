if(NOT DEFINED MPIEXEC OR NOT DEFINED MPIEXEC_NUMPROC_FLAG OR
   NOT DEFINED NRANKS OR NOT DEFINED TEST_EXE OR NOT DEFINED INPUT OR
   NOT DEFINED WORKING_DIRECTORY OR NOT DEFINED LOG OR
   NOT DEFINED DIAGNOSTIC OR NOT DEFINED CHECKER OR
   NOT DEFINED EXPECTED_COMPONENTS)
    message(FATAL_ERROR "RunSBMPrototype.cmake missing required argument")
endif()

set(_command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NRANKS})
if(DEFINED MPIEXEC_PREFLAGS AND NOT "${MPIEXEC_PREFLAGS}" STREQUAL "")
    separate_arguments(_mpi_preflags UNIX_COMMAND "${MPIEXEC_PREFLAGS}")
    list(APPEND _command ${_mpi_preflags})
endif()
list(APPEND _command ${TEST_EXE} ${INPUT})
if(DEFINED RUNTIME_OPTIONS AND NOT "${RUNTIME_OPTIONS}" STREQUAL "")
    separate_arguments(_runtime_options UNIX_COMMAND "${RUNTIME_OPTIONS}")
    list(APPEND _command ${_runtime_options})
endif()

execute_process(
    COMMAND ${_command}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${LOG}"
    ERROR_FILE "${LOG}"
    RESULT_VARIABLE simulation_result)
if(NOT simulation_result EQUAL 0)
    message(FATAL_ERROR "SBM P1 simulation failed: ${simulation_result}")
endif()

file(READ "${LOG}" simulation_log)
foreach(expected_text
        "SBM layout identity"
        "SBM P1 auxiliary state"
        "components=${EXPECTED_COMPONENTS}"
        "Coarse STEP 1 ends"
        "Coarse STEP 2 ends")
    string(FIND "${simulation_log}" "${expected_text}" expected_offset)
    if(expected_offset EQUAL -1)
        message(FATAL_ERROR "SBM P1 log is missing expected text: ${expected_text}")
    endif()
endforeach()

if(NOT EXISTS "${DIAGNOSTIC}")
    message(FATAL_ERROR "SBM P1 simulation did not produce numerical diagnostic: ${DIAGNOSTIC}")
endif()
execute_process(
    COMMAND "${CHECKER}" "${DIAGNOSTIC}"
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    RESULT_VARIABLE checker_result)
if(NOT checker_result EQUAL 0)
    message(FATAL_ERROR "SBM P1 numerical qualification failed: ${checker_result}")
endif()
