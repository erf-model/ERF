cmake_minimum_required(VERSION 3.24)

foreach(required_argument
    MPIEXEC MPIEXEC_NUMPROC_FLAG NRANKS TEST_EXE
    BASELINE_INPUT MUTANT_INPUT BASELINE_WORKING_DIRECTORY MUTANT_WORKING_DIRECTORY
    BASELINE_LOG MUTANT_LOG CHECKER TARGET_FIELD MIN_FINAL_DIFFERENCE
    MIN_BASELINE_EVOLUTION MAX_MUTANT_TO_BASELINE_EVOLUTION_RATIO)
    if(NOT DEFINED ${required_argument} OR "${${required_argument}}" STREQUAL "")
        message(FATAL_ERROR
            "RunShocMutationRegression.cmake requires ${required_argument}")
    endif()
endforeach()

function(run_mutation_case INPUT WORKING_DIRECTORY LOG OPTIONS RESULT_VARIABLE)
    set(command "${MPIEXEC}" "${MPIEXEC_NUMPROC_FLAG}" "${NRANKS}")
    if(DEFINED MPIEXEC_PREFLAGS AND NOT "${MPIEXEC_PREFLAGS}" STREQUAL "")
        separate_arguments(mpi_preflags UNIX_COMMAND "${MPIEXEC_PREFLAGS}")
        list(APPEND command ${mpi_preflags})
    endif()
    list(APPEND command "${TEST_EXE}" "${INPUT}")
    if(NOT "${OPTIONS}" STREQUAL "")
        separate_arguments(runtime_options UNIX_COMMAND "${OPTIONS}")
        list(APPEND command ${runtime_options})
    endif()

    execute_process(
        COMMAND ${command}
        WORKING_DIRECTORY "${WORKING_DIRECTORY}"
        OUTPUT_FILE "${LOG}"
        ERROR_FILE "${LOG}"
        RESULT_VARIABLE result)
    set(${RESULT_VARIABLE} "${result}" PARENT_SCOPE)
endfunction()

run_mutation_case(
    "${BASELINE_INPUT}"
    "${BASELINE_WORKING_DIRECTORY}"
    "${BASELINE_LOG}"
    "${BASELINE_OPTIONS}"
    baseline_result)
if(NOT baseline_result EQUAL 0)
    message(FATAL_ERROR "SHOC baseline simulation failed with exit code ${baseline_result}; see ${BASELINE_LOG}")
endif()

run_mutation_case(
    "${MUTANT_INPUT}"
    "${MUTANT_WORKING_DIRECTORY}"
    "${MUTANT_LOG}"
    "${MUTANT_OPTIONS}"
    mutant_result)
if(NOT mutant_result EQUAL 0)
    message(FATAL_ERROR "SHOC mutant simulation failed with exit code ${mutant_result}; see ${MUTANT_LOG}")
endif()

foreach(snapshot plt00000 plt00010 plt00020)
    if(NOT EXISTS "${BASELINE_WORKING_DIRECTORY}/${snapshot}")
        message(FATAL_ERROR "SHOC baseline output is missing: ${snapshot}")
    endif()
    if(NOT EXISTS "${MUTANT_WORKING_DIRECTORY}/${snapshot}")
        message(FATAL_ERROR "SHOC mutant output is missing: ${snapshot}")
    endif()
endforeach()

set(check_command "${MPIEXEC}" "${MPIEXEC_NUMPROC_FLAG}" "1")
if(DEFINED MPIEXEC_PREFLAGS AND NOT "${MPIEXEC_PREFLAGS}" STREQUAL "")
    separate_arguments(check_preflags UNIX_COMMAND "${MPIEXEC_PREFLAGS}")
    list(APPEND check_command ${check_preflags})
endif()
list(APPEND check_command
    "${CHECKER}"
    --field "${TARGET_FIELD}"
    --min-final-difference "${MIN_FINAL_DIFFERENCE}"
    --min-baseline-evolution "${MIN_BASELINE_EVOLUTION}"
    --max-mutant-to-baseline-evolution-ratio "${MAX_MUTANT_TO_BASELINE_EVOLUTION_RATIO}"
    --baseline-initial "${BASELINE_WORKING_DIRECTORY}/plt00000"
    --mutant-initial "${MUTANT_WORKING_DIRECTORY}/plt00000"
    --baseline-final "${BASELINE_WORKING_DIRECTORY}/plt00020"
    --mutant-final "${MUTANT_WORKING_DIRECTORY}/plt00020")
execute_process(COMMAND ${check_command} RESULT_VARIABLE checker_result)
if(NOT checker_result EQUAL 0)
    message(FATAL_ERROR "SHOC paired mutation checker failed with exit code ${checker_result}")
endif()
