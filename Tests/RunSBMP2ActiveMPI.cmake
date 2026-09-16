if(NOT DEFINED MPIEXEC OR NOT DEFINED MPIEXEC_NUMPROC_FLAG OR
   NOT DEFINED TEST_EXE OR NOT DEFINED INPUT OR NOT DEFINED WORKING_DIRECTORY OR
   NOT DEFINED LOG OR NOT DEFINED COMPARE)
    message(FATAL_ERROR "RunSBMP2ActiveMPI.cmake missing required argument")
endif()

set(_root "${WORKING_DIRECTORY}/active_mpi")
file(REMOVE_RECURSE "${_root}")
file(MAKE_DIRECTORY "${_root}")

set(_mpi_command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG})
if(DEFINED MPIEXEC_PREFLAGS AND NOT "${MPIEXEC_PREFLAGS}" STREQUAL "")
    separate_arguments(_mpi_preflags UNIX_COMMAND "${MPIEXEC_PREFLAGS}")
    list(APPEND _mpi_command ${_mpi_preflags})
endif()

foreach(_nranks IN ITEMS 1 2)
    set(_run_log "${_root}/run_${_nranks}r.log")
    set(_diagnostic "${_root}/run_${_nranks}r.composite")
    set(_command ${_mpi_command} ${_nranks} ${TEST_EXE} ${INPUT}
        erf.sbm_composite_diagnostic_file=${_diagnostic})
    execute_process(
        COMMAND ${_command}
        WORKING_DIRECTORY "${WORKING_DIRECTORY}"
        OUTPUT_FILE "${_run_log}"
        ERROR_FILE "${_run_log}"
        RESULT_VARIABLE _result)
    if(NOT _result EQUAL 0)
        message(FATAL_ERROR "active-limiter ${_nranks}-rank production run failed: ${_result}")
    endif()
    file(READ "${_diagnostic}" _diagnostic_text)
    foreach(_required IN ITEMS "format=erf-sbm-p2-composite-v1" "active_limiter=1" "passed=1")
        string(FIND "${_diagnostic_text}" "${_required}" _offset)
        if(_offset EQUAL -1)
            message(FATAL_ERROR "${_nranks}-rank active-limiter diagnostic is missing ${_required}")
        endif()
    endforeach()
endforeach()

execute_process(
    COMMAND python3 "${COMPARE}" "${_root}/run_1r.composite" "${_root}/run_2r.composite"
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${LOG}"
    ERROR_FILE "${LOG}"
    RESULT_VARIABLE _compare_result)
if(NOT _compare_result EQUAL 0)
    message(FATAL_ERROR "active-limiter MPI diagnostics differ: ${_compare_result}")
endif()
