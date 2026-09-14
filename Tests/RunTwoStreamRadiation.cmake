if(NOT DEFINED TEST_EXE OR NOT DEFINED INPUT OR
   NOT DEFINED WORKING_DIRECTORY OR NOT DEFINED SIMULATION_LOG OR
   NOT DEFINED CHECKER_LOG OR NOT DEFINED CHECKER OR NOT DEFINED PLOTFILE)
    message(FATAL_ERROR "RunTwoStreamRadiation.cmake missing required argument")
endif()

# Build the launcher prefix. A build without MPI passes an empty MPIEXEC, in
# which case the executables are run directly instead of through a launcher.
function(two_stream_launcher nranks out_var)
    if(DEFINED MPIEXEC AND NOT "${MPIEXEC}" STREQUAL "")
        set(launcher ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${nranks} ${MPIEXEC_PREFLAGS})
    else()
        set(launcher "")
    endif()
    set(${out_var} "${launcher}" PARENT_SCOPE)
endfunction()

# Echo a log to the test output so a CI failure is diagnosable without having
# to fetch the attached files.
function(two_stream_report_log label path)
    if(EXISTS "${path}")
        file(READ "${path}" contents)
        message(STATUS "---- ${label} (${path}) ----\n${contents}\n---- end ${label} ----")
    else()
        message(STATUS "---- ${label}: ${path} was never written ----")
    endif()
endfunction()

if(NOT DEFINED NRANKS OR "${NRANKS}" STREQUAL "")
    set(NRANKS 1)
endif()

# Extra command-line inputs for the simulation (optional).
if(DEFINED RUNTIME_OPTIONS AND NOT "${RUNTIME_OPTIONS}" STREQUAL "")
    separate_arguments(runtime_options NATIVE_COMMAND "${RUNTIME_OPTIONS}")
else()
    set(runtime_options "")
endif()

# The domain-mean diagnostics CSV must not depend on the decomposition: the
# NRANKS run writes one file, a 1-rank run of the same deck writes another,
# and the two must match byte for byte. Rank-local means (no MPI reduction)
# fail this whenever the columns differ, as they do over terrain.
set(diag_nranks "${WORKING_DIRECTORY}/radiation_diag_np${NRANKS}.dat")
set(diag_serial "${WORKING_DIRECTORY}/radiation_diag_np1.dat")
file(REMOVE "${diag_nranks}" "${diag_serial}")
set(diag_options erf.radiation.diag_csv_enable=true erf.radiation.diag_enable=true)

two_stream_launcher(${NRANKS} simulation_launcher)
execute_process(
    COMMAND ${simulation_launcher} ${TEST_EXE} ${INPUT} ${runtime_options}
            ${diag_options} erf.radiation.diag_file=${diag_nranks}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${SIMULATION_LOG}"
    ERROR_FILE "${SIMULATION_LOG}"
    RESULT_VARIABLE simulation_result)
if(NOT simulation_result EQUAL 0)
    two_stream_report_log("simulation log" "${SIMULATION_LOG}")
    message(FATAL_ERROR "TwoStream radiation simulation failed: ${simulation_result}")
endif()

if(NRANKS GREATER 1)
    two_stream_launcher(1 serial_launcher)
    execute_process(
        COMMAND ${serial_launcher} ${TEST_EXE} ${INPUT} ${runtime_options}
                ${diag_options} erf.radiation.diag_file=${diag_serial}
                erf.plot_int_1=-1
        WORKING_DIRECTORY "${WORKING_DIRECTORY}"
        OUTPUT_FILE "${SIMULATION_LOG}.np1"
        ERROR_FILE "${SIMULATION_LOG}.np1"
        RESULT_VARIABLE serial_result)
    if(NOT serial_result EQUAL 0)
        two_stream_report_log("1-rank simulation log" "${SIMULATION_LOG}.np1")
        message(FATAL_ERROR "TwoStream radiation 1-rank simulation failed: ${serial_result}")
    endif()
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E compare_files "${diag_nranks}" "${diag_serial}"
        RESULT_VARIABLE csv_result)
    if(NOT csv_result EQUAL 0)
        two_stream_report_log("${NRANKS}-rank diagnostics" "${diag_nranks}")
        two_stream_report_log("1-rank diagnostics" "${diag_serial}")
        message(FATAL_ERROR "TwoStream radiation diagnostics differ between ${NRANKS} ranks and 1 rank")
    endif()
endif()

two_stream_launcher(1 checker_launcher)
execute_process(
    COMMAND ${checker_launcher} ${CHECKER} ${PLOTFILE}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${CHECKER_LOG}"
    ERROR_FILE "${CHECKER_LOG}"
    RESULT_VARIABLE checker_result)
if(NOT checker_result EQUAL 0)
    two_stream_report_log("checker log" "${CHECKER_LOG}")
    message(FATAL_ERROR "TwoStream radiation column check failed: ${checker_result}")
endif()
