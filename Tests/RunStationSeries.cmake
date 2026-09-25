cmake_minimum_required(VERSION 3.24)

include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

# Station time-series regression driver.
#
# MODE = single   : run the deck once in run/, then run the checker with every
#                   check in CHECKS as it is, @RUN@ replaced by the run directory.
# MODE = analytic : run the deck once, then run the checker in analytic mode on
#                   the series STATION for every check in CHECKS.
# MODE = approach : run the deck with nudging (in on/) and with OFF_OPTIONS
#                   (in off/), then run the checker in approach mode for every
#                   check in CHECKS with on= and off= set to the station series
#                   of each run that the check names with series=<station>.
# MODE = abort    : run the deck with RUNTIME_OPTIONS and require that it
#                   aborts during start-up with EXPECTED_MESSAGE in its output.
#
# CHECKS is a list of checker argument strings separated by '|'.  In single
# mode each is the mode and key=value arguments of one checker call
# (StationSeriesCheck.cpp); in analytic and approach mode, the key=value
# arguments only.

# -DX= defines X as empty, so test for a value, not for DEFINED: an empty CHECKS
# would otherwise run the deck and report success without calling the checker
# once, and an empty EXPECTED_MESSAGE would accept any abort at all
set(required TEST_EXE INPUT WORKING_DIRECTORY MODE LOG)
if("${MODE}" STREQUAL "abort")
    list(APPEND required EXPECTED_MESSAGE)
else()
    list(APPEND required CHECKER CHECKS)
    if("${MODE}" STREQUAL "analytic")
        list(APPEND required STATION)
    elseif("${MODE}" STREQUAL "approach")
        list(APPEND required OFF_OPTIONS)
    endif()
endif()
foreach(arg IN LISTS required)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunStationSeries.cmake: ${arg} must be given and non-empty for MODE ${MODE}")
    endif()
endforeach()
if("${NRANKS}" STREQUAL "")
    set(NRANKS 1)
endif()

# On Windows the executable is named with a wildcard for the config subdirectory
# a multi-config generator picks; execute_process does not expand it
erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunStationSeries.cmake: ERF executable")

function(station_report_log label path)
    if(EXISTS "${path}")
        file(READ "${path}" contents)
        message(STATUS "---- ${label} (${path}) ----\n${contents}\n---- end ${label} ----")
    else()
        message(STATUS "---- ${label}: ${path} was never written ----")
    endif()
endfunction()

# Run the deck in dir with the extra options opts; the result goes to out_var
function(station_run dir opts log out_var)
    file(REMOVE_RECURSE "${dir}")
    file(MAKE_DIRECTORY "${dir}")
    # The deck's data files (soundings and the like), not the outputs of
    # earlier runs: every regular file except the deck and the logs
    file(GLOB inputs LIST_DIRECTORIES false "${WORKING_DIRECTORY}/*")
    list(FILTER inputs EXCLUDE REGEX "\\.(i|log|log\\..*)$")
    if(inputs)
        file(COPY ${inputs} DESTINATION "${dir}")
    endif()
    erf_mpi_launcher_command(launcher
        LAUNCHER "${MPIEXEC}"
        NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
        NRANKS ${NRANKS}
        PREFLAGS "${MPIEXEC_PREFLAGS}"
        CONTEXT "RunStationSeries.cmake")
    set(options "")
    if(NOT "${opts}" STREQUAL "")
        separate_arguments(options UNIX_COMMAND "${opts}")
    endif()
    execute_process(
        COMMAND ${launcher} ${TEST_EXE} ${INPUT} ${options} amrex.call_addr2line=0
        WORKING_DIRECTORY "${dir}"
        OUTPUT_FILE "${log}"
        ERROR_FILE "${log}"
        RESULT_VARIABLE result)
    set(${out_var} "${result}" PARENT_SCOPE)
endfunction()

function(station_check args)
    separate_arguments(check_args UNIX_COMMAND "${args}")
    execute_process(
        COMMAND ${CHECKER} ${check_args}
        WORKING_DIRECTORY "${WORKING_DIRECTORY}"
        OUTPUT_VARIABLE out
        ERROR_VARIABLE err
        RESULT_VARIABLE result)
    message(STATUS "${out}${err}")
    file(APPEND "${LOG}.checker" "${out}${err}")
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "Station-series check failed (${result}): ${CHECKER} ${args}")
    endif()
endfunction()

file(REMOVE "${LOG}.checker")

if("${MODE}" STREQUAL "abort")
    station_run("${WORKING_DIRECTORY}/abort" "${RUNTIME_OPTIONS}" "${LOG}" result)
    file(READ "${LOG}" output)
    if("${result}" STREQUAL "0")
        station_report_log("simulation log" "${LOG}")
        message(FATAL_ERROR "The run was expected to abort but succeeded")
    endif()
    string(FIND "${output}" "${EXPECTED_MESSAGE}" found)
    if(found EQUAL -1)
        station_report_log("simulation log" "${LOG}")
        message(FATAL_ERROR "The abort did not contain the expected text '${EXPECTED_MESSAGE}'")
    endif()
    string(FIND "${output}" "Coarse STEP" stepped)
    if(NOT stepped EQUAL -1)
        station_report_log("simulation log" "${LOG}")
        message(FATAL_ERROR "The run aborted only after it started time stepping")
    endif()
    message(STATUS "Aborted during start-up with: ${EXPECTED_MESSAGE}")
    return()
endif()

string(REPLACE "|" ";" check_list "${CHECKS}")

if("${MODE}" STREQUAL "single")
    station_run("${WORKING_DIRECTORY}/run" "${RUNTIME_OPTIONS}" "${LOG}" result)
    if(NOT result EQUAL 0)
        station_report_log("simulation log" "${LOG}")
        message(FATAL_ERROR "The simulation failed: ${result}")
    endif()
    foreach(check IN LISTS check_list)
        string(REPLACE "@RUN@" "${WORKING_DIRECTORY}/run" check "${check}")
        station_check("${check}")
    endforeach()
elseif("${MODE}" STREQUAL "analytic")
    station_run("${WORKING_DIRECTORY}/run" "${RUNTIME_OPTIONS}" "${LOG}" result)
    if(NOT result EQUAL 0)
        station_report_log("simulation log" "${LOG}")
        message(FATAL_ERROR "The simulation failed: ${result}")
    endif()
    foreach(check IN LISTS check_list)
        station_check("analytic file=${WORKING_DIRECTORY}/run/Output_Stations/${STATION}.dat ${check}")
    endforeach()
elseif("${MODE}" STREQUAL "approach")
    station_run("${WORKING_DIRECTORY}/on" "${RUNTIME_OPTIONS}" "${LOG}" result)
    if(NOT result EQUAL 0)
        station_report_log("simulation log" "${LOG}")
        message(FATAL_ERROR "The nudged simulation failed: ${result}")
    endif()
    station_run("${WORKING_DIRECTORY}/off" "${RUNTIME_OPTIONS} ${OFF_OPTIONS}" "${LOG}.off" result)
    if(NOT result EQUAL 0)
        station_report_log("simulation log" "${LOG}.off")
        message(FATAL_ERROR "The free simulation failed: ${result}")
    endif()
    foreach(check IN LISTS check_list)
        if(NOT check MATCHES "series=([A-Za-z0-9_]+)")
            message(FATAL_ERROR "RunStationSeries.cmake: approach check without series=: ${check}")
        endif()
        set(series "${CMAKE_MATCH_1}")
        string(REGEX REPLACE "series=[A-Za-z0-9_]+" "" check "${check}")
        station_check("approach on=${WORKING_DIRECTORY}/on/Output_Stations/${series}.dat off=${WORKING_DIRECTORY}/off/Output_Stations/${series}.dat ${check}")
    endforeach()
else()
    message(FATAL_ERROR "RunStationSeries.cmake: unknown MODE ${MODE}")
endif()
