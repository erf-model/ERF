cmake_minimum_required(VERSION 3.24)

include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")

# Observation-nudging regression driver.
#
# MODE = analytic : run the deck once, then run the checker on the station
#                   series for every check in CHECKS.
# MODE = approach : run the deck with nudging (in on/) and with OFF_OPTIONS
#                   (in off/), then run the checker for every check in CHECKS
#                   with on= and off= set to the station series of each run
#                   that the check names with series=<station>.
# MODE = single   : run the deck once, then run the checker with every check
#                   in CHECKS as it is, @RUN@ replaced by the run directory.
# MODE = abort    : run the deck with RUNTIME_OPTIONS and require that it
#                   aborts during start-up with EXPECTED_MESSAGE in its output.
#
# CHECKS is a list of checker argument strings separated by '|', each one the
# key=value arguments of one checker call (ObsNudgingCheck.cpp).

if(NOT DEFINED TEST_EXE OR NOT DEFINED INPUT OR NOT DEFINED WORKING_DIRECTORY OR
   NOT DEFINED MODE OR NOT DEFINED LOG)
    message(FATAL_ERROR "RunObsNudging.cmake missing required argument")
endif()
if(NOT DEFINED NRANKS OR "${NRANKS}" STREQUAL "")
    set(NRANKS 1)
endif()

file(GLOB test_exe_candidates "${TEST_EXE}")
list(LENGTH test_exe_candidates test_exe_count)
if(NOT test_exe_count EQUAL 1)
    message(FATAL_ERROR "The test executable pattern must resolve to exactly one file: ${TEST_EXE}")
endif()
list(GET test_exe_candidates 0 TEST_EXE)

function(obs_report_log label path)
    if(EXISTS "${path}")
        file(READ "${path}" contents)
        message(STATUS "---- ${label} (${path}) ----\n${contents}\n---- end ${label} ----")
    else()
        message(STATUS "---- ${label}: ${path} was never written ----")
    endif()
endfunction()

function(obs_launcher nranks out_var)
    erf_mpi_launcher_command(launcher
        LAUNCHER "${MPIEXEC}"
        NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
        NRANKS ${nranks}
        PREFLAGS "${MPIEXEC_PREFLAGS}"
        CONTEXT "RunObsNudging.cmake")
    set(${out_var} "${launcher}" PARENT_SCOPE)
endfunction()

# Run the deck in dir with the extra options opts; the result goes to out_var
function(obs_run dir opts log out_var)
    file(REMOVE_RECURSE "${dir}")
    file(MAKE_DIRECTORY "${dir}")
    # The deck's data files (station files, soundings), not the outputs of
    # earlier legs: every regular file except the deck and the logs
    file(GLOB inputs LIST_DIRECTORIES false "${WORKING_DIRECTORY}/*")
    list(FILTER inputs EXCLUDE REGEX "\\.(i|log|log\\..*)$")
    if(inputs)
        file(COPY ${inputs} DESTINATION "${dir}")
    endif()
    obs_launcher(${NRANKS} launcher)
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

function(obs_check args)
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
        message(FATAL_ERROR "Observation-nudging check failed: ${args}")
    endif()
endfunction()

file(REMOVE "${LOG}.checker")

if("${MODE}" STREQUAL "abort")
    obs_run("${WORKING_DIRECTORY}/abort" "${RUNTIME_OPTIONS}" "${LOG}" result)
    file(READ "${LOG}" output)
    if("${result}" STREQUAL "0")
        obs_report_log("simulation log" "${LOG}")
        message(FATAL_ERROR "The run was expected to abort but succeeded")
    endif()
    string(FIND "${output}" "${EXPECTED_MESSAGE}" found)
    if(found EQUAL -1)
        obs_report_log("simulation log" "${LOG}")
        message(FATAL_ERROR "The abort did not contain the expected text '${EXPECTED_MESSAGE}'")
    endif()
    string(FIND "${output}" "Coarse STEP" stepped)
    if(NOT stepped EQUAL -1)
        obs_report_log("simulation log" "${LOG}")
        message(FATAL_ERROR "The run aborted only after it started time stepping")
    endif()
    message(STATUS "Aborted during start-up with: ${EXPECTED_MESSAGE}")
    return()
endif()

if(NOT DEFINED CHECKER OR NOT DEFINED CHECKS)
    message(FATAL_ERROR "RunObsNudging.cmake: MODE ${MODE} needs CHECKER and CHECKS")
endif()
string(REPLACE "|" ";" check_list "${CHECKS}")

if("${MODE}" STREQUAL "analytic")
    obs_run("${WORKING_DIRECTORY}/run" "${RUNTIME_OPTIONS}" "${LOG}" result)
    if(NOT result EQUAL 0)
        obs_report_log("simulation log" "${LOG}")
        message(FATAL_ERROR "The simulation failed: ${result}")
    endif()
    foreach(check IN LISTS check_list)
        obs_check("analytic file=${WORKING_DIRECTORY}/run/Output_Stations/${STATION}.dat ${check}")
    endforeach()
elseif("${MODE}" STREQUAL "single")
    obs_run("${WORKING_DIRECTORY}/run" "${RUNTIME_OPTIONS}" "${LOG}" result)
    if(NOT result EQUAL 0)
        obs_report_log("simulation log" "${LOG}")
        message(FATAL_ERROR "The simulation failed: ${result}")
    endif()
    foreach(check IN LISTS check_list)
        string(REPLACE "@RUN@" "${WORKING_DIRECTORY}/run" check "${check}")
        obs_check("${check}")
    endforeach()
elseif("${MODE}" STREQUAL "approach")
    obs_run("${WORKING_DIRECTORY}/on" "${RUNTIME_OPTIONS}" "${LOG}" result)
    if(NOT result EQUAL 0)
        obs_report_log("simulation log" "${LOG}")
        message(FATAL_ERROR "The nudged simulation failed: ${result}")
    endif()
    obs_run("${WORKING_DIRECTORY}/off" "${RUNTIME_OPTIONS} ${OFF_OPTIONS}" "${LOG}.off" result)
    if(NOT result EQUAL 0)
        obs_report_log("simulation log" "${LOG}.off")
        message(FATAL_ERROR "The free simulation failed: ${result}")
    endif()
    foreach(check IN LISTS check_list)
        if(NOT check MATCHES "series=([A-Za-z0-9_]+)")
            message(FATAL_ERROR "RunObsNudging.cmake: approach check without series=: ${check}")
        endif()
        set(series "${CMAKE_MATCH_1}")
        string(REGEX REPLACE "series=[A-Za-z0-9_]+" "" check "${check}")
        obs_check("approach on=${WORKING_DIRECTORY}/on/Output_Stations/${series}.dat off=${WORKING_DIRECTORY}/off/Output_Stations/${series}.dat ${check}")
    endforeach()
else()
    message(FATAL_ERROR "RunObsNudging.cmake: unknown MODE ${MODE}")
endif()
