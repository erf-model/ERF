cmake_minimum_required(VERSION 3.24)

include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")

if(NOT DEFINED MPIEXEC OR NOT DEFINED MPIEXEC_NUMPROC_FLAG OR
   NOT DEFINED MPIEXEC_PREFLAGS OR NOT DEFINED TEST_EXE OR
   NOT DEFINED INPUT OR NOT DEFINED RUNTIME_OPTIONS OR
   NOT DEFINED WORKING_DIRECTORY OR NOT DEFINED LOG OR
   NOT DEFINED EXPECTED_MESSAGE OR NOT DEFINED EXPECTED_GUIDANCE_1 OR
   NOT DEFINED EXPECTED_GUIDANCE_2)
    message(FATAL_ERROR "RunShocRemovedTransportConfig.cmake missing required argument")
endif()

file(GLOB test_exe_candidates "${TEST_EXE}")
list(LENGTH test_exe_candidates test_exe_count)
if(NOT test_exe_count EQUAL 1)
    message(FATAL_ERROR
        "Native SHOC startup test executable pattern must resolve to exactly one file: ${TEST_EXE}")
endif()
list(GET test_exe_candidates 0 TEST_EXE)
if(NOT EXISTS "${INPUT}")
    message(FATAL_ERROR "Native SHOC startup test input is missing: ${INPUT}")
endif()

# MPIEXEC may be a multi-word command such as "flux run", so it is split and
# validated by the shared helper rather than used as a single argv[0].
erf_mpi_launcher_command(run_command
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunShocRemovedTransportConfig.cmake")
list(APPEND run_command "${TEST_EXE}" "${INPUT}")
if(NOT "${RUNTIME_OPTIONS}" STREQUAL "")
    separate_arguments(runtime_options UNIX_COMMAND "${RUNTIME_OPTIONS}")
    list(APPEND run_command ${runtime_options})
endif()

# These negative startup tests intentionally trigger amrex::Abort. Avoid
# expensive addr2line symbolization of the expected backtrace, particularly
# for large Debug/LTO CI executables.
list(APPEND run_command "amrex.call_addr2line=0")

execute_process(
    COMMAND ${run_command}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_VARIABLE simulation_stdout
    ERROR_VARIABLE simulation_stderr
    RESULT_VARIABLE simulation_result)

set(combined_output "${simulation_stdout}\n${simulation_stderr}")
file(WRITE "${LOG}"
    "RESULT=${simulation_result}\nINPUT=${INPUT}\n\n${combined_output}")

if("${simulation_result}" STREQUAL "0")
    message(FATAL_ERROR
        "Native SHOC removed transport option unexpectedly succeeded; see ${LOG}")
endif()

foreach(expected_text IN ITEMS "${EXPECTED_MESSAGE}" "${EXPECTED_GUIDANCE_1}" "${EXPECTED_GUIDANCE_2}")
    if(NOT "${expected_text}" STREQUAL "")
        string(FIND "${combined_output}" "${expected_text}" expected_index)
        if(expected_index EQUAL -1)
            message(FATAL_ERROR
                "Native SHOC startup failure did not contain expected text '${expected_text}'; see ${LOG}")
        endif()
    endif()
endforeach()

string(FIND "${combined_output}" "Coarse STEP" advancement_index)
if(NOT advancement_index EQUAL -1)
    message(FATAL_ERROR
        "Native SHOC removed transport option was rejected after time advancement; see ${LOG}")
endif()

string(FIND "${combined_output}" "Writing native 3D plotfile" output_index)
if(NOT output_index EQUAL -1)
    message(FATAL_ERROR
        "Native SHOC removed transport option was rejected after output; see ${LOG}")
endif()

message(STATUS "Native SHOC removed transport option rejected during startup")
