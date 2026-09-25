include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")

if(NOT DEFINED TEST_EXE OR NOT DEFINED INPUT OR
   NOT DEFINED WORKING_DIRECTORY OR NOT DEFINED SIMULATION_LOG OR
   NOT DEFINED CHECKER_LOG OR NOT DEFINED CHECKER OR NOT DEFINED PLOTFILE)
    message(FATAL_ERROR "RunTwoStreamRadiation.cmake missing required argument")
endif()

# Build the launcher prefix. A build without MPI passes an empty MPIEXEC, in
# which case the executables are run directly instead of through a launcher.
function(two_stream_launcher nranks out_var)
    erf_mpi_launcher_command(launcher
        LAUNCHER "${MPIEXEC}"
        NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
        NRANKS ${nranks}
        PREFLAGS "${MPIEXEC_PREFLAGS}"
        CONTEXT "RunTwoStreamRadiation.cmake")
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

# Levels to run the column check on, as a comma-separated string; default "0". Tests name
# them as a CMake list (CHECK_LEVELS 0 1) and add_test_two_stream_radiation joins them with
# commas on the way here -- one spelling per layer, and the comma is why:
# not just the coarse one's; each level gets the same assertions. Comma rather
# than semicolon because a semicolon in a -D argument is split when the
# COMMAND is built, which silently reduced this to level 0 only.
if(NOT DEFINED CHECK_LEVELS OR "${CHECK_LEVELS}" STREQUAL "")
    set(CHECK_LEVELS "0")
endif()
string(REPLACE "," ";" check_level_list "${CHECK_LEVELS}")
# Levels required to appear in the diagnostics CSV; defaults to the checked levels.
if(NOT DEFINED DIAG_LEVELS OR "${DIAG_LEVELS}" STREQUAL "")
    set(DIAG_LEVELS "${CHECK_LEVELS}")
endif()
string(REPLACE "," ";" diag_level_list "${DIAG_LEVELS}")
message(STATUS "TwoStream column check will run on level(s): ${check_level_list}")

# The diagnostics CSV must actually carry a row for every level that was checked.
#
# The 1-rank vs NRANKS byte-comparison above cannot see this: it passes just as well if
# every fine-level row is dropped in both runs, which is precisely the defect the per-level
# writer fixes (one shared writer deduplicates on (step, call_site, time), which every level
# reports identically). So assert the content, not only that two runs agree.
# A missing file is a failure, not a skip: the runner forces
# erf.radiation.diag_csv_enable=true above, so no file at all means the writer produced
# nothing -- a superset of the defect this block exists to catch, and exactly the "passes
# vacuously" mode the CHECK_LEVELS plumbing was fixed for.
if(NOT EXISTS "${diag_nranks}")
    message(FATAL_ERROR
        "TwoStream diagnostics CSV ${diag_nranks} was not written, though the runner forces "
        "erf.radiation.diag_csv_enable=true; the diagnostics writer produced nothing")
endif()
file(STRINGS "${diag_nranks}" diag_lines)
list(POP_FRONT diag_lines diag_header)
if(NOT diag_header MATCHES ",level$")
    message(FATAL_ERROR
        "TwoStream diagnostics CSV header does not end with the level column: ${diag_header}")
endif()
foreach(check_level IN LISTS diag_level_list)
    set(found_level FALSE)
    foreach(row IN LISTS diag_lines)
        if(row MATCHES ",${check_level}$")
            set(found_level TRUE)
            break()
        endif()
    endforeach()
    if(NOT found_level)
        message(FATAL_ERROR
            "TwoStream diagnostics CSV has no row for level ${check_level}; the per-level "
            "writer is not emitting one row set per level")
    endif()
endforeach()
message(STATUS "TwoStream diagnostics CSV carries rows for level(s): ${diag_level_list}")

two_stream_launcher(1 checker_launcher)
foreach(check_level IN LISTS check_level_list)
    set(level_log "${CHECKER_LOG}.lev${check_level}")
    execute_process(
        COMMAND ${checker_launcher} ${CHECKER} ${PLOTFILE} ${check_level}
        WORKING_DIRECTORY "${WORKING_DIRECTORY}"
        OUTPUT_FILE "${level_log}"
        ERROR_FILE "${level_log}"
        RESULT_VARIABLE checker_result)
    # Keep the un-suffixed log too, so an existing single-level test's attached
    # file name still resolves.
    configure_file("${level_log}" "${CHECKER_LOG}" COPYONLY)
    if(NOT checker_result EQUAL 0)
        two_stream_report_log("checker log (level ${check_level})" "${level_log}")
        message(FATAL_ERROR
            "TwoStream radiation column check failed at level ${check_level}: ${checker_result}")
    endif()
endforeach()
