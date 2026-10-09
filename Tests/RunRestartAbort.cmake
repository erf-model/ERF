# Run a deck to STEP_CHK writing a checkpoint there, then restart from it and require the
# restart to STOP with EXPECTED_MESSAGE in its output. For guards on the restart path, which
# add_test_abort cannot reach: that macro runs one leg from an inputs file, and there is no
# checkpoint to restart from until a first leg has written one.
#
# The two legs take separate rank counts. That is the point for the level-0 regrid guard,
# whose automatic branch fires on grids[0].size() < NProcs() -- i.e. only when the restart
# runs on more ranks than the checkpoint was written with.
#
# amrex.call_addr2line = 0 for the same reason as add_test_abort: the abort is the point of
# the test, and resolving every stack frame costs far longer than the run.
#
# -DX= defines X as empty, so test for a value, not for DEFINED
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

foreach(arg CHK_NRANKS RESTART_NRANKS TEST_EXE INPUT WORKING_DIRECTORY STEP_CHK EXPECTED_MESSAGE RUN_TIMEOUT)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunRestartAbort.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunRestartAbort.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunRestartAbort.cmake: ERF executable")

separate_arguments(common_options  UNIX_COMMAND "${COMMON_OPTIONS}")
separate_arguments(restart_options UNIX_COMMAND "${RESTART_OPTIONS}")

set(RUN_DIR "${WORKING_DIRECTORY}/run")
file(MAKE_DIRECTORY "${RUN_DIR}")

# A deck may name auxiliary inputs -- an input_sounding, a table, a terrain file --
# relatively, and the leg runs in its own subdirectory, so they have to be there too.
# Directories are skipped: an earlier run's checkpoints are not inputs.
file(GLOB _ra_aux "${WORKING_DIRECTORY}/*")
foreach(_ra_f IN LISTS _ra_aux)
    if(NOT IS_DIRECTORY "${_ra_f}")
        file(COPY "${_ra_f}" DESTINATION "${RUN_DIR}")
    endif()
endforeach()

erf_mpi_launcher_command(launch_chk
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${CHK_NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunRestartAbort.cmake")
erf_mpi_launcher_command(launch_restart
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${RESTART_NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunRestartAbort.cmake")

# checkpoint names carry the step padded to five digits
set(_s "0000${STEP_CHK}")
string(LENGTH "${_s}" _len)
math(EXPR _start "${_len} - 5")
string(SUBSTRING "${_s}" ${_start} 5 _s)
set(CHKFILE "chk${_s}")

# to the checkpoint step, writing it there
execute_process(
    COMMAND ${launch_chk} ${TEST_EXE} ${INPUT} ${common_options}
            "max_step=${STEP_CHK}" "erf.check_int=${STEP_CHK}" "erf.plot_int_1=-1"
    WORKING_DIRECTORY "${RUN_DIR}"
    OUTPUT_FILE "${RUN_DIR}/checkpoint.log"
    ERROR_FILE "${RUN_DIR}/checkpoint.log"
    TIMEOUT ${RUN_TIMEOUT}
    RESULT_VARIABLE chk_result)
if(NOT chk_result EQUAL 0)
    message(FATAL_ERROR "RunRestartAbort.cmake: the checkpoint run failed or exceeded ${RUN_TIMEOUT} s: ${chk_result} (see checkpoint.log)")
endif()
if(NOT EXISTS "${RUN_DIR}/${CHKFILE}/Header")
    message(FATAL_ERROR "RunRestartAbort.cmake: no ${CHKFILE} written by the checkpoint run")
endif()

# restart, and require it to stop
execute_process(
    COMMAND ${launch_restart} ${TEST_EXE} ${INPUT} ${common_options}
            "erf.restart=${CHKFILE}" "erf.check_int=-1" "erf.plot_int_1=-1"
            "amrex.call_addr2line=0" ${restart_options}
    WORKING_DIRECTORY "${RUN_DIR}"
    OUTPUT_FILE "${RUN_DIR}/restart.log"
    ERROR_FILE "${RUN_DIR}/restart.log"
    TIMEOUT ${RUN_TIMEOUT}
    RESULT_VARIABLE restart_result)

file(READ "${RUN_DIR}/restart.log" restart_output)

if(restart_result EQUAL 0)
    message(FATAL_ERROR
        "RunRestartAbort.cmake: the restart was expected to stop with\n"
        "  ${EXPECTED_MESSAGE}\n"
        "but it ran to completion (see restart.log)")
endif()
if(NOT restart_output MATCHES "${EXPECTED_MESSAGE}")
    message(FATAL_ERROR
        "RunRestartAbort.cmake: the restart stopped (${restart_result}) but without the "
        "expected message\n  ${EXPECTED_MESSAGE}\n"
        "A stop for the wrong reason passes nothing; see restart.log")
endif()

message(STATUS "RunRestartAbort: the restart stopped as required, matching '${EXPECTED_MESSAGE}'")
