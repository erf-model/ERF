# Run a deck straight to STEP_END, run it again to STEP_CHK with a checkpoint there,
# restart from that checkpoint to STEP_END, and require the restarted run's plotfile
# at STEP_END to equal the straight run's with fcompare. Each leg uses the forwarded
# RUN_TIMEOUT, and the enclosing CTest timeout is sized separately by the caller.
# COMMON_OPTIONS goes to all three legs; RESTART_OPTIONS goes to the restart leg only.
# CHK_NRANKS/RESTART_NRANKS let the checkpoint and restart legs run at different widths;
# REQUIRE_LEVEL0_REMAKE asserts the restart really did re-make the level-0 grids.
# -DX= defines X as empty, so test for a value, not for DEFINED
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

foreach(arg NRANKS TEST_EXE INPUT WORKING_DIRECTORY FCOMPARE STEP_CHK STEP_END RTOL ATOL RUN_TIMEOUT)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunRestartParity.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunRestartParity.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

# On Windows the executables are named with a wildcard for the config subdirectory
# a multi-config generator picks; execute_process does not expand it.
erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunRestartParity.cmake: ERF executable")
erf_resolve_executable(FCOMPARE "${FCOMPARE}" CONFIG "${CONFIG}"
    CONTEXT "RunRestartParity.cmake: fcompare")

separate_arguments(common_options   UNIX_COMMAND "${COMMON_OPTIONS}")
#
# Options for the restart leg only.  COMMON_OPTIONS cannot carry anything that changes the
# decomposition, because it goes to the straight leg as well -- that would make this a
# box-parity test instead of a restart-parity one.  Optional, so deliberately not in the
# required-argument loop above.
#
separate_arguments(restart_options  UNIX_COMMAND "${RESTART_OPTIONS}")

#
# A restart leg that re-makes the level-0 grids writes its plotfile on a different
# BoxArray, and fcompare refuses to compare those without being told to.  Opt in, so the
# other restart-parity tests keep failing if their grids ever move -- for them a changed
# decomposition is a bug, not the point of the test.
#
set(fcompare_grids "")
if(ALLOW_DIFF_GRIDS)
    set(fcompare_grids "--allow_diff_grids")
endif()

set(STRAIGHT_DIR "${WORKING_DIRECTORY}/straight")
set(RESTART_DIR  "${WORKING_DIRECTORY}/restart")
file(REMOVE_RECURSE "${STRAIGHT_DIR}" "${RESTART_DIR}")
file(MAKE_DIRECTORY "${STRAIGHT_DIR}" "${RESTART_DIR}")

# A deck may need auxiliary inputs sitting beside it -- an input_sounding, a table, a
# terrain file -- which it names relatively. Each leg runs in its own subdirectory, so
# those files have to be there too; otherwise the run aborts at start-up reading them.
# Directories are skipped: the plotfiles and checkpoints of an earlier run are not inputs.
file(GLOB _rp_aux "${WORKING_DIRECTORY}/*")
foreach(_rp_f IN LISTS _rp_aux)
    if(NOT IS_DIRECTORY "${_rp_f}")
        file(COPY "${_rp_f}" DESTINATION "${STRAIGHT_DIR}")
        file(COPY "${_rp_f}" DESTINATION "${RESTART_DIR}")
    endif()
endforeach()

# MPIEXEC may be a multi-word command such as "flux run"; the helper splits
# it, validates the program and applies MPIEXEC_PREFLAGS. An empty MPIEXEC
# yields an empty prefix, so the runs stay serial.
#
# The checkpoint leg may run on a different number of ranks from the restart, which is its
# own code path: ERF::restart re-makes the level-0 grids by itself when the checkpoint has
# fewer level-0 boxes than there are ranks, so a run continued on more ranks than it was
# written with takes the regrid branch without anyone asking for it. Default both to NRANKS
# so every existing caller is unchanged.
#
# The STRAIGHT leg runs at the restart's rank count, not the checkpoint's, so that it and the
# restart share a decomposition and the comparison isolates the restart itself.
#
set(_chk_nranks "${NRANKS}")
set(_restart_nranks "${NRANKS}")
if(NOT "${CHK_NRANKS}" STREQUAL "")
    set(_chk_nranks "${CHK_NRANKS}")
endif()
if(NOT "${RESTART_NRANKS}" STREQUAL "")
    set(_restart_nranks "${RESTART_NRANKS}")
endif()

erf_mpi_launcher_command(launch_chk
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${_chk_nranks}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunRestartParity.cmake")
erf_mpi_launcher_command(launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${_restart_nranks}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunRestartParity.cmake")
erf_mpi_launcher_command(launch_one
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunRestartParity.cmake")

# plotfile and checkpoint names carry the step padded to five digits
function(padded step out_var)
    set(_s "0000${step}")
    string(LENGTH "${_s}" _len)
    math(EXPR _start "${_len} - 5")
    string(SUBSTRING "${_s}" ${_start} 5 _s)
    set(${out_var} "${_s}" PARENT_SCOPE)
endfunction()
padded(${STEP_CHK} chk_step)
padded(${STEP_END} end_step)
set(PLTFILE "plt${end_step}")
set(CHKFILE "chk${chk_step}")

function(run_erf_with launcher dir log timeout_s)
    execute_process(
        COMMAND ${launcher} ${TEST_EXE} ${INPUT} ${common_options} ${ARGN}
        WORKING_DIRECTORY "${dir}"
        OUTPUT_FILE "${dir}/${log}"
        ERROR_FILE "${dir}/${log}"
        TIMEOUT ${timeout_s}
        RESULT_VARIABLE _result)
    if(NOT _result EQUAL 0)
        message(FATAL_ERROR "RunRestartParity.cmake: the run in ${dir} (${log}) failed or exceeded ${timeout_s} s: ${_result}")
    endif()
endfunction()

# Optional second comparison: a 2D plotfile. The 3D plotfile carries no surface
# fields, so state that is checkpointed but never restored -- a surface temperature
# reloading as its scalar default, say -- can leave plt identical while plt2d is
# wrong. Drive its cadence exactly as the 3D one is driven.
set(plot2d_end "")
set(plot2d_off "")
if(NOT "${PLT2DFILE}" STREQUAL "")
    set(plot2d_end "erf.plot2d_int_1=${STEP_END}")
    set(plot2d_off "erf.plot2d_int_1=-1")
endif()

# straight to the end, no checkpoint
run_erf_with("${launch}" "${STRAIGHT_DIR}" "simulation.log" ${RUN_TIMEOUT}
        "max_step=${STEP_END}" "erf.check_int=-1" "erf.plot_int_1=${STEP_END}"
        ${plot2d_end})
# to the checkpoint step, writing it there
run_erf_with("${launch_chk}" "${RESTART_DIR}" "checkpoint.log" ${RUN_TIMEOUT}
        "max_step=${STEP_CHK}" "erf.check_int=${STEP_CHK}" "erf.plot_int_1=-1"
        ${plot2d_off})
if(NOT EXISTS "${RESTART_DIR}/${CHKFILE}/Header")
    message(FATAL_ERROR "RunRestartParity.cmake: no ${CHKFILE} written by the checkpoint run")
endif()
# from the checkpoint to the end
run_erf_with("${launch}" "${RESTART_DIR}" "restart.log" ${RUN_TIMEOUT}
        "erf.restart=${CHKFILE}" "max_step=${STEP_END}" "erf.check_int=-1" "erf.plot_int_1=${STEP_END}"
        ${plot2d_end} ${restart_options})

#
# A regrid test that silently stops regridding still passes, because it then compares an
# ordinary restart against the straight run and those agree trivially. Require the evidence
# in the log, so the test keeps testing what its name says. ERF prints this line from
# ReadCheckpointFile whenever level 0's grids differ from the checkpoint's, regardless of
# erf.v.
#
if(REQUIRE_LEVEL0_REMAKE)
    file(READ "${RESTART_DIR}/restart.log" _restart_output)
    if(NOT _restart_output MATCHES "reading level 0 onto new grids")
        message(FATAL_ERROR
            "RunRestartParity.cmake: the restart leg was expected to read level 0 onto "
            "new grids, but its log does not say it did, so this test is no longer "
            "exercising a regrid on restart (see restart.log)")
    endif()
endif()

foreach(dir "${STRAIGHT_DIR}" "${RESTART_DIR}")
    if(NOT EXISTS "${dir}/${PLTFILE}/Header")
        message(FATAL_ERROR "RunRestartParity.cmake: no ${PLTFILE} in ${dir}")
    endif()
endforeach()

execute_process(
    COMMAND ${launch_one} ${FCOMPARE} --abort_if_not_all_found ${fcompare_grids}
            --rel_tol ${RTOL} --abs_tol ${ATOL}
            ${STRAIGHT_DIR}/${PLTFILE} ${RESTART_DIR}/${PLTFILE}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/parity.log"
    ERROR_FILE "${WORKING_DIRECTORY}/parity.log"
    RESULT_VARIABLE parity_result)
if(NOT parity_result EQUAL 0)
    message(FATAL_ERROR "RunRestartParity.cmake: the restarted run's ${PLTFILE} differs from the straight run's: ${parity_result} (see parity.log)")
endif()
message(STATUS "RunRestartParity: restart from ${CHKFILE} reproduces ${PLTFILE}")

if(NOT "${PLT2DFILE}" STREQUAL "")
    foreach(dir "${STRAIGHT_DIR}" "${RESTART_DIR}")
        if(NOT EXISTS "${dir}/${PLT2DFILE}/Header")
            message(FATAL_ERROR
                "RunRestartParity.cmake: no ${PLT2DFILE} in ${dir}; the deck must select 2D "
                "output with erf.plot2d_vars_1 for the 2D comparison to mean anything")
        endif()
    endforeach()
    execute_process(
        COMMAND ${launch_one} ${FCOMPARE} --abort_if_not_all_found ${fcompare_grids}
                --rel_tol ${RTOL} --abs_tol ${ATOL}
                ${STRAIGHT_DIR}/${PLT2DFILE} ${RESTART_DIR}/${PLT2DFILE}
        WORKING_DIRECTORY "${WORKING_DIRECTORY}"
        OUTPUT_FILE "${WORKING_DIRECTORY}/parity2d.log"
        ERROR_FILE "${WORKING_DIRECTORY}/parity2d.log"
        RESULT_VARIABLE parity2d_result)
    if(NOT parity2d_result EQUAL 0)
        message(FATAL_ERROR
            "RunRestartParity.cmake: the restarted run's ${PLT2DFILE} differs from the straight "
            "run's: ${parity2d_result} (see parity2d.log). A surface field that is written to the "
            "checkpoint but not read back looks exactly like this.")
    endif()
    message(STATUS "RunRestartParity: restart also reproduces ${PLT2DFILE}")
endif()

# Optional: a time series that the run appends to, such as a station file written by
# erf.station_names, must come out the same whether it was written in one run or in two.
# The restarted run marks the restart with a comment line the straight run does not have,
# so the comparison is of the data lines only.
if(NOT "${DATALOG}" STREQUAL "")
    function(strip_comments in_file out_file out_count)
        file(STRINGS "${in_file}" _lines)
        set(_kept "")
        foreach(_line IN LISTS _lines)
            if(NOT _line MATCHES "^#")
                list(APPEND _kept "${_line}")
            endif()
        endforeach()
        list(LENGTH _kept _n)
        string(JOIN "\n" _text ${_kept})
        file(WRITE "${out_file}" "${_text}\n")
        set(${out_count} ${_n} PARENT_SCOPE)
    endfunction()

    foreach(dir "${STRAIGHT_DIR}" "${RESTART_DIR}")
        if(NOT EXISTS "${dir}/${DATALOG}")
            message(FATAL_ERROR "RunRestartParity.cmake: no time series ${dir}/${DATALOG}")
        endif()
    endforeach()

    strip_comments("${STRAIGHT_DIR}/${DATALOG}" "${WORKING_DIRECTORY}/datalog_straight.txt" straight_rows)
    strip_comments("${RESTART_DIR}/${DATALOG}"  "${WORKING_DIRECTORY}/datalog_restart.txt"  restart_rows)
    if(straight_rows LESS 2)
        message(FATAL_ERROR "RunRestartParity.cmake: ${DATALOG} has ${straight_rows} data rows; the comparison would be trivial")
    endif()

    if("${DATALOG_SIGDIGITS}" STREQUAL "")
        set(DATALOG_SIGDIGITS 6)
    endif()
    include("${CMAKE_CURRENT_LIST_DIR}/CompareDataLogs.cmake")
    erf_compare_data_logs("${WORKING_DIRECTORY}/datalog_straight.txt"
                          "${WORKING_DIRECTORY}/datalog_restart.txt"
                          ${DATALOG_SIGDIGITS} 2 logs_agree log_message)
    if(NOT logs_agree)
        message(FATAL_ERROR "RunRestartParity.cmake: ${DATALOG} differs between the straight run "
                            "and the restarted run: ${log_message}")
    endif()
    message(STATUS "RunRestartParity: ${DATALOG} agrees (${straight_rows} rows)")
endif()
