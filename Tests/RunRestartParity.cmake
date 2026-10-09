# Run a deck straight to STEP_END, run it again to STEP_CHK with a checkpoint there,
# restart from that checkpoint to STEP_END, and require the restarted run's plotfile
# at STEP_END to equal the straight run's with fcompare. Each leg uses the forwarded
# RUN_TIMEOUT, and the enclosing CTest timeout is sized separately by the caller.
# COMMON_OPTIONS goes to all three legs; RESTART_OPTIONS goes to the restart leg only, and
# CHK_OPTIONS to the checkpoint leg only. CHK_LEG_END runs the checkpoint leg past STEP_CHK,
# so the restart replays steps that a run stopped after its last checkpoint has already written.
# BNDRY_PLANES_DIR compares the boundary-plane series (erf.bndry_output_planes_file) of the two
# runs from step BNDRY_PLANES_FIRST_STEP on; BNDRY_PLANES_STEPS lists the steps the straight run
# must have written from there; BNDRY_PLANES_STALE seeds both run directories with the time.dat
# of an earlier run first; BNDRY_PLANES_CUT_NEWLINE removes the newline that ends the checkpoint
# leg's time.dat, as a run stopped while writing its last row leaves it; BNDRY_PLANES_READ_INPUT
# names a deck, beside INPUT, that reads the restarted run's series back with erf.input_bndry_planes.
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
separate_arguments(chk_options      UNIX_COMMAND "${CHK_OPTIONS}")

# The checkpoint leg stops at STEP_CHK unless CHK_LEG_END says to run it further
set(_chk_leg_end "${STEP_CHK}")
if(NOT "${CHK_LEG_END}" STREQUAL "")
    if(NOT "${CHK_LEG_END}" MATCHES "^[0-9]+$" OR CHK_LEG_END LESS STEP_CHK)
        message(FATAL_ERROR "RunRestartParity.cmake: CHK_LEG_END (${CHK_LEG_END}) is not a step at or after STEP_CHK (${STEP_CHK})")
    endif()
    set(_chk_leg_end "${CHK_LEG_END}")
endif()

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

# BNDRY_PLANES_STALE leaves the series of an earlier, longer run in both run directories, as a
# rerun in the same directory finds it: a fresh start has to begin a new series, not append to it.
if(BNDRY_PLANES_STALE)
    if("${BNDRY_PLANES_DIR}" STREQUAL "")
        message(FATAL_ERROR "RunRestartParity.cmake: BNDRY_PLANES_STALE needs BNDRY_PLANES_DIR")
    endif()
    foreach(dir "${STRAIGHT_DIR}" "${RESTART_DIR}")
        file(WRITE "${dir}/${BNDRY_PLANES_DIR}/time.dat" "0 0\n10 0.2\n20 0.4\n30 0.6\n")
    endforeach()
endif()

# straight to the end, no checkpoint
run_erf_with("${launch}" "${STRAIGHT_DIR}" "simulation.log" ${RUN_TIMEOUT}
        "max_step=${STEP_END}" "erf.check_int=-1" "erf.plot_int_1=${STEP_END}"
        ${plot2d_end})
# to the checkpoint step, writing it there (and on to CHK_LEG_END when that is given)
run_erf_with("${launch_chk}" "${RESTART_DIR}" "checkpoint.log" ${RUN_TIMEOUT}
        "max_step=${_chk_leg_end}" "erf.check_int=${STEP_CHK}" "erf.plot_int_1=-1"
        ${plot2d_off} ${chk_options})
if(NOT EXISTS "${RESTART_DIR}/${CHKFILE}/Header")
    message(FATAL_ERROR "RunRestartParity.cmake: no ${CHKFILE} written by the checkpoint run")
endif()
if(BNDRY_PLANES_CUT_NEWLINE)
    set(_bp_cut "${RESTART_DIR}/${BNDRY_PLANES_DIR}/time.dat")
    if(NOT EXISTS "${_bp_cut}")
        message(FATAL_ERROR "RunRestartParity.cmake: BNDRY_PLANES_CUT_NEWLINE: the checkpoint leg wrote no ${_bp_cut}")
    endif()
    file(READ "${_bp_cut}" _bp_cut_text)
    string(REGEX REPLACE "\n$" "" _bp_cut_text "${_bp_cut_text}")
    file(WRITE "${_bp_cut}" "${_bp_cut_text}")
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

# Optional: the boundary planes written for a downstream run (erf.output_bndry_planes). The
# series is a directory of planes, bndry_outputNNNNN, one per output step, and a time.dat that
# lists "step time" for each; erf.input_bndry_planes reads time.dat and requires both columns to
# increase. A restart must leave the series as the straight run writes it: the same time.dat
# rows, character for character (the restart's time is the straight run's to the last bit), and
# the same planes, byte for byte. AMReX renames a plane directory it writes again to
# <name>.old.<n>, which a restart from an earlier checkpoint does by design, so those are skipped.
if(NOT "${BNDRY_PLANES_DIR}" STREQUAL "")
    set(_bp_first 0)
    if(NOT "${BNDRY_PLANES_FIRST_STEP}" STREQUAL "")
        set(_bp_first "${BNDRY_PLANES_FIRST_STEP}")
    endif()

    # The step of a time.dat row or of a plane directory name, without leading zeros
    function(bp_step_of text out_var)
        if("${text}" MATCHES "bndry_output([0-9]+)$")
            set(_digits "${CMAKE_MATCH_1}")
        elseif("${text}" MATCHES "^[ \t]*([0-9]+)[ \t]")
            set(_digits "${CMAKE_MATCH_1}")
        else()
            message(FATAL_ERROR "RunRestartParity.cmake: cannot read a step from \"${text}\"")
        endif()
        string(REGEX REPLACE "^0+([0-9])" "\\1" _step "${_digits}")
        set(${out_var} "${_step}" PARENT_SCOPE)
    endfunction()

    foreach(dir "${STRAIGHT_DIR}" "${RESTART_DIR}")
        if(NOT EXISTS "${dir}/${BNDRY_PLANES_DIR}/time.dat")
            message(FATAL_ERROR "RunRestartParity.cmake: no ${BNDRY_PLANES_DIR}/time.dat in ${dir}")
        endif()
    endforeach()

    # time.dat: the straight run's rows from the first compared step on, against all of the restart's
    file(STRINGS "${STRAIGHT_DIR}/${BNDRY_PLANES_DIR}/time.dat" _bp_rows)
    set(_bp_kept "")
    set(_bp_steps "")
    foreach(_row IN LISTS _bp_rows)
        bp_step_of("${_row}" _step)
        if(NOT _step LESS _bp_first)
            list(APPEND _bp_kept "${_row}")
            list(APPEND _bp_steps "${_step}")
        endif()
    endforeach()
    list(LENGTH _bp_kept _bp_nrows)
    if(_bp_nrows LESS 2)
        message(FATAL_ERROR "RunRestartParity.cmake: the straight run's ${BNDRY_PLANES_DIR}/time.dat has "
                            "${_bp_nrows} rows from step ${_bp_first} on; the comparison would be trivial")
    endif()
    # BNDRY_PLANES_STEPS is space separated: a ";" would split it on its way through add_test
    string(JOIN " " _bp_steps ${_bp_steps})
    if(NOT "${BNDRY_PLANES_STEPS}" STREQUAL "" AND NOT "${_bp_steps}" STREQUAL "${BNDRY_PLANES_STEPS}")
        message(FATAL_ERROR "RunRestartParity.cmake: the straight run wrote planes at steps ${_bp_steps} "
                            "from step ${_bp_first} on; expected ${BNDRY_PLANES_STEPS}")
    endif()
    string(JOIN "\n" _bp_text ${_bp_kept})
    file(STRINGS "${RESTART_DIR}/${BNDRY_PLANES_DIR}/time.dat" _bp_restart_rows)
    string(JOIN "\n" _bp_restart_text ${_bp_restart_rows})

    # erf.input_bndry_planes stops on a time.dat whose steps or times do not increase
    foreach(dir "${STRAIGHT_DIR}" "${RESTART_DIR}")
        file(STRINGS "${dir}/${BNDRY_PLANES_DIR}/time.dat" _rows)
        set(_prev_step "")
        set(_prev_time "")
        foreach(_row IN LISTS _rows)
            bp_step_of("${_row}" _step)
            string(REGEX REPLACE "^[ \t]*[0-9]+[ \t]+([^ \t]+).*$" "\\1" _time "${_row}")
            if(NOT "${_prev_step}" STREQUAL "" AND (NOT _step GREATER _prev_step OR NOT _time GREATER _prev_time))
                message(FATAL_ERROR "RunRestartParity.cmake: ${dir}/${BNDRY_PLANES_DIR}/time.dat does not increase "
                                    "at the row \"${_row}\" after step ${_prev_step} at time ${_prev_time}; "
                                    "erf.input_bndry_planes cannot read it")
            endif()
            set(_prev_step "${_step}")
            set(_prev_time "${_time}")
        endforeach()
    endforeach()

    if(NOT "${_bp_text}" STREQUAL "${_bp_restart_text}")
        message(FATAL_ERROR "RunRestartParity.cmake: ${BNDRY_PLANES_DIR}/time.dat differs between the "
                            "straight run and the restarted run\n"
                            "straight (from step ${_bp_first}):\n${_bp_text}\nrestart:\n${_bp_restart_text}")
    endif()

    # The plane directories: the same set, and the same files with the same bytes in each
    function(bp_plane_dirs root first out_var)
        file(GLOB _names RELATIVE "${root}" "${root}/bndry_output*")
        set(_kept "")
        foreach(_name IN LISTS _names)
            if(IS_DIRECTORY "${root}/${_name}" AND NOT _name MATCHES "\\.old\\.")
                bp_step_of("${_name}" _step)
                if(NOT _step LESS first)
                    list(APPEND _kept "${_name}")
                endif()
            endif()
        endforeach()
        list(SORT _kept)
        set(${out_var} "${_kept}" PARENT_SCOPE)
    endfunction()
    bp_plane_dirs("${STRAIGHT_DIR}/${BNDRY_PLANES_DIR}" ${_bp_first} _bp_straight_dirs)
    bp_plane_dirs("${RESTART_DIR}/${BNDRY_PLANES_DIR}"  ${_bp_first} _bp_restart_dirs)
    if(NOT "${_bp_straight_dirs}" STREQUAL "${_bp_restart_dirs}")
        message(FATAL_ERROR "RunRestartParity.cmake: the plane directories differ from step ${_bp_first} on: "
                            "straight has ${_bp_straight_dirs}, restart has ${_bp_restart_dirs}")
    endif()
    set(_bp_nfiles 0)
    foreach(_plane IN LISTS _bp_straight_dirs)
        file(GLOB_RECURSE _s_files RELATIVE "${STRAIGHT_DIR}/${BNDRY_PLANES_DIR}/${_plane}"
             "${STRAIGHT_DIR}/${BNDRY_PLANES_DIR}/${_plane}/*")
        file(GLOB_RECURSE _r_files RELATIVE "${RESTART_DIR}/${BNDRY_PLANES_DIR}/${_plane}"
             "${RESTART_DIR}/${BNDRY_PLANES_DIR}/${_plane}/*")
        list(SORT _s_files)
        list(SORT _r_files)
        if(NOT "${_s_files}" STREQUAL "${_r_files}")
            message(FATAL_ERROR "RunRestartParity.cmake: ${_plane} holds different files: "
                                "straight has ${_s_files}, restart has ${_r_files}")
        endif()
        foreach(_f IN LISTS _s_files)
            execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files
                "${STRAIGHT_DIR}/${BNDRY_PLANES_DIR}/${_plane}/${_f}"
                "${RESTART_DIR}/${BNDRY_PLANES_DIR}/${_plane}/${_f}"
                RESULT_VARIABLE _same)
            if(NOT _same EQUAL 0)
                message(FATAL_ERROR "RunRestartParity.cmake: ${BNDRY_PLANES_DIR}/${_plane}/${_f} differs "
                                    "between the straight run and the restarted run")
            endif()
            math(EXPR _bp_nfiles "${_bp_nfiles} + 1")
        endforeach()
    endforeach()
    # The planes must change from the first compared step to the last, or a plane written at the
    # wrong step could not be told from the right one
    list(GET _bp_straight_dirs 0 _bp_plane_first)
    list(GET _bp_straight_dirs -1 _bp_plane_last)
    file(GLOB_RECURSE _bp_data RELATIVE "${STRAIGHT_DIR}/${BNDRY_PLANES_DIR}/${_bp_plane_first}"
         "${STRAIGHT_DIR}/${BNDRY_PLANES_DIR}/${_bp_plane_first}/*_D_*")
    foreach(_f IN LISTS _bp_data)
        execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files
            "${STRAIGHT_DIR}/${BNDRY_PLANES_DIR}/${_bp_plane_first}/${_f}"
            "${STRAIGHT_DIR}/${BNDRY_PLANES_DIR}/${_bp_plane_last}/${_f}"
            RESULT_VARIABLE _same)
        if(_same EQUAL 0)
            message(FATAL_ERROR "RunRestartParity.cmake: ${_f} is the same in ${_bp_plane_first} and "
                                "${_bp_plane_last}; the plane comparison would not tell the steps apart")
        endif()
    endforeach()

    list(LENGTH _bp_straight_dirs _bp_nplanes)

    # The reader is what the series is for: run it on the restarted run's planes
    if(NOT "${BNDRY_PLANES_READ_INPUT}" STREQUAL "")
        get_filename_component(_bp_deck_dir "${INPUT}" DIRECTORY)
        set(_bp_read_dir "${WORKING_DIRECTORY}/read")
        file(REMOVE_RECURSE "${_bp_read_dir}")
        file(MAKE_DIRECTORY "${_bp_read_dir}")
        execute_process(
            COMMAND ${launch} ${TEST_EXE} "${_bp_deck_dir}/${BNDRY_PLANES_READ_INPUT}"
                    "erf.bndry_file=${RESTART_DIR}/${BNDRY_PLANES_DIR}"
            WORKING_DIRECTORY "${_bp_read_dir}"
            OUTPUT_FILE "${_bp_read_dir}/read.log"
            ERROR_FILE "${_bp_read_dir}/read.log"
            TIMEOUT ${RUN_TIMEOUT}
            RESULT_VARIABLE _bp_read_result)
        if(NOT _bp_read_result EQUAL 0)
            message(FATAL_ERROR "RunRestartParity.cmake: ${BNDRY_PLANES_READ_INPUT} could not read the restarted "
                                "run's ${BNDRY_PLANES_DIR}: ${_bp_read_result} (see read/read.log)")
        endif()
        message(STATUS "RunRestartParity: ${BNDRY_PLANES_READ_INPUT} reads the restarted run's ${BNDRY_PLANES_DIR}")
    endif()
    message(STATUS "RunRestartParity: ${BNDRY_PLANES_DIR} agrees (${_bp_nrows} time.dat rows, "
                   "${_bp_nplanes} planes, ${_bp_nfiles} files)")
endif()
