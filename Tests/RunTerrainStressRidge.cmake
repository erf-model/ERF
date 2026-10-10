# Driver for the Terrain_Stress_Ridge tests (Tests/test_files/Terrain_Stress_Ridge): stable flow
# over a steep periodic ridge on a 3 km grid with Smagorinsky2D + MRF and no numerical diffusion,
# where K_h * h^2 is far larger than K_v on the slopes (erf-model/ERF#4214).
#
# MODE = survive
#   Runs the deck with OPTIONS and requires it to finish, write PLTFILE and keep the vertical
#   velocity finite and inside (-WMAX, WMAX).  If CONTROL_OPTIONS is given, it then runs the deck
#   with those instead, and that run must start and then fail before it writes PLTFILE: the
#   control shows that the test sits where the setting under test matters.  A control that
#   survives means the test no longer tests anything, and fails the test.
#
# MODE = agree
#   Runs the deck with ON_OPTIONS and with OFF_OPTIONS (both must finish) and compares the two
#   plotfiles.  They must agree within RTOL (relative) or ATOL (absolute), which bounds a split
#   that changes the answer by more than its time discretization, and must differ at zero
#   tolerance, which shows the option is active (a misspelled option would pass agreement).
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

# -DX= defines X as empty, so test for a value, not for DEFINED
foreach(arg MODE NRANKS TEST_EXE INPUT WORKING_DIRECTORY FCOMPARE PLTFILE)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunTerrainStressRidge.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunTerrainStressRidge.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunTerrainStressRidge.cmake: ERF executable")
erf_resolve_executable(FCOMPARE "${FCOMPARE}" CONFIG "${CONFIG}"
    CONTEXT "RunTerrainStressRidge.cmake: fcompare")
if("${FEXTREMA}" STREQUAL "")
    string(REPLACE "amrex_fcompare" "amrex_fextrema" FEXTREMA "${FCOMPARE}")
endif()
erf_resolve_executable(FEXTREMA "${FEXTREMA}" CONFIG "${CONFIG}"
    CONTEXT "RunTerrainStressRidge.cmake: fextrema")

separate_arguments(options UNIX_COMMAND "${OPTIONS}")
separate_arguments(control_options UNIX_COMMAND "${CONTROL_OPTIONS}")
if("${RUN_TIMEOUT}" STREQUAL "")
    set(RUN_TIMEOUT 1200)
endif()

erf_mpi_launcher_command(launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunTerrainStressRidge.cmake")
erf_mpi_launcher_command(launch_one
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunTerrainStressRidge.cmake")

# Each leg runs in its own directory with a copy of the deck's input files (the sounding)
get_filename_component(input_dir "${INPUT}" DIRECTORY)
file(GLOB input_files LIST_DIRECTORIES false "${input_dir}/*")
function(make_leg_dir dir)
    file(REMOVE_RECURSE "${dir}")
    file(MAKE_DIRECTORY "${dir}")
    file(COPY ${input_files} DESTINATION "${dir}")
endfunction()

# run_leg(<dir> [options...]): run the deck in <dir>, fail if it does not finish
function(run_leg dir)
    make_leg_dir("${dir}")
    execute_process(
        COMMAND ${launch} ${TEST_EXE} ${INPUT} ${ARGN}
        WORKING_DIRECTORY "${dir}"
        OUTPUT_FILE "${dir}/simulation.log"
        ERROR_FILE "${dir}/simulation.log"
        TIMEOUT ${RUN_TIMEOUT}
        RESULT_VARIABLE _result)
    if(NOT _result EQUAL 0)
        message(FATAL_ERROR "RunTerrainStressRidge.cmake: the run in ${dir} failed: ${_result} (see simulation.log)")
    endif()
    if(NOT EXISTS "${dir}/${PLTFILE}/Header")
        message(FATAL_ERROR "RunTerrainStressRidge.cmake: no ${PLTFILE} in ${dir}")
    endif()
endfunction()

# check_w(<dir>): w in PLTFILE is finite and inside (-WMAX, WMAX)
function(check_w dir)
    execute_process(COMMAND ${launch_one} ${FEXTREMA} -v z_velocity ${dir}/${PLTFILE}
        OUTPUT_VARIABLE fx RESULT_VARIABLE fx_result)
    if(NOT fx_result EQUAL 0 OR NOT fx MATCHES "z_velocity[ \t]+([^ \t\n]+)[ \t]+([^ \t\n]+)")
        message(FATAL_ERROR "RunTerrainStressRidge.cmake: fextrema failed on ${dir}/${PLTFILE}: ${fx}")
    endif()
    set(wmin "${CMAKE_MATCH_1}")
    set(wmax "${CMAKE_MATCH_2}")
    # NaN compares false both ways, so require each bound inside (-WMAX, WMAX) positively
    if(NOT (wmax LESS WMAX AND wmin GREATER -${WMAX}))
        message(FATAL_ERROR "RunTerrainStressRidge.cmake: w in [${wmin}, ${wmax}] is not inside (-${WMAX}, ${WMAX})")
    endif()
    set(w_range "[${wmin}, ${wmax}]" PARENT_SCOPE)
endfunction()

if(MODE STREQUAL "survive")
    if("${WMAX}" STREQUAL "")
        message(FATAL_ERROR "RunTerrainStressRidge.cmake: WMAX must be given for MODE = survive")
    endif()
    set(dir "${WORKING_DIRECTORY}/survive")
    run_leg("${dir}" ${options})
    check_w("${dir}")

    set(control_note "no control")
    if(NOT "${CONTROL_OPTIONS}" STREQUAL "")
        set(cdir "${WORKING_DIRECTORY}/control")
        make_leg_dir("${cdir}")
        execute_process(
            COMMAND ${launch} ${TEST_EXE} ${INPUT} ${control_options}
            WORKING_DIRECTORY "${cdir}"
            OUTPUT_FILE "${cdir}/simulation.log"
            ERROR_FILE "${cdir}/simulation.log"
            TIMEOUT ${RUN_TIMEOUT}
            RESULT_VARIABLE control_result)
        file(STRINGS "${cdir}/simulation.log" control_steps REGEX "^Coarse STEP [0-9]+ ends")
        list(LENGTH control_steps n_control)
        if(control_result EQUAL 0)
            message(FATAL_ERROR "RunTerrainStressRidge.cmake: the control (${CONTROL_OPTIONS}) survived, so this time step no longer tests anything")
        endif()
        # A timeout reports a text result, not an exit code
        if(NOT control_result MATCHES "^-?[0-9]+$")
            message(FATAL_ERROR "RunTerrainStressRidge.cmake: the control did not fail, it stopped: ${control_result}")
        endif()
        # It must fail on the way, not after reaching the end (e.g. in the plotfile write)
        if(EXISTS "${cdir}/${PLTFILE}/Header")
            message(FATAL_ERROR "RunTerrainStressRidge.cmake: the control wrote ${PLTFILE}, so it failed after the run, not by instability")
        endif()
        if(n_control LESS 1)
            message(FATAL_ERROR "RunTerrainStressRidge.cmake: the control failed before its first step, not by instability (see control/simulation.log)")
        endif()
        set(control_note "the control failed after ${n_control} steps")
    endif()
    message(STATUS "RunTerrainStressRidge survive: w in ${w_range}; ${control_note}")

elseif(MODE STREQUAL "agree")
    foreach(arg ON_OPTIONS OFF_OPTIONS RTOL ATOL)
        if("${${arg}}" STREQUAL "")
            message(FATAL_ERROR "RunTerrainStressRidge.cmake: ${arg} must be given for MODE = agree")
        endif()
    endforeach()
    if("${ON_OPTIONS}" STREQUAL "${OFF_OPTIONS}")
        message(FATAL_ERROR "RunTerrainStressRidge.cmake: ON_OPTIONS and OFF_OPTIONS are the same, so the comparison would be trivial")
    endif()
    separate_arguments(on_options  UNIX_COMMAND "${ON_OPTIONS}")
    separate_arguments(off_options UNIX_COMMAND "${OFF_OPTIONS}")
    set(on_dir  "${WORKING_DIRECTORY}/option_on")
    set(off_dir "${WORKING_DIRECTORY}/option_off")
    run_leg("${on_dir}"  ${on_options})
    run_leg("${off_dir}" ${off_options})

    execute_process(
        COMMAND ${launch_one} ${FCOMPARE} --abort_if_not_all_found --rel_tol ${RTOL} --abs_tol ${ATOL}
                ${on_dir}/${PLTFILE} ${off_dir}/${PLTFILE}
        OUTPUT_FILE "${WORKING_DIRECTORY}/agree.log"
        ERROR_FILE "${WORKING_DIRECTORY}/agree.log"
        RESULT_VARIABLE agree_result)
    if(NOT agree_result EQUAL 0)
        message(FATAL_ERROR "RunTerrainStressRidge.cmake: the two runs differ by more than rel ${RTOL} / abs ${ATOL} (see agree.log)")
    endif()
    execute_process(
        COMMAND ${launch_one} ${FCOMPARE} --abort_if_not_all_found --rel_tol 0 --abs_tol 0
                ${on_dir}/${PLTFILE} ${off_dir}/${PLTFILE}
        OUTPUT_FILE "${WORKING_DIRECTORY}/identical.log"
        ERROR_FILE "${WORKING_DIRECTORY}/identical.log"
        RESULT_VARIABLE identical_result)
    if(identical_result EQUAL 0)
        message(FATAL_ERROR "RunTerrainStressRidge.cmake: the two runs are identical, so the option did nothing (see identical.log)")
    endif()
    message(STATUS "RunTerrainStressRidge agree: the runs differ, within rel ${RTOL} / abs ${ATOL}")

else()
    message(FATAL_ERROR "RunTerrainStressRidge.cmake: unknown MODE ${MODE}")
endif()
