# Driver for the Smag2D_Ridge tests (Tests/test_files/Smag2D_Ridge): stable flow over a steep
# periodic ridge on a 3 km grid with Smagorinsky2D + MRF and no numerical diffusion.
#
# MODE = steep
#   Runs the deck with OPTIONS (the WRF slope limiter) at a time step the limiter survives, and
#   requires the run to finish with a finite, bounded vertical velocity.  The start-up
#   slope-factor report must show alpha >= ALPHA_MIN.  It then runs the control, the same deck
#   with CONTROL_OPTIONS (the unlimited closure), which must start and then fail: that proves
#   the time step is beyond what the unlimited closure survives, so the test sits where the
#   limiter matters.  If a change to the terrain stress (ERF issue #4214) makes the control
#   survive, raise the time step to the new edge rather than drop the control.
#
# MODE = check
#   Runs the deck twice with the unlimited closure: once with erf.diffusive_dt_check on and a
#   small erf.diffusive_cfl, whose log must hold the diffusive Fourier-number warning, and once
#   with the check off, whose log must not.  The two plotfiles must be identical at zero
#   tolerance: the check is a diagnostic and must not change the answer.  The terrain
#   slope-factor report must appear in the first log and not in the second.
#
# MODE = limit
#   Runs the deck with an adaptive time step and erf.diffusive_dt_limit = true.  Every step
#   after the first must have DT <= the diffusive dt printed for it, and at least one step must
#   be set by that limit (DT equal to it), so the limit is shown to bind and to be obeyed.  The
#   first diffusive dt (step 2, from the diffusivities of step 1) must lie in
#   [REF_DIFFUSIVE_DT_LO, REF_DIFFUSIVE_DT_HI] (2 % about a reference run), so a wrong rate
#   (a dropped metric term, a lost density) is caught too.
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

# -DX= defines X as empty, so test for a value, not for DEFINED
foreach(arg MODE NRANKS TEST_EXE INPUT WORKING_DIRECTORY FCOMPARE PLTFILE)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunSmag2DRidge.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunSmag2DRidge.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunSmag2DRidge.cmake: ERF executable")
erf_resolve_executable(FCOMPARE "${FCOMPARE}" CONFIG "${CONFIG}"
    CONTEXT "RunSmag2DRidge.cmake: fcompare")
if("${FEXTREMA}" STREQUAL "")
    string(REPLACE "amrex_fcompare" "amrex_fextrema" FEXTREMA "${FCOMPARE}")
endif()
erf_resolve_executable(FEXTREMA "${FEXTREMA}" CONFIG "${CONFIG}"
    CONTEXT "RunSmag2DRidge.cmake: fextrema")

separate_arguments(options UNIX_COMMAND "${OPTIONS}")
separate_arguments(control_options UNIX_COMMAND "${CONTROL_OPTIONS}")
if((MODE STREQUAL "check" OR MODE STREQUAL "limit") AND "${DIFFUSIVE_CFL}" STREQUAL "")
    message(FATAL_ERROR "RunSmag2DRidge.cmake: DIFFUSIVE_CFL must be given for MODE = ${MODE}")
endif()
if(MODE STREQUAL "limit" AND ("${REF_DIFFUSIVE_DT_LO}" STREQUAL "" OR "${REF_DIFFUSIVE_DT_HI}" STREQUAL ""))
    message(FATAL_ERROR "RunSmag2DRidge.cmake: REF_DIFFUSIVE_DT_LO and REF_DIFFUSIVE_DT_HI must be given for MODE = limit")
endif()

erf_mpi_launcher_command(launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunSmag2DRidge.cmake")
erf_mpi_launcher_command(launch_one
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunSmag2DRidge.cmake")

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
        TIMEOUT 1200
        RESULT_VARIABLE _result)
    if(NOT _result EQUAL 0)
        message(FATAL_ERROR "RunSmag2DRidge.cmake: the run in ${dir} failed: ${_result} (see simulation.log)")
    endif()
    if(NOT EXISTS "${dir}/${PLTFILE}/Header")
        message(FATAL_ERROR "RunSmag2DRidge.cmake: no ${PLTFILE} in ${dir}")
    endif()
endfunction()

if(MODE STREQUAL "steep")
    foreach(arg ALPHA_MIN WMAX CONTROL_OPTIONS)
        if("${${arg}}" STREQUAL "")
            message(FATAL_ERROR "RunSmag2DRidge.cmake: ${arg} must be given for MODE = steep")
        endif()
    endforeach()
    set(dir "${WORKING_DIRECTORY}/steep")
    run_leg("${dir}" ${options})

    file(STRINGS "${dir}/simulation.log" alpha_line REGEX "Terrain slope factor alpha = h dx/dz at level 0: max")
    if(NOT alpha_line MATCHES "max ([0-9.eE+-]+),")
        message(FATAL_ERROR "RunSmag2DRidge.cmake: no terrain slope factor report in the log")
    endif()
    set(alpha_max "${CMAKE_MATCH_1}")
    if(alpha_max LESS ALPHA_MIN)
        message(FATAL_ERROR "RunSmag2DRidge.cmake: alpha_max = ${alpha_max} < ${ALPHA_MIN}: the ridge is too gentle to test the limiter")
    endif()

    execute_process(COMMAND ${launch_one} ${FEXTREMA} -v z_velocity ${dir}/${PLTFILE}
        OUTPUT_VARIABLE fx RESULT_VARIABLE fx_result)
    if(NOT fx_result EQUAL 0 OR NOT fx MATCHES "z_velocity[ \t]+([^ \t\n]+)[ \t]+([^ \t\n]+)")
        message(FATAL_ERROR "RunSmag2DRidge.cmake: fextrema failed on ${PLTFILE}: ${fx}")
    endif()
    set(wmin "${CMAKE_MATCH_1}")
    set(wmax "${CMAKE_MATCH_2}")
    # NaN compares false both ways, so require each bound inside (-WMAX, WMAX) positively
    if(NOT (wmax LESS WMAX AND wmin GREATER -${WMAX}))
        message(FATAL_ERROR "RunSmag2DRidge.cmake: w in [${wmin}, ${wmax}] is not inside (-${WMAX}, ${WMAX})")
    endif()

    # The control: the same step without the limiter must start and then fail
    set(cdir "${WORKING_DIRECTORY}/control")
    make_leg_dir("${cdir}")
    execute_process(
        COMMAND ${launch} ${TEST_EXE} ${INPUT} ${control_options}
        WORKING_DIRECTORY "${cdir}"
        OUTPUT_FILE "${cdir}/simulation.log"
        ERROR_FILE "${cdir}/simulation.log"
        TIMEOUT 1200
        RESULT_VARIABLE control_result)
    file(STRINGS "${cdir}/simulation.log" control_steps REGEX "^Coarse STEP [0-9]+ ends")
    list(LENGTH control_steps n_control)
    if(control_result EQUAL 0)
        message(FATAL_ERROR "RunSmag2DRidge.cmake: the control (${CONTROL_OPTIONS}) survived, so this time step no longer tests the limiter")
    endif()
    # A timeout reports a text result, not an exit code
    if(NOT control_result MATCHES "^-?[0-9]+$")
        message(FATAL_ERROR "RunSmag2DRidge.cmake: the control did not fail, it stopped: ${control_result}")
    endif()
    # It must fail on the way, not after reaching the end (e.g. in the plotfile write)
    if(EXISTS "${cdir}/${PLTFILE}/Header")
        message(FATAL_ERROR "RunSmag2DRidge.cmake: the control wrote ${PLTFILE}, so it failed after the run, not by instability")
    endif()
    if(n_control LESS 1)
        message(FATAL_ERROR "RunSmag2DRidge.cmake: the control failed before its first step, not by instability (see control/simulation.log)")
    endif()
    message(STATUS "RunSmag2DRidge steep: alpha_max = ${alpha_max}, w in [${wmin}, ${wmax}]; the control failed after ${n_control} steps")

elseif(MODE STREQUAL "check")
    set(on_dir  "${WORKING_DIRECTORY}/check_on")
    set(off_dir "${WORKING_DIRECTORY}/check_off")
    run_leg("${on_dir}"  ${options} erf.diffusive_dt_check=true erf.diffusive_cfl=${DIFFUSIVE_CFL})
    run_leg("${off_dir}" ${options} erf.diffusive_dt_check=false)

    set(warning "WARNING: explicit eddy diffusion at level 0 has Fourier number")
    file(STRINGS "${on_dir}/simulation.log"  on_warn  REGEX "${warning}")
    file(STRINGS "${off_dir}/simulation.log" off_warn REGEX "${warning}")
    if(on_warn STREQUAL "")
        message(FATAL_ERROR "RunSmag2DRidge.cmake: the diffusive check did not warn with erf.diffusive_cfl = ${DIFFUSIVE_CFL}")
    endif()
    if(NOT off_warn STREQUAL "")
        message(FATAL_ERROR "RunSmag2DRidge.cmake: the diffusive check warned with erf.diffusive_dt_check = false")
    endif()
    # The terrain slope-factor report belongs to the new options: present with the check on,
    # absent (nothing computed or printed) when the check, the limit and the Smagorinsky2D
    # limits are all off
    set(report "Terrain slope factor alpha = h dx/dz at level 0")
    file(STRINGS "${on_dir}/simulation.log"  on_report  REGEX "${report}")
    file(STRINGS "${off_dir}/simulation.log" off_report REGEX "${report}")
    if(on_report STREQUAL "")
        message(FATAL_ERROR "RunSmag2DRidge.cmake: no terrain slope-factor report with erf.diffusive_dt_check = true")
    endif()
    if(NOT off_report STREQUAL "")
        message(FATAL_ERROR "RunSmag2DRidge.cmake: the terrain slope-factor report printed with every new option off")
    endif()
    execute_process(
        COMMAND ${launch_one} ${FCOMPARE} --abort_if_not_all_found --rel_tol 0 --abs_tol 0
                ${on_dir}/${PLTFILE} ${off_dir}/${PLTFILE}
        OUTPUT_FILE "${WORKING_DIRECTORY}/parity.log"
        ERROR_FILE "${WORKING_DIRECTORY}/parity.log"
        RESULT_VARIABLE parity_result)
    if(NOT parity_result EQUAL 0)
        message(FATAL_ERROR "RunSmag2DRidge.cmake: the diffusive check changed ${PLTFILE} (see parity.log)")
    endif()
    message(STATUS "RunSmag2DRidge check: warned (${on_warn}) and left ${PLTFILE} unchanged")

elseif(MODE STREQUAL "limit")
    set(dir "${WORKING_DIRECTORY}/limit")
    run_leg("${dir}" ${options} erf.diffusive_dt_limit=true erf.diffusive_cfl=${DIFFUSIVE_CFL})

    file(STRINGS "${dir}/simulation.log" lines REGEX "^(Diffusive dt at level 0:|Coarse STEP [0-9]+ ends)")
    set(diff_dt "")
    set(first_diff_dt "")
    set(nsteps 0)
    set(nbinding 0)
    foreach(line IN LISTS lines)
        if(line MATCHES "^Diffusive dt at level 0:[ \t]+([^ \t]+)")
            set(diff_dt "${CMAKE_MATCH_1}")
            if(first_diff_dt STREQUAL "")
                set(first_diff_dt "${diff_dt}")
            endif()
        elseif(line MATCHES "^Coarse STEP ([0-9]+) ends.*DT = ([^ ]+)")
            set(step "${CMAKE_MATCH_1}")
            set(dt "${CMAKE_MATCH_2}")
            math(EXPR nsteps "${nsteps} + 1")
            # The diffusivities are zero before the first step, so the limit starts at step 2
            if(step GREATER 1)
                if(diff_dt STREQUAL "")
                    message(FATAL_ERROR "RunSmag2DRidge.cmake: no diffusive dt printed before step ${step}")
                endif()
                # Both are printed with 6 significant digits; a binding step prints the same text
                if(NOT dt LESS_EQUAL diff_dt)
                    message(FATAL_ERROR "RunSmag2DRidge.cmake: step ${step} took DT = ${dt} > diffusive dt ${diff_dt}")
                endif()
                if(dt STREQUAL diff_dt)
                    math(EXPR nbinding "${nbinding} + 1")
                endif()
            endif()
            set(diff_dt "")
        endif()
    endforeach()
    if(nsteps LESS 3)
        message(FATAL_ERROR "RunSmag2DRidge.cmake: only ${nsteps} steps found in the log")
    endif()
    if(nbinding EQUAL 0)
        message(FATAL_ERROR "RunSmag2DRidge.cmake: the diffusive limit never set the time step, so the test checks nothing")
    endif()
    # CMake's if() compares decimal strings as numbers; the bounds are 2 % about a
    # reference run, given by the registration
    if(first_diff_dt STREQUAL "")
        message(FATAL_ERROR "RunSmag2DRidge.cmake: no diffusive dt printed")
    endif()
    if(first_diff_dt LESS REF_DIFFUSIVE_DT_LO OR first_diff_dt GREATER REF_DIFFUSIVE_DT_HI)
        message(FATAL_ERROR "RunSmag2DRidge.cmake: first diffusive dt ${first_diff_dt} is outside [${REF_DIFFUSIVE_DT_LO}, ${REF_DIFFUSIVE_DT_HI}]")
    endif()
    message(STATUS "RunSmag2DRidge limit: ${nbinding} of ${nsteps} steps set by the diffusive limit; first diffusive dt ${first_diff_dt}")

else()
    message(FATAL_ERROR "RunSmag2DRidge.cmake: unknown MODE ${MODE}")
endif()
