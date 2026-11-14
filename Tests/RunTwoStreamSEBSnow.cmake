# Run the TwoStream_SEBSnow deck (the force-restore snow layer, ERF_ForceRestoreSnow.H) and
# check the legs with test_files/TwoStream_SEBSnow/check_two_stream_seb_snow.py:
#   bare       the deck without its erf.radiation.seb_snow_* lines: the surface without snow;
#   no_snow    the snow layer on with no snow and no frozen precipitation (Kessler has none):
#              its last 2D plotfile equals the bare leg's (fcompare);
#   melt       the deck as it is, a 50 kg/m^2 snowpack under a warm sky: the surface stays at
#              273.16 K, the snow budget closes per cell, the sublimation is the covered
#              fraction of the evaporation and the absorbed shortwave is that of the albedo
#              with snow;
#   full_melt  1 kg/m^2 of snow with the snow albedo set to the ground's: it melts in the
#              first step, which leaves the surface colder than the bare leg's by
#              L_f melt / C_s and the soil wetter by the snow's water;
#   restart    the melt leg checkpointed half-way and restarted: its last 2D plotfile equals
#              the melt leg's (the snow and its totals are checkpointed);
#   two_level  the melt leg on two levels, the second made over part of the domain half-way
#              through: the melting point and the budget hold on every cell of the finest data.
# -DX= defines X as empty, so test for a value, not for DEFINED
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

foreach(arg NRANKS TEST_EXE INPUT WORKING_DIRECTORY FEXTRACT FCOMPARE PYTHON_EXE CHECKER
            STEPS DT)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunTwoStreamSEBSnow.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunTwoStreamSEBSnow.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunTwoStreamSEBSnow.cmake: ERF executable")
erf_resolve_executable(FEXTRACT "${FEXTRACT}" CONFIG "${CONFIG}"
    CONTEXT "RunTwoStreamSEBSnow.cmake: fextract")
erf_resolve_executable(FCOMPARE "${FCOMPARE}" CONFIG "${CONFIG}"
    CONTEXT "RunTwoStreamSEBSnow.cmake: fcompare")

erf_mpi_launcher_command(launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunTwoStreamSEBSnow.cmake")

get_filename_component(input_dir "${INPUT}" DIRECTORY)

# The values the checker's closed forms use. The first four are the deck's; the ground albedo
# and the moisture layer depth are given to every leg, so that the checker does not rest on
# their defaults.
set(SWE0 50.0)
set(W_CRN 10.0)
set(ALBEDO_SNOW 0.70)
set(C_S 2.0e4)
set(ALBEDO_GROUND 0.3)
set(D_S 0.1)
set(FULL_MELT_SWE0 1.0)
set(common
    erf.radiation.surface_albedo_sw=${ALBEDO_GROUND}
    erf.radiation.seb_moisture_layer_depth_m=${D_S})
set(bare_vars "seb_t_sfc seb_q_sfc seb_lh seb_sw_net seb_rain_total seb_runoff_total precip_rain_accum")

# The deck without its snow inputs (an input cannot be taken away on the command line).
file(STRINGS "${INPUT}" deck_lines)
set(bare_deck "${WORKING_DIRECTORY}/TwoStream_SEBSnow_bare.i")
file(WRITE "${bare_deck}" "")
foreach(line IN LISTS deck_lines)
    if(NOT line MATCHES "^erf\\.radiation\\.seb_snow_")
        file(APPEND "${bare_deck}" "${line}\n")
    endif()
endforeach()

# Run DECK in WORKING_DIRECTORY/<name>, writing <log>; FRESH empties the directory first.
# Extra inputs after steps.
function(run_in name log fresh deck steps)
    set(dir "${WORKING_DIRECTORY}/${name}")
    if(fresh)
        file(REMOVE_RECURSE "${dir}")
        file(MAKE_DIRECTORY "${dir}")
    endif()
    file(COPY "${input_dir}/input_sounding" DESTINATION "${dir}")
    execute_process(
        COMMAND ${launch} ${TEST_EXE} ${deck}
                max_step=${steps} erf.fixed_dt=${DT}
                ${common}
                ${ARGN}
        WORKING_DIRECTORY "${dir}"
        OUTPUT_FILE "${dir}/${log}"
        ERROR_FILE "${dir}/${log}"
        TIMEOUT 600
        RESULT_VARIABLE _result)
    if(NOT _result EQUAL 0)
        message(FATAL_ERROR "RunTwoStreamSEBSnow.cmake: the ${name} run failed: ${_result} (see ${dir}/${log})")
    endif()
endfunction()

function(same_last_plotfile name_a name_b log)
    string(LENGTH "${STEPS}" _len)
    math(EXPR _pad "5 - ${_len}")
    string(REPEAT "0" ${_pad} _zeros)
    set(last2d "plt2d${_zeros}${STEPS}")
    execute_process(
        COMMAND ${FCOMPARE} --abort_if_not_all_found
                "${WORKING_DIRECTORY}/${name_a}/${last2d}" "${WORKING_DIRECTORY}/${name_b}/${last2d}"
        WORKING_DIRECTORY "${WORKING_DIRECTORY}"
        OUTPUT_FILE "${log}"
        ERROR_FILE "${log}"
        RESULT_VARIABLE fcompare_result)
    if(NOT fcompare_result EQUAL 0)
        file(READ "${log}" fcompare_output)
        message(FATAL_ERROR "RunTwoStreamSEBSnow.cmake: ${last2d} of the ${name_b} leg differs "
                            "from the ${name_a} leg's:\n${fcompare_output}")
    endif()
endfunction()

run_in(bare simulation.log TRUE "${bare_deck}" ${STEPS} "erf.plot2d_vars_1=${bare_vars}")

# The snow layer on, with nothing in it: the same answers as without it.
run_in(no_snow simulation.log TRUE "${INPUT}" ${STEPS} "erf.plot2d_vars_1=${bare_vars}"
       erf.radiation.seb_snow_initial_swe_kg_m2=0.0)
same_last_plotfile(bare no_snow "${WORKING_DIRECTORY}/no_snow/fcompare.log")

run_in(melt simulation.log TRUE "${INPUT}" ${STEPS})

run_in(full_melt simulation.log TRUE "${INPUT}" ${STEPS}
       erf.radiation.seb_snow_initial_swe_kg_m2=${FULL_MELT_SWE0}
       erf.radiation.seb_snow_albedo=${ALBEDO_GROUND})

# The melt leg checkpointed half-way, then restarted from that checkpoint.
math(EXPR chk_step "${STEPS} / 2")
string(LENGTH "${chk_step}" _len)
math(EXPR _pad "5 - ${_len}")
string(REPEAT "0" ${_pad} _zeros)
set(chkfile "chk${_zeros}${chk_step}")
run_in(restart checkpoint.log TRUE "${INPUT}" ${chk_step} erf.check_int=${chk_step})
run_in(restart restart.log FALSE "${INPUT}" ${STEPS} erf.restart=${chkfile})
same_last_plotfile(melt restart "${WORKING_DIRECTORY}/restart/fcompare.log")

# The melt leg on two levels. Level 1 covers the middle half of x over the whole depth (the
# balance needs whole columns) and appears half-way through, with the regrid every step and
# no subcycling, so the finest data share the coarse step.
run_in(two_level simulation.log TRUE "${INPUT}" ${STEPS}
       "geometry.prob_extent=800 400 800" "amr.n_cell=8 4 16"
       amr.max_level=1 "amr.ref_ratio_vect=2 2 1" amr.n_error_buf=0 amr.blocking_factor=2
       amr.refine_whole_domain_dir=2 erf.refinement_indicators=patch erf.patch.max_level=1
       "erf.patch.in_box_lo=200.0 0.0" "erf.patch.in_box_hi=600.0 400.0"
       erf.patch.start_time=${chk_step}.5 erf.regrid_int=1 erf.dt_ref_ratio=1)

execute_process(
    COMMAND "${PYTHON_EXE}" "${CHECKER}" --fextract "${FEXTRACT}"
            --melt-dir "${WORKING_DIRECTORY}/melt" --bare-dir "${WORKING_DIRECTORY}/bare"
            --full-melt-dir "${WORKING_DIRECTORY}/full_melt"
            --two-level-dir "${WORKING_DIRECTORY}/two_level"
            --steps ${STEPS} --dt ${DT} --swe0 ${SWE0} --full-melt-swe0 ${FULL_MELT_SWE0}
            --w-crn ${W_CRN} --albedo-snow ${ALBEDO_SNOW} --albedo-ground ${ALBEDO_GROUND}
            --c-s ${C_S} --d-s ${D_S}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/checker.log"
    ERROR_FILE "${WORKING_DIRECTORY}/checker.log"
    RESULT_VARIABLE check_result)
file(READ "${WORKING_DIRECTORY}/checker.log" check_output)
message("${check_output}")
if(NOT check_result EQUAL 0)
    message(FATAL_ERROR "RunTwoStreamSEBSnow.cmake: the checker failed: ${check_result}")
endif()
