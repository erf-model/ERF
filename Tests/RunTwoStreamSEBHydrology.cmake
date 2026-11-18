# Run the TwoStream_SEBSurfaceLayerFluxes deck with the force-restore hydrology
# (ERF_ForceRestoreHydrology.H) and check the legs with check_two_stream_seb_hydrology.py:
#   deep  the deep reservoir and drainage, no rain: the column water drains exactly toward
#         field capacity (MN96 Eq. 4) and the deep temperature relaxes exactly toward the
#         skin (NP89 Eq. 6), step by step;
#   rain  rain from Kessler on a moist sounding, with the deep reservoir and drainage, over a
#         saturated soil (q_s = q_deep = porosity): the soil takes in what the scheme
#         accumulates at the ground, water runs off, and each column's water budget closes;
#   restart  the rain leg checkpointed half-way and restarted: its last 2D plotfile equals the
#         rain leg's (the deep state and the running totals are checkpointed);
#   two_level  the rain leg on two levels, the second made over part of the domain half-way
#         through (filled from the coarse level, regridded every step, averaged down): the
#         rain and budget checks hold on every column of the finest data;
#   two_level_subcycled  the same with level 1 taking two steps per coarse step: the rain the
#         soil took in matches the accumulation on every column of the finest data. (The rain carried
#         over a step a new level's balance skips is checked by the unit test
#         RainOfASkippedStepIsCarriedIntoTheNext: that skip needs erf.interp_atmos_from_coarse.)
# -DX= defines X as empty, so test for a value, not for DEFINED
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

foreach(arg NRANKS TEST_EXE INPUT WORKING_DIRECTORY FEXTRACT FCOMPARE PYTHON_EXE CHECKER
            DEEP_STEPS RAIN_STEPS DT)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunTwoStreamSEBHydrology.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunTwoStreamSEBHydrology.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunTwoStreamSEBHydrology.cmake: ERF executable")
erf_resolve_executable(FEXTRACT "${FEXTRACT}" CONFIG "${CONFIG}"
    CONTEXT "RunTwoStreamSEBHydrology.cmake: fextract")
erf_resolve_executable(FCOMPARE "${FCOMPARE}" CONFIG "${CONFIG}"
    CONTEXT "RunTwoStreamSEBHydrology.cmake: fcompare")

erf_mpi_launcher_command(launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunTwoStreamSEBHydrology.cmake")

get_filename_component(input_dir "${INPUT}" DIRECTORY)

# The soil of both legs, in the units of the force-restore inputs. The restore period is short
# so that the deep temperature moves measurably in ten one-second steps (dt / tau = 1e-2); the
# drainage uses the 1-day period of MN96 Table 1 whatever tau is.
set(TAU 100.0)
set(D2 1.0)
set(C3_M 0.5)
set(FC 0.30)
set(POROSITY 0.45)
set(Q_DEEP0 0.40)
set(soil
    erf.radiation.seb_restore_timescale_s=${TAU}
    erf.radiation.seb_deep_prognostic_enable=true
    erf.radiation.seb_deep_layer_depth_m=${D2}
    erf.radiation.seb_drainage_enable=true
    erf.radiation.seb_drainage_c3_m=${C3_M}
    erf.radiation.seb_soil_moisture_fc=${FC}
    erf.radiation.seb_soil_porosity=${POROSITY})
set(vars "seb_t_sfc seb_q_sfc seb_lh seb_t_deep seb_q_deep seb_rain_rate seb_rain_total seb_runoff_total seb_drainage_total seb_deficit_total precip_total_accum")

# Run the deck in WORKING_DIRECTORY/<name>, writing <log>; FRESH empties the directory first.
# Extra inputs after steps.
function(run_in name log fresh sounding steps)
    set(dir "${WORKING_DIRECTORY}/${name}")
    if(fresh)
        file(REMOVE_RECURSE "${dir}")
        file(MAKE_DIRECTORY "${dir}")
    endif()
    file(COPY "${sounding}" DESTINATION "${dir}")
    get_filename_component(sounding_name "${sounding}" NAME)
    execute_process(
        COMMAND ${launch} ${TEST_EXE} ${INPUT}
                max_step=${steps} erf.fixed_dt=${DT}
                erf.input_sounding_file=${sounding_name}
                "erf.plot2d_vars_1=${vars}"
                ${soil}
                ${ARGN}
        WORKING_DIRECTORY "${dir}"
        OUTPUT_FILE "${dir}/${log}"
        ERROR_FILE "${dir}/${log}"
        TIMEOUT 600
        RESULT_VARIABLE _result)
    if(NOT _result EQUAL 0)
        message(FATAL_ERROR "RunTwoStreamSEBHydrology.cmake: the ${name} run failed: ${_result} (see ${dir}/${log})")
    endif()
endfunction()

function(run_leg name sounding steps)
    run_in(${name} simulation.log TRUE "${sounding}" ${steps} ${ARGN})
endfunction()

run_leg(deep "${input_dir}/input_sounding" ${DEEP_STEPS}
        erf.radiation.seb_q_sfc_default=0.35
        erf.radiation.seb_q_deep_default=${Q_DEEP0})

# A sounding saturated below 400 m, so that Kessler makes cloud and rain at once (the surface
# layer's surf_moist matches its surface value, as the sounding reader requires).
set(rain_sounding "${WORKING_DIRECTORY}/input_sounding_rain")
file(WRITE "${rain_sounding}"
"950.0 301.5 24.0
  10.0 300.0 24.0 5.0 0.0
 400.0 300.0 24.0 5.0 0.0
 500.0 303.0 6.0 5.0 0.0
2000.0 307.5 4.0 5.0 0.0
")
set(rain_inputs
    erf.radiation.seb_precipitation_source=microphysics
    erf.most.surf_moist=0.024
    erf.radiation.seb_q_sfc_default=${POROSITY}
    erf.radiation.seb_q_deep_default=${POROSITY})
run_leg(rain "${rain_sounding}" ${RAIN_STEPS} ${rain_inputs})

# The rain leg checkpointed half-way, then restarted from that checkpoint.
math(EXPR chk_step "${RAIN_STEPS} / 2")
string(LENGTH "${chk_step}" _len)
math(EXPR _pad "5 - ${_len}")
string(REPEAT "0" ${_pad} _zeros)
set(chkfile "chk${_zeros}${chk_step}")
run_in(restart checkpoint.log TRUE "${rain_sounding}" ${chk_step} ${rain_inputs}
       erf.check_int=${chk_step})
run_in(restart restart.log FALSE "${rain_sounding}" ${RAIN_STEPS} ${rain_inputs}
       erf.restart=${chkfile})
string(LENGTH "${RAIN_STEPS}" _len)
math(EXPR _pad "5 - ${_len}")
string(REPEAT "0" ${_pad} _zeros)
set(last2d "plt2d${_zeros}${RAIN_STEPS}")
execute_process(
    COMMAND ${FCOMPARE} --abort_if_not_all_found
            "${WORKING_DIRECTORY}/rain/${last2d}" "${WORKING_DIRECTORY}/restart/${last2d}"
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/restart/fcompare.log"
    ERROR_FILE "${WORKING_DIRECTORY}/restart/fcompare.log"
    RESULT_VARIABLE fcompare_result)
if(NOT fcompare_result EQUAL 0)
    file(READ "${WORKING_DIRECTORY}/restart/fcompare.log" fcompare_output)
    message(FATAL_ERROR "RunTwoStreamSEBHydrology.cmake: the restarted ${last2d} differs from "
                        "the straight run's:\n${fcompare_output}")
endif()

# The rain leg on two levels. Level 1 covers the middle half of x over the whole depth (the
# balance needs whole columns) and appears half-way through, with the regrid every step and
# no subcycling, so the finest data share the coarse step.
math(EXPR level_start "${RAIN_STEPS} / 2")
run_leg(two_level "${rain_sounding}" ${RAIN_STEPS} ${rain_inputs}
        "geometry.prob_extent=800 400 800" "amr.n_cell=8 4 16"
        amr.max_level=1 "amr.ref_ratio_vect=2 2 1" amr.n_error_buf=0 amr.blocking_factor=2
        amr.refine_whole_domain_dir=2 erf.refinement_indicators=patch erf.patch.max_level=1
        "erf.patch.in_box_lo=200.0 0.0" "erf.patch.in_box_hi=600.0 400.0"
        erf.patch.start_time=${level_start}.5 erf.regrid_int=1 erf.dt_ref_ratio=1)

# The same with subcycling: level 1 takes two steps per coarse step, and the coarse data under
# it are averaged down from level 1 between coarse steps.
run_leg(two_level_subcycled "${rain_sounding}" ${RAIN_STEPS} ${rain_inputs}
        "geometry.prob_extent=800 400 800" "amr.n_cell=8 4 16"
        amr.max_level=1 "amr.ref_ratio_vect=2 2 1" amr.n_error_buf=0 amr.blocking_factor=2
        amr.refine_whole_domain_dir=2 erf.refinement_indicators=patch erf.patch.max_level=1
        "erf.patch.in_box_lo=200.0 0.0" "erf.patch.in_box_hi=600.0 400.0"
        erf.patch.start_time=${level_start}.5 erf.regrid_int=1 erf.dt_ref_ratio=2)

execute_process(
    COMMAND "${PYTHON_EXE}" "${CHECKER}" --fextract "${FEXTRACT}"
            --deep-dir "${WORKING_DIRECTORY}/deep" --deep-steps ${DEEP_STEPS}
            --rain-dir "${WORKING_DIRECTORY}/rain" --rain-steps ${RAIN_STEPS}
            --two-level-dir "${WORKING_DIRECTORY}/two_level"
            --subcycled-dir "${WORKING_DIRECTORY}/two_level_subcycled"
            --dt ${DT} --tau ${TAU} --d2 ${D2} --c3-m ${C3_M} --fc ${FC}
            --porosity ${POROSITY} --q-deep0 ${Q_DEEP0}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/checker.log"
    ERROR_FILE "${WORKING_DIRECTORY}/checker.log"
    RESULT_VARIABLE check_result)
file(READ "${WORKING_DIRECTORY}/checker.log" check_output)
message("${check_output}")
if(NOT check_result EQUAL 0)
    message(FATAL_ERROR "RunTwoStreamSEBHydrology.cmake: the checker failed: ${check_result}")
endif()
