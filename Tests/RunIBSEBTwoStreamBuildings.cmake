# Runs Tests/test_files/IBSEB_TwoStreamBuildings/IBSEB_TwoStreamBuildings.i, where the
# two-stream columns and their ground balance see the buildings
# (erf.ibseb.two_stream_buildings), then check_ibseb_two_stream_buildings.py on the 2D and
# 3D plotfiles. The legs:
#   clear       the 40 m cube under a transparent sky, two steps, a checkpoint after the first;
#   restart     the clear leg restarted from that checkpoint to step 2;
#   absorbing   an absorbing sky with no scattering, one step;
#   scattering  the absorbing sky with scattering, so the shaded ground has diffuse light;
#   off         the absorbing sky with erf.ibseb.two_stream_buildings = false;
#   night       the cube and the tower under the calendar sun below the horizon (the
#               night deck of Tests/test_files/IBSEB_TwoStreamProvider), traps on;
#   two_level   the two-level deck of Tests/test_files/IBSEB_TwoStreamProvider (the tower
#               outside the refined level), one step;
#   fixed_sw    the clear leg with the ground balance on its fixed shortwave
#               (seb_use_radiation_fluxes off): footprints held, no shadow cast, a note;
#   nested      the deck of Tests/test_files/IBSEB_RefinedLevels under the two-stream
#               radiation with level 1 a patch below the domain top, whose heating rates come
#               from level 0 by interpolation, one step.
# The building files and soundings are those of Tests/test_files/IBSEB_TwoStreamProvider
# and Tests/test_files/IBSEB_RefinedLevels.
# -DX= defines X as empty, so test for a value, not for DEFINED
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

foreach(arg NRANKS TEST_EXE INPUT PROVIDER_FILES TWO_LEVEL_FILES WORKING_DIRECTORY PYTHON_EXE CHECKER PRECISION)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunIBSEBTwoStreamBuildings.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunIBSEBTwoStreamBuildings.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunIBSEBTwoStreamBuildings.cmake: ERF executable")

erf_mpi_launcher_command(launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunIBSEBTwoStreamBuildings.cmake")

# Stop with the run's own error lines (an abort message, an assertion), not only its code.
function(fail_run leg run_dir result)
    set(errors "")
    if(EXISTS "${run_dir}/simulation.log")
        file(STRINGS "${run_dir}/simulation.log" errors REGEX "[Aa]bort|[Ee]rror|[Aa]ssert|[Ee]xception|SIG")
    endif()
    list(JOIN errors "\n" errors)
    message(FATAL_ERROR "RunIBSEBTwoStreamBuildings.cmake: the ${leg} run failed (${result}):\n${errors}\n"
                        "(full log: ${run_dir}/simulation.log)")
endfunction()

# One run: <leg> directory (cleared first), deck, extra arguments.
function(run_leg leg deck)
    cmake_parse_arguments(R "" "" "ARGS" ${ARGN})
    set(run_dir "${WORKING_DIRECTORY}/${leg}")
    file(REMOVE_RECURSE "${run_dir}")
    file(MAKE_DIRECTORY "${run_dir}")
    file(COPY "${PROVIDER_FILES}/input_sounding" "${PROVIDER_FILES}/cube_40m_10m_32x32.txt"
              "${PROVIDER_FILES}/cube_and_tower_32x32.txt" "${TWO_LEVEL_FILES}/cube_and_tower_10m.txt"
         DESTINATION "${run_dir}")
    execute_process(
        COMMAND ${launch} ${TEST_EXE} ${deck} ${R_ARGS}
        WORKING_DIRECTORY "${run_dir}"
        OUTPUT_FILE "${run_dir}/simulation.log"
        ERROR_FILE "${run_dir}/simulation.log"
        TIMEOUT 600
        RESULT_VARIABLE run_result)
    if(NOT run_result EQUAL 0)
        fail_run(${leg} "${run_dir}" "${run_result}")
    endif()
endfunction()

# What the night and two-level decks need besides their own keys: the ground balance on
# the sweep's shortwave and the fields the checker reads, every step.
set(ground_args erf.radiation.seb_prognostic_enable=true erf.radiation.seb_use_radiation_fluxes=true
    erf.plot_int_1=1 erf.plot_vars_1=theta terrain_IB_mask qsrc_sw qsrc_lw
    erf.plot2d_file_1=plt2d erf.plot2d_int_1=1 erf.plot2d_vars_1=seb_sw_net ibseb_ground_shadow seb_t_sfc)

run_leg(clear ${INPUT} ARGS max_step=2 erf.check_int=1)
run_leg(restart ${INPUT} ARGS max_step=2 erf.restart=${WORKING_DIRECTORY}/clear/chk00001)
run_leg(absorbing ${INPUT} ARGS erf.radiation.tau_per_layer=0.02 erf.radiation.tau_lw_per_layer=0.3)
run_leg(scattering ${INPUT} ARGS erf.radiation.tau_per_layer=0.02 erf.radiation.tau_lw_per_layer=0.3
        erf.radiation.single_scattering_albedo=0.5)
run_leg(fixed_sw ${INPUT} ARGS max_step=2 erf.radiation.seb_use_radiation_fluxes=false
        erf.radiation.seb_sw_flux_default=300.0)
run_leg(off ${INPUT} ARGS erf.radiation.tau_per_layer=0.02 erf.radiation.tau_lw_per_layer=0.3
        erf.ibseb.two_stream_buildings=false)
run_leg(night "${PROVIDER_FILES}/IBSEB_TwoStreamProviderNight.i" ARGS max_step=1 ${ground_args})
run_leg(two_level "${PROVIDER_FILES}/IBSEB_TwoStreamProviderTwoLevel.i" ARGS max_step=1 ${ground_args})
run_leg(nested "${TWO_LEVEL_FILES}/IBSEB_RefinedLevels.i" ARGS max_step=1 erf.radiation_model=TwoStream
        erf.fixed_solar_zenith_angle=0.3420201433256688 erf.fixed_total_solar_irradiance=1000.0 erf.rad_t_sfc=300.0
        erf.radiation.tau_per_layer=0.02 erf.radiation.tau_lw_per_layer=0.3 amr.refine_whole_domain_dir=-1
        erf.city.in_box_lo=100.0 180.0 0.0 erf.city.in_box_hi=220.0 300.0 100.0
        erf.plot_int_1=1 erf.plot_vars_1=theta terrain_IB_mask qsrc_sw qsrc_lw)

# The deck's sun (cos z, zenith and azimuth of the faces' fixed sun), the top-of-atmosphere
# irradiance, the ground's albedo, the absorbing legs' optical depth per layer and their
# number of layers, the two-level deck's sun, and the precision of ERF's Real.
string(TOLOWER "${PRECISION}" precision)
execute_process(
    COMMAND "${PYTHON_EXE}" "${CHECKER}" "${WORKING_DIRECTORY}" --cosz 0.766044443118978 --zenith-deg 40
            --azimuth-deg 200 --toa 1000 --albedo 0.2 --tau 0.02 --nz 16 --two-level-zenith-deg 70
            --two-level-azimuth-deg 90 --precision ${precision}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/checker.log"
    ERROR_FILE "${WORKING_DIRECTORY}/checker.log"
    RESULT_VARIABLE check_result)
file(READ "${WORKING_DIRECTORY}/checker.log" check_output)
message("${check_output}")
if(NOT check_result EQUAL 0)
    message(FATAL_ERROR "RunIBSEBTwoStreamBuildings.cmake: the checker failed: ${check_result}")
endif()
