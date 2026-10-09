# Runs the decks of Tests/test_files/IBSEB_TwoStreamProvider, where the building faces take
# their radiation from the two-stream columns (erf.ibseb.radiation = two_stream), then
# check_ibseb_two_stream_provider.py on the face dumps and reports. The legs:
#   two_stream   the 40 m cube under a transparent sky, two steps, restarted from step 1;
#   prescribed   the same cube on the faces' own clear-sky radiation, set to that sky;
#   absorbing    an absorbing sky with no scattering (an exact answer to compare with);
#   scattering   the absorbing sky with scattering, so the sky has diffuse light;
#   night        the sun below the horizon, a tower beside the cube, traps on;
#   two_level_*  IBSEB_TwoStreamProviderTwoLevel.i (level 1 refined in z too, with the
#                columns' longwave on the faces) under a transparent and an absorbing sky.
# -DX= defines X as empty, so test for a value, not for DEFINED
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

foreach(arg NRANKS TEST_EXE INPUT TWO_LEVEL_INPUT WORKING_DIRECTORY PYTHON_EXE CHECKER PRECISION)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunIBSEBTwoStreamProvider.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunIBSEBTwoStreamProvider.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunIBSEBTwoStreamProvider.cmake: ERF executable")

erf_mpi_launcher_command(launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunIBSEBTwoStreamProvider.cmake")

get_filename_component(input_dir "${INPUT}" DIRECTORY)
get_filename_component(two_level_files "${input_dir}/../IBSEB_RefinedLevels" ABSOLUTE)

# Stop with the run's own error lines (an abort message, an assertion), not only its code.
function(fail_run leg run_dir result)
    file(STRINGS "${run_dir}/simulation.log" errors REGEX "[Aa]bort|[Ee]rror|[Aa]ssert|[Ee]xception|SIG")
    list(JOIN errors "\n" errors)
    message(FATAL_ERROR "RunIBSEBTwoStreamProvider.cmake: the ${leg} run failed (${result}):\n${errors}\n"
                        "(full log: ${run_dir}/simulation.log)")
endfunction()

# One run: <leg> directory, deck, extra arguments. The directory is cleared unless KEEP is set.
function(run_leg leg deck)
    cmake_parse_arguments(R "KEEP" "" "ARGS" ${ARGN})
    set(run_dir "${WORKING_DIRECTORY}/${leg}")
    if(NOT R_KEEP)
        file(REMOVE_RECURSE "${run_dir}")
        file(MAKE_DIRECTORY "${run_dir}/faces")
    endif()
    file(COPY "${input_dir}/input_sounding" "${input_dir}/cube_40m_10m_32x32.txt"
              "${input_dir}/cube_and_tower_32x32.txt" "${two_level_files}/cube_and_tower_10m.txt"
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

# The transparent sky, two steps with a checkpoint after each, then a restart from step 1
# into the same directory (its report appends to the same CSV, its dumps to the same files).
run_leg(two_stream ${INPUT} ARGS max_step=2 erf.check_int=1)
file(RENAME "${WORKING_DIRECTORY}/two_stream/simulation.log" "${WORKING_DIRECTORY}/two_stream/simulation_first.log")
run_leg(two_stream ${INPUT} KEEP ARGS max_step=2 erf.check_int=-1 erf.restart=chk00001)
# The prescribed provider on the transparent sky of the deck: the two-stream irradiance as
# the direct-normal one, no diffuse light, no sky longwave, and the column's ground
# (erf.radiation.surface_albedo_sw / surface_emissivity_lw, erf.rad_t_sfc).
# Restarted from step 1 too, as a restart must not report that step again with either kind
# of radiation.
set(prescribed_args erf.ibseb.radiation=prescribed erf.ibseb.sw_direct_normal=1000.0
    erf.ibseb.sw_diffuse=0.0 erf.ibseb.albedo_ground=0.2 erf.ibseb.lw_mode=fixed erf.ibseb.lw_down=0.0
    erf.ibseb.T_ground=300.0 erf.ibseb.emissivity_ground=0.95)
run_leg(prescribed ${INPUT} ARGS max_step=2 erf.check_int=1 ${prescribed_args})
file(RENAME "${WORKING_DIRECTORY}/prescribed/simulation.log" "${WORKING_DIRECTORY}/prescribed/simulation_first.log")
run_leg(prescribed ${INPUT} KEEP ARGS max_step=2 erf.check_int=-1 erf.restart=chk00001 ${prescribed_args})
# The absorbing sky: shortwave optical depth 0.02 per layer, no scattering (the default
# single-scattering albedo is zero), and some longwave depth.
run_leg(absorbing ${INPUT} ARGS max_step=2 erf.radiation.tau_per_layer=0.02 erf.radiation.tau_lw_per_layer=0.3)
# The same sky, half of its extinction scattering: the beam is the same, and the sky now
# sends diffuse light.
run_leg(scattering ${INPUT} ARGS max_step=2 erf.radiation.tau_per_layer=0.02 erf.radiation.tau_lw_per_layer=0.3
        erf.radiation.single_scattering_albedo=0.5)
# Night under the calendar sun, with the tower.
get_filename_component(night_deck "${INPUT}" NAME_WE)
run_leg(night "${input_dir}/${night_deck}Night.i" ARGS max_step=2)

# Two levels, level 1 refined in z too, the columns' longwave on the faces: transparent
# and absorbing skies (no scattering).
run_leg(two_level_clear ${TWO_LEVEL_INPUT} ARGS max_step=1 erf.radiation.tau_per_layer=0.0
        erf.radiation.tau_lw_per_layer=0.3)
run_leg(two_level_absorbing ${TWO_LEVEL_INPUT} ARGS max_step=1 erf.radiation.tau_per_layer=0.02
        erf.radiation.tau_lw_per_layer=0.3)

# The decks' grids (16 layers on level 0), their suns (cos z), the absorbing legs' depth
# per layer, and the precision of ERF's Real (the tolerances of the comparisons).
string(TOLOWER "${PRECISION}" precision)
execute_process(
    COMMAND "${PYTHON_EXE}" "${CHECKER}" "${WORKING_DIRECTORY}" --nz 16 --cosz 0.766044443118978 --tau 0.02
            --two-level-cosz 0.3420201433256688 --precision ${precision}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/checker.log"
    ERROR_FILE "${WORKING_DIRECTORY}/checker.log"
    RESULT_VARIABLE check_result)
file(READ "${WORKING_DIRECTORY}/checker.log" check_output)
message("${check_output}")
if(NOT check_result EQUAL 0)
    message(FATAL_ERROR "RunIBSEBTwoStreamProvider.cmake: the checker failed: ${check_result}")
endif()
