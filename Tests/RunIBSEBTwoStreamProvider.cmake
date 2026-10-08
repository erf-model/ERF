# Run Tests/test_files/IBSEB_TwoStreamProvider (the faces on erf.ibseb.radiation =
# two_stream) for two steps as given, under a transparent sky, and restart it from step 1;
# with the prescribed provider set to that sky; under an absorbing sky; and at night with a
# tower beside the cube (IBSEB_TwoStreamProviderNight.i). Then run
# Tests/test_files/IBSEB_RefinedLevels (a cube on level 1, a taller tower outside it) on the
# two-stream columns under a transparent and an absorbing sky. check_ibseb_two_stream_provider.py
# checks the face dumps and reports: the run-time path of the canopy forcing
# (TwoStreamRadiation::supply_canopy_forcing() and the sweep's write), the call of the
# balance after the step's radiation, the two_stream branches of
# IBFaceSet::compute_shortwave() and compute_longwave() with each face at its own height,
# and the initial report of a restart.
# -DX= defines X as empty, so test for a value, not for DEFINED
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

foreach(arg NRANKS TEST_EXE INPUT TWO_LEVEL_INPUT WORKING_DIRECTORY PYTHON_EXE CHECKER)
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

# One run: <leg> directory, deck, extra arguments. The directory is cleared unless KEEP is set.
function(run_leg leg deck)
    cmake_parse_arguments(R "KEEP" "" "ARGS" ${ARGN})
    set(run_dir "${WORKING_DIRECTORY}/${leg}")
    if(NOT R_KEEP)
        file(REMOVE_RECURSE "${run_dir}")
        file(MAKE_DIRECTORY "${run_dir}/faces")
    endif()
    file(COPY "${input_dir}/input_sounding" "${input_dir}/cube_40m_10m_32x32.txt"
              "${input_dir}/cube_and_tower_32x32.txt" DESTINATION "${run_dir}")
    execute_process(
        COMMAND ${launch} ${TEST_EXE} ${deck} ${R_ARGS}
        WORKING_DIRECTORY "${run_dir}"
        OUTPUT_FILE "${run_dir}/simulation.log"
        ERROR_FILE "${run_dir}/simulation.log"
        TIMEOUT 600
        RESULT_VARIABLE run_result)
    if(NOT run_result EQUAL 0)
        message(FATAL_ERROR "RunIBSEBTwoStreamProvider.cmake: the ${leg} run failed: ${run_result} (see ${run_dir}/simulation.log)")
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
run_leg(prescribed ${INPUT} ARGS max_step=2 erf.ibseb.radiation=prescribed erf.ibseb.sw_direct_normal=1000.0
        erf.ibseb.sw_diffuse=0.0 erf.ibseb.albedo_ground=0.2 erf.ibseb.lw_mode=fixed erf.ibseb.lw_down=0.0
        erf.ibseb.T_ground=300.0 erf.ibseb.emissivity_ground=0.95)
# The absorbing sky: shortwave optical depth 0.02 per layer, no scattering (the default
# single-scattering albedo is zero), and some longwave depth.
run_leg(absorbing ${INPUT} ARGS max_step=2 erf.radiation.tau_per_layer=0.02 erf.radiation.tau_lw_per_layer=0.3)
# Night under the calendar sun, with the tower.
get_filename_component(night_deck "${INPUT}" NAME_WE)
run_leg(night "${input_dir}/${night_deck}Night.i" ARGS max_step=2)

# The two-level deck on the two-stream columns, with the two-stream sun fixed where the
# deck's faces have theirs (70 degrees) and its own gray longwave.
get_filename_component(two_level_dir "${TWO_LEVEL_INPUT}" DIRECTORY)
set(two_level_args erf.radiation_model=TwoStream erf.fixed_solar_zenith_angle=0.3420201433256688
    erf.fixed_total_solar_irradiance=1000.0 erf.rad_t_sfc=300.0 erf.ibseb.radiation=two_stream
    erf.radiation.tau_lw_per_layer=0.3)
foreach(leg two_level_clear two_level_absorbing)
    set(run_dir "${WORKING_DIRECTORY}/${leg}")
    file(REMOVE_RECURSE "${run_dir}")
    file(MAKE_DIRECTORY "${run_dir}/faces")
    file(COPY "${two_level_dir}/input_sounding" "${two_level_dir}/cube_and_tower_10m.txt" DESTINATION "${run_dir}")
    if(leg STREQUAL "two_level_clear")
        set(tau erf.radiation.tau_per_layer=0.0)
    else()
        set(tau erf.radiation.tau_per_layer=0.02)
    endif()
    execute_process(
        COMMAND ${launch} ${TEST_EXE} ${TWO_LEVEL_INPUT} max_step=1 ${two_level_args} ${tau}
        WORKING_DIRECTORY "${run_dir}"
        OUTPUT_FILE "${run_dir}/simulation.log"
        ERROR_FILE "${run_dir}/simulation.log"
        TIMEOUT 600
        RESULT_VARIABLE run_result)
    if(NOT run_result EQUAL 0)
        message(FATAL_ERROR "RunIBSEBTwoStreamProvider.cmake: the ${leg} run failed: ${run_result} (see ${run_dir}/simulation.log)")
    endif()
endforeach()

# The decks' grids (16 layers), their suns (cos z) and the absorbing legs' depth.
execute_process(
    COMMAND "${PYTHON_EXE}" "${CHECKER}" "${WORKING_DIRECTORY}" --nz 16 --cosz 0.766044443118978 --tau 0.02
            --two-level-cosz 0.3420201433256688
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/checker.log"
    ERROR_FILE "${WORKING_DIRECTORY}/checker.log"
    RESULT_VARIABLE check_result)
file(READ "${WORKING_DIRECTORY}/checker.log" check_output)
message("${check_output}")
if(NOT check_result EQUAL 0)
    message(FATAL_ERROR "RunIBSEBTwoStreamProvider.cmake: the checker failed: ${check_result}")
endif()
