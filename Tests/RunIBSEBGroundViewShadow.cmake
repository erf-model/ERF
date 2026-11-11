# Runs Tests/test_files/IBSEB_TwoStreamProvider/IBSEB_TwoStreamProvider.i (the 40 m cube)
# with and without erf.ibseb.ground_view_shadow, where the ground reflects the direct beam
# onto a face through the sunlit part of the face's ground view only, then
# check_ibseb_ground_view_shadow.py on the face dumps. The legs:
#   prescribed_off  the faces' own clear-sky radiation (a fixed beam and diffuse light), the input off;
#   prescribed_on   the same with the input on, two steps, a checkpoint after the first;
#   restart         prescribed_on restarted from that checkpoint to step 2;
#   serial          prescribed_on on one rank;
#   low_sun         prescribed_on with the sun 18 degrees above the horizon in the west-south-west;
#   two_stream_on   the deck as it is (the two-stream columns, a transparent sky), the input on;
#   absorbing_off, absorbing_on    an absorbing sky with no scattering;
#   scattering_off, scattering_on  that sky with scattering, over a ground whose albedo for
#                   diffuse light (0.35) differs from its albedo for the beam (0.2);
#   two_level_off, two_level_on    IBSEB_TwoStreamProviderTwoLevel.i (a cube on level 1, a
#                   taller tower outside it, the sun low in the east) under the absorbing sky.
# -DX= defines X as empty, so test for a value, not for DEFINED
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

foreach(arg NRANKS TEST_EXE INPUT TWO_LEVEL_FILES WORKING_DIRECTORY PYTHON_EXE CHECKER PRECISION)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunIBSEBGroundViewShadow.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunIBSEBGroundViewShadow.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunIBSEBGroundViewShadow.cmake: ERF executable")

erf_mpi_launcher_command(launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunIBSEBGroundViewShadow.cmake")
erf_mpi_launcher_command(launch_serial
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunIBSEBGroundViewShadow.cmake")

get_filename_component(input_dir "${INPUT}" DIRECTORY)

# Stop with the run's own error lines (an abort message, an assertion), not only its code.
function(fail_run leg run_dir result)
    set(errors "")
    if(EXISTS "${run_dir}/simulation.log")
        file(STRINGS "${run_dir}/simulation.log" errors REGEX "[Aa]bort|[Ee]rror|[Aa]ssert|[Ee]xception|SIG")
    endif()
    list(JOIN errors "\n" errors)
    message(FATAL_ERROR "RunIBSEBGroundViewShadow.cmake: the ${leg} run failed (${result}):\n${errors}\n"
                        "(full log: ${run_dir}/simulation.log)")
endfunction()

# One run: <leg> directory (cleared first), extra arguments; SERIAL runs it on one rank,
# DECK names another deck than INPUT.
function(run_leg leg)
    cmake_parse_arguments(R "SERIAL" "DECK" "ARGS" ${ARGN})
    set(run_dir "${WORKING_DIRECTORY}/${leg}")
    file(REMOVE_RECURSE "${run_dir}")
    file(MAKE_DIRECTORY "${run_dir}/faces")
    file(COPY "${input_dir}/input_sounding" "${input_dir}/cube_40m_10m_32x32.txt"
              "${TWO_LEVEL_FILES}/cube_and_tower_10m.txt" DESTINATION "${run_dir}")
    if(R_SERIAL)
        set(prefix ${launch_serial})
    else()
        set(prefix ${launch})
    endif()
    if("${R_DECK}" STREQUAL "")
        set(R_DECK "${INPUT}")
    endif()
    execute_process(
        COMMAND ${prefix} ${TEST_EXE} ${R_DECK} ${R_ARGS}
        WORKING_DIRECTORY "${run_dir}"
        OUTPUT_FILE "${run_dir}/simulation.log"
        ERROR_FILE "${run_dir}/simulation.log"
        TIMEOUT 600
        RESULT_VARIABLE run_result)
    if(NOT run_result EQUAL 0)
        fail_run(${leg} "${run_dir}" "${run_result}")
    endif()
endfunction()

# The faces' own clear-sky radiation: a direct-normal irradiance of 1000 W/m2, 80 W/m2 of
# diffuse light, and the ground of the deck's columns.
set(prescribed erf.ibseb.radiation=prescribed erf.ibseb.sw_direct_normal=1000.0 erf.ibseb.sw_diffuse=80.0
    erf.ibseb.albedo_ground=0.2 erf.ibseb.lw_mode=fixed erf.ibseb.lw_down=0.0 erf.ibseb.T_ground=300.0
    erf.ibseb.emissivity_ground=0.95)
set(on erf.ibseb.ground_view_shadow=true)
set(absorbing erf.radiation.tau_per_layer=0.02 erf.radiation.tau_lw_per_layer=0.3)
set(scattering ${absorbing} erf.radiation.single_scattering_albedo=0.5
    erf.radiation.surface_albedo_sw_diffuse=0.35)

run_leg(prescribed_off ARGS ${prescribed})
run_leg(prescribed_on ARGS ${prescribed} ${on} max_step=2 erf.check_int=1)
run_leg(restart ARGS ${prescribed} ${on} max_step=2 erf.restart=${WORKING_DIRECTORY}/prescribed_on/chk00001)
run_leg(serial SERIAL ARGS ${prescribed} ${on})
# The faces' sun and the two-stream sun must be one sun: cos(72 deg) for the columns too.
run_leg(low_sun ARGS ${prescribed} ${on} erf.ibseb.sun_zenith_deg=72.0 erf.ibseb.sun_azimuth_deg=255.0
        erf.fixed_solar_zenith_angle=0.3090169943749474)
run_leg(two_stream_on ARGS ${on})
run_leg(absorbing_off ARGS ${absorbing})
run_leg(absorbing_on ARGS ${absorbing} ${on})
run_leg(scattering_off ARGS ${scattering})
run_leg(scattering_on ARGS ${scattering} ${on})
run_leg(two_level_off DECK "${input_dir}/IBSEB_TwoStreamProviderTwoLevel.i" ARGS max_step=1 ${absorbing})
run_leg(two_level_on DECK "${input_dir}/IBSEB_TwoStreamProviderTwoLevel.i" ARGS max_step=1 ${absorbing} ${on})

# The deck's grid (32 x 32 columns of 10 m, 16 layers), the view rays per face (the
# defaults of erf.ibseb.view_n_az / view_n_el), the two suns, the faces' own irradiances,
# the two-stream irradiance at the top, the ground's albedos, the legs' optical depth per
# layer and single-scattering albedo, and the precision of ERF's Real.
string(TOLOWER "${PRECISION}" precision)
execute_process(
    COMMAND "${PYTHON_EXE}" "${CHECKER}" "${WORKING_DIRECTORY}" --nx 32 --ny 32 --dx 10 --nz 16 --n-az 16 --n-el 8
            --zenith-deg 40 --azimuth-deg 200 --low-zenith-deg 72 --low-azimuth-deg 255 --dni 1000 --diffuse 80
            --toa 1000 --albedo 0.2 --albedo-diffuse 0.35 --tau 0.02 --ssa 0.5 --precision ${precision}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/checker.log"
    ERROR_FILE "${WORKING_DIRECTORY}/checker.log"
    RESULT_VARIABLE check_result)
file(READ "${WORKING_DIRECTORY}/checker.log" check_output)
message("${check_output}")
if(NOT check_result EQUAL 0)
    message(FATAL_ERROR "RunIBSEBGroundViewShadow.cmake: the checker failed: ${check_result}")
endif()
