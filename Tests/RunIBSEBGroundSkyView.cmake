# Runs Tests/test_files/IBSEB_TwoStreamBuildings/IBSEB_TwoStreamBuildings.i with
# erf.ibseb.ground_sky_view (the two-stream ground balance beside the buildings: the
# ground's view of the sky, and the walls' reflected light and longwave onto it), then
# check_ibseb_ground_sky_view.py on the 2D plotfiles and the face dumps. The legs:
#   clear_off, clear_on            the 40 m cube under a transparent sky, the input off
#                                  and on, the faces' skin temperature prognostic so that it
#                                  moves; two steps, a checkpoint after the first with it on;
#   restart                        clear_on restarted from that checkpoint to step 2;
#   serial                         clear_on on one rank;
#   scattering_off, scattering_on  an absorbing and scattering sky, the ground with a
#                                  diffuse albedo of its own, one step;
#   seam, seam_serial              a 10 m block in the periodic corner of a 64 x 64 domain in
#                                  four boxes of 32, the input on, on several ranks and on
#                                  one: with a sun floor of 20 degrees a rank's ray map
#                                  then reaches across the edge without spanning the domain;
#   two_level_off, two_level       the two-level deck of Tests/test_files/IBSEB_TwoStreamProvider
#                                  (the tower outside the refined level), the input off and on,
#                                  one step.
# 10 x 8 view rays per face, so no ray lands on a line between two columns. The building
# files and soundings are those of Tests/test_files/IBSEB_TwoStreamProvider and
# Tests/test_files/IBSEB_RefinedLevels.
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

foreach(arg NRANKS TEST_EXE INPUT PROVIDER_FILES TWO_LEVEL_FILES WORKING_DIRECTORY PYTHON_EXE CHECKER PRECISION)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunIBSEBGroundSkyView.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunIBSEBGroundSkyView.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunIBSEBGroundSkyView.cmake: ERF executable")

erf_mpi_launcher_command(launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunIBSEBGroundSkyView.cmake")

# Stop with the run's own error lines (an abort message, an assertion), not only its code.
function(fail_run leg run_dir result)
    set(errors "")
    if(EXISTS "${run_dir}/simulation.log")
        file(STRINGS "${run_dir}/simulation.log" errors REGEX "[Aa]bort|[Ee]rror|[Aa]ssert|[Ee]xception|SIG")
    endif()
    list(JOIN errors "\n" errors)
    message(FATAL_ERROR "RunIBSEBGroundSkyView.cmake: the ${leg} run failed (${result}):\n${errors}\n"
                        "(full log: ${run_dir}/simulation.log)")
endfunction()

# One run: <leg> directory (cleared first), deck, extra arguments.
function(run_leg leg deck)
    cmake_parse_arguments(R "" "" "ARGS" ${ARGN})
    set(run_dir "${WORKING_DIRECTORY}/${leg}")
    file(REMOVE_RECURSE "${run_dir}")
    file(MAKE_DIRECTORY "${run_dir}/faces")
    file(COPY "${PROVIDER_FILES}/input_sounding" "${PROVIDER_FILES}/cube_40m_10m_32x32.txt"
              "${PROVIDER_FILES}/cube_and_tower_32x32.txt" "${TWO_LEVEL_FILES}/cube_and_tower_10m.txt"
              "${WORKING_DIRECTORY}/block_10m_at_edge_64x64.txt"
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

# The face dumps and the fields the checker reads, every step.
set(common erf.ibseb.view_n_az=10 erf.ibseb.view_n_el=8 erf.ibseb.csv_int=1 erf.ibseb.dump_faces_file=faces
    erf.ibseb.dump_faces_tag_step=true
    erf.plot2d_vars_1=seb_sw_net seb_lw_net ibseb_ground_shadow ibseb_ground_wall_view seb_t_sfc)
set(on erf.ibseb.ground_sky_view=true)
set(scattering erf.radiation.tau_per_layer=0.02 erf.radiation.tau_lw_per_layer=0.3
    erf.radiation.single_scattering_albedo=0.5 erf.radiation.surface_albedo_sw_diffuse=0.35)

set(moving erf.ibseb.prognostic=true)
run_leg(clear_off ${INPUT} ARGS ${common} ${moving} max_step=2)
run_leg(clear_on ${INPUT} ARGS ${common} ${moving} ${on} max_step=2 erf.check_int=1)
run_leg(restart ${INPUT} ARGS ${common} ${moving} ${on} max_step=2 erf.restart=${WORKING_DIRECTORY}/clear_on/chk00001)
erf_mpi_launcher_command(launch_serial
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunIBSEBGroundSkyView.cmake")
set(launch_parallel ${launch})
set(launch ${launch_serial})
run_leg(serial ${INPUT} ARGS ${common} ${moving} ${on} max_step=2)
set(launch ${launch_parallel})
set(seam amr.n_cell=64 64 16 geometry.prob_extent=640 640 160 amr.max_grid_size=32
    erf.buildings_file_name=block_10m_at_edge_64x64.txt erf.ibseb.sun_min_elev_deg=20.0 erf.ibseb.debug=true)
run_leg(seam ${INPUT} ARGS ${common} ${scattering} ${on} ${seam})
set(launch ${launch_serial})
run_leg(seam_serial ${INPUT} ARGS ${common} ${scattering} ${on} ${seam})
set(launch ${launch_parallel})
run_leg(scattering_off ${INPUT} ARGS ${common} ${scattering})
run_leg(scattering_on ${INPUT} ARGS ${common} ${scattering} ${on})
set(two_level max_step=1 erf.radiation.tau_per_layer=0.02 erf.radiation.tau_lw_per_layer=0.3
    erf.radiation.seb_prognostic_enable=true erf.radiation.seb_use_radiation_fluxes=true
    erf.ibseb.view_n_az=10 erf.ibseb.view_n_el=8
    erf.plot2d_file_1=plt2d erf.plot2d_int_1=1
    erf.plot2d_vars_1=seb_sw_net seb_lw_net ibseb_ground_shadow ibseb_ground_wall_view seb_t_sfc)
run_leg(two_level_off "${PROVIDER_FILES}/IBSEB_TwoStreamProviderTwoLevel.i" ARGS ${two_level})
run_leg(two_level "${PROVIDER_FILES}/IBSEB_TwoStreamProviderTwoLevel.i" ARGS ${two_level} ${on})

# The deck's grid (32 x 32 columns of 10 m), the view rays, the ground's albedo for the
# beam and its emissivity, the scattering legs' diffuse albedo, and the precision of ERF's Real.
string(TOLOWER "${PRECISION}" precision)
execute_process(
    COMMAND "${PYTHON_EXE}" "${CHECKER}" "${WORKING_DIRECTORY}" --nx 32 --ny 32 --dx 10 --n-az 10 --n-el 8
            --albedo 0.2 --albedo-diffuse 0.35 --emissivity 0.95 --precision ${precision}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/checker.log"
    ERROR_FILE "${WORKING_DIRECTORY}/checker.log"
    RESULT_VARIABLE check_result)
file(READ "${WORKING_DIRECTORY}/checker.log" check_output)
message("${check_output}")
if(NOT check_result EQUAL 0)
    message(FATAL_ERROR "RunIBSEBGroundSkyView.cmake: the checker failed: ${check_result}")
endif()
