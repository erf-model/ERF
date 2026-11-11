# Runs Tests/test_files/IBSEB_TwoStreamBuildings/IBSEB_TwoStreamBuildings.i with
# erf.ibseb.two_stream_roofs (the two-stream column under a building starts on the
# building's roof), then check_ibseb_two_stream_roofs.py on the plotfiles and the face
# dumps. The legs:
#   clear_off, clear_on            the 40 m cube under a transparent sky, the input off
#                                  and on, the faces' skin temperature prognostic; two steps;
#   scattering_off, scattering_on  an absorbing and scattering sky; two steps, a checkpoint
#                                  after the first with the input on;
#   restart                        scattering_on restarted from that checkpoint to step 2;
#   shadow_on                      scattering_on with erf.ibseb.ground_view_shadow too;
#   serial                         scattering_on on one rank;
#   two_level_off, two_level       the two-level deck of Tests/test_files/IBSEB_TwoStreamProvider
#                                  (the tower outside the refined level), one step, traps on.
# The building files and soundings are those of Tests/test_files/IBSEB_TwoStreamProvider and
# Tests/test_files/IBSEB_RefinedLevels.
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

foreach(arg NRANKS TEST_EXE INPUT PROVIDER_FILES TWO_LEVEL_FILES WORKING_DIRECTORY PYTHON_EXE CHECKER PRECISION)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunIBSEBTwoStreamRoofs.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunIBSEBTwoStreamRoofs.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunIBSEBTwoStreamRoofs.cmake: ERF executable")

erf_mpi_launcher_command(launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunIBSEBTwoStreamRoofs.cmake")

# Stop with the run's own error lines (an abort message, an assertion), not only its code.
function(fail_run leg run_dir result)
    set(errors "")
    if(EXISTS "${run_dir}/simulation.log")
        file(STRINGS "${run_dir}/simulation.log" errors REGEX "[Aa]bort|[Ee]rror|[Aa]ssert|[Ee]xception|SIG")
    endif()
    list(JOIN errors "\n" errors)
    message(FATAL_ERROR "RunIBSEBTwoStreamRoofs.cmake: the ${leg} run failed (${result}):\n${errors}\n"
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
set(common erf.ibseb.csv_int=1 erf.ibseb.dump_faces_file=faces erf.ibseb.dump_faces_tag_step=true
    erf.ibseb.prognostic=true max_step=2
    erf.plot_int_1=1 erf.plot_vars_1=theta terrain_IB_mask qsrc_sw qsrc_lw
    erf.plot2d_vars_1=seb_sw_net seb_lw_net ibseb_ground_shadow seb_t_sfc OLR)
set(on erf.ibseb.two_stream_roofs=true)
set(scattering erf.radiation.tau_per_layer=0.02 erf.radiation.tau_lw_per_layer=0.3
    erf.radiation.single_scattering_albedo=0.5)

run_leg(clear_off ${INPUT} ARGS ${common})
run_leg(clear_on ${INPUT} ARGS ${common} ${on})
run_leg(scattering_off ${INPUT} ARGS ${common} ${scattering})
run_leg(scattering_on ${INPUT} ARGS ${common} ${scattering} ${on} erf.check_int=1)
run_leg(restart ${INPUT} ARGS ${common} ${scattering} ${on} erf.restart=${WORKING_DIRECTORY}/scattering_on/chk00001)
run_leg(shadow_on ${INPUT} ARGS ${common} ${scattering} ${on} erf.ibseb.ground_view_shadow=true)
erf_mpi_launcher_command(launch_serial
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunIBSEBTwoStreamRoofs.cmake")
set(launch_parallel ${launch})
set(launch ${launch_serial})
run_leg(serial ${INPUT} ARGS ${common} ${scattering} ${on})
set(launch ${launch_parallel})
# (A scattering sky: without scattering no light the roof reflects comes back down.)
set(two_level max_step=1 erf.radiation.tau_per_layer=0.02 erf.radiation.tau_lw_per_layer=0.3
    erf.radiation.single_scattering_albedo=0.5
    erf.radiation.seb_prognostic_enable=true erf.radiation.seb_use_radiation_fluxes=true)
run_leg(two_level_off "${PROVIDER_FILES}/IBSEB_TwoStreamProviderTwoLevel.i" ARGS ${two_level})
run_leg(two_level "${PROVIDER_FILES}/IBSEB_TwoStreamProviderTwoLevel.i" ARGS ${two_level} ${on})

# The deck's grid (16 layers; the cube is 4 cells high), its sun and sky, the ground's and
# the roofs' albedo and emissivity, and the precision of ERF's Real.
string(TOLOWER "${PRECISION}" precision)
execute_process(
    COMMAND "${PYTHON_EXE}" "${CHECKER}" "${WORKING_DIRECTORY}" --nx 32 --ny 32 --nz 16 --cosz 0.766044443118978
            --toa 1000 --tau 0.02 --ssa 0.5 --tau-lw 0.3 --albedo 0.2 --emissivity 0.95 --roof-albedo 0.3 --roof-emissivity 0.9
            --precision ${precision}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/checker.log"
    ERROR_FILE "${WORKING_DIRECTORY}/checker.log"
    RESULT_VARIABLE check_result)
file(READ "${WORKING_DIRECTORY}/checker.log" check_output)
message("${check_output}")
if(NOT check_result EQUAL 0)
    message(FATAL_ERROR "RunIBSEBTwoStreamRoofs.cmake: the checker failed: ${check_result}")
endif()
