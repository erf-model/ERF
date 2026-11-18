# Run the deck of TwoStream_SEBLandMask three times -- erf.is_land = 1, = 0, and two levels
# with water below a nest of land -- then check
# with the case's check_land_mask.py that the two-stream force-restore balance holds the
# water columns, takes the surface layer's sea surface temperature as the longwave
# boundary there, leaves them out of the means, and keeps a coarse water column under
# fine land columns at the average-down.
# -DX= defines X as empty, so test for a value, not for DEFINED
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

foreach(arg NRANKS TEST_EXE INPUT WORKING_DIRECTORY FEXTRACT PYTHON_EXE CHECKER STEPS)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunTwoStreamSEBLandMask.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunTwoStreamSEBLandMask.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunTwoStreamSEBLandMask.cmake: ERF executable")
erf_resolve_executable(FEXTRACT "${FEXTRACT}" CONFIG "${CONFIG}"
    CONTEXT "RunTwoStreamSEBLandMask.cmake: fextract")

erf_mpi_launcher_command(launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunTwoStreamSEBLandMask.cmake")

get_filename_component(input_dir "${INPUT}" DIRECTORY)
set(LAND_DIR  "${WORKING_DIRECTORY}/land")
set(WATER_DIR "${WORKING_DIRECTORY}/water")
set(NEST_DIR  "${WORKING_DIRECTORY}/nest")
file(REMOVE_RECURSE "${LAND_DIR}" "${WATER_DIR}" "${NEST_DIR}")
file(MAKE_DIRECTORY "${LAND_DIR}" "${WATER_DIR}" "${NEST_DIR}")

function(run_leg dir is_land)
    # The deck reads its sounding by a relative name.
    file(COPY "${input_dir}/input_sounding" DESTINATION "${dir}")
    execute_process(
        COMMAND ${launch} ${TEST_EXE} ${INPUT} erf.is_land=${is_land} max_step=${STEPS} ${ARGN}
        WORKING_DIRECTORY "${dir}"
        OUTPUT_FILE "${dir}/simulation.log"
        ERROR_FILE "${dir}/simulation.log"
        TIMEOUT 600
        RESULT_VARIABLE _result)
    if(NOT _result EQUAL 0)
        message(FATAL_ERROR "RunTwoStreamSEBLandMask.cmake: the erf.is_land = ${is_land} run failed: ${_result} (see ${dir}/simulation.log)")
    endif()
endfunction()

run_leg("${LAND_DIR}" 1)
run_leg("${WATER_DIR}" 0)
# Level 0 is water and a level 1 of land covers its western half from the start. No
# subcycling, so both levels are at a step's end when the 2D plotfile is written.
run_leg("${NEST_DIR}" "0 1"
    amr.max_level=1 "amr.ref_ratio_vect=2 2 1" amr.n_error_buf=0 amr.blocking_factor=2
    amr.refine_whole_domain_dir=2 erf.refinement_indicators=patch erf.patch.max_level=1
    "erf.patch.in_box_lo=0.0 0.0" "erf.patch.in_box_hi=200.0 400.0"
    erf.regrid_int=1 erf.dt_ref_ratio=1)

# The fill of level 1 from level 0 (ERF::fill_seb_from_coarse): every level 1 column is land
# over water with no land next to it, so the fill reports all of them: 4 x 8 columns of 50 m.
file(READ "${NEST_DIR}/simulation.log" nest_log)
string(FIND "${nest_log}" "ERF::fill_seb_from_coarse: level 1: 32 land column(s) lie over a coarse water column" fill_line)
if(fill_line EQUAL -1)
    message(FATAL_ERROR "RunTwoStreamSEBLandMask.cmake: the two-level run did not report 32 level 1 land columns with no land source (see ${NEST_DIR}/simulation.log)")
endif()

execute_process(
    COMMAND "${PYTHON_EXE}" "${CHECKER}"
            --fextract "${FEXTRACT}"
            --land-dir "${LAND_DIR}"
            --water-dir "${WATER_DIR}"
            --nest-dir "${NEST_DIR}"
            --steps ${STEPS}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/checker.log"
    ERROR_FILE "${WORKING_DIRECTORY}/checker.log"
    RESULT_VARIABLE check_result)
file(READ "${WORKING_DIRECTORY}/checker.log" check_output)
message("${check_output}")
if(NOT check_result EQUAL 0)
    message(FATAL_ERROR "RunTwoStreamSEBLandMask.cmake: the checker failed: ${check_result}")
endif()
