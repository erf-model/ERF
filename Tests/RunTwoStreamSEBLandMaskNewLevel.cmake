# Run the deck of TwoStream_SEBLandMaskNewLevel four times on a wrfinput land mask with a
# shore -- level 1 made by a regrid at t = 5 s, level 1 from the start, level 1 from the
# start with its own wrfinput file, and the first run restarted from a checkpoint written
# after the regrid -- then check with the case's
# check_new_level.py that the level 1 columns over coarse water hold the parent's values
# and that the fill counts no land column over water.
# -DX= defines X as empty, so test for a value, not for DEFINED
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

foreach(arg NRANKS TEST_EXE INPUT WORKING_DIRECTORY PYTHON_EXE GENERATOR CHECKER READER_DIR)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunTwoStreamSEBLandMaskNewLevel.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunTwoStreamSEBLandMaskNewLevel.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunTwoStreamSEBLandMaskNewLevel.cmake: ERF executable")

erf_mpi_launcher_command(launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunTwoStreamSEBLandMaskNewLevel.cmake")

set(REGRID_DIR  "${WORKING_DIRECTORY}/regrid")
set(SCRATCH_DIR "${WORKING_DIRECTORY}/scratch")
set(NESTFILE_DIR "${WORKING_DIRECTORY}/nestfile")
set(RESTART_DIR "${WORKING_DIRECTORY}/restart")
file(REMOVE_RECURSE "${REGRID_DIR}" "${SCRATCH_DIR}" "${NESTFILE_DIR}" "${RESTART_DIR}")
file(MAKE_DIRECTORY "${REGRID_DIR}" "${SCRATCH_DIR}" "${NESTFILE_DIR}" "${RESTART_DIR}")

# wrfinput_d01 and wrflowinp_d01, written once and copied to each leg: the deck reads
# them by a relative name.
execute_process(
    COMMAND "${PYTHON_EXE}" "${GENERATOR}"
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/generate.log"
    ERROR_FILE "${WORKING_DIRECTORY}/generate.log"
    RESULT_VARIABLE _result)
if(NOT _result EQUAL 0)
    message(FATAL_ERROR "RunTwoStreamSEBLandMaskNewLevel.cmake: make_wrfinput.py failed: ${_result} (see ${WORKING_DIRECTORY}/generate.log)")
endif()

function(run_leg dir)
    file(COPY "${WORKING_DIRECTORY}/wrfinput_d01" "${WORKING_DIRECTORY}/wrflowinp_d01"
         DESTINATION "${dir}")
    execute_process(
        COMMAND ${launch} ${TEST_EXE} ${INPUT} ${ARGN}
        WORKING_DIRECTORY "${dir}"
        OUTPUT_FILE "${dir}/simulation.log"
        ERROR_FILE "${dir}/simulation.log"
        TIMEOUT 180
        RESULT_VARIABLE _result)
    if(NOT _result EQUAL 0)
        message(FATAL_ERROR "RunTwoStreamSEBLandMaskNewLevel.cmake: the run in ${dir} failed: ${_result} (see ${dir}/simulation.log)")
    endif()
endfunction()

# Level 1 appears at the regrid after step 5 (the deck's erf.patch.start_time = 4.5 s);
# the checkpoint of step 6 holds both levels.
run_leg("${REGRID_DIR}" erf.check_int=6)
if(NOT EXISTS "${REGRID_DIR}/chk00006/Level_1")
    message(FATAL_ERROR "RunTwoStreamSEBLandMaskNewLevel.cmake: the regrid run wrote no level 1 in chk00006")
endif()
# Level 1 from the start: the patch is tagged at t = 0.
run_leg("${SCRATCH_DIR}" erf.patch.start_time=0.0)
# Level 1 from the start with its own file: it reads LANDMASK from wrfinput_d02.
file(COPY "${WORKING_DIRECTORY}/wrfinput_d02" "${WORKING_DIRECTORY}/wrflowinp_d02"
     DESTINATION "${NESTFILE_DIR}")
run_leg("${NESTFILE_DIR}" erf.patch.start_time=0.0
        erf.nc_init_file_1=wrfinput_d02 erf.nc_low_file_1=wrflowinp_d02)
run_leg("${RESTART_DIR}" "erf.restart=${REGRID_DIR}/chk00006")

execute_process(
    COMMAND "${PYTHON_EXE}" "${CHECKER}"
            --reader-dir "${READER_DIR}"
            --regrid-dir "${REGRID_DIR}"
            --scratch-dir "${SCRATCH_DIR}"
            --nestfile-dir "${NESTFILE_DIR}"
            --restart-dir "${RESTART_DIR}"
            --first-step 6
            --last-step 8
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/checker.log"
    ERROR_FILE "${WORKING_DIRECTORY}/checker.log"
    RESULT_VARIABLE check_result)
file(READ "${WORKING_DIRECTORY}/checker.log" check_output)
message("${check_output}")
if(NOT check_result EQUAL 0)
    message(FATAL_ERROR "RunTwoStreamSEBLandMaskNewLevel.cmake: the checker failed: ${check_result}")
endif()
