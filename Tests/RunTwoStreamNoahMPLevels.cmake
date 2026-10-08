# Run Tests/test_files/TwoStream_NoahMPLevels four ways and check each with
# check_two_stream_noahmp_levels.py: two-stream radiation feeding Noah-MP on two levels.
#
#   own     level 1 runs Noah-MP on its own nested land file and sweeps its own columns;
#   nested  the same land, with a level that stops below the domain top and so takes its
#           radiation (Noah-MP's forcing included) from level 0;
#   interp  no level-1 land file: level 1 takes its land state from level 0;
#   regrid  the refined region moves, which a level running Noah-MP on its own land file
#           cannot do: the run must stop with the message that says so.
#
# The level-0 land file, the sounding and Noah-MP's table are those of
# Exec/RegTests/NoahMP_Ideal (LAND_DIR).
# -DX= defines X as empty, so test for a value, not for DEFINED
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

foreach(arg NRANKS TEST_EXE INPUT WORKING_DIRECTORY LAND_DIR FEXTRACT PYTHON_EXE CHECKER)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunTwoStreamNoahMPLevels.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunTwoStreamNoahMPLevels.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunTwoStreamNoahMPLevels.cmake: ERF executable")
erf_resolve_executable(FEXTRACT "${FEXTRACT}" CONFIG "${CONFIG}"
    CONTEXT "RunTwoStreamNoahMPLevels.cmake: fextract")

erf_mpi_launcher_command(launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunTwoStreamNoahMPLevels.cmake")

get_filename_component(input_dir "${INPUT}" DIRECTORY)
get_filename_component(input_name "${INPUT}" NAME)

# Stage a leg's run directory; the driver reads every file by a relative name.
function(stage_leg dir)
    file(REMOVE_RECURSE "${dir}")
    file(MAKE_DIRECTORY "${dir}")
    file(COPY "${INPUT}" "${input_dir}/namelist.erf" "${input_dir}/wrfinput_d02"
              "${LAND_DIR}/wrfinput_d01" "${LAND_DIR}/input_sounding"
              "${LAND_DIR}/NoahmpTable.TBL"
         DESTINATION "${dir}")
endfunction()

function(run_leg leg expect_success)
    set(dir "${WORKING_DIRECTORY}/${leg}")
    execute_process(
        COMMAND ${launch} ${TEST_EXE} ${input_name} ${ARGN}
        WORKING_DIRECTORY "${dir}"
        OUTPUT_FILE "${dir}/simulation.log"
        ERROR_FILE "${dir}/simulation.log"
        TIMEOUT 600
        RESULT_VARIABLE _result)
    if(expect_success AND NOT _result EQUAL 0)
        message(FATAL_ERROR "RunTwoStreamNoahMPLevels.cmake: the ${leg} run failed: ${_result} (see ${dir}/simulation.log)")
    endif()
    set(run_result "${_result}" PARENT_SCOPE)
endfunction()

function(check_leg leg)
    set(dir "${WORKING_DIRECTORY}/${leg}")
    execute_process(
        COMMAND "${PYTHON_EXE}" "${CHECKER}" --leg ${leg} --run-dir "${dir}"
                --fextract "${FEXTRACT}"
        WORKING_DIRECTORY "${dir}"
        OUTPUT_FILE "${dir}/checker.log"
        ERROR_FILE "${dir}/checker.log"
        RESULT_VARIABLE check_result)
    file(READ "${dir}/checker.log" check_output)
    message("${check_output}")
    if(NOT check_result EQUAL 0)
        message(FATAL_ERROR "RunTwoStreamNoahMPLevels.cmake: the checker failed on the ${leg} run: ${check_result}")
    endif()
endfunction()

stage_leg("${WORKING_DIRECTORY}/own")
run_leg(own TRUE)
check_leg(own)

# A level that stops at z = 300 m: no amr.refine_whole_domain_dir, and a 3D box.
stage_leg("${WORKING_DIRECTORY}/nested")
run_leg(nested TRUE
        amr.refine_whole_domain_dir=-1
        "erf.patch.in_box_lo=250.0 250.0 0.0" "erf.patch.in_box_hi=750.0 750.0 300.0")
check_leg(nested)

stage_leg("${WORKING_DIRECTORY}/interp")
file(STRINGS "${WORKING_DIRECTORY}/interp/namelist.erf" namelist_lines)
list(FILTER namelist_lines EXCLUDE REGEX "ERF_SETUP_FILE_02")
list(JOIN namelist_lines "\n" namelist_text)
file(WRITE "${WORKING_DIRECTORY}/interp/namelist.erf" "${namelist_text}\n")
run_leg(interp TRUE)
check_leg(interp)

# The region moves at t = 0.5 s, so the regrid before step 2 rebuilds level 1.
stage_leg("${WORKING_DIRECTORY}/regrid")
run_leg(regrid FALSE
        "erf.refinement_indicators=patch moved" erf.patch.end_time=0.5
        erf.moved.max_level=1 "erf.moved.in_box_lo=0.0 0.0" "erf.moved.in_box_hi=500.0 500.0"
        erf.moved.start_time=0.5 erf.regrid_int=1)
file(READ "${WORKING_DIRECTORY}/regrid/simulation.log" regrid_log)
set(expected "Regridding level 1 would rebuild its Noah-MP land state")
if(run_result EQUAL 0)
    message(FATAL_ERROR "RunTwoStreamNoahMPLevels.cmake: the regrid run completed; it must stop when level 1, which runs Noah-MP on its own land file, is rebuilt")
endif()
string(FIND "${regrid_log}" "${expected}" found)
if(found EQUAL -1)
    message(FATAL_ERROR "RunTwoStreamNoahMPLevels.cmake: the regrid run stopped (${run_result}) without the message '${expected}' (see ${WORKING_DIRECTORY}/regrid/simulation.log)")
endif()
message("regrid: stopped with '${expected}'")
