# Run a terrain-fitted anelastic deck and check with check_terrain_poisson_divergence.py
# that the projection on LEVEL converged after every solve (divergence after the solve at
# most MAX_RATIO of the divergence before it).
#
# Arguments (all -D):
#   NRANKS TEST_EXE INPUT WORKING_DIRECTORY PYTHON_EXE CHECKER LEVEL  required
#                   (LEVEL may list several levels, space separated: each one is checked)
#   MAX_RATIO       largest allowed divergence after / before a solve (default 1e-5)
#   MAX_MEAN        largest allowed |mean| subtracted before a solve (net inflow; default no limit)
#   COMMON_OPTIONS  options added to the deck (space separated, no quotes inside)
#   MPIEXEC MPIEXEC_NUMPROC_FLAG MPIEXEC_PREFLAGS CONFIG  as the other runners
# -DX= defines X as empty, so test for a value, not for DEFINED
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

foreach(arg NRANKS TEST_EXE INPUT WORKING_DIRECTORY PYTHON_EXE CHECKER LEVEL)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunTerrainPoissonDivergence.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunTerrainPoissonDivergence.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()
if("${MAX_RATIO}" STREQUAL "")
    set(MAX_RATIO "1.0e-5")
endif()

set(mean_opt "")
if(NOT "${MAX_MEAN}" STREQUAL "")
    set(mean_opt --max-mean ${MAX_MEAN})
endif()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunTerrainPoissonDivergence.cmake: ERF executable")

erf_mpi_launcher_command(launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunTerrainPoissonDivergence.cmake")

get_filename_component(input_dir "${INPUT}" DIRECTORY)
separate_arguments(common_opts UNIX_COMMAND "${COMMON_OPTIONS}")

# Run the deck in its own directory and return the checker's exit status in <rc_var>.
function(run_leg name rc_var)
    set(dir "${WORKING_DIRECTORY}/${name}")
    file(REMOVE_RECURSE "${dir}")
    file(MAKE_DIRECTORY "${dir}")
    # The deck reads its sounding by a relative name.
    file(COPY "${input_dir}/input_sounding" DESTINATION "${dir}")
    execute_process(
        COMMAND ${launch} ${TEST_EXE} ${INPUT} ${common_opts}
        WORKING_DIRECTORY "${dir}"
        OUTPUT_FILE "${dir}/simulation.log"
        ERROR_FILE "${dir}/simulation.log"
        TIMEOUT 600
        RESULT_VARIABLE run_result)
    if(NOT run_result EQUAL 0)
        file(STRINGS "${dir}/simulation.log" abort_lines REGEX "Abort|abort|Error|error")
        message(FATAL_ERROR "RunTerrainPoissonDivergence.cmake: the run failed (${run_result}): "
                            "${abort_lines} (see ${dir}/simulation.log)")
    endif()
    separate_arguments(levels UNIX_COMMAND "${LEVEL}")
    set(worst 0)
    foreach(lev IN LISTS levels)
        execute_process(
            COMMAND "${PYTHON_EXE}" "${CHECKER}" "${dir}/simulation.log" --level ${lev} --max-ratio ${MAX_RATIO} ${mean_opt}
            OUTPUT_VARIABLE check_out ERROR_VARIABLE check_out RESULT_VARIABLE check_result)
        message("${check_out}")
        if(NOT check_result EQUAL 0)
            set(worst ${check_result})
        endif()
    endforeach()
    set(${rc_var} "${worst}" PARENT_SCOPE)
endfunction()

run_leg(run run_rc)
if(NOT run_rc EQUAL 0)
    message(FATAL_ERROR "RunTerrainPoissonDivergence.cmake: the projection at level ${LEVEL} did not converge "
                        "(checker exit ${run_rc}; see ${WORKING_DIRECTORY}/run/simulation.log)")
endif()
