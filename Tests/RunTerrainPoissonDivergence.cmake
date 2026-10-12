# Run a terrain-fitted anelastic deck and check with check_terrain_poisson_divergence.py
# that the projection on LEVEL converged after every solve (divergence after the solve at
# most MAX_RATIO of the divergence before it).  With OLD_OPTIONS the deck is run a second
# time with them, and there the checker must report a failure (exit 1): that leg runs the
# operator the fix replaced and shows the test can tell the two apart.
#
# Arguments (all -D):
#   NRANKS TEST_EXE INPUT WORKING_DIRECTORY PYTHON_EXE CHECKER LEVEL  required
#   MAX_RATIO       largest allowed divergence after / before a solve (default 1e-5)
#   COMMON_OPTIONS  options for every leg (space separated, no quotes inside)
#   NEW_OPTIONS     options of the leg that must converge
#   OLD_OPTIONS     options of the leg that must not converge (optional)
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

# Run one leg in its own directory and return the checker's exit status in <rc_var>.
function(run_leg name opts rc_var)
    set(dir "${WORKING_DIRECTORY}/${name}")
    file(REMOVE_RECURSE "${dir}")
    file(MAKE_DIRECTORY "${dir}")
    # The deck reads its sounding by a relative name.
    file(COPY "${input_dir}/input_sounding" DESTINATION "${dir}")
    separate_arguments(leg_opts UNIX_COMMAND "${opts}")
    execute_process(
        COMMAND ${launch} ${TEST_EXE} ${INPUT} ${common_opts} ${leg_opts}
        WORKING_DIRECTORY "${dir}"
        OUTPUT_FILE "${dir}/simulation.log"
        ERROR_FILE "${dir}/simulation.log"
        TIMEOUT 600
        RESULT_VARIABLE run_result)
    if(NOT run_result EQUAL 0)
        file(STRINGS "${dir}/simulation.log" abort_lines REGEX "Abort|abort|Error|error")
        message(FATAL_ERROR "RunTerrainPoissonDivergence.cmake: the ${name} run failed (${run_result}): "
                            "${abort_lines} (see ${dir}/simulation.log)")
    endif()
    execute_process(
        COMMAND "${PYTHON_EXE}" "${CHECKER}" "${dir}/simulation.log" --level ${LEVEL} --max-ratio ${MAX_RATIO}
        OUTPUT_VARIABLE check_out ERROR_VARIABLE check_out RESULT_VARIABLE check_result)
    message("${name} leg (${opts}):\n${check_out}")
    set(${rc_var} "${check_result}" PARENT_SCOPE)
endfunction()

run_leg(new "${NEW_OPTIONS}" new_rc)
if(NOT new_rc EQUAL 0)
    message(FATAL_ERROR "RunTerrainPoissonDivergence.cmake: the projection at level ${LEVEL} did not converge "
                        "with '${NEW_OPTIONS}' (checker exit ${new_rc})")
endif()

if(NOT "${OLD_OPTIONS}" STREQUAL "")
    run_leg(old "${OLD_OPTIONS}" old_rc)
    if(NOT old_rc EQUAL 1)
        message(FATAL_ERROR "RunTerrainPoissonDivergence.cmake: the leg with '${OLD_OPTIONS}' was expected "
                            "to leave a large divergence at level ${LEVEL} (checker exit 1), got exit ${old_rc}; "
                            "the deck no longer tells the two operators apart")
    endif()
endif()
