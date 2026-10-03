# Run one deck with a reference projection solver and again with each of the
# other solvers, and require the final plotfiles to agree within a tolerance.
#
# Each leg must also print a line matching its own regular expression (the
# solver's convergence report), so a leg whose option was silently ignored, and
# which would therefore trivially match the reference, fails instead.
#
# Arguments (all -D):
#   MPIEXEC, MPIEXEC_NUMPROC_FLAG, MPIEXEC_PREFLAGS, NRANKS, TEST_EXE, CONFIG
#   INPUT              the inputs file
#   WORKING_DIRECTORY  where the legs run (one subdirectory each)
#   FCOMPARE           the fcompare executable
#   PLTFILE            the plotfile to compare, e.g. plt00050
#   RUN_TIMEOUT        seconds allowed for each run
#   COMMON_OPTIONS     runtime options for every leg
#   REF_OPTIONS        runtime options of the reference leg
#   REF_REQUIRE        regular expression the reference log must match
#   LEG_OPTIONS        runtime options of the other legs, separated by '|'
#   LEG_REQUIRE        one regular expression per leg, separated by '|'
#   REL_TOL, ABS_TOL   fcompare tolerances

cmake_policy(SET CMP0053 NEW)
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

foreach(arg NRANKS TEST_EXE INPUT WORKING_DIRECTORY FCOMPARE PLTFILE RUN_TIMEOUT
            REF_OPTIONS REF_REQUIRE LEG_OPTIONS LEG_REQUIRE REL_TOL ABS_TOL)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunTerrainSolverParity.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunTerrainSolverParity.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunTerrainSolverParity.cmake: ERF executable")
erf_resolve_executable(FCOMPARE "${FCOMPARE}" CONFIG "${CONFIG}"
    CONTEXT "RunTerrainSolverParity.cmake: fcompare")

separate_arguments(common_options UNIX_COMMAND "${COMMON_OPTIONS}")
separate_arguments(ref_options    UNIX_COMMAND "${REF_OPTIONS}")

string(REPLACE "|" ";" leg_option_list  "${LEG_OPTIONS}")
string(REPLACE "|" ";" leg_require_list "${LEG_REQUIRE}")
list(LENGTH leg_option_list  n_legs)
list(LENGTH leg_require_list n_require)
if(NOT n_legs EQUAL n_require)
    message(FATAL_ERROR "RunTerrainSolverParity.cmake: LEG_OPTIONS has ${n_legs} legs but LEG_REQUIRE has ${n_require}")
endif()
foreach(leg ${leg_option_list})
    if("${leg}" STREQUAL "${REF_OPTIONS}")
        message(FATAL_ERROR "RunTerrainSolverParity.cmake: a leg is given the reference options, so its comparison would be trivial")
    endif()
endforeach()

erf_mpi_launcher_command(launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunTerrainSolverParity.cmake")
erf_mpi_launcher_command(launch_one
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunTerrainSolverParity.cmake")

# The decks name their auxiliary files (soundings, terrain) relative to the run directory
file(GLOB _deck_files LIST_DIRECTORIES false "${WORKING_DIRECTORY}/*")

function(run_leg dir require)
    file(REMOVE_RECURSE "${dir}")
    file(MAKE_DIRECTORY "${dir}")
    file(COPY ${_deck_files} DESTINATION "${dir}")
    execute_process(
        COMMAND ${launch} ${TEST_EXE} ${INPUT} ${ARGN}
        WORKING_DIRECTORY "${dir}"
        OUTPUT_FILE "${dir}/simulation.log"
        ERROR_FILE "${dir}/simulation.log"
        TIMEOUT ${RUN_TIMEOUT}
        RESULT_VARIABLE _result)
    if(NOT _result EQUAL 0)
        message(FATAL_ERROR "RunTerrainSolverParity.cmake: the run in ${dir} failed or exceeded ${RUN_TIMEOUT} s: ${_result}")
    endif()
    if(NOT EXISTS "${dir}/${PLTFILE}/Header")
        message(FATAL_ERROR "RunTerrainSolverParity.cmake: no ${PLTFILE} in ${dir}")
    endif()
    file(STRINGS "${dir}/simulation.log" _matches REGEX "${require}")
    if("${_matches}" STREQUAL "")
        message(FATAL_ERROR "RunTerrainSolverParity.cmake: the log in ${dir} never matched '${require}', so the solver under test did not run")
    endif()
endfunction()

set(REF_DIR "${WORKING_DIRECTORY}/solver_ref")
run_leg("${REF_DIR}" "${REF_REQUIRE}" ${common_options} ${ref_options})

set(ileg 0)
foreach(leg ${leg_option_list})
    list(GET leg_require_list ${ileg} require)
    separate_arguments(leg_options UNIX_COMMAND "${leg}")
    set(LEG_DIR "${WORKING_DIRECTORY}/solver_leg${ileg}")
    run_leg("${LEG_DIR}" "${require}" ${common_options} ${leg_options})

    execute_process(
        COMMAND ${launch_one} ${FCOMPARE} --abort_if_not_all_found
                --rel_tol ${REL_TOL} --abs_tol ${ABS_TOL}
                ${REF_DIR}/${PLTFILE} ${LEG_DIR}/${PLTFILE}
        WORKING_DIRECTORY "${WORKING_DIRECTORY}"
        OUTPUT_FILE "${WORKING_DIRECTORY}/parity_leg${ileg}.log"
        ERROR_FILE "${WORKING_DIRECTORY}/parity_leg${ileg}.log"
        RESULT_VARIABLE parity_result)
    if(NOT parity_result EQUAL 0)
        message(FATAL_ERROR "RunTerrainSolverParity.cmake: ${PLTFILE} of '${leg}' differs from the reference beyond rel ${REL_TOL} / abs ${ABS_TOL}: ${parity_result} (see parity_leg${ileg}.log)")
    endif()
    message(STATUS "RunTerrainSolverParity: ${PLTFILE} of '${leg}' agrees with the reference within rel ${REL_TOL} / abs ${ABS_TOL}")
    math(EXPR ileg "${ileg} + 1")
endforeach()
