# Run a multi-level anelastic deck and check, with Tests/check_terrain_multilevel.py,
# that every projection on level LEVEL converged: with erf.mg_v = 1 the log reports the
# divergence before and after each solve and the compatibility constant subtracted from
# the right-hand side, and the divergence left after a converged solve is that constant.
#
# Arguments (all -D):
#   MPIEXEC, MPIEXEC_NUMPROC_FLAG, MPIEXEC_PREFLAGS, NRANKS, TEST_EXE, CONFIG
#   INPUT, WORKING_DIRECTORY, RUN_TIMEOUT, RUNTIME_OPTIONS
#   LEVEL, MIN_SOLVES, TOL          (passed to the checker)
#   PYTHON_EXECUTABLE, CHECK_SCRIPT

cmake_policy(SET CMP0053 NEW)
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

foreach(arg NRANKS TEST_EXE INPUT WORKING_DIRECTORY RUN_TIMEOUT LEVEL MIN_SOLVES TOL PYTHON_EXECUTABLE CHECK_SCRIPT)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunTerrainMultiLevel.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunTerrainMultiLevel.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunTerrainMultiLevel.cmake: ERF executable")

separate_arguments(runtime_options UNIX_COMMAND "${RUNTIME_OPTIONS}")

erf_mpi_launcher_command(launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunTerrainMultiLevel.cmake")

# The deck names its auxiliary files relative to the run directory
file(GLOB _deck_files LIST_DIRECTORIES false "${WORKING_DIRECTORY}/*")
set(RUN_DIR "${WORKING_DIRECTORY}/run")
file(REMOVE_RECURSE "${RUN_DIR}")
file(MAKE_DIRECTORY "${RUN_DIR}")
file(COPY ${_deck_files} DESTINATION "${RUN_DIR}")

set(log "${RUN_DIR}/simulation.log")
execute_process(
    COMMAND ${launch} ${TEST_EXE} ${INPUT} ${runtime_options}
    WORKING_DIRECTORY "${RUN_DIR}"
    OUTPUT_FILE "${log}"
    ERROR_FILE "${log}"
    TIMEOUT ${RUN_TIMEOUT}
    RESULT_VARIABLE run_result)
if(NOT run_result EQUAL 0)
    message(FATAL_ERROR "RunTerrainMultiLevel.cmake: the run failed or exceeded ${RUN_TIMEOUT} s: ${run_result} (see ${log})")
endif()

execute_process(
    COMMAND ${PYTHON_EXECUTABLE} ${CHECK_SCRIPT} --level ${LEVEL} --min-solves ${MIN_SOLVES} --tol ${TOL} ${log}
    RESULT_VARIABLE check_result
    OUTPUT_VARIABLE check_output
    ERROR_VARIABLE check_output)
message(STATUS "${check_output}")
if(NOT check_result EQUAL 0)
    message(FATAL_ERROR "RunTerrainMultiLevel.cmake: the divergence check failed (${check_result})")
endif()
