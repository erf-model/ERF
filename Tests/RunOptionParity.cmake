# Run one deck twice, once with a diagnostic switched off and once with it on, and require
# the two plotfiles to be identical.  This is how a claim that some output "does not change
# the answer" is held to: a diagnostic that perturbs the solution is a bug, not a tolerance
# question, so PLTFILE is compared at zero tolerance rather than at the ERF_TEST_FCOMPARE_*
# tolerances the gold-file tests use.
#
# Both legs run on the same number of ranks with the same decomposition, so anything that
# survives is the diagnostic's own doing.
#
# REQUIRE_ON_FILE, when given, is a path relative to each run directory that the "on" leg must
# produce and the "off" leg must not.  Without it a misspelled option would leave both legs
# identical and the test would pass while checking nothing.
# -DX= defines X as empty, so test for a value, not for DEFINED
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")

foreach(arg NRANKS TEST_EXE INPUT WORKING_DIRECTORY FCOMPARE PLTFILE RUN_TIMEOUT ON_OPTIONS)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunOptionParity.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunOptionParity.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()
if("${OFF_OPTIONS}" STREQUAL "${ON_OPTIONS}")
    message(FATAL_ERROR "RunOptionParity.cmake: the two legs are given the same options, so the comparison would be trivial")
endif()

separate_arguments(common_options UNIX_COMMAND "${COMMON_OPTIONS}")
separate_arguments(off_options    UNIX_COMMAND "${OFF_OPTIONS}")
separate_arguments(on_options     UNIX_COMMAND "${ON_OPTIONS}")

set(OFF_DIR "${WORKING_DIRECTORY}/option_off")
set(ON_DIR  "${WORKING_DIRECTORY}/option_on")
file(REMOVE_RECURSE "${OFF_DIR}" "${ON_DIR}")
file(MAKE_DIRECTORY "${OFF_DIR}" "${ON_DIR}")

# MPIEXEC may be a multi-word command such as "flux run"; the helper splits
# it, validates the program and applies MPIEXEC_PREFLAGS. An empty MPIEXEC
# yields an empty prefix, so the runs stay serial.
erf_mpi_launcher_command(launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunOptionParity.cmake")
erf_mpi_launcher_command(launch_one
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunOptionParity.cmake")

function(run_leg dir log timeout_s)
    execute_process(
        COMMAND ${launch} ${TEST_EXE} ${INPUT} ${ARGN}
        WORKING_DIRECTORY "${dir}"
        OUTPUT_FILE "${dir}/${log}"
        ERROR_FILE "${dir}/${log}"
        TIMEOUT ${timeout_s}
        RESULT_VARIABLE _result)
    if(NOT _result EQUAL 0)
        message(FATAL_ERROR "RunOptionParity.cmake: the run in ${dir} (${log}) failed or exceeded ${timeout_s} s: ${_result}")
    endif()
endfunction()

run_leg("${OFF_DIR}" "simulation.log" ${RUN_TIMEOUT} ${common_options} ${off_options})
run_leg("${ON_DIR}"  "simulation.log" ${RUN_TIMEOUT} ${common_options} ${on_options})

foreach(dir "${OFF_DIR}" "${ON_DIR}")
    if(NOT EXISTS "${dir}/${PLTFILE}/Header")
        message(FATAL_ERROR "RunOptionParity.cmake: no ${PLTFILE} in ${dir}")
    endif()
endforeach()

# The comparison proves nothing unless the "on" leg really did turn the diagnostic on
if(NOT "${REQUIRE_ON_FILE}" STREQUAL "")
    if(NOT EXISTS "${ON_DIR}/${REQUIRE_ON_FILE}")
        message(FATAL_ERROR "RunOptionParity.cmake: the on leg did not write ${REQUIRE_ON_FILE}, so the diagnostic under test never ran")
    endif()
    if(EXISTS "${OFF_DIR}/${REQUIRE_ON_FILE}")
        message(FATAL_ERROR "RunOptionParity.cmake: the off leg wrote ${REQUIRE_ON_FILE}, so the diagnostic under test was not switched off")
    endif()
endif()

execute_process(
    COMMAND ${launch_one} ${FCOMPARE} --abort_if_not_all_found
            --rel_tol 0 --abs_tol 0
            ${OFF_DIR}/${PLTFILE} ${ON_DIR}/${PLTFILE}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/parity.log"
    ERROR_FILE "${WORKING_DIRECTORY}/parity.log"
    RESULT_VARIABLE parity_result)
if(NOT parity_result EQUAL 0)
    message(FATAL_ERROR "RunOptionParity.cmake: ${PLTFILE} differs between the two runs, so the diagnostic changed the answer: ${parity_result} (see parity.log)")
endif()
message(STATUS "RunOptionParity: ${PLTFILE} is unchanged by '${ON_OPTIONS}'")
