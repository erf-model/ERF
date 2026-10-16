# Run a deck and check that the domain integral of a conserved variable is the same in the first
# and the last plotfile.  The integral is the composite one amrex_fvolumesum computes (fine data
# where it exists); it is read as the "integral of VAR = value" line fvolumesum prints, and the
# two values must agree to within two units of the last of the ten significant digits it prints.
#
# fvolumesum integrates with the level's uniform cell volumes and a single refinement ratio, so the
# deck must use a mesh with detJ = 1 and the same refinement ratio in every direction.
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/CompareDataLogs.cmake")

# -DX= defines X as empty, so test for a value, not for DEFINED
foreach(arg NRANKS TEST_EXE INPUT WORKING_DIRECTORY FVOLUMESUM VAR FIRST_PLTFILE LAST_PLTFILE)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunConservation.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunConservation.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunConservation.cmake: ERF executable")
erf_resolve_executable(FVOLUMESUM "${FVOLUMESUM}" CONFIG "${CONFIG}"
    CONTEXT "RunConservation.cmake: fvolumesum")

separate_arguments(options UNIX_COMMAND "${OPTIONS}")

# Run where the test's files were copied, so that the deck finds the ones it names (an input
# sounding, say); remove the plotfiles of an earlier run first
set(RUN_DIR "${WORKING_DIRECTORY}")
file(REMOVE_RECURSE "${RUN_DIR}/${FIRST_PLTFILE}" "${RUN_DIR}/${LAST_PLTFILE}")

erf_mpi_launcher_command(launch_run
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunConservation.cmake")
erf_mpi_launcher_command(launch_one
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunConservation.cmake")

execute_process(
    COMMAND ${launch_run} ${TEST_EXE} ${INPUT} ${options}
    WORKING_DIRECTORY "${RUN_DIR}"
    OUTPUT_FILE "${RUN_DIR}/simulation.log"
    ERROR_FILE "${RUN_DIR}/simulation.log"
    RESULT_VARIABLE run_result)
if(NOT run_result EQUAL 0)
    message(FATAL_ERROR "RunConservation.cmake: simulation failed: ${run_result}")
endif()

function(integral_of plotfile out_var)
    if(NOT EXISTS "${RUN_DIR}/${plotfile}/Header")
        message(FATAL_ERROR "RunConservation.cmake: no plotfile ${RUN_DIR}/${plotfile}")
    endif()
    execute_process(
        COMMAND ${launch_one} ${FVOLUMESUM} -v ${VAR} ${RUN_DIR}/${plotfile}
        WORKING_DIRECTORY "${RUN_DIR}"
        OUTPUT_VARIABLE _out
        ERROR_VARIABLE _err
        RESULT_VARIABLE _result)
    if(NOT _result EQUAL 0)
        message(FATAL_ERROR "RunConservation.cmake: fvolumesum failed on ${plotfile}: ${_err}")
    endif()
    if(NOT _out MATCHES "integral of ${VAR} = ([-+0-9.eE]+)")
        message(FATAL_ERROR "RunConservation.cmake: cannot read the integral from: ${_out}")
    endif()
    set(${out_var} "${CMAKE_MATCH_1}" PARENT_SCOPE)
endfunction()

integral_of("${FIRST_PLTFILE}" first)
integral_of("${LAST_PLTFILE}"  last)
message(STATUS "RunConservation: integral of ${VAR} = ${first} in ${FIRST_PLTFILE}, ${last} in ${LAST_PLTFILE}")

erf_numbers_close("${first}" "${last}" 10 2 close)
if(NOT close)
    message(FATAL_ERROR "RunConservation.cmake: ${VAR} is not conserved: ${first} in ${FIRST_PLTFILE}, ${last} in ${LAST_PLTFILE}")
endif()
