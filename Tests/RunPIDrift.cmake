# Run a pseudo-incompressible deck and check that (rho theta) stays rho0 theta0.  ERF prints
# "PI DRIFT = <value>", the maximum of |rho theta - rho0 theta0| / (rho0 theta0) over the
# pseudo-incompressible levels, every erf.sum_interval steps (with erf.v >= 1).  The last value
# printed must be below 10^MAX_EXP: zero, or a number whose decimal exponent is below MAX_EXP.
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

# -DX= defines X as empty, so test for a value, not for DEFINED
foreach(arg NRANKS TEST_EXE INPUT WORKING_DIRECTORY MAX_EXP)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunPIDrift.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunPIDrift.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunPIDrift.cmake: ERF executable")

separate_arguments(options UNIX_COMMAND "${OPTIONS}")

erf_mpi_launcher_command(launch_run
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunPIDrift.cmake")

# Run where the test's files were copied, so that the deck finds the ones it names
set(log "${WORKING_DIRECTORY}/simulation.log")
execute_process(
    COMMAND ${launch_run} ${TEST_EXE} ${INPUT} ${options}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${log}"
    ERROR_FILE "${log}"
    RESULT_VARIABLE run_result)
if(NOT run_result EQUAL 0)
    message(FATAL_ERROR "RunPIDrift.cmake: simulation failed: ${run_result}")
endif()

file(STRINGS "${log}" lines REGEX "^ *PI DRIFT += ")
list(LENGTH lines nlines)
if(nlines LESS 2)
    message(FATAL_ERROR "RunPIDrift.cmake: found ${nlines} lines \"PI DRIFT = ...\" in ${log}, "
                        "need at least two (erf.v >= 1 and erf.sum_interval > 0?)")
endif()
list(GET lines -1 last_line)
if(NOT last_line MATCHES "^ *PI DRIFT += +([-+0-9.eE]+)")
    message(FATAL_ERROR "RunPIDrift.cmake: cannot read the drift from: ${last_line}")
endif()
set(drift "${CMAKE_MATCH_1}")
message(STATUS "RunPIDrift: PI DRIFT = ${drift} after ${nlines} values (limit 1e${MAX_EXP})")

if(drift MATCHES "^0(\\.0*)?$")
    return()
endif()
if(NOT drift MATCHES "^[0-9.]+[eE]([-+]?[0-9]+)$")
    message(FATAL_ERROR "RunPIDrift.cmake: PI DRIFT = ${drift} is not below 1e${MAX_EXP}")
endif()
if(NOT CMAKE_MATCH_1 LESS MAX_EXP)
    message(FATAL_ERROR "RunPIDrift.cmake: PI DRIFT = ${drift} is not below 1e${MAX_EXP}")
endif()
