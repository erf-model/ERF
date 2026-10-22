# Run a multilevel deck and check that a conserved quantity has the same composite integral at
# the start and the end of the run.  The integral is the one ERF prints every erf.sum_interval
# steps (with erf.v >= 1) as "<QUANTITY> SL/ML = <single level> <multilevel>", for instance
# QUANTITY = "MASS", "RHO SCALAR" or "RHO QTOTAL".  ERF pads the shorter names ("MASS       SL/ML"),
# and the line must start with QUANTITY, so that "MASS" does not also match "PERT MASS".  The
# multilevel value takes the fine data where it
# exists and weights each cell by its volume, detJ / (m_x m_y), so it is the quantity the
# discretization conserves on any mesh and with any refinement ratio.  The first and the last
# value printed must agree to within DIGITS significant digits, give or take UNITS units of the
# last of them; the deck must print the integral at the start of the run.
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/CompareDataLogs.cmake")

# -DX= defines X as empty, so test for a value, not for DEFINED
foreach(arg NRANKS TEST_EXE INPUT WORKING_DIRECTORY QUANTITY DIGITS UNITS)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunConservation.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunConservation.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunConservation.cmake: ERF executable")

separate_arguments(options UNIX_COMMAND "${OPTIONS}")

erf_mpi_launcher_command(launch_run
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunConservation.cmake")

# Run where the test's files were copied, so that the deck finds the ones it names (an input
# sounding or a terrain file, say)
set(log "${WORKING_DIRECTORY}/simulation.log")
execute_process(
    COMMAND ${launch_run} ${TEST_EXE} ${INPUT} ${options}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${log}"
    ERROR_FILE "${log}"
    RESULT_VARIABLE run_result)
if(NOT run_result EQUAL 0)
    message(FATAL_ERROR "RunConservation.cmake: simulation failed: ${run_result}")
endif()

file(STRINGS "${log}" lines REGEX "^ *${QUANTITY} +SL/ML = ")
list(LENGTH lines nlines)
if(nlines LESS 2)
    message(FATAL_ERROR "RunConservation.cmake: found ${nlines} lines \"${QUANTITY} SL/ML = ...\" "
                        "in ${log}, need at least two (is the deck multilevel, with erf.v >= 1 "
                        "and erf.sum_interval > 0?)")
endif()

function(multilevel_integral line out_var)
    if(NOT line MATCHES "^ *${QUANTITY} +SL/ML = +[^ ]+ +([-+0-9.eE]+)")
        message(FATAL_ERROR "RunConservation.cmake: cannot read the integral from: ${line}")
    endif()
    set(${out_var} "${CMAKE_MATCH_1}" PARENT_SCOPE)
endfunction()

list(GET lines 0 first_line)
list(GET lines -1 last_line)
multilevel_integral("${first_line}" first)
multilevel_integral("${last_line}" last)
message(STATUS "RunConservation: ${QUANTITY} = ${first} at the start, ${last} at the end "
               "(${nlines} values)")

erf_numbers_close("${first}" "${last}" ${DIGITS} ${UNITS} close)
if(NOT close)
    message(FATAL_ERROR "RunConservation.cmake: ${QUANTITY} is not conserved: ${first} at the "
                        "start, ${last} at the end")
endif()
