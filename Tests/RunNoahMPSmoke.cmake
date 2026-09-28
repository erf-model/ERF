# Run the Noah-MP smoke test.
#
# Noah-MP's driver is WRF's and reads three things from the run directory: namelist.erf
# (copied in with the test files), NoahmpTable.TBL (copied here from the Noah-MP
# submodule rather than duplicated in the repository) and a wrfinput-format land setup
# file (generated here from CDL text with ncgen). Then ERF runs and the checker reads
# Noah-MP's own outputs from the last plotfile.
#
# The exit status of the run is necessary but not sufficient: a Fortran STOP inside
# Noah-MP exits with status 0, so the checker also requires the plotfile the run should
# have reached. Stale plotfiles are removed first, so one left by an earlier run cannot
# stand in for it.
# -DX= defines X as empty, so test for a value, not for DEFINED.
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")

foreach(arg TEST_EXE INPUT WORKING_DIRECTORY NCGEN CDL TABLE FEXTREMA PYTHON_EXE CHECKER
            PLOTFILE RANGES SIMULATION_LOG CHECKER_LOG)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunNoahMPSmoke.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()

function(noahmp_report_log label path)
    if(EXISTS "${path}")
        file(READ "${path}" contents)
        message(STATUS "---- ${label} (${path}) ----\n${contents}\n---- end ${label} ----")
    else()
        message(STATUS "---- ${label}: ${path} was never written ----")
    endif()
endfunction()

# Inputs the driver opens by name.
if(NOT EXISTS "${TABLE}")
    message(FATAL_ERROR "RunNoahMPSmoke.cmake: ${TABLE} not found; is Submodules/Noah-MP checked out?")
endif()
configure_file("${TABLE}" "${WORKING_DIRECTORY}/NoahmpTable.TBL" COPYONLY)

file(REMOVE "${WORKING_DIRECTORY}/wrfinput_d01")
execute_process(
    COMMAND ${NCGEN} -o wrfinput_d01 ${CDL}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    RESULT_VARIABLE ncgen_result
    OUTPUT_VARIABLE ncgen_out
    ERROR_VARIABLE ncgen_out)
if(NOT ncgen_result EQUAL 0 OR NOT EXISTS "${WORKING_DIRECTORY}/wrfinput_d01")
    message(FATAL_ERROR "RunNoahMPSmoke.cmake: ncgen could not build wrfinput_d01 from ${CDL} "
                        "(exit ${ncgen_result}):\n${ncgen_out}")
endif()

# Clear every plotfile an earlier run may have left behind.
file(GLOB stale_plotfiles LIST_DIRECTORIES true "${WORKING_DIRECTORY}/plt*")
if(stale_plotfiles)
    file(REMOVE_RECURSE ${stale_plotfiles})
endif()

if(NOT DEFINED NRANKS OR "${NRANKS}" STREQUAL "")
    set(NRANKS 1)
endif()
erf_mpi_launcher_command(launcher
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunNoahMPSmoke.cmake")
execute_process(
    COMMAND ${launcher} ${TEST_EXE} ${INPUT}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${SIMULATION_LOG}"
    ERROR_FILE "${SIMULATION_LOG}"
    RESULT_VARIABLE simulation_result)
if(NOT simulation_result EQUAL 0)
    noahmp_report_log("simulation log" "${SIMULATION_LOG}")
    message(FATAL_ERROR "Noah-MP smoke run failed: ${simulation_result}")
endif()

string(REPLACE "," ";" range_list "${RANGES}")
set(range_args "")
foreach(r IN LISTS range_list)
    list(APPEND range_args --range ${r})
endforeach()
execute_process(
    COMMAND ${PYTHON_EXE} ${CHECKER} --plotfile ${PLOTFILE} --fextrema ${FEXTREMA} ${range_args}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${CHECKER_LOG}"
    ERROR_FILE "${CHECKER_LOG}"
    RESULT_VARIABLE checker_result)
noahmp_report_log("checker log" "${CHECKER_LOG}")
if(checker_result EQUAL 1)
    message(FATAL_ERROR "Noah-MP smoke test: an output field is outside its physical bounds")
elseif(NOT checker_result EQUAL 0)
    noahmp_report_log("simulation log" "${SIMULATION_LOG}")
    message(FATAL_ERROR "Noah-MP smoke test could not check the outputs (${checker_result}); "
                        "see the checker log -- a missing plotfile means the run stopped early")
endif()
