# Run a deck straight to STEP_END, run it again to STEP_CHK with a checkpoint there,
# restart from that checkpoint to STEP_END, and require the restarted run's plotfile
# at STEP_END to equal the straight run's with fcompare. Each leg uses the forwarded
# RUN_TIMEOUT, and the enclosing CTest timeout is sized separately by the caller.
# -DX= defines X as empty, so test for a value, not for DEFINED
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")

foreach(arg NRANKS TEST_EXE INPUT WORKING_DIRECTORY FCOMPARE STEP_CHK STEP_END RTOL ATOL RUN_TIMEOUT)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunRestartParity.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunRestartParity.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

separate_arguments(common_options   UNIX_COMMAND "${COMMON_OPTIONS}")

set(STRAIGHT_DIR "${WORKING_DIRECTORY}/straight")
set(RESTART_DIR  "${WORKING_DIRECTORY}/restart")
file(REMOVE_RECURSE "${STRAIGHT_DIR}" "${RESTART_DIR}")
file(MAKE_DIRECTORY "${STRAIGHT_DIR}" "${RESTART_DIR}")

# MPIEXEC may be a multi-word command such as "flux run"; the helper splits
# it, validates the program and applies MPIEXEC_PREFLAGS. An empty MPIEXEC
# yields an empty prefix, so the runs stay serial.
erf_mpi_launcher_command(launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunRestartParity.cmake")
erf_mpi_launcher_command(launch_one
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunRestartParity.cmake")

# plotfile and checkpoint names carry the step padded to five digits
function(padded step out_var)
    set(_s "0000${step}")
    string(LENGTH "${_s}" _len)
    math(EXPR _start "${_len} - 5")
    string(SUBSTRING "${_s}" ${_start} 5 _s)
    set(${out_var} "${_s}" PARENT_SCOPE)
endfunction()
padded(${STEP_CHK} chk_step)
padded(${STEP_END} end_step)
set(PLTFILE "plt${end_step}")
set(CHKFILE "chk${chk_step}")

function(run_erf dir log timeout_s)
    execute_process(
        COMMAND ${launch} ${TEST_EXE} ${INPUT} ${common_options} ${ARGN}
        WORKING_DIRECTORY "${dir}"
        OUTPUT_FILE "${dir}/${log}"
        ERROR_FILE "${dir}/${log}"
        TIMEOUT ${timeout_s}
        RESULT_VARIABLE _result)
    if(NOT _result EQUAL 0)
        message(FATAL_ERROR "RunRestartParity.cmake: the run in ${dir} (${log}) failed or exceeded ${timeout_s} s: ${_result}")
    endif()
endfunction()

# straight to the end, no checkpoint
run_erf("${STRAIGHT_DIR}" "simulation.log" ${RUN_TIMEOUT}
        "max_step=${STEP_END}" "erf.check_int=-1" "erf.plot_int_1=${STEP_END}")
# to the checkpoint step, writing it there
run_erf("${RESTART_DIR}" "checkpoint.log" ${RUN_TIMEOUT}
        "max_step=${STEP_CHK}" "erf.check_int=${STEP_CHK}" "erf.plot_int_1=-1")
if(NOT EXISTS "${RESTART_DIR}/${CHKFILE}/Header")
    message(FATAL_ERROR "RunRestartParity.cmake: no ${CHKFILE} written by the checkpoint run")
endif()
# from the checkpoint to the end
run_erf("${RESTART_DIR}" "restart.log" ${RUN_TIMEOUT}
        "erf.restart=${CHKFILE}" "max_step=${STEP_END}" "erf.check_int=-1" "erf.plot_int_1=${STEP_END}")

foreach(dir "${STRAIGHT_DIR}" "${RESTART_DIR}")
    if(NOT EXISTS "${dir}/${PLTFILE}/Header")
        message(FATAL_ERROR "RunRestartParity.cmake: no ${PLTFILE} in ${dir}")
    endif()
endforeach()

execute_process(
    COMMAND ${launch_one} ${FCOMPARE} --abort_if_not_all_found
            --rel_tol ${RTOL} --abs_tol ${ATOL}
            ${STRAIGHT_DIR}/${PLTFILE} ${RESTART_DIR}/${PLTFILE}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/parity.log"
    ERROR_FILE "${WORKING_DIRECTORY}/parity.log"
    RESULT_VARIABLE parity_result)
if(NOT parity_result EQUAL 0)
    message(FATAL_ERROR "RunRestartParity.cmake: the restarted run's ${PLTFILE} differs from the straight run's: ${parity_result} (see parity.log)")
endif()
message(STATUS "RunRestartParity: restart from ${CHKFILE} reproduces ${PLTFILE}")
