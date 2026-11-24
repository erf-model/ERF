# erf.init_vels_from_checkpoint on grids other than the checkpoint's (issue 5376).
# Run the deck to chk00005/plt00005 on its own boxes, then start it again from time 0 on
# other boxes (and, with more than one rank, another rank count) with
# erf.init_vels_from_checkpoint = chk00005, and stop at step 0. The velocities of the
# second run's plt00000 must equal those of plt00005 exactly; the other fields are
# re-initialized, so they are not compared.
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

# -DX= defines X as empty, so test for a value, not for DEFINED
foreach(arg NRANKS TEST_EXE INPUT WORKING_DIRECTORY FCOMPARE)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunInitVelsFromCheckpoint.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunInitVelsFromCheckpoint.cmake: ERF executable")

erf_mpi_launcher_command(launch_chk
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunInitVelsFromCheckpoint.cmake")

# The second run uses one rank fewer when it can, so the distribution differs too
set(nranks_init ${NRANKS})
if(NRANKS GREATER 1)
    math(EXPR nranks_init "${NRANKS} - 1")
endif()
erf_mpi_launcher_command(launch_init
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${nranks_init}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunInitVelsFromCheckpoint.cmake")

file(REMOVE_RECURSE "${WORKING_DIRECTORY}/chk00000" "${WORKING_DIRECTORY}/chk00005"
                    "${WORKING_DIRECTORY}/plt00000" "${WORKING_DIRECTORY}/plt00005"
                    "${WORKING_DIRECTORY}/plt_init00000")

execute_process(
    COMMAND ${launch_chk} ${TEST_EXE} ${INPUT}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/chk_run.log"
    ERROR_FILE "${WORKING_DIRECTORY}/chk_run.log"
    RESULT_VARIABLE status)
if(NOT status EQUAL 0 OR NOT EXISTS "${WORKING_DIRECTORY}/chk00005/Header")
    message(FATAL_ERROR "Checkpoint run failed (${status}) or wrote no chk00005; see chk_run.log")
endif()

execute_process(
    COMMAND ${launch_init} ${TEST_EXE} ${INPUT}
            max_step=0 amr.max_grid_size=8 erf.check_int=-1
            erf.plot_file_1=plt_init erf.init_vels_from_checkpoint=chk00005
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/init_run.log"
    ERROR_FILE "${WORKING_DIRECTORY}/init_run.log"
    RESULT_VARIABLE status)
if(NOT status EQUAL 0 OR NOT EXISTS "${WORKING_DIRECTORY}/plt_init00000/Header")
    message(FATAL_ERROR "Run with erf.init_vels_from_checkpoint failed (${status}); see init_run.log")
endif()
file(READ "${WORKING_DIRECTORY}/init_run.log" init_log)
if(NOT init_log MATCHES "Reading vels only from native checkpoint chk00005")
    message(FATAL_ERROR "The second run did not read its velocities from chk00005; see init_run.log")
endif()

# The two runs must be on different grids, or the test does not test the fix
file(READ "${WORKING_DIRECTORY}/plt00005/Level_0/Cell_H" chk_cell_h)
file(READ "${WORKING_DIRECTORY}/plt_init00000/Level_0/Cell_H" init_cell_h)
if(chk_cell_h STREQUAL init_cell_h)
    message(FATAL_ERROR "The two runs used the same grids")
endif()

execute_process(
    COMMAND "${FCOMPARE}" -a "${WORKING_DIRECTORY}/plt00005" "${WORKING_DIRECTORY}/plt_init00000"
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_VARIABLE fcompare_out
    ERROR_VARIABLE fcompare_out)
foreach(var x_velocity y_velocity z_velocity)
    if(NOT fcompare_out MATCHES "${var}[ \t]+0[ \t]+0")
        message(FATAL_ERROR "${var} read from chk00005 differs from the checkpoint run:\n${fcompare_out}")
    endif()
endforeach()
message(STATUS "Velocities read onto other grids match the checkpoint exactly")
