include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")

if(NOT DEFINED TEST_EXE OR NOT DEFINED LOG)
    message(FATAL_ERROR
        "RunCloudChamberWallDtGuardFailure.cmake missing required argument")
endif()

# In an MPI build the checker must be launched through mpiexec like every other
# Cloud Chamber checker; a bare singleton launch hangs in MPI_Init instead of
# reaching the guard.  A serial build has no MPIEXEC and is launched directly.
erf_mpi_launcher_command(launch_command
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunCloudChamberWallDtGuardFailure.cmake")
list(APPEND launch_command "${TEST_EXE}")

execute_process(
    COMMAND ${launch_command}
    RESULT_VARIABLE result
    OUTPUT_VARIABLE stdout
    ERROR_VARIABLE stderr)

set(output "${stdout}\n${stderr}")
file(WRITE "${LOG}" "${output}")

if("${result}" STREQUAL "0")
    message(FATAL_ERROR
        "Cloud Chamber fixed-dt guard unexpectedly returned success")
endif()

foreach(expected IN ITEMS
        "Cloud Chamber wall-transfer timestep violation"
        "fixed_dt="
        "wall_dt="
        "max_wall_rate=")
    string(FIND "${output}" "${expected}" position)
    if(position EQUAL -1)
        message(FATAL_ERROR
            "Missing expected fixed-dt guard diagnostic: ${expected}\n${output}")
    endif()
endforeach()
