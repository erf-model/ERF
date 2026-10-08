# Run Tests/test_files/IBSEB_TwoStreamSun (the faces on erf.ibseb.sun_mode = two_stream)
# for a few steps with a report every step, and check with check_ibseb_two_stream_sun.py
# that the faces' top-of-atmosphere shortwave and sun side match the two-stream sweep's
# diagnostics: the run-time path of ERF::ibseb_set_two_stream_sun() and the two_stream
# branch of IBFaceSet::compute_shortwave(), which the unit tests cover only as formulas.
# -DX= defines X as empty, so test for a value, not for DEFINED
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

foreach(arg NRANKS TEST_EXE INPUT WORKING_DIRECTORY PYTHON_EXE CHECKER)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunIBSEBTwoStreamSun.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunIBSEBTwoStreamSun.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunIBSEBTwoStreamSun.cmake: ERF executable")

erf_mpi_launcher_command(launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunIBSEBTwoStreamSun.cmake")

set(RUN_DIR "${WORKING_DIRECTORY}/run")
file(REMOVE_RECURSE "${RUN_DIR}")
file(MAKE_DIRECTORY "${RUN_DIR}")
get_filename_component(input_dir "${INPUT}" DIRECTORY)
file(COPY "${input_dir}/input_sounding" "${input_dir}/cube_40m_10m_32x32.txt" DESTINATION "${RUN_DIR}")

# The deck's step is 0.5 s and its sw_transmission 0.8; the checker is told both.
execute_process(
    COMMAND ${launch} ${TEST_EXE} ${INPUT} max_step=4 erf.ibseb.csv_int=1
    WORKING_DIRECTORY "${RUN_DIR}"
    OUTPUT_FILE "${RUN_DIR}/simulation.log"
    ERROR_FILE "${RUN_DIR}/simulation.log"
    TIMEOUT 600
    RESULT_VARIABLE run_result)
if(NOT run_result EQUAL 0)
    message(FATAL_ERROR "RunIBSEBTwoStreamSun.cmake: the run failed: ${run_result} (see ${RUN_DIR}/simulation.log)")
endif()

execute_process(
    COMMAND "${PYTHON_EXE}" "${CHECKER}" "${RUN_DIR}" --dt 0.5 --tau 0.8
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/checker.log"
    ERROR_FILE "${WORKING_DIRECTORY}/checker.log"
    RESULT_VARIABLE check_result)
file(READ "${WORKING_DIRECTORY}/checker.log" check_output)
message("${check_output}")
if(NOT check_result EQUAL 0)
    message(FATAL_ERROR "RunIBSEBTwoStreamSun.cmake: the checker failed: ${check_result}")
endif()
