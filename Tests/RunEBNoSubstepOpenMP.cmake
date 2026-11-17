# Issue #5137: EB redistribution must run once per stage, outside the caller's
# OpenMP region. Compare the same cut-cell problem with one and two threads.
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

foreach(arg NRANKS TEST_EXE INPUT WORKING_DIRECTORY FCOMPARE MODE RTOL ATOL THREADS)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunEBNoSubstepOpenMP.cmake: ${arg} is required")
    endif()
endforeach()
erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}")
erf_resolve_executable(FCOMPARE "${FCOMPARE}" CONFIG "${CONFIG}")
erf_mpi_launcher_command(launch
    LAUNCHER "${MPIEXEC}" NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS} PREFLAGS "${MPIEXEC_PREFLAGS}")
erf_mpi_launcher_command(launch_one
    LAUNCHER "${MPIEXEC}" NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1 PREFLAGS "${MPIEXEC_PREFLAGS}")

# Both cases use several tiles per box, so the tiled update after the
# redistribution is exercised.
if(MODE STREQUAL "Anelastic")
    # Converge the projection tightly enough that its stopping tolerance does
    # not dominate the thread-count comparison of this weakly forced flow.
    # The MLMG reductions sum in a thread-dependent order, so the runs agree
    # only to roundoff (about 1e-11 relative); RTOL/ATOL are set for that.
    set(case_options erf.timestepping_type=MidPoint erf.fixed_dt=0.01
        erf.poisson_reltol=1.e-12 erf.poisson_abstol=1.e-17
        amr.n_cell=32 4 32 amr.max_grid_size=16 fabarray.mfiter_tile_size=8 2 8)
elseif(MODE STREQUAL "Compressible")
    # dt = 4e-6 keeps the explicit acoustic CFL below 0.5 (dx = dz = 3.125e-3,
    # c ~ 347 m/s); the runs are bit-for-bit equal across thread counts.
    set(case_options erf.timestepping_type=RK3_NoSubstepping erf.fixed_dt=4.e-6
        amr.max_grid_size=32 fabarray.mfiter_tile_size=16 2 8)
else()
    message(FATAL_ERROR "Unknown mode: ${MODE}")
endif()

foreach(threads IN ITEMS 1 ${THREADS})
    set(run_dir "${WORKING_DIRECTORY}/omp_${threads}")
    file(REMOVE_RECURSE "${run_dir}")
    file(MAKE_DIRECTORY "${run_dir}")
    execute_process(
        COMMAND "${CMAKE_COMMAND}" -E env OMP_NUM_THREADS=${threads}
                OMP_DYNAMIC=FALSE OMP_THREAD_LIMIT=${THREADS} OMP_MAX_ACTIVE_LEVELS=1
                ${launch} "${TEST_EXE}" "${INPUT}" ${case_options} max_step=4
                erf.check_int=-1 erf.write_eb_surface=0
                erf.plot_file_1=plt erf.plot_int_1=4
        WORKING_DIRECTORY "${run_dir}"
        OUTPUT_FILE "${run_dir}/simulation.log"
        ERROR_FILE "${run_dir}/simulation.log"
        TIMEOUT 120
        RESULT_VARIABLE result)
    if(NOT "${result}" STREQUAL "0")
        message(FATAL_ERROR "${MODE}, ${threads} threads failed: ${result}; see ${run_dir}/simulation.log")
    endif()
    file(READ "${run_dir}/simulation.log" log)
    if(NOT log MATCHES "OMP initialized with ${threads} OMP threads")
        message(FATAL_ERROR "The run did not initialize ${threads} OpenMP threads")
    endif()
    if(NOT EXISTS "${run_dir}/plt00004/Header")
        message(FATAL_ERROR "Missing final plotfile in ${run_dir}")
    endif()
endforeach()

execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env OMP_NUM_THREADS=1
            ${launch_one} "${FCOMPARE}" --abort_if_not_all_found
            --rel_tol ${RTOL} --abs_tol ${ATOL}
            "${WORKING_DIRECTORY}/omp_1/plt00004"
            "${WORKING_DIRECTORY}/omp_${THREADS}/plt00004"
    OUTPUT_FILE "${WORKING_DIRECTORY}/parity.log"
    ERROR_FILE "${WORKING_DIRECTORY}/parity.log"
    TIMEOUT 60
    RESULT_VARIABLE result)
if(NOT "${result}" STREQUAL "0")
    message(FATAL_ERROR "${MODE} depends on OpenMP thread count: ${result}; see ${WORKING_DIRECTORY}/parity.log")
endif()
