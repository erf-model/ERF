# Have CMake discover the number of cores on the node
include(ProcessorCount)
ProcessorCount(PROCESSES)

# Shared handling of multi-word MPI launchers (e.g., "flux run").
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")

#=============================================================================
# Functions for adding tests / Categories of tests
#=============================================================================
function(resolve_test_exe TEST_DIR TEST_EXE OUT_VAR)
    if(WIN32)
        # Multi-config generators place binaries in a config subdir.
        set(${OUT_VAR} "${CMAKE_BINARY_DIR}/Exec/${TEST_DIR}/*/${TEST_EXE}.exe" PARENT_SCOPE)
    else()
        set(_exe_in_subdir "${CMAKE_BINARY_DIR}/Exec/${TEST_DIR}/${TEST_EXE}${CMAKE_EXECUTABLE_SUFFIX}")
        set(_exe_in_root  "${CMAKE_BINARY_DIR}/Exec/${TEST_EXE}${CMAKE_EXECUTABLE_SUFFIX}")
        if(EXISTS "${_exe_in_subdir}")
            set(${OUT_VAR} "${_exe_in_subdir}" PARENT_SCOPE)
        elseif(EXISTS "${_exe_in_root}")
            set(${OUT_VAR} "${_exe_in_root}" PARENT_SCOPE)
        else()
            # Keep the historical path so the error message is still informative.
            set(${OUT_VAR} "${_exe_in_subdir}" PARENT_SCOPE)
        endif()
    endif()
endfunction()

macro(setup_test)
    if(DEFINED TEST_FILES_DIR AND NOT "${TEST_FILES_DIR}" STREQUAL "")
        set(_test_source_dir_name "${TEST_FILES_DIR}")
    else()
        set(_test_source_dir_name "${TEST_NAME}")
    endif()
    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${_test_source_dir_name})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    set(PLOT_GOLD ${ERF_TEST_GOLD_FILES_DIRECTORY}/${TEST_NAME})

    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")

    if(ERF_ENABLE_MPI)
        set(NP ${ERF_TEST_NRANKS})
        set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NP} ${MPIEXEC_PREFLAGS}")
        set(MPI_FCOMP_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}")
    else()
        set(NP 1)
        unset(MPI_COMMANDS)
        unset(MPI_FCOMP_COMMANDS)
    endif()
endmacro(setup_test)

# Production contract test for native 3D plotfile names and unavailable-name warnings.
function(add_test_plotfile_header TEST_NAME TEST_DIR TEST_EXE PLTFILE)
    setup_test()

    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)
    set(header_checker "${PROJECT_SOURCE_DIR}/Tests/CheckPlotfileHeader.cmake")
    set(test_log "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log")
    set(header_file "${CURRENT_TEST_BINARY_DIR}/${PLTFILE}/Header")
    # Regression motivation: this test is launched through `sh -c`, while on
    # Windows CMAKE_COMMAND normally resides below "C:/Program Files". Keep
    # checker executable and path-bearing arguments shell-quoted.
    set(check_command
        "\"${CMAKE_COMMAND}\""
        "\"-DHEADER=${header_file}\""
        "\"-DEXPECTED_NAMES_FILE=${CURRENT_TEST_BINARY_DIR}/expected_names.txt\""
        "\"-DLOG=${test_log}\""
        "\"-DEXPECTED_UNAVAILABLE_FILE=${CURRENT_TEST_BINARY_DIR}/expected_unavailable.txt\""
        "-P"
        "\"${header_checker}\"")
    list(JOIN check_command " " check_command_string)
    set(test_input "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i")
    set(test_command sh -c
        "${MPI_COMMANDS} ${TEST_EXE} \"${test_input}\" > \"${test_log}\" 2>&1 && ${check_command_string}")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 300
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;plotfile"
        ATTACHED_FILES_ON_FAIL "${test_log};${header_file}"
    )
endfunction(add_test_plotfile_header)

# Standard regression test
function(add_test_r TEST_NAME TEST_DIR TEST_EXE PLTFILE)
    set(options )
    set(oneValueArgs "INPUT_SOUNDING" "RUNTIME_OPTIONS" "FCOMPARE_RTOL" "FCOMPARE_ATOL")
    set(multiValueArgs )
    cmake_parse_arguments(ADD_TEST_R "${options}" "${oneValueArgs}"
        "${multiValueArgs}" ${ARGN})

    setup_test()

    set(RUNTIME_OPTIONS "${ADD_TEST_R_RUNTIME_OPTIONS}")
    if(NOT "${ADD_TEST_R_INPUT_SOUNDING}" STREQUAL "")
      string(APPEND RUNTIME_OPTIONS "erf.input_sounding_file=${CURRENT_TEST_BINARY_DIR}/${ADD_TEST_R_INPUT_SOUNDING}")
    endif()

    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)

    set(_fcompare_rtol "${ERF_TEST_FCOMPARE_RTOL}")
    set(_fcompare_atol "${ERF_TEST_FCOMPARE_ATOL}")
    if(NOT "${ADD_TEST_R_FCOMPARE_RTOL}" STREQUAL "")
        set(_fcompare_rtol "${ADD_TEST_R_FCOMPARE_RTOL}")
    endif()
    if(NOT "${ADD_TEST_R_FCOMPARE_ATOL}" STREQUAL "")
        set(_fcompare_atol "${ADD_TEST_R_FCOMPARE_ATOL}")
    endif()

    set(FCOMPARE_TOLERANCE "-r ${_fcompare_rtol} --abs_tol ${_fcompare_atol}")
    set(FCOMPARE_FLAGS "--abort_if_not_all_found -a ${FCOMPARE_TOLERANCE}")
    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUNTIME_OPTIONS} > ${TEST_NAME}.log && ${MPI_FCOMP_COMMANDS} ${FCOMPARE_EXE} ${FCOMPARE_FLAGS} ${PLOT_GOLD} ${CURRENT_TEST_BINARY_DIR}/${PLTFILE}")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log"
    )
endfunction(add_test_r)

# Rotated six-wall anelastic manufactured regression. Each case keeps a
# linear theta profile stationary and checks the full field inventory.
function(add_test_anelastic_wall_diffusion TEST_NAME TEST_AXIS)
    set(TEST_FILES_DIR "${TEST_NAME}")
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)

    set(test_input "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i")
    set(test_simulation_log "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.simulation.log")
    set(test_checker_log "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.checker.log")
    add_test(${TEST_NAME} ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DINPUT=${test_input}"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DSIMULATION_LOG=${test_simulation_log}"
        "-DCHECKER_LOG=${test_checker_log}"
        "-DCHECKER=${ANELASTIC_WALL_DIFFUSION_CHECKER}"
        "-DPLOTFILE=${CURRENT_TEST_BINARY_DIR}/plt00002"
        "-DAXIS=${TEST_AXIS}"
        "-DTHETA_LO=300.0"
        "-DTHETA_HI=301.0"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunAnelasticWallDiffusion.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;anelastic;wall-diffusion"
        ATTACHED_FILES_ON_FAIL "${test_simulation_log};${test_checker_log}")
endfunction(add_test_anelastic_wall_diffusion)

# Checker-driven Cloud Chamber tests.  The short run checks the exact initial
# conserved-state correction and a bounded early buoyant response; it
# intentionally avoids a fragile turbulent gold file.
function(add_test_cloud_chamber TEST_NAME MODE)
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)
    set(test_input "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i")
    set(test_simulation_log "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.simulation.log")
    set(test_checker_log "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.checker.log")
    add_test(${TEST_NAME} ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DINPUT=${test_input}"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DSIMULATION_LOG=${test_simulation_log}"
        "-DCHECKER_LOG=${test_checker_log}"
        "-DCHECKER=${CLOUD_CHAMBER_CHECKER}"
        "-DMODE=${MODE}"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunCloudChamber.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;cloud-chamber"
        ATTACHED_FILES_ON_FAIL "${test_simulation_log};${test_checker_log}")
endfunction(add_test_cloud_chamber)

# Gold-free TwoStream radiation regression: run a short SW + LW column case
# and verify the vertical structure of qsrc_sw / qsrc_lw in the plotfile
# (surface at k = 0, cooling to space from the top layer).
function(add_test_two_stream_radiation TEST_NAME PLTFILE)
    set(oneValueArgs "RUNTIME_OPTIONS")
    cmake_parse_arguments(ADD_TEST_TSR "" "${oneValueArgs}" "" ${ARGN})
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)
    set(test_input "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i")
    set(test_simulation_log "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.simulation.log")
    set(test_checker_log "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.checker.log")
    add_test(${TEST_NAME} ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DINPUT=${test_input}"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DSIMULATION_LOG=${test_simulation_log}"
        "-DCHECKER_LOG=${test_checker_log}"
        "-DCHECKER=${TWO_STREAM_RADIATION_CHECKER}"
        "-DPLOTFILE=${CURRENT_TEST_BINARY_DIR}/${PLTFILE}"
        "-DRUNTIME_OPTIONS=${ADD_TEST_TSR_RUNTIME_OPTIONS}"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunTwoStreamRadiation.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;radiation"
        ATTACHED_FILES_ON_FAIL "${test_simulation_log};${test_checker_log}")
endfunction(add_test_two_stream_radiation)

function(add_test_cloud_chamber_parity TEST_NAME)
    set(TEST_FILES_DIR "CloudChamber_SatAdj")
    if (ARGC GREATER 1)
        set(TEST_FILES_DIR "${ARGV1}")
    endif()
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)
    add_test(${TEST_NAME} ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DINPUT=${CURRENT_TEST_BINARY_DIR}/${TEST_FILES_DIR}.i"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DCHECKER=${CLOUD_CHAMBER_CHECKER}"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunCloudChamberParity.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;cloud-chamber"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/budget_off/simulation.log;${CURRENT_TEST_BINARY_DIR}/budget_on/simulation.log;${CURRENT_TEST_BINARY_DIR}/parity.log")
endfunction(add_test_cloud_chamber_parity)

# Run the deck in test_files/<TEST_FILES_DIR> on a single box (one rank) and on a split
# BoxArray (ERF_TEST_NRANKS ranks) and compare the two plotfiles PLTFILE with fcompare.
# COMMON_OPTIONS go to both runs, REFERENCE_OPTIONS must make the grid a single box and
# SPLIT_OPTIONS give the split (the deck's own grid when empty).
function(add_test_box_parity TEST_NAME TEST_FILES_DIR PLTFILE)
    set(oneValueArgs "COMMON_OPTIONS" "REFERENCE_OPTIONS" "SPLIT_OPTIONS" "FCOMPARE_RTOL" "FCOMPARE_ATOL" "DATALOG" "DATALOG_SIGDIGITS")
    cmake_parse_arguments(ADD_TEST_BP "" "${oneValueArgs}" "" ${ARGN})
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)

    set(_fcompare_rtol "${ERF_TEST_FCOMPARE_RTOL}")
    set(_fcompare_atol "${ERF_TEST_FCOMPARE_ATOL}")
    if(NOT "${ADD_TEST_BP_FCOMPARE_RTOL}" STREQUAL "")
        set(_fcompare_rtol "${ADD_TEST_BP_FCOMPARE_RTOL}")
    endif()
    if(NOT "${ADD_TEST_BP_FCOMPARE_ATOL}" STREQUAL "")
        set(_fcompare_atol "${ADD_TEST_BP_FCOMPARE_ATOL}")
    endif()

    add_test(${TEST_NAME} ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DINPUT=${CURRENT_TEST_BINARY_DIR}/${TEST_FILES_DIR}.i"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DFCOMPARE=${FCOMPARE_EXE}"
        "-DPLTFILE=${PLTFILE}"
        "-DRTOL=${_fcompare_rtol}"
        "-DATOL=${_fcompare_atol}"
        "-DCOMMON_OPTIONS=${ADD_TEST_BP_COMMON_OPTIONS}"
        "-DREFERENCE_OPTIONS=${ADD_TEST_BP_REFERENCE_OPTIONS}"
        "-DSPLIT_OPTIONS=${ADD_TEST_BP_SPLIT_OPTIONS}"
        "-DDATALOG=${ADD_TEST_BP_DATALOG}"
        "-DDATALOG_SIGDIGITS=${ADD_TEST_BP_DATALOG_SIGDIGITS}"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunBoxParity.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;box-parity"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/one_box/simulation.log;${CURRENT_TEST_BINARY_DIR}/split/simulation.log;${CURRENT_TEST_BINARY_DIR}/parity.log")
endfunction(add_test_box_parity)

# Run the deck in test_files/<TEST_FILES_DIR> twice, once with a diagnostic off and once
# with it on, and require PLTFILE to be identical bit for bit.  This is the harness for
# "turning this output on does not change the answer": OFF_OPTIONS and ON_OPTIONS are the
# two option sets, COMMON_OPTIONS go to both, and REQUIRE_ON_FILE names a file the on leg
# must write and the off leg must not, so a misspelled option cannot pass as agreement.
function(add_test_option_parity TEST_NAME TEST_FILES_DIR PLTFILE)
    set(oneValueArgs "COMMON_OPTIONS" "OFF_OPTIONS" "ON_OPTIONS" "REQUIRE_ON_FILE" "RUN_TIMEOUT")
    cmake_parse_arguments(ADD_TEST_OP "" "${oneValueArgs}" "" ${ARGN})
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)

    set(_run_timeout 600)
    set(_ctest_timeout 600)
    if(DEFINED ADD_TEST_OP_RUN_TIMEOUT)
        set(_run_timeout "${ADD_TEST_OP_RUN_TIMEOUT}")
        math(EXPR _ctest_timeout "2 * ${_run_timeout} + 600")
    endif()

    add_test(${TEST_NAME} ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DINPUT=${CURRENT_TEST_BINARY_DIR}/${TEST_FILES_DIR}.i"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DFCOMPARE=${FCOMPARE_EXE}"
        "-DPLTFILE=${PLTFILE}"
        "-DRUN_TIMEOUT=${_run_timeout}"
        "-DCOMMON_OPTIONS=${ADD_TEST_OP_COMMON_OPTIONS}"
        "-DOFF_OPTIONS=${ADD_TEST_OP_OFF_OPTIONS}"
        "-DON_OPTIONS=${ADD_TEST_OP_ON_OPTIONS}"
        "-DREQUIRE_ON_FILE=${ADD_TEST_OP_REQUIRE_ON_FILE}"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunOptionParity.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT ${_ctest_timeout}
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;option-parity"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/option_off/simulation.log;${CURRENT_TEST_BINARY_DIR}/option_on/simulation.log;${CURRENT_TEST_BINARY_DIR}/parity.log")
endfunction(add_test_option_parity)

# The numeric log comparison add_test_box_parity's DATALOG relies on: a comparator that
# accepts everything passes every test that uses it, so it needs its own test.  Pure CMake,
# no ERF run, hence the "unit" label.
add_test(CompareDataLogs_SelfTest ${CMAKE_COMMAND}
    "-DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/test_files/CompareDataLogs_SelfTest"
    -P ${PROJECT_SOURCE_DIR}/Tests/CompareDataLogsSelfTest.cmake)
set_tests_properties(CompareDataLogs_SelfTest
    PROPERTIES
    TIMEOUT 60
    PROCESSORS 1
    LABELS "unit;box-parity")

# Restart parity: run one deck straight, then to a checkpoint and on from it, and
# require the plotfile at the end to be identical (no gold file). Every run has a
# time limit; the default stays at 600, but an explicit RUN_TIMEOUT is forwarded
# unchanged to each leg and used to size the outer CTest watchdog.
function(add_test_restart_parity TEST_NAME TEST_FILES_DIR STEP_CHK STEP_END)
    set(oneValueArgs "COMMON_OPTIONS" "FCOMPARE_RTOL" "FCOMPARE_ATOL" "RUN_TIMEOUT" "DATALOG" "DATALOG_SIGDIGITS")
    cmake_parse_arguments(ADD_TEST_RP "" "${oneValueArgs}" "" ${ARGN})
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)

    set(_fcompare_rtol "${ERF_TEST_FCOMPARE_RTOL}")
    set(_fcompare_atol "${ERF_TEST_FCOMPARE_ATOL}")
    if(NOT "${ADD_TEST_RP_FCOMPARE_RTOL}" STREQUAL "")
        set(_fcompare_rtol "${ADD_TEST_RP_FCOMPARE_RTOL}")
    endif()
    if(NOT "${ADD_TEST_RP_FCOMPARE_ATOL}" STREQUAL "")
        set(_fcompare_atol "${ADD_TEST_RP_FCOMPARE_ATOL}")
    endif()
    set(_run_timeout 600)
    set(_ctest_timeout 600)
    if(DEFINED ADD_TEST_RP_RUN_TIMEOUT)
        set(_run_timeout "${ADD_TEST_RP_RUN_TIMEOUT}")
        math(EXPR _ctest_timeout "3 * ${_run_timeout} + 600")
    endif()

    add_test(${TEST_NAME} ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DINPUT=${CURRENT_TEST_BINARY_DIR}/${TEST_FILES_DIR}.i"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DFCOMPARE=${FCOMPARE_EXE}"
        "-DSTEP_CHK=${STEP_CHK}"
        "-DSTEP_END=${STEP_END}"
        "-DRTOL=${_fcompare_rtol}"
        "-DATOL=${_fcompare_atol}"
        "-DRUN_TIMEOUT=${_run_timeout}"
        "-DCOMMON_OPTIONS=${ADD_TEST_RP_COMMON_OPTIONS}"
        "-DDATALOG=${ADD_TEST_RP_DATALOG}"
        "-DDATALOG_SIGDIGITS=${ADD_TEST_RP_DATALOG_SIGDIGITS}"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunRestartParity.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT ${_ctest_timeout}
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;restart-parity"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/straight/simulation.log;${CURRENT_TEST_BINARY_DIR}/restart/checkpoint.log;${CURRENT_TEST_BINARY_DIR}/restart/restart.log;${CURRENT_TEST_BINARY_DIR}/parity.log")
endfunction(add_test_restart_parity)

# Tiling parity: run one deck with MFIter tiling on and off and require identical
# 3D and 2D plotfiles (no gold file). Catches kernels that loop over the valid box
# while indexing per-tile work arrays. VARYING_3D / VARYING_2D list fields (space
# separated) that must take more than one value in the untiled run, so the
# agreement is not between two copies of a constant.
function(add_test_tiling_parity TEST_NAME TEST_FILES_DIR PLTFILE PLT2DFILE)
    set(options )
    set(oneValueArgs "RUNTIME_OPTIONS" "VARYING_3D" "VARYING_2D")
    set(multiValueArgs )
    cmake_parse_arguments(ADD_TEST_TP "${options}" "${oneValueArgs}"
        "${multiValueArgs}" ${ARGN})

    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)
    add_test(${TEST_NAME} ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DINPUT=${CURRENT_TEST_BINARY_DIR}/${TEST_FILES_DIR}.i"
        "-DRUNTIME_OPTIONS=${ADD_TEST_TP_RUNTIME_OPTIONS}"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DFCOMPARE=${FCOMPARE_EXE}"
        "-DFEXTREMA=${FEXTREMA_EXE}"
        "-DRTOL=${ERF_TEST_FCOMPARE_RTOL}"
        "-DATOL=${ERF_TEST_FCOMPARE_ATOL}"
        "-DPLTFILE=${PLTFILE}"
        "-DPLT2DFILE=${PLT2DFILE}"
        "-DVARYING_3D=${ADD_TEST_TP_VARYING_3D}"
        "-DVARYING_2D=${ADD_TEST_TP_VARYING_2D}"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunTilingParity.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/tiled.log;${CURRENT_TEST_BINARY_DIR}/untiled.log;${CURRENT_TEST_BINARY_DIR}/fcompare_plt.log;${CURRENT_TEST_BINARY_DIR}/fcompare_plt2d.log;${CURRENT_TEST_BINARY_DIR}/fextrema_plt.log;${CURRENT_TEST_BINARY_DIR}/fextrema_plt2d.log")
endfunction(add_test_tiling_parity)

function(add_test_cloud_chamber_budget TEST_NAME MODE SOURCE_NAME)
    set(_cloud_chamber_input_name "${SOURCE_NAME}")
    set(TEST_FILES_DIR "${SOURCE_NAME}")
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)
    add_test(${TEST_NAME} ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DINPUT=${CURRENT_TEST_BINARY_DIR}/${_cloud_chamber_input_name}.i"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DCHECKER=${CLOUD_CHAMBER_CHECKER}"
        "-DMODE=${MODE}"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunCloudChamberBudget.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;cloud-chamber"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/simulation.log;${CURRENT_TEST_BINARY_DIR}/checker.log;${CURRENT_TEST_BINARY_DIR}/cloud_chamber_budget.dat")
endfunction(add_test_cloud_chamber_budget)

# At-rest test: a hydrostatic atmosphere over terrain must stay at rest with lateral
# outflow boundaries, where the mesh is extrapolated past the domain and the base state in
# the ghost cells has to be built at the height the mesh puts them at rather than copied.
function(add_test_at_rest_terrain_outflow TEST_NAME PLTFILE TOLERANCE GRADP_TOLERANCE)
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)
    add_test(${TEST_NAME} ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DINPUT=${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DFEXTREMA=${FEXTREMA_EXE}"
        "-DPLTFILE=${PLTFILE}"
        "-DTOLERANCE=${TOLERANCE}"
        "-DGRADP_TOLERANCE=${GRADP_TOLERANCE}"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunAtRestTerrainOutflow.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/symmetry/simulation.log;${CURRENT_TEST_BINARY_DIR}/outflow/simulation.log;${CURRENT_TEST_BINARY_DIR}/at_rest.log")
endfunction(add_test_at_rest_terrain_outflow)

# Positive startup regression for the retained legacy theta/qv parser path.
# This intentionally has no physical-temperature or physical-wall keys.
function(add_test_cloud_chamber_legacy_config TEST_NAME)
    set(TEST_FILES_DIR "CloudChamber_Legacy_Config")
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)
    set(test_input "${CURRENT_TEST_BINARY_DIR}/CloudChamber_Legacy_Config.i")
    set(test_log "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log")
    set(output_directory "${CURRENT_TEST_BINARY_DIR}/legacy_plt00000")
    set(output_artifact "${output_directory}/Header")
    add_test(${TEST_NAME} ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DINPUT=${test_input}"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DLOG=${test_log}"
        "-DOUTPUT_DIRECTORY=${output_directory}"
        "-DOUTPUT_ARTIFACT=${output_artifact}"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunCloudChamberConfigSuccess.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 180
        PROCESSORS 1
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;cloud-chamber;configuration"
        ATTACHED_FILES_ON_FAIL "${test_log};${output_artifact}")
endfunction(add_test_cloud_chamber_legacy_config)
# Negative startup tests for the Native SHOC transport modes removed from the
# production input contract.  The shared fixture supplies a complete Native
# SHOC run, while the runtime option exercises the real ParmParse reader path.
function(add_test_shoc_removed_transport TEST_NAME RUNTIME_OPTION EXPECTED_MESSAGE
        EXPECTED_GUIDANCE_1 EXPECTED_GUIDANCE_2)
    set(TEST_FILES_DIR "SHOC_Stable_Clear")
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)

    set(test_input "${CURRENT_TEST_BINARY_DIR}/SHOC_Stable_Clear.i")
    set(test_log "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log")
    add_test(${TEST_NAME} ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DINPUT=${test_input}"
        "-DRUNTIME_OPTIONS=${RUNTIME_OPTION}"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DLOG=${test_log}"
        "-DEXPECTED_MESSAGE=${EXPECTED_MESSAGE}"
        "-DEXPECTED_GUIDANCE_1=${EXPECTED_GUIDANCE_1}"
        "-DEXPECTED_GUIDANCE_2=${EXPECTED_GUIDANCE_2}"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunShocRemovedTransportConfig.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 180
        PROCESSORS 1
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;shoc;configuration"
        ATTACHED_FILES_ON_FAIL "${test_log}")
endfunction(add_test_shoc_removed_transport)

add_test_shoc_removed_transport(SHOC_Removed_Scalar_Host_Diffusion
    "erf.shoc.transport_mode=host_diffusion"
    "erf.shoc.transport_mode = host_diffusion has been removed for native SHOC"
    "Use erf.shoc.transport_mode = state_update"
    "")
add_test_shoc_removed_transport(SHOC_Removed_Momentum_Host_Diffusion
    "erf.shoc.momentum_transport=host_diffusion"
    "erf.shoc.momentum_transport = host_diffusion has been removed for native SHOC"
    "state_update"
    "none")

# Production wiring regression: two dry runs differ only in xlo roughness;
# the checker requires finite output and a resolvable z0_m response.

function(add_test_cloud_chamber_neutral_momentum TEST_NAME)
    set(TEST_FILES_DIR "CloudChamber_Dry_NeutralMomentum")
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)
    add_test(${TEST_NAME} ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DINPUT=${CURRENT_TEST_BINARY_DIR}/CloudChamber_Dry_NeutralMomentum.i"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DCHECKER=${CLOUD_CHAMBER_CHECKER}"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunCloudChamberNeutralMomentum.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;cloud-chamber;neutral-roughness"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/z0_baseline/simulation.log;${CURRENT_TEST_BINARY_DIR}/z0_changed/simulation.log;${CURRENT_TEST_BINARY_DIR}/neutral_momentum_checker.log")
endfunction(add_test_cloud_chamber_neutral_momentum)

# Production wiring regression for fixed bulk aerodynamic momentum.  The
# harness changes only C_D and requires a measurable velocity response.
function(add_test_cloud_chamber_fixed_momentum TEST_NAME)
    set(TEST_FILES_DIR "CloudChamber_Dry_FixedMomentum")
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)
    add_test(${TEST_NAME} ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DINPUT=${CURRENT_TEST_BINARY_DIR}/CloudChamber_Dry_FixedMomentum.i"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DCHECKER=${CLOUD_CHAMBER_CHECKER}"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunCloudChamberFixedMomentum.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;cloud-chamber;bulk-momentum"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/cd_baseline/simulation.log;${CURRENT_TEST_BINARY_DIR}/cd_changed/simulation.log;${CURRENT_TEST_BINARY_DIR}/fixed_momentum_checker.log")
endfunction(add_test_cloud_chamber_fixed_momentum)

# Production wiring regression for all-channel horizontal MOST.  The harness
# checks wet-wall budgets in both runs and changes only horizontal momentum
# transfer to require an observable production-path momentum response.
function(add_test_cloud_chamber_most TEST_NAME)
    set(TEST_FILES_DIR "CloudChamber_SatAdj_MOSTMixedWalls")
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)
    add_test(${TEST_NAME} ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DINPUT=${CURRENT_TEST_BINARY_DIR}/CloudChamber_SatAdj_MOSTMixedWalls.i"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DCHECKER=${CLOUD_CHAMBER_CHECKER}"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunCloudChamberMOST.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;cloud-chamber;most"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/most_baseline/simulation.log;${CURRENT_TEST_BINARY_DIR}/most_changed/simulation.log;${CURRENT_TEST_BINARY_DIR}/most_momentum_checker.log")
endfunction(add_test_cloud_chamber_most)

function(add_test_cloud_chamber_fixed_dt_guard TEST_NAME)
    set(test_log "${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.log")
    add_test(NAME ${TEST_NAME} COMMAND ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DTEST_EXE=$<TARGET_FILE:erf_cloud_chamber_wall_dt_guard_check>"
        "-DLOG=${test_log}"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunCloudChamberWallDtGuardFailure.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 120
        PROCESSORS 1
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/"
        LABELS "regression;cloud-chamber;configuration"
        ATTACHED_FILES_ON_FAIL "${test_log}")
endfunction(add_test_cloud_chamber_fixed_dt_guard)

function(add_test_cloud_chamber_openmp TEST_NAME)
    set(TEST_FILES_DIR "CloudChamber_SatAdj")
    if (ARGC GREATER 1)
        set(TEST_FILES_DIR "${ARGV1}")
    endif()
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)
    add_test(${TEST_NAME} ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DINPUT=${CURRENT_TEST_BINARY_DIR}/${TEST_FILES_DIR}.i"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DCHECKER=${CLOUD_CHAMBER_CHECKER}"
        "-DCMAKE_COMMAND=${CMAKE_COMMAND}"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunCloudChamberOpenMP.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;cloud-chamber;openmp"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/omp_1_thread/simulation.log;${CURRENT_TEST_BINARY_DIR}/omp_2_threads/simulation.log;${CURRENT_TEST_BINARY_DIR}/openmp_parity.log")
endfunction(add_test_cloud_chamber_openmp)

# Native SHOC regression test.  This intentionally remains separate from
# add_test_r so existing registrations retain their exact command and
# fixture behaviour.  TEST_FILES_DIR and INPUT_FILE allow the small SHOC
# matrix to reuse a physical fixture while keeping unique binary and gold
# directories.  The checker runs before the selected gold comparison path.
function(add_test_shoc_r TEST_NAME TEST_DIR TEST_EXE PLTFILE)
    set(options SKIP_GOLD)
    set(oneValueArgs "TEST_FILES_DIR" "INPUT_FILE" "CHECK_MODE" "RUNTIME_OPTIONS" "TIMEOUT"
        "GOLD_COMPARISON" "GOLD_MODE")
    set(multiValueArgs "LABELS")
    cmake_parse_arguments(ADD_TEST_SHOC_R "${options}" "${oneValueArgs}"
        "${multiValueArgs}" ${ARGN})

    set(TEST_FILES_DIR "${ADD_TEST_SHOC_R_TEST_FILES_DIR}")
    setup_test()

    if("${ADD_TEST_SHOC_R_INPUT_FILE}" STREQUAL "")
        set(test_input "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i")
    else()
        set(test_input "${CURRENT_TEST_BINARY_DIR}/${ADD_TEST_SHOC_R_INPUT_FILE}")
    endif()

    set(RUNTIME_OPTIONS "${ADD_TEST_SHOC_R_RUNTIME_OPTIONS}")
    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)

    set(test_log "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log")
    if(ADD_TEST_SHOC_R_GOLD_COMPARISON)
        set(_shoc_gold_comparison "${ADD_TEST_SHOC_R_GOLD_COMPARISON}")
    else()
        set(_shoc_gold_comparison "fcompare")
    endif()
    add_test(${TEST_NAME} ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DINPUT=${test_input}"
        "-DRUNTIME_OPTIONS=${RUNTIME_OPTIONS}"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DLOG=${test_log}"
        "-DCHECKER=${SHOC_PLOTFILE_CHECKER}"
        "-DCHECK_MODE=${ADD_TEST_SHOC_R_CHECK_MODE}"
        "-DINITIAL=${CURRENT_TEST_BINARY_DIR}/plt00000"
        "-DMIDPOINT=${CURRENT_TEST_BINARY_DIR}/plt00010"
        "-DFINAL=${CURRENT_TEST_BINARY_DIR}/${PLTFILE}"
        "-DFCOMPARE=${FCOMPARE_EXE}"
        "-DGOLD_DIFFERENTIAL=${SHOC_GOLD_DIFFERENTIAL}"
        "-DGOLD_COMPARISON=${_shoc_gold_comparison}"
        "-DGOLD_MODE=${ADD_TEST_SHOC_R_GOLD_MODE}"
        "-DRTOL=${ERF_TEST_FCOMPARE_RTOL}"
        "-DATOL=${ERF_TEST_FCOMPARE_ATOL}"
        "-DGOLD=${PLOT_GOLD}"
        "-DSKIP_GOLD=${ADD_TEST_SHOC_R_SKIP_GOLD}"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunShocRegression.cmake)
    if(DEFINED ADD_TEST_SHOC_R_TIMEOUT)
        set(_shoc_timeout "${ADD_TEST_SHOC_R_TIMEOUT}")
    else()
        set(_shoc_timeout 600)
    endif()
    if(ADD_TEST_SHOC_R_LABELS)
        set(_shoc_labels "${ADD_TEST_SHOC_R_LABELS}")
    else()
        set(_shoc_labels "regression;shoc")
    endif()
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT ${_shoc_timeout}
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "${_shoc_labels}"
        ATTACHED_FILES_ON_FAIL "${test_log}")
endfunction(add_test_shoc_r)

function(add_test_shoc_mutation TEST_NAME MUTATION_OPTION TARGET_FIELD
        MIN_FINAL_DIFFERENCE MIN_BASELINE_EVOLUTION MAX_MUTANT_TO_BASELINE_EVOLUTION_RATIO)
    set(TEST_FILES_DIR "SHOC_Stable_Clear")
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)
    set(_baseline_dir "${CURRENT_TEST_BINARY_DIR}/baseline")
    set(_mutant_dir "${CURRENT_TEST_BINARY_DIR}/mutant")
    file(MAKE_DIRECTORY "${_baseline_dir}" "${_mutant_dir}")
    file(COPY "${CURRENT_TEST_SOURCE_DIR}/." DESTINATION "${_baseline_dir}")
    file(COPY "${CURRENT_TEST_SOURCE_DIR}/." DESTINATION "${_mutant_dir}")

    set(_baseline_log "${_baseline_dir}/${TEST_NAME}_baseline.log")
    set(_mutant_log "${_mutant_dir}/${TEST_NAME}_mutant.log")
    add_test(${TEST_NAME} ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DBASELINE_INPUT=${_baseline_dir}/SHOC_Stable_Clear.i"
        "-DMUTANT_INPUT=${_mutant_dir}/SHOC_Stable_Clear.i"
        "-DBASELINE_OPTIONS="
        "-DMUTANT_OPTIONS=${MUTATION_OPTION}"
        "-DBASELINE_WORKING_DIRECTORY=${_baseline_dir}"
        "-DMUTANT_WORKING_DIRECTORY=${_mutant_dir}"
        "-DBASELINE_LOG=${_baseline_log}"
        "-DMUTANT_LOG=${_mutant_log}"
        "-DCHECKER=${SHOC_MUTATION_DIFFERENTIAL}"
        "-DTARGET_FIELD=${TARGET_FIELD}"
        "-DMIN_FINAL_DIFFERENCE=${MIN_FINAL_DIFFERENCE}"
        "-DMIN_BASELINE_EVOLUTION=${MIN_BASELINE_EVOLUTION}"
        "-DMAX_MUTANT_TO_BASELINE_EVOLUTION_RATIO=${MAX_MUTANT_TO_BASELINE_EVOLUTION_RATIO}"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunShocMutationRegression.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;shoc;mutation"
        ATTACHED_FILES_ON_FAIL "${_baseline_log};${_mutant_log}")
endfunction(add_test_shoc_mutation)

# Negative-control regression: the true fixed_dt > dt_wall guard must fire
# and expose stable diagnostic fields for automated CI forensics.
add_test_cloud_chamber_fixed_dt_guard(CloudChamber_Bulk_FixedDtGuard)

if(ERF_ENABLE_MPI)
add_test_anelastic_wall_diffusion(AnelasticWallDiffusion_X 0)
add_test_anelastic_wall_diffusion(AnelasticWallDiffusion_Y 1)
add_test_anelastic_wall_diffusion(AnelasticWallDiffusion_Z 2)
# Same stationary state as the _X case, but with erf.anelastic_type = MidPoint so the
# vertical implicit diffusion stays on (the _X/_Y/_Z cases opt out with vert_implicit).
add_test_anelastic_wall_diffusion(AnelasticWallDiffusion_X_MidPoint 0)
add_test_cloud_chamber(CloudChamber_Dry dry)
add_test_cloud_chamber_legacy_config(CloudChamber_Legacy_Config)
add_test_cloud_chamber_neutral_momentum(CloudChamber_Dry_NeutralMomentumActivation)
add_test_cloud_chamber_fixed_momentum(CloudChamber_Dry_FixedMomentumActivation)
add_test_cloud_chamber(CloudChamber_SatAdj cloudy)
add_test_cloud_chamber_parity(CloudChamber_SatAdj_Parity)
add_test_cloud_chamber_budget(CloudChamber_SatAdj_AllDry all_dry CloudChamber_SatAdj_AllDry)
add_test_cloud_chamber_budget(CloudChamber_SatAdj_WetBudget wet_budget CloudChamber_SatAdj_WetBudget)
add_test_cloud_chamber_budget(CloudChamber_SatAdj_BulkMixedWet bulk_wet CloudChamber_SatAdj_BulkMixedWet)
add_test_cloud_chamber_budget(CloudChamber_SatAdj_NeutralWetBudget neutral_wet CloudChamber_SatAdj_NeutralWet)
add_test_cloud_chamber_budget(CloudChamber_SatAdj_MOSTWetBudget most_wet CloudChamber_SatAdj_MOSTWetBudget)
add_test_cloud_chamber_most(CloudChamber_SatAdj_MOSTMixedWalls)
if(ERF_ENABLE_OPENMP)
add_test_cloud_chamber_openmp(CloudChamber_SatAdj_OpenMP)
endif()
# ctest runs this argv directly, with no shell, so the launcher cannot be
# pasted in as one word: "flux run" has to reach ctest as two arguments.
erf_mpi_launcher_command(SHOC_DIFFERENTIAL_LAUNCHER
    LAUNCHER "${MPIEXEC_EXECUTABLE}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "SHOC microphysics differential test"
    OPTIONAL)
add_test(SHOC_Unstable_Cloud_SatAdj_vs_NoCond
    ${SHOC_DIFFERENTIAL_LAUNCHER}
    ${SHOC_MICROPHYSICS_DIFFERENTIAL}
    ${CMAKE_CURRENT_BINARY_DIR}/test_files/SHOC_Unstable_Cloud_SatAdj_Property/plt00020
    ${CMAKE_CURRENT_BINARY_DIR}/test_files/SHOC_Unstable_Cloud_NoCond_Property/plt00020)
set_tests_properties(SHOC_Unstable_Cloud_SatAdj_vs_NoCond
    PROPERTIES
    DEPENDS "SHOC_Unstable_Cloud_SatAdj_Property;SHOC_Unstable_Cloud_NoCond_Property"
    TIMEOUT 120
    PROCESSORS 1
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/"
    LABELS "regression;shoc;microphysics")
# execute_process needs mpiexec, and does not expand the executable globs used on Windows
if(NOT WIN32)
add_test_at_rest_terrain_outflow(AtRestTerrainOutflow "plt00400" 1.0e-8 0.1)
endif()
endif()

# Debug regression test with lower tolerance
function(add_test_d TEST_NAME TEST_DIR TEST_EXE PLTFILE)
    set(options )
    set(oneValueArgs "INPUT_SOUNDING" "RUNTIME_OPTIONS")
    set(multiValueArgs )
    cmake_parse_arguments(ADD_TEST_D "${options}" "${oneValueArgs}"
        "${multiValueArgs}" ${ARGN})

    setup_test()

    set(RUNTIME_OPTIONS "${ADD_TEST_D_RUNTIME_OPTIONS}")
    if(NOT "${ADD_TEST_D_INPUT_SOUNDING}" STREQUAL "")
      string(APPEND RUNTIME_OPTIONS "erf.input_sounding_file=${CURRENT_TEST_BINARY_DIR}/${ADD_TEST_D_INPUT_SOUNDING}")
    endif()

    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)
    set(FCOMPARE_TOLERANCE "-r 3.0e-9 --abs_tol 3.0e-9")
    set(FCOMPARE_FLAGS "--abort_if_not_all_found -a ${FCOMPARE_TOLERANCE}")
    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUNTIME_OPTIONS} > ${TEST_NAME}.log && ${MPI_FCOMP_COMMANDS} ${FCOMPARE_EXE} ${FCOMPARE_FLAGS} ${PLOT_GOLD} ${CURRENT_TEST_BINARY_DIR}/${PLTFILE}")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log"
    )
endfunction(add_test_d)

# Stationary test -- compare with time 0
function(add_test_0 TEST_NAME TEST_DIR TEST_EXE PLTFILE)
    set(options )
    set(oneValueArgs "INPUT_SOUNDING" "RUNTIME_OPTIONS")
    set(multiValueArgs )
    cmake_parse_arguments(ADD_TEST_0 "${options}" "${oneValueArgs}"
        "${multiValueArgs}" ${ARGN})

    setup_test()

    set(RUNTIME_OPTIONS "${ADD_TEST_0_RUNTIME_OPTIONS}")
    if(NOT "${ADD_TEST_0_INPUT_SOUNDING}" STREQUAL "")
      string(APPEND RUNTIME_OPTIONS "erf.input_sounding_file=${CURRENT_TEST_BINARY_DIR}/${ADD_TEST_0_INPUT_SOUNDING}")
    endif()

    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)
    set(FCOMPARE_TOLERANCE "-r 1e-14 --abs_tol 1.0e-14")
    set(FCOMPARE_FLAGS "-a ${FCOMPARE_TOLERANCE}")
    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i erf.input_sounding_file=${CURRENT_TEST_BINARY_DIR}/input_sounding ${RUNTIME_OPTIONS} > ${TEST_NAME}.log && ${MPI_FCOMP_COMMANDS} ${FCOMPARE_EXE} ${FCOMPARE_FLAGS} ${CURRENT_TEST_BINARY_DIR}/plt00000 ${CURRENT_TEST_BINARY_DIR}/${PLTFILE}")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log"
    )
endfunction(add_test_0)

# SDM regression test
function(add_test_sdm TEST_NAME TEST_DIR TEST_EXE PLTFILE TEST_RTOL TEST_ATOL)
    set(options )
    set(oneValueArgs "INPUT_SOUNDING" "RUNTIME_OPTIONS")
    set(multiValueArgs )
    cmake_parse_arguments(ADD_TEST_SDM "${options}" "${oneValueArgs}"
        "${multiValueArgs}" ${ARGN})

    setup_test()

    set(RUNTIME_OPTIONS "${ADD_TEST_SDM_RUNTIME_OPTIONS}")
    if(NOT "${ADD_TEST_SDM_INPUT_SOUNDING}" STREQUAL "")
      string(APPEND RUNTIME_OPTIONS "erf.input_sounding_file=${CURRENT_TEST_BINARY_DIR}/${ADD_TEST_SDM_INPUT_SOUNDING}")
    endif()

    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)

    if(ERF_SDM_SMOKE_ONLY)
        # No gold file is available for this case here, so run it to completion
        # and let assertions and aborts be the check.
        set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUNTIME_OPTIONS} > ${TEST_NAME}.log")
        set(TEST_LABELS "smoke")
    else()
        set(FCOMPARE_TOLERANCE "--rel_tol ${TEST_RTOL} --abs_tol ${TEST_ATOL}")
        set(FCOMPARE_FLAGS "--abort_if_not_all_found --allow_diff_grids ${FCOMPARE_TOLERANCE}")
        set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUNTIME_OPTIONS} > ${TEST_NAME}.log && ${MPI_FCOMP_COMMANDS} ${FCOMPARE_EXE} ${FCOMPARE_FLAGS} ${PLOT_GOLD} ${CURRENT_TEST_BINARY_DIR}/${PLTFILE}")
        set(TEST_LABELS "regression")
    endif()

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "${TEST_LABELS}"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log"
    )
endfunction(add_test_sdm)

#=============================================================================
# Regression tests
#=============================================================================

if(ERF_ENABLE_TESTS AND ERF_ENABLE_MPI)
    # The checker is a small AMReX PlotFileData consumer and is built only
    # when regression tests are enabled.  All SHOC cases use explicit
    # state-update ownership settings in their input fixture.
    # DOUBLE keeps the historical strict fcompare oracle.  SINGLE uses the
    # field-aware comparator so the shared gold remains the scientific
    # reference while represented-float roundoff is bounded per field.
    if(ERF_PRECISION STREQUAL "SINGLE")
        set(_shoc_clear_gold_comparison "field_aware")
    else()
        set(_shoc_clear_gold_comparison "fcompare")
    endif()

    add_test_shoc_r(SHOC_Stable_Clear "" "erf_exec" "plt00020"
        TEST_FILES_DIR "SHOC_Stable_Clear"
        CHECK_MODE "stable_clear"
        GOLD_COMPARISON "${_shoc_clear_gold_comparison}"
        GOLD_MODE "stable_clear"
        LABELS regression shoc
        TIMEOUT 600)
    add_test_shoc_r(SHOC_Stable_Cloud "" "erf_exec" "plt00020"
        TEST_FILES_DIR "SHOC_Stable_Cloud"
        CHECK_MODE "stable_cloud"
        GOLD_COMPARISON "field_aware"
        GOLD_MODE "stable_cloud"
        LABELS regression shoc
        TIMEOUT 600)
    add_test_shoc_r(SHOC_Unstable_Clear_BOMEX "" "erf_exec" "plt00020"
        TEST_FILES_DIR "SHOC_Unstable_Clear_BOMEX"
        CHECK_MODE "unstable_clear"
        GOLD_COMPARISON "${_shoc_clear_gold_comparison}"
        GOLD_MODE "unstable_clear"
        LABELS regression shoc
        TIMEOUT 600)
    add_test_shoc_r(SHOC_Unstable_Cloud_SatAdj "" "erf_exec" "plt00020"
        TEST_FILES_DIR "SHOC_Unstable_Cloud"
        INPUT_FILE "SHOC_Unstable_Cloud.i"
        CHECK_MODE "unstable_cloud"
        GOLD_COMPARISON "field_aware"
        GOLD_MODE "unstable_cloud"
        RUNTIME_OPTIONS "erf.moisture_model=SatAdj erf.buoyancy_type=1 "
        LABELS regression shoc microphysics
        TIMEOUT 600)
    add_test_shoc_r(SHOC_Unstable_Cloud_NoCond "" "erf_exec" "plt00020"
        TEST_FILES_DIR "SHOC_Unstable_Cloud"
        INPUT_FILE "SHOC_Unstable_Cloud.i"
        CHECK_MODE "unstable_cloud_nocond"
        GOLD_COMPARISON "field_aware"
        GOLD_MODE "unstable_cloud_nocond"
        RUNTIME_OPTIONS "erf.moisture_model=MoistNoCondensation erf.buoyancy_type=1 "
        LABELS regression shoc microphysics
        TIMEOUT 600)
    add_test_shoc_r(SHOC_Unstable_Cloud_SatAdj_Property "" "erf_exec" "plt00020"
        TEST_FILES_DIR "SHOC_Unstable_Cloud"
        INPUT_FILE "SHOC_Unstable_Cloud.i"
        CHECK_MODE "unstable_cloud"
        GOLD_COMPARISON "field_aware"
        GOLD_MODE "unstable_cloud"
        RUNTIME_OPTIONS "erf.moisture_model=SatAdj erf.buoyancy_type=1 "
        SKIP_GOLD
        LABELS regression shoc microphysics property
        TIMEOUT 600)
    add_test_shoc_r(SHOC_Unstable_Cloud_NoCond_Property "" "erf_exec" "plt00020"
        TEST_FILES_DIR "SHOC_Unstable_Cloud"
        INPUT_FILE "SHOC_Unstable_Cloud.i"
        CHECK_MODE "unstable_cloud_nocond"
        GOLD_COMPARISON "field_aware"
        GOLD_MODE "unstable_cloud_nocond"
        RUNTIME_OPTIONS "erf.moisture_model=MoistNoCondensation erf.buoyancy_type=1 "
        SKIP_GOLD
        LABELS regression shoc microphysics property
        TIMEOUT 600)
    add_test_shoc_r(SHOC_Unstable_Cloud_Kessler "" "erf_exec" "plt00020"
        TEST_FILES_DIR "SHOC_Unstable_Cloud"
        INPUT_FILE "SHOC_Unstable_Cloud_Kessler.i"
        CHECK_MODE "unstable_cloud_kessler"
        GOLD_COMPARISON "field_aware"
        GOLD_MODE "unstable_cloud_kessler"
        LABELS regression shoc microphysics
        TIMEOUT 600)
    add_test_shoc_r(SHOC_Unstable_Cloud_WSM6 "" "erf_exec" "plt00020"
        TEST_FILES_DIR "SHOC_Unstable_Cloud"
        INPUT_FILE "SHOC_Unstable_Cloud_WSM6.i"
        CHECK_MODE "unstable_cloud_wsm6"
        GOLD_COMPARISON "field_aware"
        GOLD_MODE "unstable_cloud_wsm6"
        LABELS regression shoc microphysics
        TIMEOUT 600)
    add_test_shoc_mutation(SHOC_Mutation_Disable_Tke_State_Update
        "erf.shoc.debug_disable_tke_state_update=true" rhoKE
        1.0e-3 1.0e-3 0.10)
    add_test_shoc_mutation(SHOC_Mutation_Disable_Theta_State_Update
        "erf.shoc.debug_disable_theta_state_update=true" theta
        1.0e-2 1.0e-2 0.05)
endif()

# These tests will all be built in Exec
add_test_plotfile_header(Plotfile3D_DryUnavailableSelection "" "erf_exec" "plt00000")
add_test_r(DensityCurrent                    ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(DensityCurrent_anelastic          ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(DensityCurrent_detJ2              ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(DensityCurrent_detJ2_nosub        ""  "erf_exec" "plt00020" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(DensityCurrent_detJ2_MT           ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(EkmanSpiral                       ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(IsentropicVortexStationary        ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(IsentropicVortexAdvecting         ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(IVA_NumDiff                       ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(MovingTerrain_nosub               ""  "erf_exec" "plt00020" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(MovingTerrain_sub                 ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(Terrain2Lev_STF_interp            ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(Terrain2Lev_STF_transform         ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(RayleighDamping                   ""  "erf_exec" "plt00100" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(ScalarAdvectionUniformU           ""  "erf_exec" "plt00020" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(ScalarAdvectionShearedU           ""  "erf_exec" "plt00080" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(ScalarAdvDiff_order2              ""  "erf_exec" "plt00020" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(ScalarAdvDiff_order3              ""  "erf_exec" "plt00020" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(ScalarAdvDiff_order4              ""  "erf_exec" "plt00020" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(ScalarAdvDiff_order5              ""  "erf_exec" "plt00020" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(ScalarAdvDiff_order6              ""  "erf_exec" "plt00020" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(ScalarAdvDiff_weno3               ""  "erf_exec" "plt00020" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_d(ScalarAdvDiff_weno3z              ""  "erf_exec" "plt00020" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(ScalarAdvDiff_weno5               ""  "erf_exec" "plt00020" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_d(ScalarAdvDiff_weno5z              ""  "erf_exec" "plt00020" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(ScalarAdvDiff_wenomzq3            ""  "erf_exec" "plt00020" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(ScalarDiffusionGaussian           ""  "erf_exec" "plt00020" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(ScalarDiffusionSine               ""  "erf_exec" "plt00020" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(TaylorGreenAdvecting              ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(TaylorGreenAdvectingDiffusing     ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(MSF_NoSub_IsentropicVortexAdv     ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(MSF_Sub_IsentropicVortexAdv       ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
#add_test_r(FlowInABox                       ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(ABL_MOST                          ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(ABL_MOST_IMP_DIFF                 ""  "erf_exec" "plt00010")
add_test_r(ABL_MOST_IMP_DIFF_WOA             ""  "erf_exec" "plt00010")
add_test_r(ABL_MOST_IMP_DIFF_TKE
    ""
    "erf_exec"
    "plt00010"
    FCOMPARE_ATOL "4.0e-10")
if(ERF_ENABLE_FFT)
    add_test_r(ABL_MOST_Cloudchamber         ""  "erf_exec" "plt00010")
endif()
add_test_r(ABL_MOST_SFC                      ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(ABL_MOST_SST                      ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(ABL_MYNN_PBL                      ""  "erf_exec" "plt00100" INPUT_SOUNDING "input_sounding_GABLS1" RUNTIME_OPTIONS "erf.vert_implicit=false " )
# RunTilingParity.cmake calls mpiexec and fcompare through execute_process,
# which neither drops an empty MPIEXEC nor expands the Windows exe globs.
if(ERF_ENABLE_MPI AND NOT WIN32)
  # Box, rank and tiling parity of the other closures: one box on one rank without tiling
  # against four boxes on two ranks with 8 x 8 tiles (Tests/test_files/Closure_BoxParity).
  set(_cbp_ref   "amr.max_grid_size_x=32 amr.max_grid_size_y=32 fabarray.mfiter_tile_size=1024 1024 1024")
  set(_cbp_split "fabarray.mfiter_tile_size=8 8 1024")
  set(_cbp_dir   "${CMAKE_CURRENT_BINARY_DIR}/test_files")
  set(_cbp_ke    "erf.plot_vars_1=density x_velocity y_velocity z_velocity theta KE Kmv Khv")
  foreach(_cbp IN ITEMS
      # erf.vert_implicit defaults to true, so the explicit half of the pair is the one that
      # has to be named: without erf.vert_implicit=false the two Deardorff entries are the
      # same run twice.
      "Deardorff_Explicit|erf.les_type=Deardorff erf.vert_implicit=false|sounding_dry"
      "kEqn_PBLHCap|erf.rans_type=kEqn erf.theta_ref=300 erf.most.pblh_calc=MYNN25 erf.rans_lscale_from_pblh=true|sounding_dry"
      "MYNN25|erf.pbl_type=MYNN25 erf.most.pblh_calc=MYNN25|sounding_dry"
      "MYNNEDMF|erf.pbl_type=MYNNEDMF erf.most.pblh_calc=MYNN25|sounding_dry"
      "MYJ|erf.pbl_type=MYJ erf.most.pblh_calc=MYNN25|sounding_dry"
      "NativeSHOC|erf.pbl_type=NATIVE_SHOC|sounding_dry"
      # sounding_moist is supersaturated below 150 m on purpose: over the 10 steps of a parity
      # run that is what makes Kessler condense, autoconvert and sediment, so qc and qp are
      # compared as varying fields rather than as zero against zero.  Sedimentation in
      # particular sets its substep count from a ParReduce over the whole MultiFab, which is
      # exactly the kind of global quantity a decomposition can get wrong.
      "Kessler_Smagorinsky|erf.les_type=Smagorinsky erf.Cs=0.1 erf.moisture_model=Kessler erf.plot_vars_1=density x_velocity y_velocity z_velocity theta qv qc qp Kmv|sounding_moist"
      # 32 levels starting at 19.503737114205 m and growing by 1.03 sum to the deck's
      # prob_extent[2] = 1024 m, so the stretched mesh reaches the top of the domain and the
      # geometry-derived quantities (Rayleigh damping, MOST reference height) agree with the
      # levels.  erf.initial_dz=20 overshoots to 1050.1 m and ERF prints a mismatch note.
      "Stretched_Smagorinsky|erf.les_type=Smagorinsky erf.Cs=0.1 erf.grid_stretching_ratio=1.03 erf.initial_dz=19.503737114205|sounding_dry"
      "Deardorff_Implicit|erf.les_type=Deardorff|sounding_dry")
    string(REPLACE "|" ";" _cbp_fields "${_cbp}")
    list(GET _cbp_fields 0 _cbp_name)
    list(GET _cbp_fields 1 _cbp_opts)
    list(GET _cbp_fields 2 _cbp_snd)
    set(_cbp_test "Closure_BoxParity_${_cbp_name}")
    add_test_box_parity(${_cbp_test} Closure_BoxParity "plt00010"
        COMMON_OPTIONS "${_cbp_ke} ${_cbp_opts} erf.input_sounding_file=${_cbp_dir}/${_cbp_test}/${_cbp_snd}"
        REFERENCE_OPTIONS "${_cbp_ref}"
        SPLIT_OPTIONS "${_cbp_split}"
        FCOMPARE_RTOL "1.0e-9")
  endforeach()
  if(ERF_ENABLE_FFT)
    # The anelastic MidPoint integrator; its MLMG projection does not converge on this deck,
    # so the FFT solver is used and the tests need the FFT build.  The slow scalars (k, qv)
    # are advected with the projected momentum, which used to be copied tile by tile.
    add_test_box_parity(Closure_BoxParity_Anelastic_Kessler Closure_BoxParity "plt00010"
        COMMON_OPTIONS "erf.plot_vars_1=density x_velocity y_velocity z_velocity theta qv qc qp Kmv erf.anelastic=1 erf.anelastic_type=MidPoint erf.use_fft=true erf.les_type=Smagorinsky erf.Cs=0.1 erf.moisture_model=Kessler erf.input_sounding_file=${_cbp_dir}/Closure_BoxParity_Anelastic_Kessler/sounding_moist"
        REFERENCE_OPTIONS "${_cbp_ref}"
        SPLIT_OPTIONS "${_cbp_split}"
        FCOMPARE_RTOL "1.0e-9")
    add_test_box_parity(Closure_BoxParity_Anelastic_Deardorff Closure_BoxParity "plt00010"
        COMMON_OPTIONS "${_cbp_ke} erf.anelastic=1 erf.anelastic_type=MidPoint erf.use_fft=true erf.les_type=Deardorff erf.input_sounding_file=${_cbp_dir}/Closure_BoxParity_Anelastic_Deardorff/sounding_dry"
        REFERENCE_OPTIONS "${_cbp_ref}"
        SPLIT_OPTIONS "${_cbp_split}"
        FCOMPARE_RTOL "1.0e-9")
  endif()

  # pblh (2D) and Lturb (3D) are the per-tile PBL height copied out of the
  # scheme; the deck is set up so they differ from column to column.
  add_test_tiling_parity(ABL_MRF_Tiling      ABL_MRF_Tiling "00010" "00010"
      VARYING_3D "Lturb Kmv" VARYING_2D "pblh u_star")
  add_test_tiling_parity(ABL_YSUNew_Tiling   ABL_MRF_Tiling "00010" "00010"
      RUNTIME_OPTIONS "erf.pbl_type=YSUNew erf.most.pblh_calc=YSU"
      VARYING_3D "Lturb Kmv" VARYING_2D "pblh u_star")
  # Legacy YSU aborts in unstable conditions, so cool the surface (a stronger
  # cooling than -0.02 with the 5 m/s wind stops the MOST iteration converging).
  # It covers the full-column assert only: legacy YSU never calls set_pblh, so
  # pblh is left out of the 2D plotfile rather than compared as a constant.
  add_test_tiling_parity(ABL_YSU_Tiling      ABL_MRF_Tiling "00010" "00010"
      RUNTIME_OPTIONS "erf.pbl_type=YSU erf.most.pblh_calc=YSU erf.most.surf_temp_flux=-0.02 'erf.plot2d_vars_1=u_star t_star Olen'"
      VARYING_3D "Lturb Kmv" VARYING_2D "u_star")
  # The PBLH smoothing stencil reads a column its own tile does not own, so it
  # needs its own coverage: with the stencil reading off the end of the array the
  # MRF deck differed by 24.5 m in Lturb (12%) between the tiled and untiled runs.
  # MRF and YSUNew size and fill that halo separately, so both are registered.
  add_test_tiling_parity(ABL_MRF_Tiling_Smooth    ABL_MRF_Tiling "00010" "00010"
      RUNTIME_OPTIONS "erf.enable_pblh_smoothing=true"
      VARYING_3D "Lturb Kmv" VARYING_2D "pblh u_star")
  add_test_tiling_parity(ABL_YSUNew_Tiling_Smooth ABL_MRF_Tiling "00010" "00010"
      RUNTIME_OPTIONS "erf.pbl_type=YSUNew erf.most.pblh_calc=YSU erf.enable_pblh_smoothing=true"
      VARYING_3D "Lturb Kmv" VARYING_2D "pblh u_star")
  # The immersed-boundary-aware MRF and YSUNew (erf.pbl_ib_aware) build their
  # per-column surface and work arrays on the tile work box; a cube by
  # immersed forcing makes them differ from column to column.
  add_test_tiling_parity(PBL_IBAware_MRF_Tiling    PBL_IBAware_Tiling "00010" "00010"
      VARYING_3D "Kmv" VARYING_2D "pblh u_star")
  add_test_tiling_parity(PBL_IBAware_YSUNew_Tiling PBL_IBAware_Tiling "00010" "00010"
      RUNTIME_OPTIONS "erf.pbl_type=YSUNew erf.most.pblh_calc=YSU"
      VARYING_3D "Kmv" VARYING_2D "pblh u_star")
endif()
add_test_r(ABL_InflowFile                    ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(MoistBubble                       ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
# The moist bubble with Kessler rain (the MoistBubble deck with erf.moisture_model=Kessler
# and the rain fields in the plotfile), restarted from step 4 to step 8 on its constant-dz
# mesh: the restart path set the microphysics' minimum dz only on fitted meshes, so the
# sedimentation substep count of the first restarted step came from an uninitialised
# value and the step never finished. The runner is a cmake -P script (MPI, not Windows).
if(ERF_ENABLE_MPI AND NOT WIN32)
add_test_restart_parity(MoistBubble_Kessler_Restart MoistBubble_Kessler_Restart 4 8
    COMMON_OPTIONS "erf.vert_implicit=false"
    FCOMPARE_RTOL "0.0" FCOMPARE_ATOL "0.0")
endif()
add_test_r(SquallLine_2D                     ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(SuperCell_3D                      ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
if(ERF_ENABLE_NETCDF)
  # Distributed terrain ownership and gridded forest interpolation are both
  # exercised by this one-step, two-rank regression.  The NetCDF files are
  # static fixtures so CI does not require ncgen.
  add_test_r(BellForest                       ""  "erf_exec" "plt00001"
      FCOMPARE_RTOL "2.0e-9" FCOMPARE_ATOL "2.0e-9")
endif()
if(ERF_ENABLE_PARTICLES)
  # Production regression: protect the fixed SuperDroplets water-field
  # contract against confusing the constructor sentinel with state width.
  add_test_plotfile_header(Plotfile3D_SuperDropletsSelection "" "erf_exec" "plt00000")
  add_test_r(ParticleAdvect                  ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
  add_test_r(ParticleWoA                     ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
  add_test_r(ParticleAdvect_AMR1_box         ""  "erf_exec" "plt00050" RUNTIME_OPTIONS "erf.vert_implicit=false ")
  add_test_sdm(ParticleAdvect_AMR1_pcount      ""  "erf_exec" "plt00050" 2e-8 3e-9 RUNTIME_OPTIONS "erf.vert_implicit=false ")
  # Skip AMR2_pcount for Debug/RelWithDebInfo builds with AMD GPUs (it freezes!)
  if((CMAKE_BUILD_TYPE STREQUAL "Release") OR (NOT ERF_ENABLE_HIP))
    add_test_sdm(ParticleAdvect_AMR2_pcount    ""  "erf_exec" "plt00050" 1e-7 5e-9 RUNTIME_OPTIONS "erf.vert_implicit=false ")
  endif()
endif( )
# The option name used to be misspelled (ERF_ENABLE_RRGMTP), which kept this
# test unregistered; Tests/test_files/Radiation has never existed, so it is
# registered only once someone adds the inputs.
if(ERF_ENABLE_RRTMGP AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test_files/Radiation")
  add_test_r(Radiation                       ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
endif()

# TwoStream radiation needs no external library and no gold plotfiles: the
# column-physics checker verifies the vertical structure of the heating
# rates, and the header test verifies that qsrc_sw/qsrc_lw are written.
# The column test runs through cmake -P and execute_process, which needs a
# launcher and a resolved executable path; the Windows job builds without
# MPI and resolves test executables through sh -c globs, so it is skipped
# there like the other script-driven tests.
if(ERF_ENABLE_MPI AND NOT WIN32)
  add_test_two_stream_radiation(TwoStream_ColumnHeating "plt00002")
  # Same column over a Witch-of-Agnesi hill on a terrain-fitted mesh: the
  # layer thicknesses come from the nodal heights, every column differs, and
  # the runner's 1-rank vs NRANKS comparison of the diagnostics CSV has a
  # real signal (rank-local means fail it).
  add_test_two_stream_radiation(TwoStream_ColumnHeating_Terrain "plt00002")
endif()
add_test_plotfile_header(Plotfile3D_TwoStreamHeatingSelection "" "erf_exec" "plt00000")

add_test_0(CouetteFlow_x                     "" "erf_exec" "plt00050" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_0(CouetteFlow_y                     "" "erf_exec" "plt00050" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_0(PoiseuilleFlow_x                  "" "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_0(PoiseuilleFlow_y                  "" "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_0(InitSoundingIdeal_stationary      "" "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_0(Deardorff_stationary              "" "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")

if(ERF_ENABLE_PARTICLES)
    # These tests require machine-specific gold files due to platform-dependent initial sampling.
    # Without those gold files they can still be run to completion as smoke tests, which is the
    # only coverage they get on a machine that builds with assertions enabled.
    if(ERF_TEST_ENABLE_EXTRA_SDM_TESTS OR ERF_TEST_SDM_SMOKE_GATED)
        if(NOT ERF_TEST_ENABLE_EXTRA_SDM_TESTS)
            set(ERF_SDM_SMOKE_ONLY TRUE)
        endif()
        # log-normal distribution for radius
        add_test_sdm(SDM_RICO3D_InitSampling         ""  "erf_exec"   "plt00000" 1e-14 2e-13 INPUT_SOUNDING "input_sounding" RUNTIME_OPTIONS "erf.vert_implicit=false ")
        # mass-exponential distribution for mass
        add_test_sdm(SDM_Bubble2D_Adv_InitSampling   ""  "erf_exec"   "plt00000" 1e-14 1e-14 RUNTIME_OPTIONS "erf.vert_implicit=false ")
        # per-box high-multiplicity injection (stochastic cell scatter -> platform-specific gold)
        add_test_sdm(SDM_Bubble2D_PerBoxInjection    ""  "erf_exec"   "plt00050" 5e-12 5e-12 RUNTIME_OPTIONS "erf.vert_implicit=false ")
        # passive advection of particles with injection (takes ~1200s on GitHub Windows CI)
        add_test_sdm(SDM_Bubble2D_Adv_wInjection     "" "erf_exec"  "plt00050" 5e-12 5e-12 RUNTIME_OPTIONS "erf.vert_implicit=false ")
        # INAS sampled initialization for freezing temperature
        add_test_sdm(SDM_Bubble2D_Adv_TfzINAS        ""  "erf_exec"   "plt00000" 1e-14 1e-14 RUNTIME_OPTIONS "erf.vert_implicit=false ")
        # column case to test condensation
        add_test_sdm(SDM_SineMassFlux                "" "erf_exec" "plt00050" 1e-14 1e-14 INPUT_SOUNDING "input_sounding" RUNTIME_OPTIONS "erf.vert_implicit=false ")
        # recycling
        add_test_sdm(SDM_Box3D_Recycling             "" "erf_exec"  "plt00060" 5e-13 1e-14 RUNTIME_OPTIONS "erf.vert_implicit=false ")
        # INAS immersion freezing in a 1D cooling column (Tfz sampling -> platform-specific gold)
        add_test_sdm(SDM_FreezingShaft               "" "erf_exec"  "plt00001" 1e-12 1e-12 INPUT_SOUNDING "input_sounding" RUNTIME_OPTIONS "erf.vert_implicit=false ")
        # Collision processes: stochastic pair sampling and skipped on GPU (RNG/reduction ordering differs).
        if(NOT (ERF_ENABLE_CUDA OR ERF_ENABLE_HIP OR ERF_ENABLE_SYCL))
            # ice-ice aggregation (0D box)
            add_test_sdm(SDM_Box3D_IceAgg            "" "erf_exec"  "plt04500" 1e-12 1e-12 RUNTIME_OPTIONS "erf.vert_implicit=false ")
            # warm-rain coalescence (0D box) -- one test per collection kernel.
            # These dominate the runtime and the kernels they cover already have
            # unit tests, so they are left out of the smoke pass.
            if(NOT ERF_SDM_SMOKE_ONLY)
                add_test_sdm(SDM_Box3D_Coal_Golovin      "" "erf_exec"  "plt04000" 1e-12 1e-12 RUNTIME_OPTIONS "erf.vert_implicit=false ")
                add_test_sdm(SDM_Box3D_Coal_Halls        "" "erf_exec"  "plt04000" 1e-12 1e-12 RUNTIME_OPTIONS "erf.vert_implicit=false ")
                add_test_sdm(SDM_Box3D_Coal_Longs        "" "erf_exec"  "plt04000" 1e-12 1e-12 RUNTIME_OPTIONS "erf.vert_implicit=false ")
                add_test_sdm(SDM_Box3D_Coal_Sedimentation "" "erf_exec" "plt04000" 1e-12 1e-12 RUNTIME_OPTIONS "erf.vert_implicit=false ")
            endif()
            # riming (ice collecting cloud droplets, 1D shaft)
            add_test_sdm(SDM_RimingShaft             "" "erf_exec"  "plt00400" 1e-12 1e-12 INPUT_SOUNDING "input_sounding" RUNTIME_OPTIONS "erf.vert_implicit=false ")
        endif()
        unset(ERF_SDM_SMOKE_ONLY)
    endif()

    # passive advection of particles
    add_test_sdm(SDM_Bubble2D_Adv                "" "erf_exec"  "plt00050" 1e-12 1e-12 RUNTIME_OPTIONS "erf.vert_implicit=false ")
    # super-droplets on a terrain-fitted mesh: covers the pos(2) zeta convention
    add_test_sdm(SDM_Bubble2D_WoA                "" "erf_exec"  "plt00050" 1e-12 1e-12 RUNTIME_OPTIONS "erf.vert_implicit=false ")
    # same case with MFIter tiling forced on: in-place kernels written over
    # grown tiles must still reproduce the untiled answer
    add_test_sdm(SDM_Bubble2D_Tiled              "" "erf_exec"  "plt00050" 1e-12 1e-12 RUNTIME_OPTIONS "erf.vert_implicit=false fabarray.mfiter_tile_size=8 8 8 ")
    add_test_sdm(SDM_Bubble2D_Adv_AMR1           "" "erf_exec"  "plt00050" 1e-12 1e-12 RUNTIME_OPTIONS "erf.vert_implicit=false ")
    add_test_sdm(SDM_Bubble2D_Adv_AMR2           "" "erf_exec"  "plt00025" 1e-12 1e-12 RUNTIME_OPTIONS "erf.vert_implicit=false ")
    add_test_sdm(SDM_Bubble3D_Adv                "" "erf_exec"  "plt00020" 1e-12 1e-12 RUNTIME_OPTIONS "erf.vert_implicit=false ")
    add_test_sdm(SDM_Bubble3D_Adv_AMR1           "" "erf_exec"  "plt00020" 1e-12 1e-12 RUNTIME_OPTIONS "erf.vert_implicit=false ")
    add_test_sdm(SDM_Bubble3D_Adv_AMR2           "" "erf_exec"  "plt00020" 1e-12 1e-12 RUNTIME_OPTIONS "erf.vert_implicit=false ")
    # Gold files are MPI-rank-specific (particle-to-mesh FP ordering).
    if(ERF_ENABLE_MPI)
        add_test_sdm(SDM_MoistBubble2D_AMR1      "" "erf_exec"  "plt00020" 1e-12 1e-12 RUNTIME_OPTIONS "erf.vert_implicit=false ")
        #add_test_sdm(SDM_MoistBubble2D_AMR2      "" "erf_exec" "plt00020" 1e-12 1e-12 RUNTIME_OPTIONS "erf.vert_implicit=false ")
        add_test_sdm(SDM_MoistBubble3D_AMR1      "" "erf_exec"  "plt00020" 1e-12 1e-12 RUNTIME_OPTIONS "erf.vert_implicit=false ")
        add_test_sdm(SDM_MoistBubble3D_AMR2      "" "erf_exec"  "plt00020" 1e-12 1e-12 RUNTIME_OPTIONS "erf.vert_implicit=false ")
    endif()
    # fractional injection (sub-unity per-step multiplicity accumulates to one)
    add_test_sdm(SDM_Bubble2D_FracInjection      "" "erf_exec"  "plt00050" 5e-12 5e-12 RUNTIME_OPTIONS "erf.vert_implicit=false ")
    # condensation/evaporation
    add_test_sdm(SDM_Box3D_Cond                  "" "erf_exec"  "plt00010" 2e-12 3e-13 RUNTIME_OPTIONS "erf.vert_implicit=false ")
    # ice freezing + deposition
    add_test_sdm(SDM_Box3D_IceFrzDep             "" "erf_exec"  "plt00010" 1e-14 1e-12 RUNTIME_OPTIONS "erf.vert_implicit=false ")
    if(NOT (ERF_ENABLE_HIP OR ERF_ENABLE_SYCL))
        # 1D sublimation shaft: monodisperse ice in a subsaturated column (supersedes the 0D box sublimation test)
        add_test_sdm(SDM_SublimationShaft            "" "erf_exec"  "plt00100" 1e-12 1e-12 INPUT_SOUNDING "input_sounding" RUNTIME_OPTIONS "erf.vert_implicit=false ")
    endif()
    # 1D melting layer: melting + mixed-phase fall as ice flakes descend into warmer air (supersedes the 0D box melting test)
    add_test_sdm(SDM_MeltingLayer                "" "erf_exec"  "plt00300" 1e-12 1e-12 INPUT_SOUNDING "input_sounding" RUNTIME_OPTIONS "erf.vert_implicit=false ")
    # terminal velocity
    add_test_sdm(SDM_Box3D_VTerm                 "" "erf_exec"  "plt00001" 5e-13 1e-14 RUNTIME_OPTIONS "erf.vert_implicit=false ")
    # Congestus case
    add_test_sdm(SDM_Congestus3D                 "" "erf_exec"  "plt00020" 5e-13 5e-13 INPUT_SOUNDING "input_sounding" RUNTIME_OPTIONS "erf.vert_implicit=false ")
    # RICO case
    add_test_sdm(SDM_RICO3D                      "" "erf_exec"  "plt00010" 5e-13 5e-13 INPUT_SOUNDING "input_sounding" RUNTIME_OPTIONS "erf.vert_implicit=false ")
    # multispecies setup with dummy water species
    add_test_sdm(SDM_MultiSpecies_Bubble2D       "" "erf_exec"  "plt00001" 5e-12 1e-12 RUNTIME_OPTIONS "erf.vert_implicit=false ")
endif()

#=============================================================================
# Canonical RANS cases (Exec/CanonicalTests/Canonical_RANS)
#
# Each case runs a short smoke deck and then its Python check script, which
# compares planar-averaged numbers against stated targets with tolerances.
# A clean exit alone is never the pass criterion.
#
# The decks run the anelastic projection with the FFT solver (erf.use_fft),
# which no CI configuration builds. The flat decks are therefore run here
# with the MLMG projection (erf.use_fft=false), and the terrain-fitted decks,
# whose general-terrain projection has no non-FFT path, are registered only
# when the build enables FFT (ERF_ENABLE_FFT).
#=============================================================================
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(ERF_RANS_PYTHON "${Python3_EXECUTABLE}")
else()
    set(ERF_RANS_PYTHON "python3")
endif()

function(add_test_rans TEST_NAME CASE_DIR INPUT_FILE NSTEPS CHECK_SCRIPT)
    set(options )
    set(oneValueArgs "RUNTIME_OPTIONS" "NRANKS")
    set(multiValueArgs )
    cmake_parse_arguments(ADD_TEST_RANS "${options}" "${oneValueArgs}"
        "${multiValueArgs}" ${ARGN})

    set(_rans_root ${PROJECT_SOURCE_DIR}/Exec/CanonicalTests/Canonical_RANS)
    set(CURRENT_TEST_SOURCE_DIR ${_rans_root}/${CASE_DIR})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    # shared plotfile reader and check helpers used by every check script
    file(GLOB _rans_py "${_rans_root}/*.py")
    file(COPY ${_rans_py} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")

    if(ERF_ENABLE_MPI)
        if("${ADD_TEST_RANS_NRANKS}" STREQUAL "")
            set(NP ${ERF_TEST_NRANKS})
        else()
            set(NP ${ADD_TEST_RANS_NRANKS})
        endif()
        set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NP} ${MPIEXEC_PREFLAGS}")
    else()
        set(NP 1)
        unset(MPI_COMMANDS)
    endif()

    resolve_test_exe("" "erf_exec" TEST_EXE)

    # plotfile names carry the step number padded to five digits
    set(_step "0000${NSTEPS}")
    string(LENGTH "${_step}" _len)
    math(EXPR _start "${_len} - 5")
    string(SUBSTRING "${_step}" ${_start} 5 _step)
    set(PLTFILE "plt${_step}")

    set(RUNTIME_OPTIONS "max_step=${NSTEPS} erf.plot_int_1=${NSTEPS} erf.check_int=-1 ${ADD_TEST_RANS_RUNTIME_OPTIONS}")
    set(test_log "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log")
    set(check_log "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.check.log")
    # The check script's exit code is the verdict; its table is echoed into
    # the ctest output so a failure shows the measured numbers, and the tail
    # of the run log is echoed when the executable itself exits non-zero.
    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${INPUT_FILE} ${RUNTIME_OPTIONS} > ${test_log} 2>&1 || ( tail -n 60 ${test_log} && false ) && rm -f ${CURRENT_TEST_BINARY_DIR}/CHECK_FAILED && ( ${ERF_RANS_PYTHON} ${CURRENT_TEST_BINARY_DIR}/${CHECK_SCRIPT} --smoke ${CURRENT_TEST_BINARY_DIR}/${PLTFILE} > ${check_log} 2>&1 || touch ${CURRENT_TEST_BINARY_DIR}/CHECK_FAILED ) && cat ${check_log} && test ! -f ${CURRENT_TEST_BINARY_DIR}/CHECK_FAILED")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "rans;regression"
        ATTACHED_FILES_ON_FAIL "${test_log};${check_log}"
    )
endfunction(add_test_rans)

# flat meshes: MLMG projection so the tests run in every build
add_test_rans(RANS_Neutral_ABL_Flat     Neutral_ABL_Flat     inputs_neutral     40  check_neutral.py    RUNTIME_OPTIONS "erf.use_fft=false")
add_test_rans(RANS_Stable_ABL_Flat      Stable_ABL_Flat      inputs_stable      40  check_stable.py     RUNTIME_OPTIONS "erf.use_fft=false")
add_test_rans(RANS_Convective_ABL_Flat  Convective_ABL_Flat  inputs_convective  40  check_convective.py RUNTIME_OPTIONS "erf.use_fft=false")

# Runs a case twice, with OPTIONS_A and with OPTIONS_B on top of RUNTIME_OPTIONS,
# and passes both final plotfiles to the check script.
function(add_test_rans_pair TEST_NAME CASE_DIR INPUT_FILE NSTEPS CHECK_SCRIPT)
    set(options )
    set(oneValueArgs "RUNTIME_OPTIONS" "OPTIONS_A" "OPTIONS_B" "CHECK_OPTIONS" "NRANKS")
    set(multiValueArgs )
    cmake_parse_arguments(ADD_TEST_RANS_PAIR "${options}" "${oneValueArgs}"
        "${multiValueArgs}" ${ARGN})

    set(_rans_root ${PROJECT_SOURCE_DIR}/Exec/CanonicalTests/Canonical_RANS)
    set(CURRENT_TEST_SOURCE_DIR ${_rans_root}/${CASE_DIR})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    file(GLOB _rans_py "${_rans_root}/*.py")
    file(COPY ${_rans_py} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")

    if(ERF_ENABLE_MPI)
        if("${ADD_TEST_RANS_PAIR_NRANKS}" STREQUAL "")
            set(NP ${ERF_TEST_NRANKS})
        else()
            set(NP ${ADD_TEST_RANS_PAIR_NRANKS})
        endif()
        set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NP} ${MPIEXEC_PREFLAGS}")
    else()
        set(NP 1)
        unset(MPI_COMMANDS)
    endif()

    resolve_test_exe("" "erf_exec" TEST_EXE)

    # plotfile names carry the step number padded to five digits
    set(_step "0000${NSTEPS}")
    string(LENGTH "${_step}" _len)
    math(EXPR _start "${_len} - 5")
    string(SUBSTRING "${_step}" ${_start} 5 _step)

    set(_dir ${CURRENT_TEST_BINARY_DIR})
    set(_common "max_step=${NSTEPS} erf.plot_int_1=${NSTEPS} erf.check_int=-1 ${ADD_TEST_RANS_PAIR_RUNTIME_OPTIONS}")
    set(log_a "${_dir}/${TEST_NAME}.a.log")
    set(log_b "${_dir}/${TEST_NAME}.b.log")
    set(check_log "${_dir}/${TEST_NAME}.check.log")
    # Either run's log tail is echoed if it exits non-zero; the check script's
    # exit code is the verdict and its table is echoed into the ctest output.
    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${_dir}/${INPUT_FILE} ${_common} ${ADD_TEST_RANS_PAIR_OPTIONS_A} erf.plot_file_1=${_dir}/a_plt > ${log_a} 2>&1 || ( tail -n 60 ${log_a} && false ) && ${MPI_COMMANDS} ${TEST_EXE} ${_dir}/${INPUT_FILE} ${_common} ${ADD_TEST_RANS_PAIR_OPTIONS_B} erf.plot_file_1=${_dir}/b_plt > ${log_b} 2>&1 || ( tail -n 60 ${log_b} && false ) && rm -f ${_dir}/CHECK_FAILED && ( ${ERF_RANS_PYTHON} ${_dir}/${CHECK_SCRIPT} ${ADD_TEST_RANS_PAIR_CHECK_OPTIONS} ${_dir}/a_plt${_step} ${_dir}/b_plt${_step} > ${check_log} 2>&1 || touch ${_dir}/CHECK_FAILED ) && cat ${check_log} && test ! -f ${_dir}/CHECK_FAILED")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 1800
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${_dir}/"
        LABELS "rans;regression"
        ATTACHED_FILES_ON_FAIL "${log_a};${log_b};${check_log}"
    )
endfunction(add_test_rans_pair)

# The buoyancy production of k must not depend on how theta is diffused: the
# convective case, run compressible with the implicit vertical solve (A) and
# with explicit vertical diffusion (B), must give the same KE to within the
# time-discretisation difference. The command runs through sh, so not on Windows.
if(NOT WIN32)
    add_test_rans_pair(RANS_Convective_ABL_Flat_Buoyancy_kEqn Convective_ABL_Flat inputs_convective 40 check_implicit_explicit_ke.py
        RUNTIME_OPTIONS "erf.anelastic=0 erf.use_fft=false"
        OPTIONS_A "erf.vert_implicit=true" OPTIONS_B "erf.vert_implicit=false"
        CHECK_OPTIONS "--tol 1.0e-4")
    add_test_rans_pair(RANS_Convective_ABL_Flat_Buoyancy_Deardorff Convective_ABL_Flat inputs_convective 40 check_implicit_explicit_ke.py
        RUNTIME_OPTIONS "erf.anelastic=0 erf.use_fft=false erf.rans_type=None erf.les_type=Deardorff erf.plot_vars_1=density theta KE Kmv Khv"
        OPTIONS_A "erf.vert_implicit=true" OPTIONS_B "erf.vert_implicit=false"
        CHECK_OPTIONS "--tol 3.0e-4")
endif()

# The check scripts' own pass/fail logic: kind = "range" accepted half a band
# width outside the band, so every band check was looser than it reads.
# Pure Python, no ERF run.
#
# Registered only when CMake actually found an interpreter: this is the one
# test with the "unit" label that is not a built binary, and "ctest -L unit"
# runs in the gcc, macos, ci and windows workflows. Without the guard,
# ERF_RANS_PYTHON falls back to the bare name "python3" and a configuration
# that has no such executable on PATH fails the whole unit stage on a test
# that exercises no ERF code.
if(Python3_Interpreter_FOUND)
    add_test(RANS_Checks_SelfTest ${ERF_RANS_PYTHON}
        ${PROJECT_SOURCE_DIR}/Exec/CanonicalTests/Canonical_RANS/test_rans_checks.py)
    set_tests_properties(RANS_Checks_SelfTest
        PROPERTIES
        TIMEOUT 60
        PROCESSORS 1
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/Exec/CanonicalTests/Canonical_RANS"
        LABELS "rans;unit")
endif()
if(ERF_ENABLE_FFT)
    # terrain-fitted mesh (FFT-preconditioned projection): wall distance against
    # the exact ridge distance, and the same deck flattened (prob.hmax = 1e-6)
    # against the analytic height, each with the terrain_height and Poisson paths
    add_test_rans(RANS_Neutral_Hill_2D        Neutral_Hill_2D      inputs_hill        40  check_hill.py)
    add_test_rans(RANS_Neutral_Hill_2D_Poisson Neutral_Hill_2D     inputs_hill        40  check_hill.py RUNTIME_OPTIONS "erf.wall_dist_type=poisson")
    add_test_rans(RANS_Flat_Fitted_2D         Neutral_Hill_2D      inputs_hill        40  check_flat_fitted.py RUNTIME_OPTIONS "prob.hmax=1e-6")
    add_test_rans(RANS_Flat_Fitted_2D_Poisson Neutral_Hill_2D      inputs_hill        40  check_flat_fitted.py RUNTIME_OPTIONS "prob.hmax=1e-6 erf.wall_dist_type=poisson")
    add_test_rans(RANS_Neutral_Hill_3D        Neutral_Hill_3D      inputs_hill3d      40  check_hill3d.py)
    add_test_rans(RANS_Neutral_Hill_3D_Poisson Neutral_Hill_3D     inputs_hill3d      40  check_hill3d.py RUNTIME_OPTIONS "erf.wall_dist_type=poisson")
endif()

#=============================================================================
# MOST reference height on flat stretched meshes
#
# run_most_zref.py runs one flat stretched column through the terrain-fitted,
# interpolated and no-terrain MOST lookups plus a uniform 10 m column, and
# checks u* against the log law at the reported reference height and the
# stretched column against the uniform one.  Its exit code is the verdict.
#=============================================================================
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(ERF_MOST_ZREF_PYTHON "${Python3_EXECUTABLE}")
else()
    set(ERF_MOST_ZREF_PYTHON "python3")
endif()

function(add_test_most_zref TEST_NAME)
    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${TEST_NAME})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")

    # 4x4 columns: one rank
    if(ERF_ENABLE_MPI)
        set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}")
    else()
        unset(MPI_COMMANDS)
    endif()

    resolve_test_exe("" "erf_exec" TEST_EXE)

    set(test_log "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log")
    set(test_command sh -c "rm -f ${CURRENT_TEST_BINARY_DIR}/CHECK_FAILED && ( ${ERF_MOST_ZREF_PYTHON} ${CURRENT_TEST_BINARY_DIR}/run_most_zref.py --exe ${TEST_EXE} --mpi-cmd \"${MPI_COMMANDS}\" --workdir ${CURRENT_TEST_BINARY_DIR}/runs > ${test_log} 2>&1 || touch ${CURRENT_TEST_BINARY_DIR}/CHECK_FAILED ) && cat ${test_log} && test ! -f ${CURRENT_TEST_BINARY_DIR}/CHECK_FAILED")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS 1
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression"
        ATTACHED_FILES_ON_FAIL "${test_log}"
    )
endfunction(add_test_most_zref)

add_test_most_zref(MOST_Zref_Stretched)

# Immersed-boundary surface energy balance on the faces of a height-map cube
# (prognostic skin, slab conduction, heat flux into the air), 40 steps.
add_test_r(IBSEB_Cube                        ""  "erf_exec" "plt00040")
add_test_r(PBL_IBAware_MRF_Smoothing         ""  "erf_exec" "plt00010")

#=============================================================================
# Station time-series output
#=============================================================================

# A station series is an interpolation from whichever level and whichever box
# happens to cover the point, so the two things most likely to break it are a
# change of decomposition and a restart.  Both are checked against the run that
# does it in one piece.  Center.dat is the series compared because it is the one
# that varies: the stations in the still air away from the bubble would compare
# a constant against a constant.  The series print ten significant digits, so
# that is what must agree.
add_test_box_parity(StationSampling_BoxParity StationSampling "plt00010"
    COMMON_OPTIONS ""
    REFERENCE_OPTIONS "amr.max_grid_size=1024"
    SPLIT_OPTIONS "amr.max_grid_size_x=32 amr.max_grid_size_y=2 amr.max_grid_size_z=64"
    DATALOG "Output_Stations/Center.dat"
    DATALOG_SIGDIGITS 10)

add_test_restart_parity(StationSampling_Restart StationSampling 4 10
    DATALOG "Output_Stations/Center.dat"
    DATALOG_SIGDIGITS 10)

# The docs say a run with station output turned on gives the same answer as one without,
# and the sampler is built so that it does: it asks BuildPlot3DScratch not to average the
# microphysics state down.  What it still does at every sampled step is fillpatch the state
# on every level up to the highest one a station is on, re-point the qmoist pointers, and
# fill the requested variables over whole levels, so the claim is not free and is tested
# rather than asserted.  Each test runs the same deck with the stations off and on and
# requires the plotfile to be identical bit for bit, not to a tolerance: a diagnostic that
# moves the answer at all is a bug.  The three decks cover the paths that could break it.
#
# AMR, dry: the only one of the three whose station resolves to level 1, so it is the case
# that exercises FillPatchFineLevel.  The deck's own stations are switched off with
# erf.do_station_sampling for the off leg.
add_test_option_parity(StationSampling_AnswerParity StationSampling "plt00010"
    OFF_OPTIONS "erf.do_station_sampling=false"
    ON_OPTIONS  "erf.Center.field=theta magvel vorticity_x vorticity_y vorticity_z pressure"
    REQUIRE_ON_FILE "Output_Stations/Center.dat")

# Surface layer: u_star and t_star are 2D diagnostics of the MOST path, so the on leg reads
# what the surface layer computed as well as the 3D state.
add_test_option_parity(StationSampling_AnswerParity_MOST ABL_MOST "plt00010"
    COMMON_OPTIONS "erf.vert_implicit=false"
    OFF_OPTIONS "erf.do_station_sampling=false"
    ON_OPTIONS  "erf.station_names=T erf.station_sampling_interval=1 erf.T.field=theta magvel vorticity_z pressure u_star t_star erf.T.x=500 erf.T.y=500 erf.T.height_agl=8.0 100.0"
    REQUIRE_ON_FILE "Output_Stations/T.dat")

if(ERF_ENABLE_PARTICLES)
    # Lagrangian microphysics on two levels with TwoWay coupling: the configuration in which
    # BuildPlot3DScratch would average the microphysics state down, and so the one the
    # sync_solution = false argument exists for.  It also has qmoist to re-point on every
    # level.  Unlike the SDM gold-file tests this one compares a run against itself, so it
    # needs neither the machine-specific gold files nor the flags that gate them.
    add_test_option_parity(StationSampling_AnswerParity_SDM SDM_MoistBubble2D_AMR1 "plt00020"
        COMMON_OPTIONS "erf.vert_implicit=false"
        OFF_OPTIONS "erf.do_station_sampling=false"
        ON_OPTIONS  "erf.station_names=T erf.station_sampling_interval=1 erf.T.field=theta magvel vorticity_z qv qc qrain pressure erf.T.x=10000 erf.T.y=200 erf.T.height_agl=500.0 2000.0"
        REQUIRE_ON_FILE "Output_Stations/T.dat")
endif()

#=============================================================================
# Performance tests
#=============================================================================
