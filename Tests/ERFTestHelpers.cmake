# General CTest helpers: resolve_test_exe, setup_test, add_test_plotfile_header, add_test_r,
# add_test_restart_abort and add_test_abort. Tests/CTestList.cmake includes this file, and
# a project that adds ERF with add_subdirectory can include it too and register its own
# tests the way ERF does:
#
#   add_subdirectory(Submodules/ERF ERF)
#   set(ERF_TEST_GOLD_FILES_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/GoldFiles")
#   include(Submodules/ERF/Tests/ERFTestHelpers.cmake)
#   add_test_r(MyCase "" erf_exec plt00010)
#
# A parent that builds its own executable against ERF (erf_setup_executable) names that
# target instead of erf_exec: add_test_r(MyCase "" my_exec plt00010). add_test_abort and
# add_test_restart_abort take no executable argument; they run ${ERF_TEST_EXE} when the
# caller sets that normal variable, erf_exec otherwise.
#
# Each test reads its files from ${CMAKE_CURRENT_SOURCE_DIR}/test_files/<test name> of the
# directory that calls the helper (or test_files/${TEST_FILES_DIR} when the caller sets
# TEST_FILES_DIR), and add_test_r compares against ${ERF_TEST_GOLD_FILES_DIRECTORY}/<test
# name>. A normal variable of that name takes the place of ERF's cache entry in the calling
# directory only; set it after add_subdirectory, or ERF's own tests also use it. The ERF
# executables and checker scripts are found through ERF_BINARY_DIR and ERF_SOURCE_DIR,
# which project(ERF) sets as cache entries, so they name ERF's folders from any directory
# (a parent project that is itself named ERF would shadow them). Include the file in each
# directory that calls the helpers: ERF_SH_MPI_LAUNCHER is a normal variable of the
# including directory. add_test_r needs FCOMPARE_EXE (ERF_ENABLE_FCOMPARE=ON or
# ERF_ENABLE_TESTS=ON), ERF_TEST_GOLD_FILES_DIRECTORY, ERF_TEST_FCOMPARE_RTOL and _ATOL
# (cache entries with ERF_ENABLE_TESTS=ON, or normal variables of the caller) and, with
# ERF_ENABLE_MPI, ERF_TEST_NRANKS; it stops the configure when one of them is empty.

if(NOT DEFINED ERF_BINARY_DIR OR NOT DEFINED ERF_SOURCE_DIR)
    message(FATAL_ERROR "ERFTestHelpers.cmake: include it after project(ERF) or "
                        "add_subdirectory(<ERF source tree>)")
endif()

# Shared handling of multi-word MPI launchers (e.g., "flux run").
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")

# The launcher for the test commands handed to `sh -c`. Each launcher token is quoted,
# so a path with spaces ("C:/Program Files/Microsoft MPI/Bin/mpiexec.exe") stays one
# word and a multi-word launcher ("flux run") stays several.
set(ERF_SH_MPI_LAUNCHER "")
if(ERF_ENABLE_MPI)
    erf_mpi_launcher_tokens(_mpi_launcher_tokens "${MPIEXEC_EXECUTABLE}")
    foreach(_token IN LISTS _mpi_launcher_tokens)
        string(APPEND ERF_SH_MPI_LAUNCHER "\"${_token}\" ")
    endforeach()
endif()

#=============================================================================
# Functions for adding tests
#=============================================================================
function(resolve_test_exe TEST_DIR TEST_EXE OUT_VAR)
    # An executable target a parent project defined (erf_setup_executable) is built in the
    # parent's binary folders, outside ERF's: name it by its target file. Every caller hands
    # the result to add_test, which evaluates the generator expression.
    if(TARGET "${TEST_EXE}")
        get_target_property(_exe_binary_dir "${TEST_EXE}" BINARY_DIR)
        file(RELATIVE_PATH _exe_rel "${ERF_BINARY_DIR}" "${_exe_binary_dir}")
        if(_exe_rel MATCHES "^\\.\\./" OR IS_ABSOLUTE "${_exe_rel}")
            set(${OUT_VAR} "$<TARGET_FILE:${TEST_EXE}>" PARENT_SCOPE)
            return()
        endif()
    endif()
    get_property(_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
    if(WIN32 AND _multi_config)
        # Multi-config generators place binaries in a config subdir, and which config is
        # built is not known here.  The tests that launch through `sh -c` let the shell
        # expand the wildcard; the ones driven through `cmake -P` call execute_process,
        # which never invokes a shell, and expand it with erf_resolve_executable
        # (Tests/ResolveExecutable.cmake) out of the CONFIG those tests are given.
        set(${OUT_VAR} "${ERF_BINARY_DIR}/Exec/${TEST_DIR}/*/${TEST_EXE}.exe" PARENT_SCOPE)
    else()
        set(_exe_in_subdir "${ERF_BINARY_DIR}/Exec/${TEST_DIR}/${TEST_EXE}${CMAKE_EXECUTABLE_SUFFIX}")
        set(_exe_in_root  "${ERF_BINARY_DIR}/Exec/${TEST_EXE}${CMAKE_EXECUTABLE_SUFFIX}")
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
        set(MPI_COMMANDS "${ERF_SH_MPI_LAUNCHER}${MPIEXEC_NUMPROC_FLAG} ${NP} ${MPIEXEC_PREFLAGS}")
        set(MPI_FCOMP_COMMANDS "${ERF_SH_MPI_LAUNCHER}${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}")
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
    set(header_checker "${ERF_SOURCE_DIR}/Tests/CheckPlotfileHeader.cmake")
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

    add_test(NAME ${TEST_NAME} COMMAND ${test_command})
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

    # Inputs a parent project may not have set; ERF's own configure always sets them.
    set(_required FCOMPARE_EXE ERF_TEST_GOLD_FILES_DIRECTORY)
    if("${ADD_TEST_R_FCOMPARE_RTOL}" STREQUAL "")
        list(APPEND _required ERF_TEST_FCOMPARE_RTOL)
    endif()
    if("${ADD_TEST_R_FCOMPARE_ATOL}" STREQUAL "")
        list(APPEND _required ERF_TEST_FCOMPARE_ATOL)
    endif()
    if(ERF_ENABLE_MPI)
        list(APPEND _required ERF_TEST_NRANKS)
    endif()
    foreach(_var IN LISTS _required)
        if("${${_var}}" STREQUAL "")
            message(FATAL_ERROR "add_test_r(${TEST_NAME}): ${_var} is empty; see the header of "
                                "Tests/ERFTestHelpers.cmake")
        endif()
    endforeach()

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

    add_test(NAME ${TEST_NAME} COMMAND ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log"
    )
endfunction(add_test_r)

# Restart abort: run a deck to a checkpoint, restart from it, and require the restart to stop
# with EXPECTED_MESSAGE. For guards that only a restart can reach, which add_test_abort cannot
# test because it runs a single leg and no checkpoint exists yet.  CHK_NRANKS and
# RESTART_NRANKS are separate so a test can reproduce a restart on more ranks than the
# checkpoint was written with, which is its own code path.
function(add_test_restart_abort TEST_NAME TEST_FILES_DIR STEP_CHK EXPECTED_MESSAGE)
    set(oneValueArgs "COMMON_OPTIONS" "RESTART_OPTIONS" "CHK_NRANKS" "RESTART_NRANKS" "RUN_TIMEOUT")
    cmake_parse_arguments(ADD_TEST_RA "" "${oneValueArgs}" "" ${ARGN})
    setup_test()
    if(DEFINED ERF_TEST_EXE AND NOT "${ERF_TEST_EXE}" STREQUAL "")
        resolve_test_exe("" "${ERF_TEST_EXE}" TEST_EXE)
    else()
        resolve_test_exe("" "erf_exec" TEST_EXE)
    endif()

    set(_chk_nranks "${NP}")
    set(_restart_nranks "${NP}")
    if(NOT "${ADD_TEST_RA_CHK_NRANKS}" STREQUAL "")
        set(_chk_nranks "${ADD_TEST_RA_CHK_NRANKS}")
    endif()
    if(NOT "${ADD_TEST_RA_RESTART_NRANKS}" STREQUAL "")
        set(_restart_nranks "${ADD_TEST_RA_RESTART_NRANKS}")
    endif()
    set(_run_timeout 600)
    set(_ctest_timeout 600)
    if(DEFINED ADD_TEST_RA_RUN_TIMEOUT)
        set(_run_timeout "${ADD_TEST_RA_RUN_TIMEOUT}")
        math(EXPR _ctest_timeout "2 * ${_run_timeout} + 600")
    endif()
    # The watchdog has to cover both legs, and the processor reservation the wider leg.
    set(_procs "${_chk_nranks}")
    if(_restart_nranks GREATER _procs)
        set(_procs "${_restart_nranks}")
    endif()

    add_test(NAME ${TEST_NAME} COMMAND ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DCHK_NRANKS=${_chk_nranks}"
        "-DRESTART_NRANKS=${_restart_nranks}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DCONFIG=$<CONFIG>"
        "-DINPUT=${CURRENT_TEST_BINARY_DIR}/${TEST_FILES_DIR}.i"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DSTEP_CHK=${STEP_CHK}"
        "-DEXPECTED_MESSAGE=${EXPECTED_MESSAGE}"
        "-DRUN_TIMEOUT=${_run_timeout}"
        "-DCOMMON_OPTIONS=${ADD_TEST_RA_COMMON_OPTIONS}"
        "-DRESTART_OPTIONS=${ADD_TEST_RA_RESTART_OPTIONS}"
        -P ${ERF_SOURCE_DIR}/Tests/RunRestartAbort.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT ${_ctest_timeout}
        PROCESSORS ${_procs}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;restart-parity"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/run/checkpoint.log;${CURRENT_TEST_BINARY_DIR}/run/restart.log")
endfunction(add_test_restart_abort)

# Start-up check: copy SOURCE_DIR, run INPUT_FILE on one rank with RUNTIME_OPTIONS
# that break a start-up requirement, and pass when the run stops with
# EXPECTED_MESSAGE in its output. The run is meant to abort, so its exit status is
# dropped by the pipe into tee (a ';' here would split the CMake command list).
#
# amrex.call_addr2line = 0 because the abort is the point of the test: AMReX's
# SIGABRT handler runs addr2line once per stack frame, and on a build with
# debug info that is about 55 s per test -- far longer than the run itself
# (measured: the abort message and "See Backtrace.0 file for details" 56 s
# apart in CI). Backtrace.0 is still written, with raw addresses, for a test
# that stops somewhere unexpected; `addr2line -Cpfie <exe> <address>` resolves
# them by hand.
function(add_test_abort TEST_NAME SOURCE_DIR INPUT_FILE EXPECTED_MESSAGE RUNTIME_OPTIONS)
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    file(GLOB TEST_FILES "${SOURCE_DIR}/*")
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")

    if(ERF_ENABLE_MPI)
        set(MPI_COMMANDS "${ERF_SH_MPI_LAUNCHER}${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}")
    else()
        unset(MPI_COMMANDS)
    endif()

    if(DEFINED ERF_TEST_EXE AND NOT "${ERF_TEST_EXE}" STREQUAL "")
        resolve_test_exe("" "${ERF_TEST_EXE}" TEST_EXE)
    else()
        resolve_test_exe("" "erf_exec" TEST_EXE)
    endif()

    set(test_log "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log")
    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${INPUT_FILE} max_step=1 erf.plot_int_1=-1 erf.plot_int_2=-1 erf.check_int=-1 amrex.call_addr2line=0 ${RUNTIME_OPTIONS} 2>&1 | tee ${test_log}")

    add_test(NAME ${TEST_NAME} COMMAND ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS 1
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression"
        PASS_REGULAR_EXPRESSION "${EXPECTED_MESSAGE}"
        ATTACHED_FILES_ON_FAIL "${test_log}"
    )
endfunction(add_test_abort)
