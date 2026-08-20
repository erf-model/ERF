# Have CMake discover the number of cores on the node
include(ProcessorCount)
ProcessorCount(PROCESSES)

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
        TIMEOUT 5400
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
        -DMPIEXEC=${MPIEXEC_EXECUTABLE}
        -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
        -DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}
        -DNRANKS=${NP}
        -DTEST_EXE=${TEST_EXE}
        -DINPUT=${test_input}
        -DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}
        -DSIMULATION_LOG=${test_simulation_log}
        -DCHECKER_LOG=${test_checker_log}
        -DCHECKER=${ANELASTIC_WALL_DIFFUSION_CHECKER}
        -DPLOTFILE=${CURRENT_TEST_BINARY_DIR}/plt00002
        -DAXIS=${TEST_AXIS}
        -DTHETA_LO=300.0
        -DTHETA_HI=301.0
        -P ${PROJECT_SOURCE_DIR}/Tests/RunAnelasticWallDiffusion.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 900
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;anelastic;wall-diffusion"
        ATTACHED_FILES_ON_FAIL "${test_simulation_log};${test_checker_log}")
endfunction(add_test_anelastic_wall_diffusion)

# Checker-driven Stage 1 Cloud Chamber tests.  The short run checks the exact
# initial conserved-state correction and a bounded early buoyant response;
# it intentionally avoids a fragile turbulent gold file.
function(add_test_cloud_chamber TEST_NAME MODE)
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)
    set(test_input "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i")
    set(test_simulation_log "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.simulation.log")
    set(test_checker_log "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.checker.log")
    add_test(${TEST_NAME} ${CMAKE_COMMAND}
        -DMPIEXEC=${MPIEXEC_EXECUTABLE}
        -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
        -DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}
        -DNRANKS=${NP}
        -DTEST_EXE=${TEST_EXE}
        -DINPUT=${test_input}
        -DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}
        -DSIMULATION_LOG=${test_simulation_log}
        -DCHECKER_LOG=${test_checker_log}
        -DCHECKER=${CLOUD_CHAMBER_CHECKER}
        -DMODE=${MODE}
        -P ${PROJECT_SOURCE_DIR}/Tests/RunCloudChamber.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 900
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;cloud-chamber"
        ATTACHED_FILES_ON_FAIL "${test_simulation_log};${test_checker_log}")
endfunction(add_test_cloud_chamber)

function(add_test_cloud_chamber_parity TEST_NAME)
    set(TEST_FILES_DIR "CloudChamber_SatAdj")
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)
    add_test(${TEST_NAME} ${CMAKE_COMMAND}
        -DMPIEXEC=${MPIEXEC_EXECUTABLE}
        -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
        -DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}
        -DNRANKS=${NP}
        -DTEST_EXE=${TEST_EXE}
        -DINPUT=${CURRENT_TEST_BINARY_DIR}/CloudChamber_SatAdj.i
        -DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}
        -DCHECKER=${CLOUD_CHAMBER_CHECKER}
        -P ${PROJECT_SOURCE_DIR}/Tests/RunCloudChamberParity.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 1200
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;cloud-chamber"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/budget_off/simulation.log;${CURRENT_TEST_BINARY_DIR}/budget_on/simulation.log;${CURRENT_TEST_BINARY_DIR}/parity.log")
endfunction(add_test_cloud_chamber_parity)

function(add_test_cloud_chamber_budget TEST_NAME MODE SOURCE_NAME)
    set(_cloud_chamber_input_name "${SOURCE_NAME}")
    set(TEST_FILES_DIR "${SOURCE_NAME}")
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)
    add_test(${TEST_NAME} ${CMAKE_COMMAND}
        -DMPIEXEC=${MPIEXEC_EXECUTABLE}
        -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
        -DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}
        -DNRANKS=${NP}
        -DTEST_EXE=${TEST_EXE}
        -DINPUT=${CURRENT_TEST_BINARY_DIR}/${_cloud_chamber_input_name}.i
        -DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}
        -DCHECKER=${CLOUD_CHAMBER_CHECKER}
        -DMODE=${MODE}
        -P ${PROJECT_SOURCE_DIR}/Tests/RunCloudChamberBudget.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 900
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;cloud-chamber"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/simulation.log;${CURRENT_TEST_BINARY_DIR}/checker.log;${CURRENT_TEST_BINARY_DIR}/cloud_chamber_budget.dat")
endfunction(add_test_cloud_chamber_budget)

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
        -DMPIEXEC=${MPIEXEC_EXECUTABLE}
        -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
        -DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}
        -DTEST_EXE=${TEST_EXE}
        -DINPUT=${test_input}
        -DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}
        -DLOG=${test_log}
        -DOUTPUT_DIRECTORY=${output_directory}
        -DOUTPUT_ARTIFACT=${output_artifact}
        -P ${PROJECT_SOURCE_DIR}/Tests/RunCloudChamberConfigSuccess.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 180
        PROCESSORS 1
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;cloud-chamber;configuration"
        ATTACHED_FILES_ON_FAIL "${test_log};${output_artifact}")
endfunction(add_test_cloud_chamber_legacy_config)

function(add_test_cloud_chamber_openmp TEST_NAME)
    set(TEST_FILES_DIR "CloudChamber_SatAdj")
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)
    add_test(${TEST_NAME} ${CMAKE_COMMAND}
        -DMPIEXEC=${MPIEXEC_EXECUTABLE}
        -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
        -DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}
        -DNRANKS=${NP}
        -DTEST_EXE=${TEST_EXE}
        -DINPUT=${CURRENT_TEST_BINARY_DIR}/CloudChamber_SatAdj.i
        -DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}
        -DCHECKER=${CLOUD_CHAMBER_CHECKER}
        -DCMAKE_COMMAND=${CMAKE_COMMAND}
        -P ${PROJECT_SOURCE_DIR}/Tests/RunCloudChamberOpenMP.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 1200
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
        -DMPIEXEC=${MPIEXEC_EXECUTABLE}
        -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
        -DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}
        -DNRANKS=${NP}
        -DTEST_EXE=${TEST_EXE}
        -DINPUT=${test_input}
        -DRUNTIME_OPTIONS=${RUNTIME_OPTIONS}
        -DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}
        -DLOG=${test_log}
        -DCHECKER=${SHOC_PLOTFILE_CHECKER}
        -DCHECK_MODE=${ADD_TEST_SHOC_R_CHECK_MODE}
        -DINITIAL=${CURRENT_TEST_BINARY_DIR}/plt00000
        -DMIDPOINT=${CURRENT_TEST_BINARY_DIR}/plt00010
        -DFINAL=${CURRENT_TEST_BINARY_DIR}/${PLTFILE}
        -DFCOMPARE=${FCOMPARE_EXE}
        -DGOLD_DIFFERENTIAL=${SHOC_GOLD_DIFFERENTIAL}
        -DGOLD_COMPARISON=${_shoc_gold_comparison}
        -DGOLD_MODE=${ADD_TEST_SHOC_R_GOLD_MODE}
        -DRTOL=${ERF_TEST_FCOMPARE_RTOL}
        -DATOL=${ERF_TEST_FCOMPARE_ATOL}
        -DGOLD=${PLOT_GOLD}
        -DSKIP_GOLD=${ADD_TEST_SHOC_R_SKIP_GOLD}
        -P ${PROJECT_SOURCE_DIR}/Tests/RunShocRegression.cmake)
    if(ADD_TEST_SHOC_R_TIMEOUT)
        set(_shoc_timeout "${ADD_TEST_SHOC_R_TIMEOUT}")
    else()
        set(_shoc_timeout 900)
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
        -DMPIEXEC=${MPIEXEC_EXECUTABLE}
        -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
        -DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}
        -DNRANKS=${NP}
        -DTEST_EXE=${TEST_EXE}
        -DBASELINE_INPUT=${_baseline_dir}/SHOC_Stable_Clear.i
        -DMUTANT_INPUT=${_mutant_dir}/SHOC_Stable_Clear.i
        -DBASELINE_OPTIONS=
        -DMUTANT_OPTIONS=${MUTATION_OPTION}
        -DBASELINE_WORKING_DIRECTORY=${_baseline_dir}
        -DMUTANT_WORKING_DIRECTORY=${_mutant_dir}
        -DBASELINE_LOG=${_baseline_log}
        -DMUTANT_LOG=${_mutant_log}
        -DCHECKER=${SHOC_MUTATION_DIFFERENTIAL}
        -DTARGET_FIELD=${TARGET_FIELD}
        -DMIN_FINAL_DIFFERENCE=${MIN_FINAL_DIFFERENCE}
        -DMIN_BASELINE_EVOLUTION=${MIN_BASELINE_EVOLUTION}
        -DMAX_MUTANT_TO_BASELINE_EVOLUTION_RATIO=${MAX_MUTANT_TO_BASELINE_EVOLUTION_RATIO}
        -P ${PROJECT_SOURCE_DIR}/Tests/RunShocMutationRegression.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 900
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;shoc;mutation"
        ATTACHED_FILES_ON_FAIL "${_baseline_log};${_mutant_log}")
endfunction(add_test_shoc_mutation)

if(ERF_ENABLE_MPI)
add_test_anelastic_wall_diffusion(AnelasticWallDiffusion_X 0)
add_test_anelastic_wall_diffusion(AnelasticWallDiffusion_Y 1)
add_test_anelastic_wall_diffusion(AnelasticWallDiffusion_Z 2)
add_test_cloud_chamber(CloudChamber_Dry dry)
add_test_cloud_chamber_legacy_config(CloudChamber_Legacy_Config)
add_test_cloud_chamber_budget(CloudChamber_Dry_ThermalBudget thermal_budget CloudChamber_Dry)
add_test_cloud_chamber(CloudChamber_SatAdj cloudy)
add_test_cloud_chamber_parity(CloudChamber_SatAdj_Parity)
add_test_cloud_chamber_budget(CloudChamber_SatAdj_AllDry all_dry CloudChamber_SatAdj_AllDry)
add_test_cloud_chamber_budget(CloudChamber_SatAdj_WetBudget wet_budget CloudChamber_SatAdj_WetBudget)
if(ERF_ENABLE_OPENMP)
add_test_cloud_chamber_openmp(CloudChamber_SatAdj_OpenMP)
endif()
add_test(SHOC_Unstable_Cloud_SatAdj_vs_NoCond
    ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
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
        TIMEOUT 5400
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
        TIMEOUT 5400
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
        TIMEOUT 5400
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
    add_test_shoc_r(SHOC_Stable_Clear "" "erf_exec" "plt00020"
        TEST_FILES_DIR "SHOC_Stable_Clear"
        CHECK_MODE "stable_clear"
        GOLD_COMPARISON "fcompare"
        LABELS regression shoc
        TIMEOUT 900)
    add_test_shoc_r(SHOC_Stable_Cloud "" "erf_exec" "plt00020"
        TEST_FILES_DIR "SHOC_Stable_Cloud"
        CHECK_MODE "stable_cloud"
        GOLD_COMPARISON "field_aware"
        GOLD_MODE "stable_cloud"
        LABELS regression shoc
        TIMEOUT 900)
    add_test_shoc_r(SHOC_Unstable_Clear_BOMEX "" "erf_exec" "plt00020"
        TEST_FILES_DIR "SHOC_Unstable_Clear_BOMEX"
        CHECK_MODE "unstable_clear"
        GOLD_COMPARISON "fcompare"
        LABELS regression shoc
        TIMEOUT 900)
    add_test_shoc_r(SHOC_Unstable_Cloud_SatAdj "" "erf_exec" "plt00020"
        TEST_FILES_DIR "SHOC_Unstable_Cloud"
        INPUT_FILE "SHOC_Unstable_Cloud.i"
        CHECK_MODE "unstable_cloud"
        GOLD_COMPARISON "field_aware"
        GOLD_MODE "unstable_cloud"
        RUNTIME_OPTIONS "erf.moisture_model=SatAdj erf.buoyancy_type=1 "
        LABELS regression shoc microphysics
        TIMEOUT 900)
    add_test_shoc_r(SHOC_Unstable_Cloud_NoCond "" "erf_exec" "plt00020"
        TEST_FILES_DIR "SHOC_Unstable_Cloud"
        INPUT_FILE "SHOC_Unstable_Cloud.i"
        CHECK_MODE "unstable_cloud_nocond"
        GOLD_COMPARISON "field_aware"
        GOLD_MODE "unstable_cloud_nocond"
        RUNTIME_OPTIONS "erf.moisture_model=MoistNoCondensation erf.buoyancy_type=1 "
        LABELS regression shoc microphysics
        TIMEOUT 900)
    add_test_shoc_r(SHOC_Unstable_Cloud_SatAdj_Property "" "erf_exec" "plt00020"
        TEST_FILES_DIR "SHOC_Unstable_Cloud"
        INPUT_FILE "SHOC_Unstable_Cloud.i"
        CHECK_MODE "unstable_cloud"
        GOLD_COMPARISON "field_aware"
        GOLD_MODE "unstable_cloud"
        RUNTIME_OPTIONS "erf.moisture_model=SatAdj erf.buoyancy_type=1 "
        SKIP_GOLD
        LABELS regression shoc microphysics property
        TIMEOUT 900)
    add_test_shoc_r(SHOC_Unstable_Cloud_NoCond_Property "" "erf_exec" "plt00020"
        TEST_FILES_DIR "SHOC_Unstable_Cloud"
        INPUT_FILE "SHOC_Unstable_Cloud.i"
        CHECK_MODE "unstable_cloud_nocond"
        GOLD_COMPARISON "field_aware"
        GOLD_MODE "unstable_cloud_nocond"
        RUNTIME_OPTIONS "erf.moisture_model=MoistNoCondensation erf.buoyancy_type=1 "
        SKIP_GOLD
        LABELS regression shoc microphysics property
        TIMEOUT 900)
    add_test_shoc_r(SHOC_Unstable_Cloud_Kessler "" "erf_exec" "plt00020"
        TEST_FILES_DIR "SHOC_Unstable_Cloud"
        INPUT_FILE "SHOC_Unstable_Cloud_Kessler.i"
        CHECK_MODE "unstable_cloud_kessler"
        GOLD_COMPARISON "field_aware"
        GOLD_MODE "unstable_cloud_kessler"
        LABELS regression shoc microphysics
        TIMEOUT 900)
    add_test_shoc_r(SHOC_Unstable_Cloud_WSM6 "" "erf_exec" "plt00020"
        TEST_FILES_DIR "SHOC_Unstable_Cloud"
        INPUT_FILE "SHOC_Unstable_Cloud_WSM6.i"
        CHECK_MODE "unstable_cloud_wsm6"
        GOLD_COMPARISON "field_aware"
        GOLD_MODE "unstable_cloud_wsm6"
        LABELS regression shoc microphysics
        TIMEOUT 900)
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
add_test_r(ABL_MYNN_PBL                      ""  "erf_exec" "plt00100" INPUT_SOUNDING "input_sounding_GABLS1" RUNTIME_OPTIONS "erf.vert_implicit=false " )
add_test_r(ABL_InflowFile                    ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(MoistBubble                       ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(SquallLine_2D                     ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_r(SuperCell_3D                      ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
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
if(ERF_ENABLE_RRGMTP)
  add_test_r(Radiation                       ""  "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
endif()

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
    # passive advection of particles with injection
    add_test_sdm(SDM_Bubble2D_Adv_wInjection     "" "erf_exec"  "plt00050" 5e-12 5e-12 RUNTIME_OPTIONS "erf.vert_implicit=false ")
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
# Performance tests
#=============================================================================
