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
        # Multi-config generators place binaries in a config subdir, and which config is
        # built is not known here.  The tests that launch through `sh -c` let the shell
        # expand the wildcard; the ones driven through `cmake -P` call execute_process,
        # which never invokes a shell, and expand it with erf_resolve_executable
        # (Tests/ResolveExecutable.cmake) out of the CONFIG those tests are given.
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
    set(oneValueArgs "RUNTIME_OPTIONS" "SEB_PARITY_PLOTFILE" "SEB_PARITY_TOL" "SEB_EVOLVED_FROM")
    # CHECK_LEVELS is multi-value: as a one-value arg CMake's list semantics
    # split "0;1" into two arguments and only the first was ever seen, so the
    # fine level went unchecked and the test passed vacuously.
    set(multiValueArgs "CHECK_LEVELS" "DIAG_LEVELS" "SEB_CREATED_FROM" "SEB_REGRIDDED_FROM")
    cmake_parse_arguments(ADD_TEST_TSR "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    # Join with a comma, not a semicolon: a semicolon inside a -D argument is
    # split again when the COMMAND is built. The runner splits on the comma.
    set(tsr_check_levels "0")
    if(DEFINED ADD_TEST_TSR_CHECK_LEVELS)
        string(JOIN "," tsr_check_levels ${ADD_TEST_TSR_CHECK_LEVELS})
    endif()
    # Levels that must appear in the diagnostics CSV. Usually the same list, but a nested
    # patch is checked in the plotfile while writing no CSV row of its own: it never sweeps,
    # so it has no flux diagnostics to report.
    #
    # DEFINED, not truthiness: if(<var>) treats the string "0" as false, so a list of just
    # level 0 would silently fall back to the default.
    set(tsr_diag_levels "${tsr_check_levels}")
    if(DEFINED ADD_TEST_TSR_DIAG_LEVELS)
        string(JOIN "," tsr_diag_levels ${ADD_TEST_TSR_DIAG_LEVELS})
    endif()
    setup_test()
    # Optional second check: the prognostic surface state a refined run keeps on every
    # level must satisfy the relation average_down establishes -- each coarse cell holding
    # the mean of the fine cells above it. Only the decks that switch the prognostic SEB
    # on ask for this.
    set(tsr_seb_parity "")
    if(DEFINED ADD_TEST_TSR_SEB_PARITY_PLOTFILE AND NOT "${ADD_TEST_TSR_SEB_PARITY_PLOTFILE}" STREQUAL "")
        set(tsr_seb_parity "${CURRENT_TEST_BINARY_DIR}/${ADD_TEST_TSR_SEB_PARITY_PLOTFILE}")
    endif()
    set(tsr_seb_tol "1.0e-8")
    if(DEFINED ADD_TEST_TSR_SEB_PARITY_TOL AND NOT "${ADD_TEST_TSR_SEB_PARITY_TOL}" STREQUAL "")
        set(tsr_seb_tol "${ADD_TEST_TSR_SEB_PARITY_TOL}")
    endif()
    # A shallow nest never sweeps, so there is no fine solution to compare against; what
    # must hold is that the coarse level kept evolving underneath it.
    set(tsr_seb_evolved "")
    if(DEFINED ADD_TEST_TSR_SEB_EVOLVED_FROM AND NOT "${ADD_TEST_TSR_SEB_EVOLVED_FROM}" STREQUAL "")
        set(tsr_seb_evolved "${ADD_TEST_TSR_SEB_EVOLVED_FROM}")
    endif()
    # A level created mid-run: the three 2D plotfiles (initial, and the last two before
    # the level existed) the checker extrapolates the parent's surface from. Joined with a
    # comma for the same reason as CHECK_LEVELS.
    set(tsr_seb_created "")
    if(DEFINED ADD_TEST_TSR_SEB_CREATED_FROM)
        set(tsr_seb_created_paths "")
        foreach(created_plt ${ADD_TEST_TSR_SEB_CREATED_FROM})
            list(APPEND tsr_seb_created_paths "${CURRENT_TEST_BINARY_DIR}/${created_plt}")
        endforeach()
        string(JOIN "," tsr_seb_created ${tsr_seb_created_paths})
    endif()
    # A regrid that moves an existing fine level: the same three plotfiles, with the fine
    # level present in the last two on the grids the regrid replaces.
    set(tsr_seb_regridded "")
    if(DEFINED ADD_TEST_TSR_SEB_REGRIDDED_FROM)
        set(tsr_seb_regridded_paths "")
        foreach(regridded_plt ${ADD_TEST_TSR_SEB_REGRIDDED_FROM})
            list(APPEND tsr_seb_regridded_paths "${CURRENT_TEST_BINARY_DIR}/${regridded_plt}")
        endforeach()
        string(JOIN "," tsr_seb_regridded ${tsr_seb_regridded_paths})
    endif()
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
        "-DCHECK_LEVELS=${tsr_check_levels}"
        "-DDIAG_LEVELS=${tsr_diag_levels}"
        "-DSEB_PARITY_PLOTFILE=${tsr_seb_parity}"
        "-DSEB_PARITY_TOL=${tsr_seb_tol}"
        "-DSEB_EVOLVED_FROM=${tsr_seb_evolved}"
        "-DSEB_CREATED_FROM=${tsr_seb_created}"
        "-DSEB_REGRIDDED_FROM=${tsr_seb_regridded}"
        "-DSEB_PARITY_CHECKER=${TWO_STREAM_SEB_PARITY_CHECKER}"
        "-DFEXTRACT=${FEXTRACT_EXE}"
        "-DPYTHON_EXE=${ERF_TEST_PYTHON}"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunTwoStreamRadiation.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;radiation"
        ATTACHED_FILES_ON_FAIL "${test_simulation_log};${test_checker_log}")
endfunction(add_test_two_stream_radiation)

# Run a two-stream deck with the surface energy balance's H and LE taken from the
# surface layer and from the scalar defaults, and check the balance removed the surface
# layer's fluxes from the ground (Tests/check_two_stream_seb_flux_source.py).
function(add_test_two_stream_seb_flux_source TEST_NAME)
    set(oneValueArgs "CHECKER_OPTIONS" "DT")
    cmake_parse_arguments(ADD_TEST_SEBFS "" "${oneValueArgs}" "" ${ARGN})
    set(_sebfs_dt "1.0")
    if(DEFINED ADD_TEST_SEBFS_DT AND NOT "${ADD_TEST_SEBFS_DT}" STREQUAL "")
        set(_sebfs_dt "${ADD_TEST_SEBFS_DT}")
    endif()
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)
    add_test(NAME ${TEST_NAME} COMMAND ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DCONFIG=$<CONFIG>"
        "-DINPUT=${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DFEXTRACT=${FEXTRACT_EXE}"
        "-DPYTHON_EXE=${ERF_TEST_PYTHON}"
        "-DCHECKER=${TWO_STREAM_SEB_FLUX_SOURCE_CHECKER}"
        "-DSTEPS=10"
        "-DDT=${_sebfs_dt}"
        "-DHEAT_CAPACITY=2.0e4"
        "-DCHECKER_OPTIONS=${ADD_TEST_SEBFS_CHECKER_OPTIONS}"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunTwoStreamSEBFluxSource.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 1200
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;radiation"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/surface_layer/simulation.log;${CURRENT_TEST_BINARY_DIR}/defaults/simulation.log;${CURRENT_TEST_BINARY_DIR}/two_way/simulation.log;${CURRENT_TEST_BINARY_DIR}/checker.log")
endfunction(add_test_two_stream_seb_flux_source)

# Two-stream radiation feeding Noah-MP on two levels: a fine level that runs Noah-MP on its
# own nested land file and sweeps its own columns, the same land under a nested patch that
# takes its radiation from level 0, a fine level without a land file, and a regrid that must
# stop (Tests/RunTwoStreamNoahMPLevels.cmake, Tests/check_two_stream_noahmp_levels.py).
function(add_test_two_stream_noahmp_levels TEST_NAME)
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)
    add_test(NAME ${TEST_NAME} COMMAND ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DCONFIG=$<CONFIG>"
        "-DINPUT=${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DLAND_DIR=${PROJECT_SOURCE_DIR}/Exec/RegTests/NoahMP_Ideal"
        "-DFEXTRACT=${FEXTRACT_EXE}"
        "-DPYTHON_EXE=${ERF_TEST_PYTHON}"
        "-DCHECKER=${TWO_STREAM_NOAHMP_LEVELS_CHECKER}"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunTwoStreamNoahMPLevels.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 1200
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;radiation;noahmp"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/own/simulation.log;${CURRENT_TEST_BINARY_DIR}/own/checker.log;${CURRENT_TEST_BINARY_DIR}/nested/simulation.log;${CURRENT_TEST_BINARY_DIR}/nested/checker.log;${CURRENT_TEST_BINARY_DIR}/interp/simulation.log;${CURRENT_TEST_BINARY_DIR}/interp/checker.log;${CURRENT_TEST_BINARY_DIR}/regrid/simulation.log")
endfunction(add_test_two_stream_noahmp_levels)

# The two-stream balance's skin and soil moisture both coupled to the surface layer, on the
# deck of TwoStream_SEBSurfaceLayerFluxes in eight legs: the linear soil-water factor at the
# wilting point, at field capacity and in between; that soil as a Noah-MP soil type, bare
# (the soil resistance and pore humidity), with grassland at fraction 0 (which must match
# bare) and bare at the wilting point (which must not evaporate); under
# grassland (canopy and soil resistances); and a one-step leg without erf.most.z0, whose
# job_info must record the tables' roughness
# (Tests/RunTwoStreamSEBMoisture.cmake, Tests/check_two_stream_seb_moisture.py).
function(add_test_two_stream_seb_moisture TEST_NAME)
    set(TEST_FILES_DIR "TwoStream_SEBSurfaceLayerFluxes")
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)
    add_test(NAME ${TEST_NAME} COMMAND ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DCONFIG=$<CONFIG>"
        "-DINPUT=${CURRENT_TEST_BINARY_DIR}/TwoStream_SEBSurfaceLayerFluxes.i"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DFEXTRACT=${FEXTRACT_EXE}"
        "-DPYTHON_EXE=${ERF_TEST_PYTHON}"
        "-DCHECKER=${TWO_STREAM_SEB_MOISTURE_CHECKER}"
        "-DSTEPS=10"
        "-DDT=1.0"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunTwoStreamSEBMoisture.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 1200
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;radiation"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/dry/simulation.log;${CURRENT_TEST_BINARY_DIR}/wet/simulation.log;${CURRENT_TEST_BINARY_DIR}/mid/simulation.log;${CURRENT_TEST_BINARY_DIR}/veg/simulation.log;${CURRENT_TEST_BINARY_DIR}/bare/simulation.log;${CURRENT_TEST_BINARY_DIR}/bare_f0/simulation.log;${CURRENT_TEST_BINARY_DIR}/bare_dry/simulation.log;${CURRENT_TEST_BINARY_DIR}/tables/simulation.log;${CURRENT_TEST_BINARY_DIR}/checker.log")
endfunction(add_test_two_stream_seb_moisture)

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

    add_test(NAME ${TEST_NAME} COMMAND ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DCONFIG=$<CONFIG>"
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

    add_test(NAME ${TEST_NAME} COMMAND ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DCONFIG=$<CONFIG>"
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

# Stable flow over a steep ridge with Smagorinsky2D (Tests/test_files/Smag2D_Ridge), driven by
# Tests/RunSmag2DRidge.cmake in one of three modes: "steep" (the WRF slope limiter at a time
# step the unlimited closure does not survive), "check" (the diffusive time-step warning fires
# and leaves the answer unchanged) and "limit" (erf.diffusive_dt_limit sets and bounds dt).
function(add_test_smag2d_ridge TEST_NAME MODE PLTFILE)
    set(oneValueArgs "OPTIONS" "CONTROL_OPTIONS" "ALPHA_MIN" "WMAX" "DIFFUSIVE_CFL"
                     "REF_DIFFUSIVE_DT_LO" "REF_DIFFUSIVE_DT_HI")
    cmake_parse_arguments(ADD_TEST_SR "" "${oneValueArgs}" "" ${ARGN})
    set(TEST_FILES_DIR Smag2D_Ridge)
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)

    add_test(NAME ${TEST_NAME} COMMAND ${CMAKE_COMMAND}
        "-DMODE=${MODE}"
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DCONFIG=$<CONFIG>"
        "-DINPUT=${CURRENT_TEST_BINARY_DIR}/Smag2D_Ridge.i"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DFCOMPARE=${FCOMPARE_EXE}"
        "-DFEXTREMA=${FEXTREMA_EXE}"
        "-DPLTFILE=${PLTFILE}"
        "-DOPTIONS=${ADD_TEST_SR_OPTIONS}"
        "-DCONTROL_OPTIONS=${ADD_TEST_SR_CONTROL_OPTIONS}"
        "-DALPHA_MIN=${ADD_TEST_SR_ALPHA_MIN}"
        "-DWMAX=${ADD_TEST_SR_WMAX}"
        "-DDIFFUSIVE_CFL=${ADD_TEST_SR_DIFFUSIVE_CFL}"
        "-DREF_DIFFUSIVE_DT_LO=${ADD_TEST_SR_REF_DIFFUSIVE_DT_LO}"
        "-DREF_DIFFUSIVE_DT_HI=${ADD_TEST_SR_REF_DIFFUSIVE_DT_HI}"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunSmag2DRidge.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 2400
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/steep/simulation.log;${CURRENT_TEST_BINARY_DIR}/control/simulation.log;${CURRENT_TEST_BINARY_DIR}/check_on/simulation.log;${CURRENT_TEST_BINARY_DIR}/check_off/simulation.log;${CURRENT_TEST_BINARY_DIR}/parity.log;${CURRENT_TEST_BINARY_DIR}/limit/simulation.log")
endfunction(add_test_smag2d_ridge)

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

# The wildcard expansion the same scripts rely on to find erf_exec and fcompare in a
# multi-config build tree.  A resolution that picks the wrong binary, or none at all, only
# shows up on Windows, and there as a regression test that fails before it starts, so it is
# tested here on every platform.  Pure CMake, no ERF run, hence the "unit" label.
add_test(ResolveExecutable_SelfTest ${CMAKE_COMMAND}
    "-DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/test_files/ResolveExecutable_SelfTest"
    -P ${PROJECT_SOURCE_DIR}/Tests/ResolveExecutableSelfTest.cmake)
set_tests_properties(ResolveExecutable_SelfTest
    PROPERTIES
    TIMEOUT 60
    PROCESSORS 1
    LABELS "unit")

# The nightly case Exec/RegTests/NoahMP_Ideal keeps a copy of Noah-MP's parameter table so it
# runs in place; fail as soon as that copy stops matching the submodule, so the nightly
# reference never silently tests stale parameters. A file comparison: no Noah-MP build needed.
# Registered only when the submodule is checked out, which every CI job does.
set(ERF_NOAHMP_TABLE "${PROJECT_SOURCE_DIR}/Submodules/Noah-MP/parameters/NoahmpTable.TBL")
if(EXISTS "${ERF_NOAHMP_TABLE}")
  add_test(NoahMP_Ideal_TableMatchesSubmodule ${CMAKE_COMMAND}
      "-DSUBMODULE_TABLE=${ERF_NOAHMP_TABLE}"
      "-DCOPY_TABLE=${PROJECT_SOURCE_DIR}/Exec/RegTests/NoahMP_Ideal/NoahmpTable.TBL"
      -P ${PROJECT_SOURCE_DIR}/Tests/CheckNoahmpTableCopy.cmake)
  set_tests_properties(NoahMP_Ideal_TableMatchesSubmodule
      PROPERTIES
      TIMEOUT 60
      PROCESSORS 1
      LABELS "unit;noahmp")

  # The two-stream balance's copies of Noah-MP's soil and vegetation parameters
  # (erf.radiation.seb_soil_type, seb_vegetation_type) against the same table. Plain Python, no ERF run: every build with the submodule has it.
  if(ERF_TEST_PYTHON)
    add_test(NAME NoahMPSoilTable_MatchesSubmodule
        COMMAND ${ERF_TEST_PYTHON} ${PROJECT_SOURCE_DIR}/Tests/check_noahmp_soil_table.py
                --table ${ERF_NOAHMP_TABLE}
                --header ${PROJECT_SOURCE_DIR}/Source/Radiation/TwoStream/ERF_NoahMPSoilTable.H
                --vegetation-header ${PROJECT_SOURCE_DIR}/Source/Radiation/TwoStream/ERF_NoahMPVegetationTable.H)
    set_tests_properties(NoahMPSoilTable_MatchesSubmodule
        PROPERTIES
        TIMEOUT 60
        PROCESSORS 1
        LABELS "unit;radiation")
  endif()
endif()

# Restart parity: run one deck straight, then to a checkpoint and on from it, and
# require the plotfile at the end to be identical (no gold file). Every run has a
# time limit; the default stays at 600, but an explicit RUN_TIMEOUT is forwarded
# unchanged to each leg and used to size the outer CTest watchdog.
# COMMON_OPTIONS reaches all three legs; RESTART_OPTIONS reaches the restart leg alone,
# which is where anything that changes the decomposition has to go.  ALLOW_DIFF_GRIDS
# lets fcompare compare plotfiles written on different BoxArrays, which a restart leg
# that re-makes the level-0 grids needs and no other restart test should want.
# CHK_NRANKS/RESTART_NRANKS give the checkpoint and restart legs separate rank counts, for
# the case of continuing a run on more ranks than it was written with; the straight leg
# follows the restart's count so the comparison isolates the restart, not the decomposition.
# REQUIRE_LEVEL0_REMAKE makes a regrid test fail if the restart stopped regridding, which
# would otherwise leave it silently comparing an ordinary restart and passing.
function(add_test_restart_parity TEST_NAME TEST_FILES_DIR STEP_CHK STEP_END)
    set(oneValueArgs "COMMON_OPTIONS" "RESTART_OPTIONS" "CHK_NRANKS" "RESTART_NRANKS" "FCOMPARE_RTOL" "FCOMPARE_ATOL" "RUN_TIMEOUT" "DATALOG" "DATALOG_SIGDIGITS" "PLT2DFILE")
    cmake_parse_arguments(ADD_TEST_RP "ALLOW_DIFF_GRIDS;REQUIRE_LEVEL0_REMAKE" "${oneValueArgs}" "" ${ARGN})
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

    add_test(NAME ${TEST_NAME} COMMAND ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DCONFIG=$<CONFIG>"
        "-DINPUT=${CURRENT_TEST_BINARY_DIR}/${TEST_FILES_DIR}.i"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DFCOMPARE=${FCOMPARE_EXE}"
        "-DSTEP_CHK=${STEP_CHK}"
        "-DSTEP_END=${STEP_END}"
        "-DRTOL=${_fcompare_rtol}"
        "-DATOL=${_fcompare_atol}"
        "-DRUN_TIMEOUT=${_run_timeout}"
        "-DCOMMON_OPTIONS=${ADD_TEST_RP_COMMON_OPTIONS}"
        "-DRESTART_OPTIONS=${ADD_TEST_RP_RESTART_OPTIONS}"
        "-DALLOW_DIFF_GRIDS=${ADD_TEST_RP_ALLOW_DIFF_GRIDS}"
        "-DREQUIRE_LEVEL0_REMAKE=${ADD_TEST_RP_REQUIRE_LEVEL0_REMAKE}"
        "-DCHK_NRANKS=${ADD_TEST_RP_CHK_NRANKS}"
        "-DRESTART_NRANKS=${ADD_TEST_RP_RESTART_NRANKS}"
        "-DDATALOG=${ADD_TEST_RP_DATALOG}"
        "-DDATALOG_SIGDIGITS=${ADD_TEST_RP_DATALOG_SIGDIGITS}"
        "-DPLT2DFILE=${ADD_TEST_RP_PLT2DFILE}"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunRestartParity.cmake)
    # The reservation has to cover the widest leg, which need not be NP.
    set(_procs "${NP}")
    foreach(_n "${ADD_TEST_RP_CHK_NRANKS}" "${ADD_TEST_RP_RESTART_NRANKS}")
        if(NOT "${_n}" STREQUAL "" AND _n GREATER _procs)
            set(_procs "${_n}")
        endif()
    endforeach()

    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT ${_ctest_timeout}
        PROCESSORS ${_procs}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;restart-parity"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/straight/simulation.log;${CURRENT_TEST_BINARY_DIR}/restart/checkpoint.log;${CURRENT_TEST_BINARY_DIR}/restart/restart.log;${CURRENT_TEST_BINARY_DIR}/parity.log")
endfunction(add_test_restart_parity)

# Restart abort: run a deck to a checkpoint, restart from it, and require the restart to stop
# with EXPECTED_MESSAGE. For guards that only a restart can reach, which add_test_abort cannot
# test because it runs a single leg and no checkpoint exists yet.  CHK_NRANKS and
# RESTART_NRANKS are separate so a test can reproduce a restart on more ranks than the
# checkpoint was written with, which is its own code path.
function(add_test_restart_abort TEST_NAME TEST_FILES_DIR STEP_CHK EXPECTED_MESSAGE)
    set(oneValueArgs "COMMON_OPTIONS" "RESTART_OPTIONS" "CHK_NRANKS" "RESTART_NRANKS" "RUN_TIMEOUT")
    cmake_parse_arguments(ADD_TEST_RA "" "${oneValueArgs}" "" ${ARGN})
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)

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
        -P ${PROJECT_SOURCE_DIR}/Tests/RunRestartAbort.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT ${_ctest_timeout}
        PROCESSORS ${_procs}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;restart-parity"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/run/checkpoint.log;${CURRENT_TEST_BINARY_DIR}/run/restart.log")
endfunction(add_test_restart_abort)

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

# Field-bounds test: run one deck and require the extrema of a plotfile variable to stay
# within [LO, HI]. The three inflow cases keep a 300 K box at 300 K: a primitive theta from
# a dirichlet_file is multiplied by the ghost density the face prescribes on an Inflow face
# (InflowThetaDensity) and on an inflow_outflow face (InflowOutflowThetaFile), and only the
# face that read the file uses it (InflowThetaFileOtherFace, a second inflow face with its
# own density and theta).
function(add_test_field_bounds TEST_NAME PLTFILE VARIABLE LO HI)
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
        "-DVARIABLE=${VARIABLE}"
        "-DLO=${LO}"
        "-DHI=${HI}"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunFieldBounds.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/simulation.log")
endfunction(add_test_field_bounds)

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
add_test_field_bounds(InflowThetaDensity "plt00010" theta 299.999 300.001)
add_test_field_bounds(InflowOutflowThetaFile "plt00010" theta 299.999 300.001)
add_test_field_bounds(InflowThetaFileOtherFace "plt00010" theta 299.999 300.001)
# Input sponge with a refined patch that does not reach the domain top: ran outside the patch's boxes at start-up
add_test_field_bounds(InputSponge_FinePatch "plt00004" x_velocity 5.9 8.5)
# A refined patch below the inversion takes the coarse level's PBL height (~513 m), not its own
add_test_field_bounds(PBLH_FinePatch "plt2d00004" pblh 450.0 600.0)
# A refined patch aloft (256-768 m) has no column at the ground: it takes level 0's height (~461 m), not zero
add_test_field_bounds(PBLH_PatchAloft "plt2d00004" pblh 400.0 600.0)
# The same patch on one box and on 8 x 8 columns: the coarse heights the patch takes, and the
# length cap and eddy viscosity they set, must not depend on the decomposition
add_test_box_parity(PBLH_FinePatch_BoxParity PBLH_FinePatch "plt00004"
    COMMON_OPTIONS "erf.input_sounding_file=${CMAKE_CURRENT_BINARY_DIR}/test_files/PBLH_FinePatch_BoxParity/sounding_inversion"
    REFERENCE_OPTIONS "amr.max_grid_size=1024"
    SPLIT_OPTIONS "amr.max_grid_size_x=8 amr.max_grid_size_y=8 amr.max_grid_size_z=64")
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

# Regression test for land surface models
function(add_test_lsm TEST_NAME TEST_DIR TEST_EXE)
    set(options )
    set(oneValueArgs "INPUT_SOUNDING" "RUNTIME_OPTIONS")
    set(multiValueArgs PLTFILES "EXTRA_FILES" "LABELS")
    cmake_parse_arguments(ADD_TEST_LSM "${options}" "${oneValueArgs}"
        "${multiValueArgs}" ${ARGN})

    # Check additional external files before creating the test directory.
    foreach(EXTRA_FILE IN LISTS ADD_TEST_LSM_EXTRA_FILES)
        if(NOT EXISTS "${EXTRA_FILE}")
            message(WARNING
                "Skipping LSM test '${TEST_NAME}': extra file does not exist: "
                "'${EXTRA_FILE}'")
            return()
        endif()
    endforeach()

    setup_test()

    set(RUNTIME_OPTIONS "${ADD_TEST_LSM_RUNTIME_OPTIONS}")
    if(NOT "${ADD_TEST_LSM_INPUT_SOUNDING}" STREQUAL "")
      string(APPEND RUNTIME_OPTIONS "erf.input_sounding_file=${CURRENT_TEST_BINARY_DIR}/${ADD_TEST_LSM_INPUT_SOUNDING}")
    endif()

    # Copy any additional external files needed to the test directory
    foreach(EXTRA_FILE IN LISTS ADD_TEST_LSM_EXTRA_FILES)
        message(DEBUG " -- Copying extra file '${EXTRA_FILE}' to test directory '${CURRENT_TEST_BINARY_DIR}'")
        file(COPY "${EXTRA_FILE}" DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    endforeach()

    if (ADD_TEST_LSM_LABELS)
        set(test_labels "")
        foreach(LABEL ${ADD_TEST_LSM_LABELS})
            list(APPEND test_labels "${LABEL}")
        endforeach()
    else()
        set(test_labels "regression")
    endif()

    if(WIN32)
        set(TEST_EXE "${CMAKE_BINARY_DIR}/Exec/${TEST_DIR}/*/${TEST_EXE}.exe")
    else()
        set(TEST_EXE "${CMAKE_BINARY_DIR}/Exec/${TEST_DIR}/${TEST_EXE}")
    endif()

    set(FCOMPARE_TOLERANCE "-r ${ERF_TEST_FCOMPARE_RTOL} --abs_tol ${ERF_TEST_FCOMPARE_ATOL}")
    set(FCOMPARE_FLAGS "--abort_if_not_all_found -a ${FCOMPARE_TOLERANCE}")

    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUNTIME_OPTIONS} > ${TEST_NAME}.log")
    # These tests are gated on external input data, so their reference plotfiles are not carried
    # in Tests/ERFGoldFiles either.  Compare against whichever ones the configured gold directory
    # actually has, and run the rest to completion as smoke tests rather than failing every run on
    # a reference that is not there -- fcompare is invoked with --abort_if_not_all_found.
    foreach(PLTFILE ${ADD_TEST_LSM_PLTFILES})
        if(EXISTS "${PLOT_GOLD}/${PLTFILE}")
            set(test_command "${test_command} && ${MPI_FCOMP_COMMANDS} ${FCOMPARE_EXE} ${FCOMPARE_FLAGS} ${PLOT_GOLD}/${PLTFILE} ${CURRENT_TEST_BINARY_DIR}/${PLTFILE}")
        else()
            message(STATUS
                " -- LSM test '${TEST_NAME}': no gold file '${PLOT_GOLD}/${PLTFILE}', "
                "running without a plotfile comparison")
        endif()
    endforeach()
    message(DEBUG "TEST COMMAND FOR '${TEST_NAME}': ${test_command}")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 5400
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "${test_labels}"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log"
    )
endfunction(add_test_lsm)

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
  # MYNN25 carries a prognostic TKE, and its source terms read the stage state with a
  # VERTICAL stencil: AddTurbKESources forms d(theta_v)/dz from k-1 and k+1 while
  # erf_slow_rhs_post writes that same state (new_cons = cur_cons) later in the same
  # MFIter iteration.  Tiled, a tile's write landed before the next tile's read, so a
  # cell beside a tile boundary saw the updated state on one side and the stage-entry
  # state on the other and the answer depended on fabarray.mfiter_tile_size.  The
  # harness tiles 1024000 8 8, which splits in z -- the direction this stencil reaches.
  # The Closure_BoxParity entry for MYNN25 below does not cover it: that one splits
  # 8 8 1024, so it tiles in x and y only and a k+-1 hazard is invisible to it.
  # KE is added to the 3D plotfile because it is the field that breaks first: it moved
  # by ~2e-3 relative after a single step, which the MYNN length scales then amplified
  # into Kmv/Khv and from there into the solution.
  add_test_tiling_parity(ABL_MYNN25_Tiling   ABL_MRF_Tiling "00010" "00010"
      RUNTIME_OPTIONS "erf.pbl_type=MYNN25 erf.most.pblh_calc=MYNN25 'erf.plot_vars_1=density x_velocity y_velocity z_velocity pressure theta Kmv Khv Lturb KE'"
      VARYING_3D "Lturb Kmv KE" VARYING_2D "pblh u_star")
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
  # Two levels. The same column physics must hold on the fine level, which runs
  # its own sweep: CHECK_LEVELS 0 1 runs the vertical-structure assertions on
  # both, so a fine level left at the allocation's zero heating fails the
  # "qsrc_sw is zero everywhere" check. The refinement patch is tagged (not an
  # explicit erf.boxN), so amr.refine_whole_domain_dir = 2 is what makes it span
  # z -- which is also the remediation the model's abort recommends.
  add_test_two_stream_radiation(TwoStream_ColumnHeating_TwoLevel "plt00002"
                                CHECK_LEVELS 0 1)
endif()

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
        set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}")
    else()
        unset(MPI_COMMANDS)
    endif()

    resolve_test_exe("" "erf_exec" TEST_EXE)

    set(test_log "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log")
    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${INPUT_FILE} max_step=1 erf.plot_int_1=-1 erf.plot_int_2=-1 erf.check_int=-1 amrex.call_addr2line=0 ${RUNTIME_OPTIONS} 2>&1 | tee ${test_log}")

    add_test(${TEST_NAME} ${test_command})
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

if(ERF_ENABLE_MPI AND NOT WIN32)
  # A shallow nest -- a fine level that stops below the domain top -- has no complete
  # column, so the sweep cannot run on it. That is a supported configuration, not an
  # error: advance_radiation interpolates the level's heating rates and fluxes from its
  # parent, the same route RRTMGP takes for a nested patch. Same deck as the two-level
  # case with the grid guarantee switched off, so the patch really is shallow (it comes
  # out as k = 0..7 of a 32-cell domain). The checker detects the nested level from the
  # data and asserts what remains meaningful there -- finite, non-negative, and not
  # identically zero, which is exactly what fails if the interpolation never happened.
  add_test_two_stream_radiation(TwoStream_NestedPatch "plt00002"
                                CHECK_LEVELS 0 1
                                DIAG_LEVELS 0)

  # The prognostic surface energy balance now runs on a refined hierarchy: every level
  # evolves its own force-restore surface temperature, and ERF averages t_sfc and q_sfc
  # down so the levels agree about the ground they share. (This replaces the abort test for
  # the old single-level restriction, which this branch lifts.)
  #
  # The parity checker asserts the relation average_down establishes exactly -- a coarse
  # cell holds the mean of the fine cells above it -- so its tolerance is a round-off
  # bound. Dropping the average-down moves the worst cell to ~3e-5 K, four orders of
  # magnitude above it, from the first output onward. Domain means would NOT catch this:
  # average_down is mean-preserving, so the two levels' means agree either way.
  if(ERF_TEST_PYTHON)
    add_test_two_stream_radiation(TwoStream_PrognosticSEBMultiLevel "plt00006"
                                  CHECK_LEVELS 0 1
                                  SEB_PARITY_PLOTFILE "plt2d00006")
  endif()

  # The complement of the case above: a SHALLOW nest, whose boxes do not span the domain
  # in z. Such a level never sweeps, so its surface state stays frozen at what
  # fill_seb_from_coarse wrote, and averaging that down would pin level 0's surface under
  # the patch at its level-creation value -- a regression against the old level-0-only
  # behaviour. post_timestep skips the transfer for such a level; this asserts every cell
  # of level 0 moved away from erf.rad_t_sfc. With the guard removed those cells sit at
  # exactly 300.0 instead of 300.0476.
  if(ERF_TEST_PYTHON)
    add_test_two_stream_radiation(TwoStream_PrognosticSEBShallowNest "plt00006"
                                  CHECK_LEVELS 0
                                  DIAG_LEVELS 0
                                  SEB_PARITY_PLOTFILE "plt2d00006"
                                  SEB_EVOLVED_FROM "300.0")
  endif()

  # A level created MID-RUN must start from the surface its parent has reached.
  # TwoStream_PrognosticSEBMultiLevel cannot see that: its fine level exists from t = 0, when
  # both levels are uniform at erf.rad_t_sfc and interpolating from the parent is a no-op.
  # Here the tagging switches on only after level 0 has run alone for 10 steps, so the regrid
  # at step 11 builds level 1 over a surface that has drifted ~0.08 K. The checker asserts
  # level 1 averages to the parent's surface extrapolated one step on, to a quarter of one
  # step's change; without ERF::fill_seb_from_coarse at level creation it misses by the whole
  # drift (7.9e-2 against a 2.0e-3 tolerance, where the transfer gives 7.9e-5).
  #
  # The runner's 1-rank vs 2-rank comparison matters here too: it is what showed the
  # interpolation reading a periodic halo that fill_seb_from_coarse had clamped, which left
  # the new level's edge columns dependent on the box layout.
  if(ERF_TEST_PYTHON)
    add_test_two_stream_radiation(TwoStream_PrognosticSEBLateLevel "plt00011"
                                  CHECK_LEVELS 0 1
                                  SEB_PARITY_PLOTFILE "plt2d00011"
                                  SEB_CREATED_FROM "plt2d00000" "plt2d00009" "plt2d00010")
  endif()

  # A regrid that MOVES an existing fine level must keep the surface it has evolved. The
  # cases above tag a fixed region, so the fine BoxArray never changes and RemakeLevel never
  # runs. Here the refinement box moves at step 11, keeping fine x-cells 4-7, adding 8-11 and
  # dropping 0-3. The checker asserts the kept cells carry on from the fine surface (block
  # means and sub-coarse deviations, each extrapolated one step) and the added cells start
  # from the parent. Each piece of RemakeLevel's restore fails its own assertion when
  # removed: the whole block (means off by 7.9e-2 against 2.0e-3), only the copy of the
  # retained values (deviations off by 2.2e-4 against 3.8e-5), only the interpolation from
  # the parent (added cells off by 7.9e-2 against 2.0e-3). The correct run sits at 7e-6,
  # 3e-6 and 4e-5.
  if(ERF_TEST_PYTHON)
    add_test_two_stream_radiation(TwoStream_PrognosticSEBRegrid "plt00011"
                                  CHECK_LEVELS 0 1
                                  SEB_PARITY_PLOTFILE "plt2d00011"
                                  SEB_REGRIDDED_FROM "plt2d00000" "plt2d00009" "plt2d00010")
  endif()

  # The force-restore surface state is checkpointed per level. A level whose copy is
  # never written, or is written and never read back, restarts from the erf.rad_t_sfc
  # scalar instead of the surface the run had reached -- and nothing in the 3D plotfile
  # would show it, because the 3D plotfile carries no surface fields. So this case
  # selects the surface state as 2D output and compares plt2d as well as plt; without
  # PLT2DFILE the check would pass on a surface temperature that reset to its default.
  add_test_restart_parity(TwoStream_PrognosticSEB_Restart TwoStream_PrognosticSEBRestart 3 6
                          PLT2DFILE "plt2d00006")

  # The prognostic surface energy balance must remove from the ground the sensible and
  # latent heat the surface layer puts into the air. The deck runs twice, with
  # erf.radiation.seb_turbulent_flux_source = surface_layer and = defaults; the checker
  # asserts the balance's seb_hfx/seb_lh equal the surface layer's sensible_heat_flux/
  # latent_heat_flux at every step (to 1e-10 relative, with |H| and |LE| above 1 W/m^2),
  # that the defaults leg keeps the constants, and that the skin ends cooler by
  # sum(dt (H + LE)) / C_s to 5 %. With the balance reading the defaults in both legs
  # (the behaviour before this test) seb_hfx is 0 against a sensible_heat_flux of tens
  # of W/m^2, and the two skins agree.
  if(ERF_TEST_PYTHON)
    add_test_two_stream_seb_flux_source(TwoStream_SEBSurfaceLayerFluxes)
    # The same on two levels over a ridge, level 1 created at step 7 over the middle half:
    # every level and column checked, with per-column surface-layer fluxes (H spread
    # across the coarse columns at least 0.5 W/m^2, so the columns are told apart). No
    # subcycling, so a 0.5 s step keeps the fine level's acoustic substeps stable.
    add_test_two_stream_seb_flux_source(TwoStream_SEBSurfaceLayerFluxesMultiLevel
                                        DT 0.5
                                        CHECKER_OPTIONS "--multilevel --min-spread 0.5")
    # The same deck with the soil moisture coupled as well: the surface mixing ratio blends
    # q_sat and the air's by beta, vegetation lowers the evaporation, and the balance's soil
    # loses the water the air gains.
    add_test_two_stream_seb_moisture(TwoStream_SEBSoilMoisture)
  endif()

  # Two-stream radiation feeding Noah-MP on two levels (see the function above). Noah-MP
  # needs a parallel NetCDF build, which no CI job has, so this runs where one exists.
  if(ERF_ENABLE_NOAHMP AND ERF_TEST_PYTHON)
    add_test_two_stream_noahmp_levels(TwoStream_NoahMPLevels)
  endif()

  # With seb_turbulent_flux_source = surface_layer (the default) the balance takes H from
  # the surface layer wherever its flux field exists -- including an adiabatic surface layer,
  # whose flux is zero -- so a nonzero erf.radiation.seb_hfx_default in the deck is not used.
  # That must be said at start-up rather than happen silently. One step of the deck above.
  add_test_abort(TwoStream_SEBDefaultReplacedWarning
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_files/TwoStream_SEBSurfaceLayerFluxes
                 TwoStream_SEBSurfaceLayerFluxes.i
                 "seb_hfx_default = 10 is not used"
                 "erf.radiation.seb_hfx_default=10")
endif()
add_test_plotfile_header(Plotfile3D_TwoStreamHeatingSelection "" "erf_exec" "plt00000")

# The two-stream radiation tests stay registered but do not run: DISABLED keeps
# them listed by ctest (reported as "Not Run (Disabled)") without executing them.
# Remove this block to re-enable them. TwoStream_SEBSoilMoisture is not on the list:
# its eight legs of 10 steps on a 4 x 4 x 16 grid take about a minute in Debug.
foreach(_two_stream_test IN ITEMS
    TwoStream_ColumnHeating
    TwoStream_ColumnHeating_Terrain
    TwoStream_ColumnHeating_TwoLevel
    TwoStream_NoahMPLevels
    TwoStream_NestedPatch
    TwoStream_PrognosticSEBMultiLevel
    TwoStream_PrognosticSEBShallowNest
    TwoStream_PrognosticSEBLateLevel
    TwoStream_PrognosticSEBRegrid
    TwoStream_PrognosticSEB_Restart
    TwoStream_SEBSurfaceLayerFluxes
    TwoStream_SEBSurfaceLayerFluxesMultiLevel
    TwoStream_SEBDefaultReplacedWarning
    Plotfile3D_TwoStreamHeatingSelection)
  if(TEST ${_two_stream_test})
    set_tests_properties(${_two_stream_test} PROPERTIES DISABLED TRUE)
  endif()
endforeach()

add_test_0(CouetteFlow_x                     "" "erf_exec" "plt00050" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_0(CouetteFlow_y                     "" "erf_exec" "plt00050" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_0(PoiseuilleFlow_x                  "" "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_0(PoiseuilleFlow_y                  "" "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_0(InitSoundingIdeal_stationary      "" "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")
add_test_0(Deardorff_stationary              "" "erf_exec" "plt00010" RUNTIME_OPTIONS "erf.vert_implicit=false ")

# LSM tests are gated because they require large input files (wrfinput, etc)
# and can take several hours to run.
if(ERF_TEST_ENABLE_EXTRA_LSM_TESTS)
    # CASS case with external-driven radiation fluxes to SLM (no RRTMGP).
    # The listed plotfiles are compared only if the configured gold directory has them.
    add_test_lsm(SLM_CASS_SAMRadiation            "" "erf_exec"
                                                  LABELS "slm"
                                                  EXTRA_FILES "${CMAKE_SOURCE_DIR}/Tests/test_files/SLM_CASS_SAMRadiation/sounding_cass_interpolated"
                                                              "${CMAKE_SOURCE_DIR}/Tests/test_files/SLM_CASS_SAMRadiation/lsf_cass"
                                                              "${ERF_TEST_EXTRA_FILES_DIRECTORY}/CASS_32x32x156_50m_50m_1s_rad_coszrs_combined.nc"
                                                  PLTFILES "plt34500"
                                                           "plt_lsm_34500"
                                                           "plt_lsm_2D_34500")

    # LBA case using RRTMGP radiation
    add_test_lsm(SLM_LBA_RRTMGP                   "" "erf_exec"
                                                  LABELS "slm" "manual"
                                                  EXTRA_FILES "${CMAKE_SOURCE_DIR}/Tests/test_files/SLM_LBA_RRTMGP/snd_lba"
                                                              "${ERF_TEST_EXTRA_FILES_DIRECTORY}/rrtmgp-gas-sw-g112.nc"
                                                              "${ERF_TEST_EXTRA_FILES_DIRECTORY}/rrtmgp-gas-lw-g128.nc"
                                                              "${ERF_TEST_EXTRA_FILES_DIRECTORY}/rrtmgp-cloud-optics-coeffs-sw.nc"
                                                              "${ERF_TEST_EXTRA_FILES_DIRECTORY}/rrtmgp-cloud-optics-coeffs-lw.nc")

    # AWAKEN case testing the fully coupled real pathway
    add_test_lsm(SLM_AWAKEN                       "" "erf_exec"
                                                  LABELS "slm" "manual"
                                                  EXTRA_FILES "${ERF_TEST_EXTRA_FILES_DIRECTORY}/SLM_AWAKEN/wrfinput_d01"
                                                              "${ERF_TEST_EXTRA_FILES_DIRECTORY}/SLM_AWAKEN/wrfbdy_d01"
                                                              "${ERF_TEST_EXTRA_FILES_DIRECTORY}/rrtmgp-gas-sw-g112.nc"
                                                              "${ERF_TEST_EXTRA_FILES_DIRECTORY}/rrtmgp-gas-lw-g128.nc"
                                                              "${ERF_TEST_EXTRA_FILES_DIRECTORY}/rrtmgp-cloud-optics-coeffs-sw.nc"
                                                              "${ERF_TEST_EXTRA_FILES_DIRECTORY}/rrtmgp-cloud-optics-coeffs-lw.nc")
endif()

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

# The balance on two levels (Tests/RunIBSEBRefinedLevels.cmake): a cube on level 1 and a
# tower outside it, whose shadow the cube's level-1 faces must find through level 0's
# column map; the level-1 face dumps, and the refinement-box script's verdicts on the
# deck's box and on one cut through the cube. Then the start-up aborts: a refined level
# whose edge cuts a building (one level 0 resolves, and one only level 1 does), the faces on their own sun under the two-stream radiation
# (sun_mode = solar, or a fixed two-stream sun the faces do not share), and
# sun_mode = two_stream without it; IBSEB_TwoStreamSunRun runs that deck and checks the sun.
function(add_test_ibseb_refined_levels TEST_NAME)
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)
    add_test(NAME ${TEST_NAME} COMMAND ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DCONFIG=$<CONFIG>"
        "-DINPUT=${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DPYTHON_EXE=${ERF_TEST_PYTHON}"
        "-DCHECKER=${CMAKE_CURRENT_SOURCE_DIR}/check_ibseb_refined_levels.py"
        "-DBOX_SCRIPT=${PROJECT_SOURCE_DIR}/Exec/CanonicalTests/SEB/ibseb_refinement_box.py"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunIBSEBRefinedLevels.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 1200
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/cube_only/simulation.log;${CURRENT_TEST_BINARY_DIR}/both/simulation.log;${CURRENT_TEST_BINARY_DIR}/checker.log")
endfunction(add_test_ibseb_refined_levels)

# The faces on the two-stream sun in a run (Tests/RunIBSEBTwoStreamSun.cmake).
function(add_test_ibseb_two_stream_sun TEST_NAME TEST_FILES_DIR)
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)
    add_test(NAME ${TEST_NAME} COMMAND ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DCONFIG=$<CONFIG>"
        "-DINPUT=${CURRENT_TEST_BINARY_DIR}/${TEST_FILES_DIR}.i"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DPYTHON_EXE=${ERF_TEST_PYTHON}"
        "-DCHECKER=${CMAKE_CURRENT_SOURCE_DIR}/check_ibseb_two_stream_sun.py"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunIBSEBTwoStreamSun.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 1200
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/run/simulation.log;${CURRENT_TEST_BINARY_DIR}/checker.log")
endfunction(add_test_ibseb_two_stream_sun)

if(ERF_ENABLE_MPI AND NOT WIN32)
  if(NOT "${ERF_TEST_PYTHON}" STREQUAL "")
    add_test_ibseb_refined_levels(IBSEB_RefinedLevels)
    add_test_ibseb_two_stream_sun(IBSEB_TwoStreamSunRun IBSEB_TwoStreamSun)
  endif()
  add_test_abort(IBSEB_RefinedLevelCutsBuilding
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_files/IBSEB_RefinedLevels
                 IBSEB_RefinedLevels.i
                 "lie on the edge of level 1"
                 "erf.city.in_box_lo=160.0 180.0")
  # A 4 m block across level 1's edge that only level 1 resolves (5 m cells there, 10 m
  # below): the level-0 check cannot see it, IBFaceSet::build() must.
  add_test_abort(IBSEB_RefinedLevelCutsLowBuilding
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_files/IBSEB_RefinedLevels
                 IBSEB_RefinedLevels.i
                 "cell faces against solid cells outside its grids"
                 "amr.ref_ratio_vect=2 2 2 erf.buildings_file_name=cube_and_low_block_10m.txt")
  # A report file from an earlier version (20 columns, no building_level0): a restart
  # appending to it must stop rather than write rows its header does not describe.
  add_test_abort(IBSEB_ReportOldColumns
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_files/IBSEB_RefinedLevels
                 IBSEB_RefinedLevels.i
                 "has another set of columns than this version writes"
                 "erf.ibseb.csv_file=old_ibseb_buildings.csv")
  add_test_abort(IBSEB_TwoStreamSunSolar
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_files/IBSEB_TwoStreamSun
                 IBSEB_TwoStreamSun.i
                 "the faces need erf.ibseb.sun_mode = two_stream"
                 "erf.ibseb.sun_mode=solar")
  add_test_abort(IBSEB_TwoStreamSunWithoutTwoStream
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_files/IBSEB_TwoStreamSun
                 IBSEB_TwoStreamSun.i
                 "sun_mode = two_stream needs erf.radiation_model = TwoStream"
                 "erf.radiation_model=None")
  add_test_abort(IBSEB_TwoStreamSunNoDate
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_files/IBSEB_TwoStreamSun
                 IBSEB_TwoStreamSunNoDate.i
                 "sun_mode = two_stream: the two-stream sun follows the calendar and no start date"
                 "")
  add_test_abort(IBSEB_TwoStreamSunSolarInputs
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_files/IBSEB_TwoStreamSun
                 IBSEB_TwoStreamSun.i
                 "latitude_deg is not used with erf.ibseb.sun_mode = two_stream"
                 "erf.ibseb.latitude_deg=40.0")
  add_test_abort(IBSEB_TwoStreamSunFixed
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_files/IBSEB_TwoStreamSun
                 IBSEB_TwoStreamSun.i
                 "the two-stream sun is fixed"
                 "erf.fixed_solar_zenith_angle=0.5")
endif()

# The cube through a checkpoint at step 17, to step 40. A restart rebuilds the
# immersed forcing's blanking rather than reading it back; it rebuilt it without
# clearing the almost-fluid cells (eb2.small_volfrac), so the restarted run forced
# cells the straight run leaves alone (2e-4 in terrain_IB_mask, 6e-5 m/s in u, 8e-6 K
# in the face skin temperatures). The deck plots no velocities, so they are added
# here. Each test launches the deck three times, so both pass a RUN_TIMEOUT: 600 s
# is the whole budget of the single 40-step IBSEB_Cube run (and half of what the
# two-level IBSEB_RefinedLevels test gets for two runs of two steps, the view
# factors being what costs), and the default watchdog would cut the legs off
# looking like a parity failure rather than a timeout.
# IBSEB_RefinedLevels_Restart does the same on the two levels of the
# IBSEB_RefinedLevels deck (a cube on level 1, a tower outside it); its face dumps go
# to a plain file name, since the deck's faces/ directory does not exist in the
# runner's legs. It compares to 1e-8 relative, not zero: the restart leg writes a
# plotfile at the restart step and the straight leg does not, and on two levels
# writing a plotfile changes the solution at round-off (issue 4224; 5e-15 relative in theta
# without buildings, 1.7e-10 relative in w here by step 20); with the plotfiles at
# the same steps in all legs the restart is bit-exact. terrain_IB_mask still shows
# the uncleared blanking (2.6e-3); 1e-8 relative is about 3e-6 K on the skin temperatures.
# The runner is a cmake -P script (MPI, not Windows).
if(ERF_ENABLE_MPI AND NOT WIN32)
add_test_restart_parity(IBSEB_Cube_Restart IBSEB_Cube 17 40
    COMMON_OPTIONS "erf.plot_vars_1=density x_velocity y_velocity z_velocity theta terrain_IB_mask ibseb_nfaces ibseb_tskin ibseb_sw_abs ibseb_lw_net ibseb_H ibseb_G"
    RUN_TIMEOUT 900
    FCOMPARE_RTOL "0.0" FCOMPARE_ATOL "0.0")
add_test_restart_parity(IBSEB_RefinedLevels_Restart IBSEB_RefinedLevels 7 20
    COMMON_OPTIONS "erf.ibseb.dump_faces_file=faces_set erf.plot_vars_1=density x_velocity y_velocity z_velocity theta terrain_IB_mask ibseb_nfaces ibseb_tskin ibseb_sw_abs ibseb_lw_net ibseb_H ibseb_G"
    RUN_TIMEOUT 1200
    FCOMPARE_RTOL "1.0e-8" FCOMPARE_ATOL "0.0")
endif()
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
# moves the answer at all is a bug.  The two decks below cover the dry AMR path and the
# surface-layer path.
#
# AMR, dry: the one of the two whose station resolves to level 1, so it is the case
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

#=============================================================================
# Input sounding on refined levels (#4143): each level samples the sounding, and
# the large-scale forcing profiles, at its own cell centres
#=============================================================================
# The two sounding files differ only in their surface line, which lies below every cell
# centre, so the initial states must be identical on both levels.
add_test_option_parity(InputSounding_FineLevelInit InputSoundingFineLevels "plt00000"
    OFF_OPTIONS "erf.input_sounding_file=../input_sounding_sfc300"
    ON_OPTIONS  "erf.input_sounding_file=../input_sounding_sfc303")

# Nudging towards the sounding the run started from adds exactly zero on every level:
# theta here (the wind is not nudged, since the sheared wind is advected by the w that the
# theta profile sets going) ...
add_test_option_parity(InputSounding_FineLevelNudging InputSoundingFineLevels "plt00001"
    COMMON_OPTIONS "max_step=1 erf.nudging_u=false"
    OFF_OPTIONS "erf.nudging_from_input_sounding=false"
    ON_OPTIONS  "erf.nudging_from_input_sounding=true erf.tau_nudging=0.5")

# ... and u, v through the momentum sources, over a constant theta so that nothing moves.
add_test_option_parity(InputSounding_FineLevelWindNudging InputSoundingFineLevels "plt00001"
    COMMON_OPTIONS "max_step=1 erf.input_sounding_file=../input_sounding_wind"
    OFF_OPTIONS "erf.nudging_from_input_sounding=false"
    ON_OPTIONS  "erf.nudging_from_input_sounding=true erf.tau_nudging=0.5")

# Large-scale forcing that relaxes the wind towards the sounding's own wind, with no
# tendencies or subsidence, adds exactly zero on every level, refined in z included.
add_test_option_parity(InputSounding_FineLevelLSF InputSoundingFineLevels "plt00001"
    COMMON_OPTIONS "max_step=1 erf.input_sounding_file=../input_sounding_wind"
    OFF_OPTIONS "erf.large_scale_forcing=false"
    ON_OPTIONS  "erf.nudging_from_input_sounding=true erf.large_scale_forcing=true erf.large_scale_forcing_file=../lsf_zero_tendency erf.forcing_timescale=0.5")

#=============================================================================
# Terrain: decomposition and station output over a hill
#=============================================================================

# Run a deck and check its station series (see Tests/RunStationSeries.cmake):
# MODE single runs it once and hands CHECKS, separated by '|', to the checker;
# analytic compares the series STATION with a closed-form answer, approach
# compares it between runs with and without nudging (OFF_OPTIONS), and abort
# requires a start-up abort containing EXPECTED_MESSAGE.
function(add_test_station_series TEST_NAME TEST_FILES_DIR MODE)
    set(oneValueArgs "RUNTIME_OPTIONS" "CHECKS" "OFF_OPTIONS" "STATION" "EXPECTED_MESSAGE" "LABELS")
    # PARSE_ARGV takes the arguments from ARGV verbatim, so a value that itself
    # holds a ';' (LABELS "regression;obs-nudging") stays one argument.  Passing
    # an unquoted ${ARGN} instead would split it and drop all but the first word.
    cmake_parse_arguments(PARSE_ARGV 3 ADD_TEST_SS "" "${oneValueArgs}" "")
    setup_test()
    resolve_test_exe("" "erf_exec" TEST_EXE)
    set(test_log "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log")
    set(_nranks ${NP})
    if("${MODE}" STREQUAL "abort")
        set(_nranks 1)
    endif()
    set(_labels "regression;station")
    if(NOT "${ADD_TEST_SS_LABELS}" STREQUAL "")
        set(_labels "${ADD_TEST_SS_LABELS}")
    endif()
    # The checker is named by its target file, which is exact under every
    # generator; the ERF executable is resolved by the runner (it may carry a
    # wildcard for the config subdirectory on Windows)
    add_test(NAME ${TEST_NAME} COMMAND ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${_nranks}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DCONFIG=$<CONFIG>"
        "-DINPUT=${CURRENT_TEST_BINARY_DIR}/${TEST_FILES_DIR}.i"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DLOG=${test_log}"
        "-DMODE=${MODE}"
        "-DCHECKER=$<TARGET_FILE:erf_station_series_check>"
        "-DCHECKS=${ADD_TEST_SS_CHECKS}"
        "-DSTATION=${ADD_TEST_SS_STATION}"
        "-DRUNTIME_OPTIONS=${ADD_TEST_SS_RUNTIME_OPTIONS}"
        "-DOFF_OPTIONS=${ADD_TEST_SS_OFF_OPTIONS}"
        "-DEXPECTED_MESSAGE=${ADD_TEST_SS_EXPECTED_MESSAGE}"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunStationSeries.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS ${_nranks}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "${_labels}"
        ATTACHED_FILES_ON_FAIL "${test_log};${test_log}.off;${test_log}.checker")
endfunction(add_test_station_series)

# The parity drivers run each leg in a subdirectory, so the sounding is named by
# absolute path.
function(terrain_hill_files TEST_NAME OUT_VAR)
    set(${OUT_VAR} "erf.input_sounding_file=${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME}/input_sounding" PARENT_SCOPE)
endfunction()

# Two levels over a 100 m hill on a terrain-fitted mesh.  The zero-gradient
# condition below the mesh is corrected by the terrain slope times a lateral
# gradient that each box can only take one-sided in its outermost ghost
# columns, so the copies of a ghost cell below the ground differed between
# boxes, and the refined level, which interpolates from coarse ghost cells,
# moved with the decomposition (x-velocity 9.4e-6 apart on level 1 after 20
# steps).  The copies are now made to agree (BelowGroundGhostSync), and the
# answer must not depend on the decomposition.
terrain_hill_files(Terrain2Lev_Hill_BoxParity _hill_files)
add_test_box_parity(Terrain2Lev_Hill_BoxParity TerrainHill "plt00020"
    COMMON_OPTIONS "${_hill_files}"
    REFERENCE_OPTIONS "amr.max_grid_size=1024"
    SPLIT_OPTIONS "amr.max_grid_size_x=8 amr.max_grid_size_y=8 amr.max_grid_size_z=64"
    DATALOG "Output_Stations/gate.dat"
    DATALOG_SIGDIGITS 10)

# The same with the refined level starting 100 m up, so that none of its boxes
# reaches the ground: the synchronisation has no cells below the ground to own
# on that level and must leave it alone (a first version copied from boxes
# that do not reach the bottom, and the run aborted in Debug).
terrain_hill_files(Terrain2Lev_HillAloft_BoxParity _hill_files)
add_test_box_parity(Terrain2Lev_HillAloft_BoxParity TerrainHill "plt00020"
    COMMON_OPTIONS "${_hill_files} erf.box1.in_box_lo=400.0 200.0 100.0"
    REFERENCE_OPTIONS "amr.max_grid_size=1024"
    SPLIT_OPTIONS "amr.max_grid_size_x=8 amr.max_grid_size_y=8 amr.max_grid_size_z=64")

# Station output over terrain: a height above the local terrain is measured
# from the ground at the station, so a series asked for 40 m above the terrain
# on the flank of the hill must be the series asked for at 93.08 m above z = 0
# (the ground is 53.08 m up there).  They agree to 1.4e-5 m/s on the fitted
# mesh and 1.2e-4 m/s over the immersed hill, and to 3e-6 K.  The previous
# code fails both: it interpolated the ground from cell averages between cell
# centres, 3.4 m too low on the flank (0.047 m/s away on the fitted mesh), and
# with immersed-forcing terrain it measured the height from the flat bottom of
# the mesh, so the series was the flow inside the hill (3.8 m/s away).
add_test_station_series(StationSampling_FittedTerrain TerrainHill single
    RUNTIME_OPTIONS "amr.max_level=0 erf.station_names=mast mastabs erf.mastabs.field=x_velocity y_velocity theta erf.mastabs.x=700.0 erf.mastabs.y=400.0 erf.mastabs.height_abs=93.08"
    CHECKS "equal a=@RUN@/Output_Stations/mast.dat:2 b=@RUN@/Output_Stations/mastabs.dat:2 tol=0.005|equal a=@RUN@/Output_Stations/mast.dat:4 b=@RUN@/Output_Stations/mastabs.dat:4 tol=0.001")

add_test_station_series(StationSampling_ImmersedTerrain TerrainHill single
    RUNTIME_OPTIONS "amr.max_level=0 erf.terrain_type=ImmersedForcing erf.immersed_forcing_substep=true eb2.small_volfrac=0.005 erf.station_names=mast mastabs erf.mastabs.field=x_velocity y_velocity theta erf.mastabs.x=700.0 erf.mastabs.y=400.0 erf.mastabs.height_abs=93.08"
    CHECKS "equal a=@RUN@/Output_Stations/mast.dat:2 b=@RUN@/Output_Stations/mastabs.dat:2 tol=0.005|equal a=@RUN@/Output_Stations/mast.dat:4 b=@RUN@/Output_Stations/mastabs.dat:4 tol=0.001")

# The same hill as immersed-forcing terrain through a checkpoint at step 7, to step
# 20, on one level and on the deck's two. A restart rebuilt the blanking without
# clearing the almost-fluid cells of the hill's tails (eb2.small_volfrac), so the
# wall law acted there after the restart only: 0.2 m/s in u at step 20 on one level.
# The one-level leg compares at zero tolerance. The two-level one compares to 1e-8
# relative, for the same reason as IBSEB_RefinedLevels_Restart above: the restart leg
# writes a plotfile at the restart step and the straight leg does not, and on two
# levels writing a plotfile moves the solution at round-off (issue 4224).
# The runner is a cmake -P script (MPI, not Windows).
if(ERF_ENABLE_MPI AND NOT WIN32)
add_test_restart_parity(ImmersedTerrain_Hill_Restart TerrainHill 7 20
    COMMON_OPTIONS "amr.max_level=0 erf.terrain_type=ImmersedForcing erf.immersed_forcing_substep=true eb2.small_volfrac=0.005 erf.plot_vars_1=density x_velocity y_velocity z_velocity theta terrain_IB_mask"
    FCOMPARE_RTOL "0.0" FCOMPARE_ATOL "0.0")
add_test_restart_parity(ImmersedTerrain_Hill_TwoLevel_Restart TerrainHill 7 20
    COMMON_OPTIONS "erf.terrain_type=ImmersedForcing erf.immersed_forcing_substep=true eb2.small_volfrac=0.005 erf.plot_vars_1=density x_velocity y_velocity z_velocity theta terrain_IB_mask"
    FCOMPARE_RTOL "1.0e-8" FCOMPARE_ATOL "0.0")
endif()

# A restart that re-makes the level-0 grids (erf.regrid_level_0_on_restart, and the
# same branch ERF::restart takes on its own when the checkpoint has fewer level-0
# boxes than there are ranks).  The remake used to copy the old state onto the new
# grids offering its ghost cells as a source, and ReadCheckpointFile leaves those at
# bogus_large_value, so valid cells of the new grids came out of the remake holding
# 1e150 and the first estTimeStep trapped on the cast of fixed_dt/dt_sub_max (issue
# 4225).  RESTART_OPTIONS, not COMMON_OPTIONS: max_grid_size must move on the restart
# leg alone, or this becomes a box-parity test.  max_grid_size_z is left spanning the
# column because the vertical diffusion is implicit and define_column_kextent refuses
# a split column.  Terrain2Lev_Hill_BoxParity already holds this deck decomposition
# independent under the same 8/8/64 split; both legs here in fact come out bit-for-bit
# equal, flat mesh and terrain-fitted alike, so they are held at zero tolerance.
# ALLOW_DIFF_GRIDS because the restart leg writes its plotfile on the re-made grids.
if(ERF_ENABLE_MPI AND NOT WIN32)
add_test_restart_parity(TerrainHill_RegridOnRestart TerrainHill 7 20
    COMMON_OPTIONS  "amr.max_level=0 erf.terrain_type=None"
    RESTART_OPTIONS "erf.regrid_level_0_on_restart=1 amr.max_grid_size_x=8 amr.max_grid_size_y=8 amr.max_grid_size_z=64"
    ALLOW_DIFF_GRIDS REQUIRE_LEVEL0_REMAKE FCOMPARE_RTOL "0.0" FCOMPARE_ATOL "0.0")
add_test_restart_parity(TerrainHill_RegridOnRestart_Fitted TerrainHill 7 20
    COMMON_OPTIONS  "amr.max_level=0"
    RESTART_OPTIONS "erf.regrid_level_0_on_restart=1 amr.max_grid_size_x=8 amr.max_grid_size_y=8 amr.max_grid_size_z=64"
    ALLOW_DIFF_GRIDS REQUIRE_LEVEL0_REMAKE FCOMPARE_RTOL "0.0" FCOMPARE_ATOL "0.0")
endif()

# The automatic branch, end to end and with nothing set. ERF::restart re-makes the level-0
# grids on its own when the checkpoint has fewer level-0 boxes than there are ranks, and a
# fresh start gives level 0 exactly NProcs() boxes, so writing the checkpoint on one rank and
# restarting on two is enough to take it. That is what a user gets by continuing a run on more
# ranks than it was written with -- no flag, no warning -- and before #4240 it silently
# corrupted the state rather than failing. Nothing else covers a rank count that changes
# across the checkpoint: every other restart test runs all three legs at the same width.
#
# The straight leg runs at the restart's rank count, so the comparison isolates the restart
# rather than the decomposition. This deck carries none of the state the level-0 remake drops,
# so it is required to be bit-for-bit exact; a deck that does carry such state is refused by
# the guard instead, which the RefusesOnMoreRanks case below asserts.
if(ERF_ENABLE_MPI AND NOT WIN32)
add_test_restart_parity(TerrainHill_RestartOnMoreRanks TerrainHill 7 20
    COMMON_OPTIONS "amr.max_level=0 erf.terrain_type=None"
    CHK_NRANKS 1 RESTART_NRANKS 2
    ALLOW_DIFF_GRIDS REQUIRE_LEVEL0_REMAKE FCOMPARE_RTOL "0.0" FCOMPARE_ATOL "0.0")
endif()

# The state that only the checkpoint carries. These used to be refused, because the level-0
# remake rebuilt the level through init_stuff and dropped them; now the checkpoint is read
# straight onto the new grids, so each is simply required to survive. Three categories the
# TerrainHill deck can switch on: the microphysics accumulators, the velocity time averages
# and the interval means. The first runs on the automatic branch -- checkpoint on one rank,
# restart on two, nothing set -- because that is how a user reaches this without asking.
if(ERF_ENABLE_MPI AND NOT WIN32)
# PLT2DFILE because rain_accum is a surface field: the 3D plotfile carries none, so without
# it this would compare density and velocity -- which the remake never lost -- and pass while
# the accumulator silently restarted from zero.
add_test_restart_parity(TerrainHill_RegridOnRestart_Moisture TerrainHill 7 20
    COMMON_OPTIONS "amr.max_level=0 erf.terrain_type=None erf.moisture_model=Kessler erf.plot2d_vars_1=precip_rain_accum"
    CHK_NRANKS 1 RESTART_NRANKS 2
    PLT2DFILE "plt2d_1_00020"
    ALLOW_DIFF_GRIDS REQUIRE_LEVEL0_REMAKE FCOMPARE_RTOL "0.0" FCOMPARE_ATOL "0.0")
add_test_restart_parity(TerrainHill_RegridOnRestart_TimeAvg TerrainHill 7 20
    COMMON_OPTIONS  "amr.max_level=0 erf.terrain_type=None erf.time_avg_vel=true erf.plot_vars_1=density x_velocity y_velocity theta u_t_avg v_t_avg umag_t_avg"
    RESTART_OPTIONS "erf.regrid_level_0_on_restart=1 amr.max_grid_size_x=8 amr.max_grid_size_y=8 amr.max_grid_size_z=64"
    ALLOW_DIFF_GRIDS REQUIRE_LEVEL0_REMAKE FCOMPARE_RTOL "0.0" FCOMPARE_ATOL "0.0")
# Not zero tolerance, and not because of the regrid: the interval means are already not
# restart-exact without one. A plain same-rank restart of this deck reproduces theta_mean to
# 1 ulp (1.9e-16 relative) and w_mean to 6.5e-23 absolute, on a w_mean that is itself ~1e-21,
# and the same two numbers appear whether or not level 0 is re-made. Two decompositions run
# straight through, with no restart at all, agree exactly, so it is the restart and not the
# decomposition. That is issue 4243; the bound here is set just above what it costs so this
# test still fails if the means are actually lost, which would be O(1). Put it back to zero
# when 4243 is fixed.
add_test_restart_parity(TerrainHill_RegridOnRestart_IntervalMeans TerrainHill 7 20
    COMMON_OPTIONS  "amr.max_level=0 erf.terrain_type=None erf.compute_mean_vars=true erf.plot_vars_1=density x_velocity theta u_mean v_mean w_mean theta_mean"
    RESTART_OPTIONS "erf.regrid_level_0_on_restart=1 amr.max_grid_size_x=8 amr.max_grid_size_y=8 amr.max_grid_size_z=64"
    ALLOW_DIFF_GRIDS REQUIRE_LEVEL0_REMAKE FCOMPARE_RTOL "1.0e-12" FCOMPARE_ATOL "1.0e-20")
endif()

# SLM across a restart. Nothing covered this: the two SLM decks in the tree are gated behind
# ERF_TEST_ENABLE_EXTRA_LSM_TESTS and need input files this repository does not carry, so an
# SLM restart was never exercised at all -- which is how the failure reported on issue 4225
# (SLM stopping on t_canop > tfriz after a restart) got in.
#
# PLT2DFILE is plt_lsm_2D, not an erf.plot2d file: SLM keeps prognostic state ERF does not
# expose through lsm_data -- t_canop, t_skin, t_ground_skin, t_cas, q_cas, mw, mws,
# wet_canop -- and plt_lsm_2D is the only output that carries it. erf.plot_lsm in the deck
# writes it alongside plotfile 1, so the runner's existing plot_int_1 cadence drives it. The
# 3D comparison alone would see the atmosphere and pass while the canopy state silently
# restarted from its initialization value.
if(ERF_ENABLE_MPI AND NOT WIN32)
add_test_restart_parity(SLM_Restart SLM_Restart 7 20
    PLT2DFILE "plt_lsm_2D_00020"
    FCOMPARE_RTOL "0.0" FCOMPARE_ATOL "0.0")
endif()

# The same deck across a level-0 regrid, which ReadCheckpointFile still refuses: SLM does its
# own checkpoint I/O and reads onto the grids the run uses, so it cannot yet take grids that
# moved. Checkpoint on one rank and restart on two to reach the automatic branch without
# setting anything. When SLM's reader is taught to redistribute, this becomes an
# add_test_restart_parity case like the one above and the category leaves the guard.
if(ERF_ENABLE_MPI AND NOT WIN32)
add_test_restart_abort(SLM_RegridOnRestart_Refused SLM_Restart 7
    "a land-surface model"
    CHK_NRANKS 1 RESTART_NRANKS 2)
endif()

# The same for terrain carried by an embedded boundary: the mesh is flat there
# too, so the ground under the station is the surface the EB was built from and
# not the bottom of the mesh.  The hill is 50 m up at the station, so 40 m above
# the terrain is 90 m above z = 0.  Measuring from the mesh instead puts the
# station at 40 m, inside the hill, where the velocity is held at zero and the
# checker reports a constant series.
#add_test_station_series(StationSampling_EBTerrain HillEB single
#    RUNTIME_OPTIONS "erf.station_names=mast mastabs erf.mastabs.field=x_velocity theta erf.mastabs.x=400.0 erf.mastabs.y=10.0 erf.mastabs.height_abs=90.0"
#    CHECKS "equal a=@RUN@/Output_Stations/mast.dat:2 b=@RUN@/Output_Stations/mastabs.dat:2 tol=0.001")

#=============================================================================
# Observation nudging
#=============================================================================

# The obs-nudging tests run through add_test_station_series (above), in its
# analytic, approach, single and abort modes.

# A uniform flow relaxing to a target linear in time: u, v and theta at three
# heights of a probe far from the station must follow the closed form (see the
# deck).  The error is a few 1e-6 at dt / tau = 1/20; 1e-4 still fails a rate,
# time interpolation or target that is wrong by a percent.
add_test_station_series(ObsNudging_Uniform ObsNudging_Uniform analytic
    STATION "probe"
    LABELS "regression;obs-nudging"
    CHECKS "col=2 phi0=4.0 a=6.0 b=0.01 tau=20.0 tol=1.0e-4|col=7 phi0=1.0 a=-1.0 b=0.0 tau=20.0 tol=1.0e-4|col=13 phi0=300.0 a=301.0 b=0.005 tau=20.0 tol=1.0e-4")

# The same flow with the band mean +- sigma: u and theta start below the band
# and v above it, so each relaxes to the near edge, mean - sigma (u: 6 - 0.5,
# theta: 301 - 0.2) or mean + sigma (v: -1 + 0.25), and never enters it.
add_test_station_series(ObsNudging_Uniform_SigmaBand ObsNudging_Uniform analytic
    RUNTIME_OPTIONS "erf.obs_nudging.sigma_factor=1.0"
    STATION "probe"
    LABELS "regression;obs-nudging"
    CHECKS "col=2 phi0=4.0 a=5.5 b=0.01 tau=20.0 tol=1.0e-4|col=7 phi0=1.0 a=-0.75 b=0.0 tau=20.0 tol=1.0e-4|col=13 phi0=300.0 a=300.8 b=0.005 tau=20.0 tol=1.0e-4")

# Over a hill on a terrain-fitted mesh, with a refined level from the ground up:
# after 10 s the nudged run must be closer than the free one to the measurements
# interpolated to that time (mast: u 7.0167, v 0.5083, theta 300.5083 at 40 m;
# lidar gate at 120 m: u 6.4083, w 0.4), by the factors below.  The ratios
# measured when the test was written are 0.53, 0.60, 0.36, 0.50 and 0.71.
add_test_station_series(ObsNudging_Hill ObsNudging_Hill approach
    LABELS "regression;obs-nudging"
    OFF_OPTIONS "erf.nudging_from_observations=false"
    CHECKS "series=mast col=2 target=7.0167 factor=0.7|series=mast col=3 target=0.5083 factor=0.75|series=mast col=4 target=300.5083 factor=0.6|series=gate col=2 target=6.4083 factor=0.7|series=gate col=3 target=0.4 factor=0.85")

# The parity drivers run each leg in a subdirectory, so the deck's data files
# are named by absolute path.
function(obs_nudging_hill_files TEST_NAME OUT_VAR)
    set(_d "${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME}")
    set(${OUT_VAR} "erf.input_sounding_file=${_d}/input_sounding erf.obs_nudging.mast.file=${_d}/mast.txt erf.obs_nudging.lidar.file=${_d}/lidar.txt" PARENT_SCOPE)
endfunction()

# The nudging reads a face's neighbours (the density, the mesh and the terrain
# under it), the terrain is gathered from the boxes that touch the ground, and
# the two copies of a face on a periodic boundary must see the same distance to
# every station, so the answer is checked across a change of decomposition on
# both levels, with the station series compared as well as the plotfile.
obs_nudging_hill_files(ObsNudging_Hill_BoxParity _obs_files)
add_test_box_parity(ObsNudging_Hill_BoxParity ObsNudging_Hill "plt00020"
    COMMON_OPTIONS "${_obs_files}"
    REFERENCE_OPTIONS "amr.max_grid_size=1024"
    SPLIT_OPTIONS "amr.max_grid_size_x=8 amr.max_grid_size_y=8 amr.max_grid_size_z=64"
    DATALOG "Output_Stations/gate.dat"
    DATALOG_SIGDIGITS 10)

# The same hill as an immersed boundary in a flat mesh: the station heights are
# measured from the immersed terrain surface, 53 m up at the stations, so the
# nudging reaches the station series (placed at the absolute heights 93.08 and
# 173.08 m) only if the terrain under the stations is found.  The ratios
# measured when the test was written are 0.72, 0.74, 0.52, 0.46 and 0.79.
add_test_station_series(ObsNudging_HillIF ObsNudging_HillIF approach
    LABELS "regression;obs-nudging"
    OFF_OPTIONS "erf.nudging_from_observations=false"
    CHECKS "series=mast col=2 target=7.0167 factor=0.85|series=mast col=3 target=0.5083 factor=0.85|series=mast col=4 target=300.5083 factor=0.7|series=gate col=2 target=6.4083 factor=0.65|series=gate col=3 target=0.4 factor=0.9")

function(obs_nudging_hillif_files TEST_NAME OUT_VAR)
    set(_d "${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME}")
    set(${OUT_VAR} "erf.input_sounding_file=${_d}/input_sounding erf.obs_nudging.mast.file=${_d}/mast.txt erf.obs_nudging.lidar.file=${_d}/lidar.txt" PARENT_SCOPE)
endfunction()
obs_nudging_hillif_files(ObsNudging_HillIF_BoxParity _obs_files)
add_test_box_parity(ObsNudging_HillIF_BoxParity ObsNudging_HillIF "plt00020"
    COMMON_OPTIONS "${_obs_files}"
    REFERENCE_OPTIONS "amr.max_grid_size=1024"
    SPLIT_OPTIONS "amr.max_grid_size_x=8 amr.max_grid_size_y=8 amr.max_grid_size_z=64"
    DATALOG "Output_Stations/gate.dat"
    DATALOG_SIGDIGITS 10)

# Every input the nudging checks must be refused at start-up, before the first
# step, with a message that names it.  Each test changes one input of a deck
# that otherwise runs.
foreach(_case IN ITEMS
        "NoTau|StationSampling|erf.nudging_from_observations=true|needs erf.obs_nudging.tau"
        "NoStations|StationSampling|erf.nudging_from_observations=true erf.obs_nudging.tau=10|needs erf.obs_nudging.stations"
        "EBTerrain|StationSampling|erf.nudging_from_observations=true erf.terrain_type=EB|does not support erf.terrain_type = EB"
        "NegativeTau|ObsNudging_Uniform|erf.obs_nudging.tau=-1|erf.obs_nudging.tau must be positive"
        "ZeroRadius|ObsNudging_Uniform|erf.obs_nudging.horizontal_radius=0|horizontal_radius must be positive"
        "NegativeSigma|ObsNudging_Uniform|erf.obs_nudging.sigma_factor=-1|sigma_factor must not be negative"
        "NothingNudged|ObsNudging_Uniform|erf.obs_nudging.nudge_wind=false erf.obs_nudging.nudge_w=false erf.obs_nudging.nudge_theta=false|nothing would be nudged"
        "NoUsableStation|ObsNudging_Uniform|erf.obs_nudging.nudge_wind=false erf.obs_nudging.nudge_theta=false|no station measures any of the quantities"
        "EpochNoDate|ObsNudging_Uniform|erf.obs_nudging.time_type=epoch|needs the run to know its start date"
        "BadTimeType|ObsNudging_Uniform|erf.obs_nudging.time_type=utc|must be elapsed or epoch"
        "BadWindFrame|ObsNudging_Uniform|erf.obs_nudging.mast.wind_frame=north|wind_frame must be earth or grid"
        "BadHeightRef|ObsNudging_Uniform|erf.obs_nudging.mast.height_ref=asl|height_ref must be agl or msl"
        "MissingFile|ObsNudging_Uniform|erf.obs_nudging.mast.file=no_such_station.txt|cannot open 'no_such_station.txt'"
        "BadColumn|ObsNudging_Uniform|erf.obs_nudging.mast.file=bad_columns_station.txt|unknown column 'pressure'"
        "HeightMismatch|ObsNudging_Uniform|erf.obs_nudging.mast.file=mismatched_heights_station.txt|does not match the heights of the first time"
        "XYAndLatLon|ObsNudging_Uniform|erf.obs_nudging.mast.lat=40.0 erf.obs_nudging.mast.long=-105.0|give either lat/long or x/y, not both"
        "LatLonNoArrays|ObsNudging_Uniform|erf.obs_nudging.stations=geo erf.obs_nudging.geo.file=uniform_station.txt erf.obs_nudging.geo.lat=40.0 erf.obs_nudging.geo.long=-105.0|lat/long needs a run with latitude/longitude arrays"
        "OutsideDomain|ObsNudging_Uniform|erf.obs_nudging.mast.x=5000.0|is outside the problem domain"
        "LevelOffGround|ObsNudging_Hill|erf.box1.in_box_lo=400.0 200.0 100.0|do not reach the bottom of the domain")
    string(REPLACE "|" ";" _fields "${_case}")
    list(GET _fields 0 _name)
    list(GET _fields 1 _deck)
    list(GET _fields 2 _options)
    list(GET _fields 3 _message)
    add_test_station_series(ObsNudging_Abort_${_name} ${_deck} abort
        LABELS "regression;obs-nudging"
        RUNTIME_OPTIONS "${_options}"
        EXPECTED_MESSAGE "${_message}")
endforeach()

# The targets are interpolated in time from the run time and the terrain is
# rebuilt from the grids, so a restart must continue the run exactly.
obs_nudging_hill_files(ObsNudging_Hill_Restart _obs_files)
add_test_restart_parity(ObsNudging_Hill_Restart ObsNudging_Hill 10 20
    COMMON_OPTIONS "${_obs_files}"
    DATALOG "Output_Stations/mast.dat"
    DATALOG_SIGDIGITS 10)

#=============================================================================
# Smagorinsky2D on terrain: WRF smag2d_km limits and the diffusive time-step check
#=============================================================================
# Steep ridge (alpha about 19.5 in the first cells), 3 km grid, Smagorinsky2D + MRF, no numerical
# diffusion.  Measured with erf_exec (Release, 2 ranks, 6 simulated hours per run): the unlimited
# closure survives dt = 45 s and fails from 46 s; with erf.smag2d_slope_limiter it survives to
# 70 s and fails at 75 s, the same edge as with no LES at all.  The test runs the limiter at
# 65 s, and the control (no limiter) at 65 s must fail (it does at step 2).
add_test_smag2d_ridge(Smag2D_Ridge_SteepLimiter steep "plt00060"
    OPTIONS "erf.smag2d_slope_limiter=true erf.fixed_dt=65 erf.fixed_mri_dt_ratio=40 max_step=60 erf.plot_int_1=60"
    CONTROL_OPTIONS "erf.fixed_dt=65 erf.fixed_mri_dt_ratio=40 max_step=60 erf.plot_int_1=60"
    ALPHA_MIN 15
    WMAX 15)
# The opt-in diffusive check warns on this deck and leaves the answer unchanged, with an
# adaptive dt so that a check that moved dt would show.
add_test_smag2d_ridge(Smag2D_Ridge_DiffusiveCheck check "plt00010"
    OPTIONS "erf.fixed_dt=-1 max_step=10 erf.plot_int_1=10"
    DIFFUSIVE_CFL 0.3)
# erf.diffusive_dt_limit sets an adaptive dt and is never exceeded; the first diffusive dt
# (3.6881 s in the reference run, Release, 1 and 2 ranks) must match to 2 %.
add_test_smag2d_ridge(Smag2D_Ridge_DiffusiveLimit limit "plt00010"
    OPTIONS "erf.fixed_dt=-1 max_step=10 erf.plot_int_1=10"
    DIFFUSIVE_CFL 0.2
    REF_DIFFUSIVE_DT_LO 3.6143
    REF_DIFFUSIVE_DT_HI 3.7619)

# Each new input is checked at start-up and aborts naming itself.  add_test_abort runs through
# `sh -c ... | tee`, so these sit under the same guard as its other use.
if(ERF_ENABLE_MPI AND NOT WIN32)
set(_smag2d_dir ${CMAKE_CURRENT_SOURCE_DIR}/test_files/Smag2D_Ridge)
add_test_abort(Smag2D_Abort_SlopeWithout2D ${_smag2d_dir} Smag2D_Ridge.i
               "apply only to erf.les_type = Smagorinsky2D"
               "erf.les_type=Smagorinsky erf.Cs=0.1 erf.pbl_type=None erf.smag2d_slope_limiter=true")
add_test_abort(Smag2D_Abort_CapWithout2D ${_smag2d_dir} Smag2D_Ridge.i
               "apply only to erf.les_type = Smagorinsky2D"
               "erf.les_type=None erf.smag2d_kh_cap=10")
add_test_abort(Smag2D_Abort_CapNegative ${_smag2d_dir} Smag2D_Ridge.i
               "erf.smag2d_kh_cap must be >= 0"
               "erf.smag2d_kh_cap=-1")
add_test_abort(Smag2D_Abort_SlopeWithoutTerrain ${_smag2d_dir} Smag2D_Ridge.i
               "erf.smag2d_slope_limiter requires erf.terrain_type"
               "erf.terrain_type=None erf.grid_stretching_ratio=0 prob.custom_terrain_type=None erf.smag2d_slope_limiter=true")
add_test_abort(Smag2D_Abort_CapWithEB ${CMAKE_CURRENT_SOURCE_DIR}/test_files/HillEB HillEB.i
               "erf.smag2d_kh_cap is not supported with erf.terrain_type = EB"
               "erf.les_type=Smagorinsky2D erf.Cs=0.1 erf.smag2d_kh_cap=10")
add_test_abort(DiffusiveDt_Abort_CflRange ${_smag2d_dir} Smag2D_Ridge.i
               "erf.diffusive_cfl must be in"
               "erf.diffusive_cfl=1.5")
add_test_abort(DiffusiveDt_Abort_CflUnused ${_smag2d_dir} Smag2D_Ridge.i
               "erf.diffusive_cfl is used only with"
               "erf.diffusive_dt_check=false erf.diffusive_cfl=0.3")
add_test_abort(DiffusiveDt_Abort_LimitFixedDt ${_smag2d_dir} Smag2D_Ridge.i
               "erf.diffusive_dt_limit = true cannot change erf.fixed_dt"
               "erf.diffusive_dt_limit=true")
add_test_abort(DiffusiveDt_Abort_LimitNoClosure ${_smag2d_dir} Smag2D_Ridge.i
               "erf.diffusive_dt_limit = true needs an eddy-diffusivity closure"
               "erf.diffusive_dt_limit=true erf.fixed_dt=-1 erf.les_type=None erf.pbl_type=None")
endif()

#=============================================================================
# Performance tests
#=============================================================================
