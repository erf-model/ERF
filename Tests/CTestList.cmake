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
    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${TEST_NAME})
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

# Standard regression test
function(add_test_r TEST_NAME TEST_DIR TEST_EXE PLTFILE)
    set(options )
    set(oneValueArgs "INPUT_SOUNDING" "RUNTIME_OPTIONS")
    set(multiValueArgs )
    cmake_parse_arguments(ADD_TEST_R "${options}" "${oneValueArgs}"
        "${multiValueArgs}" ${ARGN})

    setup_test()

    set(RUNTIME_OPTIONS "${ADD_TEST_R_RUNTIME_OPTIONS}")
    if(NOT "${ADD_TEST_R_INPUT_SOUNDING}" STREQUAL "")
      string(APPEND RUNTIME_OPTIONS "erf.input_sounding_file=${CURRENT_TEST_BINARY_DIR}/${ADD_TEST_R_INPUT_SOUNDING}")
    endif()

    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)

    set(FCOMPARE_TOLERANCE "-r ${ERF_TEST_FCOMPARE_RTOL} --abs_tol ${ERF_TEST_FCOMPARE_ATOL}")
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

# Debug regression test with lower tolerance
function(add_test_d TEST_NAME TEST_DIR TEST_EXE PLTFILE)
    setup_test()

    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)
    set(FCOMPARE_TOLERANCE "-r 3.0e-9 --abs_tol 3.0e-9")
    set(FCOMPARE_FLAGS "--abort_if_not_all_found -a ${FCOMPARE_TOLERANCE}")
    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i > ${TEST_NAME}.log && ${MPI_FCOMP_COMMANDS} ${FCOMPARE_EXE} ${FCOMPARE_FLAGS} ${PLOT_GOLD} ${CURRENT_TEST_BINARY_DIR}/${PLTFILE}")

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
    setup_test()

    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)
    set(FCOMPARE_TOLERANCE "-r 1e-14 --abs_tol 1.0e-14")
    set(FCOMPARE_FLAGS "-a ${FCOMPARE_TOLERANCE}")
    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i erf.input_sounding_file=${CURRENT_TEST_BINARY_DIR}/input_sounding > ${TEST_NAME}.log && ${MPI_FCOMP_COMMANDS} ${FCOMPARE_EXE} ${FCOMPARE_FLAGS} ${CURRENT_TEST_BINARY_DIR}/plt00000 ${CURRENT_TEST_BINARY_DIR}/${PLTFILE}")

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

    set(FCOMPARE_TOLERANCE "--rel_tol ${TEST_RTOL} --abs_tol ${TEST_ATOL}")
    set(FCOMPARE_FLAGS "--abort_if_not_all_found --allow_diff_grids ${FCOMPARE_TOLERANCE}")
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
endfunction(add_test_sdm)

#=============================================================================
# Regression tests
#=============================================================================

# These tests will all be built in Exec
add_test_r(DensityCurrent                    ""  "erf_exec" "plt00010")
add_test_r(DensityCurrent_anelastic          ""  "erf_exec" "plt00010")
add_test_r(DensityCurrent_detJ2              ""  "erf_exec" "plt00010")
add_test_r(DensityCurrent_detJ2_nosub        ""  "erf_exec" "plt00020")
add_test_r(DensityCurrent_detJ2_MT           ""  "erf_exec" "plt00010")
add_test_r(EkmanSpiral                       ""  "erf_exec" "plt00010")
add_test_r(IsentropicVortexStationary        ""  "erf_exec" "plt00010")
add_test_r(IsentropicVortexAdvecting         ""  "erf_exec" "plt00010")
add_test_r(IVA_NumDiff                       ""  "erf_exec" "plt00010")
add_test_r(MovingTerrain_nosub               ""  "erf_exec" "plt00020")
add_test_r(MovingTerrain_sub                 ""  "erf_exec" "plt00010")
add_test_r(RayleighDamping                   ""  "erf_exec" "plt00100")
add_test_r(ScalarAdvectionUniformU           ""  "erf_exec" "plt00020")
add_test_r(ScalarAdvectionShearedU           ""  "erf_exec" "plt00080")
add_test_r(ScalarAdvDiff_order2              ""  "erf_exec" "plt00020")
add_test_r(ScalarAdvDiff_order3              ""  "erf_exec" "plt00020")
add_test_r(ScalarAdvDiff_order4              ""  "erf_exec" "plt00020")
add_test_r(ScalarAdvDiff_order5              ""  "erf_exec" "plt00020")
add_test_r(ScalarAdvDiff_order6              ""  "erf_exec" "plt00020")
add_test_r(ScalarAdvDiff_weno3               ""  "erf_exec" "plt00020")
add_test_d(ScalarAdvDiff_weno3z              ""  "erf_exec" "plt00020")
add_test_r(ScalarAdvDiff_weno5               ""  "erf_exec" "plt00020")
add_test_d(ScalarAdvDiff_weno5z              ""  "erf_exec" "plt00020")
add_test_r(ScalarAdvDiff_wenomzq3            ""  "erf_exec" "plt00020")
add_test_r(ScalarDiffusionGaussian           ""  "erf_exec" "plt00020")
add_test_r(ScalarDiffusionSine               ""  "erf_exec" "plt00020")
add_test_r(TaylorGreenAdvecting              ""  "erf_exec" "plt00010")
add_test_r(TaylorGreenAdvectingDiffusing     ""  "erf_exec" "plt00010")
add_test_r(MSF_NoSub_IsentropicVortexAdv     ""  "erf_exec" "plt00010")
add_test_r(MSF_Sub_IsentropicVortexAdv       ""  "erf_exec" "plt00010")
#add_test_r(FlowInABox                       ""  "erf_exec" "plt00010")
add_test_r(ABL_MOST                          ""  "erf_exec" "plt00010")
add_test_r(ABL_MOST_IMP_DIFF                 ""  "erf_exec" "plt00010")
add_test_r(ABL_MYNN_PBL                      ""  "erf_exec" "plt00100" INPUT_SOUNDING "input_sounding_GABLS1")
add_test_r(ABL_InflowFile                    ""  "erf_exec" "plt00010")
add_test_r(MoistBubble                       ""  "erf_exec" "plt00010")
add_test_r(SquallLine_2D                     ""  "erf_exec" "plt00010")
add_test_r(SuperCell_3D                      ""  "erf_exec" "plt00010")
if(ERF_ENABLE_PARTICLES)
  add_test_r(ParticleAdvect                  ""  "erf_exec" "plt00010")
  add_test_r(ParticleWoA                     ""  "erf_exec" "plt00010")
endif( )
if(ERF_ENABLE_RRGMTP)
  add_test_r(Radiation                       ""  "erf_exec" "plt00010")
endif()

add_test_0(CouetteFlow_x                     "" "erf_exec" "plt00050")
add_test_0(CouetteFlow_y                     "" "erf_exec" "plt00050")
add_test_0(PoiseuilleFlow_x                  "" "erf_exec" "plt00010")
add_test_0(PoiseuilleFlow_y                  "" "erf_exec" "plt00010")
add_test_0(InitSoundingIdeal_stationary      "" "erf_exec" "plt00010")
add_test_0(Deardorff_stationary              "" "erf_exec" "plt00010")

if(ERF_ENABLE_PARTICLES)
    # These tests require machine-specific gold files due to platform-dependent initial sampling
    if(ERF_TEST_ENABLE_EXTRA_SDM_TESTS)
        # log-normal distribution for radius
        add_test_sdm(SDM_RICO3D_InitSampling         ""  "erf_exec"   "plt00000" 1e-14 2e-13 INPUT_SOUNDING "input_sounding")
        # mass-exponential distribution for mass
        add_test_sdm(SDM_Bubble2D_Adv_InitSampling   ""  "erf_exec"   "plt00000" 1e-14 1e-14)
        # column case to test condensation
        add_test_sdm(SDM_SineMassFlux                "" "erf_exec" "plt00050" 1e-14 1e-14 INPUT_SOUNDING "input_sounding")
    endif()

    # passive advection of particles
    add_test_sdm(SDM_Bubble2D_Adv                "" "erf_exec"  "plt00050" 1e-12 1e-12)
    # passive advection of particles with injection
    add_test_sdm(SDM_Bubble2D_Adv_wInjection     "" "erf_exec"  "plt00050" 5e-12 5e-12)
    # condensation/evaporation
    add_test_sdm(SDM_Box3D_Cond                  "" "erf_exec"  "plt00010" 2e-12 3e-13)
    # terminal velocity
    add_test_sdm(SDM_Box3D_VTerm                 "" "erf_exec"  "plt00001" 5e-13 1e-14)
    # recycling
    add_test_sdm(SDM_Box3D_Recycling             "" "erf_exec"  "plt00020" 5e-13 1e-14)
    # Congestus case
    add_test_sdm(SDM_Congestus3D                 "" "erf_exec"  "plt00020" 5e-13 5e-13 INPUT_SOUNDING "input_sounding")
    # RICO case
    add_test_sdm(SDM_RICO3D                      "" "erf_exec"  "plt00010" 5e-13 5e-13 INPUT_SOUNDING "input_sounding")
    # multispecies setup with dummy water species
    add_test_sdm(SDM_MultiSpecies_Bubble2D       "" "erf_exec"  "plt00001" 5e-12 1e-12)
endif()

#=============================================================================
# Performance tests
#=============================================================================
