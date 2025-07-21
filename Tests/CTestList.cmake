# Have CMake discover the number of cores on the node
include(ProcessorCount)
ProcessorCount(PROCESSES)

#=============================================================================
# Functions for adding tests / Categories of tests
#=============================================================================
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

    if(WIN32)
        set(TEST_EXE "${CMAKE_BINARY_DIR}/Exec/${TEST_DIR}/*/${TEST_EXE}.exe")
    else()
        set(TEST_EXE "${CMAKE_BINARY_DIR}/Exec/${TEST_DIR}/${TEST_EXE}")
    endif()

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

    if(WIN32)
        set(TEST_EXE "${CMAKE_BINARY_DIR}/Exec/${TEST_DIR}/*/${TEST_EXE}.exe")
    else()
        set(TEST_EXE "${CMAKE_BINARY_DIR}/Exec/${TEST_DIR}/${TEST_EXE}")
    endif()
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

    if(WIN32)
        set(TEST_EXE "${CMAKE_BINARY_DIR}/Exec/${TEST_DIR}/*/${TEST_EXE}.exe")
    else()
        set(TEST_EXE "${CMAKE_BINARY_DIR}/Exec/${TEST_DIR}/${TEST_EXE}")
    endif()
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

#=============================================================================
# Regression tests
#=============================================================================

add_test_r(DensityCurrent                    "DryRegTests/DensityCurrent" "erf_density_current" "plt00010")
add_test_r(DensityCurrent_anelastic          "DryRegTests/DensityCurrent" "erf_density_current" "plt00010")
add_test_r(DensityCurrent_detJ2              "DryRegTests/DensityCurrent" "erf_density_current" "plt00010")
add_test_r(DensityCurrent_detJ2_nosub        "DryRegTests/DensityCurrent" "erf_density_current" "plt00020")
add_test_r(DensityCurrent_detJ2_MT           "DryRegTests/DensityCurrent" "erf_density_current" "plt00010")
add_test_r(EkmanSpiral                       "DryRegTests/EkmanSpiral" "erf_ekman_spiral" "plt00010")
#add_test_r(FlowInABox                        "DevTests/FlowInABox"     "erf_flow_in_a_box" "plt00010")
add_test_r(IsentropicVortexStationary        "DryRegTests/IsentropicVortex" "erf_isentropic_vortex" "plt00010")
add_test_r(IsentropicVortexAdvecting         "DryRegTests/IsentropicVortex" "erf_isentropic_vortex" "plt00010")
add_test_r(IVA_NumDiff                       "DryRegTests/IsentropicVortex" "erf_isentropic_vortex" "plt00010")
add_test_r(MovingTerrain_nosub               "DevTests/MovingTerrain" "erf_moving_terrain"   "plt00020")
add_test_r(MovingTerrain_sub                 "DevTests/MovingTerrain" "erf_moving_terrain"   "plt00010")
add_test_r(RayleighDamping                   "DryRegTests/ScalarAdvDiff" "erf_scalar_advdiff" "plt00100")
add_test_r(ScalarAdvectionUniformU           "DryRegTests/ScalarAdvDiff" "erf_scalar_advdiff" "plt00020")
add_test_r(ScalarAdvectionShearedU           "DryRegTests/ScalarAdvDiff" "erf_scalar_advdiff" "plt00080")
add_test_r(ScalarAdvDiff_order2              "DryRegTests/ScalarAdvDiff" "erf_scalar_advdiff" "plt00020")
add_test_r(ScalarAdvDiff_order3              "DryRegTests/ScalarAdvDiff" "erf_scalar_advdiff" "plt00020")
add_test_r(ScalarAdvDiff_order4              "DryRegTests/ScalarAdvDiff" "erf_scalar_advdiff" "plt00020")
add_test_r(ScalarAdvDiff_order5              "DryRegTests/ScalarAdvDiff" "erf_scalar_advdiff" "plt00020")
add_test_r(ScalarAdvDiff_order6              "DryRegTests/ScalarAdvDiff" "erf_scalar_advdiff" "plt00020")
add_test_r(ScalarAdvDiff_weno3               "DryRegTests/ScalarAdvDiff" "erf_scalar_advdiff" "plt00020")
add_test_d(ScalarAdvDiff_weno3z              "DryRegTests/ScalarAdvDiff" "erf_scalar_advdiff" "plt00020")
add_test_r(ScalarAdvDiff_weno5               "DryRegTests/ScalarAdvDiff" "erf_scalar_advdiff" "plt00020")
add_test_d(ScalarAdvDiff_weno5z              "DryRegTests/ScalarAdvDiff" "erf_scalar_advdiff" "plt00020")
add_test_r(ScalarAdvDiff_wenomzq3            "DryRegTests/ScalarAdvDiff" "erf_scalar_advdiff" "plt00020")
add_test_r(ScalarDiffusionGaussian           "DryRegTests/ScalarAdvDiff" "erf_scalar_advdiff" "plt00020")
add_test_r(ScalarDiffusionSine               "DryRegTests/ScalarAdvDiff" "erf_scalar_advdiff" "plt00020")
add_test_r(TaylorGreenAdvecting              "DryRegTests/TaylorGreenVortex" "erf_taylor_green" "plt00010")
add_test_r(TaylorGreenAdvectingDiffusing     "DryRegTests/TaylorGreenVortex" "erf_taylor_green" "plt00010")
add_test_r(MSF_NoSub_IsentropicVortexAdv     "DryRegTests/IsentropicVortex" "erf_isentropic_vortex" "plt00010")
add_test_r(MSF_Sub_IsentropicVortexAdv       "DryRegTests/IsentropicVortex" "erf_isentropic_vortex" "plt00010")
add_test_r(ABL_MOST                          "ABL" "erf_abl" "plt00010")
add_test_r(ABL_MYNN_PBL                      "ABL" "erf_abl" "plt00100" INPUT_SOUNDING "input_sounding_GABLS1")
add_test_r(ABL_InflowFile                    "ABL" "erf_abl" "plt00010")
add_test_r(MoistBubble                       "MoistRegTests/Bubble" "erf_bubble" "plt00010")
add_test_r(SquallLine_2D                     "MoistRegTests/SquallLine_2D" "erf_squallline" "plt00010")
add_test_r(SuperCell_3D                      "MoistRegTests/SuperCell_3D" "erf_supercell"   "plt00010")
if(ERF_ENABLE_PARTICLES)
  add_test_r(ParticleAdvect                  "DryRegTests/ParticleAdvection" "erf_particles_advect" "plt00010")
  add_test_r(ParticleWoA                     "DryRegTests/ParticlesOverWoA" "erf_particles_over_woa" "plt00010")
endif()
if(ERF_ENABLE_RRGMTP)
  add_test_r(Radiation                       "DevTests/Radiation" "erf_radiation" "plt00010")
endif()

add_test_0(CouetteFlow_x                     "DryRegTests/Couette_Poiseuille" "erf_couette_poiseuille" "plt00050")
add_test_0(CouetteFlow_y                     "DryRegTests/Couette_Poiseuille" "erf_couette_poiseuille" "plt00050")
add_test_0(PoiseuilleFlow_x                  "DryRegTests/Couette_Poiseuille" "erf_couette_poiseuille" "plt00010")
add_test_0(PoiseuilleFlow_y                  "DryRegTests/Couette_Poiseuille" "erf_couette_poiseuille" "plt00010")
add_test_0(InitSoundingIdeal_stationary      "ABL" "erf_abl" "plt00010")
add_test_0(Deardorff_stationary              "ABL" "erf_abl" "plt00010")

#=============================================================================
# Performance tests
#=============================================================================
