
# Have CMake discover the number of cores on the node
include(ProcessorCount)
ProcessorCount(PROCESSES)

set(FCOMPARE_GOLD_FILES_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ERFGoldFiles)

#=============================================================================
# Functions for adding tests / Categories of tests
#=============================================================================
macro(setup_test)
    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${TEST_NAME})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    set(PLOT_GOLD ${FCOMPARE_GOLD_FILES_DIRECTORY}/${TEST_NAME})

    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")

    if(ERF_ENABLE_MPI)
        set(NP 4)
        set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NP} ${MPIEXEC_PREFLAGS}")
    else()
        set(NP 1)
        unset(MPI_COMMANDS)
    endif()

    # Set some default runtime options for all tests in this category
    # set(RUNTIME_OPTIONS "time.max_step=10 amr.plot_file=plt time.plot_interval=10 amrex.throw_exception=1 amrex.signal_handling=0")
    # set(RUNTIME_OPTIONS "max_step=10 amr.plot_file=plt amr.checkpoint_files_output=0 amr.plot_files_output=1 amrex.signal_handling=0")    

endmacro(setup_test)

# Standard regression test
function(add_test_r TEST_NAME TEST_EXE PLTFILE)
    setup_test()

    set(TEST_EXE ${CMAKE_BINARY_DIR}/Exec/${TEST_EXE})
    set(FCOMPARE_TOLERANCE "-r 1e-12 --abs_tol 1.0e-12")
    set(FCOMPARE_FLAGS "-a ${FCOMPARE_TOLERANCE}")
    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUNTIME_OPTIONS} > ${TEST_NAME}.log && ${FCOMPARE_EXE} ${FCOMPARE_FLAGS} ${PLOT_GOLD} ${CURRENT_TEST_BINARY_DIR}/${PLTFILE}")

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

# Stationary test -- compare with time 0
function(add_test_0 TEST_NAME TEST_EXE PLTFILE)
    setup_test()

    set(TEST_EXE ${CMAKE_BINARY_DIR}/Exec/${TEST_EXE})
    set(FCOMPARE_TOLERANCE "-r 1e-14 --abs_tol 1.0e-14")
    set(FCOMPARE_FLAGS "-a ${FCOMPARE_TOLERANCE}")
    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i erf.input_sounding_file=${CURRENT_TEST_BINARY_DIR}/input_sounding ${RUNTIME_OPTIONS} > ${TEST_NAME}.log && ${FCOMPARE_EXE} ${FCOMPARE_FLAGS} ${CURRENT_TEST_BINARY_DIR}/plt00000 ${CURRENT_TEST_BINARY_DIR}/${PLTFILE}")

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

# Standard unit test
function(add_test_u TEST_NAME)
    setup_test()
    add_test(${TEST_NAME} sh -c "${MPI_COMMANDS} ${CMAKE_BINARY_DIR}/${amr_wind_unit_test_exe_name}")
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 500
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "unit"
    )
endfunction(add_test_u)

#=============================================================================
# Unit tests
#=============================================================================
# add_test_u(unit_tests)

#=============================================================================
# Regression tests
#=============================================================================
#add_test_r(Bubble_DensityCurrent             "Bubble/bubble" "plt00010")
add_test_r(CouetteFlow                       "RegTests/CouetteFlow/erf_couette_flow" "plt00050")
add_test_r(DensityCurrent                    "RegTests/DensityCurrent/density_current" "plt00010")
add_test_r(DensityCurrent_detJ2              "RegTests/DensityCurrent/density_current" "plt00010")
add_test_r(DensityCurrent_detJ2_nosub        "RegTests/DensityCurrent/density_current" "plt00020")
add_test_r(DensityCurrent_detJ2_MT           "RegTests/DensityCurrent/density_current" "plt00010")
add_test_r(EkmanSpiral                       "RegTests/EkmanSpiral_custom/ekman_spiral_custom" "plt00010")
add_test_r(IsentropicVortexStationary        "RegTests/IsentropicVortex/erf_isentropic_vortex" "plt00010")
add_test_r(IsentropicVortexAdvecting         "RegTests/IsentropicVortex/erf_isentropic_vortex" "plt00010")
add_test_r(MovingTerrain_nosub               "DevTests/MovingTerrain/moving_terrain"   "plt00020")
add_test_r(MovingTerrain_sub                 "DevTests/MovingTerrain/moving_terrain"   "plt00010")
add_test_r(PoiseuilleFlow                    "RegTests/PoiseuilleFlow/erf_poiseuille_flow" "plt00010")
add_test_r(RayleighDamping                   "RegTests/ScalarAdvDiff/erf_scalar_advdiff" "plt00100")
add_test_r(ScalarAdvectionUniformU           "RegTests/ScalarAdvDiff/erf_scalar_advdiff" "plt00020")
add_test_r(ScalarAdvectionShearedU           "RegTests/ScalarAdvDiff/erf_scalar_advdiff" "plt00080")
add_test_r(ScalarAdvDiff_order2              "RegTests/ScalarAdvDiff/erf_scalar_advdiff" "plt00020")
add_test_r(ScalarAdvDiff_order3              "RegTests/ScalarAdvDiff/erf_scalar_advdiff" "plt00020")
add_test_r(ScalarAdvDiff_order4              "RegTests/ScalarAdvDiff/erf_scalar_advdiff" "plt00020")
add_test_r(ScalarAdvDiff_order5              "RegTests/ScalarAdvDiff/erf_scalar_advdiff" "plt00020")
add_test_r(ScalarAdvDiff_order6              "RegTests/ScalarAdvDiff/erf_scalar_advdiff" "plt00020")
add_test_r(ScalarAdvDiff_weno3               "RegTests/ScalarAdvDiff/erf_scalar_advdiff" "plt00020")
add_test_r(ScalarAdvDiff_weno3z              "RegTests/ScalarAdvDiff/erf_scalar_advdiff" "plt00020")
add_test_r(ScalarAdvDiff_weno5               "RegTests/ScalarAdvDiff/erf_scalar_advdiff" "plt00020")
add_test_r(ScalarAdvDiff_weno5z              "RegTests/ScalarAdvDiff/erf_scalar_advdiff" "plt00020")
add_test_r(ScalarAdvDiff_wenomzq3            "RegTests/ScalarAdvDiff/erf_scalar_advdiff" "plt00020")
add_test_r(ScalarDiffusionGaussian           "RegTests/ScalarAdvDiff/erf_scalar_advdiff" "plt00020")
add_test_r(ScalarDiffusionSine               "RegTests/ScalarAdvDiff/erf_scalar_advdiff" "plt00020")
add_test_r(TaylorGreenAdvecting              "RegTests/TaylorGreenVortex/taylor_green" "plt00010")
add_test_r(TaylorGreenAdvectingDiffusing     "RegTests/TaylorGreenVortex/taylor_green" "plt00010")
add_test_r(MSF_NoSub_IsentropicVortexAdv     "RegTests/IsentropicVortex/erf_isentropic_vortex" "plt00010")
add_test_r(MSF_Sub_IsentropicVortexAdv       "RegTests/IsentropicVortex/erf_isentropic_vortex" "plt00010")

add_test_0(Deardorff_stationary              "ABL/erf_abl" "plt00010")

#=============================================================================
# Performance tests
#=============================================================================

