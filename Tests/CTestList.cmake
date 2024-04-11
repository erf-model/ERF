
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
        set(NP ${ERF_TEST_NRANKS})
        set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NP} ${MPIEXEC_PREFLAGS}")
        set(MPI_FCOMP_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}")
    else()
        set(NP 1)
        unset(MPI_COMMANDS)
        unset(MPI_FCOMP_COMMANDS)
    endif()

    # Set some default runtime options for all tests in this category
    # set(RUNTIME_OPTIONS "time.max_step=10 amr.plot_file=plt time.plot_interval=10 amrex.throw_exception=1 amrex.signal_handling=0")
    # set(RUNTIME_OPTIONS "max_step=10 amr.plot_file=plt amr.checkpoint_files_output=0 amr.plot_files_output=1 amrex.signal_handling=0")    

endmacro(setup_test)

# Standard regression test
function(add_test_r TEST_NAME TEST_EXE PLTFILE)
    setup_test()

    set(TEST_EXE ${CMAKE_BINARY_DIR}/Exec/${TEST_EXE})
    set(FCOMPARE_TOLERANCE "-r 2e-10 --abs_tol 2.0e-10")
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

#=============================================================================
# Regression tests
#=============================================================================
if(WIN32)
#add_test_r(Bubble_DensityCurrent             "Bubble/bubble.exe" "plt00010")
add_test_r(CouetteFlow                       "RegTests/Couette_Poiseuille/*/erf_couette_poiseuille.exe" "plt00050")
add_test_r(DensityCurrent                    "RegTests/DensityCurrent/*/erf_density_current.exe" "plt00010")
add_test_r(DensityCurrent_detJ2              "RegTests/DensityCurrent/*/erf_density_current.exe" "plt00010")
add_test_r(DensityCurrent_detJ2_nosub        "RegTests/DensityCurrent/*/erf_density_current.exe" "plt00020")
add_test_r(DensityCurrent_detJ2_MT           "RegTests/DensityCurrent/*/erf_density_current.exe" "plt00010")
add_test_r(EkmanSpiral                       "RegTests/EkmanSpiral/*/erf_ekman_spiral.exe" "plt00010")
add_test_r(IsentropicVortexStationary        "RegTests/IsentropicVortex/*/erf_isentropic_vortex.exe" "plt00010")
add_test_r(IsentropicVortexAdvecting         "RegTests/IsentropicVortex/*/erf_isentropic_vortex.exe" "plt00010")
add_test_r(IVA_NumDiff                       "RegTests/IsentropicVortex/*/erf_isentropic_vortex.exe" "plt00010")
add_test_r(MovingTerrain_nosub               "DevTests/MovingTerrain/*/erf_moving_terrain.exe"   "plt00020")
add_test_r(MovingTerrain_sub                 "DevTests/MovingTerrain/*/erf_moving_terrain.exe"   "plt00010")
add_test_r(PoiseuilleFlow                    "RegTests/Couette_Poiseuille/*/erf_couette_poiseuille.exe" "plt00010")
add_test_r(RayleighDamping                   "RegTests/ScalarAdvDiff/*/erf_scalar_advdiff.exe" "plt00100")
add_test_r(ScalarAdvectionUniformU           "RegTests/ScalarAdvDiff/*/erf_scalar_advdiff.exe" "plt00020")
add_test_r(ScalarAdvectionShearedU           "RegTests/ScalarAdvDiff/*/erf_scalar_advdiff.exe" "plt00080")
add_test_r(ScalarAdvDiff_order2              "RegTests/ScalarAdvDiff/*/erf_scalar_advdiff.exe" "plt00020")
add_test_r(ScalarAdvDiff_order3              "RegTests/ScalarAdvDiff/*/erf_scalar_advdiff.exe" "plt00020")
add_test_r(ScalarAdvDiff_order4              "RegTests/ScalarAdvDiff/*/erf_scalar_advdiff.exe" "plt00020")
add_test_r(ScalarAdvDiff_order5              "RegTests/ScalarAdvDiff/*/erf_scalar_advdiff.exe" "plt00020")
add_test_r(ScalarAdvDiff_order6              "RegTests/ScalarAdvDiff/*/erf_scalar_advdiff.exe" "plt00020")
add_test_r(ScalarAdvDiff_weno3               "RegTests/ScalarAdvDiff/*/erf_scalar_advdiff.exe" "plt00020")
add_test_r(ScalarAdvDiff_weno3z              "RegTests/ScalarAdvDiff/*/erf_scalar_advdiff.exe" "plt00020")
add_test_r(ScalarAdvDiff_weno5               "RegTests/ScalarAdvDiff/*/erf_scalar_advdiff.exe" "plt00020")
add_test_r(ScalarAdvDiff_weno5z              "RegTests/ScalarAdvDiff/*/erf_scalar_advdiff.exe" "plt00020")
add_test_r(ScalarAdvDiff_wenomzq3            "RegTests/ScalarAdvDiff/*/erf_scalar_advdiff.exe" "plt00020")
add_test_r(ScalarDiffusionGaussian           "RegTests/ScalarAdvDiff/*/erf_scalar_advdiff.exe" "plt00020")
add_test_r(ScalarDiffusionSine               "RegTests/ScalarAdvDiff/*/erf_scalar_advdiff.exe" "plt00020")
add_test_r(TaylorGreenAdvecting              "RegTests/TaylorGreenVortex/*/erf_taylor_green.exe" "plt00010")
add_test_r(TaylorGreenAdvectingDiffusing     "RegTests/TaylorGreenVortex/*/erf_taylor_green.exe" "plt00010")
add_test_r(MSF_NoSub_IsentropicVortexAdv     "RegTests/IsentropicVortex/*/erf_isentropic_vortex.exe" "plt00010")
add_test_r(MSF_Sub_IsentropicVortexAdv       "RegTests/IsentropicVortex/*/erf_isentropic_vortex.exe" "plt00010")
add_test_r(ABL_MOST                          "ABL/*/erf_abl.exe" "plt00010")
add_test_r(MoistBubble                       "RegTests/Bubble/*/erf_bubble.exe" "plt00010")
    
add_test_0(Deardorff_stationary              "ABL/*/erf_abl.exe" "plt00010")

else()
#add_test_r(Bubble_DensityCurrent             "Bubble/bubble" "plt00010")
add_test_r(CouetteFlow                       "RegTests/Couette_Poiseuille/erf_couette_poiseuille" "plt00050")
add_test_r(DensityCurrent                    "RegTests/DensityCurrent/erf_density_current" "plt00010")
add_test_r(DensityCurrent_detJ2              "RegTests/DensityCurrent/erf_density_current" "plt00010")
add_test_r(DensityCurrent_detJ2_nosub        "RegTests/DensityCurrent/erf_density_current" "plt00020")
add_test_r(DensityCurrent_detJ2_MT           "RegTests/DensityCurrent/erf_density_current" "plt00010")
add_test_r(EkmanSpiral                       "RegTests/EkmanSpiral/erf_ekman_spiral" "plt00010")
add_test_r(IsentropicVortexStationary        "RegTests/IsentropicVortex/erf_isentropic_vortex" "plt00010")
add_test_r(IsentropicVortexAdvecting         "RegTests/IsentropicVortex/erf_isentropic_vortex" "plt00010")
add_test_r(IVA_NumDiff                       "RegTests/IsentropicVortex/erf_isentropic_vortex" "plt00010")
add_test_r(MovingTerrain_nosub               "DevTests/MovingTerrain/erf_moving_terrain"   "plt00020")
add_test_r(MovingTerrain_sub                 "DevTests/MovingTerrain/erf_moving_terrain"   "plt00010")
add_test_r(PoiseuilleFlow                    "RegTests/Couette_Poiseuille/erf_couette_poiseuille" "plt00010")
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
add_test_r(TaylorGreenAdvecting              "RegTests/TaylorGreenVortex/erf_taylor_green" "plt00010")
add_test_r(TaylorGreenAdvectingDiffusing     "RegTests/TaylorGreenVortex/erf_taylor_green" "plt00010")
add_test_r(MSF_NoSub_IsentropicVortexAdv     "RegTests/IsentropicVortex/erf_isentropic_vortex" "plt00010")
add_test_r(MSF_Sub_IsentropicVortexAdv       "RegTests/IsentropicVortex/erf_isentropic_vortex" "plt00010")
add_test_r(ABL_MOST                          "ABL/erf_abl" "plt00010")
add_test_r(MoistBubble                       "RegTests/Bubble/erf_bubble" "plt00010")

add_test_0(InitSoundingIdeal_stationary      "ABL/erf_abl" "plt00010")
add_test_0(Deardorff_stationary              "ABL/erf_abl" "plt00010")
endif()
#=============================================================================
# Performance tests
#=============================================================================

