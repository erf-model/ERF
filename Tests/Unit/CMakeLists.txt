set(BUILD_GMOCK OFF CACHE BOOL
    "Build Google Mock targets for the vendored GoogleTest dependency")
set(INSTALL_GTEST OFF CACHE BOOL
    "Install targets from the vendored GoogleTest dependency")
if(MSVC)
  set(gtest_force_shared_crt ON CACHE BOOL
      "Use the shared CRT when building the vendored GoogleTest dependency")
endif()

add_subdirectory(${PROJECT_SOURCE_DIR}/Submodules/googletest
                 ${PROJECT_BINARY_DIR}/Submodules/googletest
                 EXCLUDE_FROM_ALL)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|IntelLLVM")
  foreach(_erf_gtest_target gtest gtest_main)
    if(TARGET ${_erf_gtest_target})
      target_compile_options(${_erf_gtest_target}
        PRIVATE
          -Wno-unknown-warning-option
      )
    endif()
  endforeach()
endif()

set(erf_lib_name erf_srclib)
set(erf_exe_name erf_unit_tests)

add_executable(${erf_exe_name})
target_sources(${erf_exe_name}
  PRIVATE
    ${PROJECT_SOURCE_DIR}/Exec/ERF_Prob.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ERF_GTestMain.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ERF_GTestHelloWorld.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ERF_GTestRuntimeStubs.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/EOS/ERF_GTestEOSScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/EOS/ERF_GTestEOSKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/Interpolation/ERF_GTestInterpolationScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/Interpolation/ERF_GTestInterpolationKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/Interpolation/ERF_GTestInterpolationPolicyKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/Interpolation/ERF_GTestInterpolationWENOScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/Interpolation/ERF_GTestInterpolationWENOKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/Microphysics/ERF_GTestMicrophysicsUtilsScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils/Microphysics/ERF_GTestMicrophysicsUtilsKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SatAdj/ERF_GTestSatAdjThermo.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SatAdj/ERF_GTestSatAdjCell.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/SatAdj/ERF_GTestSatAdjMultiFab.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/Kessler/ERF_GTestKesslerScalar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/Kessler/ERF_GTestKesslerKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Microphysics/Kessler/ERF_GTestKesslerPhysicalProperties.cpp
)
build_erf_exe(${erf_exe_name})
target_link_libraries(${erf_exe_name} PRIVATE GTest::gtest)

include(GoogleTest)
gtest_discover_tests(${erf_exe_name}
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS "unit"
)
