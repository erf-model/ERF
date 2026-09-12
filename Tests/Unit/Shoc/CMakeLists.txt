if(NOT TARGET GTest::gtest)
  message(FATAL_ERROR
    "erf_shoc_unit_tests expects the vendored googletest target from "
    "Tests/Unit/CMakeLists.txt. Configure unit tests through the top-level "
    "ERF test build so Submodules/googletest is added first.")
endif()

include(${PROJECT_SOURCE_DIR}/CMake/BuildERFExe.cmake)
include(${PROJECT_SOURCE_DIR}/CMake/SetERFCompileFlags.cmake)

add_executable(erf_shoc_unit_tests
  ${PROJECT_SOURCE_DIR}/Exec/ERF_Prob.cpp
  ERF_ShocUnitMain.cpp
  ${PROJECT_SOURCE_DIR}/Tests/Unit/ERF_GTestRuntimeStubs.cpp
  ERF_ShocStructureTests.cpp
  ERF_ShocTkeTests.cpp
  ERF_ShocMomentsTests.cpp
  ERF_ShocPDFTests.cpp
  ERF_ShocImplicitTests.cpp
  ERF_ShocDriverTests.cpp
  ERF_ShocPhysicalPropertyTests.cpp
  ERF_ShocThermoTests.cpp
  ERF_ShocRadiationCouplingTests.cpp
)

build_erf_exe(erf_shoc_unit_tests NO_ERF_MAIN)

target_include_directories(erf_shoc_unit_tests PRIVATE
  ${CMAKE_SOURCE_DIR}/Source
  ${CMAKE_SOURCE_DIR}/Source/DataStructs
  ${CMAKE_SOURCE_DIR}/Source/PBL
  ${CMAKE_SOURCE_DIR}/Source/PBL/Shoc
  ${CMAKE_SOURCE_DIR}/Source/Utils
  ${CMAKE_SOURCE_DIR}/Source/PhysicsInterfaces/Radiation
)

target_compile_definitions(erf_shoc_unit_tests PRIVATE
  ERF_SHOC_UNIT_FIXTURE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/fixtures"
  ERF_HAS_SHOC_FAMILY
)

target_link_libraries(erf_shoc_unit_tests PRIVATE GTest::gtest)

erf_register_gtest_suite(erf_shoc_unit_tests
  LABELS unit shoc
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DISCOVERY_TIMEOUT 60
  EXECUTION_TIMEOUT 900
)
