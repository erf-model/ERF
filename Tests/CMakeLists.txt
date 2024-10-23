# Additional testing options
set(ERF_TEST_NRANKS 2 CACHE STRING  "Number of MPI ranks to use for each test")
set(ERF_TEST_FCOMPARE_RTOL "2.0e-10" CACHE STRING "fcompare relative tolerance")
set(ERF_TEST_FCOMPARE_ATOL "2.0e-10" CACHE STRING "fcompare absolute tolerance")

message(STATUS "ERF testing configuration summary:")
message(STATUS "   Number of ranks               = ${ERF_TEST_NRANKS}")
message(STATUS "   fcompare executable           = ${FCOMPARE_EXE}")
message(STATUS "   comparison relative tolerance = ${ERF_TEST_FCOMPARE_RTOL}")
message(STATUS "   comparison absolute tolerance = ${ERF_TEST_FCOMPARE_ATOL}")

include(${CMAKE_CURRENT_SOURCE_DIR}/CTestList.cmake)
