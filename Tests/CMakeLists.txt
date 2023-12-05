
set(FCOMPARE_EXE ${CMAKE_BINARY_DIR}/Submodules/AMReX/Tools/Plotfile/amrex_fcompare CACHE INTERNAL "Path to fcompare executable for regression tests")
set(ERF_TEST_NRANKS 4 CACHE STRING  "Number of MPI ranks to use for each test")
include(${CMAKE_CURRENT_SOURCE_DIR}/CTestList.cmake)
