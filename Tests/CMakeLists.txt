
set(FCOMPARE_EXE ${CMAKE_BINARY_DIR}/Submodules/AMReX/Tools/Plotfile/amrex_fcompare CACHE INTERNAL "Path to fcompare executable for regression tests")
include(${CMAKE_CURRENT_SOURCE_DIR}/CTestList.cmake)
