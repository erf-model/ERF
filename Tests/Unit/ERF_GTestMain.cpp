#include <AMReX.H>

#include <cstdio>
#include <cstdlib>
#include <gtest/gtest.h>

int
main (int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    // When invoked with --gtest_list_tests (gtest_discover_tests prints the
    // test list this way and CMake consumes it from stdout), skip
    // amrex::Initialize so we don't hang the MPI/GPU runtime, and use
    // std::_Exit to skip static destructors that can block on Cray MPICH
    // atexit handlers. Flush stdio first because _Exit bypasses normal stream
    // cleanup.
    if (::testing::GTEST_FLAG(list_tests))
    {
        const int list_result = RUN_ALL_TESTS();
        std::fflush(nullptr);
        std::_Exit(list_result);
    }

    amrex::Initialize(argc, argv);
    const int result = RUN_ALL_TESTS();
    amrex::Finalize();

    return result;
}
