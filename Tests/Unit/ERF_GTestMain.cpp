#include <AMReX.H>
#include <gtest/gtest.h>

int
main (int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    // When invoked with --gtest_list_tests (gtest_discover_tests prints the
    // test list this way), skip amrex::Initialize so we don't hang the MPI/GPU
    // runtime, and use std::_Exit to skip static destructors that can block on
    // Cray MPICH atexit handlers.
    if (::testing::GTEST_FLAG(list_tests))
    {
        const int list_result = RUN_ALL_TESTS();
        std::_Exit(list_result);
    }

    amrex::Initialize(argc, argv);
    const int result = RUN_ALL_TESTS();
    amrex::Finalize();

    return result;
}
