#include <AMReX.H>
#include <gtest/gtest.h>

int
main (int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    amrex::Initialize(argc, argv);
    const int result = RUN_ALL_TESTS();
    amrex::Finalize();

    return result;
}
