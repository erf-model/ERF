#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#include <cstdlib>

TEST(ParallelEnvironment, HelloWorld)
{
    const int my_rank = amrex::ParallelDescriptor::MyProc();
    const int n_ranks = amrex::ParallelDescriptor::NProcs();

    EXPECT_GE(my_rank, 0);
    EXPECT_LT(my_rank, n_ranks);
    EXPECT_GE(n_ranks, 1);

    if (const char* expected = std::getenv("ERF_PARALLEL_TEST_NPROCS")) {
        char* end = nullptr;
        const long expected_nprocs = std::strtol(expected, &end, 10);
        ASSERT_NE(end, expected);
        ASSERT_EQ(*end, '\0');
        EXPECT_GT(expected_nprocs, 0);
        EXPECT_EQ(n_ranks, expected_nprocs);
    }

    int local_val = 1;
    amrex::ParallelDescriptor::ReduceIntSum(local_val);

    EXPECT_EQ(local_val, n_ranks);
}