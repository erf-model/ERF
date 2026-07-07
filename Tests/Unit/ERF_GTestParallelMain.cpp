#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ccse-mpi.H>

#include <gtest/gtest.h>

#include <cstdlib>

int
main (int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    int provided = MPI_THREAD_SINGLE;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    amrex::Initialize(argc, argv, true, MPI_COMM_WORLD);

    auto& listeners = ::testing::UnitTest::GetInstance()->listeners();

    if (rank != 0 && std::getenv("ERF_GTEST_VERBOSE_RANKS") == nullptr) {
        delete listeners.Release(listeners.default_result_printer());
    }

    const int local_result = RUN_ALL_TESTS();

    int global_result = local_result;
    MPI_Allreduce(MPI_IN_PLACE, &global_result, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);

    amrex::Finalize();
    MPI_Finalize();

    return global_result;
}