#include <mpi.h>

#include <cstdlib>

#include <gtest/gtest.h>

#include "ERF_GTestMPIListener.H"

TEST(ERFMPIListenerInfrastructure, NonRootFailureIsReported)
{
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 1) {
        ADD_FAILURE() << "intentional non-root failure";
    }
}

int main (int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    int provided = MPI_THREAD_SINGLE;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    (void) provided;

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    auto& listeners = ::testing::UnitTest::GetInstance()->listeners();
    if (rank != 0 && std::getenv("ERF_GTEST_VERBOSE_RANKS") == nullptr) {
        delete listeners.Release(listeners.default_result_printer());
    }
    auto* listener = new erf_gtest::MPIFailureListener(rank);
    listeners.Append(listener);

    const int local_result = RUN_ALL_TESTS();
    erf_gtest::gather_and_report_failures(MPI_COMM_WORLD, *listener);

    int global_result = local_result;
    MPI_Allreduce(MPI_IN_PLACE, &global_result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Finalize();
    return global_result;
}
