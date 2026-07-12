#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ccse-mpi.H>

#include <gtest/gtest.h>

#include <cstdlib>
#include <string>
#include <vector>

#include "ERF_GTestMPIListener.H"

namespace {

bool gtest_list_requested (const int argc, char** argv)
{
    for (int index = 1; index < argc; ++index) {
        if (std::string(argv[index]) == "--gtest_list_tests" ||
            std::string(argv[index]).rfind("--gtest_list_tests=", 0) == 0) {
            return true;
        }
    }
    return false;
}

void make_mpi_xml_path_unique (std::vector<std::string>& arguments,
                               const int rank)
{
    const char* nprocs = std::getenv("ERF_PARALLEL_TEST_NPROCS");
    if (nprocs == nullptr) {
        return;
    }

    for (std::string& argument : arguments) {
        if (argument.rfind("--gtest_output=", 0) != 0) {
            continue;
        }
        const std::string output = argument.substr(std::string("--gtest_output=").size());
        const std::size_t separator = output.find(':');
        if (separator == std::string::npos) {
            continue;
        }
        argument = "--gtest_output=" + output.substr(0, separator + 1) +
                   output.substr(separator + 1) + ".np" + nprocs + ".rank" +
                   std::to_string(rank) + ".xml";
        return;
    }
}

} // namespace

int
main (int argc, char** argv)
{
    const std::vector<std::string> original_arguments(argv, argv + argc);

    // Discovery must not initialize MPI or AMReX. This mirrors the serial
    // main and avoids launcher/runtime initialization during CTest discovery.
    if (gtest_list_requested(argc, argv)) {
        ::testing::InitGoogleTest(&argc, argv);
        const int list_result = RUN_ALL_TESTS();
        std::fflush(nullptr);
        std::_Exit(list_result);
    }

    int provided = MPI_THREAD_SINGLE;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    (void) provided;

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // MPI implementations are permitted to rewrite argc/argv during
    // initialization. Keep the GoogleTest arguments captured before MPI so
    // XML output and filters are not lost.
    std::vector<std::string> argument_storage = original_arguments;
    make_mpi_xml_path_unique(argument_storage, rank);
    std::vector<char*> argument_values;
    argument_values.reserve(argument_storage.size() + 1);
    for (std::string& argument : argument_storage) {
        argument_values.push_back(argument.data());
    }
    argument_values.push_back(nullptr);
    int gtest_argc = static_cast<int>(argument_storage.size());
    char** gtest_argv = argument_values.data();
    ::testing::InitGoogleTest(&gtest_argc, gtest_argv);

    amrex::Initialize(gtest_argc, gtest_argv, true, MPI_COMM_WORLD);

    auto& listeners = ::testing::UnitTest::GetInstance()->listeners();

    if (rank != 0 && std::getenv("ERF_GTEST_VERBOSE_RANKS") == nullptr) {
        delete listeners.Release(listeners.default_result_printer());
    }
    auto* listener = new erf_gtest::MPIFailureListener(rank);
    listeners.Append(listener);

    const int local_result = RUN_ALL_TESTS();

    erf_gtest::gather_and_report_failures(MPI_COMM_WORLD, *listener);

    int global_result = local_result;
    MPI_Allreduce(MPI_IN_PLACE, &global_result, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);

    amrex::Finalize();
    MPI_Finalize();

    return global_result;
}
