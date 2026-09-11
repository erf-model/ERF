#ifndef ERF_GTEST_MPI_LISTENER_H_
#define ERF_GTEST_MPI_LISTENER_H_

#include <mpi.h>

#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

namespace erf_gtest {

struct MPIFailureRecord {
    int rank{0};
    std::string suite;
    std::string test;
    std::string file;
    int line{0};
    bool fatal{false};
    std::string message;
};

class MPIFailureListener final : public ::testing::EmptyTestEventListener {
public:
    explicit MPIFailureListener (const int rank) : m_rank(rank) {}

    void OnTestStart (const ::testing::TestInfo& test_info) override
    {
        // Do not query UnitTest::current_test_info() from OnTestPartResult:
        // GoogleTest may hold its internal event mutex while invoking the
        // callback, and querying the current test can recursively acquire it.
        m_suite = test_info.test_suite_name();
        m_test = test_info.name();
    }

    void OnTestPartResult (const ::testing::TestPartResult& result) override
    {
        if (!result.failed()) {
            return;
        }

        MPIFailureRecord record;
        record.rank = m_rank;
        record.suite = m_suite;
        record.test = m_test;
        record.file = result.file_name() != nullptr ? result.file_name() : "unknown";
        record.line = result.line_number();
        record.fatal = result.fatally_failed();
        record.message = result.summary();
        m_failures.push_back(std::move(record));
    }

    const std::vector<MPIFailureRecord>& failures () const { return m_failures; }

private:
    int m_rank;
    std::string m_suite;
    std::string m_test;
    std::vector<MPIFailureRecord> m_failures;
};

inline std::string serialize_failures (const std::vector<MPIFailureRecord>& failures)
{
    std::ostringstream output;
    for (const auto& failure : failures) {
        output << "rank=" << failure.rank << '\n'
               << "suite=" << failure.suite << '\n'
               << "test=" << failure.test << '\n'
               << "file=" << failure.file << '\n'
               << "line=" << failure.line << '\n'
               << "fatal=" << (failure.fatal ? "true" : "false") << '\n'
               << "message=" << failure.message << '\n'
               << "---\n";
    }
    return output.str();
}

inline void gather_and_report_failures (const MPI_Comm comm,
                                         const MPIFailureListener& listener)
{
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const std::string local_payload = serialize_failures(listener.failures());
    const int local_size = static_cast<int>(local_payload.size());
    std::vector<int> sizes(rank == 0 ? size : 0);
    MPI_Gather(&local_size, 1, MPI_INT,
               rank == 0 ? sizes.data() : nullptr, 1, MPI_INT, 0, comm);

    std::vector<int> displacements(rank == 0 ? size : 0);
    int total_size = 0;
    if (rank == 0) {
        for (int source = 0; source < size; ++source) {
            displacements[source] = total_size;
            total_size += sizes[source];
        }
    }

    std::vector<char> all_payload(rank == 0 ? static_cast<std::size_t>(total_size) : 0);
    MPI_Gatherv(local_payload.data(), local_size, MPI_CHAR,
                rank == 0 && total_size > 0 ? all_payload.data() : nullptr,
                rank == 0 ? sizes.data() : nullptr,
                rank == 0 ? displacements.data() : nullptr,
                MPI_CHAR, 0, comm);

    if (rank == 0 && total_size > 0) {
        std::cout << "\n=== ERF GoogleTest MPI failure report ===\n";
        for (int source = 0; source < size; ++source) {
            if (sizes[source] == 0 || source == 0) {
                continue;
            }
            std::cout << "--- failures collected from rank " << source << " ---\n";
            std::cout.write(all_payload.data() + displacements[source], sizes[source]);
        }
        std::cout << "=== end ERF GoogleTest MPI failure report ===\n";
        std::cout.flush();
    }
}

} // namespace erf_gtest

#endif
