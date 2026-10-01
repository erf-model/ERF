#ifndef ERF_GTEST_TEMP_DIR_H_
#define ERF_GTEST_TEMP_DIR_H_

#include <filesystem>
#include <random>
#include <string>

#include <AMReX_INT.H>
#include <AMReX_ParallelDescriptor.H>

/**
 * @brief A fresh scratch path under the system temporary directory, drawn anew on every
 * call (two calls with the same stem give two different paths), and the same on every
 * rank of the calling process. Store the result; do not call again expecting the same
 * directory.
 *
 * ctest -j runs several test processes at once -- each discovered test, and
 * erf_unit_tests_shuffle_repeat, which runs all of them again -- and two build trees
 * share one temporary directory. A fixed name there is shared too: one process removes
 * or rewrites the checkpoint another is still reading, and the reader aborts. The tag is
 * drawn on the I/O rank and broadcast, so a collective write still sees one path.
 */
inline std::filesystem::path erf_gtest_temp_path (const std::string& stem)
{
    amrex::Long tag = 0;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::random_device rd;
        const unsigned long long bits =
            (static_cast<unsigned long long>(rd()) << 32) ^ static_cast<unsigned long long>(rd());
        tag = static_cast<amrex::Long>(bits >> 1);  // non-negative
    }
    amrex::ParallelDescriptor::Bcast(&tag, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    return std::filesystem::temp_directory_path() / (stem + "_" + std::to_string(tag));
}

#endif // ERF_GTEST_TEMP_DIR_H_
