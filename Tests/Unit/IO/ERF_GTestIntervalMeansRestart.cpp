#include <filesystem>
#include <string>
#include <utility>

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_VisMF.H>

#include <gtest/gtest.h>

namespace {

constexpr const char* restart_directory = "erf_unit_interval_means_restart";

void
fill_spatial_field (amrex::MultiFab& field)
{
    for (amrex::MFIter mfi(field); mfi.isValid(); ++mfi) {
        const auto box = mfi.validbox();
        auto array = field.array(mfi);
        amrex::ParallelFor(box, field.nComp(),
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                array(i, j, k, n) = amrex::Real(1000 * n + 100 * k + 10 * j + i);
            });
    }
    amrex::Gpu::streamSynchronize();
}

amrex::DistributionMapping
make_mapping (const amrex::BoxArray& boxes, bool reverse)
{
    amrex::Vector<int> processors(boxes.size());
    const int nprocs = amrex::ParallelDescriptor::NProcs();
    for (int ibox = 0; ibox < boxes.size(); ++ibox) {
        const int slot = ibox % nprocs;
        processors[ibox] = reverse ? (nprocs - 1 - slot) : slot;
    }
    return amrex::DistributionMapping(std::move(processors));
}

} // namespace

// VisMF reads the checkpoint BoxArray and redistributes its data into the
// caller's DistributionMapping. This is the operation used by ERF's interval-
// mean restart path when the MPI decomposition changes.
TEST(IntervalMeansRestart, RedistributesSpatiallyVaryingComponents)
{
    const std::filesystem::path directory(restart_directory);
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::filesystem::remove_all(directory);
        std::filesystem::create_directories(directory);
    }
    amrex::ParallelDescriptor::Barrier();

    const amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(7, 7, 3));
    amrex::BoxArray boxes(domain);
    boxes.maxSize(2);
    const auto write_dm = make_mapping(boxes, false);
    const auto read_dm = make_mapping(boxes, true);
    const std::string field_name = (directory / "IntervalMeans").string();

    amrex::MultiFab written(boxes, write_dm, 2, 0);
    fill_spatial_field(written);
    amrex::VisMF::Write(written, field_name);

    amrex::MultiFab restored(boxes, read_dm, 2, 0);
    amrex::VisMF::Read(restored, field_name);

    amrex::MultiFab expected(boxes, read_dm, 2, 0);
    fill_spatial_field(expected);
    amrex::MultiFab::Subtract(restored, expected, 0, 0, 2, 0);
    EXPECT_EQ(restored.norm0(0, 1, false), amrex::Real(0.0));

    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::filesystem::remove_all(directory);
    }
    amrex::ParallelDescriptor::Barrier();
}
