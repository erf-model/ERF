#include <cmath>
#include <cstdint>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Gpu.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#include "ERF_GTestHashRngCommon.H"

using namespace hash_rng_test;

namespace hash_rng_test {

struct LocatedMismatch
{
    int rank{-1};
    int i{-1};
    int j{-1};
    int k{-1};
    amrex::Real actual{amrex::Real(0.0)};
    amrex::Real expected{amrex::Real(0.0)};
    amrex::Real absolute_error{amrex::Real(0.0)};
};

void PrintTo (const LocatedMismatch& mismatch, std::ostream* os)
{
    *os << "LocatedMismatch{rank=" << mismatch.rank
        << ", cell=(" << mismatch.i << "," << mismatch.j << "," << mismatch.k << ")"
        << ", actual=" << mismatch.actual
        << ", expected=" << mismatch.expected
        << ", absolute_error=" << mismatch.absolute_error << "}";
}

} // namespace hash_rng_test

namespace {

enum class LayoutKind { SingleBox, FixedMaxSize, RankAdaptive };

struct LayoutCase
{
    std::string name;
    LayoutKind kind;
    amrex::IntVect max_size;
};

std::vector<LayoutCase> make_layout_cases ()
{
    return {
        {"single-box", LayoutKind::SingleBox, amrex::IntVect(AMREX_D_DECL(23, 19, 5))},
        {"maxsize-7x5x3", LayoutKind::FixedMaxSize, amrex::IntVect(AMREX_D_DECL(7, 5, 3))},
        {"maxsize-4x11x2", LayoutKind::FixedMaxSize, amrex::IntVect(AMREX_D_DECL(4, 11, 2))},
        {"rank-adaptive", LayoutKind::RankAdaptive, amrex::IntVect(AMREX_D_DECL(1, 1, 1))}
    };
}

// The harness launches this binary at every rank count in
// ERF_PARALLEL_TEST_NRANKS, and on a cluster that becomes `srun -n` / `flux run
// -n` with a value we do not choose.  The fixed layouts above deliberately keep
// awkward box boundaries, but "single-box" has exactly one box, so above one rank
// most ranks own nothing and the per-rank coverage check is vacuous for them.
// This layout chops the domain until it holds at least as many boxes as there are
// ranks, so whatever -n the launcher picks, every rank owns work and the
// invariance assertion is exercised on all of them.
amrex::BoxArray make_rank_adaptive_boxarray (const amrex::Box& domain,
                                             const int nprocs)
{
    amrex::IntVect max_size = domain.size();
    amrex::BoxArray ba(domain);
    ba.maxSize(max_size);
    while (ba.size() < nprocs) {
        int longest = 0;
        for (int dir = 1; dir < AMREX_SPACEDIM; ++dir) {
            if (max_size[dir] > max_size[longest]) {
                longest = dir;
            }
        }
        if (max_size[longest] <= 1) {
            break; // domain cannot be chopped any finer
        }
        max_size[longest] = (max_size[longest] + 1) / 2;
        ba = amrex::BoxArray(domain);
        ba.maxSize(max_size);
    }
    return ba;
}

amrex::BoxArray make_layout_boxarray (const amrex::Box& domain,
                                      const LayoutCase& layout,
                                      const int nprocs)
{
    if (layout.kind == LayoutKind::RankAdaptive) {
        return make_rank_adaptive_boxarray(domain, nprocs);
    }
    amrex::BoxArray ba(domain);
    if (layout.kind == LayoutKind::FixedMaxSize) {
        ba.maxSize(layout.max_size);
    }
    return ba;
}

std::string max_size_string (const LayoutCase& layout)
{
    if (layout.kind == LayoutKind::SingleBox) {
        return "domain";
    }
    if (layout.kind == LayoutKind::RankAdaptive) {
        return "adaptive";
    }
    std::ostringstream os;
    os << "(" << layout.max_size[0] << ","
       << layout.max_size[1] << ","
       << layout.max_size[2] << ")";
    return os.str();
}

// Keep the extended device lambda outside GoogleTest's private TestBody so
// CUDA extended-lambda diagnostics do not reject the test.
void fill_hash_field (amrex::MultiFab& field)
{
    for (amrex::MFIter mfi(field, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        const auto values = field.array(mfi);
        amrex::ParallelFor(
            bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                values(i,j,k) = erf_hash_rng::hash_uniform(
                    i, j, k, kTemperatureComp, 2, kTestSeed);
            });
    }
    amrex::Gpu::streamSynchronize();
}

amrex::MultiFab make_host_mirror (const amrex::MultiFab& source)
{
    amrex::MultiFab host(source.boxArray(),
                         source.DistributionMap(),
                         source.nComp(),
                         source.nGrowVect(),
                         amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
    amrex::MultiFab::Copy(host, source, 0, 0, source.nComp(), source.nGrowVect());
    amrex::Gpu::streamSynchronize();
    return host;
}

} // namespace

// Motivation: every value is keyed only by its global cell index, component,
// level, and seed. Distributed ownership and awkward box boundaries must not
// change the field at any MPI rank count registered by the test harness.
TEST(HashRngParallel, DistributedLayoutsMatchIndexReference)
{
    // Deliberately awkward extents so that no max_size below divides the domain
    // evenly.  No amrex::Geometry is constructed: the hash is keyed purely on the
    // global cell index, so the test needs no physical coordinates -- and the
    // Box-only Geometry constructor would read geometry.prob_hi/prob_extent from
    // ParmParse, which the unit-test harness does not set.
    const amrex::Box domain(amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
                            amrex::IntVect(AMREX_D_DECL(22, 18, 4)));
    const int rank = amrex::ParallelDescriptor::MyProc();
    const int nprocs = amrex::ParallelDescriptor::NProcs();

    for (const LayoutCase& layout : make_layout_cases()) {
        const amrex::BoxArray ba = make_layout_boxarray(domain, layout, nprocs);
        const amrex::DistributionMapping dm(ba);
        amrex::MultiFab field(ba, dm, 1, 0);
        fill_hash_field(field);
        const amrex::MultiFab field_host = make_host_mirror(field);

        amrex::Long local_checked = 0;
        amrex::Long local_mismatches = 0;
        LocatedMismatch local_worst;
        for (amrex::MFIter mfi(field_host); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.validbox();
            const auto values = field_host.const_array(mfi);
            for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
                for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                    for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                        ++local_checked;
                        const amrex::Real actual = values(i,j,k);
                        const amrex::Real expected = erf_hash_rng::hash_uniform(
                            i, j, k, kTemperatureComp, 2, kTestSeed);
                        if (actual != expected) {
                            ++local_mismatches;
                            const amrex::Real error = std::abs(actual - expected);
                            if (error > local_worst.absolute_error) {
                                local_worst = LocatedMismatch{
                                    rank, i, j, k, actual, expected, error};
                            }
                        }
                    }
                }
            }
        }

        amrex::Long global_checked = local_checked;
        amrex::Long global_mismatches = local_mismatches;
        amrex::ParallelDescriptor::ReduceLongSum(global_checked);
        amrex::ParallelDescriptor::ReduceLongSum(global_mismatches);

        // Anti-vacuity guard.  The global reduction below proves the whole domain
        // was covered, but says nothing about whether THIS rank checked anything:
        // a layout with fewer boxes than ranks leaves ranks idle, and at a large
        // -n the test would report success while verifying nothing on most of
        // them.  The rank-adaptive layout exists to prevent that, so assert both
        // that its chopping actually reached the rank count and that this rank
        // consequently owns work.  The sole precondition is the domain's own
        // capacity -- deliberately NOT "enough boxes were produced", which would
        // make the guard true only when it is already trivially satisfied.
        if (layout.kind == LayoutKind::RankAdaptive &&
            static_cast<amrex::Long>(nprocs) <= domain.numPts()) {
            EXPECT_GE(ba.size(), static_cast<amrex::Long>(nprocs))
                << "layout=" << layout.name
                << " rank=" << rank
                << " nprocs=" << nprocs
                << " nboxes=" << ba.size();
            EXPECT_GT(local_checked, amrex::Long(0))
                << "layout=" << layout.name
                << " rank=" << rank
                << " nprocs=" << nprocs
                << " nboxes=" << ba.size();
        }

        EXPECT_EQ(global_checked, domain.numPts())
            << "layout=" << layout.name
            << " max_size=" << max_size_string(layout)
            << " rank=" << rank
            << " nprocs=" << nprocs
            << " local_checked=" << local_checked;
        EXPECT_EQ(global_mismatches, amrex::Long(0))
            << "layout=" << layout.name
            << " max_size=" << max_size_string(layout)
            << " rank=" << rank
            << " nprocs=" << nprocs
            << " local_mismatches=" << local_mismatches
            << " local_worst=" << ::testing::PrintToString(local_worst);
    }
}
