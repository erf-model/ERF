#include <cmath>
#include <initializer_list>
#include <string>
#include <utility>

#include <AMReX_BoxArray.H>
#include <AMReX_BoxList.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>

#include <gtest/gtest.h>

#include "ERF_Utils.H"

namespace {

// An 8 x 8 x 24 cell domain, 800 m x 800 m x 1200 m, over a 60 m bump
constexpr int nx = 8;
constexpr int ny = 8;
constexpr int nz = 24;

// Terrain-following node heights: the bump flattens out towards the top of the domain
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real z_node (int i, int j, int k) noexcept
{
    const amrex::Real xi = amrex::Real(i) - amrex::Real(0.5)*amrex::Real(nx);
    const amrex::Real yj = amrex::Real(j) - amrex::Real(0.5)*amrex::Real(ny);
    const amrex::Real h  = amrex::Real(60.0) * std::exp(-(xi*xi + yj*yj) / amrex::Real(8.0));
    return h + (amrex::Real(1200.0) - h) * amrex::Real(k) / amrex::Real(nz);
}

struct ColumnState
{
    amrex::MultiFab rho;
    amrex::MultiFab theta;
    amrex::MultiFab qv;
    amrex::MultiFab qt;
    amrex::MultiFab z_nd;
};

// A moist, stably stratified state on the given boxes whose density is several percent out
// of hydrostatic balance, varying in k, so that where an integration starts matters
ColumnState make_state (const amrex::BoxArray& ba)
{
    const amrex::DistributionMapping dm(ba);
    ColumnState s;
    s.rho.define  (ba, dm, 1, 0);
    s.theta.define(ba, dm, 1, 0);
    s.qv.define   (ba, dm, 1, 0);
    s.qt.define   (ba, dm, 1, 0);
    s.z_nd.define (amrex::convert(ba, amrex::IntVect(1,1,1)), dm, 1, 0);

    for (amrex::MFIter mfi(s.rho); mfi.isValid(); ++mfi) {
        const auto z   = s.z_nd.array(mfi);
        const auto rho = s.rho.array(mfi);
        const auto th  = s.theta.array(mfi);
        const auto qv  = s.qv.array(mfi);
        const auto qt  = s.qt.array(mfi);
        amrex::ParallelFor(s.z_nd[mfi].box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            z(i,j,k) = z_node(i,j,k);
        });
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real zc = amrex::Real(0.125) *
                (z_node(i,j,k  ) + z_node(i+1,j,k  ) + z_node(i,j+1,k  ) + z_node(i+1,j+1,k  ) +
                 z_node(i,j,k+1) + z_node(i+1,j,k+1) + z_node(i,j+1,k+1) + z_node(i+1,j+1,k+1));
            th (i,j,k) = amrex::Real(300.0) + amrex::Real(0.003)*zc
                       + amrex::Real(0.2)*std::cos(amrex::Real(0.7)*amrex::Real(i) + amrex::Real(0.3)*amrex::Real(j));
            qv (i,j,k) = amrex::Real(0.008) * std::exp(-zc / amrex::Real(1500.0));
            qt (i,j,k) = qv(i,j,k) + amrex::Real(0.0005);
            rho(i,j,k) = amrex::Real(1.2) * std::exp(-zc / amrex::Real(9000.0))
                       * (amrex::Real(1.0) + amrex::Real(0.03)*std::sin(amrex::Real(0.9)*amrex::Real(k) + amrex::Real(0.4)*amrex::Real(i)));
        });
    }
    return s;
}

amrex::Geometry make_geom ()
{
    const amrex::Box domain(amrex::IntVect(0,0,0), amrex::IntVect(nx-1,ny-1,nz-1));
    const amrex::RealBox rb(amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0),
                            amrex::Real(800.0), amrex::Real(800.0), amrex::Real(1200.0));
    return amrex::Geometry(domain, rb, amrex::CoordSys::cartesian, {1, 1, 0});
}

void rebalance (ColumnState& s, bool maintain_Th, bool use_sfc)
{
    const amrex::Geometry geom = make_geom();
    rebalance_columns(s.rho, s.theta, s.qv, s.qt, &s.z_nd, geom, maintain_Th, use_sfc);
}

// Largest |other - ref| over the cells of ref, relative to the largest |ref|
amrex::Real max_rel_difference (const amrex::MultiFab& ref, const amrex::MultiFab& other)
{
    amrex::MultiFab diff(ref.boxArray(), ref.DistributionMap(), 1, 0);
    diff.ParallelCopy(other, 0, 0, 1);
    amrex::MultiFab::Subtract(diff, ref, 0, 0, 1, 0);
    return diff.norminf() / ref.norminf();
}

amrex::BoxArray chopped (const amrex::Box& region, const amrex::IntVect& max_size)
{
    amrex::BoxArray ba(region);
    ba.maxSize(max_size);
    return ba;
}

// Boxes stacked in z whose footprints differ between the lower and the upper boxes
amrex::BoxArray staggered_footprints ()
{
    amrex::BoxList bl;
    bl.push_back(amrex::Box(amrex::IntVect(0,0, 0), amrex::IntVect(3,7, 9)));
    bl.push_back(amrex::Box(amrex::IntVect(4,0, 0), amrex::IntVect(7,7, 9)));
    bl.push_back(amrex::Box(amrex::IntVect(0,0,10), amrex::IntVect(7,3,23)));
    bl.push_back(amrex::Box(amrex::IntVect(0,4,10), amrex::IntVect(7,7,17)));
    bl.push_back(amrex::Box(amrex::IntVect(0,4,18), amrex::IntVect(7,7,23)));
    return amrex::BoxArray(std::move(bl));
}

// The same arithmetic runs whatever the boxes, so the results agree to round-off
amrex::Real tolerance ()
{
    return (sizeof(amrex::Real) == 8) ? amrex::Real(1.0e-13) : amrex::Real(1.0e-6);
}

void expect_layouts_match_one_box (const amrex::Box& region, bool maintain_Th, bool use_sfc)
{
    ColumnState ref = make_state(amrex::BoxArray(region));
    rebalance(ref, maintain_Th, use_sfc);

    // The rebalance must have changed the density, or the comparison proves nothing
    ColumnState untouched = make_state(amrex::BoxArray(region));
    ASSERT_GT(max_rel_difference(ref.rho, untouched.rho), amrex::Real(1.0e-3));

    const std::pair<const char*, amrex::BoxArray> layouts[] = {
        {"4 columns",                    chopped(region, amrex::IntVect(4,4,nz))},
        {"3 boxes per column",           chopped(region, amrex::IntVect(8,8,8))},
        {"4 columns of 3 boxes",         chopped(region, amrex::IntVect(4,4,8))},
        {"4 columns of 5 uneven boxes",  chopped(region, amrex::IntVect(4,4,5))},
    };

    for (const auto& layout : layouts) {
        SCOPED_TRACE(layout.first);
        ColumnState split = make_state(layout.second);
        rebalance(split, maintain_Th, use_sfc);
        EXPECT_LE(max_rel_difference(ref.rho,   split.rho),   tolerance());
        EXPECT_LE(max_rel_difference(ref.theta, split.theta), tolerance());
    }

    if (region == make_geom().Domain()) {
        SCOPED_TRACE("staggered footprints");
        ColumnState split = make_state(staggered_footprints());
        rebalance(split, maintain_Th, use_sfc);
        EXPECT_LE(max_rel_difference(ref.rho,   split.rho),   tolerance());
        EXPECT_LE(max_rel_difference(ref.theta, split.theta), tolerance());
    }
}

} // namespace

// A column split into boxes stacked in z continues its integration across the split, so the
// rebalanced state is the one a single box gives; each box used to start afresh from its own
// lowest, unbalanced cell
TEST(RebalanceColumnsParallel, SplitInZMatchesOneBox)
{
    for (const bool maintain_Th : {true, false}) {
        SCOPED_TRACE(std::string("maintain_Th ") + (maintain_Th ? "true" : "false"));
        expect_layouts_match_one_box(make_geom().Domain(), maintain_Th, false);
    }
}

// Seeded from p_0 at the surface, a split column continues from the box below
TEST(RebalanceColumnsParallel, SurfaceSeededSplitInZMatchesOneBox)
{
    for (const bool maintain_Th : {true, false}) {
        SCOPED_TRACE(std::string("maintain_Th ") + (maintain_Th ? "true" : "false"));
        expect_layouts_match_one_box(make_geom().Domain(), maintain_Th, true);
    }
}

// Boxes that do not reach the bottom of the domain, as on a refined patch aloft: the lowest
// boxes start from their own lowest cell, and the boxes stacked on them continue from there
TEST(RebalanceColumnsParallel, PatchAloftMatchesOneBox)
{
    const amrex::Box aloft(amrex::IntVect(0,0,8), amrex::IntVect(nx-1,ny-1,nz-1));
    expect_layouts_match_one_box(aloft, true, false);
}
