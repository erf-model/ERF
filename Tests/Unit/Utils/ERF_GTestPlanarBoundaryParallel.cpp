#include <cmath>
#include <initializer_list>
#include <string>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Gpu.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Periodicity.H>
#include <AMReX_Reduce.H>

#include <gtest/gtest.h>

#include "ERF_PlanarBoundary.H"

namespace {

// A 16 x 16 x 12 cell domain, periodic in x and y as in the surface-layer cases
constexpr int nx = 16;
constexpr int ny = 16;
constexpr int nz = 12;

// Stands in for the bogus_large_value the planar fields are allocated with
constexpr amrex::Real placeholder = amrex::Real(1.0e30);

struct PlanarLayout
{
    amrex::BoxArray ba3d;
    amrex::BoxArray ba2d;
    amrex::DistributionMapping dm;
};

// The domain chopped into boxes of at most max_size cells, keeping the boxes whose lowest
// cell is at or above kmin, and its z-collapse, built the way SurfaceLayer and MOSTAverage
// build their planar BoxArrays
PlanarLayout make_layout (const amrex::IntVect& max_size, int kmin = 0)
{
    amrex::BoxArray chopped(amrex::Box(amrex::IntVect(0,0,0), amrex::IntVect(nx-1,ny-1,nz-1)));
    chopped.maxSize(max_size);

    amrex::BoxList bl3d;
    for (int ib = 0; ib < static_cast<int>(chopped.size()); ++ib) {
        if (chopped[ib].smallEnd(2) >= kmin) { bl3d.push_back(chopped[ib]); }
    }

    PlanarLayout layout;
    layout.ba3d = amrex::BoxArray(std::move(bl3d));
    layout.dm   = amrex::DistributionMapping(layout.ba3d);

    amrex::BoxList bl2d = layout.ba3d.boxList();
    for (auto& b : bl2d) { b.setRange(2,0); }
    layout.ba2d = amrex::BoxArray(std::move(bl2d));
    return layout;
}

// Value at (i,j) in component n, periodic in x and y
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real expected_value (int i, int j, int n, amrex::Real offset) noexcept
{
    const int iw = ((i % nx) + nx) % nx;
    const int jw = ((j % ny) + ny) % ny;
    return offset + amrex::Real(iw + 100*jw + 10000*n);
}

// Write the expected value into the valid region of the copies that belong to the 3D
// boxes touching the surface, the only copies the surface layer computes
void set_surface_copies (amrex::MultiFab& mf, const PlanarLayout& layout, amrex::Real offset)
{
    const int ncomp = mf.nComp();
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        if (layout.ba3d[mfi.index()].smallEnd(2) != 0) { continue; }
        const auto arr = mf.array(mfi);
        amrex::ParallelFor(mfi.validbox(), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            arr(i,j,k,n) = expected_value(i,j,n,offset);
        });
    }
}

// Largest |value - expected| over the valid region and ghost cells of every box on every rank
amrex::Real max_error (const amrex::MultiFab& mf, amrex::Real offset)
{
    amrex::ReduceOps<amrex::ReduceOpMax> reduce_op;
    amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    const int ncomp = mf.nComp();
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const auto arr = mf.const_array(mfi);
        reduce_op.eval(mfi.fabbox(), ncomp, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) -> ReduceTuple
            {
                return std::abs(arr(i,j,k,n) - expected_value(i,j,n,offset));
            });
    }
    amrex::Real err = amrex::get<0>(reduce_data.value(reduce_op));
    amrex::ParallelDescriptor::ReduceRealMax(err);
    return err;
}

// Largest |value - placeholder| over the valid region and ghost cells of every box on every rank
amrex::Real max_change_from_placeholder (const amrex::MultiFab& mf)
{
    amrex::ReduceOps<amrex::ReduceOpMax> reduce_op;
    amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    const int ncomp = mf.nComp();
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const auto arr = mf.const_array(mfi);
        reduce_op.eval(mfi.fabbox(), ncomp, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) -> ReduceTuple
            {
                return std::abs(arr(i,j,k,n) - placeholder);
            });
    }
    amrex::Real change = amrex::get<0>(reduce_data.value(reduce_op));
    amrex::ParallelDescriptor::ReduceRealMax(change);
    return change;
}

} // namespace

// The fills below only copy values that are small integers or half-integers, which are
// exact in single and double precision, so every comparison is exact.

TEST(PlanarBoundaryParallel, SurfaceCopiesAreTheLowestBoxOfEachColumn)
{
    // Four columns of three boxes
    const PlanarLayout layout = make_layout(amrex::IntVect(8,8,4));
    ASSERT_EQ(static_cast<int>(layout.ba3d.size()), 12);

    PlanarBoundary bndry;
    bndry.define(layout.ba3d, layout.ba2d, layout.dm, 0);

    const amrex::BoxArray& ba_sfc = bndry.surface_boxes();
    const amrex::Vector<int>& src = bndry.surface_index();
    ASSERT_EQ(static_cast<int>(ba_sfc.size()), 4);
    ASSERT_EQ(static_cast<int>(src.size()), 4);
    EXPECT_TRUE(ba_sfc.isDisjoint());

    for (int is = 0; is < 4; ++is) {
        SCOPED_TRACE("surface copy " + std::to_string(is));
        EXPECT_EQ(layout.ba3d[src[is]].smallEnd(2), 0);
        EXPECT_EQ(ba_sfc[is], layout.ba2d[src[is]]);
    }
}

// On a BoxArray split in z, every copy of a planar field, valid region and ghost cells,
// holds the surface copy's value after the fill, for cell- and face-centered fields and
// for several components.  A second fill with new surface values reuses the buffer and
// must hand back the new values.
TEST(PlanarBoundaryParallel, FillsEveryCopyFromTheSurfaceCopy)
{
    const PlanarLayout layout = make_layout(amrex::IntVect(8,8,4));
    const amrex::Periodicity period(amrex::IntVect(nx,ny,0));

    PlanarBoundary bndry;
    bndry.define(layout.ba3d, layout.ba2d, layout.dm, 0);

    const amrex::IntVect ng(2,2,0);
    for (const amrex::IntVect& ixtype : {amrex::IntVect(0,0,0), amrex::IntVect(1,0,0), amrex::IntVect(0,1,0)}) {
        for (const int ncomp : {1, 3}) {
            SCOPED_TRACE("index type (" + std::to_string(ixtype[0]) + "," + std::to_string(ixtype[1]) +
                         "), ncomp " + std::to_string(ncomp));

            amrex::MultiFab mf(amrex::convert(layout.ba2d, ixtype), layout.dm, ncomp, ng);
            mf.setVal(placeholder);

            set_surface_copies(mf, layout, amrex::Real(0.0));
            bndry.fill(mf, period);
            EXPECT_EQ(max_error(mf, amrex::Real(0.0)), amrex::Real(0.0));

            set_surface_copies(mf, layout, amrex::Real(0.5));
            bndry.fill(mf, period);
            EXPECT_EQ(max_error(mf, amrex::Real(0.5)), amrex::Real(0.0));
        }
    }
}

// Without a split in z every planar box is a surface copy and the fill is a FillBoundary
TEST(PlanarBoundaryParallel, WithoutSplitFillsLikeFillBoundary)
{
    const PlanarLayout layout = make_layout(amrex::IntVect(8,8,nz));
    ASSERT_EQ(static_cast<int>(layout.ba3d.size()), 4);
    const amrex::Periodicity period(amrex::IntVect(nx,ny,0));

    PlanarBoundary bndry;
    bndry.define(layout.ba3d, layout.ba2d, layout.dm, 0);
    EXPECT_EQ(static_cast<int>(bndry.surface_boxes().size()), 4);

    amrex::MultiFab mf(layout.ba2d, layout.dm, 2, amrex::IntVect(2,2,0));
    mf.setVal(placeholder);
    set_surface_copies(mf, layout, amrex::Real(0.0));
    bndry.fill(mf, period);
    EXPECT_EQ(max_error(mf, amrex::Real(0.0)), amrex::Real(0.0));
}

// A level none of whose boxes reaches the surface has no computed copy, and the fill
// leaves the field as it is
TEST(PlanarBoundaryParallel, LevelAboveTheSurfaceIsLeftAsItIs)
{
    const PlanarLayout layout = make_layout(amrex::IntVect(8,8,4), 4);
    ASSERT_EQ(static_cast<int>(layout.ba3d.size()), 8);
    const amrex::Periodicity period(amrex::IntVect(nx,ny,0));

    PlanarBoundary bndry;
    bndry.define(layout.ba3d, layout.ba2d, layout.dm, 0);
    EXPECT_EQ(static_cast<int>(bndry.surface_boxes().size()), 0);

    amrex::MultiFab mf(layout.ba2d, layout.dm, 1, amrex::IntVect(2,2,0));
    mf.setVal(placeholder);
    bndry.fill(mf, period);
    EXPECT_EQ(max_change_from_placeholder(mf), amrex::Real(0.0));
}
