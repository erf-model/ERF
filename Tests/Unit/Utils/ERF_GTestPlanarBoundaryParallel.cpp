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

// Give every surface copy its own constant, so that neighbouring surface boxes disagree at
// any face they share.  Component n is offset by n.  Returns the constants, indexed as
// surface_boxes(), so that a caller can say which box a filled value came from.
amrex::Vector<amrex::Real> set_distinct_surface_copies (amrex::MultiFab& mf,
                                                        const PlanarBoundary& bndry)
{
    const amrex::Vector<int>& src = bndry.surface_index();
    const int nsfc = static_cast<int>(src.size());

    amrex::Vector<amrex::Real> vals(nsfc);
    for (int is = 0; is < nsfc; ++is) { vals[is] = amrex::Real(is + 1); }

    const int ncomp = mf.nComp();
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        int is = -1;
        for (int s = 0; s < nsfc; ++s) { if (src[s] == mfi.index()) { is = s; break; } }
        if (is < 0) { continue; }  // not a surface copy, so the surface layer never computes it
        const auto arr = mf.array(mfi);
        const amrex::Real v = vals[is];
        amrex::ParallelFor(mfi.validbox(), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            arr(i,j,k,n) = v + amrex::Real(n);
        });
    }
    return vals;
}

// The surface box a point belongs to: the lowest-index box whose converted box contains it,
// counting the lateral periodic images.  At a face shared by two surface boxes that is the
// one the fill must hand back, in every copy -- it is the precedence OverrideSync gives, and
// it depends on neither the decomposition nor the rank count.
int owning_surface_box (const amrex::BoxArray& ba_sfc, const amrex::IntVect& ixtype, int i, int j)
{
    for (int is = 0; is < static_cast<int>(ba_sfc.size()); ++is) {
        const amrex::Box b = amrex::convert(ba_sfc[is], ixtype);
        for (int si = -1; si <= 1; ++si) {
            for (int sj = -1; sj <= 1; ++sj) {
                if (b.contains(amrex::IntVect(i + si*nx, j + sj*ny, 0))) { return is; }
            }
        }
    }
    return -1;
}

// Largest |value - the owning surface box's constant| over the valid region and ghost cells
// of every copy on every rank.  Non-zero means some copy disagrees with the owner, which is
// what a shared face supplied by the wrong neighbour looks like.
amrex::Real max_error_against_owner (const amrex::MultiFab& mf,
                                     const PlanarBoundary& bndry,
                                     const amrex::Vector<amrex::Real>& vals,
                                     const amrex::IntVect& ixtype)
{
    // Tabulate the owner's constant over every (i,j) a box or its ghost cells can reach
    const amrex::BoxArray& ba_sfc = bndry.surface_boxes();
    const int pad  = 4;
    const int ilo  = -pad,  jlo = -pad;
    const int ni   = nx + 1 + 2*pad;
    const int nj   = ny + 1 + 2*pad;

    amrex::Vector<amrex::Real> h_val(static_cast<std::size_t>(ni)*nj);
    for (int j = 0; j < nj; ++j) {
        for (int i = 0; i < ni; ++i) {
            const int is = owning_surface_box(ba_sfc, ixtype, ilo+i, jlo+j);
            EXPECT_GE(is, 0) << "no surface box owns (" << ilo+i << "," << jlo+j << ")";
            h_val[static_cast<std::size_t>(i) + static_cast<std::size_t>(ni)*j] =
                (is >= 0) ? vals[is] : amrex::Real(0.0);
        }
    }
    amrex::Gpu::DeviceVector<amrex::Real> d_val(h_val.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_val.begin(), h_val.end(), d_val.begin());
    amrex::Gpu::streamSynchronize();
    const amrex::Real* owner = d_val.data();

    amrex::ReduceOps<amrex::ReduceOpMax> reduce_op;
    amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    const int ncomp = mf.nComp();
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const auto arr = mf.const_array(mfi);
        reduce_op.eval(mfi.fabbox(), ncomp, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) -> ReduceTuple
            {
                const amrex::Real e = owner[(i-ilo) + ni*(j-jlo)] + amrex::Real(n);
                return std::abs(arr(i,j,k,n) - e);
            });
    }
    amrex::Real err = amrex::get<0>(reduce_data.value(reduce_op));
    amrex::ParallelDescriptor::ReduceRealMax(err);
    return err;
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

// A face-centered planar field has faces that two surface boxes both hold -- the MOST
// velocity averages are the real case.  Both boxes compute the same value there today, so a
// fill that takes such a face from either one looks right; this test breaks that tie by
// giving every surface copy a different constant, and then requires the filled field to hand
// back the owning box's constant in *every* copy, valid region and ghost cells alike.
//
// That pins down two things a uniform value cannot see: that all the copies of a shared face
// agree with each other, and that which box supplies it is fixed rather than left to the copy
// order.  The test runs on 1 and on 2 ranks, so a winner that changed with the rank count
// would fail here.
TEST(PlanarBoundaryParallel, SharedFacesComeFromOneSurfaceBox)
{
    const PlanarLayout layout = make_layout(amrex::IntVect(8,8,4));
    const amrex::Periodicity period(amrex::IntVect(nx,ny,0));

    PlanarBoundary bndry;
    bndry.define(layout.ba3d, layout.ba2d, layout.dm, 0);
    ASSERT_EQ(static_cast<int>(bndry.surface_boxes().size()), 4);

    const amrex::IntVect ng(2,2,0);
    for (const amrex::IntVect& ixtype : {amrex::IntVect(0,0,0), amrex::IntVect(1,0,0), amrex::IntVect(0,1,0)}) {
        for (const int ncomp : {1, 3}) {
            SCOPED_TRACE("index type (" + std::to_string(ixtype[0]) + "," + std::to_string(ixtype[1]) +
                         "), ncomp " + std::to_string(ncomp));

            amrex::MultiFab mf(amrex::convert(layout.ba2d, ixtype), layout.dm, ncomp, ng);
            mf.setVal(placeholder);

            const amrex::Vector<amrex::Real> vals = set_distinct_surface_copies(mf, bndry);
            bndry.fill(mf, period);
            EXPECT_EQ(max_error_against_owner(mf, bndry, vals, ixtype), amrex::Real(0.0));
        }
    }
}
