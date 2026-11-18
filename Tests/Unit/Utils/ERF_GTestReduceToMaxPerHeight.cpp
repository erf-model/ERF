#include <memory>

#include <AMReX_BoxArray.H>
#include <AMReX_BoxList.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_MultiFab.H>

#include <gtest/gtest.h>

#include "ERF_ParFunctions.H"

// Motivation: the input sponge takes, for every k of a level, the largest cell-centre height at
// that k (reduce_to_max_per_height over z_phys_cc).  The reduction read every k from every box,
// so a box that does not span the column -- a refined patch near the ground, or grids split in
// z -- was read outside its data: garbage heights, or a crash.  Each box must contribute only
// the k it holds, and a k that no box holds must come back as the nominal height of the level,
// so that no caller is handed a value that is not a height.

using namespace amrex;

namespace {

constexpr int nx = 16;
constexpr int ny = 8;
constexpr int nz = 8;

// A value that depends on the index only, with a known maximum per k
Real
cell_value (int i, int j, int k)
{
    return Real(100.)*k + Real(2.)*i + j;
}

// Nominal face heights of the column, unlike any cell value
Vector<Real>
nominal_faces ()
{
    Vector<Real> zs(nz+1);
    for (int k = 0; k <= nz; ++k) { zs[k] = Real(-10.) - Real(4.)*k*k; }
    return zs;
}

std::unique_ptr<MultiFab>
make_field (const BoxArray& ba)
{
    auto mf = std::make_unique<MultiFab>(ba, DistributionMapping(ba), 1, 2);
    // ghost cells hold a value larger than any valid one: the reduction must not read them
    mf->setVal(Real(1.0e6));
    for (MFIter mfi(*mf); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        auto a = mf->array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a(i,j,k) = Real(100.)*k + Real(2.)*i + j;
        });
    }
    return mf;
}

} // namespace

// Boxes that stop at k = 3 of an 8-cell column, the left one shorter: k = 0..3 get the exact
// maximum over the boxes that hold them, k = 4..7 are not covered and get the nominal height
// of the cell centre.
TEST(ReduceToMaxPerHeight, BoxesShorterThanTheColumn)
{
    BoxList bl;
    bl.push_back(Box(IntVect(0, 0, 0), IntVect(7, ny-1, 3)));
    bl.push_back(Box(IntVect(8, 0, 0), IntVect(nx-1, ny-1, 2)));
    const BoxArray ba(bl);
    auto mf = make_field(ba);

    const Vector<Real> zs = nominal_faces();
    Vector<Real> v(nz, Real(-7.0));
    reduce_to_max_per_height(v, mf, zs);

    for (int k = 0; k < nz; ++k) {
        if (k <= 2) {
            EXPECT_EQ(v[k], cell_value(nx-1, ny-1, k)) << "k = " << k;
        } else if (k == 3) {
            EXPECT_EQ(v[k], cell_value(7, ny-1, k)) << "k = " << k;   // only the left box holds k = 3
        } else {
            EXPECT_EQ(v[k], Real(0.5)*(zs[k] + zs[k+1])) << "k = " << k;
        }
    }
}

// A patch that reaches neither the ground nor the top, on a field that is nodal in z (the
// subsidence profile reads z_phys_nd): the faces k = 2..5 get the exact maximum, the faces
// below and above get the nominal face height, k = 0 included (a profile may anchor on it).
TEST(ReduceToMaxPerHeight, NodalPatchOffTheGround)
{
    const BoxArray ba_cc(Box(IntVect(4, 0, 2), IntVect(11, ny-1, 4)));
    const BoxArray ba = amrex::convert(ba_cc, IntVect(0, 0, 1));
    auto mf = make_field(ba);

    const Vector<Real> zs = nominal_faces();
    Vector<Real> v(nz+1, Real(-7.0));
    reduce_to_max_per_height(v, mf, zs);

    for (int k = 0; k <= nz; ++k) {
        if (k >= 2 && k <= 5) {
            EXPECT_EQ(v[k], cell_value(11, ny-1, k)) << "k = " << k;
        } else {
            EXPECT_EQ(v[k], zs[k]) << "k = " << k;
        }
    }
}

// Boxes stacked in z (grids split in z): every k is held by exactly one box.
TEST(ReduceToMaxPerHeight, BoxesStackedInZ)
{
    BoxList bl;
    bl.push_back(Box(IntVect(0, 0, 0), IntVect(nx-1, ny-1, 2)));
    bl.push_back(Box(IntVect(0, 0, 3), IntVect(nx-1, ny-1, nz-1)));
    const BoxArray ba(bl);
    auto mf = make_field(ba);

    Vector<Real> v(nz, Real(-7.0));
    reduce_to_max_per_height(v, mf, nominal_faces());

    for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(v[k], cell_value(nx-1, ny-1, k)) << "k = " << k;
    }
}
