#include <cmath>
#include <limits>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <gtest/gtest.h>

#include "ERF_Utils.H"

// Motivation (erf-model/ERF#4920): when the surface lies on a coarse cell's mid-plane, the coarse
// cell is half solid over fine cells that are all solid and all fluid, and the two levels disagree
// on the surface along the lateral coarse/fine faces (the case kynema-sgf also excluded); the
// start-up warning (erf.if_cf_mismatch_warning) must count those coarse cells along the c/f faces
// and only there.

using namespace amrex;

namespace {

constexpr int ncx = 8, ncy = 2, ncz = 8;   // coarse cells, 20 m
constexpr Real dxc = 20.;
constexpr int fi_lo = 4, fi_hi = 11;       // fine patch in x (coarse cells 2-5)
constexpr int fk_hi = 9;                   // fine patch from the ground to 100 m
const double tol = 1000. * std::numeric_limits<Real>::epsilon();

struct Levels {
    Geometry cgeom, fgeom;
    BoxArray cba, fba;
    DistributionMapping cdm, fdm;
    Levels () {
        RealBox rb({0., 0., 0.}, {ncx*dxc, ncy*dxc, ncz*dxc});
        Array<int,3> per{1, 1, 0};
        Box cdom(IntVect(0), IntVect(ncx-1, ncy-1, ncz-1));
        cgeom.define(cdom, rb, 0, per);
        fgeom.define(amrex::refine(cdom, 2), rb, 0, per);
        cba = BoxArray(cdom);
        fba = BoxArray(Box(IntVect(fi_lo, 0, 0), IntVect(fi_hi, 2*ncy-1, fk_hi)));
        cdm = DistributionMapping(cba);
        fdm = DistributionMapping(fba);
    }
};

// fine solid fraction: solid below the fine cell ksolid, fraction bpart in it, fluid above
void fill_tblank (MultiFab& tb, int ksolid, Real bpart)
{
    for (MFIter mfi(tb); mfi.isValid(); ++mfi) {
        auto const& a = tb.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a(i,j,k) = (k < ksolid) ? Real(1.) : ((k == ksolid) ? bpart : Real(0.));
        });
    }
    Gpu::streamSynchronize();
}

} // namespace

TEST(ImmersedCFFaces, MismatchCountsSurfaceOnTheCoarseMidPlaneAlongCFFacesOnly)
{
    Levels L;
    MultiFab tb(L.fba, L.fdm, 1, 0);

    // surface at 50 m, the mid-plane of the coarse cell 40-60 m: the fine cells under it are all
    // solid (40-50 m) and all fluid (50-60 m). Along the two c/f faces: 2 coarse columns x 2 in y.
    // The covered coarse cells away from the c/f faces (x = 3, 4) disagree too but are not counted.
    fill_tblank(tb, 5, Real(0.));
    auto [n1, worst1] = if_cf_surface_mismatch(tb, L.fba, L.cba, L.cdm, L.cgeom, IntVect(2), Real(0.75));
    EXPECT_EQ(n1, 4);
    EXPECT_NEAR(worst1, 1., tol);

    // surface at 45 m: the fine cell 40-50 m is half solid, the spread 0.5 is below the threshold
    fill_tblank(tb, 4, Real(0.5));
    auto [n2, worst2] = if_cf_surface_mismatch(tb, L.fba, L.cba, L.cdm, L.cgeom, IntVect(2), Real(0.75));
    EXPECT_EQ(n2, 0);
    EXPECT_NEAR(worst2, 0.5, tol);
}
