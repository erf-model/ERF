#include <cmath>
#include <limits>
#include <vector>

#include <AMReX_BCRec.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PhysBCFunct.H>

#include <gtest/gtest.h>

#include "ERF_FillPatcher.H"

// Motivation (erf-model/ERF#4920): on a lateral coarse/fine face the fine faces are filled from the
// coarse face with a linear slope that does not know the wall. The coarse value next to the wall is a point value at
// the centroid of its fluid part, so the fine face at the wall got too much flow (+15 % over flat
// ground, +9 % over immersed terrain at a refined patch's inflow edge). erf.cf_loglaw_fill shares
// the flux of the first coarse face above the wall among its fine sub-faces by ln(1 + z/z0) at
// their heights above the wall. It must keep the flux of every coarse face (the c/f faces are
// averaged down, and moving flux between coarse faces drives a circulation beside the fine box),
// close the solid sub-faces of immersed terrain, and leave the faces above the wall face alone.
// The solid-face closure alone (erf.if_cf_fill = fluid) must close every solid sub-face and move
// its flux to the fluid sub-faces of the same coarse face, or up the column when it has none.

using namespace amrex;

namespace {

constexpr int ncx = 8, ncy = 2, ncz = 8;   // coarse cells, 20 m
constexpr Real dxc = 20.;
constexpr int fi_lo = 4, fi_hi = 11;       // fine patch in x (coarse cells 2-5)
constexpr int fk_hi = 9;                   // fine patch from the ground to 100 m
constexpr Real z0 = 0.1;
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

// coarse x-momentum: a function of height only (u_of_k), ghost cells included
void fill_coarse (MultiFab& mf, std::vector<Real> const& u_of_k)
{
    Gpu::DeviceVector<Real> du(u_of_k.size());
    Gpu::copy(Gpu::hostToDevice, u_of_k.begin(), u_of_k.end(), du.begin());
    Real const* p = du.data();
    const int nk = static_cast<int>(u_of_k.size());
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        auto const& a = mf.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // u_of_k[0] is the ghost row below the ground
            const int kk = amrex::min(amrex::max(k + 1, 0), nk - 1);
            a(i,j,k) = p[kk];
        });
    }
    Gpu::streamSynchronize();
}

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
}

// fill the fine c/f x-faces; returns the fine x-momentum copied to host memory
MultiFab fill_fine (Levels const& L, MultiFab const& crse, MultiFab const* tb, bool loglaw)
{
    BoxArray fbax = convert(L.fba, IntVect(1,0,0));
    ERFFillPatcher fp(fbax, L.fdm, L.fgeom, convert(L.cba, IntVect(1,0,0)), L.cdm, L.cgeom,
                      0, 0, 1, &face_cons_linear_interp);
    fp.SetSolidFraction(tb);
    fp.SetLogLawFill(loglaw, z0, tb);
    fp.RegisterCoarseData({&crse, &crse}, {0., 1.});

    MultiFab fine(fbax, L.fdm, 1, 0);
    fine.setVal(Real(-999.));
    PhysBCFunctNoOp null_bc;
    Vector<BCRec> bcs(1);
    fp.FillSet(fine, 0., null_bc, bcs);

    MultiFab host(fbax, L.fdm, 1, 0, MFInfo().SetArena(The_Pinned_Arena()));
    MultiFab::Copy(host, fine, 0, 0, 1, 0);
    Gpu::streamSynchronize();
    return host;
}

// fine x-face value on the inflow c/f face (i = fi_lo)
Real at (MultiFab const& host, int j, int k)
{
    for (MFIter mfi(host); mfi.isValid(); ++mfi) {
        if (mfi.validbox().contains(IntVect(fi_lo, j, k))) { return host.const_array(mfi)(fi_lo, j, k); }
    }
    return std::numeric_limits<Real>::quiet_NaN();
}

Real ln1p (Real z) { return std::log1p(z / z0); }

} // namespace

TEST(CFLogLawFill, FlatSurfaceKeepsTheCoarseFluxAndSharesItByTheLogLaw)
{
    Levels L;
    MultiFab crse(convert(L.cba, IntVect(1,0,0)), L.cdm, 1, 1);
    // ghost row (zero gradient, as the surface-layer bottom fills it), then the coarse cells
    fill_coarse(crse, {5.674, 5.674, 7.043, 7.8, 8.3, 8.7, 9.0, 9.2, 9.4, 9.5});

    MultiFab lin = fill_fine(L, crse, nullptr, false);
    MultiFab log = fill_fine(L, crse, nullptr, true);

    for (int j = 0; j < 2*ncy; ++j) {
        // the flux of the coarse face at the wall is kept
        EXPECT_NEAR(at(log,j,0) + at(log,j,1), at(lin,j,0) + at(lin,j,1), tol * at(lin,j,1));
        // its two fine faces, 5 and 15 m above the ground, share it by ln(1 + z/z0)
        EXPECT_NEAR(at(log,j,0) / at(log,j,1), ln1p(5.) / ln1p(15.), tol);
        // the linear fill does not (else this test shows nothing)
        EXPECT_GT(std::abs(at(lin,j,0) / at(lin,j,1) - ln1p(5.) / ln1p(15.)), 0.05);
        // the faces above the wall face are left alone
        for (int k = 2; k <= fk_hi; ++k) { EXPECT_EQ(at(log,j,k), at(lin,j,k)); }
    }
}

TEST(CFLogLawFill, ImmersedSurfaceClosesTheSolidAndSharesTheWallFace)
{
    Levels L;
    MultiFab crse(convert(L.cba, IntVect(1,0,0)), L.cdm, 1, 1);
    // surface at 45 m: coarse cells 0-40 m solid, 40-60 m a quarter solid
    fill_coarse(crse, {0., 0., 0., 5.3, 6.93, 7.78, 8.35, 8.7, 9.0, 9.2});
    // fine: 0-40 m solid, the 40-50 m cell half solid
    MultiFab tb(L.fba, L.fdm, 1, 2);
    fill_tblank(tb, 4, Real(0.5));

    MultiFab lin = fill_fine(L, crse, &tb, false);
    MultiFab log = fill_fine(L, crse, &tb, true);

    for (int j = 0; j < 2*ncy; ++j) {
        Real sum_lin = 0., sum_log = 0.;
        for (int k = 0; k <= 5; ++k) { sum_lin += at(lin,j,k); sum_log += at(log,j,k); }
        // the solid sub-faces carry nothing
        for (int k = 0; k < 4; ++k) { EXPECT_EQ(at(log,j,k), Real(0.)); }
        // the flux of the solid coarse faces and of the wall coarse face is kept
        EXPECT_NEAR(sum_log, sum_lin, tol * sum_lin);
        // the fine faces of the wall coarse face, 2.5 and 10 m above the surface (the fluid centroid
        // of the half-solid face, and the face above), share it by ln(1 + z/z0)
        EXPECT_NEAR(at(log,j,4) / at(log,j,5), ln1p(2.5) / ln1p(10.), tol);
        // (the linear fill gives 0.733 here against 0.706)
        EXPECT_GT(std::abs(at(lin,j,4) / at(lin,j,5) - ln1p(2.5) / ln1p(10.)), 0.01);
        for (int k = 6; k <= fk_hi; ++k) { EXPECT_EQ(at(log,j,k), at(lin,j,k)); }
    }
}

TEST(CFLogLawFill, SolidClosureMovesTheFluxWithinTheCoarseFaceAndUpTheColumn)
{
    Levels L;
    MultiFab crse(convert(L.cba, IntVect(1,0,0)), L.cdm, 1, 1);
    // the coarse faces below 60 m carry flow too (the coarse projection puts flux through solid
    // faces, which the immersed forcing only removes afterwards)
    fill_coarse(crse, {1., 1., 2., 3., 6., 7., 8., 8.5, 9., 9.2});
    // surface at 50 m, the mid-plane of the coarse cell 40-60 m: fine cells 0-50 m solid
    MultiFab tb(L.fba, L.fdm, 1, 2);
    fill_tblank(tb, 5, Real(0.));

    MultiFab lin = fill_fine(L, crse, nullptr, false);
    MultiFab cls = fill_fine(L, crse, &tb, false);

    for (int j = 0; j < 2*ncy; ++j) {
        // the solid sub-faces (0-50 m) carry nothing
        for (int k = 0; k < 5; ++k) { EXPECT_EQ(at(cls,j,k), Real(0.)); }
        // the coarse faces 0-40 m have no fluid sub-face: their flux is carried up to the coarse
        // face 40-60 m, whose one fluid sub-face (50-60 m) takes it all with its own
        Real carried = 0.;
        for (int k = 0; k < 6; ++k) { carried += at(lin,j,k); }
        EXPECT_NEAR(at(cls,j,5), carried, tol * carried);
        EXPECT_GT(at(cls,j,5), at(lin,j,5));   // (else the closure moved nothing)
        // the faces above are left alone
        for (int k = 6; k <= fk_hi; ++k) { EXPECT_EQ(at(cls,j,k), at(lin,j,k)); }
    }
}
