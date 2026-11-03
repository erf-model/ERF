#include <cmath>
#include <limits>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFab.H>

#include <gtest/gtest.h>

#include "ERF_FillPatcher.H"
#include "ERF_Utils.H"

// Motivation (erf-model/ERF#4920): two multi-level immersed-forcing checks at the coarse/fine faces.
// (1) When the surface lies on a coarse cell's mid-plane, the coarse cell is half solid over fine
// cells that are all solid and all fluid, and the two levels disagree on the surface (the case
// kynema-sgf also excluded); the start-up warning must count those coarse cells along the c/f faces
// and only there. (2) The coarse level zeroes its fully solid faces after its projection, so the
// fluxes filled into a fine patch's c/f faces need not sum to zero; the singular fine projection
// then spread the imbalance as a uniform divergence, a spurious heat source. The balance must make
// the net c/f flux zero, taking it off the fluid c/f faces only, and leave the solid faces and the
// faces on a non-periodic side of the domain alone.

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

// x-momentum: inflow value on the face i = fi_lo, outflow value on i = fi_hi+1, 1 on the solid
// sub-faces (k < 4) of both, 0 elsewhere; z-momentum w_top on the top faces; y-momentum 0
void fill_mom (MultiFab& xm, MultiFab& ym, MultiFab& zm, Real u_in, Real u_out, Real w_top)
{
    for (MFIter mfi(xm); mfi.isValid(); ++mfi) {
        auto const& a = xm.array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const bool edge = (i == fi_lo || i == fi_hi + 1);
            a(i,j,k) = !edge ? Real(0.) : ((k < 4) ? Real(1.) : ((i == fi_lo) ? u_in : u_out));
        });
    }
    ym.setVal(Real(0.));
    for (MFIter mfi(zm); mfi.isValid(); ++mfi) {
        auto const& a = zm.array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a(i,j,k) = (k == fk_hi + 1) ? w_top : Real(0.);
        });
    }
    Gpu::streamSynchronize();
}

// host copy
MultiFab on_host (MultiFab const& mf)
{
    MultiFab h(mf.boxArray(), mf.DistributionMap(), 1, 0, MFInfo().SetArena(The_Pinned_Arena()));
    MultiFab::Copy(h, mf, 0, 0, 1, 0);
    Gpu::streamSynchronize();
    return h;
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

TEST(ImmersedCFFaces, BalanceZeroesTheNetFluxOnTheFluidFaces)
{
    Levels L;
    // the fill masks of the three face directions, as ERF builds them for two-way coupling
    ERFFillPatcher fx(convert(L.fba, IntVect(1,0,0)), L.fdm, L.fgeom, convert(L.cba, IntVect(1,0,0)), L.cdm, L.cgeom,
                      0, 0, 1, &face_cons_linear_interp);
    ERFFillPatcher fy(convert(L.fba, IntVect(0,1,0)), L.fdm, L.fgeom, convert(L.cba, IntVect(0,1,0)), L.cdm, L.cgeom,
                      0, 0, 1, &face_cons_linear_interp);
    ERFFillPatcher fz(convert(L.fba, IntVect(0,0,1)), L.fdm, L.fgeom, convert(L.cba, IntVect(0,0,1)), L.cdm, L.cgeom,
                      0, 0, 1, &face_cons_linear_interp);

    MultiFab tb(L.fba, L.fdm, 1, 1);
    fill_tblank(tb, 4, Real(0.5));                 // surface at 45 m
    MultiFab xm(convert(L.fba, IntVect(1,0,0)), L.fdm, 1, 0);
    MultiFab ym(convert(L.fba, IntVect(0,1,0)), L.fdm, 1, 0);
    MultiFab zm(convert(L.fba, IntVect(0,0,1)), L.fdm, 1, 0);
    const Real u_in = 5., u_out = 4., w_top = 0.2;
    fill_mom(xm, ym, zm, u_in, u_out, w_top);

    Vector<BoxArray> subdomains{L.fba};
    if_balance_cf_fluxes(L.fgeom, subdomains, L.fba, {&xm, &ym, &zm},
                         {fx.GetMask(), fy.GetMask(), fz.GetMask()},
                         {fx.GetSetMaskVal(), fy.GetSetMaskVal(), fz.GetSetMaskVal()}, tb);

    const MultiFab hx = on_host(xm), hy = on_host(ym), hz = on_host(zm);
    const Real ax = (dxc/2) * (dxc/2), az = (dxc/2) * (dxc/2);
    Real net = 0.;
    for (int j = 0; j < 2*ncy; ++j) {
        for (int k = 0; k <= fk_hi; ++k) {
            const Real in  = hx[0](IntVect(fi_lo,   j, k));
            const Real out = hx[0](IntVect(fi_hi+1, j, k));
            net += (in - out) * ax;
            if (k < 4) {   // solid sub-faces are left alone
                EXPECT_EQ(in,  Real(1.));
                EXPECT_EQ(out, Real(1.));
            }
            // interior faces are not c/f faces
            EXPECT_EQ(hx[0](IntVect(fi_lo+3, j, k)), Real(0.));
        }
    }
    for (int j = 0; j < 2*ncy; ++j) {
        for (int i = fi_lo; i <= fi_hi; ++i) {
            net -= hz[0](IntVect(i, j, fk_hi+1)) * az;
            // the bottom faces are on the (non-periodic) domain boundary
            EXPECT_EQ(hz[0](IntVect(i, j, 0)), Real(0.));
        }
    }
    // the y faces at the periodic sides are fine/fine faces (the patch spans y): not c/f, untouched
    for (int k = 0; k <= fk_hi; ++k) { EXPECT_EQ(hy[0](IntVect(fi_lo, 0, k)), Real(0.)); }

    const Real net_before = (2*ncy) * 6 * (u_in - u_out) * ax - (2*ncy) * (fi_hi - fi_lo + 1) * w_top * az;
    ASSERT_GT(std::abs(net_before), 100.);
    EXPECT_NEAR(net, 0., tol * std::abs(net_before));
}
