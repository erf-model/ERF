#include <cmath>
#include <limits>

#include <AMReX_BCRec.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Geometry.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFab.H>

#include <gtest/gtest.h>

#include "ERF_Constants.H"
#include "ERF_LandWaterInterp.H"

// Motivation (erf-model/ERF#4390): a finer level that takes its land state from a coarser one
// interpolated the land fields across the coastline, so the slope of a coastal land cell was
// set partly by its water neighbour (soil wetness 1 over water in SLM) and the fine cells inside
// it got a spurious land/water gradient. LandWaterConsLinear interpolates land from land and
// water from water only, as WRF's masked nest interpolation does. It must keep each coarse
// cell's mean, leave every cell whose neighbours are all of its own type exactly as
// cell_cons_interp has it, and give no slope toward a neighbour outside the domain.

using namespace amrex;

namespace {

constexpr int nc = 8;        // coarse cells in x and y (one in z)
constexpr int ncoast = 4;    // coarse cells i < ncoast are land, the rest water
const IntVect ratio(2, 2, 1);

struct Levels {
    Geometry cgeom, fgeom;
    BoxArray cba, fba;
    DistributionMapping cdm, fdm;
    explicit Levels (bool periodic_x) {
        RealBox rb({0., 0., 0.}, {Real(nc), Real(nc), 1.});
        Array<int,3> per{int(periodic_x), 1, 0};
        Box cdom(IntVect(0), IntVect(nc-1, nc-1, 0));
        cgeom.define(cdom, rb, 0, per);
        fgeom.define(amrex::refine(cdom, ratio), rb, 0, per);
        cba = BoxArray(cdom);
        fba = BoxArray(amrex::refine(cdom, ratio));
        cdm = DistributionMapping(cba);
        fdm = DistributionMapping(fba);
    }
};

Real land_value (int i, int j) { return Real(0.3) + Real(0.05)*i + Real(0.02)*j; }

// Coarse field with its land/water indicator (1 water, 0 land), one ghost cell, both filled
// across periodic boundaries; the indicator is NaN outside a non-periodic domain
MultiFab make_coarse (const Levels& L, Real water_value, Real land_offset)
{
    MultiFab host(L.cba, L.cdm, 2, 1, MFInfo().SetArena(The_Pinned_Arena()));
    for (MFIter mfi(host); mfi.isValid(); ++mfi) {
        const auto a = host.array(mfi);
        const Box b = mfi.fabbox();
        for (int j = b.smallEnd(1); j <= b.bigEnd(1); ++j) {
        for (int i = b.smallEnd(0); i <= b.bigEnd(0); ++i) {
            const int iw = ((i % nc) + nc) % nc;   // the periodic image
            const bool water = iw >= ncoast;
            a(i,j,0,0) = water ? water_value : land_value(iw, ((j % nc) + nc) % nc) + land_offset;
            a(i,j,0,1) = water ? Real(1.) : Real(0.);
        }}
    }
    host.setDomainBndry(std::numeric_limits<Real>::quiet_NaN(), 1, 1, L.cgeom);
    MultiFab crse(L.cba, L.cdm, 2, 1);
    MultiFab::Copy(crse, host, 0, 0, 2, 1);
    return crse;
}

// Fine field interpolated from crse (with interp; a single component unless masked), on host
MultiFab interpolate (const Levels& L, const MultiFab& crse, Interpolater* interp, int ncomp)
{
    MultiFab fine(L.fba, L.fdm, ncomp, 0);
    fine.setVal(Real(0.));
    Vector<BCRec> bcs(2);
    for (auto& bc : bcs) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            const int t = L.cgeom.isPeriodic(d) ? BCType::int_dir : BCType::foextrap;
            bc.setLo(d, t); bc.setHi(d, t);
        }
    }
    InterpFromCoarseLevel(fine, IntVect(0), IntVect(0), crse, 0, 0, ncomp,
                          L.cgeom, L.fgeom, ratio, interp, bcs, 0);
    MultiFab host(L.fba, L.fdm, ncomp, 0, MFInfo().SetArena(The_Pinned_Arena()));
    MultiFab::Copy(host, fine, 0, 0, ncomp, 0);
    Gpu::streamSynchronize();
    return host;
}

Real at (const MultiFab& host, int i, int j, int n = 0)
{
    for (MFIter mfi(host); mfi.isValid(); ++mfi) {
        if (mfi.validbox().contains(IntVect(i, j, 0))) { return host.const_array(mfi)(i, j, 0, n); }
    }
    return std::numeric_limits<Real>::quiet_NaN();
}

bool coarse_is_water (int ic) { return ((ic % nc) + nc) % nc >= ncoast; }

// As NOAHMP::interp_from_lev0 does: the field alone (water holding lsm_undefined), its
// indicator the validity of each value, interpolated with interp_land_water; on host
MultiFab interpolate_by_validity (const Levels& L, Real land_offset)
{
    MultiFab src(L.cba, L.cdm, 1, 1);
    MultiFab::Copy(src, make_coarse(L, lsm_undefined, land_offset), 0, 0, 1, 1);
    MultiFab valid(L.cba, L.cdm, 1, 1);
    for (MFIter mfi(src); mfi.isValid(); ++mfi) {
        const auto v = src.const_array(mfi);
        const auto w = valid.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            w(i,j,k) = is_valid_lsm_value(v(i,j,k)) ? Real(0.0) : Real(1.0);
        });
    }
    MultiFab fine(L.fba, L.fdm, 1, 0);
    fine.setVal(Real(0.));
    Vector<BCRec> bcs(2);
    for (auto& bc : bcs) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            const int t = L.cgeom.isPeriodic(d) ? BCType::int_dir : BCType::foextrap;
            bc.setLo(d, t); bc.setHi(d, t);
        }
    }
    interp_land_water(fine, src, valid, L.cgeom, L.fgeom, ratio, bcs, 0);
    MultiFab host(L.fba, L.fdm, 1, 0, MFInfo().SetArena(The_Pinned_Arena()));
    MultiFab::Copy(host, fine, 0, 0, 1, 0);
    Gpu::streamSynchronize();
    return host;
}

} // namespace

TEST(LandWaterInterp, LandAndWaterDoNotSeeEachOther)
{
    Levels L(true);
    LandWaterConsLinear interp;
    const MultiFab a = interpolate(L, make_coarse(L, Real(1.0), Real(0.)), &interp, 2);
    const MultiFab b = interpolate(L, make_coarse(L, Real(5.0), Real(0.)), &interp, 2);   // other water
    const MultiFab c = interpolate(L, make_coarse(L, Real(1.0), Real(0.3)), &interp, 2);  // other land
    for (int j = 0; j < nc*ratio[1]; ++j) {
        for (int i = 0; i < nc*ratio[0]; ++i) {
            if (coarse_is_water(i / ratio[0])) {
                EXPECT_EQ(at(a, i, j), at(c, i, j)) << "water at i=" << i << " j=" << j;
                EXPECT_EQ(at(a, i, j, 1), Real(1.));
            } else {
                EXPECT_EQ(at(a, i, j), at(b, i, j)) << "land at i=" << i << " j=" << j;
                EXPECT_EQ(at(a, i, j, 1), Real(0.));
            }
        }
    }
}

TEST(LandWaterInterp, KeepsCoarseMeansAndMatchesCellConsAwayFromTheCoast)
{
    Levels L(true);
    LandWaterConsLinear interp;
    const MultiFab crse = make_coarse(L, Real(1.0), Real(0.));
    const MultiFab masked = interpolate(L, crse, &interp, 2);
    const MultiFab plain  = interpolate(L, crse, &cell_cons_interp, 1);
    for (int jc = 0; jc < nc; ++jc) {
        for (int ic = 0; ic < nc; ++ic) {
            Real sum = Real(0.);
            for (int jf = 0; jf < ratio[1]; ++jf) {
            for (int iff = 0; iff < ratio[0]; ++iff) {
                const int i = ic*ratio[0] + iff, j = jc*ratio[1] + jf;
                sum += at(masked, i, j);
                // all x neighbours of the same type (y is uniform in type): exactly cell_cons_interp
                if (coarse_is_water(ic-1) == coarse_is_water(ic) &&
                    coarse_is_water(ic+1) == coarse_is_water(ic)) {
                    EXPECT_EQ(at(masked, i, j), at(plain, i, j)) << "i=" << i << " j=" << j;
                }
            }}
            const Real parent = coarse_is_water(ic) ? Real(1.0) : land_value(ic, jc);
            EXPECT_NEAR(sum / Real(ratio[0]*ratio[1]), parent, Real(100.)*std::numeric_limits<Real>::epsilon()) << "ic=" << ic << " jc=" << jc;
        }
    }
    // The coastal land cell (ic = ncoast-1) has no x slope, and keeps its y slope
    const int ic = ncoast - 1, jc = 3;
    EXPECT_EQ(at(masked, ic*2, jc*2), at(masked, ic*2+1, jc*2));
    EXPECT_GT(at(masked, ic*2, jc*2+1), at(masked, ic*2, jc*2));
    // ... where cell_cons_interp had the water value steepen its x slope
    EXPECT_NE(at(plain, ic*2, jc*2), at(plain, ic*2+1, jc*2));
}

TEST(LandWaterInterp, NoSlopeTowardOutsideTheDomain)
{
    Levels L(false);   // not periodic in x: the indicator is NaN beyond i = 0 and i = nc-1
    LandWaterConsLinear interp;
    const MultiFab f = interpolate(L, make_coarse(L, Real(1.0), Real(0.)), &interp, 2);
    for (int j = 0; j < nc*ratio[1]; ++j) {
        EXPECT_EQ(at(f, 0, j), at(f, 1, j)) << "j=" << j;   // land cell at the domain edge
        EXPECT_TRUE(std::isfinite(at(f, 0, j)));
    }
    // an interior land cell still has its x slope
    EXPECT_LT(at(f, 2, 4), at(f, 3, 4));
}

// Motivation: Noah-MP holds the lsm_undefined sentinel over water and its fine levels take
// their land state from level 0 every step. The sentinel must come through exactly (consumers
// reject it with is_valid_lsm_value) and must not enter a land cell's slope.
TEST(LandWaterInterp, SentinelComesThroughExactlyAndStaysOutOfLandSlopes)
{
    Levels L(true);
    const MultiFab f = interpolate_by_validity(L, Real(0.));
    for (int j = 0; j < nc*ratio[1]; ++j) {
        for (int i = 0; i < nc*ratio[0]; ++i) {
            if (coarse_is_water(i / ratio[0])) {
                EXPECT_EQ(at(f, i, j), lsm_undefined) << "i=" << i << " j=" << j;
            } else {
                EXPECT_TRUE(is_valid_lsm_value(at(f, i, j))) << "i=" << i << " j=" << j;
                EXPECT_LT(at(f, i, j), Real(1.)) << "i=" << i << " j=" << j;
            }
        }
    }
    // the coastal land cell has no slope toward the sentinel
    const int ic = ncoast - 1, jc = 3;
    EXPECT_EQ(at(f, ic*2, jc*2), at(f, ic*2+1, jc*2));
}
