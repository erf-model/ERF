#include <cmath>
#include <limits>

#include <AMReX_BCRec.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Geometry.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFab.H>

#include <gtest/gtest.h>

#include "ERF_Constants.H"
#include "ERF_TypedInterp.H"

// Motivation (erf-model/ERF#4390): a finer level that takes its land state from a coarser one
// interpolated the land fields across the coastline, so the slope of a coastal land cell was
// set partly by its water neighbour (soil wetness 1 over water in SLM) and the fine cells inside
// it got a spurious land/water gradient. TypedConsLinear interpolates land from land and
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
// indicator the validity of each value, interpolated with interp_by_type; on host
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
    interp_by_type(fine, src, valid, L.cgeom, L.fgeom, ratio, bcs, 0);
    MultiFab host(L.fba, L.fdm, 1, 0, MFInfo().SetArena(The_Pinned_Arena()));
    MultiFab::Copy(host, fine, 0, 0, 1, 0);
    Gpu::streamSynchronize();
    return host;
}

} // namespace

TEST(TypedInterp, LandAndWaterDoNotSeeEachOther)
{
    Levels L(true);
    TypedConsLinear interp;
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

TEST(TypedInterp, KeepsCoarseMeansAndMatchesCellConsAwayFromTheCoast)
{
    Levels L(true);
    TypedConsLinear interp;
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

TEST(TypedInterp, NoSlopeTowardOutsideTheDomain)
{
    Levels L(false);   // not periodic in x: the indicator is NaN beyond i = 0 and i = nc-1
    TypedConsLinear interp;
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
TEST(TypedInterp, SentinelComesThroughExactlyAndStaysOutOfLandSlopes)
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

namespace {

// interp_by_type of src (any number of fields) with indicator, from a BC vector of one entry
MultiFab interpolate_by_type (const Levels& L, const MultiFab& src, const MultiFab& indicator)
{
    MultiFab fine(L.fba, L.fdm, src.nComp(), 0);
    fine.setVal(Real(0.));
    Vector<BCRec> bcs(1);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const int t = L.cgeom.isPeriodic(d) ? BCType::int_dir : BCType::foextrap;
        bcs[0].setLo(d, t); bcs[0].setHi(d, t);
    }
    interp_by_type(fine, src, indicator, L.cgeom, L.fgeom, ratio, bcs, 0);
    MultiFab host(L.fba, L.fdm, src.nComp(), 0, MFInfo().SetArena(The_Pinned_Arena()));
    MultiFab::Copy(host, fine, 0, 0, src.nComp(), 0);
    Gpu::streamSynchronize();
    return host;
}

// The validity of each value of each field of src as its indicator (0 valid, 1 lsm_undefined);
// device lambdas cannot sit in a TEST body, which gtest makes a private member function
MultiFab make_validity (const Levels& L, const MultiFab& src)
{
    const int ncomp = src.nComp();
    MultiFab valid(L.cba, L.cdm, ncomp, 1);
    for (MFIter mfi(src); mfi.isValid(); ++mfi) {
        const auto v = src.const_array(mfi);
        const auto w = valid.array(mfi);
        ParallelFor(mfi.fabbox(), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            w(i,j,k,n) = is_valid_lsm_value(v(i,j,k,n)) ? Real(0.0) : Real(1.0);
        });
    }
    return valid;
}

// The sentinel in column ic (and its periodic images) of crse, typed as land/water type + 2
void put_sentinel_column (MultiFab& crse, int ic)
{
    for (MFIter mfi(crse); mfi.isValid(); ++mfi) {
        const auto a = crse.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (((i % nc) + nc) % nc == ic) {
                a(i,j,k,0) = lsm_undefined;
                a(i,j,k,1) += Real(2.0);
            }
        });
    }
}

} // namespace

// Motivation: Noah-MP packs all its fields into one MultiFab, each with its own validity as
// indicator, and interpolates them in one call to move them in one ParallelCopy. Each field must
// come out exactly as it does on its own, whether it has its own indicator or shares one, and
// the call must not read past the one BCRec it is handed.
TEST(TypedInterp, PackedFieldsMatchTheirSingleFieldInterpolation)
{
    Levels L(false);
    constexpr int nf = 3;
    const MultiFab c0 = make_coarse(L, Real(1.0), Real(0.));
    const MultiFab c1 = make_coarse(L, lsm_undefined, Real(0.1));
    const MultiFab c2 = make_coarse(L, Real(5.0), Real(0.3));
    const MultiFab* cs[nf] = {&c0, &c1, &c2};

    MultiFab src(L.cba, L.cdm, nf, 1);
    for (int n = 0; n < nf; ++n) { MultiFab::Copy(src, *cs[n], 0, n, 1, 1); }

    // Own indicators: the validity of each value (only field 1 holds the sentinel)
    const MultiFab valid = make_validity(L, src);
    // Shared indicator: land/water
    MultiFab water(L.cba, L.cdm, 1, 1);
    MultiFab::Copy(water, c0, 1, 0, 1, 1);

    const MultiFab own    = interpolate_by_type(L, src, valid);
    const MultiFab shared = interpolate_by_type(L, src, water);

    for (int n = 0; n < nf; ++n) {
        MultiFab src_n(L.cba, L.cdm, 1, 1), valid_n(L.cba, L.cdm, 1, 1);
        MultiFab::Copy(src_n, src, n, 0, 1, 1);
        MultiFab::Copy(valid_n, valid, n, 0, 1, 1);
        const MultiFab own_n    = interpolate_by_type(L, src_n, valid_n);
        const MultiFab shared_n = interpolate_by_type(L, src_n, water);
        for (int j = 0; j < nc*ratio[1]; ++j) {
            for (int i = 0; i < nc*ratio[0]; ++i) {
                EXPECT_EQ(at(own, i, j, n), at(own_n, i, j)) << "n=" << n << " i=" << i << " j=" << j;
                EXPECT_EQ(at(shared, i, j, n), at(shared_n, i, j)) << "n=" << n << " i=" << i << " j=" << j;
            }
        }
    }
    // field 0 has no sentinel, so with its own indicator it is cell_cons_interp's, coast and all;
    // with the land/water one the coastal land cell has no x slope
    EXPECT_NE(at(own, 2*ncoast-2, 4, 0), at(own, 2*ncoast-1, 4, 0));
    EXPECT_EQ(at(shared, 2*ncoast-2, 4, 0), at(shared, 2*ncoast-1, 4, 0));
}

// Motivation: SLM interpolates a flux it has not computed yet (lsm_undefined) on land, where
// SurfaceLayer falls back to MOST. Typed by land/water alone, such a cell would lend its
// neighbours a slope toward the sentinel and be blended with them; as a type of its own (land/
// water type + 2) the sentinel comes through exactly and the land around it takes no slope
// toward it. The packed interface takes the field and its indicator as built by the caller.
TEST(TypedInterp, SentinelOnLandIsATypeOfItsOwn)
{
    Levels L(true);
    constexpr int ic_sentinel = 1;   // a land column holding the sentinel
    MultiFab crse = make_coarse(L, Real(1.0), Real(0.));
    put_sentinel_column(crse, ic_sentinel);
    MultiFab fine(L.fba, L.fdm, 2, 0);
    fine.setVal(Real(0.));
    BCRec bc;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { bc.setLo(d, BCType::int_dir); bc.setHi(d, BCType::int_dir); }
    interp_by_type_packed(fine, crse, 1, L.cgeom, L.fgeom, ratio, bc);
    MultiFab f(L.fba, L.fdm, 2, 0, MFInfo().SetArena(The_Pinned_Arena()));
    MultiFab::Copy(f, fine, 0, 0, 2, 0);
    Gpu::streamSynchronize();

    for (int j = 0; j < nc*ratio[1]; ++j) {
        for (int i = 0; i < nc*ratio[0]; ++i) {
            const int ic = i / ratio[0];
            if (ic == ic_sentinel) {
                EXPECT_EQ(at(f, i, j), lsm_undefined) << "i=" << i << " j=" << j;
            } else if (coarse_is_water(ic)) {
                EXPECT_EQ(at(f, i, j), Real(1.0)) << "i=" << i << " j=" << j;
            } else {
                EXPECT_TRUE(is_valid_lsm_value(at(f, i, j))) << "i=" << i << " j=" << j;
                EXPECT_LT(at(f, i, j), Real(1.)) << "i=" << i << " j=" << j;
            }
        }
    }
    // the land cells beside the sentinel column have no x slope
    for (int ic : {ic_sentinel-1, ic_sentinel+1}) {
        EXPECT_EQ(at(f, ic*2, 6), at(f, ic*2+1, 6)) << "ic=" << ic;
    }
}
