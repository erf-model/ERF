#include <cmath>
#include <limits>

#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Vector.H>

#include <gtest/gtest.h>

#include "ERF_TerrainMetrics.H"

// Motivation: erf.initial_dz without erf.grid_stretching_ratio asks init_zlevels to solve for
// the ratio, so that an idealized or metgrid case can name its first layer and its domain top
// and get a grid that closes on both.  An explicit ratio must keep building the grid it always
// has, whatever top that reaches.

using namespace amrex;

namespace {

constexpr int nz = 32;

Geometry
make_geom (Real zlo, Real zhi)
{
    const Box domain(IntVect(0,0,0), IntVect(7,7,nz-1));
    const RealBox rb(Real(0.), Real(0.), zlo, Real(800.), Real(800.), zhi);
    return Geometry(domain, rb, CoordSys::cartesian, {1, 1, 0});
}

struct ZLevels
{
    Vector<Real> zlevels;
    Vector<Real> dz;
};

ZLevels
build (Real zhi, Real ratio, Real zlo, Real dz0)
{
    Vector<Geometry> geom{make_geom(zlo, zhi)};
    Vector<IntVect>  ref_ratio;

    Vector<Vector<Real>> zlevels_stag(1);
    Vector<Vector<Real>> stretched_dz_h(1);
    Vector<Gpu::DeviceVector<Real>> stretched_dz_d(1);

    init_zlevels(zlevels_stag, stretched_dz_h, stretched_dz_d, geom, ref_ratio,
                 ratio, dz0);

    return {zlevels_stag[0], stretched_dz_h[0]};
}

Real
allowance (Real scale)
{
    return Real(16.0) * nz * std::numeric_limits<Real>::epsilon() * std::abs(scale);
}

} // namespace

// Motivation: the solved ratio has to close the column on prob_hi, start one initial_dz above
// the surface, and give back the ratio an explicit run of the same grid would have used.
TEST(InitZLevels, SolvesTheRatioThatClosesOnProbHi)
{
    // 32 layers from 19.503737114205 m growing by 1.03 sum to 1024 m
    const Real dz0  = Real(19.503737114205);
    const Real ztop = Real(1024.0);
    const auto z = build(ztop, Real(0.0), Real(0.0), dz0);

    ASSERT_EQ(static_cast<int>(z.zlevels.size()), nz+1);
    ASSERT_EQ(static_cast<int>(z.dz.size()), nz);
    EXPECT_EQ(z.zlevels[0], Real(0.0));
    EXPECT_NEAR(z.zlevels[nz], ztop, allowance(ztop));
    EXPECT_NEAR(z.dz[0], dz0, allowance(ztop));

    for (int k(0); k<nz; ++k) {
        EXPECT_GT(z.dz[k], Real(0.0)) << "at k = " << k;
        EXPECT_NEAR(z.dz[k], z.zlevels[k+1] - z.zlevels[k], allowance(ztop)) << "at k = " << k;
    }
    // The top layer is assigned z_top exactly rather than by accumulation
    for (int k(1); k<nz-1; ++k) {
        EXPECT_NEAR(z.dz[k]/z.dz[k-1], Real(1.03), Real(1.e-6)) << "at k = " << k;
    }
}

// Motivation: the grid is built above geometry.prob_lo[2], so the column it has to fill is
// prob_hi - prob_lo rather than prob_hi.
TEST(InitZLevels, SolvedRatioStartsAtProbLo)
{
    const Real zlo  = Real(100.0);
    const Real ztop = Real(1124.0);
    const auto z = build(ztop, Real(0.0), zlo, Real(19.503737114205));

    EXPECT_EQ(z.zlevels[0], zlo);
    EXPECT_NEAR(z.zlevels[nz], ztop, allowance(ztop));
    for (int k(1); k<nz-1; ++k) {
        EXPECT_NEAR(z.dz[k]/z.dz[k-1], Real(1.03), Real(1.e-6)) << "at k = " << k;
    }
}

// Motivation: with an explicit ratio the grid also starts at geometry.prob_lo[2], which may be
// below zero (#4442); a uniform ratio from prob_lo = -1 must give back the uniform levels.
TEST(InitZLevels, AnExplicitRatioStartsAtProbLo)
{
    const Real zlo = Real(-1.0);
    const Real dz0 = Real(2.0) / nz;
    const auto z = build(Real(1.0), Real(1.0), zlo, dz0);

    EXPECT_EQ(z.zlevels[0], zlo);
    for (int k(0); k<=nz; ++k) {
        EXPECT_NEAR(z.zlevels[k], zlo + k*dz0, allowance(Real(1.0))) << "at k = " << k;
    }
}

// Motivation: an explicit ratio is a request for that grid, not for one that closes on
// prob_hi, so it must still end wherever its layers do.
TEST(InitZLevels, AnExplicitRatioIsNotSolvedFor)
{
    const Real dz0  = Real(20.0);
    const auto z = build(Real(1024.0), Real(1.03), Real(0.0), dz0);

    Real top = Real(0.0);
    Real dz  = dz0;
    for (int k(0); k<nz; ++k) { top += dz; dz *= Real(1.03); }

    EXPECT_NEAR(z.zlevels[nz], top, allowance(top));
    EXPECT_GT(z.zlevels[nz], Real(1024.0));
}

// Motivation: without initial_dz or a ratio the grid stays the uniform one from the geometry.
TEST(InitZLevels, NeitherGivenIsUniform)
{
    const auto z = build(Real(1024.0), Real(0.0), Real(0.0), Real(0.0));
    for (int k(0); k<nz; ++k) {
        EXPECT_NEAR(z.dz[k], Real(32.0), allowance(Real(1024.0))) << "at k = " << k;
    }
}
