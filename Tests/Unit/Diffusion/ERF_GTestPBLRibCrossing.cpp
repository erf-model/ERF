// Unit tests for the bulk Richardson number scan that sets the YSUNew PBL height
// (find_rib_crossing, ERF_PBLRibCrossing.H; erf-model/ERF#4458), on synthetic columns
// whose answer is known in closed form:
//  * Rib linear in height on a stretched column: the interpolated height is Ribcr / slope
//    whatever the level spacing, and kpbl is the first level at or above Ribcr;
//  * a crossing at the top cell khi is interpolated, not sent to the floor;
//  * no crossing up to khi gives the top of the scan, z(khi), capped by pblh_max (WRF brint = 1);
//  * the surface cell and the next level both at or above Ribcr give the floor pblh_min, finite
//    (the interpolation used to divide 0 by 0 there); a surface cell above Ribcr with the next
//    level below it does not stop the scan, as in WRF, where the surface value only seeds it;
//  * the bounds apply in the order min then max, so pblh_min wins when pblh_max < pblh_min;
//  * the same answers from inside a ParallelFor, so the helper is usable in device code.

#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>

#include <ERF_PBLRibCrossing.H>

#include <gtest/gtest.h>

#include <cmath>
#include <utility>
#include <vector>

using namespace amrex;

namespace {

constexpr Real Ribcr = Real(0.25);
constexpr Real tol   = (sizeof(Real) == 8) ? Real(1.0e-12) : Real(1.0e-4);

// Ten levels, k = 0..9, on a stretched column z_k = 2 + 3 k + 0.5 k^2 (m).
constexpr int nz  = 10;
constexpr int khi = nz - 1;

Real z_stretched (int k) { return Real(2.0) + Real(3.0) * k + Real(0.5) * k * k; }

// Level k of an array, callable from host and device code (a plain host lambda handed to the
// host-device find_rib_crossing would draw nvcc's host-from-host-device warning).
struct LevelOf {
    const Real* v;
    AMREX_GPU_HOST_DEVICE Real operator() (int k) const noexcept { return v[k]; }
};

struct Column {
    std::vector<Real> rib;
    std::vector<Real> z;
    explicit Column (std::vector<Real> r) : rib(std::move(r)), z(nz) {
        for (int k = 0; k < nz; ++k) { z[k] = z_stretched(k); }
    }
    RibCrossing scan (int ksrf, Real pblh_min, Real pblh_max) const {
        return find_rib_crossing(ksrf, khi, Ribcr, LevelOf{rib.data()}, LevelOf{z.data()},
                                 pblh_min, pblh_max);
    }
};

} // namespace

// Rib = s z: the crossing height is Ribcr / s for any spacing, between levels 4 (z = 22)
// and 5 (z = 29.5) for s = 0.01, so kpbl = 5 and pblh = 25.
TEST(PBLRibCrossing, LinearRibOnStretchedColumn)
{
    std::vector<Real> r(nz);
    for (int k = 0; k < nz; ++k) { r[k] = Real(0.01) * z_stretched(k); }
    const RibCrossing c = Column(r).scan(0, Real(10.0), Real(1000.0));
    EXPECT_TRUE(c.crossed);
    EXPECT_EQ(c.kpbl, 5);
    EXPECT_NEAR(c.pblh, Real(25.0), tol * Real(25.0));
}

// The column starts at ksrf = 2 (an immersed solid below): the scan and the
// interpolation use levels 2 and up only, and the answer is unchanged.
TEST(PBLRibCrossing, LinearRibAboveImmersedSolid)
{
    std::vector<Real> r(nz);
    for (int k = 0; k < nz; ++k) { r[k] = Real(0.01) * z_stretched(k); }
    r[0] = Real(5.0); r[1] = Real(5.0);   // inside the solid: must not be read as a crossing
    const RibCrossing c = Column(r).scan(2, Real(10.0), Real(1000.0));
    EXPECT_TRUE(c.crossed);
    EXPECT_EQ(c.kpbl, 5);
    EXPECT_NEAR(c.pblh, Real(25.0), tol * Real(25.0));
}

// The first crossing is at the top cell: Rib = 0 up to level 8 (z = 58) and 0.5 at level 9
// (z = 69.5), so pblh = 58 + 11.5 * 0.25 / 0.5 = 63.75. The scan used to send this column
// to pblh_min (10).
TEST(PBLRibCrossing, CrossingAtTopCellIsInterpolated)
{
    std::vector<Real> r(nz, Real(0.0));
    r[khi] = Real(0.5);
    const RibCrossing c = Column(r).scan(0, Real(10.0), Real(1000.0));
    EXPECT_TRUE(c.crossed);
    EXPECT_EQ(c.kpbl, khi);
    EXPECT_NEAR(c.pblh, Real(63.75), tol * Real(63.75));
}

// Rib stays below Ribcr in every cell: the height is the top of the scan, z(khi) = 69.5,
// and kpbl = khi; with the 0.9 z(khi) cap it is 62.55. The scan used to give pblh_min (10).
TEST(PBLRibCrossing, NoCrossingGivesTopOfScan)
{
    std::vector<Real> r(nz, Real(-0.5));
    const Column col(r);
    const RibCrossing open = col.scan(0, Real(10.0), Real(1000.0));
    EXPECT_FALSE(open.crossed);
    EXPECT_EQ(open.kpbl, khi);
    EXPECT_NEAR(open.pblh, z_stretched(khi), tol * z_stretched(khi));

    const Real cap = Real(0.9) * z_stretched(khi);
    const RibCrossing capped = col.scan(0, Real(10.0), cap);
    EXPECT_FALSE(capped.crossed);
    EXPECT_EQ(capped.kpbl, khi);
    EXPECT_NEAR(capped.pblh, cap, tol * cap);
}

// The surface cell and the level above it are at or above Ribcr: the height is z(ksrf), raised
// to pblh_min, and kpbl = ksrf (WRF brint = 0).  The interpolation used to run with both ends on
// the surface cell, 0 / 0.
TEST(PBLRibCrossing, SurfaceCellAboveCriticalGivesFloor)
{
    std::vector<Real> r(nz, Real(1.0));
    const RibCrossing c = Column(r).scan(0, Real(10.0), Real(1000.0));
    EXPECT_TRUE(c.crossed);
    EXPECT_EQ(c.kpbl, 0);
    EXPECT_TRUE(std::isfinite(c.pblh));
    EXPECT_EQ(c.pblh, Real(10.0));

    // Exactly at Ribcr counts as a crossing too.
    r[0] = Ribcr;
    const RibCrossing at = Column(r).scan(0, Real(10.0), Real(1000.0));
    EXPECT_EQ(at.kpbl, 0);
    EXPECT_EQ(at.pblh, Real(10.0));
}

// The surface cell is above Ribcr but level 1 is below it: the scan goes on, as WRF's does, to the
// crossing between level 4 (z = 22, Rib = 0) and level 5 (z = 29.5, Rib = 0.5): pblh = 25.75.
TEST(PBLRibCrossing, SurfaceCellAboveCriticalDoesNotStopTheScan)
{
    std::vector<Real> r(nz, Real(0.0));
    r[0] = Real(1.0);
    for (int k = 5; k < nz; ++k) { r[k] = Real(0.5); }
    const RibCrossing c = Column(r).scan(0, Real(10.0), Real(1000.0));
    EXPECT_TRUE(c.crossed);
    EXPECT_EQ(c.kpbl, 5);
    EXPECT_NEAR(c.pblh, Real(25.75), tol * Real(25.75));
}

// A crossing exactly at a level gives that level's height.
TEST(PBLRibCrossing, RibEqualToCriticalAtALevel)
{
    std::vector<Real> r(nz, Real(0.0));
    r[3] = Ribcr;
    const RibCrossing c = Column(r).scan(0, Real(1.0), Real(1000.0));
    EXPECT_EQ(c.kpbl, 3);
    EXPECT_NEAR(c.pblh, z_stretched(3), tol * z_stretched(3));
}

// The bounds: an interpolated height below pblh_min is raised to it, one above pblh_max is
// lowered to it, and pblh_min wins when pblh_max < pblh_min (a column only a few cells deep).
// A column whose surface cell is the top cell (ksrf = khi) has no level to scan.
TEST(PBLRibCrossing, BoundsAndOrder)
{
    std::vector<Real> r(nz);
    for (int k = 0; k < nz; ++k) { r[k] = Real(0.01) * z_stretched(k); }
    const Column col(r);
    EXPECT_EQ(col.scan(0, Real(30.0), Real(1000.0)).pblh, Real(30.0));
    EXPECT_EQ(col.scan(0, Real(10.0), Real(20.0)).pblh,   Real(20.0));
    EXPECT_EQ(col.scan(0, Real(30.0), Real(20.0)).pblh,   Real(30.0));

    std::vector<Real> low(nz, Real(0.0));
    const RibCrossing top = Column(low).scan(khi, Real(75.0), Real(0.9) * z_stretched(khi));
    EXPECT_FALSE(top.crossed);
    EXPECT_EQ(top.kpbl, khi);
    EXPECT_EQ(top.pblh, Real(75.0));
}

namespace {

// Run the scan on ncol columns of nz levels inside a ParallelFor (a GPU kernel in a GPU
// build).  The launch sits at namespace scope: nvcc rejects an extended device lambda in
// the private TestBody of a TEST.
void scan_on_device (const std::vector<Real>& h_rib, const std::vector<Real>& h_z, int ncol,
                     std::vector<Real>& h, std::vector<int>& kp)
{
    Gpu::DeviceVector<Real> d_rib(h_rib.size()), d_z(h_z.size()), d_h(ncol);
    Gpu::DeviceVector<int>  d_k(ncol);
    Gpu::copy(Gpu::hostToDevice, h_rib.begin(), h_rib.end(), d_rib.begin());
    Gpu::copy(Gpu::hostToDevice, h_z.begin(), h_z.end(), d_z.begin());
    const Real* rib = d_rib.data();
    const Real* z   = d_z.data();
    Real* hout = d_h.data();
    int*  kout = d_k.data();
    const int  nlev  = nz;
    const int  ktop  = khi;
    const Real rcrit = Ribcr;
    ParallelFor(ncol, [=] AMREX_GPU_DEVICE (int n) noexcept
    {
        auto rib_at = [=] (int k) { return rib[n * nlev + k]; };
        auto z_at   = [=] (int k) { return z[k]; };
        const RibCrossing c = find_rib_crossing(0, ktop, rcrit, rib_at, z_at, Real(10.0), Real(1000.0));
        hout[n] = c.pblh;
        kout[n] = c.kpbl;
    });
    h.resize(ncol);
    kp.resize(ncol);
    Gpu::copy(Gpu::deviceToHost, d_h.begin(), d_h.end(), h.begin());
    Gpu::copy(Gpu::deviceToHost, d_k.begin(), d_k.end(), kp.begin());
    Gpu::streamSynchronize();
}

} // namespace

// The same three answers from inside a ParallelFor.
TEST(PBLRibCrossing, DeviceKernel)
{
    constexpr int ncol = 3;
    std::vector<Real> h_rib(ncol * nz), h_z(nz);
    for (int k = 0; k < nz; ++k) {
        h_z[k] = z_stretched(k);
        h_rib[0 * nz + k] = Real(0.01) * z_stretched(k);       // interior crossing: 25
        h_rib[1 * nz + k] = (k == khi) ? Real(0.5) : Real(0);  // crossing at khi: 63.75
        h_rib[2 * nz + k] = Real(-0.5);                         // no crossing: 69.5
    }
    std::vector<Real> h;
    std::vector<int>  kp;
    scan_on_device(h_rib, h_z, ncol, h, kp);
    EXPECT_NEAR(h[0], Real(25.0),  tol * Real(25.0));
    EXPECT_NEAR(h[1], Real(63.75), tol * Real(63.75));
    EXPECT_NEAR(h[2], Real(69.5),  tol * Real(69.5));
    EXPECT_EQ(kp[0], 5);
    EXPECT_EQ(kp[1], khi);
    EXPECT_EQ(kp[2], khi);
}
