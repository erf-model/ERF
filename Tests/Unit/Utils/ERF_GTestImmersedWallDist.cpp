#include <cmath>
#include <limits>
#include <vector>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#include "ERF_Utils.H"

// Motivation (erf-model/ERF#4134): with terrain immersed forcing the mesh is flat, and the RANS wall
// distance is the height above the bottom of the domain, inside the solid. With the fraction-stress
// wall law (erf.if_wall_form = fraction_stress) immersed_wall_dist measures from the surface the
// wall law uses, h = z_bot + beta dz of the wall cell (the partial cell, or the fluid cell on a solid
// cell), with the wall cell at its fluid centroid, so that the length scale and the wall law see the
// same geometry.

using namespace amrex;

namespace {

constexpr int  nx = 8;
constexpr int  ny = 4;
constexpr int  nz = 12;
constexpr Real small_volfrac = Real(0.005);

// Node heights z_k of a flat mesh: stretched, so that no cell height is a round number
std::vector<Real>
stretched_levels ()
{
    std::vector<Real> z(nz+1);
    for (int k = 0; k <= nz; ++k) { z[k] = Real(10.)*k + Real(0.5)*k*k; }
    return z;
}

struct Fields
{
    MultiFab z_nd, z_cc, beta, wdist;
};

// Flat mesh from the node heights, solid fraction beta(i,k) from a table (the same in every j),
// on a BoxArray split in x and y but holding whole columns
Fields
make_fields (const std::vector<Real>& z_levels, const std::vector<Real>& beta_table)
{
    const Box domain(IntVect(0,0,0), IntVect(nx-1,ny-1,nz-1));
    BoxArray ba(domain);
    ba.maxSize(IntVect(2,2,nz));
    DistributionMapping dm(ba);
    BoxArray ba_nd(ba); ba_nd.surroundingNodes();

    Fields f{MultiFab(ba_nd, dm, 1, 1), MultiFab(ba, dm, 1, 1), MultiFab(ba, dm, 1, 1), MultiFab(ba, dm, 1, 1)};

    Gpu::DeviceVector<Real> zl(z_levels.size()), bt(beta_table.size());
    Gpu::copy(Gpu::hostToDevice, z_levels.begin(), z_levels.end(), zl.begin());
    Gpu::copy(Gpu::hostToDevice, beta_table.begin(), beta_table.end(), bt.begin());
    const Real* zp = zl.data();
    const Real* bp = bt.data();
    const int k_top = nz;
    const int n_x   = nx;

    for (MFIter mfi(f.z_nd); mfi.isValid(); ++mfi) {
        auto const z = f.z_nd.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            z(i,j,k) = zp[amrex::max(0, amrex::min(k, k_top))];
        });
    }
    for (MFIter mfi(f.z_cc); mfi.isValid(); ++mfi) {
        auto const zc = f.z_cc.array(mfi);
        auto const b  = f.beta.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int kk = amrex::max(0, amrex::min(k, k_top-1));
            const int ii = amrex::max(0, amrex::min(i, n_x-1));
            zc(i,j,k) = Real(0.5) * (zp[kk] + zp[kk+1]);
            b(i,j,k)  = bp[ii*k_top + kk];
        });
    }
    f.wdist.setVal(Real(-1.));
    Gpu::streamSynchronize();

    immersed_wall_dist(f.wdist, f.beta, f.z_nd, f.z_cc, domain, small_volfrac);
    Gpu::streamSynchronize();
    return f;
}

// The wall distance of every valid cell, on the host, indexed (i, j, k) -> (k*ny + j)*nx + i
std::vector<Real>
to_host (const MultiFab& wdist)
{
    const Box domain(IntVect(0,0,0), IntVect(nx-1,ny-1,nz-1));
    MultiFab one(BoxArray(domain), DistributionMapping(BoxArray(domain)), 1, 0);
    one.ParallelCopy(wdist, 0, 0, 1);
    std::vector<Real> out(static_cast<std::size_t>(nx*ny*nz), Real(-2.));
    for (MFIter mfi(one); mfi.isValid(); ++mfi) {
        Gpu::copy(Gpu::deviceToHost, one[mfi].dataPtr(), one[mfi].dataPtr() + one[mfi].size(), out.begin());
    }
    Gpu::streamSynchronize();
    return out;
}

Real d_at (const std::vector<Real>& d, int i, int j, int k) { return d[static_cast<std::size_t>((k*ny + j)*nx + i)]; }

Real
tol_for (Real z_max)
{
    return Real(8.) * std::numeric_limits<Real>::epsilon() * z_max;
}

// Fine-level boxes that do not hold whole columns (split at k = 4, dz = 10): columns 0-1 have the
// surface at 65 m, so their lower box is wholly solid; columns 2-3 at 25 m, so their upper box is
// wholly fluid and does not reach the bottom. Those two boxes must take the coarser level's wall
// height (given here as 3 m above the local one, so that the test sees which one is used); the
// other two see the surface and measure from it. Returns wdist and the wall height of each box.
struct SplitResult { std::vector<Real> d; std::vector<Real> hw; };

SplitResult
split_columns ()
{
    constexpr int sx = 4, sy = 2, sz = 12;
    const Real dz = Real(10.);
    const Box domain(IntVect(0,0,0), IntVect(sx-1,sy-1,sz-1));
    BoxList bl;
    for (int i0 = 0; i0 < sx; i0 += 2) {
        bl.push_back(Box(IntVect(i0,0,0), IntVect(i0+1,sy-1,3)));
        bl.push_back(Box(IntVect(i0,0,4), IntVect(i0+1,sy-1,sz-1)));
    }
    BoxArray ba(std::move(bl));
    DistributionMapping dm(ba);
    BoxArray ba_nd(ba); ba_nd.surroundingNodes();
    BoxList bl2d = ba.boxList();
    for (Box& b : bl2d) { b.setRange(2, 0); }
    BoxArray ba2d(std::move(bl2d));

    MultiFab z_nd(ba_nd, dm, 1, 1), z_cc(ba, dm, 1, 1), beta(ba, dm, 1, 1), wdist(ba, dm, 1, 1);
    MultiFab hw(ba2d, dm, 1, 0), hc(ba2d, dm, 1, 0);
    for (MFIter mfi(beta); mfi.isValid(); ++mfi) {
        auto const zn = z_nd.array(mfi);
        auto const zc = z_cc.array(mfi);
        auto const b  = beta.array(mfi);
        auto const c  = hc.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Real h = (i < 2) ? Real(65.) : Real(25.);
            zc(i,j,k) = dz * (k + Real(0.5));
            b(i,j,k)  = amrex::min(Real(1.), amrex::max(Real(0.), (h - dz*k) / dz));
        });
        ParallelFor(hc[mfi].box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            c(i,j,k) = ((i < 2) ? Real(65.) : Real(25.)) + Real(3.);
        });
        ParallelFor(z_nd[mfi].box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            zn(i,j,k) = dz * k;
        });
    }
    wdist.setVal(Real(-1.));
    hw.setVal(Real(-1.));
    Gpu::streamSynchronize();
    immersed_wall_dist(wdist, beta, z_nd, z_cc, domain, small_volfrac, &hw, &hc);
    Gpu::streamSynchronize();

    SplitResult r;
    MultiFab one(BoxArray(domain), DistributionMapping(BoxArray(domain)), 1, 0);
    one.ParallelCopy(wdist, 0, 0, 1);
    r.d.assign(static_cast<std::size_t>(sx*sy*sz), Real(-2.));
    for (MFIter mfi(one); mfi.isValid(); ++mfi) {
        Gpu::copy(Gpu::deviceToHost, one[mfi].dataPtr(), one[mfi].dataPtr() + one[mfi].size(), r.d.begin());
    }
    // wall height of each box at its (i = i0, j = 0) column, box order: (cols 0-1, low), (cols 0-1,
    // high), (cols 2-3, low), (cols 2-3, high)
    r.hw.assign(ba.size(), Real(-2.));
    for (MFIter mfi(hw); mfi.isValid(); ++mfi) {
        Real v = Real(-2.);
        Gpu::copy(Gpu::deviceToHost, hw[mfi].dataPtr(), hw[mfi].dataPtr() + 1, &v);
        r.hw[mfi.index()] = v;
    }
    Gpu::streamSynchronize();
    ParallelDescriptor::ReduceRealMax(r.hw.data(), static_cast<int>(r.hw.size()));
    return r;
}

} // namespace

// Every cell against a column-by-column reference on a stretched mesh: no solid, a face-aligned
// surface, partial cells, a fraction below and exactly at small_volfrac, and a solid block
// floating above fluid (the topmost solid cell wins). Solid cells get their distance to the wall.
TEST(ImmersedWallDist, MatchesTheWallCellColumnByColumn)
{
    const std::vector<Real> z = stretched_levels();
    std::vector<Real> beta(static_cast<std::size_t>(nx*nz), Real(0.));
    auto set = [&] (int i, int k, Real b) { beta[i*nz + k] = b; };
    // column i: expected surface height h and wall cell kw (the partial cell, or the fluid cell on
    // the solid; d = (z[kw+1] - h) / 2 there, |z_cc - h| elsewhere)
    std::vector<Real> h(nx);
    std::vector<int>  kw(nx);
    auto dz3 = z[4] - z[3];
    h[0] = z[0]; kw[0] = 0;                                                                  // no solid: the bottom of the domain
    for (int k = 0; k < 3; ++k) { set(1, k, Real(1.)); }  h[1] = z[3];                   kw[1] = 3;  // face-aligned surface at z[3]
    for (int k = 0; k < 3; ++k) { set(2, k, Real(1.)); }  set(2, 3, Real(0.25)); h[2] = z[3] + Real(0.25)*dz3; kw[2] = 3;
    for (int k = 0; k < 3; ++k) { set(3, k, Real(1.)); }  set(3, 3, Real(0.75)); h[3] = z[3] + Real(0.75)*dz3; kw[3] = 3;
    for (int k = 0; k < 3; ++k) { set(4, k, Real(1.)); }  set(4, 3, Real(0.004)); h[4] = z[3]; kw[4] = 3;      // sliver: fluid
    for (int k = 0; k < 3; ++k) { set(5, k, Real(1.)); }  set(5, 3, small_volfrac); h[5] = z[3] + small_volfrac*dz3; kw[5] = 3; // at the threshold: partial
    for (int k = 0; k < 5; ++k) { set(6, k, Real(1.)); }  set(6, 5, Real(0.5));  h[6] = z[5] + Real(0.5)*(z[6]-z[5]); kw[6] = 5;
    for (int k = 0; k < 2; ++k) { set(7, k, Real(1.)); }  set(7, 4, Real(1.));   h[7] = z[5]; kw[7] = 5;      // floating block: the top solid cell wins

    const Fields f = make_fields(z, beta);
    const std::vector<Real> d = to_host(f.wdist);
    const Real tol = tol_for(z[nz]);

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                SCOPED_TRACE(testing::Message() << "i = " << i << ", j = " << j << ", k = " << k);
                const Real z_cc = Real(0.5) * (z[k] + z[k+1]);
                const Real expect = (k == kw[i]) ? Real(0.5) * (z[k+1] - h[i]) : std::abs(z_cc - h[i]);
                EXPECT_NEAR(d_at(d, i, j, k), expect, tol);
            }
        }
    }
}

// A fine level whose boxes do not hold whole columns: boxes that see the surface measure from it,
// the others (wholly solid, or above the surface without reaching the bottom) from the coarser
// level's wall height
TEST(ImmersedWallDist, FineBoxesWithoutTheSurfaceUseTheCoarseWallHeight)
{
    constexpr int sx = 4, sy = 2, sz = 12;
    const Real dz = Real(10.);
    const Real tol = tol_for(dz*sz);
    const SplitResult r = split_columns();
    for (int i = 0; i < sx; ++i) {
    for (int j = 0; j < sy; ++j) {
    for (int k = 0; k < sz; ++k) {
        SCOPED_TRACE(testing::Message() << "i = " << i << ", k = " << k);
        const bool left   = (i < 2);
        const Real h      = left ? Real(65.) : Real(25.);
        const bool sees   = left ? (k >= 4) : (k < 4);
        const Real z_cc   = dz * (k + Real(0.5));
        Real expect;
        if (!sees) {
            expect = std::abs(z_cc - (h + Real(3.)));
        } else {
            const int kw = static_cast<int>(std::floor(h / dz));   // partial cell: the wall cell
            expect = (k == kw) ? Real(0.5) * (dz * (kw + 1) - h) : std::abs(z_cc - h);
        }
        EXPECT_NEAR(r.d[static_cast<std::size_t>((k*sy + j)*sx + i)], expect, tol);
    }}}
    EXPECT_NEAR(r.hw[0], Real(68.), tol);   // columns 0-1, low box: inside the solid
    EXPECT_NEAR(r.hw[1], Real(65.), tol);
    EXPECT_NEAR(r.hw[2], Real(25.), tol);
    EXPECT_NEAR(r.hw[3], Real(28.), tol);   // columns 2-3, high box: above the surface
}
