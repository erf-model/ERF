#include <cmath>
#include <limits>

#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <gtest/gtest.h>

#include "ERF_NodalReconstruction.H"

// Motivation: reconstruct_nodal_height_slice ran its BaseFab copies and fills as
// copy<RunOn::Host>/setVal<RunOn::Host>, which on SYCL builds abort
// unconditionally ("HOST_DEVICE disabled for Intel"): a WRFInput/Metgrid init
// with the default avg_grid_faces_to_nodes = false died at init on Aurora.
// These tests pin the behavior the replacement host loops must keep, on both
// solver paths (exact de-averaging, and the CG fallback that rough terrain
// forces), plus the slice-into-level fill.

using namespace amrex;

namespace {

constexpr int nx = 16;
constexpr int ny = 12;
constexpr int nz = 4;

// Nodal heights of an analytic field: 100 + 10 sin(3x) + 5xy.  Its exact avg4
// average gives the w-face slice.  Face averages do not determine a unique
// nodal field (checkerboard modes lie in the null space of the averaging
// operator), so the reconstruction is not required to return this field; it is
// the informational reference in DeviationFromAnalytic below.
Real
analytic_nodal (Real x, Real y)
{
    return Real(100.) + Real(10.)*std::sin(Real(3.)*x) + Real(5.)*x*y;
}

// The w-face slice whose exact avg4 de-averaging is the analytic nodal field.
FArrayBox
make_smooth_slice (const Box& wface_box)
{
    FArrayBox z_w(wface_box, 1, The_Managed_Arena());
    const Real dx = Real(1.)/nx;
    const Real dy = Real(1.)/ny;
    auto const arr = z_w.array();
    LoopOnCpu(wface_box, [=] (int i, int j, int k) noexcept
    {
        const Real x = i*dx, y = j*dy;
        const Real nd_lo = analytic_nodal(x, y);
        const Real nd_hi = analytic_nodal(x + dx, y);
        const Real nd_lh = analytic_nodal(x, y + dy);
        const Real nd_hh = analytic_nodal(x + dx, y + dy);
        arr(i,j,0) = Real(0.25)*(nd_lo + nd_hi + nd_lh + nd_hh);
    });
    return z_w;
}

Real
avg4_identity_error (const FArrayBox& nodal, const FArrayBox& wface,
                     const Box& wface_box)
{
    auto const g = nodal.const_array();
    auto const w = wface.const_array();
    Real worst = 0;
    LoopOnCpu(wface_box, [&] (int i, int j, int k) noexcept
    {
        const Real a4 = Real(0.25)*(g(i,j,0) + g(i+1,j,0) + g(i,j+1,0) + g(i+1,j+1,0));
        worst = amrex::max(worst, std::abs(a4 - w(i,j,0)));
    });
    return worst;
}

} // namespace

// The exact de-averaging path on a smooth slice: the averaging identity must
// close to round-off.  The reported tolerance is 4x machine epsilon scaled by
// the field magnitude (the identity combines four nodal values and a
// pre-averaged face value).
TEST(NodalReconstruction, ExactDeAveragingIdentity)
{
    const Box domain(IntVect(0,0,0), IntVect(nx-1, ny-1, nz-1));
    const RealBox rb(Real(0.), Real(0.), Real(0.), Real(1.), Real(1.), Real(1.));
    const Geometry geom(domain, rb, CoordSys::cartesian, {0, 0, 0});
    Box wface_box = convert(domain, IntVect(0,0,1));
    wface_box.setSmall(2, 0);
    wface_box.setBig(2, 0);

    const FArrayBox z_w = make_smooth_slice(wface_box);
    const FArrayBox got = reconstruct_nodal_height_slice(wface_box, geom, z_w, 0, "test z-faces");

    const Real tol = Real(16.)*std::numeric_limits<Real>::epsilon()*Real(200.);
    EXPECT_LT(avg4_identity_error(got, z_w, wface_box), tol);
}

// A second call on the same input must reproduce the first exactly: the solve
// is deterministic, so a caller that reconstructs twice (e.g. two levels, or a
// restart) gets the same terrain.
TEST(NodalReconstruction, DeterministicAcrossCalls)
{
    const Box domain(IntVect(0,0,0), IntVect(nx-1, ny-1, nz-1));
    const RealBox rb(Real(0.), Real(0.), Real(0.), Real(1.), Real(1.), Real(1.));
    const Geometry geom(domain, rb, CoordSys::cartesian, {0, 0, 0});
    Box wface_box = convert(domain, IntVect(0,0,1));
    wface_box.setSmall(2, 0);
    wface_box.setBig(2, 0);

    const FArrayBox z_w = make_smooth_slice(wface_box);
    const FArrayBox got1 = reconstruct_nodal_height_slice(wface_box, geom, z_w, 0, "test z-faces");
    const FArrayBox got2 = reconstruct_nodal_height_slice(wface_box, geom, z_w, 0, "test z-faces");

    auto const a1 = got1.const_array();
    auto const a2 = got2.const_array();
    Real worst = 0;
    LoopOnCpu(got1.box(), [&] (int i, int j, int k) noexcept
    {
        worst = amrex::max(worst, std::abs(a1(i,j,k) - a2(i,j,k)));
    });
    EXPECT_EQ(worst, Real(0.0));
}

// A checkered-sawtooth slice is rough enough to force the CG fallback.  The
// fallback is regularized, so only a usable field is asserted: finite
// everywhere and the averaging identity within the loose bound the solver's
// own residual tolerance implies.
TEST(NodalReconstruction, RoughSliceFallsBackToCG)
{
    const Box domain(IntVect(0,0,0), IntVect(nx-1, ny-1, nz-1));
    const RealBox rb(Real(0.), Real(0.), Real(0.), Real(1.), Real(1.), Real(1.));
    const Geometry geom(domain, rb, CoordSys::cartesian, {0, 0, 0});
    Box wface_box = convert(domain, IntVect(0,0,1));
    wface_box.setSmall(2, 0);
    wface_box.setBig(2, 0);

    FArrayBox rough(wface_box, 1, The_Managed_Arena());
    {
        auto const arr = rough.array();
        const Real dx = Real(1.)/nx;
        const Real dy = Real(1.)/ny;
        LoopOnCpu(wface_box, [=] (int i, int j, int k) noexcept
        {
            // Steep checkered sawtooth: direct interpolation is badly wrong and
            // the exact de-averaging overshoots the source range, so the
            // solver takes the CG fallback.  Deterministic values, no RNG, so
            // the test is reproducible under any decomposition.
            const Real x = i*dx, y = j*dy;
            arr(i,j,0) = Real(500.) + Real(200.)*x*(Real(1.) - y)
                       + ((i + j) % 2 == 0 ? Real(400.) : Real(-400.));
        });
    }

    const FArrayBox got = reconstruct_nodal_height_slice(wface_box, geom, rough, 0, "rough z-faces");

    auto const g = got.const_array();
    Real range_lo = std::numeric_limits<Real>::max();
    Real range_hi = std::numeric_limits<Real>::lowest();
    LoopOnCpu(got.box(), [&] (int i, int j, int k) noexcept
    {
        const Real v = g(i,j,k);
        ASSERT_TRUE(std::isfinite(v)) << "at (" << i << ", " << j << ")";
        range_lo = amrex::min(range_lo, v);
        range_hi = amrex::max(range_hi, v);
    });
    // The nodal field must stay of the order of the input (a diverged solve
    // would produce garbage well outside it).
    EXPECT_LT(range_hi, Real(500.) + Real(200.) + Real(400.) + Real(100.));
    EXPECT_GT(range_lo, Real(500.) - Real(200.) - Real(400.) - Real(100.));

    const Real tol = Real(16.)*std::numeric_limits<Real>::epsilon()*Real(200.);
    EXPECT_LT(avg4_identity_error(got, rough, wface_box), tol);
}

// fill_nodal_level_from_slice must copy the slice into the level's k = 0 nodal
// plane everywhere, clamping lateral ghost and out-of-domain nodes to the
// nearest valid node.  The MultiFab mirrors the real allocation in
// ERF_MakeNewArrays.cpp: surroundingNodes boxes with one ghost cell.
TEST(NodalReconstruction, FillNodalLevelFromSliceClampsGhosts)
{
    const Box domain(IntVect(0,0,0), IntVect(nx-1, ny-1, nz-1));
    const RealBox rb(Real(0.), Real(0.), Real(0.), Real(1.), Real(1.), Real(1.));
    const Geometry geom(domain, rb, CoordSys::cartesian, {0, 0, 0});
    Box wface_box = convert(domain, IntVect(0,0,1));
    wface_box.setSmall(2, 0);
    wface_box.setBig(2, 0);

    const FArrayBox z_w = make_smooth_slice(wface_box);
    const FArrayBox got = reconstruct_nodal_height_slice(wface_box, geom, z_w, 0, "test z-faces");

    BoxArray ba(domain);
    ba.maxSize(8);
    ba.surroundingNodes();
    MultiFab z_phys_nd(ba, DistributionMapping(ba), 1, 1);
    fill_nodal_level_from_slice(z_phys_nd, 0, got);

    const Box slb = got.box();
    const int ilo = slb.smallEnd(0), ihi = slb.bigEnd(0);
    const int jlo = slb.smallEnd(1), jhi = slb.bigEnd(1);
    auto const sarr = got.const_array();
    for (MFIter mfi(z_phys_nd); mfi.isValid(); ++mfi) {
        auto const zarr = z_phys_nd.const_array(mfi);
        // The fill writes only the k = 0 plane; compare just that plane.
        const Box cbx = makeSlab(mfi.growntilebox(), 2, 0);
        LoopOnCpu(cbx, [&] (int i, int j, int k) noexcept
        {
            const int ii = amrex::max(ilo, amrex::min(i, ihi));
            const int jj = amrex::max(jlo, amrex::min(j, jhi));
            ASSERT_EQ(zarr(i,j,k), sarr(ii,jj,0))
                << "at (" << i << ", " << j << ")";
        });
    }
}
