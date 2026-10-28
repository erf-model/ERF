#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>

#include "ERF_MetgridUtils.H"

// A named namespace: the device lambda below needs an enclosing function with external
// linkage under nvcc
namespace metgrid_force_sfc_test {

// ERF levels: nodes 100, 180, 260, ... m, so the cell centres are 140, 220, 300, ... m
constexpr int kNewCells = 12;

#ifdef AMREX_USE_FLOAT
constexpr amrex::Real kTol = amrex::Real(1.0e-4);
#else
constexpr amrex::Real kTol = amrex::Real(1.0e-12);
#endif
constexpr amrex::Real kSfc = amrex::Real(100.0);

// A data profile that is not linear in z, so dropping an origin level changes the result
amrex::Real profile (amrex::Real z) { return (z / amrex::Real(100.0)) * (z / amrex::Real(100.0)); }

// Interpolate one origin column (surface first, then the other origin levels in the order
// given) onto the ERF levels with linear interpolation and erf.metgrid_force_sfc_k = force.
// interpolate_column_metgrid is a device function, so it runs in a one-cell ParallelFor on
// FABs that the host can also read and write (pinned memory in a GPU build).
std::vector<amrex::Real> interpolate (const std::vector<amrex::Real>& orig_z, int force,
                                      bool retain_sfc = false, bool extra_top = true)
{
    amrex::Arena* arena = amrex::The_Arena();
#ifdef AMREX_USE_GPU
    arena = amrex::The_Pinned_Arena();
#endif

    // Optionally an extra origin level 1000 m above the last, which the ERF column does not reach
    const int nz_orig = static_cast<int>(orig_z.size());
    const int nz_box  = extra_top ? nz_orig + 1 : nz_orig;
    const amrex::Box obox(amrex::IntVect(0,0,0), amrex::IntVect(0,0,nz_box-1));
    amrex::FArrayBox oz(obox, 1, arena), od(obox, 1, arena);
    for (int k = 0; k < nz_box; ++k) {
        const amrex::Real z = (k < nz_orig) ? orig_z[k] : orig_z.back() + amrex::Real(1000.0);
        oz(amrex::IntVect(0,0,k)) = z;
        od(amrex::IntVect(0,0,k)) = profile(z);
    }

    const amrex::Box nbox(amrex::IntVect(0,0,0), amrex::IntVect(1,1,kNewCells));
    amrex::FArrayBox nz(nbox, 1, arena);
    for (int k = 0; k <= kNewCells; ++k) {
        for (int j = 0; j <= 1; ++j) {
            for (int i = 0; i <= 1; ++i) {
                nz(amrex::IntVect(i,j,k)) = kSfc + amrex::Real(80.0) * static_cast<amrex::Real>(k);
            }
        }
    }
    const amrex::Box dbox(amrex::IntVect(0,0,0), amrex::IntVect(0,0,kNewCells-1));
    amrex::FArrayBox nd(dbox, 1, arena);
    nd.setVal<amrex::RunOn::Host>(amrex::Real(0.0));

    const auto oz_a = oz.const_array();
    const auto od_a = od.const_array();
    const auto nz_a = nz.const_array();
    const auto nd_a = nd.array();
    amrex::ParallelFor(amrex::Box(amrex::IntVect(0,0,0), amrex::IntVect(0,0,0)),
    [=] AMREX_GPU_DEVICE (int, int, int) noexcept
    {
        interpolate_column_metgrid(/*use_below_sfc=*/false, /*use_sfc=*/true, /*exp_interp=*/false,
                                   retain_sfc, /*proximity=*/amrex::Real(0.0), /*order=*/1,
                                   force, 0, 0, kNewCells-1, 0, 0, 'U', 'M', oz_a, od_a, nz_a, nd_a);
    });
    amrex::Gpu::streamSynchronize();

    std::vector<amrex::Real> out(kNewCells);
    for (int k = 0; k < kNewCells; ++k) { out[k] = nd(amrex::IntVect(0,0,k)); }
    return out;
}

// Origin levels above the surface: 250 m lies between the second (220 m) and the third
// (300 m) ERF cell centre, which is where the two branches used to disagree
const std::vector<amrex::Real> kAbove = {amrex::Real(250.0), amrex::Real(400.0), amrex::Real(600.0),
                                         amrex::Real(900.0), amrex::Real(1300.0), amrex::Real(1800.0)};

std::vector<amrex::Real> column_without_below_sfc ()
{
    std::vector<amrex::Real> z = {kSfc};
    z.insert(z.end(), kAbove.begin(), kAbove.end());
    return z;
}

std::vector<amrex::Real> column_with_below_sfc ()
{
    // An origin level below the surface (50 m) sends the interpolation down its other branch
    std::vector<amrex::Real> z = {kSfc, amrex::Real(50.0)};
    z.insert(z.end(), kAbove.begin(), kAbove.end());
    return z;
}

} // namespace metgrid_force_sfc_test

using namespace metgrid_force_sfc_test;

// Motivation: erf.metgrid_force_sfc_k = N uses the origin surface through N ERF levels, so
// origin levels between the surface and the centre of ERF level N-1 (0-based) are dropped, as
// in WRF's vert_interp. Here N = 2: the 250 m level lies above the 220 m centre of level 1 and
// must be kept, so level 1 is interpolated between the surface (100 m) and 250 m. The branch
// for a column whose surface is the lowest origin level compared with level N instead and
// dropped the 250 m level (erf-model/ERF#4607).
TEST(MetgridForceSfc, KeepsTheFirstOriginLevelAboveLevelNMinus1)
{
    const auto v = interpolate(column_without_below_sfc(), 2);
    const amrex::Real w = (amrex::Real(220.0) - kSfc) / (amrex::Real(250.0) - kSfc);
    EXPECT_NEAR(v[1], profile(kSfc) + w * (profile(amrex::Real(250.0)) - profile(kSfc)), kTol * profile(amrex::Real(250.0)));
}

// Motivation: with metgrid_use_below_sfc = false an origin level below the surface is not
// used at all, so adding one must not change the result: both branches of the column
// interpolation must drop the same levels for every metgrid_force_sfc_k, up to the number
// of ERF levels (the largest value erf.metgrid_force_sfc_k may take).
TEST(MetgridForceSfc, BothBranchesAgree)
{
    for (int force = 0; force <= kNewCells; ++force) {
        const auto a = interpolate(column_without_below_sfc(), force);
        const auto b = interpolate(column_with_below_sfc(), force);
        for (int k = 0; k < kNewCells; ++k) {
            EXPECT_NEAR(a[k], b[k], kTol * (amrex::Real(1.0) + std::abs(a[k]))) << "metgrid_force_sfc_k = " << force << ", level " << k;
        }
    }
}

// Motivation: erf.metgrid_retain_sfc sets the lowest ERF level to the origin surface value.
// With an origin level below the surface, the first of the levels as ordered is that
// below-ground level, not the surface, and it was taken instead.
TEST(MetgridForceSfc, RetainSfcTakesTheSurfaceValue)
{
    const auto v = interpolate(column_with_below_sfc(), 0, /*retain_sfc=*/true);
    EXPECT_NEAR(v[0], profile(kSfc), kTol * profile(kSfc));
}

// Motivation: the top origin level must be used.  The column interpolation read origin
// levels 0 .. ubound-1, so an ERF level between the two highest origin levels was
// extrapolated from those below them.
TEST(MetgridForceSfc, UsesTheTopOriginLevel)
{
    // The top ERF cell centre, 1020 m, lies between the two highest origin levels
    const std::vector<amrex::Real> z = {kSfc, amrex::Real(300.0), amrex::Real(600.0), amrex::Real(1100.0)};
    const auto v = interpolate(z, 0, /*retain_sfc=*/false, /*extra_top=*/false);
    const amrex::Real z_top = amrex::Real(1020.0);
    const amrex::Real w = (z_top - amrex::Real(600.0)) / (amrex::Real(1100.0) - amrex::Real(600.0));
    const amrex::Real expected = profile(amrex::Real(600.0))
                               + w * (profile(amrex::Real(1100.0)) - profile(amrex::Real(600.0)));
    EXPECT_NEAR(v[kNewCells-1], expected, kTol * expected);
}
