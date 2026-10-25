#include <gtest/gtest.h>

#include <vector>

#include <AMReX_FArrayBox.H>

#include "ERF_MetgridUtils.H"

namespace {

// ERF levels: nodes 100, 180, 260, ... m, so the cell centres are 140, 220, 300, ... m
constexpr int kNewCells = 12;
constexpr amrex::Real kSfc = 100.0;

// A data profile that is not linear in z, so dropping an origin level changes the result
amrex::Real profile (amrex::Real z) { return (z / 100.0) * (z / 100.0); }

// Interpolate one origin column (surface first, then the other origin levels in the order
// given) onto the ERF levels with linear interpolation and erf.metgrid_force_sfc_k = force.
std::vector<amrex::Real> interpolate (const std::vector<amrex::Real>& orig_z, int force)
{
    // interpolate_column_metgrid reads origin levels 0 .. ubound-1, so add one unused level
    const int nz_orig = static_cast<int>(orig_z.size());
    const amrex::Box obox(amrex::IntVect(0,0,0), amrex::IntVect(0,0,nz_orig));
    amrex::FArrayBox oz(obox, 1), od(obox, 1);
    for (int k = 0; k <= nz_orig; ++k) {
        const amrex::Real z = (k < nz_orig) ? orig_z[k] : orig_z.back() + 1000.0;
        oz(amrex::IntVect(0,0,k)) = z;
        od(amrex::IntVect(0,0,k)) = profile(z);
    }

    const amrex::Box nbox(amrex::IntVect(0,0,0), amrex::IntVect(1,1,kNewCells));
    amrex::FArrayBox nz(nbox, 1);
    for (int k = 0; k <= kNewCells; ++k) {
        for (int j = 0; j <= 1; ++j) {
            for (int i = 0; i <= 1; ++i) {
                nz(amrex::IntVect(i,j,k)) = kSfc + 80.0 * k;
            }
        }
    }
    const amrex::Box dbox(amrex::IntVect(0,0,0), amrex::IntVect(0,0,kNewCells-1));
    amrex::FArrayBox nd(dbox, 1);
    nd.setVal<amrex::RunOn::Host>(0.0);

    interpolate_column_metgrid(/*use_below_sfc=*/false, /*use_sfc=*/true, /*exp_interp=*/false,
                               /*retain_sfc=*/false, /*proximity=*/0.0, /*order=*/1, force,
                               0, 0, kNewCells-1, 0, 0, 'U', 'M',
                               oz.const_array(), od.const_array(), nz.const_array(), nd.array());

    std::vector<amrex::Real> out(kNewCells);
    for (int k = 0; k < kNewCells; ++k) { out[k] = nd(amrex::IntVect(0,0,k)); }
    return out;
}

// Origin levels above the surface: 250 m lies between the second (220 m) and the third
// (300 m) ERF cell centre, which is where the two branches used to disagree
const std::vector<amrex::Real> kAbove = {250.0, 400.0, 600.0, 900.0, 1300.0, 1800.0};

std::vector<amrex::Real> column_without_below_sfc ()
{
    std::vector<amrex::Real> z = {kSfc};
    z.insert(z.end(), kAbove.begin(), kAbove.end());
    return z;
}

std::vector<amrex::Real> column_with_below_sfc ()
{
    // An origin level below the surface (50 m) sends the interpolation down its other branch
    std::vector<amrex::Real> z = {kSfc, 50.0};
    z.insert(z.end(), kAbove.begin(), kAbove.end());
    return z;
}

} // namespace

// Motivation: erf.metgrid_force_sfc_k = N uses the origin surface through N ERF levels, so
// origin levels between the surface and the centre of ERF level N-1 (0-based) are dropped, as
// in WRF's vert_interp. Here N = 2: the 250 m level lies above the 220 m centre of level 1 and
// must be kept, so level 1 is interpolated between the surface (100 m) and 250 m. The branch
// for a column whose surface is the lowest origin level compared with level N instead and
// dropped the 250 m level (erf-model/ERF#4607).
TEST(MetgridForceSfc, KeepsTheFirstOriginLevelAboveLevelNMinus1)
{
    const auto v = interpolate(column_without_below_sfc(), 2);
    const amrex::Real w = (220.0 - kSfc) / (250.0 - kSfc);
    EXPECT_NEAR(v[1], profile(kSfc) + w * (profile(250.0) - profile(kSfc)), 1.0e-12);
}

// Motivation: with metgrid_use_below_sfc = false an origin level below the surface is not
// used at all, so adding one must not change the result: both branches of the column
// interpolation must drop the same levels for every metgrid_force_sfc_k.
TEST(MetgridForceSfc, BothBranchesAgree)
{
    for (int force = 0; force <= 4; ++force) {
        const auto a = interpolate(column_without_below_sfc(), force);
        const auto b = interpolate(column_with_below_sfc(), force);
        for (int k = 0; k < kNewCells; ++k) {
            EXPECT_NEAR(a[k], b[k], 1.0e-12) << "metgrid_force_sfc_k = " << force << ", level " << k;
        }
    }
}
