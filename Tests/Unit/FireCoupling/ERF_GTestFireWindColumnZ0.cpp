// Unit tests for the column roughness of the fire wind sampler (ERF_FireWindExtract.H,
// fill_fire_wind_from_interpolation with col_z0).
//
// A 4 x 4 atmospheric level on flat ground carries a wind that is uniform in height but differs
// from column to column, and a roughness length that differs from column to column. With a
// first cell center at 10 m above ground, the wind at 6.1 m comes from the neutral log profile
// below the first center, U(6.1) = U ln(6.1 / z0) / ln(10 / z0), and with a sampling height of
// 30 m (between centers, where the wind is U) the resampled wind is U ln(6.1 / z0) / ln(30 / z0).
// Both must use the roughness of each column. The expected values are written out from these
// closed forms, not from the sampler.

#include <ERF_FireWindExtract.H>
#include <ERF_FireGrid.H>

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>

using namespace amrex;

namespace {

constexpr int  nx    = 4;
constexpr int  ny    = 4;
constexpr int  nz    = 6;
constexpr int  C     = 3;               // fire cells per atmospheric cell
constexpr Real dz    = Real(20.0);      // first cell center at 10 m
constexpr Real z_ref = Real(6.1);
constexpr Real z0_scalar = Real(0.1);   // the scalar roughness an engine passes
// Relative tolerance: a few dozen roundings of the working precision
constexpr Real rel_tol = Real(64) * std::numeric_limits<Real>::epsilon();

// Column wind: cell-center u of column i is u_a + u_b (i + 1/2); v is fixed
constexpr Real u_a = Real(4.0);
constexpr Real u_b = Real(0.5);
constexpr Real v_c = Real(-1.5);

// Column roughness [m], 0.02 to 1.02 over the columns and their ghost ring
Real z0_of (int i, int j) { return Real(0.02) + Real(0.05) * Real(i + 1) + Real(0.15) * Real(j + 1); }
// The column the roughness is read from: across a non-periodic face of the region
// (columns lo to lo + n - 1) the edge column, across a periodic face the periodic image
int z0_column (int i, int lo, int n, bool periodic)
{
    if (periodic) { return lo + (((i - lo) % n) + n) % n; }
    return std::max(lo, std::min(i, lo + n - 1));
}
Real u_of (int i) { return u_a + u_b * (Real(i) + Real(0.5)); }

struct Sampler {
    BoxArray ba_atm, ba_atm_2d;
    DistributionMapping dm;
    FireGrid fg;
    MultiFab xvel, yvel, zcc, z0col;
    MultiFab fire_wind, extract_z, surf_z, col_ground;

    bool periodic = false;

    // lo: lower corner of the atmospheric boxes, as on a refined level (fg.atm_lo);
    // periodic: the region spans a domain periodic in x and y
    explicit Sampler (IntVect lo = IntVect(0, 0, 0), bool per = false)
        : periodic(per)
    {
        ba_atm    = BoxArray(Box(lo, lo + IntVect(nx - 1, ny - 1, nz - 1)));
        ba_atm_2d = BoxArray(Box(lo, lo + IntVect(nx - 1, ny - 1, 0)));
        dm = DistributionMapping(ba_atm);

        fg.C  = C;
        fg.atm_lo = lo;
        fg.ba = BoxArray(Box(IntVect(0, 0, 0), IntVect(C * nx - 1, C * ny - 1, 0)));
        fg.dm = dm;
        {
            const Box fire_dom(IntVect(0, 0, 0), IntVect(C * nx - 1, C * ny - 1, 0));
            const RealBox rb({AMREX_D_DECL(Real(0.0), Real(0.0), Real(0.0))},
                             {AMREX_D_DECL(Real(C * nx), Real(C * ny), Real(1.0))});
            const int ip = periodic ? 1 : 0;
            fg.geom = Geometry(fire_dom, rb, CoordSys::cartesian, {AMREX_D_DECL(ip, ip, 0)});
        }

        const MFInfo host_info = MFInfo().SetArena(The_Managed_Arena());
        xvel.define(convert(ba_atm, IntVect(1, 0, 0)), dm, 1, 2, host_info);
        yvel.define(convert(ba_atm, IntVect(0, 1, 0)), dm, 1, 2, host_info);
        zcc.define(ba_atm, dm, 1, 2, host_info);
        z0col.define(ba_atm_2d, dm, 1, IntVect(1, 1, 0), host_info);
        fire_wind.define(fg.ba, dm, 2, 0, host_info);
        extract_z.define(fg.ba, dm, 1, 0, host_info);
        surf_z.define(fg.ba, dm, 1, 0, host_info);
        col_ground.define(fg.ba, dm, 4, 0, host_info);

        for (MFIter mfi(xvel); mfi.isValid(); ++mfi) {
            const auto u = xvel.array(mfi);
            const Box gb = mfi.fabbox();
            // Face i between centers i-1 and i: the face average gives u_of(i) at center i
            for (int k = gb.smallEnd(2); k <= gb.bigEnd(2); ++k) {
            for (int j = gb.smallEnd(1); j <= gb.bigEnd(1); ++j) {
            for (int i = gb.smallEnd(0); i <= gb.bigEnd(0); ++i) {
                u(i, j, k) = u_a + u_b * Real(i);
            }}}
        }
        yvel.setVal(v_c);
        for (MFIter mfi(zcc); mfi.isValid(); ++mfi) {
            const auto z = zcc.array(mfi);
            const Box gb = mfi.fabbox();
            for (int k = gb.smallEnd(2); k <= gb.bigEnd(2); ++k) {
            for (int j = gb.smallEnd(1); j <= gb.bigEnd(1); ++j) {
            for (int i = gb.smallEnd(0); i <= gb.bigEnd(0); ++i) {
                z(i, j, k) = (Real(k) + Real(0.5)) * dz;
            }}}
        }
        for (MFIter mfi(z0col); mfi.isValid(); ++mfi) {
            const auto z0 = z0col.array(mfi);
            const Box gb = mfi.fabbox();
            for (int j = gb.smallEnd(1); j <= gb.bigEnd(1); ++j) {
            for (int i = gb.smallEnd(0); i <= gb.bigEnd(0); ++i) {
                z0(i, j, 0) = periodic ? z0_of(z0_column(i, lo[0], nx, true),
                                               z0_column(j, lo[1], ny, true))
                                       : z0_of(i, j);
            }}
        }
        surf_z.setVal(0.0);
        col_ground.setVal(0.0);
    }

    void run (int wind_interp, Real z_sample, bool below_first_log, const MultiFab* z0)
    {
        fill_fire_wind_from_interpolation(fire_wind, extract_z, xvel, yvel, zcc, surf_z, col_ground,
                                          fg, z_ref, nz, nullptr, wind_interp, nullptr, nullptr,
                                          z_sample, z0_scalar, below_first_log, z0);
        Gpu::streamSynchronize();
    }
};

// The closed-form factor of column (i, j) at z_ref, with the roughness of the column it is read from
Real column_factor (int i, int j, Real z_sample, const IntVect& lo, bool periodic)
{
    const Real z0 = z0_of(z0_column(i, lo[0], nx, periodic), z0_column(j, lo[1], ny, periodic));
    if (z_sample > Real(0.0)) { return std::log(z_ref / z0) / std::log(z_sample / z0); }
    return std::log(z_ref / z0) / std::log(Real(0.5) * dz / z0);
}

// Expected wind of fire cell (i_f, j_f): nearest column, or the bilinear blend of four;
// lo is the atmospheric corner of the fire grid
void expected (int i_f, int j_f, int wind_interp, Real z_sample, const IntVect& lo, bool periodic,
               Real& u, Real& v)
{
    if (wind_interp == 0) {
        const int i = i_f / C + lo[0];
        const int j = j_f / C + lo[1];
        const Real f = column_factor(i, j, z_sample, lo, periodic);
        u = f * u_of(i);
        v = f * v_c;
        return;
    }
    const Real gx = (Real(i_f) + Real(0.5)) / Real(C) - Real(0.5);
    const Real gy = (Real(j_f) + Real(0.5)) / Real(C) - Real(0.5);
    const int i0l = static_cast<int>(std::floor(gx));
    const int j0l = static_cast<int>(std::floor(gy));
    const Real wx = gx - Real(i0l);
    const Real wy = gy - Real(j0l);
    const int i0 = i0l + lo[0];
    const int j0 = j0l + lo[1];
    u = Real(0.0);
    v = Real(0.0);
    for (int dj = 0; dj <= 1; ++dj) {
        for (int di = 0; di <= 1; ++di) {
            const Real w = (di ? wx : Real(1.0) - wx) * (dj ? wy : Real(1.0) - wy);
            const Real f = column_factor(i0 + di, j0 + dj, z_sample, lo, periodic);
            u += w * f * u_of(i0 + di);
            v += w * f * v_c;
        }
    }
}

void check_against_closed_form (Sampler& s, int wind_interp, Real z_sample)
{
    const IntVect lo = s.fg.atm_lo;
    int ncells = 0;
    for (MFIter mfi(s.fire_wind); mfi.isValid(); ++mfi) {
        const auto fw = s.fire_wind.const_array(mfi);
        const Box bx = mfi.validbox();
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
        for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
            Real u_ex, v_ex;
            expected(i, j, wind_interp, z_sample, lo, s.periodic, u_ex, v_ex);
            EXPECT_NEAR(fw(i, j, 0, 0), u_ex, rel_tol * std::abs(u_ex))
                << "u, fire cell (" << i << ", " << j << ")";
            EXPECT_NEAR(fw(i, j, 0, 1), v_ex, rel_tol * std::abs(v_ex))
                << "v, fire cell (" << i << ", " << j << ")";
            ++ncells;
        }}
    }
    ParallelDescriptor::ReduceIntSum(ncells);
    EXPECT_EQ(ncells, C * nx * C * ny);
}

} // namespace

// Below the first center, nearest column: each fire cell takes its column's own log profile
TEST(FireWindColumnZ0, BelowFirstCenterNearestUsesColumnRoughness)
{
    Sampler s;
    s.run(0, Real(0.0), true, &s.z0col);
    check_against_closed_form(s, 0, Real(0.0));
}

// Below the first center, bilinear: the four columns are each scaled by their own roughness
TEST(FireWindColumnZ0, BelowFirstCenterBilinearUsesColumnRoughness)
{
    Sampler s;
    s.run(1, Real(0.0), true, &s.z0col);
    check_against_closed_form(s, 1, Real(0.0));
}

// Resampling from 30 m: the factor ln(z_ref / z0) / ln(z_sample / z0) is taken per column
TEST(FireWindColumnZ0, ResampleUsesColumnRoughness)
{
    Sampler s;
    s.run(1, Real(30.0), false, &s.z0col);
    check_against_closed_form(s, 1, Real(30.0));
    s.run(0, Real(30.0), false, &s.z0col);
    check_against_closed_form(s, 0, Real(30.0));
}

// Without col_z0 the scalar roughness applies everywhere; a column field equal to it everywhere
// gives the same wind (to rounding on the resampled path, where the factor moves inside the blend)
TEST(FireWindColumnZ0, UniformColumnFieldMatchesScalar)
{
    for (const Real z_sample : {Real(0.0), Real(30.0)}) {
        for (const int interp : {0, 1}) {
            Sampler s;
            s.run(interp, z_sample, true, nullptr);
            MultiFab scalar_wind(s.fire_wind.boxArray(), s.dm, 2, 0);
            MultiFab::Copy(scalar_wind, s.fire_wind, 0, 0, 2, 0);

            s.z0col.setVal(z0_scalar);
            s.run(interp, z_sample, true, &s.z0col);
            MultiFab::Subtract(scalar_wind, s.fire_wind, 0, 0, 2, 0);
            const Real tol = (z_sample > Real(0.0)) ? rel_tol * Real(10.0) : Real(0.0);
            EXPECT_LE(scalar_wind.norm0(0), tol) << "interp " << interp << ", z_sample " << z_sample;
            EXPECT_LE(scalar_wind.norm0(1), tol) << "interp " << interp << ", z_sample " << z_sample;
        }
    }
}

// With the clamp below the first center (no log profile) and no resampling the column field is
// not used: the wind is the column wind as it is
TEST(FireWindColumnZ0, ClampIgnoresRoughness)
{
    Sampler s;
    s.run(0, Real(0.0), false, &s.z0col);
    for (MFIter mfi(s.fire_wind); mfi.isValid(); ++mfi) {
        const auto fw = s.fire_wind.const_array(mfi);
        const Box bx = mfi.validbox();
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
        for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
            EXPECT_NEAR(fw(i, j, 0, 0), u_of(i / C), rel_tol * u_of(i / C));
            EXPECT_EQ(fw(i, j, 0, 1), v_c);
        }}
    }
}

// On a refined level the fire grid starts at atmospheric cell atm_lo: the roughness must come
// from the atmospheric column (fire index / C + atm_lo), not from the fire-local index
TEST(FireWindColumnZ0, RefinedLevelOffsetUsesAtmosphericColumn)
{
    Sampler s(IntVect(2, 1, 0));
    for (const int interp : {0, 1}) {
        s.run(interp, Real(0.0), true, &s.z0col);
        check_against_closed_form(s, interp, Real(0.0));
        s.run(interp, Real(30.0), false, &s.z0col);
        check_against_closed_form(s, interp, Real(30.0));
    }
}

// Sampling at 8 m, below the 10 m first center, and resampling to z_ref: the two factors of each
// column combine to ln(z_ref / z0) / ln(z1 / z0) with that column's z0
TEST(FireWindColumnZ0, ResampleBelowFirstCenterCombines)
{
    Sampler s;
    for (const int interp : {0, 1}) {
        s.run(interp, Real(8.0), true, &s.z0col);
        check_against_closed_form(s, interp, Real(0.0));
    }
}

// The range check: positive, and e z0 no higher than the lowest height of the log profiles
TEST(FireWindColumnZ0, RangeCheck)
{
    Sampler s;
    // Largest valid z0 is 0.82 m (the ghost ring is not checked): e * 0.82 = 2.229 m
    EXPECT_TRUE(fire_column_z0_range_error(s.z0col, z_ref).empty());
    EXPECT_TRUE(fire_column_z0_range_error(s.z0col, Real(2.24)).empty());
    EXPECT_FALSE(fire_column_z0_range_error(s.z0col, Real(2.22)).empty());

    // A zero roughness in one column
    for (MFIter mfi(s.z0col); mfi.isValid(); ++mfi) {
        if (mfi.validbox().contains(IntVect(1, 2, 0))) { s.z0col[mfi](IntVect(1, 2, 0)) = Real(0.0); }
    }
    EXPECT_FALSE(fire_column_z0_range_error(s.z0col, z_ref).empty());
}

// Non-periodic edges: the outer ghost ring of the roughness is left unfilled (NaN), as
// FillBoundary leaves it outside a non-periodic face. Bilinear edge cells reach that ring;
// the sampler must take the edge column's roughness and never read the ring
TEST(FireWindColumnZ0, NonPeriodicEdgeIgnoresGhostRing)
{
    for (const IntVect lo : {IntVect(0, 0, 0), IntVect(2, 1, 0)}) {
        Sampler s(lo);
        const Real qnan = std::numeric_limits<Real>::quiet_NaN();
        for (MFIter mfi(s.z0col); mfi.isValid(); ++mfi) {
            const auto z0 = s.z0col.array(mfi);
            const Box gb = mfi.fabbox();
            const Box vb = mfi.validbox();
            for (int j = gb.smallEnd(1); j <= gb.bigEnd(1); ++j) {
            for (int i = gb.smallEnd(0); i <= gb.bigEnd(0); ++i) {
                if (!vb.contains(IntVect(i, j, 0))) { z0(i, j, 0) = qnan; }
            }}
        }
        for (const Real z_sample : {Real(0.0), Real(30.0)}) {
            s.run(1, z_sample, z_sample == Real(0.0), &s.z0col);
            check_against_closed_form(s, 1, z_sample);
        }
    }
}

// Periodic in x and y: the ghost ring holds the periodic images (as FillBoundary fills it),
// and an edge column's neighbour across the face uses the image's roughness
TEST(FireWindColumnZ0, PeriodicEdgeUsesImageRoughness)
{
    Sampler s(IntVect(0, 0, 0), true);
    for (const Real z_sample : {Real(0.0), Real(30.0)}) {
        s.run(1, z_sample, z_sample == Real(0.0), &s.z0col);
        check_against_closed_form(s, 1, z_sample);
    }
}

// No log profile (clamp below the first center, no resampling): the answer does not depend on
// the roughness, so a field that fails the range check must not abort the run
TEST(FireWindColumnZ0, RangeCheckSkippedWithoutLogProfile)
{
    Sampler s;
    s.z0col.setVal(Real(50.0));
    ASSERT_FALSE(fire_column_z0_range_error(s.z0col, z_ref).empty());
    s.run(0, Real(0.0), false, &s.z0col);
    for (MFIter mfi(s.fire_wind); mfi.isValid(); ++mfi) {
        const auto fw = s.fire_wind.const_array(mfi);
        const Box bx = mfi.validbox();
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
        for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
            EXPECT_NEAR(fw(i, j, 0, 0), u_of(i / C), rel_tol * u_of(i / C));
        }}
    }
}
