// The 10 m wind diagnostics: the speed of the surface-layer momentum profile at 10 m, in
// the direction of the first-cell wind, turned earth-relative with the grid orientation.
#include <cmath>

#include <gtest/gtest.h>

#include <AMReX_MultiFab.H>

#include "Diagnostics/ERF_NearSurfaceDiagnostics.H"

namespace
{

using amrex::Real;
using near_surface_diagnostics::WindSources;
using near_surface_diagnostics::fill_wind_10m;

constexpr Real missing = Real(-999.0);

// One column of two cells in x on a three-cell-deep grid, with its 2D surface fields
struct WindCase
{
    amrex::BoxArray ba{amrex::Box(amrex::IntVect(0,0,0), amrex::IntVect(1,0,2))};
    amrex::BoxArray ba2d{amrex::Box(amrex::IntVect(0,0,0), amrex::IntVect(1,0,0))};
    amrex::DistributionMapping dm{ba};
    amrex::MultiFab xvel{amrex::convert(ba, amrex::IntVect(1,0,0)), dm, 1, 0};
    amrex::MultiFab yvel{amrex::convert(ba, amrex::IntVect(0,1,0)), dm, 1, 0};
    amrex::MultiFab ustar{ba2d, dm, 1, 0};
    amrex::MultiFab z0{ba2d, dm, 1, 0};
    amrex::MultiFab olen{ba2d, dm, 1, 0};
    amrex::MultiFab cosa{ba2d, dm, 1, 0};
    amrex::MultiFab sina{ba2d, dm, 1, 0};
    amrex::MultiFab dst{ba2d, dm, 4, 0};

    WindCase (Real u, Real v, Real u_star, Real roughness, Real obukhov_length)
    {
        xvel.setVal(u);
        yvel.setVal(v);
        ustar.setVal(u_star);
        z0.setVal(roughness);
        olen.setVal(obukhov_length);
        cosa.setVal(Real(1.0));
        sina.setVal(Real(0.0));
        dst.setVal(Real(0.0));
    }

    WindSources sources (bool with_rotation = false)
    {
        WindSources s;
        s.u_star = &ustar;
        s.roughness_height = &z0;
        s.obukhov_length = &olen;
        s.xvel = &xvel;
        s.yvel = &yvel;
        if (with_rotation) {
            s.cos_alpha = &cosa;
            s.sin_alpha = &sina;
        }
        s.klo = 0;
        return s;
    }

    Real value (int comp) { return dst.min(comp); }
};

} // namespace

// Motivation: in near-neutral conditions the 10 m speed is u_star/kappa ln(10/z0), and
// it must point along the first-cell wind.
TEST(NearSurfaceWind, NeutralSpeedAlongTheFirstCellWind)
{
    const Real u_star = Real(0.4);
    const Real z0 = Real(0.1);
    WindCase c(Real(3.0), Real(4.0), u_star, z0, Real(1.0e10));
    fill_wind_10m(c.dst, 0, 1, 2, 3, c.sources());

    const Real speed = u_star / KAPPA * std::log(Real(10.0) / z0);
    const Real tol = Real(1.0e-5) * speed;
    EXPECT_NEAR(c.value(0), Real(0.6) * speed, tol);
    EXPECT_NEAR(c.value(1), Real(0.8) * speed, tol);
    // Without a grid orientation the earth-relative winds are the grid-relative ones
    EXPECT_NEAR(c.value(2), Real(0.6) * speed, tol);
    EXPECT_NEAR(c.value(3), Real(0.8) * speed, tol);
}

// Motivation: the earth-relative winds turn the grid-relative ones by the grid
// orientation, here a grid whose x axis points north.
TEST(NearSurfaceWind, EarthRelativeWindsTurnWithTheGrid)
{
    WindCase c(Real(1.0), Real(0.0), Real(0.3), Real(0.05), Real(1.0e10));
    c.cosa.setVal(Real(0.0));
    c.sina.setVal(Real(1.0));
    fill_wind_10m(c.dst, 0, 1, 2, 3, c.sources(true));

    const Real speed = Real(0.3) / KAPPA * std::log(Real(10.0) / Real(0.05));
    const Real tol = Real(1.0e-5) * speed;
    EXPECT_NEAR(c.value(0), speed, tol);
    EXPECT_NEAR(c.value(1), Real(0.0), tol);
    EXPECT_NEAR(c.value(2), Real(0.0), tol);  // u_east
    EXPECT_NEAR(c.value(3), speed, tol);      // v_north
}

// Motivation: the stability correction enters through psi_m, with 10/L capped at the
// same bound as the 2 m profile.
TEST(NearSurfaceWind, StableProfileUsesTheMomentumCorrection)
{
    const similarity_funs sfuns{};
    const Real u_star = Real(0.2);
    const Real z0 = Real(0.05);
    const Real obukhov_length = Real(50.0);
    WindCase c(Real(2.0), Real(0.0), u_star, z0, obukhov_length);
    fill_wind_10m(c.dst, 0, -1, -1, -1, c.sources());

    const Real speed = u_star / KAPPA *
        (std::log(Real(10.0) / z0) - sfuns.calc_psi_m(Real(10.0) / obukhov_length));
    EXPECT_NEAR(c.value(0), speed, Real(1.0e-5) * speed);
}

// Motivation: with no surface layer, or a roughness height at or above 10 m, there is
// no 10 m wind to give, and every requested component is the missing value.
TEST(NearSurfaceWind, MissingWithoutAnEvaluableProfile)
{
    {
        WindCase c(Real(3.0), Real(4.0), Real(0.4), Real(12.0), Real(1.0e10));
        fill_wind_10m(c.dst, 0, 1, 2, 3, c.sources());
        for (int n = 0; n < 4; ++n) { EXPECT_EQ(c.value(n), missing); }
    }
    {
        WindCase c(Real(3.0), Real(4.0), Real(0.4), Real(0.1), Real(1.0e10));
        WindSources s = c.sources();
        s.u_star = nullptr;
        fill_wind_10m(c.dst, 0, 1, 2, 3, s);
        for (int n = 0; n < 4; ++n) { EXPECT_EQ(c.value(n), missing); }
    }
}
