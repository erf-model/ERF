// The stability coordinate z/L used by the 2 m MOST profile is bounded on the
// stable side, so a very small positive Obukhov length cannot drive the
// diagnostic to an unphysical value.
#include <cmath>
#include <limits>

#include <gtest/gtest.h>

#include <AMReX_FArrayBox.H>

#include "Diagnostics/ERF_NearSurfaceDiagnostics.H"
#include "ERF_MOSTStress.H"

namespace
{

using amrex::Real;
using near_surface_diagnostics::MostProfile;
using near_surface_diagnostics::MostState;
using near_surface_diagnostics::evaluate_most_factor;
using near_surface_diagnostics::evaluate_most_profile;
using near_surface_diagnostics::most_profile_zeta_max;

MostState state_with (Real roughness, Real obukhov_length)
{
    return MostState{Real(270.0), Real(0.56), Real(0.002), Real(0.0),
                     roughness, obukhov_length};
}

} // namespace

// Motivation: with z0 at its floor and L of a few centimetres, the unbounded
// stable correction gave a 2 m temperature of several hundred K.
TEST(NearSurfaceStability, TinyPositiveObukhovLengthGivesABoundedFactor)
{
    const similarity_funs sfuns{};
    const Real z0 = Real(1.0e-7);
    const Real height = Real(2.0);
    Real factor = Real(0.0);
    ASSERT_TRUE(evaluate_most_factor(state_with(z0, Real(0.04)), height, factor));
    const Real expected = std::log(height / z0) - sfuns.calc_psi_h(most_profile_zeta_max());
    EXPECT_NEAR(factor, expected, Real(1.0e-4));

    MostProfile profile{};
    ASSERT_TRUE(evaluate_most_profile(state_with(z0, Real(0.04)), height, profile));
    EXPECT_LT(profile.potential_temperature, Real(400.0));
}

// Motivation: the bound must not change any cell whose z/L is below it.
TEST(NearSurfaceStability, ModerateStabilityIsUnchanged)
{
    const similarity_funs sfuns{};
    const Real z0 = Real(0.05);
    const Real height = Real(2.0);
    const Real obukhov_length = Real(20.0);
    Real factor = Real(0.0);
    ASSERT_TRUE(evaluate_most_factor(state_with(z0, obukhov_length), height, factor));
    EXPECT_NEAR(factor,
                std::log(height / z0) - sfuns.calc_psi_h(height / obukhov_length),
                Real(1.0e-6));
}

// Motivation: the unstable branch is not bounded here and must stay as it was.
TEST(NearSurfaceStability, UnstableBranchIsUnchanged)
{
    const similarity_funs sfuns{};
    const Real z0 = Real(0.05);
    const Real height = Real(2.0);
    const Real obukhov_length = Real(-5.0);
    Real factor = Real(0.0);
    ASSERT_TRUE(evaluate_most_factor(state_with(z0, obukhov_length), height, factor));
    EXPECT_NEAR(factor,
                std::log(height / z0) - sfuns.calc_psi_h(height / obukhov_length),
                Real(1.0e-6));
}

// Motivation: the factor must not decrease as stability grows, and must reach
// its bound.
TEST(NearSurfaceStability, FactorIsMonotoneAndSaturatesAtTheBound)
{
    const Real z0 = Real(0.01);
    const Real height = Real(2.0);
    Real previous = -std::numeric_limits<Real>::infinity();
    Real saturated = Real(0.0);
    for (const Real obukhov_length : {Real(1000.0), Real(10.0), Real(2.0), Real(1.0),
                                      Real(0.5), Real(0.04)}) {
        Real factor = Real(0.0);
        ASSERT_TRUE(evaluate_most_factor(state_with(z0, obukhov_length), height, factor));
        EXPECT_GE(factor, previous);
        previous = factor;
        saturated = factor;
    }
    Real at_bound = Real(0.0);
    ASSERT_TRUE(evaluate_most_factor(state_with(z0, Real(2.0)), height, at_bound));
    EXPECT_EQ(saturated, at_bound);
}

// erf-model/ERF#4765: with a prescribed upward heat flux, psi_h grows without bound as z/L -> -inf and
// passes ln(zref / z0), so the diagnosed surface temperature falls below the air and the Moeng local
// flux, which divides by theta_mean - theta_surf, is singular. The crossing is at z/L = -8.4 for
// zref / z0 = 40 and -3.9 for zref / z0 = 20, so a fixed bound on z/L cannot cover both.
#ifndef AMREX_USE_GPU
namespace
{

struct FluxResult { Real t_surf; Real u_star; Real t_mean; };

FluxResult prescribed_flux_surface (Real flux, Real wind, Real z0, Real zref)
{
    const amrex::Box bx(amrex::IntVect(0), amrex::IntVect(0));
    amrex::FArrayBox zref_f(bx), z0_f(bx), umm_f(bx), tm_f(bx), qv_f(bx), us_f(bx), ts_f(bx),
                     qs_f(bx), tstar_f(bx), qstar_f(bx), ol_f(bx), pblh_f(bx);
    zref_f.setVal<amrex::RunOn::Host>(zref);
    z0_f.setVal<amrex::RunOn::Host>(z0);
    umm_f.setVal<amrex::RunOn::Host>(wind);
    tm_f.setVal<amrex::RunOn::Host>(Real(300.0));
    qv_f.setVal<amrex::RunOn::Host>(Real(0.0));
    us_f.setVal<amrex::RunOn::Host>(bogus_large_value);
    ts_f.setVal<amrex::RunOn::Host>(Real(300.0));
    qs_f.setVal<amrex::RunOn::Host>(Real(0.0));
    tstar_f.setVal<amrex::RunOn::Host>(Real(0.0));
    qstar_f.setVal<amrex::RunOn::Host>(Real(0.0));
    ol_f.setVal<amrex::RunOn::Host>(Real(0.0));
    pblh_f.setVal<amrex::RunOn::Host>(Real(1000.0));

    surface_flux sf(flux, Real(0.0), false);
    sf.iterate_flux(0, 0, 0, 100,
                    zref_f.const_array(), z0_f.const_array(), umm_f.const_array(),
                    tm_f.const_array(), tm_f.const_array(), qv_f.const_array(),
                    us_f.array(), amrex::Array4<Real>{}, tstar_f.array(), qstar_f.array(),
                    ts_f.array(), qs_f.array(), ol_f.array(), pblh_f.array(),
                    amrex::Array4<Real>{}, amrex::Array4<Real>{}, amrex::Array4<Real>{});
    return FluxResult{ts_f.array()(0,0,0), us_f.array()(0,0,0), Real(300.0)};
}

} // namespace

TEST(MostPrescribedFlux, SurfaceIsWarmerThanAirForAnUpwardFluxAtTheReproducerValues)
{
    // zref / z0 = 40: the case of #4765
    const FluxResult r = prescribed_flux_surface(Real(0.1408), Real(0.1), Real(0.05), Real(2.0));
    EXPECT_GT(r.t_surf, r.t_mean);
    EXPECT_GT(r.u_star, Real(0.0));
}

TEST(MostPrescribedFlux, SurfaceIsWarmerThanAirWhereAFixedZetaBoundWouldNotHold)
{
    // zref / z0 = 20 (default z0 at 2 m): the bracket crossed zero at z/L = -3.9, inside a bound of -5
    const FluxResult r = prescribed_flux_surface(Real(0.1408), Real(0.1), Real(0.1), Real(2.0));
    EXPECT_GT(r.t_surf, r.t_mean);
    EXPECT_GT(r.u_star, Real(0.0));
}

TEST(MostPrescribedFlux, SurfaceIsWarmerThanAirAtARoughSmallRatio)
{
    // zref / z0 = 10 (10 m, 1 m): the bracket crossed zero at z/L = -1.7
    const FluxResult r = prescribed_flux_surface(Real(0.1408), Real(0.1), Real(1.0), Real(10.0));
    EXPECT_GT(r.t_surf, r.t_mean);
    EXPECT_GT(r.u_star, Real(0.0));
}

TEST(MostPrescribedFlux, ModerateInstabilityIsUnchangedByTheCap)
{
    // Convective deck values (zref 10 m, z0 0.1 m, F = 0.24 K m/s, 5 m/s): z/L of order -1, far from the cap
    const FluxResult r = prescribed_flux_surface(Real(0.24), Real(5.0), Real(0.1), Real(10.0));
    EXPECT_GT(r.t_surf, r.t_mean);
    EXPECT_LT(r.t_surf - r.t_mean, Real(5.0));
}

TEST(MostPrescribedFlux, PsiCapKeepsTheBracketAtLeastATenthOfTheLog)
{
    const Real zref = Real(2.0), z0 = Real(0.1);
    const Real C = std::log(zref / z0);
    EXPECT_NEAR(most_capped_psi(Real(100.0), zref, z0), Real(0.9) * C, Real(1e-6));
    EXPECT_EQ(most_capped_psi(Real(1.0), zref, z0), Real(1.0));      // below the cap: bit for bit
    EXPECT_EQ(most_capped_psi(Real(-3.0), zref, z0), Real(-3.0));    // stable side untouched
}
#endif
