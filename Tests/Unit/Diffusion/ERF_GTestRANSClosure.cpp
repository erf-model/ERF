// Unit tests for the Axell & Liungman (2001) one-equation RANS closure
// relations in ERF_RANSClosure.H: the same inline functions that
// ComputeTurbulentViscosityRANS evaluates on the device.

#include <ERF_RANSClosure.H>

#include <gtest/gtest.h>

#include <cmath>
#include <limits>

using amrex::Real;

namespace {

constexpr Real Cmu0    = Real(0.5562);
constexpr Real Cb      = Real(0.35);
constexpr Real Rt_crit = Real(-1.0);
constexpr Real Rt_min  = Real(-3.0);
constexpr Real kappa   = Real(0.41);

constexpr Real Cmu0_pow3 = Cmu0 * Cmu0 * Cmu0;
constexpr Real inv_Cb_sq = Real(1.0) / (Cb * Cb);

Real tol (Real factor = Real(1e3))
{
    return factor * std::numeric_limits<Real>::epsilon();
}

} // namespace

TEST(RANSClosure, NeutralStabilityFunctionsReduceToCmu0)
{
    EXPECT_NEAR(AL01::cmu(Real(0), Cmu0), Cmu0, tol());
    EXPECT_NEAR(AL01::cmu_prime(Real(0), Cmu0), Cmu0, tol());
}

TEST(RANSClosure, StabilityFunctionsMonotoneAndFiniteOnAdmissibleRange)
{
    // From the strongest smoothed unstable value up to the stable limit
    // Cb^2 / Cmu0^6 (about 4.1) the functions are finite, positive and cmu'
    // decreases with Rt (more stable air mixes scalars less).
    const Real Rt_stable_max = Cb * Cb / (Cmu0_pow3 * Cmu0_pow3);
    Real prev_cmu_prime = std::numeric_limits<Real>::max();
    for (int i = 0; i <= 400; ++i) {
        Real Rt = Rt_min + (Rt_stable_max - Rt_min) * Real(i) / Real(400);
        Real c  = AL01::cmu(Rt, Cmu0);
        Real cp = AL01::cmu_prime(Rt, Cmu0);
        ASSERT_TRUE(std::isfinite(c))  << "cmu at Rt=" << Rt;
        ASSERT_TRUE(std::isfinite(cp)) << "cmu' at Rt=" << Rt;
        EXPECT_GT(c,  Real(0)) << "Rt=" << Rt;
        EXPECT_GT(cp, Real(0)) << "Rt=" << Rt;
        EXPECT_LE(cp, prev_cmu_prime + tol()) << "cmu' not monotone at Rt=" << Rt;
        prev_cmu_prime = cp;
    }
}

TEST(RANSClosure, SmoothingIsIdentityAboveRtCritAndContinuousAtIt)
{
    EXPECT_EQ(AL01::smooth_Rt(Real(0.5), Rt_crit, Rt_min), Real(0.5));
    EXPECT_EQ(AL01::smooth_Rt(Rt_crit, Rt_crit, Rt_min), Rt_crit);
    const Real below = Rt_crit - Real(1e-6);
    EXPECT_NEAR(AL01::smooth_Rt(below, Rt_crit, Rt_min), Rt_crit, Real(1e-5));
}

TEST(RANSClosure, SmoothingIsBoundedBelowByRtMinForAnyRt)
{
    // The mapped value must stay in (Rt_min, Rt_crit] for every unstable Rt,
    // including the enormous magnitudes reached when k sits at its floor
    // under a strong unstable N^2, where the formula cancels catastrophically.
    const Real samples[] = {Real(-1.5), Real(-3), Real(-10), Real(-1e3),
                            Real(-1e8), Real(-1e16), Real(-1e30)};
    Real prev = Rt_crit;
    for (Real Rt : samples) {
        Real s = AL01::smooth_Rt(Rt, Rt_crit, Rt_min);
        ASSERT_TRUE(std::isfinite(s)) << "Rt=" << Rt;
        EXPECT_GE(s, Rt_min) << "Rt=" << Rt;
        EXPECT_LE(s, Rt_crit) << "Rt=" << Rt;
        EXPECT_LE(s, prev + tol()) << "smoothing not monotone at Rt=" << Rt;
        prev = s;
        // the stability functions stay finite and positive on the mapped value
        EXPECT_GT(AL01::cmu(s, Cmu0), Real(0)) << "Rt=" << Rt;
        EXPECT_GT(AL01::cmu_prime(s, Cmu0), Real(0)) << "Rt=" << Rt;
    }
    EXPECT_GT(Rt_min, AL01::Rt_min_lower_bound);
}

TEST(RANSClosure, GeometricLengthCapIsHarmonic)
{
    const Real l_g_max = Real(30);
    // far below the cap the length is kappa (z + z0) to first order
    Real small = kappa * Real(1.0);
    EXPECT_NEAR(AL01::geom_length(small, l_g_max), small, Real(0.02) * small);
    // the cap is never exceeded and is approached from below
    EXPECT_LT(AL01::geom_length(Real(1e6), l_g_max), l_g_max);
    EXPECT_GT(AL01::geom_length(Real(1e6), l_g_max), Real(0.999) * l_g_max);
    // equal arguments give half the cap
    EXPECT_NEAR(AL01::geom_length(l_g_max, l_g_max), Real(0.5) * l_g_max, tol());
}

TEST(RANSClosure, TurbulentLengthLimits)
{
    const Real eps = std::numeric_limits<Real>::epsilon();
    const Real l_g = Real(10);
    const Real tke = Real(0.5);

    auto length = [&](Real N2, Real k) {
        return AL01::turb_length(l_g, N2, k, Cmu0_pow3, inv_Cb_sq, Rt_crit, Rt_min, eps);
    };

    // neutral: the geometric length
    EXPECT_EQ(length(Real(0), tke), l_g);

    // stable: shorter than l_g and approaching Cb sqrt(k / N^2) for strong N^2
    const Real N2_weak = Real(1e-5);
    EXPECT_LT(length(N2_weak, tke), l_g);
    const Real N2_strong = Real(1.0);
    Real l_strong = length(N2_strong, tke);
    EXPECT_NEAR(l_strong, Cb * std::sqrt(tke / N2_strong), Real(0.01) * l_strong);

    // unstable: longer than l_g, finite, and never above the bound set by
    // Rt_min, however strong the convection or small the TKE
    const Real bound = AL01::unstable_length_bound(l_g, Cmu0_pow3, inv_Cb_sq, Rt_min);
    EXPECT_NEAR(bound / l_g, std::sqrt(Real(1) + Cmu0_pow3 * Cmu0_pow3 * inv_Cb_sq * Real(3)), tol());
    Real prev = l_g;
    for (Real N2 : {-N2_weak, Real(-1e-3), Real(-1.0), Real(-1e3)}) {
        for (Real k : {tke, Real(1e-6), Real(1e-16)}) {
            Real l = length(N2, k);
            ASSERT_TRUE(std::isfinite(l)) << "N2=" << N2 << " k=" << k;
            EXPECT_GT(l, l_g) << "N2=" << N2 << " k=" << k;
            EXPECT_LE(l, bound * (Real(1) + tol())) << "N2=" << N2 << " k=" << k;
        }
        // stronger convection never shortens the length
        Real l = length(N2, tke);
        EXPECT_GE(l, prev - tol());
        prev = l;
    }
    // moderate convection: Eq. 28 above Rt_crit, unsmoothed
    const Real Rt_mod = Real(-0.5);
    const Real N2_mod = Rt_mod * tke * Cmu0_pow3 * Cmu0_pow3 / (l_g * l_g);
    EXPECT_NEAR(length(N2_mod, tke), l_g * std::sqrt(Real(1) - Cmu0_pow3 * Cmu0_pow3 * inv_Cb_sq * Rt_mod), tol() * l_g);
}

TEST(RANSClosure, DissipationAndRichardsonAreConsistent)
{
    // Rt from richardson() equals k^2 N^2 / eps^2 with eps from dissipation()
    const Real rho = Real(1.1), tke = Real(0.4), length = Real(7), N2 = Real(2e-4);
    Real eps_v = AL01::dissipation(rho, Cmu0_pow3, tke, length) / rho;   // per unit mass
    Real Rt_direct = tke * tke * N2 / (eps_v * eps_v);
    // relative tolerance aware of the build precision (1e-12 double, 1e-5 float)
    const Real rtol = (sizeof(Real) == 8) ? Real(1e-12) : Real(1e-5);
    EXPECT_NEAR(AL01::richardson(length, N2, tke, Cmu0_pow3), Rt_direct, rtol * std::abs(Rt_direct));
}
