// Unit tests for the Axell & Liungman (2001) one-equation RANS closure
// relations in ERF_RANSClosure.H: the same inline functions that
// ComputeTurbulentViscosityRANS evaluates on the device.

#include <ERF_RANSClosure.H>
#include <ERF_EddyViscosity.H>

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <vector>

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

TEST(RANSClosure, RtMinLowerBoundSitsJustAboveThePoleOfCmu)
{
    // The pole of cmu is the larger root of 1 + 0.308 Rt + 0.00837 Rt^2,
    // -3.5986889850 (computed here in double, independently of the header).
    // The bound must lie above it, and not far above.
    const double b = 0.308, c = 0.00837;
    const double pole = (-b + std::sqrt(b*b - 4.0*c)) / (2.0*c);
    EXPECT_GT(double(AL01::Rt_min_lower_bound), pole);
    EXPECT_LT(double(AL01::Rt_min_lower_bound) - pole, 1e-6);
    // the pole of cmu' (-1/0.277) lies below, so the bound covers both
    EXPECT_LT(-1.0 / 0.277, pole);

    // With the smallest Rt_min the check accepts, the smoothed Rt keeps cmu and
    // cmu' finite and positive however unstable the air. Near the pole the
    // denominator of cmu is about 4e-8, below float round-off, so this part
    // runs in double precision only; -3.598 (cmu about 980) is checked in both.
    const Real samples[] = {Real(-1.5), Real(-10), Real(-1e4), Real(-1e8), Real(-1e30)};
    std::vector<Real> rt_mins = {Real(-3.598)};
    if constexpr (sizeof(Real) == sizeof(double)) {
        rt_mins.push_back(std::nextafter(AL01::Rt_min_lower_bound, Real(0)));
    }
    for (Real rt_min : rt_mins) {
        for (Real Rt : samples) {
            const Real s = AL01::smooth_Rt(Rt, Rt_crit, rt_min);
            const Real cm = AL01::cmu(s, Cmu0);
            const Real cp = AL01::cmu_prime(s, Cmu0);
            ASSERT_TRUE(std::isfinite(cm)) << "Rt_min=" << rt_min << " Rt=" << Rt;
            EXPECT_GT(cm, Real(0)) << "Rt_min=" << rt_min << " Rt=" << Rt;
            EXPECT_GT(cp, Real(0)) << "Rt_min=" << rt_min << " Rt=" << Rt;
        }
    }

    // The old bound (-3.6) accepted Rt_min = -3.599, below the pole: a strongly
    // unstable cell (raw Rt below about -2.2e4) then gave a negative cmu.
    if constexpr (sizeof(Real) == sizeof(double)) {
        const Real s = AL01::smooth_Rt(Real(-1e8), Rt_crit, Real(-3.599));
        EXPECT_LT(AL01::cmu(s, Cmu0), Real(0));
        EXPECT_LT(Real(-3.599), AL01::Rt_min_lower_bound);
    }
}

TEST(RANSClosure, RtMinAboveTheZeroOfTheCmuNumerator)
{
    // cmu = 0 where Cmu0 + 0.108 Rt = 0; below that (and above the pole) it is
    // negative, which the pole bound alone does not exclude when Cmu0 < 0.3887.
    EXPECT_NEAR(AL01::cmu(AL01::cmu_numerator_root(Cmu0), Cmu0), Real(0), tol());
    EXPECT_NEAR(AL01::cmu_numerator_root(Cmu0), Real(-5.15), Real(1e-4));
    // the default Cmu0: the pole is the binding limit
    EXPECT_LT(AL01::cmu_numerator_root(Cmu0), AL01::Rt_min_lower_bound);

    // Cmu0 = 0.3 with the default Rt_min = -3: the numerator's zero (-2.78)
    // lies above Rt_min, so the start-up check must reject the pair, and a
    // strongly unstable cell would otherwise get a negative cmu
    const Real cmu0_low = Real(0.3);
    EXPECT_GT(AL01::cmu_numerator_root(cmu0_low), Rt_min);
    const Real s = AL01::smooth_Rt(Real(-1e4), Rt_crit, Rt_min);
    EXPECT_LT(AL01::cmu(s, cmu0_low), Real(0));

    // just above the zero, cmu and cmu' are positive for every smoothed Rt
    const Real rt_min_ok = AL01::cmu_numerator_root(cmu0_low) + Real(1e-3);
    const Real samples[] = {Real(-1.5), Real(-10), Real(-1e4), Real(-1e8), Real(-1e30)};
    for (Real Rt : samples) {
        const Real t = AL01::smooth_Rt(Rt, Real(-1), rt_min_ok);
        EXPECT_GT(AL01::cmu(t, cmu0_low), Real(0)) << "Rt=" << Rt;
        EXPECT_GT(AL01::cmu_prime(t, cmu0_low), Real(0)) << "Rt=" << Rt;
    }
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

// Realizability cap nu_t <= C_r k / S (erf.rans_realizability)
TEST(RANSClosure, RealizabilityCapMatchesClosedForm)
{
    const Real Cr = Cmu0 * Cmu0;                    // 0.3094, the equilibrium value
    const Real tke = Real(0.5), S = Real(0.2);
    const Real cap = Cr * tke / S;                  // 0.7735
    Real ratio = Real(-1);

    // above the cap: capped to C_r k / S, ratio = nut / cap > 1
    const Real nut_hi = Real(2) * cap;
    EXPECT_NEAR(AL01::realizable_nut(nut_hi, tke, S, Cr, ratio), cap, tol() * cap);
    EXPECT_NEAR(ratio, Real(2), tol());

    // below the cap: unchanged, ratio < 1
    const Real nut_lo = Real(0.25) * cap;
    EXPECT_EQ(AL01::realizable_nut(nut_lo, tke, S, Cr, ratio), nut_lo);
    EXPECT_NEAR(ratio, Real(0.25), tol());

    // S enters linearly (not S^2): doubling S halves the cap
    EXPECT_NEAR(AL01::realizable_nut(nut_hi, tke, Real(2) * S, Cr, ratio), Real(0.5) * cap, tol() * cap);
    EXPECT_NEAR(ratio, Real(4), tol());
}

TEST(RANSClosure, RealizabilityCapIsNeutralInEquilibriumShear)
{
    // Constant-stress layer: k = u*^2 / Cmu0^2, S = u* / (kappa z), nu_t = kappa u* z
    // (AL01 Eqs. 11, 16 and 19 with P = eps), so nu_t S / k = Cmu0^2 exactly and the
    // cap with C_r = Cmu0^2 returns nu_t unchanged (ratio = 1).
    const Real ustar = Real(0.3), z = Real(10);
    const Real tke = ustar * ustar / (Cmu0 * Cmu0);
    const Real S   = ustar / (kappa * z);
    // nu_t from the closure functions themselves (neutral: Rt = 0, l = kappa z)
    const Real nut = AL01::cmu(Real(0), Cmu0) * std::sqrt(tke) * (kappa * z);
    EXPECT_NEAR(nut, kappa * ustar * z, tol() * nut);
    // and P = eps there: nu_t S^2 against dissipation() per unit mass
    EXPECT_NEAR(nut * S * S, AL01::dissipation(Real(1), Cmu0_pow3, tke, kappa * z), tol() * nut * S * S);
    Real ratio = Real(0);
    EXPECT_NEAR(AL01::realizable_nut(nut, tke, S, Cmu0 * Cmu0, ratio), nut, tol() * nut);
    EXPECT_NEAR(ratio, Real(1), tol());
    // a slightly larger C_r leaves it alone, a slightly smaller one trims it
    EXPECT_EQ(AL01::realizable_nut(nut, tke, S, Real(1.01) * Cmu0 * Cmu0, ratio), nut);
    EXPECT_LT(AL01::realizable_nut(nut, tke, S, Real(0.99) * Cmu0 * Cmu0, ratio), nut);
}

TEST(RANSClosure, RealizabilityCapZeroStrainAndTinyK)
{
    // (C_r = 0, the cap off, never reaches the function: the kernel skips the cap
    // and the strain read, and writes cap_ratio = 0 itself)
    const Real tke = Real(0.5), nut = Real(3);
    Real ratio = Real(7);
    // S = 0 (uniform flow): no cap and a finite ratio, no division by zero
    Real out = AL01::realizable_nut(nut, tke, Real(0), Real(0.309), ratio);
    EXPECT_EQ(out, nut);
    ASSERT_TRUE(std::isfinite(ratio));
    EXPECT_GE(ratio, Real(0));
    EXPECT_LT(ratio, Real(1));
    // k at a tiny floor with finite S: the cap goes to zero with k
    out = AL01::realizable_nut(nut, Real(1e-16), Real(1), Real(0.309), ratio);
    EXPECT_NEAR(out, Real(0.309e-16), Real(1e-20));
}

// The strain magnitude the cap divides by: S = sqrt(2 ComputeSmnSmn), the
// cell-centred average of the face/edge strains the kernel reads.
namespace {

// Fill the six strain arrays of a one-cell box from the given uniform values
// (tau12/13/23 are nodal on their edges) and return sqrt(2 SmnSmn) at the cell.
Real strain_magnitude (Real t11, Real t22, Real t33, Real t12, Real t13, Real t23,
                       int node13 = -1)
{
    using namespace amrex;
    const Box cc(IntVect(0), IntVect(0));
    FArrayBox f11(cc, 1, The_Pinned_Arena()), f22(cc, 1, The_Pinned_Arena()), f33(cc, 1, The_Pinned_Arena());
    FArrayBox f12(convert(cc, IntVect(1,1,0)), 1, The_Pinned_Arena());
    FArrayBox f13(convert(cc, IntVect(1,0,1)), 1, The_Pinned_Arena());
    FArrayBox f23(convert(cc, IntVect(0,1,1)), 1, The_Pinned_Arena());
    f11.setVal<RunOn::Host>(t11); f22.setVal<RunOn::Host>(t22); f33.setVal<RunOn::Host>(t33);
    f12.setVal<RunOn::Host>(t12); f13.setVal<RunOn::Host>(t13); f23.setVal<RunOn::Host>(t23);
    if (node13 >= 0) {
        // only one of the four tau13 edges carries the value: the average is a quarter
        f13.setVal<RunOn::Host>(Real(0));
        const IntVect nodes[4] = {IntVect(0,0,0), IntVect(1,0,0), IntVect(0,0,1), IntVect(1,0,1)};
        f13.array()(nodes[node13][0], nodes[node13][1], nodes[node13][2]) = t13;
    }
    FArrayBox out(cc, 1, The_Pinned_Arena());
    const Array4<Real const> a11 = f11.const_array(), a22 = f22.const_array(), a33 = f33.const_array();
    const Array4<Real const> a12 = f12.const_array(), a13 = f13.const_array(), a23 = f23.const_array();
    const Array4<Real> o = out.array();
    ParallelFor(cc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        o(i,j,k) = std::sqrt(Real(2) * ComputeSmnSmn(i,j,k,a11,a22,a33,a12,a13,a23));
    });
    Gpu::streamSynchronize();
    return out.array()(0,0,0);
}

} // namespace

TEST(RANSClosure, StrainMagnitudeForTheCap)
{
    // simple shear u = a z: S_13 = a/2 on every edge, every other component 0, S = |a|
    const Real a = Real(0.8);
    EXPECT_NEAR(strain_magnitude(0, 0, 0, 0, Real(0.5) * a, 0), a, tol());
    // plane strain S_11 = b, S_33 = -b: 2 S_ij S_ij = 4 b^2, S = 2 |b|
    const Real b = Real(0.3);
    EXPECT_NEAR(strain_magnitude(b, 0, -b, 0, 0, 0), Real(2) * b, tol());
    // the edge strains are averaged over four edges: one edge of four at a/2 gives a/4
    EXPECT_NEAR(strain_magnitude(0, 0, 0, 0, Real(0.5) * a, 0, 2), Real(0.25) * a, tol());
    // all three edge strains at c: S_ij S_ij = 2 (3 c^2), S = sqrt(12) c
    const Real c = Real(0.2);
    EXPECT_NEAR(strain_magnitude(0, 0, 0, c, c, c), c * std::sqrt(Real(12)), tol());
}
