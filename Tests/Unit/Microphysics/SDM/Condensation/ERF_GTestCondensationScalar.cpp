#include <gtest/gtest.h>

#include "ERF_GTestCondensationCommon.H"

using namespace condensation_test;

namespace {

struct ConstantRadiusSquaredRate
{
    float rate;

    AMREX_GPU_HOST_DEVICE
    float rhs_func (float, float, float, float, float, float, float, float, float) const noexcept
    {
        return rate;
    }

    AMREX_GPU_HOST_DEVICE
    float rhs_jac (float, float, float, float, float, float, float, float) const noexcept
    {
        return 0;
    }
};

} // namespace

// Motivation: ERF's warm-phase super-droplet growth integrates dR^2/dt with the
// Koehler curvature and solute terms (Source/Particles/ERF_SuperDropletPCMassChange.H,
// SDMassChangeUtils_LV). This checks the diffusional r^2-law, the growth/evaporation
// sign, the Kelvin and Raoult limits, Koehler equilibrium and activation, the
// time-integrated growth, and the ventilation factor -- all on the host.
TEST(CondensationScalar, GrowthPropertiesHoldOnHost)
{
    for (int i = 0; i < kNumMetrics; ++i) {
        EXPECT_LE(condensation_metric(i), kTol) << metric_label(i);
    }
}

// Motivation: A Newton tolerance below float resolution can reject resolvable
// short steps and suppress physical exchange. Independently bracketed implicit
// roots for both signs distinguish that loss from an unresolvably weak fixture.
TEST(CondensationScalar, FloatBackwardEulerRetainsResolvableShortStepExchange)
{
    using RT = float;
    const auto material = water_material();
    const dRsqdt<RT> ode{RT(lcond), RT(therco), RT(R_v), RT(kRhoWater)};
    const RT dt = RT(0.00025);
    const RT initial = RT(1.e-10);
    const RT temperature = RT(kT);
    const RT es = RT(e_sat_water(temperature));
    const RT diffusion = RT(material.coeffMolecularDiffusion(temperature, kP));
    const RT curvature = RT(material.coeffCurv(temperature));
    const RT solute = RT(material.coeffVPSolute());
    const NewtonSolver<dRsqdt<RT>, RT> newton{
        ode, RT(1.e-6), RT(0), RT(1.e-12), 10};
    for (const RT saturation : {RT(0.2), RT(1.5)}) {
        TI<dRsqdt<RT>, NewtonSolver<dRsqdt<RT>, RT>, RT> ti{
            ode, newton, dt, RT(100), saturation, temperature, es, diffusion, curvature,
            solute, RT(0), RT(1000), RT(1.e-40), RT(1.e-3), RT(1.e-6),
            false, false, false, RT(kVentA1), RT(kVentB1), RT(kVentA2), RT(kVentB2), RT(kVentFcap)};
        RT actual = initial;
        bool success = false;
        ti.be(actual, success);
        ASSERT_TRUE(success);

        // Bracket the implicit equation independently of Newton's convergence test.
        const dRsqdt<double> reference_ode{ode.L, ode.K, ode.Rv, ode.rho_l};
        double low = 0.5 * initial, high = 1.5 * initial;
        for (int n = 0; n < 100; ++n) {
            const double mid = 0.5 * (low + high);
            const double residual = mid - initial - dt * reference_ode.rhs_func(
                mid, saturation, temperature, es, diffusion, curvature, solute, 0., 1.);
            if (residual < 0) { low = mid; }
            else { high = mid; }
        }
        const double expected = 0.5 * (low + high);
        const double roundoff = std::numeric_limits<RT>::epsilon() * initial;
        ASSERT_GT(std::abs(expected - initial), 128 * roundoff);
        EXPECT_NEAR(actual, expected, 8 * roundoff) << "S=" << saturation;
    }
}

// Motivation: At small positive R^2, float residual/step squares can vanish and
// suppress resolvable corrections. A manufactured constant-rate equation
// isolates both stopping paths against an analytic root, including a zero
// initial guess, independently of the small-radius particle-physics Jacobian.
TEST(CondensationScalar, FloatNewtonResolvesSmallPositiveRadiusSquared)
{
    using RT = float;
    constexpr RT cases[][2] = {
        {RT(1.e-18), RT(1.e-23)},
        {RT(1.e-18), RT(-1.e-23)},
        {RT(0), RT(1.e-23)}
    };
    for (const RT mu : {RT(1), RT(1.e12)}) {
        for (const auto& point : cases) {
            const RT initial = point[0];
            const ConstantRadiusSquaredRate ode{mu * point[1]};
            SCOPED_TRACE(::testing::Message() << "mu=" << mu << " initial=" << initial
                                             << " rate=" << ode.rate);
            const NewtonSolver<ConstantRadiusSquaredRate, RT> newton{
                ode, RT(1.e-6), RT(0), RT(1.e-12), 10};
            RT actual = initial;
            RT rhs = mu * initial;
            RT absolute_residual = 0, relative_residual = 0;
            bool success = false;
            newton(actual, rhs, mu, RT(0), RT(0), RT(0), RT(0), RT(0),
                   RT(0), RT(0), RT(1), absolute_residual, relative_residual, success);
            ASSERT_TRUE(success);

            const double expected = (double(rhs) + double(ode.rate)) / double(mu);
            const double roundoff = std::numeric_limits<RT>::epsilon()
                                  * std::max(double(initial), expected);
            ASSERT_GT(std::abs(expected - initial), 64 * roundoff);
            EXPECT_NEAR(actual, expected, 2 * roundoff);
            EXPECT_LE(std::abs(double(mu) * actual - double(rhs) - double(ode.rate)),
                      2 * double(mu) * roundoff);
        }
    }
}
