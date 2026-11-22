// Unit tests for the divergence damping coefficient of the acoustic substeps
// (ERF_TimestepUtils.H): the coefficient is erf.beta_d over the horizontal acoustic Courant
// number of the substep, so the damping diffusivity beta c^2 dtau does not depend on the
// number of substeps (erf-model/ERF#5003).

#include <ERF_TimestepUtils.H>

#include <gtest/gtest.h>

using amrex::Real;

namespace {
    // c = 340 m/s on a 20 m by 20 m grid
    const Real c_snd = Real(340.0);
    const amrex::GpuArray<Real, AMREX_SPACEDIM> dxinv = {Real(0.05), Real(0.05), Real(0.1)};
}

TEST(DivergenceDampingCoeff, HorizAcousticInverseDt)
{
    // 3-D: c (1/dx + 1/dy) = 340 * 0.1 = 34 1/s
    EXPECT_NEAR(Compute_InvDt_HorizAcoustic(c_snd, Real(1.0), Real(1.0), dxinv, 8, 8), Real(34.0), Real(1.e-4));
    // 2-D runs keep their own direction only: 17 1/s
    EXPECT_NEAR(Compute_InvDt_HorizAcoustic(c_snd, Real(1.0), Real(1.0), dxinv, 8, 1), Real(17.0), Real(1.e-4));
    EXPECT_NEAR(Compute_InvDt_HorizAcoustic(c_snd, Real(1.0), Real(1.0), dxinv, 1, 8), Real(17.0), Real(1.e-4));
    // map factors scale the spacing of their own direction: 340 * (2*0.05 + 0.5*0.05) = 42.5 1/s
    EXPECT_NEAR(Compute_InvDt_HorizAcoustic(c_snd, Real(2.0), Real(0.5), dxinv, 8, 8), Real(42.5), Real(1.e-4));
}

TEST(DivergenceDampingCoeff, CoefficientAtCourantNumberOne)
{
    // dtau inv_dt = 1: the coefficient is erf.beta_d
    EXPECT_NEAR(DivergenceDampingCoeff(Real(0.1), Real(1.0)/Real(34.0), Real(34.0)), Real(0.1), Real(1.e-6));
}

TEST(DivergenceDampingCoeff, InverseInCourantNumber)
{
    // Courant number 0.5 -> 0.2; 0.25 -> 0.4; 2 -> 0.05
    EXPECT_NEAR(DivergenceDampingCoeff(Real(0.1), Real(0.5) /Real(34.0), Real(34.0)), Real(0.2),  Real(1.e-6));
    EXPECT_NEAR(DivergenceDampingCoeff(Real(0.1), Real(0.25)/Real(34.0), Real(34.0)), Real(0.4),  Real(1.e-6));
    EXPECT_NEAR(DivergenceDampingCoeff(Real(0.1), Real(2.0) /Real(34.0), Real(34.0)), Real(0.05), Real(1.e-6));
}

TEST(DivergenceDampingCoeff, DiffusivityIndependentOfSubstepCount)
{
    // beta dtau, to which the damping diffusivity is proportional, is beta_d / inv_dt for
    // every number of substeps of a 1 s step
    const Real inv_dt = Compute_InvDt_HorizAcoustic(c_snd, Real(1.0), Real(1.0), dxinv, 8, 8);
    for (int n : {4, 6, 38, 72, 144, 288}) {
        const Real dtau = Real(1.0) / static_cast<Real>(n);
        const Real beta = DivergenceDampingCoeff(Real(0.1), dtau, inv_dt);
        EXPECT_NEAR(beta * dtau, Real(0.1) / Real(34.0), Real(1.e-7)) << "n = " << n;
    }
}

TEST(DivergenceDampingCoeff, NoHorizontalLimit)
{
    // a single column or an empty level stores zero: the coefficient is erf.beta_d itself
    EXPECT_EQ(DivergenceDampingCoeff(Real(0.1), Real(0.01), Real(0.0)), Real(0.1));
    // no damping asked for: none applied
    EXPECT_EQ(DivergenceDampingCoeff(Real(0.0), Real(0.01), Real(34.0)), Real(0.0));
}
