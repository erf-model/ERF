#include "ERF_IBSEBWallFunction.H"

#include <gtest/gtest.h>

#include <cmath>

using amrex::Real;

namespace {

const Real tol = (sizeof(Real) == 8) ? Real(1.0e-12) : Real(1.0e-5);

// a^2 for a 5 m reference height over a 0.01 m roughness, the wall function's
// half cell on a 10 m mesh.
Real a2_5m ()
{
    const Real lnm = std::log(Real(5.0) / Real(0.01));
    return Real(0.41) * Real(0.41) / (lnm * lnm);
}

}

TEST(IBSEBLouis, NeutralFactorsAreOne)
{
    Real Fm = 0, Fh = 0;
    ibseb::louis_factors(Real(0.0), Real(500.0), a2_5m(), Fm, Fh);
    EXPECT_NEAR(Fm, Real(1.0), tol);
    EXPECT_NEAR(Fh, Real(1.0), tol);
}

TEST(IBSEBLouis, StableMatchesClosedForm)
{
    for (Real Ri : {Real(0.01), Real(0.1), Real(0.5), Real(2.0)}) {
        Real Fm = 0, Fh = 0;
        ibseb::louis_factors(Ri, Real(500.0), a2_5m(), Fm, Fh);
        const Real q = Real(1.0) + Real(4.7) * Ri;
        EXPECT_NEAR(Fm, Real(1.0) / (q * q), tol) << "Ri = " << Ri;
        EXPECT_NEAR(Fh, Fm, tol) << "Ri = " << Ri;
    }
}

TEST(IBSEBLouis, UnstableMatchesClosedForm)
{
    const Real a2 = a2_5m();
    for (Real Ri : {Real(-0.01), Real(-0.3), Real(-1.0), Real(-5.0)}) {
        Real Fm = 0, Fh = 0;
        ibseb::louis_factors(Ri, Real(500.0), a2, Fm, Fh);
        const Real s = std::sqrt(-Ri * Real(500.0));
        EXPECT_NEAR(Fm, Real(1.0) - Real(9.4) * Ri / (Real(1.0) + Real(7.4) * a2 * Real(9.4) * s), tol) << "Ri = " << Ri;
        EXPECT_NEAR(Fh, Real(1.0) - Real(9.4) * Ri / (Real(1.0) + Real(5.3) * a2 * Real(9.4) * s), tol) << "Ri = " << Ri;
    }
}

TEST(IBSEBLouis, MonotoneAndContinuousThroughNeutral)
{
    const Real a2 = a2_5m();
    Real Fm_prev = 1.0e30, Fh_prev = 1.0e30;
    for (int n = -400; n <= 400; ++n) {
        const Real Ri = Real(n) * Real(0.005);
        Real Fm = 0, Fh = 0;
        ibseb::louis_factors(Ri, Real(500.0), a2, Fm, Fh);
        EXPECT_TRUE(std::isfinite(Fm) && std::isfinite(Fh));
        EXPECT_GT(Fm, Real(0.0));
        EXPECT_LE(Fm, Fm_prev + tol);   // mixing weakens as Ri grows
        EXPECT_LE(Fh, Fh_prev + tol);
        // Heat is mixed at least as strongly as momentum when unstable (C*_h < C*_m).
        if (Ri < 0) { EXPECT_GE(Fh, Fm - tol); }
        Fm_prev = Fm; Fh_prev = Fh;
    }
    Real Fm_m = 0, Fh_m = 0, Fm_p = 0, Fh_p = 0;
    ibseb::louis_factors(Real(-1.0e-6), Real(500.0), a2, Fm_m, Fh_m);
    ibseb::louis_factors(Real( 1.0e-6), Real(500.0), a2, Fm_p, Fh_p);
    EXPECT_NEAR(Fm_m, Fm_p, Real(1.0e-4));
    EXPECT_NEAR(Fh_m, Fh_p, Real(1.0e-4));
}
