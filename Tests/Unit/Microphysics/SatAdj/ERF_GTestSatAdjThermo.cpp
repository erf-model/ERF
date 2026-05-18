#include <array>

#include <gtest/gtest.h>

#include "ERF_GTestSatAdjCommon.H"

using namespace satadj_test;

TEST(SatAdjThermo, SaturationDerivativeConsistency)
{
    const std::array<amrex::Real, 5> temperatures = {
        amrex::Real(240.0), amrex::Real(260.0), amrex::Real(273.15), amrex::Real(290.0), amrex::Real(310.0)};
    const std::array<amrex::Real, 3> pressures = {
        amrex::Real(500.0), amrex::Real(800.0), amrex::Real(1000.0)};

    for (const amrex::Real pres_mbar : pressures) {
        amrex::Real qsat_prev = -amrex::Real(1.0);

        for (const amrex::Real tabs : temperatures) {
            amrex::Real qsat_local;
            amrex::Real dqsat_local;
            erf_qsatw(tabs, pres_mbar, qsat_local);
            erf_dtqsatw(tabs, pres_mbar, dqsat_local);

            const amrex::Real esat = erf_esatw(tabs);
            const amrex::Real expected = Rd_on_Rv * erf_dtesatw(tabs) * pres_mbar /
                                         ((pres_mbar - esat) * (pres_mbar - esat));

            EXPECT_NEAR(dqsat_local, expected, scaled_tol(expected, kThermoTolFactor));
            EXPECT_GT(dqsat_local, amrex::Real(0.0));

            if (qsat_prev > amrex::Real(0.0)) {
                EXPECT_GT(qsat_local, qsat_prev);
            }
            qsat_prev = qsat_local;
        }
    }

    for (const amrex::Real tabs : temperatures) {
        amrex::Real qsat_prev = qsat(tabs, pressures.front());

        for (int ip = 1; ip < static_cast<int>(pressures.size()); ++ip) {
            const amrex::Real qsat_local = qsat(tabs, pressures[ip]);
            EXPECT_LT(qsat_local, qsat_prev);
            qsat_prev = qsat_local;
        }
    }
}