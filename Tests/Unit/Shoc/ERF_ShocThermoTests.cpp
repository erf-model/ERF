#include "ERF_ShocThermoUtils.H"

#include <gtest/gtest.h>

namespace
{
using namespace amrex::literals;

void
expect_round_trip (amrex::Real theta,
                   amrex::Real exner,
                   amrex::Real qc,
                   amrex::Real qi)
{
    const amrex::Real thetal = shoc::thetal_from_theta(theta, qc, qi, exner);
    const amrex::Real tabs = shoc::temperature_from_thetal(thetal, qc, qi, exner);
    const amrex::Real theta_back = tabs / exner;

    EXPECT_NEAR(theta_back, theta, 1.0e-12_rt);
    EXPECT_NEAR(shoc::thetal_from_theta(theta_back, qc, qi, exner), thetal, 1.0e-12_rt);
    EXPECT_GT(tabs, shoc::constants::min_temp());
}
} // namespace

TEST(ShocThermo, CondensateLatentHeatIsPhaseAware)
{
    const struct {
        amrex::Real qc;
        amrex::Real qi;
        amrex::Real expected;
    } cases[] = {
        { 1.0e-3_rt, 0.0_rt, L_v * 1.0e-3_rt / Cp_d },
        { 0.0_rt, 1.0e-3_rt, shoc::latent_sublimation() * 1.0e-3_rt / Cp_d },
        { 0.7e-3_rt, 0.3e-3_rt,
          (L_v * 0.7e-3_rt + shoc::latent_sublimation() * 0.3e-3_rt) / Cp_d },
    };

    for (const auto& tc : cases) {
        EXPECT_NEAR(shoc::condensate_latent_heat(tc.qc, tc.qi) / Cp_d,
                    tc.expected, 1.0e-15_rt);
    }
}

TEST(ShocThermo, ThetaAndThetalRoundTripIsConsistent)
{
    const amrex::Real exner = 0.9_rt;
    expect_round_trip(302.0_rt, exner, 1.0e-3_rt, 0.0_rt);
    expect_round_trip(302.0_rt, exner, 0.0_rt, 1.0e-3_rt);
    expect_round_trip(302.0_rt, exner, 0.7e-3_rt, 0.3e-3_rt);
}

TEST(ShocThermo, ColdNoIceSeedReconstructionDoesNotCreateIce)
{
    const amrex::Real thetal = 260.0_rt;
    const amrex::Real qw = 0.010_rt;
    const amrex::Real exner = 1.0_rt;
    const amrex::Real pdf_ql = 0.001_rt;
    const amrex::Real qi_seed = 0.0_rt;

    amrex::Real tabs = 0.0_rt;
    amrex::Real qv = 0.0_rt;
    amrex::Real qc = 0.0_rt;
    amrex::Real qi = 0.0_rt;

    shoc::reconstruct_pdf_state(thetal, qw, exner, qi_seed, pdf_ql,
                                tabs, qv, qc, qi);

    EXPECT_LT(shoc::temperature_from_thetal(thetal, pdf_ql, 0.0_rt, exner),
              shoc::constants::freezing_temp());
    EXPECT_NEAR(qc, pdf_ql, 1.0e-14_rt);
    EXPECT_NEAR(qi, 0.0_rt, 1.0e-14_rt);
    EXPECT_NEAR(qv, qw - pdf_ql, 1.0e-14_rt);
    EXPECT_NEAR(tabs, thetal + L_v * pdf_ql / Cp_d, 1.0e-12_rt);
    EXPECT_NEAR(tabs, shoc::temperature_from_thetal(thetal, pdf_ql, qi_seed, exner),
                1.0e-12_rt);
    EXPECT_LT(tabs, shoc::constants::freezing_temp());
}

TEST(ShocThermo, ExistingIceSeedContributesSublimationHeatWithoutNewIceSource)
{
    const amrex::Real thetal = 260.0_rt;
    const amrex::Real qw = 0.010_rt;
    const amrex::Real exner = 1.0_rt;
    const amrex::Real pdf_ql = 0.001_rt;
    const amrex::Real qi_seed = 1.0e-6_rt;

    amrex::Real tabs = 0.0_rt;
    amrex::Real qv = 0.0_rt;
    amrex::Real qc = 0.0_rt;
    amrex::Real qi = 0.0_rt;

    shoc::reconstruct_pdf_state(thetal, qw, exner, qi_seed, pdf_ql,
                                tabs, qv, qc, qi);

    EXPECT_LT(shoc::temperature_from_thetal(thetal, pdf_ql, qi_seed, exner),
              shoc::constants::freezing_temp());
    EXPECT_NEAR(qc, pdf_ql, 1.0e-14_rt);
    EXPECT_NEAR(qi, qi_seed, 1.0e-14_rt);
    EXPECT_NEAR(qv, qw - qc - qi, 1.0e-14_rt);
    EXPECT_NEAR(tabs, shoc::temperature_from_thetal(thetal, qc, qi, exner),
                1.0e-12_rt);
    EXPECT_NEAR(tabs,
                thetal + (L_v * pdf_ql + shoc::latent_sublimation() * qi_seed) / Cp_d,
                1.0e-12_rt);
}

TEST(ShocThermo, WarmReconstructionRemainsLiquidWhenIceSeedIsZero)
{
    const amrex::Real thetal = 285.0_rt;
    const amrex::Real qw = 0.010_rt;
    const amrex::Real exner = 1.0_rt;
    const amrex::Real pdf_ql = 0.001_rt;
    const amrex::Real qi_seed = 0.0_rt;

    amrex::Real tabs = 0.0_rt;
    amrex::Real qv = 0.0_rt;
    amrex::Real qc = 0.0_rt;
    amrex::Real qi = 0.0_rt;

    shoc::reconstruct_pdf_state(thetal, qw, exner, qi_seed, pdf_ql,
                                tabs, qv, qc, qi);

    EXPECT_GT(shoc::temperature_from_thetal(thetal, pdf_ql, 0.0_rt, exner),
              shoc::constants::freezing_temp());
    EXPECT_NEAR(qc, pdf_ql, 1.0e-14_rt);
    EXPECT_NEAR(qi, 0.0_rt, 1.0e-14_rt);
    EXPECT_NEAR(qv, qw - pdf_ql, 1.0e-14_rt);
    EXPECT_NEAR(tabs, thetal + L_v * pdf_ql / Cp_d, 1.0e-12_rt);
}

TEST(ShocThermo, MoistEnergyDistinguishesLiquidAndIce)
{
    const amrex::Real tabs = 280.0_rt;
    const amrex::Real qv = 0.008_rt;
    const amrex::Real liquid_energy = shoc::moist_energy(tabs, 0.0_rt, qv,
                                                         1.0e-3_rt, 0.0_rt,
                                                         0.0_rt, 0.0_rt, 0.0_rt);
    const amrex::Real ice_energy = shoc::moist_energy(tabs, 0.0_rt, qv,
                                                      0.0_rt, 1.0e-3_rt,
                                                      0.0_rt, 0.0_rt, 0.0_rt);

    EXPECT_NEAR(liquid_energy - ice_energy,
                shoc::constants::latent_ice() * 1.0e-3_rt,
                1.0e-12_rt);
}
