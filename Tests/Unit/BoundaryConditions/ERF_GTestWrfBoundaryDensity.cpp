#include <gtest/gtest.h>

namespace {

struct BoundaryTargets {
    double rho;
    double rhotheta;
    double rhoqv;
};

BoundaryTargets make_targets (double rho_erf, double rho_wrf,
                               double theta_wrf, double qv_wrf,
                               bool use_wrf_density)
{
    const double rho = use_wrf_density ? rho_wrf : rho_erf;
    return {rho, rho * theta_wrf, rho * qv_wrf};
}

double davies_density_rhs (double rho_erf, double rho_wrf,
                           double nudge_factor, double dt)
{
    return (rho_wrf - rho_erf) / (nudge_factor * dt);
}

} // namespace

TEST(WrfBoundaryDensity, ConservedTargetsUseWrfDensityWithoutSquaring)
{
    const auto target = make_targets(1.0, 1.2, 300.0, 0.01, true);
    EXPECT_DOUBLE_EQ(target.rho, 1.2);
    EXPECT_DOUBLE_EQ(target.rhotheta, 360.0);
    EXPECT_DOUBLE_EQ(target.rhoqv, 0.012);
}

TEST(WrfBoundaryDensity, DaviesRelaxationUsesDensityFactor)
{
    EXPECT_DOUBLE_EQ(davies_density_rhs(1.0, 1.2, 50.0, 2.0), 0.002);
}

TEST(WrfBoundaryDensity, TimeInterpolationUsesOneConsistentDensity)
{
    const double alpha = 0.5;
    const double rho = (1.0 - alpha) * 1.0 + alpha * 1.4;
    const auto target = make_targets(1.0, rho, 300.0, 0.01, true);
    EXPECT_DOUBLE_EQ(target.rho, 1.2);
    EXPECT_DOUBLE_EQ(target.rhotheta, 360.0);
    EXPECT_DOUBLE_EQ(target.rhoqv, 0.012);
}

TEST(WrfBoundaryDensity, DisabledModeRetainsErfDensity)
{
    const auto target = make_targets(1.0, 1.2, 300.0, 0.01, false);
    EXPECT_DOUBLE_EQ(target.rho, 1.0);
    EXPECT_DOUBLE_EQ(target.rhotheta, 300.0);
    EXPECT_DOUBLE_EQ(target.rhoqv, 0.01);
}
