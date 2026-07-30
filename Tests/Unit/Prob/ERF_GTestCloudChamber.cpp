#include <AMReX_Array.H>

#include <cmath>

#include <gtest/gtest.h>

#include "../../../Source/Prob/ERF_CloudChamber.H"

using amrex::GpuArray;
using amrex::Real;

namespace {

// Motivation: the chamber must use a decomposition-independent analytic
// initializer; these tests catch coordinate normalization and endpoint
// regressions without launching a full ERF simulation.
TEST(CloudChamberProfile, LinearProfileAndZeroAmplitude)
{
    erf_cloud_chamber::Config config;
    config.prob_lo = {Real(0.0), Real(0.0), Real(0.0)};
    config.prob_hi = {Real(2.0), Real(2.0), Real(1.0)};
    config.theta_bottom = Real(299.0);
    config.theta_top = Real(280.0);
    config.theta_perturbation_amplitude = Real(0.0);

    EXPECT_DOUBLE_EQ(erf_cloud_chamber::theta_at(config, Real(0.25),
                                                  Real(0.75), Real(0.0)),
                     Real(299.0));
    EXPECT_DOUBLE_EQ(erf_cloud_chamber::theta_at(config, Real(0.25),
                                                  Real(0.75), Real(1.0)),
                     Real(280.0));
    EXPECT_DOUBLE_EQ(erf_cloud_chamber::theta_at(config, Real(0.25),
                                                  Real(0.75), Real(0.5)),
                     Real(289.5));
}

TEST(CloudChamberProfile, PerturbationIsBoundedAndZeroOnVerticalFaces)
{
    erf_cloud_chamber::Config config;
    config.prob_lo = {Real(-1.0), Real(2.0), Real(4.0)};
    config.prob_hi = {Real(1.0), Real(6.0), Real(5.0)};
    config.theta_bottom = Real(300.0);
    config.theta_top = Real(301.0);
    config.theta_perturbation_amplitude = Real(0.2);
    const GpuArray<Real, AMREX_SPACEDIM> length = {Real(2.0), Real(4.0), Real(1.0)};

    const Real low = erf_cloud_chamber::deterministic_perturbation(
        Real(0.0), Real(4.0), Real(4.0), config.prob_lo, length, Real(0.2));
    const Real high = erf_cloud_chamber::deterministic_perturbation(
        Real(0.0), Real(4.0), Real(5.0), config.prob_lo, length, Real(0.2));
    const Real interior = erf_cloud_chamber::deterministic_perturbation(
        Real(0.0), Real(4.0), Real(4.5), config.prob_lo, length, Real(0.2));

    EXPECT_NEAR(low, Real(0.0), Real(1.0e-15));
    EXPECT_NEAR(high, Real(0.0), Real(1.0e-15));
    EXPECT_LE(std::abs(interior), Real(0.2));
}

} // namespace
