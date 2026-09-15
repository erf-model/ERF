#include <cmath>
#include <limits>

#include <gtest/gtest.h>

#include <ERF_Constants.H>
#include <ERF_RRTMGP_SurfaceTemperature.H>

namespace {

constexpr amrex::Real kSurfacePressure = amrex::Real(0.9) * p_0;
constexpr amrex::Real kSurfaceTheta = amrex::Real(300.0);
constexpr amrex::Real kDefaultTemperature = amrex::Real(280.0);

} // namespace

// Motivation: SurfaceLayer stores potential temperature, but RRTMGP's
// surface boundary inputs are absolute temperature. At a pressure below p0,
// using theta directly changes both the surface emission and the bottom
// layer boundary. The fallback and its LSM writeback must use the same
// Exner-converted temperature.
TEST(RRTMGP_SurfaceTemperature, SurfaceLayerFallbackConvertsThetaForAllBoundaryOutputs)
{
    const amrex::Real expected_temperature =
        kSurfaceTheta * std::pow(kSurfacePressure / p_0, RdoCp);
    const amrex::Real expected_tolerance =
        amrex::Real(64.0) * std::numeric_limits<amrex::Real>::epsilon() *
        expected_temperature;
    ASSERT_NE(expected_temperature, kSurfaceTheta);

    amrex::Real t_sfc = -1.0;
    amrex::Real t_lev_bottom = -2.0;
    amrex::Real lsm_t_sfc = lsm_undefined;
    rrtmgp::resolve_surface_temperature(
        true, true, false, lsm_t_sfc,
        true, kSurfaceTheta, kSurfacePressure, kDefaultTemperature,
        t_sfc, t_lev_bottom, &lsm_t_sfc);

    EXPECT_NEAR(t_sfc, expected_temperature, expected_tolerance);
    EXPECT_EQ(t_lev_bottom, t_sfc);
    EXPECT_NEAR(lsm_t_sfc, expected_temperature, expected_tolerance);
}

// Motivation: LSM t_sfc is already an absolute-temperature contract. The
// SurfaceLayer conversion must be limited to the fallback source, so a valid
// LSM value must pass through unchanged even when SurfaceLayer is available.
TEST(RRTMGP_SurfaceTemperature, ValidLsmTemperatureIsNotExnerConverted)
{
    const amrex::Real valid_lsm_temperature = amrex::Real(282.0);
    amrex::Real t_sfc = -1.0;
    amrex::Real t_lev_bottom = -2.0;
    amrex::Real lsm_t_sfc = valid_lsm_temperature;
    rrtmgp::resolve_surface_temperature(
        true, true, true, valid_lsm_temperature,
        true, kSurfaceTheta, kSurfacePressure, kDefaultTemperature,
        t_sfc, t_lev_bottom, &lsm_t_sfc);

    EXPECT_EQ(t_sfc, valid_lsm_temperature);
    EXPECT_EQ(t_lev_bottom, valid_lsm_temperature);
    EXPECT_EQ(lsm_t_sfc, valid_lsm_temperature);
}

// Motivation: the default RRTMGP surface temperature is already absolute
// temperature. It must remain unchanged when no SurfaceLayer or valid LSM
// value supplies the boundary condition.
TEST(RRTMGP_SurfaceTemperature, DefaultTemperatureIsAlreadyAbsolute)
{
    amrex::Real t_sfc = -1.0;
    amrex::Real t_lev_bottom = -2.0;
    rrtmgp::resolve_surface_temperature(
        true, false, false, 0.0,
        false, 0.0, kSurfacePressure, kDefaultTemperature,
        t_sfc, t_lev_bottom, nullptr);

    EXPECT_EQ(t_sfc, kDefaultTemperature);
    EXPECT_EQ(t_lev_bottom, kDefaultTemperature);
}
