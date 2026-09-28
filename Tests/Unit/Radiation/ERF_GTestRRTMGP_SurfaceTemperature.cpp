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

TEST(RRTMGP_SurfaceTemperature, ConfiguredRdOcpControlsInverseExner)
{
    const amrex::Real pressure = amrex::Real(0.73) * p_0;
    const amrex::Real custom_rdOcp = amrex::Real(0.31);
    const amrex::Real custom_inverse = rrtmgp::inverse_exner(pressure, custom_rdOcp);
    const amrex::Real expected_inverse = one / getExnergivenP(pressure, custom_rdOcp);
    const amrex::Real default_inverse = one / getExnergivenP(pressure, RdoCp);
    const amrex::Real tolerance = amrex::Real(64.0) *
        std::numeric_limits<amrex::Real>::epsilon() * std::abs(expected_inverse);

    EXPECT_NEAR(custom_inverse, expected_inverse, tolerance);
    EXPECT_GT(std::abs(custom_inverse - default_inverse), amrex::Real(100.0) * tolerance);
}

// Motivation: SurfaceLayer stores potential temperature, but RRTMGP's
// surface boundary inputs are absolute temperature. At a pressure below p0,
// using theta directly changes both the surface emission and the bottom
// layer boundary. The fallback must use the Exner-converted temperature
// without modifying the provider-owned input.
TEST(RRTMGP_SurfaceTemperature, SurfaceLayerFallbackConvertsThetaReadOnly)
{
    const amrex::Real expected_temperature =
        kSurfaceTheta * std::pow(kSurfacePressure / p_0, RdoCp);
    const amrex::Real expected_tolerance =
        amrex::Real(64.0) * std::numeric_limits<amrex::Real>::epsilon() *
        expected_temperature;
    ASSERT_NE(expected_temperature, kSurfaceTheta);

    amrex::Real t_sfc = -1.0;
    const amrex::Real lsm_t_sfc = lsm_undefined;
    rrtmgp::resolve_surface_temperature(
        true, true, false, lsm_t_sfc,
        true, kSurfaceTheta, kSurfacePressure, RdoCp, kDefaultTemperature,
        t_sfc);

    EXPECT_NEAR(t_sfc, expected_temperature, expected_tolerance);
    EXPECT_EQ(lsm_t_sfc, lsm_undefined);
}

// Motivation: water columns do not use an available LSM surface temperature
// in this fallback chain. SurfaceLayer theta must take precedence, be
// converted at the physical surface pressure diagnosed from the lowest
// atmospheric cell, without modifying the available LSM field.
TEST(RRTMGP_SurfaceTemperature, SurfaceLayerFallbackTakesPrecedenceOverWaterLsmTemperature)
{
    const amrex::Real expected_temperature =
        kSurfaceTheta * std::pow(kSurfacePressure / p_0, RdoCp);
    const amrex::Real valid_lsm_temperature = amrex::Real(282.0);
    const amrex::Real expected_tolerance =
        amrex::Real(64.0) * std::numeric_limits<amrex::Real>::epsilon() *
        expected_temperature;

    amrex::Real t_sfc = -1.0;
    const amrex::Real lsm_t_sfc = valid_lsm_temperature;
    rrtmgp::resolve_surface_temperature(
        false, true, true, valid_lsm_temperature,
        true, kSurfaceTheta, kSurfacePressure, RdoCp, kDefaultTemperature,
        t_sfc);

    EXPECT_NEAR(t_sfc, expected_temperature, expected_tolerance);
    EXPECT_NE(t_sfc, valid_lsm_temperature);
    EXPECT_EQ(lsm_t_sfc, valid_lsm_temperature);
}

// Motivation: LSM t_sfc is already an absolute-temperature contract. The
// SurfaceLayer conversion must be limited to the fallback source, so a valid
// LSM value must pass through unchanged even when SurfaceLayer is available.
TEST(RRTMGP_SurfaceTemperature, ValidLsmTemperatureIsNotExnerConverted)
{
    const amrex::Real valid_lsm_temperature = amrex::Real(282.0);
    amrex::Real t_sfc = -1.0;
    const amrex::Real lsm_t_sfc = valid_lsm_temperature;
    rrtmgp::resolve_surface_temperature(
        true, true, true, valid_lsm_temperature,
        true, kSurfaceTheta, kSurfacePressure, RdoCp, kDefaultTemperature,
        t_sfc);

    EXPECT_EQ(t_sfc, valid_lsm_temperature);
    EXPECT_EQ(lsm_t_sfc, valid_lsm_temperature);
}

// Motivation: the default RRTMGP surface temperature is already absolute
// temperature. It must remain unchanged when no SurfaceLayer or valid LSM
// value supplies the boundary condition.
TEST(RRTMGP_SurfaceTemperature, DefaultTemperatureIsAlreadyAbsolute)
{
    amrex::Real t_sfc = -1.0;
    rrtmgp::resolve_surface_temperature(
        true, false, false, 0.0,
        false, 0.0, kSurfacePressure, RdoCp, kDefaultTemperature,
        t_sfc);

    EXPECT_EQ(t_sfc, kDefaultTemperature);
}

// Motivation: a malformed lowest-cell state must not turn the SurfaceLayer
// fallback into a NaN RRTMGP boundary condition. The absolute fallback remains
// the only valid contract when Exner conversion cannot be evaluated.
TEST(RRTMGP_SurfaceTemperature, InvalidSurfacePressureUsesAbsoluteFallback)
{
    amrex::Real t_sfc = -1.0;
    const amrex::Real lsm_t_sfc = lsm_undefined;
    rrtmgp::resolve_surface_temperature(
        true, true, false, lsm_undefined,
        true, kSurfaceTheta, amrex::Real(-1.0), RdoCp, kDefaultTemperature,
        t_sfc);

    EXPECT_EQ(t_sfc, kDefaultTemperature);
    EXPECT_EQ(lsm_t_sfc, lsm_undefined);
}
