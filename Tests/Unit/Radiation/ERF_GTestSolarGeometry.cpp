#include <gtest/gtest.h>
#include <cmath>

#include <ERF_Constants.H>
#include <ERF_SolarGeometry.H>

// Motivation: compute_solar_azimuth_angle had sin and cos of the zenith
// angle swapped, which put cos(azimuth) outside [-1, 1] at solar noon on
// the solstice at 40N; the function was unused, so nothing caught it. These
// pin the standard formula: the sun is due south at solar noon in the
// northern mid-latitudes, east of south in the morning, west of south in
// the afternoon, and the azimuth is undefined (returned as 0) at the pole.

namespace {

constexpr amrex::Real kDeg = PI / 180.0;

amrex::Real zenith_from (amrex::Real lat_rad, amrex::Real decl_rad, amrex::Real hour_rad)
{
    const amrex::Real cosz = std::sin(lat_rad) * std::sin(decl_rad) +
                             std::cos(lat_rad) * std::cos(decl_rad) * std::cos(hour_rad);
    return std::acos(cosz);
}

} // namespace

TEST(SolarGeometry, AzimuthIsDueSouthAtSolarNoonInNorthernMidLatitudes)
{
    const amrex::Real lat = 40.0;
    const amrex::Real decl = 23.44 * kDeg;   // June solstice
    const amrex::Real hour = 0.0;
    const amrex::Real zen = zenith_from(lat * kDeg, decl, hour);
    const amrex::Real az = compute_solar_azimuth_angle(lat, decl, hour, zen);
    EXPECT_NEAR(az, PI, 1.0e-6);
}

TEST(SolarGeometry, AzimuthMovesFromEastOfSouthToWestOfSouth)
{
    const amrex::Real lat = 40.0;
    const amrex::Real decl = 0.0;            // equinox
    const amrex::Real morning = -45.0 * kDeg; // 09:00 solar time
    const amrex::Real evening = +45.0 * kDeg; // 15:00 solar time
    const amrex::Real az_am = compute_solar_azimuth_angle(lat, decl, morning, zenith_from(lat * kDeg, decl, morning));
    const amrex::Real az_pm = compute_solar_azimuth_angle(lat, decl, evening, zenith_from(lat * kDeg, decl, evening));
    EXPECT_GT(az_am, 0.5 * PI);
    EXPECT_LT(az_am, PI);
    EXPECT_GT(az_pm, PI);
    EXPECT_LT(az_pm, 1.5 * PI);
    // Symmetric about south
    EXPECT_NEAR(PI - az_am, az_pm - PI, 1.0e-6);
}

TEST(SolarGeometry, AzimuthMatchesTheClosedFormAtEquinoxSunrise)
{
    // At the equinox the sun rises due east (azimuth pi/2) everywhere
    // except the poles: hour angle -90 deg, zenith 90 deg.
    const amrex::Real lat = 40.0;
    const amrex::Real decl = 0.0;
    const amrex::Real hour = -90.0 * kDeg;
    const amrex::Real zen = zenith_from(lat * kDeg, decl, hour);
    EXPECT_NEAR(zen, 0.5 * PI, 1.0e-9);
    const amrex::Real az = compute_solar_azimuth_angle(lat, decl, hour, zen);
    EXPECT_NEAR(az, 0.5 * PI, 1.0e-6);
}

TEST(SolarGeometry, AzimuthIsDefinedAsZeroAtThePole)
{
    const amrex::Real decl = 23.44 * kDeg;
    const amrex::Real hour = 30.0 * kDeg;
    const amrex::Real zen = zenith_from(90.0 * kDeg, decl, hour);
    EXPECT_DOUBLE_EQ(compute_solar_azimuth_angle(90.0, decl, hour, zen), 0.0);
}
