#include <cmath>

#include <gtest/gtest.h>

#include <ERF_Constants.H>
#include <ERF_OrbCosZenith.H>
#include <ERF_RadStruct.H>
#include <ERF_TwoStreamRadiation.H>
#include <ERF_IBSEBSolar.H>

// erf.ibseb.sun_mode = two_stream: the building faces must see the sun the
// two-stream columns see. The faces build their sun from a declination and an
// hour angle (ibseb::solar_zenith / solar_azimuth); the columns take
// orbital_cos_zenith_instant() of the same calendar day. These tests pin the
// hour angle that makes the two the same sun, the date the shared
// two_stream_sun_date() reads, and the side of the sky the sun is on.

namespace {

double faces_cos_zenith (double calday, double lat_deg, double lon_deg, double declin)
{
    const amrex::Real ha = ibseb::orbital_hour_angle(calday, lon_deg * PI / 180.0);
    return std::cos(ibseb::solar_zenith(lat_deg, declin, ha));
}

} // namespace

// The faces' zenith from orbital_hour_angle() is the columns' zenith, at every
// hour of the day, at several sites (both hemispheres, both sides of Greenwich)
// and declinations. A hour angle off by the equation of time, the error the
// prescribed Spencer sun makes against this formula, moves the cosine by up to
// 0.07 in this range.
TEST(IBSEBTwoStreamSun, FacesZenithIsTheColumnsZenith)
{
    const double tol = (sizeof(amrex::Real) == 8) ? 1.0e-12 : 1.0e-5;
    for (const double lat : {40.0, -33.9, 0.0, 64.8}) {
        for (const double lon : {-100.0, 151.2, 0.0, -179.5}) {
            for (const double declin : {0.2894, -0.4091, 0.0}) {
                for (int hour = 0; hour < 24; ++hour) {
                    const double calday = 218.0 + (hour + 0.25) / 24.0;
                    const double columns = orbital_cos_zenith_instant(calday, lat * PI / 180.0, lon * PI / 180.0, declin);
                    EXPECT_NEAR(faces_cos_zenith(calday, lat, lon, declin), columns, tol)
                        << "lat " << lat << " lon " << lon << " declin " << declin << " hour " << hour;
                }
            }
        }
    }
}

// The hour angle is negative before local solar noon, zero at it (12:00 UTC
// at Greenwich, 18:40 UTC at 100 W), and wrapped to [-pi, pi); the faces then
// put a morning sun in the east (azimuth clockwise from north in (0, 180)) and
// an afternoon sun in the west, as the prescribed sun does.
TEST(IBSEBTwoStreamSun, HourAngleAndAzimuthFollowLocalSolarTime)
{
    const double tol = (sizeof(amrex::Real) == 8) ? 1.0e-12 : 1.0e-5;
    EXPECT_NEAR(ibseb::orbital_hour_angle(100.5, 0.0), 0.0, tol);
    EXPECT_NEAR(ibseb::orbital_hour_angle(100.0 + (12.0 + 100.0 / 15.0) / 24.0, -100.0 * PI / 180.0), 0.0, 10 * tol);
    // 15:00 UTC at 100 W is 08:20 local solar time: 3 h 40 min before noon.
    const amrex::Real h_morning = ibseb::orbital_hour_angle(218.625, -100.0 * PI / 180.0);
    EXPECT_NEAR(h_morning, -(3.0 + 40.0 / 60.0) * 15.0 * PI / 180.0, 10 * tol);
    for (int hour = 0; hour < 48; ++hour) {
        const amrex::Real h = ibseb::orbital_hour_angle(1.0 + hour / 48.0, 2.0);
        EXPECT_GE(h, -PI);
        EXPECT_LT(h, PI);
    }
    const double declin = 0.2894;   // early August
    const amrex::Real z_m = ibseb::solar_zenith(40.0, declin, h_morning);
    const amrex::Real a_m = ibseb::solar_azimuth(40.0, declin, h_morning, z_m) * 180.0 / PI;
    EXPECT_GT(a_m, 0.0);
    EXPECT_LT(a_m, 180.0);
    const amrex::Real h_afternoon = -h_morning;
    const amrex::Real z_a = ibseb::solar_zenith(40.0, declin, h_afternoon);
    const amrex::Real a_a = ibseb::solar_azimuth(40.0, declin, h_afternoon, z_a) * 180.0 / PI;
    EXPECT_GT(a_a, 180.0);
    EXPECT_LT(a_a, 360.0);
    EXPECT_NEAR(z_m, z_a, tol);
}

// The date the faces and the columns share: 2024-08-05 15:00 UTC (a leap year)
// is calendar day 218.625, the declination is early August's +16.7 degrees,
// and the distance factor is below one (the Earth is near aphelion in July).
// The orbit is cached by year (the erf.rad_orbital_* overrides, which the
// shared routine passes on, are not exercised here).
TEST(IBSEBTwoStreamSun, SharedSunDateOfStartDatetime)
{
    RadChoice rc;
    TwoStreamRadiation::OrbitalCache orbit;
    // 2024-08-05 15:00:00 UTC
    const double epoch = 1722870000.0;
    const TwoStreamSunDate sun = two_stream_sun_date(rc, orbit, epoch);
    EXPECT_NEAR(sun.calday, 218.625, 1.0e-9);
    EXPECT_NEAR(sun.declin * 180.0 / PI, 16.7, 0.3);
    EXPECT_GT(sun.eccf, 0.96);
    EXPECT_LT(sun.eccf, 1.0);
    EXPECT_EQ(orbit.year, 2024);
    // The same call an hour later moves the calendar day by 1/24 and keeps the cached orbit.
    const TwoStreamSunDate later = two_stream_sun_date(rc, orbit, epoch + 3600.0);
    EXPECT_NEAR(later.calday - sun.calday, 1.0 / 24.0, 1.0e-9);
    EXPECT_EQ(orbit.year, 2024);
}
