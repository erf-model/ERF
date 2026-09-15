#include <cmath>

#include <gtest/gtest.h>

#include <ERF_Constants.H>
#include <ERF_OrbCosZenith.H>

// The two-stream model places the sun with the orbital code the RRTMGP
// interface uses (ERF_OrbCosZenith.H). These tests pin the two helpers the
// two-stream path added to that header and the quantities it takes from the
// existing routines:
//   1. orbital_calday: 1.0 at 00:00 UTC on 1 January, leap-aware after
//      February, fraction of the day from the seconds.
//   2. orbital_cos_zenith_instant (device-callable) is the formula
//      orbital_cos_zenith uses without an averaging interval, and gives the
//      expected sun at noon and midnight on the equator.
//   3. orbital_params + orbital_decl give the June-solstice declination and
//      the Earth-Sun distance factor at perihelion and aphelion, which set the
//      irradiance when erf.fixed_total_solar_irradiance is not given.

namespace {
constexpr double kDeg = PI / 180.0;
}

TEST(OrbitalGeometry, CalendarDayCountsFromNewYearAndKnowsLeapYears)
{
    EXPECT_DOUBLE_EQ(orbital_calday(2021, 1, 1, 0), 1.0);
    EXPECT_DOUBLE_EQ(orbital_calday(2021, 1, 1, 43200), 1.5);
    EXPECT_DOUBLE_EQ(orbital_calday(2021, 3, 1, 0), 60.0);
    EXPECT_DOUBLE_EQ(orbital_calday(2020, 3, 1, 0), 61.0);      // leap year
    EXPECT_DOUBLE_EQ(orbital_calday(2020, 2, 28, 0), 59.0);     // before the leap day
    EXPECT_DOUBLE_EQ(orbital_calday(2100, 3, 1, 0), 60.0);      // century, not leap
    EXPECT_DOUBLE_EQ(orbital_calday(2000, 3, 1, 0), 61.0);      // 400-year rule
    EXPECT_DOUBLE_EQ(orbital_calday(2021, 6, 21, 43200), 172.5);
    EXPECT_DOUBLE_EQ(orbital_calday(2021, 12, 31, 86399), 365.0 + 86399.0 / 86400.0);
}

TEST(OrbitalGeometry, InstantaneousCosZenithMatchesTheHostRoutine)
{
    const double lats[] = {-60.0 * kDeg, 0.0, 40.0 * kDeg, 75.0 * kDeg};
    const double lons[] = {-120.0 * kDeg, 0.0, 90.0 * kDeg};
    const double declins[] = {-23.44 * kDeg, 0.0, 23.44 * kDeg};
    const double jdays[] = {1.0, 80.25, 172.5, 300.9};
    for (double lat : lats) {
        for (double lon : lons) {
            for (double declin : declins) {
                for (double jday : jdays) {
                    // orbital_cos_zenith takes lvalue references.
                    double jd = jday, la = lat, lo = lon, de = declin;
                    const double host = orbital_cos_zenith(jd, la, lo, de);
                    const double inst = orbital_cos_zenith_instant(jday, lat, lon, declin);
                    EXPECT_DOUBLE_EQ(inst, host) << "lat " << lat << " lon " << lon
                                                 << " declin " << declin << " jday " << jday;
                }
            }
        }
    }
}

TEST(OrbitalGeometry, EquatorialSunAtNoonAndMidnight)
{
    // Zero declination, longitude 0: overhead at 12:00 UTC, antipodal at 00:00.
    EXPECT_NEAR(orbital_cos_zenith_instant(80.5, 0.0, 0.0, 0.0), 1.0, 1.0e-12);
    EXPECT_NEAR(orbital_cos_zenith_instant(80.0, 0.0, 0.0, 0.0), -1.0, 1.0e-12);
    // 90 degrees east sees noon six hours earlier (06:00 UTC).
    EXPECT_NEAR(orbital_cos_zenith_instant(80.25, 0.0, 90.0 * kDeg, 0.0), 1.0, 1.0e-12);
    // At 40N on the June solstice noon: cos(40 - 23.44 degrees).
    const double declin = 23.44 * kDeg;
    EXPECT_NEAR(orbital_cos_zenith_instant(172.5, 40.0 * kDeg, 0.0, declin),
                std::cos(40.0 * kDeg - declin), 1.0e-12);
}

TEST(OrbitalGeometry, DeclinationAndEarthSunDistanceFactorOfTheDate)
{
    // The Berger (1978) parameters of 2021 from the year alone, as the
    // driver forms them when erf.rad_orbital_* are left unset.
    int year = 2021;
    double eccen = -9999.0, obliq = -9999.0, mvelp = -9999.0;
    double obliqr = 0.0, lambm0 = 0.0, mvelpp = 0.0;
    orbital_params(year, eccen, obliq, mvelp, obliqr, lambm0, mvelpp);
    EXPECT_NEAR(eccen, 0.0167, 5.0e-4);
    EXPECT_NEAR(obliq, 23.44, 0.05);

    double delta = 0.0, eccf = 0.0;
    double calday = orbital_calday(2021, 6, 21, 43200);
    orbital_decl(calday, eccen, mvelpp, lambm0, obliqr, delta, eccf);
    EXPECT_NEAR(delta / kDeg, 23.44, 0.1);      // June solstice
    EXPECT_NEAR(eccf, 0.967, 2.0e-3);           // near aphelion

    calday = orbital_calday(2021, 1, 3, 0);
    orbital_decl(calday, eccen, mvelpp, lambm0, obliqr, delta, eccf);
    EXPECT_NEAR(eccf, 1.034, 2.0e-3);           // near perihelion
    EXPECT_NEAR(delta / kDeg, -22.9, 0.3);

    calday = orbital_calday(2021, 3, 20, 43200);
    orbital_decl(calday, eccen, mvelpp, lambm0, obliqr, delta, eccf);
    EXPECT_NEAR(delta / kDeg, 0.0, 0.5);        // equinox (the Berger series is good to a few tenths of a degree)
}
