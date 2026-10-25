#include <cmath>

#include <gtest/gtest.h>

#include <ERF_IBSEBSolar.H>

// ibseb::solar_declination() and ibseb::equation_of_time() are Spencer's
// (1971) Fourier series, with the coefficients of Iqbal (1983), eqs. 1.3.1
// and 1.4.1 (also pvlib's declination_spencer71 and equation_of_time_spencer71):
//   delta = 0.006918 - 0.399912 cos(g) + 0.070257 sin(g)
//         - 0.006758 cos(2g) + 0.000907 sin(2g)
//         - 0.002697 cos(3g) + 0.00148  sin(3g)
//   E     = 229.18 (0.000075 + 0.001868 cos(g) - 0.032077 sin(g)
//         - 0.014615 cos(2g) - 0.040849 sin(2g))          [minutes]
// with Spencer's day angle g = 2 pi (day - 1) / 365. The values below were
// evaluated once from these formulas, outside ERF. Until erf-model/ERF#4657
// the third harmonic of delta was -0.00248 cos(3g) + 0.00031 sin(3g), and
// until erf-model/ERF#4656 the day angle divided by 365.25.

TEST(IBSEBSolarDeclination, MatchesSpencerSeries)
{
    const double tol = (sizeof(amrex::Real) == 8) ? 1.0e-9 : 1.0e-5;
    const struct { double day; double delta; } ref[] = {
        {  1.0, -0.402449},
        { 35.0, -0.2874653347989078},
        { 80.0, -0.0011505915019576882},
        {172.0,  0.40931542032971796},
        {266.0,  0.004342022176745865},
        {355.0, -0.40875419804612767},
    };
    for (const auto& r : ref) {
        EXPECT_NEAR(ibseb::solar_declination(r.day), r.delta, tol) << "day " << r.day;
    }
}

TEST(IBSEBSolarDeclination, EquationOfTimeMatchesSpencerSeries)
{
    // minutes; the series reaches about 16 min, so the float tolerance is
    // relative to that
    const double tol = (sizeof(amrex::Real) == 8) ? 1.0e-8 : 2.0e-4;
    const struct { double day; double eot; } ref[] = {
        {  1.0,  -2.90416896},
        { 35.0, -13.61411487319927},
        { 80.0,  -7.858093122027174},
        {172.0,  -1.3282368002224763},
        {266.0,   7.641989853169503},
        {355.0,   2.17052594886478},
    };
    for (const auto& r : ref) {
        EXPECT_NEAR(ibseb::equation_of_time(r.day), r.eot, tol) << "day " << r.day;
    }
}

// On Spencer's 365-day angle, day 366 (Dec 31 of a leap year) closes the
// cycle: the declination, the equation of time and the Earth-Sun distance
// come back to their day-1 values. With 365.25 in two of the three they did
// not (the declination was off by 3.2e-4 rad, the equation of time by 0.11 min).
TEST(IBSEBSolarDeclination, OneDayAngleForTheWholeSeries)
{
    const double tol = (sizeof(amrex::Real) == 8) ? 1.0e-12 : 1.0e-6;
    EXPECT_NEAR(ibseb::spencer_day_angle(366.0), 2.0 * PI, tol);
    EXPECT_NEAR(ibseb::solar_declination(366.0), ibseb::solar_declination(1.0), tol);
    EXPECT_NEAR(ibseb::equation_of_time(366.0) / 16.0, ibseb::equation_of_time(1.0) / 16.0, tol);
    EXPECT_NEAR(ibseb::earth_sun_distance_factor(366.0), ibseb::earth_sun_distance_factor(1.0), tol);
}

// At the June solstice (day 172) the declination is the obliquity of the
// ecliptic; an ephemeris gives 23.437 deg for 2025-2028. Spencer's series
// gives 23.452 deg there, the old third harmonic 23.403 deg. The band is
// narrower than the series' own accuracy (about 0.03 deg), so this is a
// check on the coefficients, not a bound on the physics.
TEST(IBSEBSolarDeclination, JuneSolsticeIsTheObliquity)
{
    const double deg = ibseb::solar_declination(172.0) * 180.0 / PI;
    EXPECT_NEAR(deg, 23.44, 0.02);
}
