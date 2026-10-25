#include <cmath>

#include <gtest/gtest.h>

#include <ERF_IBSEBSolar.H>

// ibseb::solar_declination() is Spencer's (1971) Fourier series, with the
// coefficients of Iqbal (1983), eq. 1.3.1:
//   delta = 0.006918 - 0.399912 cos(g) + 0.070257 sin(g)
//         - 0.006758 cos(2g) + 0.000907 sin(2g)
//         - 0.002697 cos(3g) + 0.00148  sin(3g),   g = 2 pi (day - 1) / 365.25.
// Spencer's day angle divides by 365; ERF divides by 365.25, and the values
// below use ERF's day angle. They were evaluated once from the published
// coefficients, outside ERF. Until erf-model/ERF#4657 the third harmonic was
// -0.00248 cos(3g) + 0.00031 sin(3g), which moves delta by 2e-4 to 1.2e-3 rad
// on these days.

TEST(IBSEBSolarDeclination, MatchesSpencerSeries)
{
    const double tol = (sizeof(amrex::Real) == 8) ? 1.0e-9 : 1.0e-5;
    const struct { double day; double delta; } ref[] = {
        {  1.0, -0.402449},
        { 35.0, -0.2875854048491383},
        { 80.0, -0.001523377106028515},
        {172.0,  0.4093003372479783},
        {266.0,  0.00557348280093537},
        {355.0, -0.4087076730501455},
    };
    for (const auto& r : ref) {
        EXPECT_NEAR(ibseb::solar_declination(r.day), r.delta, tol) << "day " << r.day;
    }
}

// At the June solstice (day 172) the declination is the obliquity of the
// ecliptic; an ephemeris gives 23.437 deg for 2025-2028. Spencer's series
// gives 23.451 deg there, the old third harmonic 23.403 deg. The band is
// narrower than the series' own accuracy (about 0.03 deg), so this is a
// check on the coefficients, not a bound on the physics.
TEST(IBSEBSolarDeclination, JuneSolsticeIsTheObliquity)
{
    const double deg = ibseb::solar_declination(172.0) * 180.0 / PI;
    EXPECT_NEAR(deg, 23.44, 0.02);
}
