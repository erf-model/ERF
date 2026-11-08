#include <cmath>

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <gtest/gtest.h>

#include <ERF_IBFaceSet.H>
#include <ERF_IBSEBParams.H>
#include <ERF_IBSEBSolar.H>

// ibseb::solar_declination() and ibseb::equation_of_time() are Spencer's
// (1971) Fourier series, with the coefficients of Iqbal (1983), eqs. 1.3.1
// and 1.4.1 (pvlib's declination_spencer71 is the same; its
// equation_of_time_spencer71 uses 0.0000075 for the constant term, 0.93 s lower):
//   delta = 0.006918 - 0.399912 cos(g) + 0.070257 sin(g)
//         - 0.006758 cos(2g) + 0.000907 sin(2g)
//         - 0.002697 cos(3g) + 0.00148  sin(3g)
//   E     = 229.18 (0.000075 + 0.001868 cos(g) - 0.032077 sin(g)
//         - 0.014615 cos(2g) - 0.040849 sin(2g))          [minutes]
// with Spencer's day angle g = 2 pi (day - 1) / 365. The values below were
// evaluated once from these formulas, outside ERF. Until erf-model/ERF#4657
// the third harmonic of delta was -0.00248 cos(3g) + 0.00031 sin(3g), and
// until erf-model/ERF#4656 the day angle divided by 365.25.

namespace {

// Mark the cells of building as building (blank = 1). A free function: nvcc
// rejects an extended device lambda inside a gtest body, which is a private
// member function of the generated test class.
void mark_building (amrex::MultiFab& blank, const amrex::Box& building)
{
    for (amrex::MFIter mfi(blank); mfi.isValid(); ++mfi) {
        auto const& a = blank.array(mfi);
        const amrex::Box bld = building & mfi.validbox();
        amrex::ParallelFor(bld, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { a(i, j, k) = 1.0; });
    }
}

} // namespace

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
// gives 23.452 deg there, the old third harmonic 23.403 deg. Spencer quotes
// 0.0006 rad (0.035 deg) for his fit; for a given year and hour the error
// against an ephemeris reaches about 0.5 deg, since the day carries neither.
// The band is a check on the coefficients, not a bound on the physics.
TEST(IBSEBSolarDeclination, JuneSolsticeIsTheObliquity)
{
    const double deg = ibseb::solar_declination(172.0) * 180.0 / PI;
    EXPECT_NEAR(deg, 23.44, 0.02);
}

// The date advances at 00:00 UTC: a run that starts on day 172 at 18 UTC is
// on day 173 after 6 h and on day 174 after 30 h. Before this, the start
// day's sun repeated however long the run (a 3-day run near an equinox
// drifted about 0.8 deg in declination).
TEST(IBSEBSolarDeclination, DateAdvancesAtMidnightUTC)
{
    EXPECT_EQ(ibseb::day_of_year_at(amrex::Real(172), amrex::Real(18 * 3600)), amrex::Real(172));
    EXPECT_EQ(ibseb::day_of_year_at(amrex::Real(172), amrex::Real(86399)), amrex::Real(172));
    EXPECT_EQ(ibseb::day_of_year_at(amrex::Real(172), amrex::Real(86400)), amrex::Real(173));
    EXPECT_EQ(ibseb::day_of_year_at(amrex::Real(172), amrex::Real(18 * 3600 + 30 * 3600)), amrex::Real(174));
    // the year wraps through the 365-day series: day 365 + 1 day = day 1's sun
    const double tol = (sizeof(amrex::Real) == 8) ? 1.0e-12 : 1.0e-6;
    const amrex::Real d = ibseb::day_of_year_at(amrex::Real(365), amrex::Real(86400 + 3600));
    EXPECT_NEAR(ibseb::solar_declination(d), ibseb::solar_declination(amrex::Real(1)), tol);
}

// The face set's solar sun takes the advancing date: a run that starts on day
// 80 and is 24 h in sees the sun of a run that starts on day 81 at the same
// UTC time. Near the equinox one day moves the declination by about 0.4 deg,
// so a face set that kept the start day would differ by that much.
TEST(IBSEBSolarDeclination, FaceSetSunFollowsTheDate)
{
    using namespace amrex;
    const RealBox rb({0.0, 0.0, 0.0}, {80.0, 80.0, 40.0});
    const Box dom(IntVect(0, 0, 0), IntVect(3, 3, 3));
    const Geometry geom(dom, rb, 0, {1, 1, 0});
    const BoxArray ba(dom);
    const DistributionMapping dm(ba);
    MultiFab blank(ba, dm, 1, 1);
    blank.setVal(0.0);
    mark_building(blank, Box(IntVect(1, 1, 0), IntVect(2, 2, 1)));
    blank.FillBoundary(geom.periodicity());

    auto sun_at = [&] (Real day0, Real time) {
        IBSEBParams params;
        params.enable = true;
        params.sun_mode = "solar";
        params.latitude_deg = 40.0;
        params.longitude_deg = -105.0;
        params.day_of_year = day0;
        params.time_zero_utc_s = 19.0 * 3600.0;
        params.view_n_az = 4;
        params.view_n_el = 2;
        IBFaceSet faces(params, 0);
        faces.build(blank, geom);
        faces.release_labels();
        faces.compute_view_fractions();
        faces.compute_shortwave(time);
        return faces.sun();
    };
    const auto a = sun_at(Real(80), Real(86400));
    const auto b = sun_at(Real(81), Real(0));
    const auto c = sun_at(Real(80), Real(0));
    const Real tol = (sizeof(Real) == 8) ? Real(1.e-12) : Real(1.e-5);
    EXPECT_NEAR(a.zenith, b.zenith, tol);
    EXPECT_NEAR(a.azimuth, b.azimuth, tol);
    EXPECT_NEAR(a.dni, b.dni, Real(1e3) * tol);
    // and the date matters: the start day's sun is about 0.4 deg away
    EXPECT_GT(std::abs(a.zenith - c.zenith) * Real(180.0 / PI), Real(0.2));
}
