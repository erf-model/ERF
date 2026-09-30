#include "ERF_StationSampler.H"

#include <cmath>

#include <gtest/gtest.h>

// The local inverse of the latitude/longitude map is the one piece of new
// geometry the station sampler carries.  It is exercised end to end only by a
// run initialized from a WRF or metgrid file, which the regression suite has no
// deck for, so the solve itself is checked here.

namespace {

// The solve is exact up to rounding: 1e-12 in double, a few float epsilons in single
constexpr amrex::Real tol = (sizeof(amrex::Real) == 8) ? amrex::Real(1.0e-12) : amrex::Real(1.0e-6);

// A map whose gradients are (dlat_di, dlon_di) and (dlat_dj, dlon_dj): given an
// offset in cells, the lat/lon difference it produces.
void
forward (amrex::Real di, amrex::Real dj,
         amrex::Real dlat_di, amrex::Real dlon_di,
         amrex::Real dlat_dj, amrex::Real dlon_dj,
         amrex::Real& dlat, amrex::Real& dlon)
{
    dlat = dlat_di*di + dlat_dj*dj;
    dlon = dlon_di*di + dlon_dj*dj;
}

} // namespace

// An axis-aligned map: latitude runs with j, longitude with i, which is what a
// mass-point lat/lon array looks like away from the poles on a regular grid.
TEST(StationSamplerLatLonInverse, AxisAligned)
{
    const amrex::Real dlat_di = amrex::Real(0.0),  dlon_di = amrex::Real(0.01);
    const amrex::Real dlat_dj = amrex::Real(0.009), dlon_dj = amrex::Real(0.0);

    amrex::Real di = 0, dj = 0;
    ASSERT_TRUE(invert_latlon_offset(amrex::Real(0.0045), amrex::Real(-0.003),
                                     dlat_di, dlon_di, dlat_dj, dlon_dj, di, dj));
    EXPECT_NEAR(di, amrex::Real(-0.3), tol);
    EXPECT_NEAR(dj, amrex::Real( 0.5), tol);
}

// A rotated, sheared map, as a Lambert conformal grid gives away from its
// standard parallel: the solve must still return the offset that generated the
// lat/lon difference.
TEST(StationSamplerLatLonInverse, RoundTripsARotatedMap)
{
    const amrex::Real dlat_di = amrex::Real(0.0021), dlon_di = amrex::Real(0.0097);
    const amrex::Real dlat_dj = amrex::Real(0.0088), dlon_dj = amrex::Real(-0.0019);

    for (const amrex::Real di_in : {amrex::Real(-0.5), amrex::Real(0.0), amrex::Real(0.37)}) {
        for (const amrex::Real dj_in : {amrex::Real(-0.25), amrex::Real(0.5)}) {
            amrex::Real dlat = 0, dlon = 0;
            forward(di_in, dj_in, dlat_di, dlon_di, dlat_dj, dlon_dj, dlat, dlon);

            amrex::Real di = 0, dj = 0;
            ASSERT_TRUE(invert_latlon_offset(dlat, dlon, dlat_di, dlon_di, dlat_dj, dlon_dj, di, dj));
            EXPECT_NEAR(di, di_in, tol);
            EXPECT_NEAR(dj, dj_in, tol);
        }
    }
}

// Parallel gradients mean the two index directions move along the same line on
// the sphere, so a lat/lon pair does not determine a cell offset.  The solve has
// to say so rather than divide by a near-zero determinant.
TEST(StationSamplerLatLonInverse, RejectsADegenerateMap)
{
    amrex::Real di = 0, dj = 0;
    EXPECT_FALSE(invert_latlon_offset(amrex::Real(0.001), amrex::Real(0.002),
                                      amrex::Real(0.01), amrex::Real(0.02),
                                      amrex::Real(0.02), amrex::Real(0.04),
                                      di, dj));
    // A map with no gradient at all is degenerate in the same way
    EXPECT_FALSE(invert_latlon_offset(amrex::Real(0.001), amrex::Real(0.002),
                                      amrex::Real(0.0), amrex::Real(0.0),
                                      amrex::Real(0.0), amrex::Real(0.0),
                                      di, dj));
}
