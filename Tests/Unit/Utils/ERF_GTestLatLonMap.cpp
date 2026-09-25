#include <cmath>

#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Math.H>

#include <gtest/gtest.h>

#include "ERF_LatLonMap.H"

// Placing a latitude/longitude on the grid (ERF_LatLonMap), as the station
// time-series output does: the nearest mass point of the lat/lon arrays and one
// local linear solve, with longitude differences taken the short way round.

namespace {

using amrex::Real;

// A grid 100 m a cell, rotated by alpha counterclockwise from east, around
// (40 N, 105 W), with lat/lon from the local tangent plane
void
fill_rotated_latlon (amrex::FArrayBox& fab, Real alpha)
{
    const Real R = Real(6371000.0), lat0 = Real(40.0), lon0 = Real(-105.0);
    const Real deg = Real(180.0)/amrex::Math::pi<Real>();
    const auto a = fab.array();
    amrex::LoopOnCpu(fab.box(), [&](int i, int j, int) {
        const Real x = (Real(i) + Real(0.5))*Real(100.0);
        const Real y = (Real(j) + Real(0.5))*Real(100.0);
        const Real east  = x*std::cos(alpha) - y*std::sin(alpha);
        const Real north = x*std::sin(alpha) + y*std::cos(alpha);
        a(i,j,0,0) = lat0 + north/R*deg;
        a(i,j,0,1) = lon0 + east/(R*std::cos(lat0/deg))*deg;
    });
}

} // namespace

TEST(LatLonMap, PlacesAPointOnARotatedGrid)
{
    const amrex::Box dom(amrex::IntVect(0,0,0), amrex::IntVect(39,39,0));
    amrex::FArrayBox fab(dom, 2, amrex::The_Pinned_Arena());
    const Real alpha = amrex::Math::pi<Real>()/Real(6.0);
    fill_rotated_latlon(fab, alpha);

    const amrex::GpuArray<Real,AMREX_SPACEDIM> problo{Real(0.0), Real(0.0), Real(0.0)};
    const amrex::GpuArray<Real,AMREX_SPACEDIM> dx{Real(100.0), Real(100.0), Real(100.0)};

    // The lat/lon of x = 1234, y = 2345 in the rotated plane
    const Real R = Real(6371000.0), deg = Real(180.0)/amrex::Math::pi<Real>();
    const Real x = Real(1234.0), y = Real(2345.0);
    const Real lat = Real(40.0) + (x*std::sin(alpha) + y*std::cos(alpha))/R*deg;
    const Real lon = Real(-105.0) + (x*std::cos(alpha) - y*std::sin(alpha))/(R*std::cos(Real(40.0)/deg))*deg;

    LatLonLocation loc;
    ASSERT_EQ(locate_latlon_on_grid(fab.const_array(), dom, problo, dx, lat, lon, loc), LatLonStatus::Ok);
    // A float holds a latitude of 40 degrees to about 4e-6 degrees, 0.4 m
    const Real xtol = (sizeof(Real) == 8) ? Real(0.5) : Real(5.0);
    EXPECT_NEAR(loc.x, x, xtol);
    EXPECT_NEAR(loc.y, y, xtol);

    // A point well outside the grid
    EXPECT_EQ(locate_latlon_on_grid(fab.const_array(), dom, problo, dx, lat + Real(1.0), lon, loc),
              LatLonStatus::TooFar);
}

TEST(LatLonMap, PlacesAPointOnAGridAcrossTheAntimeridian)
{
    // 40 x 40 cells of 0.005 degrees in longitude centred on 180 E, given as a
    // WRF file gives them: 179.9 ... 180 and then -179.995 ... -179.9
    const amrex::Box dom(amrex::IntVect(0,0,0), amrex::IntVect(39,39,0));
    amrex::FArrayBox fab(dom, 2, amrex::The_Pinned_Arena());
    const auto a = fab.array();
    amrex::LoopOnCpu(dom, [&](int i, int j, int) {
        Real lon = Real(179.9) + (Real(i) + Real(0.5))*Real(0.005);
        if (lon > Real(180.0)) { lon -= Real(360.0); }
        a(i,j,0,0) = Real(-17.0) + (Real(j) + Real(0.5))*Real(0.005);
        a(i,j,0,1) = lon;
    });
    const amrex::GpuArray<Real,AMREX_SPACEDIM> problo{Real(0.0), Real(0.0), Real(0.0)};
    const amrex::GpuArray<Real,AMREX_SPACEDIM> dx{Real(500.0), Real(500.0), Real(100.0)};

    // x = 20.3 cells, y = 12.6 cells: longitude -179.99835 on the far side of 180
    LatLonLocation loc;
    ASSERT_EQ(locate_latlon_on_grid(fab.const_array(), dom, problo, dx,
                                    Real(-17.0) + Real(12.6)*Real(0.005),
                                    Real(179.9) + Real(20.3)*Real(0.005) - Real(360.0), loc),
              LatLonStatus::Ok);
    // A float holds a longitude near 180 to 1e-5 degrees, about 1 m here
    const Real xtol = (sizeof(Real) == 8) ? Real(1.0e-3) : Real(5.0);
    EXPECT_NEAR(loc.x, Real(20.3)*Real(500.0), xtol);
    EXPECT_NEAR(loc.y, Real(12.6)*Real(500.0), xtol);
}

TEST(LatLonMap, ReportsALongitudeInRangeAcrossTheAntimeridian)
{
    // 40 x 40 cells of 0.005 degrees whose last column sits at 179.9975: a grid
    // that stops just short of the antimeridian and so needs no wrap of its own
    const amrex::Box dom(amrex::IntVect(0,0,0), amrex::IntVect(39,39,0));
    amrex::FArrayBox fab(dom, 2, amrex::The_Pinned_Arena());
    const auto a = fab.array();
    amrex::LoopOnCpu(dom, [&](int i, int j, int) {
        a(i,j,0,0) = Real(-17.0) + (Real(j) + Real(0.5))*Real(0.005);
        a(i,j,0,1) = Real(179.8) + (Real(i) + Real(0.5))*Real(0.005);
    });
    const amrex::GpuArray<Real,AMREX_SPACEDIM> problo{Real(0.0), Real(0.0), Real(0.0)};
    const amrex::GpuArray<Real,AMREX_SPACEDIM> dx{Real(500.0), Real(500.0), Real(100.0)};

    // 0.7 of a cell past the last column: 180.001 degrees east, which is the
    // same place as -179.999
    const Real req_lon = Real(179.8) + Real(40.2)*Real(0.005) - Real(360.0);
    LatLonLocation loc;
    ASSERT_EQ(locate_latlon_on_grid(fab.const_array(), dom, problo, dx,
                                    Real(-17.0) + Real(12.5)*Real(0.005), req_lon, loc),
              LatLonStatus::Ok);

    const Real xtol   = (sizeof(Real) == 8) ? Real(1.0e-3) : Real(5.0);
    const Real lontol = (sizeof(Real) == 8) ? Real(1.0e-9) : Real(1.0e-4);
    EXPECT_NEAR(loc.x, Real(40.2)*Real(500.0), xtol);
    // The nearest grid point is on the near side of 180 and the resolved point
    // on the far side, so the sum runs past +180; the longitude reported for
    // the station has to name a real place
    EXPECT_LE(loc.lon, Real(180.0));
    EXPECT_GE(loc.lon, Real(-180.0));
    EXPECT_NEAR(loc.lon, req_lon, lontol);
}

TEST(LatLonMap, LongitudeDifferencesWrapAtTheAntimeridian)
{
    const Real tol = (sizeof(Real) == 8) ? Real(1.0e-12) : Real(1.0e-4);
    EXPECT_NEAR(wrap_longitude_difference(Real(359.0)),  Real(-1.0), tol);
    EXPECT_NEAR(wrap_longitude_difference(Real(-359.0)), Real( 1.0), tol);
    EXPECT_NEAR(wrap_longitude_difference(Real(10.0)),   Real(10.0), tol);
}
