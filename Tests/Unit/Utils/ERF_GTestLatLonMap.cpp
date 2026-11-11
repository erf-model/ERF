#include <algorithm>
#include <cmath>

#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#include "ERF_EarthRelativeWind.H"
#include "ERF_LatLonMap.H"

// Placing a latitude/longitude on the grid (ERF_LatLonMap), as the station
// time-series output does: the nearest mass point of the lat/lon arrays and one
// local linear solve, with longitude differences taken the short way round; and
// the orientation of the grid there, which the observation nudging uses.

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

// The orientation of the grid, which rotates an earth-relative wind into the
// grid frame for the observation nudging: the i direction of the rotated grid
// is alpha counterclockwise from east, to within the turn of the local east
// direction across the grid (the latitude changes by 0.04 degrees, which turns
// it by a few 1e-4 radians).
TEST(LatLonMap, MeasuresTheGridRotation)
{
    const amrex::Box dom(amrex::IntVect(0,0,0), amrex::IntVect(39,39,0));
    amrex::FArrayBox fab(dom, 2, amrex::The_Pinned_Arena());
    const Real alpha = amrex::Math::pi<Real>()/Real(6.0);
    fill_rotated_latlon(fab, alpha);

    Real c = 0, s = 0;
    grid_rotation_from_latlon(fab.const_array(), dom, 12, 23, c, s);
    EXPECT_NEAR(c, std::cos(alpha), Real(1.0e-3));
    EXPECT_NEAR(s, std::sin(alpha), Real(1.0e-3));

    // The placement reports the rotation at the point it placed
    const amrex::GpuArray<Real,AMREX_SPACEDIM> problo{Real(0.0), Real(0.0), Real(0.0)};
    const amrex::GpuArray<Real,AMREX_SPACEDIM> dx{Real(100.0), Real(100.0), Real(100.0)};
    LatLonLocation loc;
    ASSERT_EQ(locate_latlon_on_grid(fab.const_array(), dom, problo, dx,
                                    fab.const_array()(12,23,0,0), fab.const_array()(12,23,0,1), loc),
              LatLonStatus::Ok);
    EXPECT_NEAR(loc.cos_alpha, std::cos(alpha), Real(1.0e-3));
    EXPECT_NEAR(loc.sin_alpha, std::sin(alpha), Real(1.0e-3));

    // A grid aligned with east and north has no rotation
    fill_rotated_latlon(fab, Real(0.0));
    grid_rotation_from_latlon(fab.const_array(), dom, 12, 23, c, s);
    EXPECT_NEAR(c, Real(1.0), Real(1.0e-6));
    EXPECT_NEAR(s, Real(0.0), Real(1.0e-6));
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

// The earth-relative winds of the 2D sampled-level fields and the stations
// (#4496): a grid-relative wind rotated by the grid orientation, its speed, and
// the direction it comes from, clockwise from north.
TEST(EarthRelativeWind, RotatesAndGivesTheMeteorologicalDirection)
{
    const Real tol = (sizeof(Real) == 8) ? Real(1.0e-12) : Real(1.0e-5);
    const Real alpha = amrex::Math::pi<Real>()/Real(6.0);

    // 10 m/s along a grid x axis 30 degrees north of east: (8.66, 5) earth-relative
    Real u_e = 0, v_n = 0;
    erf_earth_wind::grid_to_earth(Real(10.0), Real(0.0), std::cos(alpha), std::sin(alpha), u_e, v_n);
    EXPECT_NEAR(u_e, Real(10.0)*std::cos(alpha), tol*Real(10.0));
    EXPECT_NEAR(v_n, Real(5.0), tol*Real(10.0));
    // and along the grid y axis: (-5, 8.66)
    erf_earth_wind::grid_to_earth(Real(0.0), Real(10.0), std::cos(alpha), std::sin(alpha), u_e, v_n);
    EXPECT_NEAR(u_e, Real(-5.0), tol*Real(10.0));
    EXPECT_NEAR(v_n, Real(10.0)*std::cos(alpha), tol*Real(10.0));
    EXPECT_NEAR(erf_earth_wind::horizontal_speed(u_e, v_n), Real(10.0), tol*Real(10.0));

    // Wind from the west blows east (270), from the south north (180), from the
    // east west (90), from the north south (0, not 360)
    const Real missing = Real(-999.0);
    EXPECT_NEAR(erf_earth_wind::meteorological_direction(Real( 5.0), Real( 0.0), missing), Real(270.0), tol*Real(360.0));
    EXPECT_NEAR(erf_earth_wind::meteorological_direction(Real( 0.0), Real( 5.0), missing), Real(180.0), tol*Real(360.0));
    EXPECT_NEAR(erf_earth_wind::meteorological_direction(Real(-5.0), Real( 0.0), missing), Real( 90.0), tol*Real(360.0));
    EXPECT_NEAR(erf_earth_wind::meteorological_direction(Real( 0.0), Real(-5.0), missing), Real(  0.0), tol*Real(360.0));
    // From the south-west
    EXPECT_NEAR(erf_earth_wind::meteorological_direction(Real( 1.0), Real( 1.0), missing), Real(225.0), tol*Real(360.0));
    // A calm wind has no direction
    EXPECT_EQ(erf_earth_wind::meteorological_direction(Real(0.0), Real(0.0), missing), missing);
}

// The rotation the 2D plotfile uses (fill_grid_rotation, one value per mass point on
// the level's own boxes) must be the rotation the stations and the nudging use
// (grid_rotation_from_latlon on the whole array), at every point and in particular
// at the box seams, where the centred difference reads a neighbouring box.
TEST(LatLonMap, FillsTheGridRotationOnEveryBoxAsAtOnePoint)
{
    const amrex::Box dom(amrex::IntVect(0,0,0), amrex::IntVect(39,39,0));
    const Real alpha = amrex::Math::pi<Real>()/Real(6.0);

    // The reference: the whole array on one FAB
    amrex::FArrayBox fab(dom, 2, amrex::The_Pinned_Arena());
    fill_rotated_latlon(fab, alpha);

    // The same lat/lon on 25 boxes of 8 x 8, and unequal ones for the output
    amrex::BoxArray ba(dom);
    ba.maxSize(8);
    const amrex::DistributionMapping dm(ba);
    // One ghost cell, filled as the WRF reader fills it: from the whole array,
    // clamped to the domain
    amrex::MultiFab lat(ba, dm, 1, amrex::IntVect(1,1,0), amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
    amrex::MultiFab lon(ba, dm, 1, amrex::IntVect(1,1,0), amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
    const auto ref = fab.const_array();
    for (amrex::MFIter mfi(lat); mfi.isValid(); ++mfi) {
        const auto la = lat.array(mfi);
        const auto lo = lon.array(mfi);
        amrex::LoopOnCpu(mfi.fabbox(), [&](int i, int j, int k) {
            const int ii = std::min(std::max(i, dom.smallEnd(0)), dom.bigEnd(0));
            const int jj = std::min(std::max(j, dom.smallEnd(1)), dom.bigEnd(1));
            la(i,j,k) = ref(ii,jj,0,0);
            lo(i,j,k) = ref(ii,jj,0,1);
        });
    }

    amrex::BoxArray ba_out(dom);
    ba_out.maxSize(amrex::IntVect(13,7,1));
    const amrex::DistributionMapping dm_out(ba_out);
    amrex::MultiFab c(ba_out, dm_out, 1, 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
    amrex::MultiFab s(ba_out, dm_out, 1, 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));

    const amrex::RealBox rb({0.0, 0.0, 0.0}, {4000.0, 4000.0, 100.0});
    const amrex::Array<int,AMREX_SPACEDIM> periodic{0, 0, 0};
    const amrex::Geometry geom(amrex::Box(amrex::IntVect(0,0,0), amrex::IntVect(39,39,0)),
                               rb, amrex::CoordSys::cartesian, periodic);
    fill_grid_rotation(c, s, lat, lon, geom);
    amrex::Gpu::streamSynchronize();

    Real worst = 0;
    for (amrex::MFIter mfi(c); mfi.isValid(); ++mfi) {
        const auto ca = c.const_array(mfi);
        const auto sa = s.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int) {
            Real cr = 0, sr = 0;
            grid_rotation_from_latlon(ref, dom, i, j, cr, sr);
            worst = std::max(worst, std::max(std::abs(ca(i,j,0) - cr), std::abs(sa(i,j,0) - sr)));
        });
    }
    amrex::ParallelDescriptor::ReduceRealMax(worst);
    EXPECT_EQ(worst, Real(0.0));
    EXPECT_NEAR(c.max(0), std::cos(alpha), Real(1.0e-3));
    EXPECT_NEAR(s.min(0), std::sin(alpha), Real(1.0e-3));
}

// On a refined level the boxes cover only a patch of the domain. A neighbour across
// the edge of the patch, inside the domain, is held only in a ghost cell of lat/lon
// (ERF fills them from the coarser level). The rotation at the patch edge must use
// it, not stop at the patch: it must equal the rotation from the whole array.
TEST(LatLonMap, FillsTheGridRotationOnAPatchThatDoesNotCoverTheDomain)
{
    const amrex::Box dom(amrex::IntVect(0,0,0), amrex::IntVect(39,39,0));
    const Real alpha = amrex::Math::pi<Real>()/Real(6.0);
    amrex::FArrayBox fab(dom, 2, amrex::The_Pinned_Arena());
    fill_rotated_latlon(fab, alpha);
    const auto ref = fab.const_array();

    // A patch away from every domain edge, as two boxes, with one ghost cell
    // holding the values of the cells around it
    const amrex::Box patch(amrex::IntVect(10,5,0), amrex::IntVect(29,34,0));
    amrex::BoxArray ba(patch);
    ba.maxSize(amrex::IntVect(10,30,1));
    const amrex::DistributionMapping dm(ba);
    const amrex::MFInfo pinned = amrex::MFInfo().SetArena(amrex::The_Pinned_Arena());
    amrex::MultiFab lat(ba, dm, 1, amrex::IntVect(1,1,0), pinned);
    amrex::MultiFab lon(ba, dm, 1, amrex::IntVect(1,1,0), pinned);
    for (amrex::MFIter mfi(lat); mfi.isValid(); ++mfi) {
        const auto la = lat.array(mfi);
        const auto lo = lon.array(mfi);
        amrex::LoopOnCpu(mfi.fabbox(), [&](int i, int j, int k) {
            la(i,j,k) = ref(i,j,0,0);
            lo(i,j,k) = ref(i,j,0,1);
        });
    }

    amrex::MultiFab c(ba, dm, 1, 0, pinned);
    amrex::MultiFab s(ba, dm, 1, 0, pinned);
    const amrex::RealBox rb({0.0, 0.0, 0.0}, {4000.0, 4000.0, 100.0});
    const amrex::Array<int,AMREX_SPACEDIM> periodic{0, 0, 0};
    const amrex::Geometry geom(dom, rb, amrex::CoordSys::cartesian, periodic);
    fill_grid_rotation(c, s, lat, lon, geom);
    amrex::Gpu::streamSynchronize();

    Real worst = 0;
    for (amrex::MFIter mfi(c); mfi.isValid(); ++mfi) {
        const auto ca = c.const_array(mfi);
        const auto sa = s.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int) {
            Real cr = 0, sr = 0;
            grid_rotation_from_latlon(ref, dom, i, j, cr, sr);
            worst = std::max(worst, std::max(std::abs(ca(i,j,0) - cr), std::abs(sa(i,j,0) - sr)));
        });
    }
    amrex::ParallelDescriptor::ReduceRealMax(worst);
    EXPECT_EQ(worst, Real(0.0));
}
