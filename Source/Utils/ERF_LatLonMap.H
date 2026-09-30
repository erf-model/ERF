/**
 * \file ERF_LatLonMap.H
 */
#ifndef ERF_LATLONMAP_H_
#define ERF_LATLONMAP_H_

#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <algorithm>
#include <cmath>
#include <limits>

//
// Placing a point given by latitude and longitude on the model grid.
//
// The model has no map projection of its own: a run initialized from a WRF or
// metgrid file carries the latitude and longitude of every level-0 mass point
// (lat_m / lon_m), whatever projection produced them.  A requested point is
// found by searching those arrays for the nearest mass point and then inverting
// the map locally with one linear solve, which is exact for a map that is
// linear over a cell and well under a cell width in error for any smooth one.
// Longitude differences are taken the short way round, so a grid that crosses
// the antimeridian is handled like any other.
//
// The same arrays give the orientation of the grid: the direction in which the
// i index grows, measured counterclockwise from true east, is the angle alpha
// by which an earth-relative (east, north) wind has to be rotated to become a
// grid-relative (u, v) wind.
//

//
// Invert the latitude/longitude map locally.  Given the difference between the
// requested lat/lon and the lat/lon of a grid point, and the gradients of
// latitude and longitude with respect to the two horizontal index directions at
// that point, return the offset in cells from the grid point to the request.
//
// Returns false when the two gradients are parallel and the map is therefore not
// locally invertible.  "Parallel" is judged relative to the gradients' size and
// to the precision: in single precision the determinant of an exactly
// degenerate map comes out of rounding at about 1e-8 of its scale squared, far
// above a fixed 1e-12, so the threshold is at least 64 machine epsilons (which
// leaves double precision at 1e-12).
//
inline bool
invert_latlon_offset (amrex::Real dlat, amrex::Real dlon,
                      amrex::Real dlat_di, amrex::Real dlon_di,
                      amrex::Real dlat_dj, amrex::Real dlon_dj,
                      amrex::Real& di, amrex::Real& dj)
{
    const amrex::Real scale = std::max(std::max(std::abs(dlat_di), std::abs(dlon_di)),
                                       std::max(std::abs(dlat_dj), std::abs(dlon_dj)));
    const amrex::Real det = dlat_di*dlon_dj - dlat_dj*dlon_di;
    const amrex::Real rel = std::max(amrex::Real(1.0e-12),
                                     amrex::Real(64) * std::numeric_limits<amrex::Real>::epsilon());
    if (scale <= amrex::Real(0.0) || std::abs(det) <= rel*scale*scale) {
        return false;
    }
    di = ( dlat*dlon_dj - dlon*dlat_dj) / det;
    dj = (-dlat*dlon_di + dlon*dlat_di) / det;
    return true;
}

//
// A longitude difference in degrees, brought into [-180, 180] so that two
// points either side of the antimeridian are close rather than 360 degrees apart.
//
inline amrex::Real
wrap_longitude_difference (amrex::Real dlon)
{
    if (dlon >  amrex::Real(180.0)) { return dlon - amrex::Real(360.0); }
    if (dlon < -amrex::Real(180.0)) { return dlon + amrex::Real(360.0); }
    return dlon;
}

//
// The result of placing a latitude/longitude on the grid.
//
struct LatLonLocation
{
    // Domain coordinates of the point
    amrex::Real x = amrex::Real(0.0);
    amrex::Real y = amrex::Real(0.0);

    // The lat/lon actually reached by (x,y): equals the request to within the
    // accuracy of the local inverse
    amrex::Real lat = amrex::Real(0.0);
    amrex::Real lon = amrex::Real(0.0);

    // The nearest mass point, reported when the request cannot be placed
    amrex::Real near_lat = amrex::Real(0.0);
    amrex::Real near_lon = amrex::Real(0.0);

    // Orientation of the grid x axis, counterclockwise from true east
    amrex::Real cos_alpha = amrex::Real(1.0);
    amrex::Real sin_alpha = amrex::Real(0.0);
};

enum class LatLonStatus {
    Ok,         // placed
    Degenerate, // the lat/lon gradients are parallel near the point
    TooFar      // more than one cell from the nearest mass point: outside the domain
};

//
// The orientation of the grid at mass point (i,j) of a lat/lon array ll
// (component 0 latitude, component 1 longitude, both in degrees) over the 2D
// domain dom: the direction of increasing i, counterclockwise from true east.
//
void grid_rotation_from_latlon (const amrex::Array4<const amrex::Real>& ll,
                                const amrex::Box& dom, int i, int j,
                                amrex::Real& cos_alpha, amrex::Real& sin_alpha);

//
// Place (req_lat, req_lon) on the grid described by the lat/lon array ll over
// the 2D domain dom, whose cell (i,j) is centred at problo + (i+1/2) dx.  Host
// code; the whole array must be addressable.
//
LatLonStatus locate_latlon_on_grid (const amrex::Array4<const amrex::Real>& ll,
                                    const amrex::Box& dom,
                                    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& problo,
                                    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx,
                                    amrex::Real req_lat, amrex::Real req_lon,
                                    LatLonLocation& loc);

// Copy lat and lon (2D, one component each) into the two components of latlon,
// all three on the same BoxArray and DistributionMapping
void pack_latlon (amrex::MultiFab& latlon, const amrex::MultiFab& lat, const amrex::MultiFab& lon);

//
// The level-0 lat/lon arrays gathered onto the I/O rank, with the two
// questions asked of them answered there and broadcast.  Every member function
// is collective: call it on every rank with the same arguments.
//
class LatLonMap
{
public:
    // lat and lon are the level-0 mass-point arrays (lat_m[0], lon_m[0]),
    // defined on ba2d/dm; geom0 is the level-0 geometry.
    LatLonMap (const amrex::MultiFab& lat, const amrex::MultiFab& lon,
               const amrex::BoxArray& ba2d, const amrex::DistributionMapping& dm,
               const amrex::Geometry& geom0);

    // Place a latitude/longitude on the grid
    LatLonStatus locate (amrex::Real req_lat, amrex::Real req_lon, LatLonLocation& loc) const;

    // Orientation of the grid at domain coordinates (x, y), from the cell that contains it
    void rotation_at (amrex::Real x, amrex::Real y,
                      amrex::Real& cos_alpha, amrex::Real& sin_alpha) const;

private:
    amrex::MultiFab m_ll;   // one box on the I/O rank, lat and lon
    amrex::Box m_dom;       // the 2D level-0 domain
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_problo;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_dx;
};

#endif
