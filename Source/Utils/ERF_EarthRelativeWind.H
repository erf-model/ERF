/**
 * \file ERF_EarthRelativeWind.H
 */
#ifndef ERF_EARTH_RELATIVE_WIND_H_
#define ERF_EARTH_RELATIVE_WIND_H_

#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include "ERF_NumericalConstants.H"

#include <cmath>

//
// Earth-relative winds from grid-relative ones, shared by the 2D sampled-level
// fields and the station time series so that the two cannot disagree.
//
// (cos_alpha, sin_alpha) is the direction of the grid x axis (increasing i),
// counterclockwise from true east, as LatLonMap and grid_rotation_from_latlon
// return it.  A grid-relative (u, v) wind is then
//
//     u_east  = u cos_alpha - v sin_alpha
//     v_north = u sin_alpha + v cos_alpha
//
namespace erf_earth_wind {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void grid_to_earth (amrex::Real u, amrex::Real v,
                    amrex::Real cos_alpha, amrex::Real sin_alpha,
                    amrex::Real& u_east, amrex::Real& v_north) noexcept
{
    u_east  = u * cos_alpha - v * sin_alpha;
    v_north = u * sin_alpha + v * cos_alpha;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real horizontal_speed (amrex::Real u_east, amrex::Real v_north) noexcept
{
    return std::sqrt(u_east * u_east + v_north * v_north);
}

//
// Meteorological wind direction in degrees, clockwise from north, of where the
// wind comes from, in [0, 360).  A calm wind (speed <= 1e-12 m/s) has no
// direction, and missing_value is returned; the threshold catches an exactly
// zero wind, not a light one.
//
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real meteorological_direction (amrex::Real u_east, amrex::Real v_north,
                                      amrex::Real missing_value) noexcept
{
    if (horizontal_speed(u_east, v_north) <= amrex::Real(1.0e-12)) {
        return missing_value;
    }
    amrex::Real direction = amrex::Real(270.0)
        - std::atan2(v_north, u_east) * amrex::Real(180.0) / PI;
    while (direction < amrex::Real(0.0)) {
        direction += amrex::Real(360.0);
    }
    while (direction >= amrex::Real(360.0)) {
        direction -= amrex::Real(360.0);
    }
    return direction;
}

//
// The orientation of the grid at mass point (i,j) of a lat/lon array ll
// (component 0 latitude, component 1 longitude, both in degrees): the
// direction of increasing i, counterclockwise from true east.  The centred
// difference in i collapses to one side at ilo and ihi, the x extent of the
// domain; a domain one cell wide has no i direction, and gives (1, 0).
//
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void grid_rotation_at (const amrex::Array4<const amrex::Real>& ll, int ilo, int ihi,
                       int i, int j, amrex::Real& cos_alpha, amrex::Real& sin_alpha) noexcept
{
    const int im = (i - 1 > ilo) ? i - 1 : ilo;
    const int ip = (i + 1 < ihi) ? i + 1 : ihi;
    const amrex::Real inv_di = (ip > im) ? amrex::Real(1.0) / amrex::Real(ip - im)
                                         : amrex::Real(0.0);

    // Displacement on the sphere per cell in i: east and north components, in
    // degrees of latitude (the common factor, the Earth's radius, cancels).
    // The longitude difference is taken the short way round the antimeridian.
    amrex::Real dlon = ll(ip,j,0,1) - ll(im,j,0,1);
    if (dlon >  amrex::Real(180.0)) { dlon -= amrex::Real(360.0); }
    if (dlon < -amrex::Real(180.0)) { dlon += amrex::Real(360.0); }
    const amrex::Real lat   = ll(i,j,0,0);
    const amrex::Real north = (ll(ip,j,0,0) - ll(im,j,0,0)) * inv_di;
    const amrex::Real east  = dlon * inv_di * std::cos(lat * PI / amrex::Real(180.0));
    const amrex::Real norm  = std::sqrt(east*east + north*north);

    if (norm > amrex::Real(0.0)) {
        cos_alpha = east  / norm;
        sin_alpha = north / norm;
    } else {
        cos_alpha = amrex::Real(1.0);
        sin_alpha = amrex::Real(0.0);
    }
}

} // namespace erf_earth_wind

#endif
