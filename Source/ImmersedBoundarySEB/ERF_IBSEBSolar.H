#ifndef ERF_IBSEB_SOLAR_H
#define ERF_IBSEB_SOLAR_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Array4.H>
#include <limits>
#include <ERF_NumericalConstants.H>
#include <cmath>

/**
 * \file ERF_IBSEBSolar.H
 * \brief Sun position, clear-sky shortwave, and the ray cast that decides
 *        whether a building face is in shadow.
 *
 * The solar geometry (Spencer's declination and equation of time, the hour
 * angle from UTC and the site longitude, the zenith, the clockwise-from-north
 * azimuth, the Earth-Sun distance) is kept here: the radiation models use
 * the CESM orbital formulation of ERF_OrbCosZenith.H, which gives no
 * azimuth and has no equation of time; orbital_hour_angle() turns its
 * calendar day into the hour angle these functions take, for
 * erf.ibseb.sun_mode = two_stream. The clear-sky formulas follow Bird (1984) as used by the SLUCM
 * branch. The ray cast is new here.
 *
 * **Sun vector convention.** ``ibseb_sun_vector()`` returns the unit vector
 * from a surface toward the sun, with azimuth measured clockwise from north:
 * east is +x, north is +y, so ``s = (sin z sin a, sin z cos a, cos z)`` with
 * z the zenith angle and a the azimuth. A face with outward normal n receives
 * direct beam only when ``n . s > 0``.
 *
 * References: Spencer (1971); Duffie and Beckman (1991); Iqbal (1983);
 * Bird and Hulstrom (1981).
 */
namespace ibseb {

/**
 * Day angle [rad] of Spencer's (1971) series, 2 pi (day - 1) / 365. The
 * series' coefficients were fitted with this day angle; the declination, the
 * equation of time and the Earth-Sun distance below all use it, so day 366
 * gives the same values as day 1.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real spencer_day_angle (amrex::Real day_of_year)
{
    return 2.0 * PI * (day_of_year - 1.0) / 365.0;
}

/**
 * Day of year at a UTC time: the start day plus the whole days elapsed since
 * 00:00 UTC of that day (t_utc_s may exceed 86400 s in a multi-day run). The
 * Spencer series repeat every 365 days, so the result needs no wrapping.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real day_of_year_at (amrex::Real day0, amrex::Real t_utc_s)
{
    return day0 + std::floor(t_utc_s / 86400.0);
}

/** Solar declination [rad] from the day of year (Spencer 1971). */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real solar_declination (amrex::Real day_of_year)
{
    const amrex::Real gamma = spencer_day_angle(day_of_year);
    return 0.006918
        - 0.399912 * std::cos(gamma) + 0.070257 * std::sin(gamma)
        - 0.006758 * std::cos(2.0 * gamma) + 0.000907 * std::sin(2.0 * gamma)
        - 0.002697 * std::cos(3.0 * gamma) + 0.00148 * std::sin(3.0 * gamma);
}

/** Equation of time [minutes] from the day of year (Spencer 1971). */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real equation_of_time (amrex::Real day_of_year)
{
    const amrex::Real gamma = spencer_day_angle(day_of_year);
    return 229.18 * (0.000075 + 0.001868 * std::cos(gamma) - 0.032077 * std::sin(gamma)
                     - 0.014615 * std::cos(2.0 * gamma) - 0.040849 * std::sin(2.0 * gamma));
}

/**
 * Hour angle [rad] from the UTC time of day and the site: solar time is UTC
 * plus one hour per 15 degrees of longitude (east positive) plus the equation
 * of time, and the hour angle is 15 degrees per hour from solar noon,
 * negative in the morning. The time zone plays no part (it would cancel).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real solar_hour_angle (amrex::Real time_utc_s, amrex::Real longitude_deg, amrex::Real day_of_year)
{
    amrex::Real t_utc_h = std::fmod(time_utc_s / 3600.0, 24.0);
    if (t_utc_h < 0.0) { t_utc_h += 24.0; }
    const amrex::Real t_solar_h = t_utc_h + longitude_deg / 15.0 + equation_of_time(day_of_year) / 60.0;
    return 15.0 * (t_solar_h - 12.0) * PI / 180.0;
}

/**
 * Hour angle [rad] of the orbital formula the radiation models use
 * (orbital_cos_zenith_instant() of ERF_OrbCosZenith.H): its
 * ``cos z = sin(lat) sin(decl) - cos(lat) cos(decl) cos(2 pi frac(calday) + lon)``
 * is the usual form with the hour angle ``2 pi frac(calday) + lon - pi``. The
 * fraction of the calendar day is the UTC time of day and there is no equation
 * of time, so local solar time is UTC plus one hour per 15 degrees of
 * longitude. Wrapped to [-pi, pi), negative in the morning.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real orbital_hour_angle (double calday, amrex::Real longitude_rad)
{
    const double two_pi = 2.0 * double(PI);
    double h = two_pi * (calday - std::floor(calday)) + double(longitude_rad) - double(PI);
    h = std::fmod(h + double(PI), two_pi);
    if (h < 0.0) { h += two_pi; }
    return static_cast<amrex::Real>(h - double(PI));
}

/** Zenith angle [rad] from the latitude, declination and hour angle. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real solar_zenith (amrex::Real latitude_deg, amrex::Real decl, amrex::Real hour_angle)
{
    const amrex::Real lat = latitude_deg * PI / 180.0;
    amrex::Real cz = std::sin(lat) * std::sin(decl) + std::cos(lat) * std::cos(decl) * std::cos(hour_angle);
    if (cz < -1.0) { cz = -1.0; }
    if (cz >  1.0) { cz =  1.0; }
    return std::acos(cz);
}

/**
 * Azimuth [rad] clockwise from north (east is pi/2) from the latitude,
 * declination, hour angle and zenith; north when the sun is at the zenith or
 * the site is at a pole, where the azimuth is undefined.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real solar_azimuth (amrex::Real latitude_deg, amrex::Real decl,
                           amrex::Real hour_angle, amrex::Real zenith)
{
    const amrex::Real lat = latitude_deg * PI / 180.0;
    const amrex::Real sz = std::sin(zenith);
    const amrex::Real cz = std::cos(zenith);
    if (sz < 1.0e-6) { return 0.0; }
    const amrex::Real cos_lat = std::cos(lat);
    if (std::abs(cos_lat) < 1.0e-6) { return 0.0; }
    // The hour angle is negative before solar noon, when the sun is east of
    // south, hence the sign for the clockwise-from-north convention.
    const amrex::Real sin_az = -std::sin(hour_angle) * std::cos(decl) / sz;
    const amrex::Real cos_az = (std::sin(decl) - cz * std::sin(lat)) / (sz * cos_lat);
    amrex::Real az = std::atan2(sin_az, cos_az);
    if (az < 0.0) { az += 2.0 * PI; }
    return az;
}

/**
 * Earth-Sun distance factor (d0/d)^2 that scales the solar constant for the
 * day of year (Spencer 1971): about 1.034 at perihelion in early January,
 * 0.967 at aphelion in early July, annual mean 1.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real earth_sun_distance_factor (amrex::Real day_of_year)
{
    const amrex::Real G = spencer_day_angle(day_of_year);
    return 1.000110 + 0.034221 * std::cos(G) + 0.001280 * std::sin(G)
         + 0.000719 * std::cos(2.0 * G) + 0.000077 * std::sin(2.0 * G);
}

/** Unit vector toward the sun from zenith and azimuth (east +x, north +y). */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void sun_vector (amrex::Real zenith, amrex::Real azimuth,
                 amrex::Real& sx, amrex::Real& sy, amrex::Real& sz)
{
    const amrex::Real s = std::sin(zenith);
    sx = s * std::sin(azimuth);
    sy = s * std::cos(azimuth);
    sz = std::cos(zenith);
}

/**
 * Clear-sky direct-normal irradiance [W/m2] (Bird form): the solar constant,
 * corrected for the Earth-Sun distance, attenuated by ``tau^(1/cos z)``.
 * Zero when the sun is below the horizon.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real clear_sky_dni (amrex::Real cos_zenith, amrex::Real S0, amrex::Real tau,
                           amrex::Real distance_factor)
{
    if (cos_zenith <= 1.0e-3) return 0.0;
    return S0 * distance_factor * std::pow(tau, 1.0 / cos_zenith);
}

/**
 * Clear-sky diffuse irradiance on a horizontal surface [W/m2]: a fraction
 * ``k_d`` of what the direct beam lost to the atmosphere reaches the ground
 * as diffuse light (Liu and Jordan form).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real clear_sky_diffuse_h (amrex::Real cos_zenith, amrex::Real S0, amrex::Real tau,
                                 amrex::Real distance_factor, amrex::Real k_d)
{
    if (cos_zenith <= 1.0e-3) return 0.0;
    const amrex::Real tr = std::pow(tau, 1.0 / cos_zenith);
    return k_d * S0 * distance_factor * cos_zenith * (1.0 - tr);
}

/**
 * The column tops a ray cast on one rank sees: a k-collapsed map over this
 * rank's patch of the level (its boxes grown by the halo of
 * IBFaceSet::build_ray_map()), with the domain and the grid it is drawn on.
 * An entry is the domain-relative k of the highest solid cell of its column,
 * -1 where the column is fluid; column_top() turns it into a height. The
 * patch holds the columns of other ranks and the periodic images within the
 * halo, so the walk needs no communication. Where ``wrap_x`` (``wrap_y``) is
 * set the patch spans the whole periodic domain in x (y) and a ray that
 * leaves it comes back in on the other side; otherwise a ray that leaves the
 * patch has left the map: the far field is open ground (see locate()). A
 * ray is followed for as many columns as the halo's reach holds
 * (``max_steps``), which is every column a ray no shallower than the map's
 * floor can cross before it clears the tallest column or reaches the
 * ground; a shallower ray stops there too, so a building beyond the halo's
 * reach never shades, whether the map wraps or not. A level ray is followed
 * for two rounds of the map and then counts as sky.
 */
struct RayMap
{
    amrex::Array4<const int> top;       ///< Column tops over the patch (any k index; use ``top.begin[2]``)
    int  ilo = 0, jlo = 0;              ///< Low corner of the domain, in cells
    int  nx = 0, ny = 0;                ///< Columns of the domain in x and y
    bool wrap_x = false, wrap_y = false;///< The patch spans the periodic domain in x / y
    amrex::Real x_lo = 0.0, y_lo = 0.0; ///< Physical position of the domain's low corner [m]
    amrex::Real dx = 0.0, dy = 0.0, dz = 0.0;  ///< Cell sizes [m]
    amrex::Real z_ground = 0.0;         ///< Ground height [m]
    amrex::Real z_max = 0.0;            ///< Tallest column the map holds [m]; a ray above it meets nothing
    int  max_steps = 0;                 ///< Columns a ray at the map's floor elevation crosses before it clears z_max (the halo's reach in x plus y)

    /** Whether the map holds column ``(ci, cj)`` (wrapped into the patch where the patch spans the domain). */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool locate (int& ci, int& cj) const noexcept
    {
        if (ci < top.begin[0] || ci >= top.end[0]) {
            if (!wrap_x) { return false; }
            ci = top.begin[0] + ((ci - top.begin[0]) % nx + nx) % nx;
        }
        if (cj < top.begin[1] || cj >= top.end[1]) {
            if (!wrap_y) { return false; }
            cj = top.begin[1] + ((cj - top.begin[1]) % ny + ny) % ny;
        }
        return true;
    }

    /** Top height of column ``(ci, cj)`` of the map (the ground where it is fluid). */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real column_top (int ci, int cj) const noexcept
    {
        const int kt = top(ci, cj, top.begin[2]);
        return (kt < 0) ? z_ground : z_ground + (kt + 1) * dz;
    }
};

/**
 * Halo of a ray map, in cells: how far a ray no shallower than ``sz_min``
 * (the sine of its elevation) can travel horizontally before it rises above
 * ``height`` or falls that far, ``height / tan(elevation)``, rounded up and
 * with one cell for the column it starts in. At least one cell, so a face
 * always finds the column of its solid neighbour; a level ray would reach
 * for ever, which the caller caps at the domain.
 */
inline int ray_halo_cells (amrex::Real height, amrex::Real dx, amrex::Real sz_min)
{
    if (height <= 0.0) { return 1; }
    const amrex::Real cz = std::sqrt(amrex::max(amrex::Real(0.0), amrex::Real(1.0) - sz_min * sz_min));
    if (sz_min <= 0.0 || cz <= 0.0) { return std::numeric_limits<int>::max() / 4; }
    const amrex::Real reach = height * cz / sz_min / dx;
    if (reach >= amrex::Real(std::numeric_limits<int>::max() / 8)) { return std::numeric_limits<int>::max() / 4; }
    return amrex::max(1, static_cast<int>(std::ceil(reach)) + 1);
}

/**
 * Sine of the elevation of the shallowest ray of the hemisphere sampling
 * (hemisphere_direction()) with ``n_az`` azimuths and ``n_el`` elevations:
 * on a roof the ray furthest from the normal, ``cos(theta_max)``; on a wall,
 * where a ray's sine of elevation is ``sin(theta) sin(phi)``, the ray
 * nearest the normal at the azimuth nearest the horizontal,
 * ``sin(theta_min) sin(pi / n_az)``. With an odd ``n_az`` a wall has level
 * rays (``phi = pi``), which the walk follows for two rounds of the map.
 */
inline amrex::Real hemisphere_min_sz (int n_az, int n_el)
{
    const amrex::Real u_min = 0.5 / n_el, u_max = (n_el - 0.5) / n_el;
    const amrex::Real roof  = std::sqrt(1.0 - u_max);
    const amrex::Real wall  = std::sqrt(u_min) * std::sin(PI / n_az);
    return amrex::min(roof, wall);
}

/**
 * Whether the ray from a face toward the sun hits a building.
 *
 * Buildings are the column tops of the rank's ray map (RayMap), so a ray is
 * blocked wherever its height on entering a column is below that column's
 * top: buildings stand on the ground, so the lowest point of the ray inside
 * a column is at its entry, and the height only grows along a ray toward a
 * sun above the horizon. The walk over the columns is a 2D digital
 * differential analyser. It stops when the ray rises above the tallest
 * column or leaves the map (the far field is open ground; in a periodic
 * direction the map spans, the ray wraps round instead).
 *
 * The start point is nudged a little along the ray so a face does not
 * test its own column at its own height: a wall starts in its fluid cell's
 * column, whose top is below the face (no overhangs), and a roof starts a
 * hair above its own top.
 *
 * @return true when blocked (the face is in shadow).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool ray_blocked (amrex::Real x0, amrex::Real y0, amrex::Real z0,
                  amrex::Real sx, amrex::Real sy, amrex::Real sz, const RayMap& m)
{
    if (sz <= 0.0) return true;              // sun at or below the horizon
    if (z0 >= m.z_max) return false;         // nothing stands this high
    const amrex::Real dx = m.dx, dy = m.dy;
    const amrex::Real eps = 1.0e-6 * amrex::min(dx, dy);
    const amrex::Real t0  = 1.0e-3 * amrex::min(dx, dy);
    amrex::Real x = x0 + t0 * sx, y = y0 + t0 * sy, z = z0 + t0 * sz;
    int ci = m.ilo + static_cast<int>(std::floor((x - m.x_lo) / dx));
    int cj = m.jlo + static_cast<int>(std::floor((y - m.y_lo) / dy));
    const int stepx = (sx > 0.0) ? 1 : ((sx < 0.0) ? -1 : 0);
    const int stepy = (sy > 0.0) ? 1 : ((sy < 0.0) ? -1 : 0);
    const amrex::Real big = 1.0e30;
    const amrex::Real tdx = (stepx != 0) ? dx / std::abs(sx) : big;
    const amrex::Real tdy = (stepy != 0) ? dy / std::abs(sy) : big;
    // Distance along the ray to the next column boundary in x and in y.
    amrex::Real tmx = big, tmy = big;
    if (stepx > 0)      tmx = ((m.x_lo + (ci - m.ilo + 1) * dx) - x) / sx;
    else if (stepx < 0) tmx = ((m.x_lo + (ci - m.ilo) * dx) - x) / sx;
    if (stepy > 0)      tmy = ((m.y_lo + (cj - m.jlo + 1) * dy) - y) / sy;
    else if (stepy < 0) tmy = ((m.y_lo + (cj - m.jlo) * dy) - y) / sy;
    // Columns crossed before the ray passes the tallest one, plus a margin,
    // and never more than the halo's reach (see RayMap).
    const amrex::Real t_max = (m.z_max - z) / sz;
    const amrex::Real steps = amrex::min(t_max * (std::abs(sx) / dx + std::abs(sy) / dy), amrex::Real(m.max_steps));
    const int max_iter = static_cast<int>(amrex::min(amrex::Real(1.0e9), steps)) + 8;
    for (int it = 0; it < max_iter; ++it) {
        // Column (ci, cj) with the ray entering at height z.
        if (!m.locate(ci, cj)) return false;
        if (z + eps < m.column_top(ci, cj)) return true;
        // Step to the next column.
        amrex::Real t;
        if (tmx < tmy) { t = tmx; tmx += tdx; ci += stepx; }
        else           { t = tmy; tmy += tdy; cj += stepy; }
        z = z0 + t0 * sz + t * sz;
        if (z >= m.z_max) return false;
    }
    return false;
}

/// What a ray from a face ends on.
enum RayHit : int { RAY_SKY = 0, RAY_GROUND = 1, RAY_BUILDING = 2 };

/**
 * Where a ray from a face ends: sky, ground or a building, for any
 * direction (the hemisphere sampling of the view fractions). The same column walk as
 * ray_blocked(), with the height test made direction-aware:
 *
 * - rising or level rays (``sz >= 0``) are blocked by a column whose top is
 *   above the entry height; a level ray that leaves the map, or has crossed
 *   the domain twice where the map wraps, ends at the horizon, which counts as sky;
 * - falling rays (``sz < 0``) are blocked by a solid column whose top is
 *   above the exit height (the ray descends inside the column) and reach
 *   the ground where the exit height is at or below ``z_ground``.
 *
 * Solid columns are those with a top above the ground. Leaving the map
 * counts as sky for a rising ray and as ground for a falling one: the far
 * field is open ground. Where the map spans a periodic direction the ray
 * wraps round instead.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int ray_hit (amrex::Real x0, amrex::Real y0, amrex::Real z0,
             amrex::Real sx, amrex::Real sy, amrex::Real sz, const RayMap& m)
{
    const amrex::Real dx = m.dx, dy = m.dy;
    const amrex::Real eps = 1.0e-6 * amrex::min(dx, dy);
    if (sz > 0.0 && z0 >= m.z_max) return RAY_SKY;
    const amrex::Real t0 = 1.0e-3 * amrex::min(dx, dy);
    amrex::Real x = x0 + t0 * sx, y = y0 + t0 * sy, z = z0 + t0 * sz;
    int ci = m.ilo + static_cast<int>(std::floor((x - m.x_lo) / dx));
    int cj = m.jlo + static_cast<int>(std::floor((y - m.y_lo) / dy));
    const int stepx = (sx > 0.0) ? 1 : ((sx < 0.0) ? -1 : 0);
    const int stepy = (sy > 0.0) ? 1 : ((sy < 0.0) ? -1 : 0);
    const amrex::Real big = 1.0e30;
    const amrex::Real tdx = (stepx != 0) ? dx / std::abs(sx) : big;
    const amrex::Real tdy = (stepy != 0) ? dy / std::abs(sy) : big;
    amrex::Real tmx = big, tmy = big;
    if (stepx > 0)      tmx = ((m.x_lo + (ci - m.ilo + 1) * dx) - x) / sx;
    else if (stepx < 0) tmx = ((m.x_lo + (ci - m.ilo) * dx) - x) / sx;
    if (stepy > 0)      tmy = ((m.y_lo + (cj - m.jlo + 1) * dy) - y) / sy;
    else if (stepy < 0) tmy = ((m.y_lo + (cj - m.jlo) * dy) - y) / sy;
    if (stepx == 0 && stepy == 0) {          // vertical ray
        if (sz > 0.0) return RAY_SKY;
        // straight down from a wall face never happens (the face is vertical);
        // from a roof the hemisphere points up. Treat as ground for safety.
        return RAY_GROUND;
    }
    // Columns crossed before the ray passes the tallest column or reaches
    // the ground, plus a margin, and never more than the halo's reach (see
    // RayMap). A ray within 1e-9 of level (phi = pi of an odd n_az gives sz
    // of 1e-16) is a level ray and gets two rounds of the map.
    const bool level = (std::abs(sz) < 1.0e-9);
    const amrex::Real t_max = (sz > 0.0) ? (m.z_max - z) / sz : ((sz < 0.0) ? (z - m.z_ground) / (-sz) : big);
    const amrex::Real rounds = amrex::Real(2 * ((m.top.end[0] - m.top.begin[0]) + (m.top.end[1] - m.top.begin[1])));
    const amrex::Real steps = level ? rounds
                            : amrex::min(t_max * (std::abs(sx) / dx + std::abs(sy) / dy), amrex::Real(m.max_steps));
    const int max_iter = static_cast<int>(amrex::min(amrex::Real(1.0e9), steps)) + 8;
    for (int it = 0; it < max_iter; ++it) {
        if (!m.locate(ci, cj)) return (sz >= 0.0) ? RAY_SKY : RAY_GROUND;
        const amrex::Real top = m.column_top(ci, cj);
        const bool solid = (top > m.z_ground + eps);
        const amrex::Real t_next = amrex::min(tmx, tmy);
        const amrex::Real z_out  = z0 + t0 * sz + t_next * sz;
        if (sz >= 0.0) {
            if (solid && z + eps < top) return RAY_BUILDING;
        } else {
            if (solid && z + eps < top) return RAY_BUILDING;     // entered below the top
            if (solid && z_out < top - eps) return RAY_BUILDING; // descends into it
            if (z_out <= m.z_ground + eps) return RAY_GROUND;
        }
        if (tmx < tmy) { tmx += tdx; ci += stepx; }
        else           { tmy += tdy; cj += stepy; }
        z = z_out;
        if (sz > 0.0 && z >= m.z_max) return RAY_SKY;
    }
    return (sz >= 0.0) ? RAY_SKY : RAY_GROUND;
}

/**
 * Direction of the ``(ia, ie)`` sample of a cosine-weighted hemisphere
 * around the outward normal of a face with direction ``dir`` and normal
 * sign ``nsign`` (+1 along the axis, -1 against it). Stratified in the
 * azimuth ``phi = 2 pi (ia + 1/2) / n_az`` and in ``u = (ie + 1/2) / n_el``
 * with ``theta = asin(sqrt(u))`` from the normal, so every sample carries
 * the same weight and the counts are view factors.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void hemisphere_direction (int dir, int nsign, int ia, int ie, int n_az, int n_el,
                           amrex::Real& sx, amrex::Real& sy, amrex::Real& sz)
{
    const amrex::Real u     = (ie + 0.5) / n_el;
    const amrex::Real theta = std::asin(std::sqrt(u));
    const amrex::Real phi   = 2.0 * PI * (ia + 0.5) / n_az;
    const amrex::Real cn = std::cos(theta), ct = std::sin(theta) * std::cos(phi), st = std::sin(theta) * std::sin(phi);
    // Normal and two tangents for each axis.
    amrex::Real n[3] = {0.0, 0.0, 0.0}, t1[3] = {0.0, 0.0, 0.0}, t2[3] = {0.0, 0.0, 0.0};
    n[dir] = nsign;
    if (dir == 0)      { t1[1] = 1.0; t2[2] = 1.0; }
    else if (dir == 1) { t1[0] = 1.0; t2[2] = 1.0; }
    else               { t1[0] = 1.0; t2[1] = 1.0; }
    sx = cn * n[0] + ct * t1[0] + st * t2[0];
    sy = cn * n[1] + ct * t1[1] + st * t2[1];
    sz = cn * n[2] + ct * t1[2] + st * t2[2];
}

} // namespace ibseb

#endif
