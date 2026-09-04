#ifndef ERF_IBSEB_SOLAR_H
#define ERF_IBSEB_SOLAR_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <cmath>

/**
 * \file ERF_IBSEBSolar.H
 * \brief Sun position, clear-sky shortwave, and the ray cast that decides
 *        whether a building face is in shadow.
 *
 * The solar geometry (declination, equation of time, hour angle, zenith,
 * azimuth, Earth-Sun distance) is copied from the ERF-Radiation branch's
 * ERF_SolarGeometry.H so that this branch carries no dependency on it;
 * when that branch is merged the two should be unified. The clear-sky
 * formulas follow Bird (1984) as used by the SLUCM branch. The ray cast is
 * new here.
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real clamp_day_of_year (amrex::Real doy)
{
    if (!(doy == doy)) return 172.0;
    if (doy < 1.0) return 1.0;
    if (doy > 366.0) return 366.0;
    return doy;
}

/** Solar declination [rad] from the day of year (Spencer 1971). */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real solar_declination (amrex::Real day_of_year)
{
    const amrex::Real g = 2.0 * M_PI * (clamp_day_of_year(day_of_year) - 1.0) / 365.25;
    return 0.006918
        - 0.399912 * std::cos(g)       + 0.070257 * std::sin(g)
        - 0.006758 * std::cos(2.0 * g) + 0.000907 * std::sin(2.0 * g)
        - 0.002480 * std::cos(3.0 * g) + 0.000310 * std::sin(3.0 * g);
}

/** Equation of time [minutes] (Spencer 1971). */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real equation_of_time_min (amrex::Real day_of_year)
{
    const amrex::Real g = 2.0 * M_PI * (clamp_day_of_year(day_of_year) - 1.0) / 365.25;
    return 229.18 * (0.000075 + 0.001868 * std::cos(g) - 0.032077 * std::sin(g)
                     - 0.014615 * std::cos(2.0 * g) - 0.040849 * std::sin(2.0 * g));
}

/**
 * Solar hour angle [rad], zero at solar noon, negative in the morning.
 * @param time_utc_s   Seconds after 00:00 UTC of the day.
 * @param longitude_deg Site longitude, east positive.
 * @param tz_offset_h  Time-zone offset from UTC [h]; the longitude correction
 *                     is taken relative to that zone's standard meridian.
 * @param day_of_year  For the equation of time.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real solar_hour_angle (amrex::Real time_utc_s, amrex::Real longitude_deg,
                              amrex::Real tz_offset_h, amrex::Real day_of_year)
{
    amrex::Real t_utc_h = std::fmod(time_utc_s / 3600.0, 24.0);
    if (t_utc_h < 0.0) t_utc_h += 24.0;
    const amrex::Real std_meridian = 15.0 * tz_offset_h;
    const amrex::Real t_solar_h = t_utc_h + tz_offset_h + (longitude_deg - std_meridian) / 15.0
                                + equation_of_time_min(day_of_year) / 60.0;
    return 15.0 * (t_solar_h - 12.0) * M_PI / 180.0;
}

/** Solar zenith angle [rad] from latitude, declination and hour angle. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real solar_zenith (amrex::Real latitude_deg, amrex::Real decl, amrex::Real hour_angle)
{
    const amrex::Real lat = latitude_deg * M_PI / 180.0;
    amrex::Real cz = std::sin(lat) * std::sin(decl) + std::cos(lat) * std::cos(decl) * std::cos(hour_angle);
    cz = amrex::max(amrex::Real(-1.0), amrex::min(amrex::Real(1.0), cz));
    return std::acos(cz);
}

/** Solar azimuth [rad], clockwise from north, in [0, 2 pi). */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real solar_azimuth (amrex::Real latitude_deg, amrex::Real decl,
                           amrex::Real hour_angle, amrex::Real zenith)
{
    const amrex::Real lat = latitude_deg * M_PI / 180.0;
    const amrex::Real sz = std::sin(zenith), cz = std::cos(zenith);
    if (sz < 1.0e-6) return 0.0;
    const amrex::Real sin_az = std::sin(hour_angle) * std::cos(decl) / sz;
    const amrex::Real cos_az = (std::sin(decl) - cz * std::sin(lat)) / (sz * std::cos(lat));
    amrex::Real az = std::atan2(sin_az, cos_az);
    if (az < 0.0) az += 2.0 * M_PI;
    return az;
}

/** Earth-Sun distance correction factor (Spencer 1971), about 0.97 to 1.03. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real earth_sun_distance_factor (amrex::Real day_of_year)
{
    const amrex::Real g = 2.0 * M_PI * (clamp_day_of_year(day_of_year) - 1.0) / 365.0;
    return 1.000110 + 0.034221 * std::cos(g) + 0.001280 * std::sin(g)
         + 0.000719 * std::cos(2.0 * g) + 0.000077 * std::sin(2.0 * g);
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
 * Whether the ray from a face toward the sun hits a building.
 *
 * Buildings are represented by the top height of every column of the level
 * (``col_top``, replicated on every rank, ``nx * ny`` values indexed
 * ``ci * ny + cj``), so a ray is blocked wherever its height on entering a
 * column is below that column's top: buildings stand on the ground, so the
 * lowest point of the ray inside a column is at its entry, and the height
 * only grows along a ray toward a sun above the horizon. The walk over the
 * columns is a 2D digital differential analyser. It stops when the ray
 * rises above the tallest column, leaves the domain in a non-periodic
 * direction, or has travelled ``max_path``.
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
                  amrex::Real sx, amrex::Real sy, amrex::Real sz,
                  const amrex::Real* col_top, int nx, int ny,
                  amrex::Real x_lo, amrex::Real y_lo, amrex::Real dx, amrex::Real dy,
                  bool per_x, bool per_y, amrex::Real z_max, amrex::Real max_path)
{
    if (sz <= 0.0) return true;              // sun at or below the horizon
    if (z0 >= z_max) return false;           // nothing stands this high
    const amrex::Real eps = 1.0e-6 * amrex::min(dx, dy);
    const amrex::Real t0  = 1.0e-3 * amrex::min(dx, dy);
    amrex::Real x = x0 + t0 * sx, y = y0 + t0 * sy, z = z0 + t0 * sz;
    int ci = static_cast<int>(std::floor((x - x_lo) / dx));
    int cj = static_cast<int>(std::floor((y - y_lo) / dy));
    const int stepx = (sx > 0.0) ? 1 : ((sx < 0.0) ? -1 : 0);
    const int stepy = (sy > 0.0) ? 1 : ((sy < 0.0) ? -1 : 0);
    const amrex::Real big = 1.0e30;
    const amrex::Real tdx = (stepx != 0) ? dx / std::abs(sx) : big;
    const amrex::Real tdy = (stepy != 0) ? dy / std::abs(sy) : big;
    // Distance along the ray to the next column boundary in x and in y.
    amrex::Real tmx = big, tmy = big;
    if (stepx > 0)      tmx = ((x_lo + (ci + 1) * dx) - x) / sx;
    else if (stepx < 0) tmx = ((x_lo + ci * dx) - x) / sx;
    if (stepy > 0)      tmy = ((y_lo + (cj + 1) * dy) - y) / sy;
    else if (stepy < 0) tmy = ((y_lo + cj * dy) - y) / sy;
    amrex::Real t = 0.0;
    const int max_iter = 4 * (nx + ny) + 8;
    for (int it = 0; it < max_iter; ++it) {
        // Column (ci, cj) with the ray entering at height z.
        if (ci < 0 || ci >= nx) {
            if (!per_x) return false;
            ci = (ci % nx + nx) % nx;
        }
        if (cj < 0 || cj >= ny) {
            if (!per_y) return false;
            cj = (cj % ny + ny) % ny;
        }
        if (z + eps < col_top[ci * ny + cj]) return true;
        // Step to the next column.
        if (tmx < tmy) { t = tmx; tmx += tdx; ci += stepx; }
        else           { t = tmy; tmy += tdy; cj += stepy; }
        z = z0 + t0 * sz + t * sz;
        if (z >= z_max) return false;
        if (t > max_path) return false;
    }
    return false;
}

/// What a ray from a face ends on.
enum RayHit : int { RAY_SKY = 0, RAY_GROUND = 1, RAY_BUILDING = 2 };

/**
 * Where a ray from a face ends: sky, ground or a building, for any
 * direction (phase 3 hemisphere sampling). The same column walk as
 * ray_blocked(), with the height test made direction-aware:
 *
 * - rising or level rays (``sz >= 0``) are blocked by a column whose top is
 *   above the entry height; a level ray that crosses ``max_path`` without a
 *   hit ends at the horizon, which counts as sky;
 * - falling rays (``sz < 0``) are blocked by a solid column whose top is
 *   above the exit height (the ray descends inside the column) and reach
 *   the ground where the exit height is at or below ``z_ground``.
 *
 * Solid columns are those with a top above the ground. The domain is
 * wrapped in periodic directions; leaving it in a non-periodic direction
 * counts as sky for a rising ray and as ground for a falling one.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int ray_hit (amrex::Real x0, amrex::Real y0, amrex::Real z0,
             amrex::Real sx, amrex::Real sy, amrex::Real sz,
             const amrex::Real* col_top, int nx, int ny,
             amrex::Real x_lo, amrex::Real y_lo, amrex::Real dx, amrex::Real dy,
             bool per_x, bool per_y, amrex::Real z_ground, amrex::Real z_max,
             amrex::Real max_path)
{
    const amrex::Real eps = 1.0e-6 * amrex::min(dx, dy);
    if (sz > 0.0 && z0 >= z_max) return RAY_SKY;
    const amrex::Real t0 = 1.0e-3 * amrex::min(dx, dy);
    amrex::Real x = x0 + t0 * sx, y = y0 + t0 * sy, z = z0 + t0 * sz;
    int ci = static_cast<int>(std::floor((x - x_lo) / dx));
    int cj = static_cast<int>(std::floor((y - y_lo) / dy));
    const int stepx = (sx > 0.0) ? 1 : ((sx < 0.0) ? -1 : 0);
    const int stepy = (sy > 0.0) ? 1 : ((sy < 0.0) ? -1 : 0);
    const amrex::Real big = 1.0e30;
    const amrex::Real tdx = (stepx != 0) ? dx / std::abs(sx) : big;
    const amrex::Real tdy = (stepy != 0) ? dy / std::abs(sy) : big;
    amrex::Real tmx = big, tmy = big;
    if (stepx > 0)      tmx = ((x_lo + (ci + 1) * dx) - x) / sx;
    else if (stepx < 0) tmx = ((x_lo + ci * dx) - x) / sx;
    if (stepy > 0)      tmy = ((y_lo + (cj + 1) * dy) - y) / sy;
    else if (stepy < 0) tmy = ((y_lo + cj * dy) - y) / sy;
    if (stepx == 0 && stepy == 0) {          // vertical ray
        if (sz > 0.0) return RAY_SKY;
        // straight down from a wall face never happens (the face is vertical);
        // from a roof the hemisphere points up. Treat as ground for safety.
        return RAY_GROUND;
    }
    amrex::Real t = 0.0;
    const int max_iter = 4 * (nx + ny) + 8;
    for (int it = 0; it < max_iter; ++it) {
        if (ci < 0 || ci >= nx) {
            if (!per_x) return (sz >= 0.0) ? RAY_SKY : RAY_GROUND;
            ci = (ci % nx + nx) % nx;
        }
        if (cj < 0 || cj >= ny) {
            if (!per_y) return (sz >= 0.0) ? RAY_SKY : RAY_GROUND;
            cj = (cj % ny + ny) % ny;
        }
        const amrex::Real top = col_top[ci * ny + cj];
        const bool solid = (top > z_ground + eps);
        const amrex::Real t_next = amrex::min(tmx, tmy);
        const amrex::Real z_out  = z0 + t0 * sz + t_next * sz;
        if (sz >= 0.0) {
            if (solid && z + eps < top) return RAY_BUILDING;
        } else {
            if (solid && z + eps < top) return RAY_BUILDING;     // entered below the top
            if (solid && z_out < top - eps) return RAY_BUILDING; // descends into it
            if (z_out <= z_ground + eps) return RAY_GROUND;
        }
        if (tmx < tmy) { t = tmx; tmx += tdx; ci += stepx; }
        else           { t = tmy; tmy += tdy; cj += stepy; }
        z = z_out;
        if (sz > 0.0 && z >= z_max) return RAY_SKY;
        if (t > max_path) return (sz >= 0.0) ? RAY_SKY : RAY_GROUND;
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
    const amrex::Real phi   = 2.0 * M_PI * (ia + 0.5) / n_az;
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
