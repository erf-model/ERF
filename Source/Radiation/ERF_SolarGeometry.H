#ifndef ERF_SOLAR_GEOMETRY_H_
#define ERF_SOLAR_GEOMETRY_H_

#include <AMReX_GpuControl.H>
#include <AMReX_REAL.H>
#include <cmath>

/**
 * @file ERF_SolarGeometry.H
 * @brief Solar geometry calculations for time-varying solar position (Phase 16).
 *
 * Provides GPU-safe helper functions to compute solar position from:
 * - Simulation time (UTC seconds from day start or absolute time)
 * - Site latitude and longitude
 * - Day-of-year (for solar declination)
 * - Time zone offset (for local solar time)
 *
 * All functions are AMREX_GPU_DEVICE AMREX_FORCE_INLINE for GPU safety.
 * Input validation includes guards against NaN, Inf, and out-of-range values
 * following the clamp_finite() conventions established in prior phases.
 *
 * References:
 * -----------
 * - Spencer, J. W., 1971: Fourier series representation of the position of the sun.
 *   Search, 2(5), 172-172.
 * - Duffie, J. A., and W. A. Beckman, 1991: Solar Engineering of Thermal Processes.
 *   John Wiley & Sons.
 * - Iqbal, M., 1983: An Introduction to Solar Radiation.
 *   Academic Press.
 */

/**
 * @brief Clamp latitude to valid range [-90, 90] degrees.
 * @param[in] lat_deg Input latitude in degrees.
 * @return Clamped latitude value.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real clamp_latitude(amrex::Real lat_deg)
{
    if (!std::isfinite(lat_deg)) return 0.0;  // Default to equator if invalid
    if (lat_deg < -90.0) return -90.0;
    if (lat_deg > 90.0) return 90.0;
    return lat_deg;
}

/**
 * @brief Clamp longitude to valid range [-180, 180] degrees.
 * @param[in] lon_deg Input longitude in degrees.
 * @return Clamped longitude value (wrapped to [-180, 180]).
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real clamp_longitude(amrex::Real lon_deg)
{
    if (!std::isfinite(lon_deg)) return 0.0;  // Default to prime meridian if invalid
    
    // Normalize to [-180, 180]
    amrex::Real lon = lon_deg;
    while (lon > 180.0) lon -= 360.0;
    while (lon < -180.0) lon += 360.0;
    return lon;
}

/**
 * @brief Clamp day-of-year to valid range [1, 366].
 * @param[in] doy Input day-of-year (1-based).
 * @return Clamped day-of-year value.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real clamp_day_of_year(amrex::Real doy)
{
    if (!std::isfinite(doy)) return 172.0;  // Default to summer solstice if invalid
    if (doy < 1.0) return 1.0;
    if (doy > 366.0) return 366.0;
    return doy;
}

/**
 * @brief Compute solar declination angle from day-of-year.
 *
 * Uses Spencer's (1971) Fourier-series approximation:
 * declination [rad] = sum of Fourier terms with period 365.25 days.
 *
 * @param[in] day_of_year Day-of-year (1-based, 1-366).
 * @return Solar declination angle [radians].
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_solar_declination(amrex::Real day_of_year)
{
    // Clamp input
    day_of_year = clamp_day_of_year(day_of_year);
    
    // Fractional year in radians (365.25 day average year)
    amrex::Real gamma = 2.0 * M_PI * (day_of_year - 1.0) / 365.25;
    
    // Spencer's Fourier approximation (truncated to 4 terms for speed)
    // Full formula includes eqtime and declination; here we extract declination only
    amrex::Real decl = 0.006918
        - 0.399912 * std::cos(gamma) + 0.070257 * std::sin(gamma)
        - 0.006758 * std::cos(2.0 * gamma) + 0.000907 * std::sin(2.0 * gamma)
        - 0.00248 * std::cos(3.0 * gamma) + 0.00031 * std::sin(3.0 * gamma);
    
    return decl;
}

/**
 * @brief Compute equation of time correction [minutes].
 *
 * Uses Spencer's (1971) Fourier-series approximation to account for
 * the difference between solar time and mean solar time.
 *
 * @param[in] day_of_year Day-of-year (1-based, 1-366).
 * @return Equation of time [minutes].
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_equation_of_time(amrex::Real day_of_year)
{
    // Clamp input
    day_of_year = clamp_day_of_year(day_of_year);
    
    // Fractional year in radians
    amrex::Real gamma = 2.0 * M_PI * (day_of_year - 1.0) / 365.25;
    
    // Spencer's Fourier approximation for equation of time [minutes]
    amrex::Real eqtime = 229.18
        * (0.000075 + 0.001868 * std::cos(gamma) - 0.032077 * std::sin(gamma)
           - 0.014615 * std::cos(2.0 * gamma) - 0.040849 * std::sin(2.0 * gamma));
    
    return eqtime;
}

/**
 * @brief Compute solar hour angle from simulation time and site parameters.
 *
 * The hour angle is the angular displacement of the sun from the local meridian,
 * computed as:
 *   hour_angle [deg] = 15 * (solar_time_hours - 12)
 * where:
 *   solar_time_hours = time_utc_hours + time_zone_offset_hours + equation_of_time_hours
 *
 * @param[in] time_utc_seconds UTC time [seconds] from start of simulation day (00:00 UTC).
 *                             Assumes time wraps modulo 86400 s (24 hours).
 * @param[in] longitude_deg Site longitude [degrees, -180 to 180].
 * @param[in] time_zone_offset_hours Time zone offset from UTC [hours].
 * @param[in] day_of_year Day-of-year (1-based, 1-366) for equation-of-time correction.
 * @return Solar hour angle [radians].
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_solar_hour_angle(
    amrex::Real time_utc_seconds,
    amrex::Real longitude_deg,
    amrex::Real time_zone_offset_hours,
    amrex::Real day_of_year)
{
    // Validate inputs
    if (!std::isfinite(time_utc_seconds)) return 0.0;
    if (!std::isfinite(longitude_deg)) longitude_deg = 0.0;
    if (!std::isfinite(time_zone_offset_hours)) time_zone_offset_hours = 0.0;
    
    // Clamp longitude to [-180, 180]
    longitude_deg = clamp_longitude(longitude_deg);
    
    // Convert UTC time to hours within day [0, 24)
    amrex::Real time_utc_hours = std::fmod(time_utc_seconds / 3600.0, 24.0);
    if (time_utc_hours < 0.0) time_utc_hours += 24.0;
    
    // Add time zone offset to get local time
    amrex::Real time_local_hours = time_utc_hours + time_zone_offset_hours;
    
    // Account for longitude: solar time differs from local time due to longitude offset
    // Standard meridian offset: each 15 degrees of longitude = 1 hour of time
    // Longitude convention: east of Greenwich is positive, west is negative
    // For solar time, we want solar noon when sun crosses the local meridian
    amrex::Real time_solar_hours = time_local_hours + longitude_deg / 15.0;
    
    // Apply equation of time correction [minutes -> hours]
    amrex::Real eqtime_hours = compute_equation_of_time(day_of_year) / 60.0;
    time_solar_hours += eqtime_hours;
    
    // Compute hour angle: 15 degrees per hour (360 deg / 24 hours)
    // Hour angle = 0 at solar noon, negative in morning, positive in afternoon
    amrex::Real hour_angle_deg = 15.0 * (time_solar_hours - 12.0);
    
    // Convert to radians
    amrex::Real hour_angle_rad = hour_angle_deg * M_PI / 180.0;
    
    return hour_angle_rad;
}

/**
 * @brief Compute solar zenith angle from declination, hour angle, and latitude.
 *
 * Standard spherical-trigonometry formula:
 *   cos(zenith) = sin(lat) * sin(decl) + cos(lat) * cos(decl) * cos(hour_angle)
 *
 * @param[in] latitude_deg Site latitude [degrees, -90 to 90, positive north].
 * @param[in] declination_rad Solar declination [radians].
 * @param[in] hour_angle_rad Solar hour angle [radians].
 * @return Solar zenith angle [radians, 0 = sun directly overhead, pi/2 = horizon].
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_solar_zenith_angle(
    amrex::Real latitude_deg,
    amrex::Real declination_rad,
    amrex::Real hour_angle_rad)
{
    // Clamp latitude to [-90, 90]
    latitude_deg = clamp_latitude(latitude_deg);
    
    // Validate trig inputs
    if (!std::isfinite(declination_rad)) declination_rad = 0.0;
    if (!std::isfinite(hour_angle_rad)) hour_angle_rad = 0.0;
    
    // Convert latitude to radians
    amrex::Real latitude_rad = latitude_deg * M_PI / 180.0;
    
    // Compute cosine of zenith angle
    amrex::Real cos_zenith = std::sin(latitude_rad) * std::sin(declination_rad)
        + std::cos(latitude_rad) * std::cos(declination_rad) * std::cos(hour_angle_rad);
    
    // Clamp to physically valid range [-1, 1] (handles numerical precision issues)
    if (cos_zenith < -1.0) cos_zenith = -1.0;
    if (cos_zenith > 1.0) cos_zenith = 1.0;
    
    // Compute zenith angle from cosine
    // Use std::acos for GPU safety (part of std:: math on GPU)
    amrex::Real zenith_rad = std::acos(cos_zenith);
    
    return zenith_rad;
}

/**
 * @brief Compute solar azimuth angle from declination, hour angle, latitude, and zenith.
 *
 * Standard spherical-trigonometry formula (accounting for all quadrants):
 *   sin(azimuth) = sin(hour_angle) * cos(decl) / sin(zenith)
 *   cos(azimuth) = (sin(decl) - sin(zenith_angle) * sin(lat)) / (cos(zenith_angle) * cos(lat))
 *
 * Azimuth convention: 0° = north, 90° = east, 180° = south, 270° = west.
 *
 * @param[in] latitude_deg Site latitude [degrees, -90 to 90].
 * @param[in] declination_rad Solar declination [radians].
 * @param[in] hour_angle_rad Solar hour angle [radians].
 * @param[in] zenith_rad Solar zenith angle [radians] (computed from declination/hour_angle/lat).
 * @return Solar azimuth angle [radians, 0 = north, pi/2 = east, pi = south, 3pi/2 = west].
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_solar_azimuth_angle(
    amrex::Real latitude_deg,
    amrex::Real declination_rad,
    amrex::Real hour_angle_rad,
    amrex::Real zenith_rad)
{
    // Clamp latitude
    latitude_deg = clamp_latitude(latitude_deg);
    
    // Validate trig inputs
    if (!std::isfinite(declination_rad)) declination_rad = 0.0;
    if (!std::isfinite(hour_angle_rad)) hour_angle_rad = 0.0;
    if (!std::isfinite(zenith_rad)) zenith_rad = M_PI / 2.0;  // Default to horizon if invalid
    
    // Convert latitude to radians
    amrex::Real latitude_rad = latitude_deg * M_PI / 180.0;
    
    // Compute sine and cosine of zenith for azimuth formula
    amrex::Real sin_zenith = std::sin(zenith_rad);
    amrex::Real cos_zenith = std::cos(zenith_rad);
    
    // Guard against division by zero when sun is directly overhead (zenith ≈ 0)
    // In this case, azimuth is undefined; we return 0 (north) as convention
    if (sin_zenith < 1.0e-6) return 0.0;
    
    // Compute numerator and denominator for sine/cosine of azimuth
    amrex::Real sin_azimuth = std::sin(hour_angle_rad) * std::cos(declination_rad) / sin_zenith;
    amrex::Real cos_azimuth = (std::sin(declination_rad) - sin_zenith * std::sin(latitude_rad))
        / (cos_zenith * std::cos(latitude_rad));
    
    // Use atan2 to get full-quadrant azimuth angle [-pi, pi]
    // Then adjust to [0, 2*pi] (north = 0, east = pi/2)
    amrex::Real azimuth_rad = std::atan2(sin_azimuth, cos_azimuth);
    
    // Normalize to [0, 2*pi]
    if (azimuth_rad < 0.0) azimuth_rad += 2.0 * M_PI;
    
    return azimuth_rad;
}

/**
 * @brief Compute cosine of solar zenith angle (convenience wrapper for flux calculations).
 *
 * This is the primary interface for radiation flux calculations:
 * Shortwave direct-beam flux = S0 * cos_zenith * transmission_factor
 *
 * @param[in] time_utc_seconds UTC time [seconds] from start of simulation day.
 * @param[in] latitude_deg Site latitude [degrees].
 * @param[in] longitude_deg Site longitude [degrees].
 * @param[in] day_of_year Day-of-year (1-based).
 * @param[in] time_zone_offset_hours Time zone offset from UTC [hours].
 * @return cos(zenith_angle) in range [-1, 1], where values ≤ 0 indicate sun below horizon.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_cos_zenith_angle(
    amrex::Real time_utc_seconds,
    amrex::Real latitude_deg,
    amrex::Real longitude_deg,
    amrex::Real day_of_year,
    amrex::Real time_zone_offset_hours = 0.0)
{
    // Compute solar declination
    amrex::Real declination = compute_solar_declination(day_of_year);
    
    // Compute solar hour angle
    amrex::Real hour_angle = compute_solar_hour_angle(
        time_utc_seconds, longitude_deg, time_zone_offset_hours, day_of_year);
    
    // Compute zenith angle
    amrex::Real zenith = compute_solar_zenith_angle(latitude_deg, declination, hour_angle);
    
    // Return cosine (negative values mean sun below horizon)
    amrex::Real cos_zenith = std::cos(zenith);
    
    return cos_zenith;
}

#endif // ERF_SOLAR_GEOMETRY_H_
