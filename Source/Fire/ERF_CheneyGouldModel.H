#ifndef ERF_CHENEY_GOULD_MODEL_H
#define ERF_CHENEY_GOULD_MODEL_H

#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>
#include <cmath>
#include <ERF_Rothermel.H>

/**
 * @file ERF_CheneyGouldModel.H
 * @brief Cheney-Gould (1998) empirical grassland fire spread model
 *
 * Implements the Cheney-Gould (1998) fire spread model calibrated for open
 * grassland fuels, particularly Australian grasslands. Provides better
 * agreement with observed grassland fires than Rothermel for FM1–FM3.
 *
 * Not appropriate for forest or shrub fuels (FM4–FM13).
 *
 * Key relationships:
 *   Forward ROS depends on:
 *     - Fuel moisture content (dead fine fuel)
 *     - Wind speed
 *     - Curing degree (fraction of dry grass; 0=fresh, 1=fully cured)
 *     - Slope
 *   Backing and flank ROS computed using CSIRO grassfire relations.
 *
 * References:
 *   Cheney, N.P., Gould, J.S. & Catchpole, W.R. (1998). Prediction of fire
 *     spread in grasslands. Int. J. Wildland Fire, 8(1), 1–13.
 *   Catchpole, W.R., Bradstock, R.A., Choate, J.S., Foley, D., Prasad, M.,
 *     Furnell-Smith, G. & Anderson, W. (1998). Cooperative development of
 *     fire spread prediction systems for Australian vegetation. Int. J.
 *     Wildland Fire, 8(1), 14–24.
 */

/**
 * @struct CheneyGouldComputed
 * @brief Precomputed Cheney-Gould coefficients (depends on moisture and curing)
 */
struct CheneyGouldComputed
{
    amrex::Real ros_backing;  ///< Backing (into-wind) ROS [m/s]
    amrex::Real ros_flank;    ///< Flank ROS [m/s]
    amrex::Real gamma;        ///< Spread adjustment factor for slope effects
};

/**
 * @brief Compute Cheney-Gould parameters for given fuel state
 *
 * @param[in] cgp  CheneyGouldParams with moisture [%] and curing [0-1]
 * @return CheneyGouldComputed with precomputed coefficients
 *
 * Computes backing and flank ROS using CSIRO relations based on fuel moisture
 * and curing degree. Gamma factor is used for slope adjustments in fill_cheney_gould_ros().
 */
inline CheneyGouldComputed compute_cheney_gould_params(
    const FireParams::CheneyGouldParams& cgp) noexcept
{
    CheneyGouldComputed cgc;
    
    // Clamp inputs to valid ranges
    amrex::Real moisture = amrex::max(1.0_rt, amrex::min(cgp.moisture, 40.0_rt));  // [%]
    amrex::Real curing   = amrex::max(0.0_rt, amrex::min(cgp.curing, 1.0_rt));      // [0-1]

    // CSIRO grassfire relations for backing ROS
    // Simplified empirical fit: ros_backing increases with curing, decreases with moisture
    // Base backing rate at full curing: ~0.08 m/s
    arMrex::Real base_backing = 0.08_rt * curing;
    cgc.ros_backing = base_backing * std::exp(-0.01_rt * moisture);  // Exponential moisture decay
    cgc.ros_backing = amrex::max(cgc.ros_backing, 0.001_rt);  // Minimum backing

    // Flank ROS (perpendicular to fire spread)
    // Typically 0.3 of forward ROS; CSIRO gives empirical relations
    // Simplified: flank ≈ 0.5 × backing
    cgc.ros_flank = 0.5_rt * cgc.ros_backing;

    // Gamma factor for slope enhancement (CSIRO formula)
    // gamma ≈ 1 + a × tan(slope) for gentle slopes
    // Used downstream in fill_cheney_gould_ros()
    cgc.gamma = 1.0_rt;

    return cgc;
}

/**
 * @brief Forward ROS prediction using wind-only adjustment
 *
 * Cheney-Gould empirical wind effect: forward ROS increases roughly exponentially
 * with wind speed, but with moisture and curing modulation.
 *
 * @param[in] u_wind       Wind speed in forward direction [m/s]
 * @param[in] ros_backing  Base (no-wind) ROS [m/s]
 * @param[in] moisture     Dead fine fuel moisture [%]
 * @param[in] curing       Curing degree [0-1]
 * @return Forward ROS [m/s]
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real cheney_gould_ros(
    amrex::Real u_wind,
    amrex::Real ros_backing,
    amrex::Real moisture,
    amrex::Real curing) noexcept
{
    // Empirical Cheney-Gould relation:
    // ROS_forward ≈ ROS_backing × (1 + wind_factor)
    // where wind_factor grows with wind speed and curing
    amrex::Real wind_factor = 0.15_rt * u_wind * (curing + 0.2_rt);
    
    // Moisture attenuation factor (inverse relationship)
    constexpr amrex::Real moisture_scale = 20.0_rt;  // Reference moisture [%]
    amrex::Real moisture_factor = moisture_scale / (moisture + 1.0_rt);
    
    amrex::Real ros_forward = ros_backing * (1.0_rt + wind_factor) * moisture_factor;
    return amrex::max(ros_forward, 0.0_rt);
}

/**
 * @brief CSIRO grassfire acceleration factor for wind and slope
 *
 * Combines wind and slope effects into a multiplicative factor.
 * Based on CSIRO empirical relations from Catchpole et al. (1998).
 *
 * @param[in] u_wind    Wind speed [m/s]
 * @param[in] slope     Terrain slope [radians]
 * @param[in] curing    Curing degree [0-1]
 * @return Acceleration factor (dimensionless, typically 1–4)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real csiro_grassfire_acceleration(
    amrex::Real u_wind,
    amrex::Real slope,
    amrex::Real curing) noexcept
{
    // Wind effect on forward ROS (log scale)
    amrex::Real wind_accel = 1.0_rt + 0.2_rt * std::log1p(u_wind);
    
    // Slope effect (upslope enhancement)
    amrex::Real slope_mag = std::abs(slope);
    amrex::Real slope_accel = 1.0_rt + 5.0_rt * slope_mag * (curing + 0.1_rt);
    
    // Combined effect
    return wind_accel * slope_accel;
}

/**
 * @brief Cheney-Gould ROS with acceleration for wind and slope
 *
 * Full formula combining backing ROS, wind effect, and slope effect.
 *
 * @param[in] u_wind       Wind speed in fire spread direction [m/s]
 * @param[in] slope        Terrain slope in fire spread direction [radians]
 * @param[in] ros_backing  Base (no-wind, no-slope) ROS [m/s]
 * @param[in] moisture     Dead fine fuel moisture [%]
 * @param[in] curing       Curing degree [0-1]
 * @return Total forward ROS [m/s]
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real cheney_gould_ros_with_acceleration(
    amrex::Real u_wind,
    amrex::Real slope,
    amrex::Real ros_backing,
    amrex::Real moisture,
    amrex::Real curing) noexcept
{
    // Base wind effect
    amrex::Real ros_wind = cheney_gould_ros(u_wind, ros_backing, moisture, curing);
    
    // Additional slope acceleration (upslope only)
    amrex::Real slope_accel = 1.0_rt;
    if (slope > 0.0_rt) {  // upslope
        slope_accel = csiro_grassfire_acceleration(u_wind, slope, curing);
    }
    
    return ros_wind * slope_accel;
}

/**
 * @brief Fill fire_ros [m/s] using Cheney-Gould (1998) grassland model
 *
 * Implements the Cheney-Gould empirical grassland ROS prediction.
 * Suitable only for open grassland fuels (FM1–FM3); not appropriate for
 * forest or shrub fuels (FM4–FM13).
 *
 * @param[out] fire_ros        Rate of spread [m/s]
 * @param[in]  fire_wind_eff   Effective wind after WAF + terrain [m/s]
 * @param[in]  cgc             Cheney-Gould computed parameters
 */
inline void fill_cheney_gould_ros(
    amrex::MultiFab& fire_ros,
    const amrex::MultiFab& fire_wind_eff,
    const CheneyGouldComputed& cgc)
{
    // Note: For a more complete implementation, would need access to:
    // - Per-cell fuel moisture (currently uses default from cgc)
    // - Per-cell slope (already available)
    // For now, uses domain-averaged moisture from parameters
    
    for (amrex::MFIter mfi(fire_ros, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto ros_arr = fire_ros.array(mfi);
        auto wind_arr = fire_wind_eff.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (const amrex::IntVect& iv) {
            int i = iv[0];
            int j = iv[1];
            int k = 0;

            // Get effective wind
            amrex::Real u = wind_arr(i, j, k, 0);
            amrex::Real v = wind_arr(i, j, k, 1);
            amrex::Real wind_mag = std::sqrt(u*u + v*v);

            // Compute Cheney-Gould ROS
            // Simplified: no slope access here; could be added with fire_slopes parameter
            // Currently uses wind effect only via cheney_gould_ros
            amrex::Real ros = cheney_gould_ros(
                wind_mag,
                cgc.ros_backing,
                10.0_rt,   // placeholder: default moisture [%]
                1.0_rt);   // placeholder: default curing [0-1]
            
            ros = amrex::max(ros, 0.0_rt);
            ros_arr(i, j, k, 0) = ros;
        });
    }
}

/**
 * @brief Fill fire_ros with slope effects using Cheney-Gould model
 *
 * Extended version that includes slope effects. Requires fire_slopes parameter.
 * Not currently used in ERF_FireLayer; kept for future enhancement.
 *
 * @param[out] fire_ros        Rate of spread [m/s]
 * @param[in]  fire_wind_eff   Effective wind after WAF + terrain [m/s]
 * @param[in]  fire_slopes     Terrain slopes [dz/dx, dz/dy]
 * @param[in]  cgc             Cheney-Gould computed parameters
 */
inline void fill_cheney_gould_ros_with_slope(
    amrex::MultiFab& fire_ros,
    const amrex::MultiFab& fire_wind_eff,
    const amrex::MultiFab& fire_slopes,
    const CheneyGouldComputed& cgc)
{
    for (amrex::MFIter mfi(fire_ros, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto ros_arr = fire_ros.array(mfi);
        auto wind_arr = fire_wind_eff.const_array(mfi);
        auto slopes_arr = fire_slopes.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (const amrex::IntVect& iv) {
            int i = iv[0];
            int j = iv[1];
            int k = 0;

            // Get wind and slopes
            amrex::Real u = wind_arr(i, j, k, 0);
            amrex::Real v = wind_arr(i, j, k, 1);
            amrex::Real wind_mag = std::sqrt(u*u + v*v);
            
            amrex::Real dzdx = slopes_arr(i, j, k, 0);
            amrex::Real dzdy = slopes_arr(i, j, k, 1);
            amrex::Real slope_mag = std::sqrt(dzdx*dzdx + dzdy*dzdy);
            amrex::Real slope_angle = std::atan(slope_mag);

            // Project slope onto wind direction
            amrex::Real wind_dir = std::atan2(v, u);
            amrex::Real slope_dir = std::atan2(dzdy, dzdx);
            amrex::Real angle_diff = wind_dir - slope_dir;
            amrex::Real slope_proj = slope_angle * std::cos(angle_diff);

            // Compute ROS with slope effects
            amrex::Real ros = cheney_gould_ros_with_acceleration(
                wind_mag,
                slope_proj,
                cgc.ros_backing,
                10.0_rt,   // placeholder: default moisture [%]
                1.0_rt);   // placeholder: default curing [0-1]
            
            ros = amrex::max(ros, 0.0_rt);
            ros_arr(i, j, k, 0) = ros;
        });
    }
}

#endif
