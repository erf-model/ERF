#ifndef ERF_MACARTHUR_MODEL_H
#define ERF_MACARTHUR_MODEL_H

#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>
#include <cmath>

/**
 * @file ERF_MacArthurModel.H
 * @brief MacArthur (1966) Australian fire spread model
 *
 * Implements the MacArthur (1966) Mark 5 Forest Fire Danger Meter ROS formula.
 * Used in WRF-SFIRE as ibeh=0 (Clark et al. 2004).
 *
 * Formula:
 *   R [m/s] = 0.18 * (exp(0.8424 * max(U, 0)) - 1) + 0.18
 *           = 0.18 * exp(0.8424 * max(U, 0))    [with wind]
 *   R_backing = 0.18 m/s                         [no-wind backing rate]
 *
 * This is a simple empirical wind-only formula suitable for Australian
 * grassland and open forest fuels. Not appropriate for complex terrain or
 * high-fuel-complexity environments.
 *
 * References:
 *   Clark, T.L., Coen, J.L. & Latham, D. (2004). Description of a coupled
 *     atmosphere-fire model. Int. J. Wildland Fire, 13, 49–63.
 *   McArthur, A.G. (1966). Weather and Grassland Fire Behaviour. Leaflet 100,
 *     Forestry and Timber Bureau, Canberra.
 */

/**
 * @brief MacArthur (1966) ROS formula for given wind speed
 *
 * Computes rate of spread using the simple exponential wind formula.
 *
 * @param[in] wind_speed_ms  Wind speed in the fire spread direction [m/s]
 * @return ROS [m/s]
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real macarthur_ros(amrex::Real wind_speed_ms) noexcept
{
    constexpr amrex::Real backing = 0.18;
    const amrex::Real wind_contribution = backing * (std::exp(0.8424 * amrex::max(wind_speed_ms, 0.0)) - 1.0);
    return backing + wind_contribution;
}

/**
 * @brief Fill fire_ros [m/s] using MacArthur (1966) formula
 *
 * Simple wind-driven ROS model for grassland and open forest.
 *
 * @param[out] fire_ros        Rate of spread [m/s]
 * @param[in]  fire_wind_eff   Effective wind after WAF + terrain [m/s]
 */
inline void fill_macarthur_ros(
    amrex::MultiFab& fire_ros,
    const amrex::MultiFab& fire_wind_eff)
{
    for (amrex::MFIter mfi(fire_ros, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto ros_arr = fire_ros.array(mfi);
        auto wind_arr = fire_wind_eff.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (const amrex::IntVect& iv) {
            int i = iv[0];
            int j = iv[1];
            int k = 0;

            // Get effective wind magnitude
            amrex::Real u = wind_arr(i, j, k, 0);
            amrex::Real v = wind_arr(i, j, k, 1);
            amrex::Real wind_mag = std::sqrt(u*u + v*v);

            // Compute MacArthur ROS
            amrex::Real ros = macarthur_ros(wind_mag);
            
            ros_arr(i, j, k, 0) = ros;
        });
    }
}

#endif
