#ifndef ERF_METGRID_LAND_UTILS_H_
#define ERF_METGRID_LAND_UTILS_H_

/**
 * \file ERF_MetgridLandUtils.H
 *
 * Interpolation in depth and in time used to build SLM's initial state from a met_em file
 * (ERF::init_lsm_from_metgrid).
 */
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <string>

namespace metgrid_land {

// The Noah soil layers real.exe writes to wrfinput (DZS) for a 4-layer LSM, which SLM's
// initial state from met_em files is interpolated to, and so the layers SLM takes then
constexpr int nsoil_noah = 4;
constexpr amrex::Real dzs_noah[nsoil_noah] = {amrex::Real(0.1), amrex::Real(0.3),
                                              amrex::Real(0.6), amrex::Real(1.0)};

/**
 * Linear interpolation in depth from the source layer midpoints to z, held constant
 * beyond the shallowest and deepest midpoints.
 *
 * @param z_src  depths of the source layer midpoints, increasing
 * @param v_src  values at those depths
 * @param z      depth to interpolate to
 */
inline amrex::Real
interp_in_depth (const amrex::Vector<amrex::Real>& z_src,
                 const amrex::Vector<amrex::Real>& v_src,
                 amrex::Real z)
{
    const int n = z_src.size();
    if (z <= z_src[0])   { return v_src[0]; }
    if (z >= z_src[n-1]) { return v_src[n-1]; }
    int m = 0;
    while (z > z_src[m+1]) { ++m; }
    const amrex::Real w = (z - z_src[m]) / (z_src[m+1] - z_src[m]);
    return (amrex::Real(1.0) - w) * v_src[m] + w * v_src[m+1];
}

/**
 * The two months whose mid-month values bracket the given date, and the weight of the
 * second, as in real.exe's monthly_interp_to_date: each monthly value is taken to apply on
 * the 15th of its month, and values in between are interpolated linearly in time.
 *
 * @param date "YYYY-MM-DD_hh:mm:ss"
 * @param m0   the earlier month (0 = January)
 * @param m1   the later month
 * @param w1   the weight of the later month
 */
inline void
monthly_interp_weights (const std::string& date, int& m0, int& m1, amrex::Real& w1)
{
    const int year  = std::stoi(date.substr(0,4));
    const int month = std::stoi(date.substr(5,2));
    const int day   = std::stoi(date.substr(8,2));
    const int hour  = std::stoi(date.substr(11,2));

    const bool leap = (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
    const int days_in_month[12] = {31, leap ? 29 : 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    // Day of the year (starting from 0) of the date and of the middle of each month
    amrex::Real mid[12];
    int start = 0;
    for (int m = 0; m < 12; ++m) {
        mid[m] = amrex::Real(start + 14);
        start += days_in_month[m];
    }
    const int days_in_year = start;
    int doy = day - 1;
    for (int m = 0; m < month-1; ++m) { doy += days_in_month[m]; }
    const amrex::Real t = amrex::Real(doy) + amrex::Real(hour) / amrex::Real(24.0);

    if (t < mid[0] || t >= mid[11]) {
        // Between mid-December and mid-January, wrapping around the end of the year
        m0 = 11;
        m1 = 0;
        const amrex::Real t0 = mid[11];
        const amrex::Real t1 = mid[0] + amrex::Real(days_in_year);
        const amrex::Real tt = (t < mid[0]) ? t + amrex::Real(days_in_year) : t;
        w1 = (tt - t0) / (t1 - t0);
    } else {
        m0 = 0;
        while (t >= mid[m0+1]) { ++m0; }
        m1 = m0 + 1;
        w1 = (t - mid[m0]) / (mid[m1] - mid[m0]);
    }
}

} // namespace metgrid_land

#endif
