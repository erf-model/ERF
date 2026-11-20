#ifndef ERF_INFLOW_PROFILE_H_
#define ERF_INFLOW_PROFILE_H_

/**
 * Terrain-following inflow profiles for Inflow and InflowOutflow faces
 * (xlo.inflow_profile = "log_law" or "file").
 *
 * A profile tabulates u, v and, when given, theta and tke against height above
 * the local ground. Each boundary cell looks the profile up at its own height
 * above the terrain beneath it, so the inflow follows the terrain along a
 * sloping face, as WindNinja's inlets do. The terrain on the boundary is left
 * as it is.
 */

#include <istream>
#include <string>

#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

struct InflowProfile
{
    bool active    = false;
    bool has_theta = false;
    bool has_tke   = false;
    std::string source;   ///< "log_law" or the profile file name

    /// Host tables against height above the ground (strictly increasing)
    amrex::Vector<amrex::Real> z, u, v, theta, tke;

    /// Device copies used by the boundary fills
    amrex::Gpu::DeviceVector<amrex::Real> z_d, u_d, v_d, theta_d, tke_d;

    void copy_to_device ();
};

/** The neutral log law of inflow_log_law.* */
struct InflowLogLaw
{
    amrex::Real speed      = amrex::Real(-1.0);   ///< wind speed at `height` [m/s] (required)
    amrex::Real height     = amrex::Real(10.0);   ///< reference height above the ground [m]
    amrex::Real direction  = amrex::Real(270.0);  ///< direction the wind blows from [deg, meteorological]
    amrex::Real z0         = amrex::Real(0.1);    ///< roughness length [m]
    amrex::Real max_speed  = amrex::Real(-1.0);   ///< cap on the speed; <= 0 means none
    amrex::Real tke_zscale = amrex::Real(700.0);  ///< tke tapers linearly to zero at u* tke_zscale
    amrex::Real Cmu0       = amrex::Real(0.5562); ///< tke at the ground is u*^2 / Cmu0^2
};

/** Parse a profile table; returns an empty string on success, otherwise the error. */
std::string parse_inflow_profile (std::istream& is, InflowProfile& prof);

/** Read a profile file, aborting with the file name on any error. */
void read_inflow_profile_file (const std::string& fname, InflowProfile& prof);

/** Tabulate the log law from the ground to ztop; returns an empty string on success. */
std::string make_log_law_inflow_profile (const InflowLogLaw& p, amrex::Real ztop, InflowProfile& prof);

/** Linear interpolation in a strictly increasing table, held constant beyond both ends. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
inflow_profile_value (const amrex::Real* z, const amrex::Real* f, const int n, const amrex::Real zq) noexcept
{
    if (zq <= z[0])   { return f[0]; }
    if (zq >= z[n-1]) { return f[n-1]; }
    int lo = 0;
    int hi = n - 1;
    while (hi - lo > 1) {
        const int mid = (lo + hi) / 2;
        if (z[mid] <= zq) { lo = mid; } else { hi = mid; }
    }
    const amrex::Real w = (zq - z[lo]) / (z[hi] - z[lo]);
    return f[lo] + w * (f[hi] - f[lo]);
}

/**
 * Mean height above the ground of the nodes [ilo,ihi] x [jlo,jhi] between
 * levels k and k+1 of the terrain-fitted mesh: a cell centre (two nodes in each
 * lateral direction), an x-face (ilo == ihi) or a y-face (jlo == jhi). The
 * ground is the same nodes on level klo. Indices are clamped to the array.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
height_above_ground (const amrex::Array4<const amrex::Real>& z_nd,
                     int ilo, int ihi, int jlo, int jhi, const int k, const int klo) noexcept
{
    const auto lo = amrex::lbound(z_nd);
    const auto hi = amrex::ubound(z_nd);
    ilo = amrex::max(lo.x, amrex::min(ilo, hi.x));
    ihi = amrex::max(lo.x, amrex::min(ihi, hi.x));
    jlo = amrex::max(lo.y, amrex::min(jlo, hi.y));
    jhi = amrex::max(lo.y, amrex::min(jhi, hi.y));
    const int kc = amrex::max(lo.z, amrex::min(k,   hi.z - 1));
    const int kg = amrex::max(lo.z, amrex::min(klo, hi.z));
    amrex::Real sum = amrex::Real(0.0);
    int count = 0;
    for (int j = jlo; j <= jhi; ++j) {
        for (int i = ilo; i <= ihi; ++i) {
            sum += amrex::Real(0.5) * (z_nd(i,j,kc) + z_nd(i,j,kc+1)) - z_nd(i,j,kg);
            ++count;
        }
    }
    return sum / static_cast<amrex::Real>(count);
}

#endif
