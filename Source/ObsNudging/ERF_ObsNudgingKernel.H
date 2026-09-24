/**
 * \file ERF_ObsNudgingKernel.H
 */
#ifndef ERF_OBSNUDGINGKERNEL_H_
#define ERF_OBSNUDGINGKERNEL_H_

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

#include <cmath>

namespace obs_nudging {

//
// Everything the nudging tendency reads, as plain device pointers and numbers
// so that it can be captured by value in a device lambda.
//
// The measured profiles of every station are stored back to back: component c
// of station s has count[c*nstations + s] heights, starting at offset
// [c*nstations + s] in z, mean and sigma.  A count of zero means the station
// does not measure c at this moment (or is inactive), and it is skipped.
//
struct ObsNudgingView
{
    int nstations = 0;
    const amrex::Real* x = nullptr;         // station position, domain coordinates
    const amrex::Real* y = nullptr;
    const int*         agl = nullptr;       // 1: heights above the local terrain
    const int*         offset = nullptr;
    const int*         count = nullptr;
    const amrex::Real* z = nullptr;
    const amrex::Real* mean = nullptr;
    const amrex::Real* sigma = nullptr;

    amrex::Real inv_rh2 = amrex::Real(0.0);      // 1 / R_h^2
    amrex::Real inv_rz2 = amrex::Real(0.0);      // 1 / R_z^2
    amrex::Real qmax = amrex::Real(0.0);         // stations with q above this are skipped
    amrex::Real sigma_factor = amrex::Real(0.0); // alpha
    amrex::Real inv_tau = amrex::Real(0.0);      // 1 / tau
    amrex::Real max_rate = amrex::Real(0.0);     // 1 / dt

    // Domain length in a periodic horizontal direction, 0 in one that is not
    amrex::Real period_x = amrex::Real(0.0);
    amrex::Real period_y = amrex::Real(0.0);
};

//
// The separation d of two points along a direction of period L: the nearest
// periodic image when L > 0, so that the two copies of a face on a periodic
// boundary are equally far from every station.
//
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
periodic_separation (amrex::Real d, amrex::Real L) noexcept
{
    if (L > amrex::Real(0.0)) {
        d -= L * amrex::Math::round(d / L);
    }
    return d;
}

//
// The nudging tendency of one component, d(phi)/dt, at a point.
//
// Station s has the weight
//
//     w_s = exp(-q_s / 4),   q_s = r_s^2 / R_h^2 + d_s^2 / R_z^2,
//
// where r_s is the horizontal distance to the station (to its nearest image in
// a periodic direction) and d_s the vertical distance outside the measured
// height range (zero inside it).  Its target is the measured mean, interpolated
// linearly in height (held at the end values outside the range) and widened to
// the band mean +- alpha sigma: the local value is left alone inside the band
// and pulled to its nearer edge outside it,
//
//     t_s = min(max(phi, mean - alpha sigma), mean + alpha sigma).
//
// The stations are combined into one relaxation,
//
//     S = -min(W / tau, 1 / dt) (phi - tbar),
//     W = min(1, sum_s w_s),   tbar = sum_s w_s t_s / sum_s w_s,
//
// so that a point never relaxes faster than 1/tau however many stations are
// near it, and never overshoots its target within a step.
//
// x, y and z are the coordinates of the point, z_ground the terrain height
// under it.
//
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
obs_nudging_tendency (const ObsNudgingView& v, int comp,
                      amrex::Real x, amrex::Real y, amrex::Real z, amrex::Real z_ground,
                      amrex::Real phi) noexcept
{
    amrex::Real sum_w  = amrex::Real(0.0);
    amrex::Real sum_wt = amrex::Real(0.0);

    for (int s = 0; s < v.nstations; ++s)
    {
        const int idx = comp*v.nstations + s;
        const int n = v.count[idx];
        if (n == 0) { continue; }
        const int o = v.offset[idx];

        const amrex::Real h   = (v.agl[s] != 0) ? (z - z_ground) : z;
        const amrex::Real zlo = v.z[o];
        const amrex::Real zhi = v.z[o+n-1];
        const amrex::Real d   = (h < zlo) ? (zlo - h) : ((h > zhi) ? (h - zhi) : amrex::Real(0.0));
        const amrex::Real dxs = periodic_separation(x - v.x[s], v.period_x);
        const amrex::Real dys = periodic_separation(y - v.y[s], v.period_y);
        const amrex::Real q   = (dxs*dxs + dys*dys) * v.inv_rh2 + d*d * v.inv_rz2;
        if (q > v.qmax) { continue; }
        const amrex::Real w = std::exp(amrex::Real(-0.25) * q);

        amrex::Real mean, sig;
        if (h <= zlo) {
            mean = v.mean[o];
            sig  = v.sigma[o];
        } else if (h >= zhi) {
            mean = v.mean[o+n-1];
            sig  = v.sigma[o+n-1];
        } else {
            int k = o;
            while (v.z[k+1] < h) { ++k; }
            const amrex::Real f = (h - v.z[k]) / (v.z[k+1] - v.z[k]);
            mean = (amrex::Real(1.0) - f) * v.mean[k]  + f * v.mean[k+1];
            sig  = (amrex::Real(1.0) - f) * v.sigma[k] + f * v.sigma[k+1];
        }

        const amrex::Real band   = v.sigma_factor * sig;
        const amrex::Real target = amrex::min(amrex::max(phi, mean - band), mean + band);

        sum_w  += w;
        sum_wt += w * target;
    }

    if (sum_w <= amrex::Real(0.0)) { return amrex::Real(0.0); }

    const amrex::Real tbar = sum_wt / sum_w;
    const amrex::Real rate = amrex::min(amrex::min(sum_w, amrex::Real(1.0)) * v.inv_tau, v.max_rate);
    return -rate * (phi - tbar);
}

} // namespace obs_nudging

#endif
