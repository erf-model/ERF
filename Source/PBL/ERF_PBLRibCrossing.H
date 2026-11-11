#ifndef ERF_PBL_RIB_CROSSING_H_
#define ERF_PBL_RIB_CROSSING_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Algorithm.H>

/**
 * Result of a bulk Richardson number scan up one column.
 */
struct RibCrossing {
    int          kpbl;     //!< level where Rib first reaches Ribcr; khi when it never does
    amrex::Real  pblh;     //!< PBL height, kept in [pblh_min, pblh_max]
    bool         crossed;  //!< Rib reached Ribcr at or below khi
};

/**
 * Find the PBL height where the bulk Richardson number first reaches Ribcr,
 * scanning one column from its surface cell ksrf up to the top cell khi
 * (YSU, Hong et al. 2006; WRF bl_ysu.F90, the "brint" interpolation).
 *
 * - The scan looks for the first level k in (ksrf, khi] with Rib >= Ribcr.
 *   The surface value Rib(ksrf) only seeds the lower end of the first
 *   interval, as WRF's br does: a surface cell above Ribcr does not stop the
 *   scan by itself.
 * - If Rib(ksrf) and Rib(ksrf+1) are both at or above Ribcr, the height is
 *   z(ksrf) and kpbl = ksrf (WRF: brint = 0); with the bounds below this is
 *   pblh_min.
 * - Otherwise the height is interpolated linearly in Rib between k-1, the last
 *   level below Ribcr, and k. The two levels bracket Ribcr, so the Rib
 *   difference is positive.
 * - No crossing up to khi: the height is the top of the scan, z(khi), and
 *   kpbl = khi (WRF: brint = 1). The pblh_max cap then applies.
 *
 * Every height is then kept in [pblh_min, pblh_max] in that order, so
 * pblh_min wins when a shallow column has pblh_max < pblh_min.
 *
 * @param[in] ksrf      first fluid cell of the column
 * @param[in] khi       top cell of the scan
 * @param[in] Ribcr     critical bulk Richardson number
 * @param[in] rib_at    callable, rib_at(k) gives the bulk Richardson number at level k
 * @param[in] z_at      callable, z_at(k) gives the height of level k
 * @param[in] pblh_min  lower bound on the height
 * @param[in] pblh_max  upper bound on the height
 */
template <typename RibAt, typename ZAt>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
RibCrossing
find_rib_crossing (int ksrf, int khi, amrex::Real Ribcr,
                   RibAt const& rib_at, ZAt const& z_at,
                   amrex::Real pblh_min, amrex::Real pblh_max) noexcept
{
    amrex::Real rib_lo = rib_at(ksrf);
    amrex::Real z_lo   = z_at(ksrf);
    for (int k = ksrf + 1; k <= khi; ++k) {
        const amrex::Real rib_k = rib_at(k);
        const amrex::Real z_k   = z_at(k);
        if (rib_k >= Ribcr) {
            if (rib_lo >= Ribcr) {
                // only possible for k = ksrf + 1: both surface levels at or above Ribcr
                return RibCrossing{ksrf, amrex::max(amrex::min(z_lo, pblh_max), pblh_min), true};
            }
            const amrex::Real pblh = z_lo + (z_k - z_lo) / (rib_k - rib_lo) * (Ribcr - rib_lo);
            return RibCrossing{k, amrex::max(amrex::min(pblh, pblh_max), pblh_min), true};
        }
        rib_lo = rib_k;
        z_lo   = z_k;
    }
    return RibCrossing{khi, amrex::max(amrex::min(z_lo, pblh_max), pblh_min), false};
}

#endif
