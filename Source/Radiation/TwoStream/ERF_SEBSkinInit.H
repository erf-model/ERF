#ifndef ERF_SEB_SKIN_INIT_H_
#define ERF_SEB_SKIN_INIT_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_INT.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_iMultiFab.H>

/**
 * @file ERF_SEBSkinInit.H
 * @brief The first value of the two-stream force-restore balance's skin temperature from
 * the skin temperature of the initial file (TSK of a wrfinput file, SKINTEMP of a metgrid
 * file), with erf.radiation.seb_init_skin_from_file.
 *
 * Without the input every column starts at the one value erf.rad_t_sfc. With it a land
 * column starts at the file's value. A water column keeps erf.rad_t_sfc: the balance does
 * not run there and the sweep takes the surface layer's sea surface temperature.
 */

/**
 * @brief True when @p temperature [K] can start the skin: finite and inside the range the
 * balance holds the skin in, [@p t_lo, @p t_hi] (erf.radiation.seb_prognostic_t_min_k and
 * seb_prognostic_t_max_k).
 *
 * A value outside the range would be moved to the nearest bound by the first update of the
 * balance, so it is refused here.
 */
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool
seb_skin_init_usable (const amrex::Real temperature, const amrex::Real t_lo,
                      const amrex::Real t_hi) noexcept
{
    return amrex::Math::isfinite(temperature) && temperature >= t_lo && temperature <= t_hi;
}

/** @brief What seb_init_skin_from_file did on the boxes of this rank. */
struct SEBSkinInitReport
{
    amrex::Long n_set = 0;      //!< land columns that took the file's value
    amrex::Long n_refused = 0;  //!< land columns whose file value is not usable (left as they were)
    amrex::Real t_min = 0;      //!< lowest value set [K]; meaningful when n_set > 0
    amrex::Real t_max = 0;      //!< highest value set [K]; meaningful when n_set > 0
};

/**
 * @brief Copy @p file_skin [K] into the land columns of @p t_sfc [K], valid columns only.
 *
 * @param[in,out] t_sfc     the balance's skin temperature of one level
 * @param[in]     file_skin the file's skin temperature on the same grids
 * @param[in]     land_mask ERF::lmask_lev of the level (0 is water); nullptr means all land
 * @param[in]     t_lo      lowest usable value [K]
 * @param[in]     t_hi      highest usable value [K]
 * @return the counts and range of this rank. The caller sums them over the ranks, stops
 *         when a column was refused, and fills the ghost columns of @p t_sfc.
 */
[[nodiscard]] SEBSkinInitReport
seb_init_skin_from_file (amrex::MultiFab& t_sfc, const amrex::MultiFab& file_skin,
                         const amrex::iMultiFab* land_mask, amrex::Real t_lo, amrex::Real t_hi);

#endif
