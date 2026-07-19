#ifndef ERF_FIRE_DUST_COUPLING_H
#define ERF_FIRE_DUST_COUPLING_H

#if defined(ERF_ENABLE_FIRE) && defined(ERF_USE_DUST)

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

/**
 * @brief Fire-Dust coupling interactions.
 *
 * Handles interactions between fire and dust modules, particularly
 * burned area effects on dust emission.
 *
 * The fire level-set phi lives on the fire sub-grid (grid_ratio times
 * finer than the atmosphere). The dust crust index lives on the dust
 * sub-grid (which may have a different grid_ratio). We map burned
 * status from fire to dust using physical (x,y) coordinates so that
 * the two grids do not need to be aligned.
 */
struct FireDustCoupling {
    bool        enabled                   = false;
    amrex::Real post_fire_crust_reduction = 0.8;  // fraction of crust removed in burned cells

    // Fire level-set phi on the fire sub-grid (phi < 0 => burned).
    const amrex::MultiFab* fire_phi_mf  = nullptr;
    // Geometry of the fire sub-grid (needed for physical-coordinate mapping).
    amrex::Geometry        geom_fire;

    /**
     * @brief Reduce crust_index in burned cells, lowering u*_t and raising dust emission.
     *
     * For each dust cell we compute its physical centre (x,y), find the
     * containing fire cell, and if fire_phi < 0 there we apply
     *   crust *= (1 - post_fire_crust_reduction),  crust = max(crust, 0).
     *
     * @param[in,out] dust_crust_index  Dust layer crust index MultiFab (dust grid)
     * @param[in]     geom_dust         Geometry of the dust sub-grid
     */
    void apply_burned_area_to_crust(amrex::MultiFab&       dust_crust_index,
                                    const amrex::Geometry& geom_dust) const;

    // ── Phase 2: fire outflow wind → dust u* ──────────────────────────────────
    // When fire_wind_to_dust = true, fire_wind_eff (computed by FireLayer each
    // timestep) is averaged over the fire-to-dust grid ratio and converted to
    // a friction velocity via log-law. The maximum of background MRF u* and the
    // fire-derived u* is written back to dust_ustar_in. Cells immediately
    // downwind of the fire perimeter see elevated u* above the MRF baseline,
    // triggering or enhancing dust emission that would not occur from background
    // wind alone.
    bool        fire_wind_to_dust  = true;
    amrex::Real fire_wind_z0       = 0.1;  // surface roughness for log-law [m]
    amrex::Real fire_wind_zref     = 6.1;  // reference height matching fire wind_ref_ht [m]

    // Pointer to fire effective wind field (u, v) on fire sub-grid.
    // Set from FireLayer::get_wind_eff() each timestep before calling
    // apply_fire_wind_to_dust_ustar.
    const amrex::MultiFab* fire_wind_eff_mf = nullptr;
    // Geometry of the fire sub-grid (needed for physical-coordinate mapping,
    // same pattern as geom_fire used by apply_burned_area_to_crust).
    amrex::Geometry geom_fire_wind;

    /**
     * @brief Interaction 2: fire outflow wind raises effective u* on dust grid.
     *
     * For each dust cell, averages fire_wind_eff over the C×C fire cells
     * covering that dust cell, derives u* via log-law
     * (u* = U * kappa / ln(zref/z0)), and takes the maximum with the
     * background MRF u* already in dust_ustar_in.
     *
     * MPI safety: fire_wind_eff_mf must already be copied into a scratch
     * MultiFab on the dust BoxArray/DM by the caller (ERF.cpp) before this
     * function is called — same pattern as LNG_MPI_SKILLS Rule B3. This
     * function receives the pre-copied scratch MultiFab as fire_wind_scratch.
     *
     * @param[in,out] dust_ustar_in      u* on dust grid [m/s] — modified in place
     * @param[in]     fire_wind_scratch  fire_wind_eff already ParallelCopied onto
     *                                   the dust BoxArray/DM (2-comp: u, v) [m/s]
     * @param[in]     geom_dust          Geometry of the dust sub-grid
     * @param[in]     z0                 Surface roughness [m]
     * @param[in]     zref               Reference height [m]
     * @param[in]     C                  Fire-to-dust grid ratio (fire cells per dust cell)
     */
    void apply_fire_wind_to_dust_ustar(
       amrex::MultiFab&       dust_ustar_in,
       const amrex::MultiFab& fire_wind_scratch,
       const amrex::Geometry& geom_dust,
       amrex::Real            z0,
       amrex::Real            zref,
       int                    C) const;
};

#endif
#endif
