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
};

#endif
#endif
