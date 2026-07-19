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
 */
struct FireDustCoupling {
    bool        enabled                  = false;
    amrex::Real post_fire_crust_reduction = 0.8;   // fraction of crust removed in burned cells
    const amrex::MultiFab* fire_phi_mf   = nullptr; // phi < 0 = burned

    /**
     * @brief Reduce crust_index in burned cells, lowering u*_t and raising dust emission.
     *
     * Implementation: GPU kernel over dust grid cells; where fire_phi < 0,
     * apply crust *= (1 - reduction); crust = max(crust, 0).
     *
     * @param[in,out] dust_crust_index  Dust layer crust index MultiFab
     * @param[in]     geom_dust         Dust grid geometry
     */
    void apply_burned_area_to_crust(amrex::MultiFab& dust_crust_index,
                                    const amrex::Geometry& geom_dust) const;
};

#endif
#endif
