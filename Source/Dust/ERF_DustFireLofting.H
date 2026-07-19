#ifndef ERF_DUST_FIRE_LOFTING_H
#define ERF_DUST_FIRE_LOFTING_H

#ifdef ERF_USE_DUST

#include <AMReX_MultiFab.H>

/**
 * @brief Enhance dust emission flux where fire heat flux drives convective lofting.
 *
 * Applies: flux *= (1 + k_loft * max(Q_fire - Q_threshold, 0) / Q_ref)
 * capped at (1 + k_loft). No effect when Q_fire <= Q_threshold.
 *
 * MPI safety: fire_heat_scratch must already be ParallelCopied onto the
 * dust BoxArray/DM by the caller before this function is called.
 *
 * @param[in,out] dust_emission_flux  Emission flux on dust grid [kg/m2/s]
 * @param[in]     fire_heat_scratch   Fire heat flux on dust BA [W/m2]
 * @param[in]     n_bins              Number of active size bins
 * @param[in]     k_loft              Enhancement factor (dimensionless)
 * @param[in]     Q_threshold         Minimum heat flux to trigger lofting [W/m2]
 * @param[in]     Q_ref               Reference heat flux for normalisation [W/m2]
 * @param[in]     dust_debug          Print debug output if true
 * @param[in]     step                Current step number for debug output
 */
void apply_fire_lofting_to_emission_flux(
    amrex::MultiFab&       dust_emission_flux,
    const amrex::MultiFab& fire_heat_scratch,
    int                    n_bins,
    amrex::Real            k_loft,
    amrex::Real            Q_threshold,
    amrex::Real            Q_ref,
    bool                   dust_debug,
    int                    step);

#endif // ERF_USE_DUST
#endif // ERF_DUST_FIRE_LOFTING_H
