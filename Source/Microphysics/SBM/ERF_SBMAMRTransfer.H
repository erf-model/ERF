#ifndef ERF_SBM_AMR_TRANSFER_H_
#define ERF_SBM_AMR_TRANSFER_H_

#include "ERF_SBMLayout.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <string>

namespace erf_sbm {

/** A spectrum, carrier, and mapped measure tied to explicit semantic times. */
struct SBMAMRStateView {
    const amrex::MultiFab* spectrum{nullptr};
    double spectrum_time{0.0};
    const amrex::MultiFab* dry_air_density{nullptr};
    int density_component{-1};
    double density_time{0.0};
    const amrex::MultiFab* mapped_measure{nullptr};
    int measure_component{-1};
    double measure_time{0.0};
};

/**
 * Form an H-conservative coarse candidate from matching-time spectral views.
 * On success, overwrite the destination's valid cells with the complete,
 * admitted result. On failure (including a rejection on another MPI rank),
 * leave every destination valid cell and all source views unchanged. This
 * operator does not itself commit authoritative SBM state.
 */
bool RestrictMappedSpectrum (const SBMLayout& layout,
                             const SBMAMRStateView& fine,
                             const SBMAMRStateView& coarse,
                             const amrex::IntVect& ratio,
                             int coarse_level,
                             amrex::MultiFab& coarse_candidate,
                             std::string& diagnostic);

/**
 * Form a same-time piecewise-constant carrier-relative fine candidate. On
 * success, overwrite the destination's valid cells with the complete, admitted
 * result. On failure (including a rejection on another MPI rank), leave every
 * destination valid cell and all source views unchanged. This operator does
 * not itself commit authoritative SBM state.
 */
bool ProlongCarrierRelativeSpectrum (const SBMLayout& layout,
                                     const SBMAMRStateView& coarse,
                                     const SBMAMRStateView& fine_target,
                                     const amrex::Geometry& coarse_geometry,
                                     const amrex::Geometry& fine_geometry,
                                     const amrex::IntVect& ratio,
                                     int fine_level,
                                     amrex::MultiFab& fine_candidate,
                                     std::string& diagnostic);

} // namespace erf_sbm

#endif
