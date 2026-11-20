#ifndef ERF_SBM_AMR_QUALIFICATION_H_
#define ERF_SBM_AMR_QUALIFICATION_H_

// Test-only support for the fixed two-level periodic SBM AMR qualification.
// Compiled only with ERF_ENABLE_SBM_AMR_QUALIFICATION; the ERF members it
// implements are declared in ERF.H.

#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>

struct SolverChoice;

namespace erf_sbm {

struct SBMAMRStageParent;
struct SBMLayoutSpec;

namespace qualification {

// Abort unless the inputs describe the supported fixed two-level fixture.
void validate_configuration (const SolverChoice& choice, int max_level);

// Attach the number-carried tracer used by the 2M qualification decks.
void add_layout_properties (SBMLayoutSpec& spec, const SolverChoice& choice);

// Compare received split-layout ghost values with their owning valid FABs.
void verify_split_ghost_values (const amrex::MultiFab& destination,
                                const amrex::MultiFab& fine_source_state,
                                const amrex::MultiFab& fine_source_conserved,
                                const amrex::MultiFab& coarse_old_spectrum,
                                const amrex::MultiFab& coarse_new_spectrum,
                                const SBMAMRStageParent& parent,
                                const amrex::Geometry& fine_geometry,
                                double theta,
                                int stage,
                                double requested_time);

// First locally owned cell marked 1 in a coarse-only ghost mask.
[[nodiscard]] amrex::IntVect first_coarse_only_ghost (const amrex::MultiFab& mask);

} // namespace qualification
} // namespace erf_sbm

#endif
