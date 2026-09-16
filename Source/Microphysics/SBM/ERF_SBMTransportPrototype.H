#ifndef ERF_SBM_TRANSPORT_PROTOTYPE_H_
#define ERF_SBM_TRANSPORT_PROTOTYPE_H_

#include "ERF_AuxiliaryStageContext.H"
#include "ERF_AuxiliaryFaceTransfer.H"
#include "ERF_AuxiliaryStateManager.H"
#include "ERF_SBMLayout.H"
#include "ERF_SBMFCT.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <cstddef>
#include <vector>

namespace erf_sbm {

using HostState = std::vector<amrex::Real>;

/** Exact reduced-state recurrence used by the production auxiliary hook. */
void update_stage(const ::erf_auxiliary::StageContext& context,
                  const HostState& old_state,
                  const HostState& evaluation_state,
                  const HostState& rhs,
                  HostState& output_state);

/** Return the accepted full-step ledger for a sequence of stage face fluxes. */
HostState accepted_ledger(const ::erf_auxiliary::StageContext& context,
                          const std::vector<HostState>& stage_fluxes,
                          double full_step);

/**
 * First-order donor transport using ERF's already-constructed dry-air carrier
 * mass fluxes.  This is intentionally restricted to the qualified static
 * Cartesian P1/P2 auxiliary contract.
 */
void advance_stage(::erf_auxiliary::AuxiliaryStateManager& manager,
                   const SBMLayout& layout,
                   const ::erf_auxiliary::StageContext& context,
                   const amrex::MultiFab& rho_evaluation,
                   amrex::MultiFab& core_state,
                   const amrex::MultiFab& carrier_x,
                   const amrex::MultiFab& carrier_y,
                   const amrex::MultiFab& carrier_z,
                   const amrex::Geometry& geometry,
                   ::erf_auxiliary::AuxiliaryFaceTransfer& stage_flux,
                   TransportMethod method = TransportMethod::DonorCell,
                   int level = 0,
                   amrex::Real diffusion_coefficient = amrex::Real(0.0),
                   int chunk_size = 0,
                   amrex::Real* minimum_accepted_limiter = nullptr);

/** Fail-closed finite/nonnegative check for authoritative spectral state. */
void validate_nonnegative_state(const amrex::MultiFab& state, int ncomp,
                                const char* context = "SBM auxiliary state");

/** Validate physical storage and two-moment endpoint realizability. */
void validate_admissible_state(::erf_auxiliary::AuxiliaryStateManager& manager,
                               const SBMLayout& layout, int level);

/**
 * Return the allocation-bound temporary payload for one production grouped
 * FCT chunk.  The estimate includes grown ratio/budget FABs and the three
 * chunk-local face candidates; authoritative state and accepted ledgers are
 * intentionally excluded because they are persistent state, not temporary
 * chunk work.
 */
[[nodiscard]] std::size_t grouped_fct_peak_working_bytes(
    const SBMLayout& layout, const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm, int chunk_size);

} // namespace erf_sbm

#endif
