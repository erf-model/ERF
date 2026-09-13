#ifndef ERF_SBM_TRANSPORT_PROTOTYPE_H_
#define ERF_SBM_TRANSPORT_PROTOTYPE_H_

#include "ERF_AuxiliaryStageContext.H"
#include "ERF_AuxiliaryStateManager.H"
#include "ERF_SBMLayout.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

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
 * mass fluxes.  This is intentionally restricted to static Cartesian P1.
 */
void advance_stage(::erf_auxiliary::AuxiliaryStateManager& manager,
                   const SBMLayout& layout,
                   const ::erf_auxiliary::StageContext& context,
                   const amrex::MultiFab& rho_evaluation,
                   amrex::MultiFab& core_state,
                   const amrex::MultiFab& carrier_x,
                   const amrex::MultiFab& carrier_y,
                   const amrex::MultiFab& carrier_z,
                   const amrex::Geometry& geometry);

} // namespace erf_sbm

#endif
