#ifndef ERF_RESOLVED_WALL_FLUX_H_
#define ERF_RESOLVED_WALL_FLUX_H_

#include "ERF_CloudChamberWallFlux.H"

namespace erf_resolved_wall_flux {

// Compatibility wrapper for callers that only provide the legacy resolved
// wall arguments. Production Cloud Chamber calls the generalized evaluator
// with stage-specific staggered velocity fields.
inline void apply (
    const amrex::Box& bx, const amrex::Box& domain, const int quantity,
    const int flux_comp, const amrex::Array4<const amrex::Real>& state,
    const amrex::Array4<const amrex::Real>& prim,
    const amrex::Array4<const amrex::Real>& base_state,
    const amrex::Array4<amrex::Real>& rhs,
    const amrex::Array4<amrex::Real>& xflux,
    const amrex::Array4<amrex::Real>& yflux,
    const amrex::Array4<amrex::Real>& zflux,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx_inv,
    const erf_wall_thermodynamics::Boundary& walls,
    const amrex::Real alpha_T, const amrex::Real alpha_C,
    const amrex::Real rdOcp)
{
    // This compatibility entry point has no velocity views. Keep its
    // precondition explicit so a future bulk wall cannot dereference an empty
    // Array4 through the resolved-only API.
    AMREX_ALWAYS_ASSERT(!erf_cloud_chamber_wall_flux::has_velocity_dependent_scalar_wall(walls));
    erf_cloud_chamber_wall_flux::apply(
        bx, domain, quantity, flux_comp, state, prim, base_state,
        amrex::Array4<const amrex::Real>{},
        amrex::Array4<const amrex::Real>{},
        amrex::Array4<const amrex::Real>{},
        rhs, xflux, yflux, zflux, dx_inv, walls, alpha_T, alpha_C, rdOcp);
}

} // namespace erf_resolved_wall_flux

#endif
