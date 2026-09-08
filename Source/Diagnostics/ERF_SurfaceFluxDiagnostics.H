#ifndef ERF_SURFACE_FLUX_DIAGNOSTICS_H_
#define ERF_SURFACE_FLUX_DIAGNOSTICS_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include "ERF_Constants.H"

namespace surface_flux_diagnostics
{

/**
 * Convert applied conservative surface fluxes to W m^-2 for 2D plotfile output.
 *
 * `SFS_hfx3` and `SFS_q1fx3` are already conservative applied fluxes. Noah-MP
 * LSM kinematic fluxes are converted to conservative fluxes in
 * `SurfaceLayer::compute_SurfaceLayer_bcs`, so this layer must not multiply by
 * density. The output sign convention follows ERF's lower-boundary flux
 * convention. These helpers exist so the 2D writer does not duplicate unit
 * semantics.
 *
 * Maintenance note:
 * Do not add source-specific Noah-MP/MOST branches here. If this helper ever
 * appears to need a density factor, the applied SurfaceLayer flux contract has
 * likely been broken upstream.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
sensible_heat_flux_wm2_from_rhotheta_flux (amrex::Real rhotheta_flux) noexcept
{
    return Cp_d * rhotheta_flux;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
latent_heat_flux_wm2_from_rhoqv_flux (amrex::Real rhoqv_flux) noexcept
{
    return L_v * rhoqv_flux;
}

} // namespace surface_flux_diagnostics

#endif
