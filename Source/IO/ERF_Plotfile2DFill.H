/**
 * \file ERF_Plotfile2DFill.H
 */
#ifndef ERF_PLOTFILE2DFILL_H_
#define ERF_PLOTFILE2DFILL_H_

#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

// Mechanical fill helpers for ERF 2D plotfile assembly. These helpers copy,
// fill, and convert existing diagnostic sources into output components.
//
// Preconditions:
//   - dst is the already-allocated 2D plotfile MultiFab for one AMR level.
//   - source MultiFabs use a compatible BoxArray/DistributionMapping layout.
//   - src_k is interpreted in the source MultiFab's vertical index space.
//
// Non-goals:
//   - Do not add diagnostic science here.
//   - Do not encode runtime source selection here.
//   - Do not use this layer as a registry for all future diagnostics.
//
// Column reductions, such as integrated_qv and scheme-aware condensed
// water-path diagnostics, need separate helpers because they reduce over k
// rather than copying one k level. Future pressure-surface diagnostics need
// interpolation-specific helpers.

namespace plotfile2d
{

// Fill one destination component with a scalar value. Used for documented
// missing values and zero fills.
void fill_component_with_value (amrex::MultiFab& dst,
                                int dst_comp,
                                amrex::Real value);

// Copy one source component at one source k level into one 2D destination
// component. The caller owns source selection and chooses src_k.
void fill_component_from_klevel (amrex::MultiFab& dst,
                                 int dst_comp,
                                 const amrex::MultiFab& src,
                                 int src_k,
                                 int src_comp = 0);

// Copy from src when present; otherwise fill the destination component with the
// documented missing value. This keeps optional diagnostics branch-local in the
// writer.
void fill_component_from_klevel_or_value (amrex::MultiFab& dst,
                                          int dst_comp,
                                          const amrex::MultiFab* src,
                                          int src_k,
                                          amrex::Real missing_value,
                                          int src_comp = 0);

// Noah-MP and other land-surface providers use internal sentinels that must
// not enter a public plotfile. Valid zero and negative physical values pass
// through unchanged.
void fill_land_surface_component_from_klevel_or_missing (amrex::MultiFab& dst,
                                                         int dst_comp,
                                                         const amrex::MultiFab* src,
                                                         int src_k,
                                                         amrex::Real missing_value);

// Convert the selected conservative sensible-flux source to W m^-2 for the
// public sensible_heat_flux diagnostic. Unit semantics remain in
// surface_flux_diagnostics.
void fill_sensible_heat_flux_from_klevel_or_missing (amrex::MultiFab& dst,
                                                     int dst_comp,
                                                     const amrex::MultiFab* src,
                                                     int src_k,
                                                     amrex::Real missing_value);

// Convert the selected conservative moisture-flux source to W m^-2 for the
// public latent_heat_flux diagnostic. Unit semantics remain in
// surface_flux_diagnostics.
void fill_latent_heat_flux_from_klevel_or_missing (amrex::MultiFab& dst,
                                                   int dst_comp,
                                                   const amrex::MultiFab* src,
                                                   int src_k,
                                                   amrex::Real missing_value);

} // namespace plotfile2d

#endif
