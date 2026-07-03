#ifndef ERF_PLOTFILE2DFILL_H_
#define ERF_PLOTFILE2DFILL_H_

#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

// Mechanical fill helpers for ERF 2D plotfile assembly. These helpers copy,
// fill, and convert existing diagnostic sources into output components. They
// must not own science formulas beyond simple unit conversions delegated to
// existing diagnostic helpers.

namespace plotfile2d
{

void fill_component_with_value (amrex::MultiFab& dst,
                                int dst_comp,
                                amrex::Real value);

void fill_component_from_klevel (amrex::MultiFab& dst,
                                 int dst_comp,
                                 const amrex::MultiFab& src,
                                 int src_k,
                                 int src_comp = 0);

void fill_component_from_klevel_or_value (amrex::MultiFab& dst,
                                          int dst_comp,
                                          const amrex::MultiFab* src,
                                          int src_k,
                                          amrex::Real missing_value,
                                          int src_comp = 0);

void fill_sensible_heat_flux_from_klevel_or_missing (amrex::MultiFab& dst,
                                                     int dst_comp,
                                                     const amrex::MultiFab* src,
                                                     int src_k,
                                                     amrex::Real missing_value);

void fill_latent_heat_flux_from_klevel_or_missing (amrex::MultiFab& dst,
                                                   int dst_comp,
                                                   const amrex::MultiFab* src,
                                                   int src_k,
                                                   amrex::Real missing_value);

} // namespace plotfile2d

#endif
