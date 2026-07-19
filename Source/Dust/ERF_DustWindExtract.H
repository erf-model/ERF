/**
 * @file ERF_DustWindExtract.H
 * @brief Wind and surface field extraction from ERF 3D solver to the 2D dust grid.
 *
 * The wind interpolation algorithm copies fill_fire_wind_from_interpolation
 * from Source/Fire/ERF_FireWindExtract.cpp, substituting DustGrid for FireGrid
 * and grid_ratio for C. Per-fuel and WAF parameters are not applicable to dust.
 *
 * SurfaceLayer accessors used:
 *   get_u_star(0)  — friction velocity [m/s]
 *   get_t_surf(0)  — surface skin temperature [K]
 *   get_pblh(0)    — PBL height [m]
 *
 * References:
 *   Source/Fire/ERF_FireWindExtract.H (algorithm source):
 *     https://github.com/hgopalan/ERF/blob/a2ca3d1a7cea08a480c9d10e93ffdfab1b24b24e/
 *     Source/Fire/ERF_FireWindExtract.H
 *   Hong & Pan (1996). Mon. Wea. Rev., 124, 2322.
 *     https://doi.org/10.1175/1520-0493(1996)124<2322:NBLVDI>2.0.CO;2
 *   Marticorena & Bergametti (1995). J. Geophys. Res., 100, 16415.
 *     https://doi.org/10.1029/95JD00690
 */

#ifndef ERF_DUST_WIND_EXTRACT_H
#define ERF_DUST_WIND_EXTRACT_H

#ifdef ERF_USE_DUST

#include <AMReX_MultiFab.H>
#include <ERF_DustGrid.H>

// Vertical wind interpolation to zref. Copy of fill_fire_wind_from_interpolation
// with FireGrid->DustGrid and fg.C->dg.grid_ratio.
// Steps per dust cell: map to atm column, find k bracket at z_surf+zref,
// average face-staggered u/v to cell centres, linearly interpolate to z_target.
void fill_dust_wind_from_interpolation(
    amrex::MultiFab&       dust_wind_ref,  // ncomp=2, u and v [m/s]
    const amrex::MultiFab& xvel_mf,
    const amrex::MultiFab& yvel_mf,
    const amrex::MultiFab& z_phys_cc_mf,
    const DustGrid&        dg,
    amrex::Real            zref,
    int                    nz);

// Map u* from SurfaceLayer::get_u_star(0) to dust grid (i_a=i/C, j_a=j/C, k=0).
void fill_dust_ustar_from_surface_layer(
    amrex::MultiFab&       dust_ustar_in,
    const amrex::MultiFab& ustar_atm,
    const DustGrid&        dg);

// Generic: map any scalar atm k=0 field to dust grid.
// Used for T_sfc (get_t_surf) and PBLH (get_pblh).
void fill_dust_scalar_from_atm(
    amrex::MultiFab&       dust_field,
    const amrex::MultiFab& atm_field,
    const DustGrid&        dg);

// Compute u* from wind magnitude using simple log-profile: u* = κ * U / ln(z_ref/z0)
// @param dust_ustar_in [out] Friction velocity on dust grid [m/s]
// @param dust_wind_ref [in]  Wind vector (u,v) at reference height on dust grid [m/s]
// @param z_ref        [in]  Reference height [m]
// @param z0           [in]  Surface roughness length [m]
void compute_dust_ustar_from_wind(
    amrex::MultiFab&       dust_ustar_in,
    const amrex::MultiFab& dust_wind_ref,
    amrex::Real            z_ref,
    amrex::Real            z0);

#endif // ERF_USE_DUST

#endif // ERF_DUST_WIND_EXTRACT_H
