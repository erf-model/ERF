/**
 * @file ERF_DustAtmCoupling.H
 * @brief One-way dust-to-atmosphere coupling: inject dust emission flux
 *        into the ERF 3D conserved state at k=0.
 *
 * Structural analog of ERF_FireAtmCoupling.H from Source/Fire/:
 *   coarsen_fire_flux_to_atm  ->  coarsen_dust_flux_to_atm
 *   apply_fire_tendency_to_cc_source  ->  apply_dust_tendency_to_cc_source
 *
 * The fire module distributes heat flux vertically via an exponential
 * profile (Mandel et al. 2011). Dust injection is a surface-only source:
 * the emission flux is added to the k=0 layer and transported by the
 * dycore advection in subsequent steps.
 *
 * Injection formula at k=0:
 *   d(RhoDust)/dt = F_dust / dz_k0
 * where F_dust [kg/m^2/s] is the coarsened total emission flux and
 * dz_k0 [m] is the depth of the first atmospheric cell.
 *
 * One-step explicit lag: flux from time n is injected at step n+1,
 * before advance_dycore(). Consistent with ERF_FireAtmCoupling.H.
 *
 * Scalar slot: RhoAdv_comp (first passive advected scalar). Phase 11
 * extends to n_size_bins slots when DustParams::transport_bins_separately
 * is true.
 *
 * References:
 *   Mandel, J. et al. (2011). Math. Comput. Simul., 79(3), 584-606.
 *     (WRF-SFIRE fire_tendency coupling approach)
 *   Marticorena, B., & Bergametti, G. (1995). J. Geophys. Res., 100, 16415.
 *     https://doi.org/10.1029/95JD00690
 */

#ifndef ERF_DUST_ATM_COUPLING_H
#define ERF_DUST_ATM_COUPLING_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_MFIter.H>

// Coarsen dust_emission_flux from dust grid to atm level-0 grid.
// Uses amrex::average_down with refinement ratio (grid_ratio, grid_ratio, 1).
// Direct analog of coarsen_fire_flux_to_atm in ERF_FireAtmCoupling.H.
inline void coarsen_dust_flux_to_atm(
    amrex::MultiFab&       Q_dust_atm,     // output: ncomp=1, on atm grid
    const amrex::MultiFab& Q_dust,         // input:  dust_emission_flux summed over bins
    const amrex::Geometry& geom_dust,
    const amrex::Geometry& geom_atm,
    int                    grid_ratio)
{
    // Sum over bins is done before calling this function.
    // grid_ratio = 1: direct copy (no averaging needed).
    amrex::average_down(Q_dust, Q_dust_atm, geom_dust, geom_atm,
                        0, 1, amrex::IntVect(grid_ratio, grid_ratio, 1));
}

// Apply dust emission flux tendency to cc_source at k=0.
// All levels above k=0 receive zero tendency from this function.
// Transport away from k=0 is handled by the dycore advection.
//
// Tendency: d(RhoAdv_comp)/dt|_{k=0} += F_dust / dz_k0
//
// Optional dust_debug prints diagnostics after the kernel:
//   [DUST COUPLING] F_dust_max=... kg/m^2/s  RhoDust_tend_max=... kg/m^3/s
void apply_dust_tendency_to_cc_source(
    amrex::MultiFab&       cc_source,
    const amrex::MultiFab& Q_dust_atm,    // coarsened emission flux [kg/m^2/s]
    const amrex::MultiFab& z_phys_cc,     // cell-centre terrain height [m]
    const amrex::Geometry& geom_atm,
    int                    dust_scalar_comp, // scalar component index in cc_source
    amrex::Real            feedback,         // [0,1]; 0 = off, 1 = full coupling
    bool                   dust_debug = false);

#endif // ERF_DUST_ATM_COUPLING_H
