#ifndef _INTEGRATION_H_
#define _INTEGRATION_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_InterpFaceRegister.H>
#include "DataStruct.H"
#include "IndexDefines.H"

using namespace amrex;

// cons: cell-centered state
// U,V,W: face-based velocity
// source: cell-centered source term

namespace IntVar {
    enum {
        cons = 0,
        xmom,
        ymom,
        zmom,
        xflux,
        yflux,
        zflux
    };
}

void MomentumToVelocity(MultiFab& xvel_out, MultiFab& yvel_out, MultiFab& zvel_out,
                        const MultiFab& cons_in,
                        const MultiFab& xmom_in, const MultiFab& ymom_in, const MultiFab& zmom_in,
                        int ngrow, int l_spatial_order);

void VelocityToMomentum(const MultiFab& xvel_in , const MultiFab& yvel_in, const MultiFab& zvel_in,
                        const MultiFab& cons_in,
                        MultiFab& xmom_out, MultiFab& ymom_out, MultiFab& zmom_out,
                        int l_spatial_order);

/**
 * Update the cell-centered state and face-based velocity using a time integrator
 *
 * \param[in]  S_old  state on cell centers
 * \param[in]  U_old  x-velocity on x-faces
 * \param[in]  V_old  y-velocity on y-faces
 * \param[in]  W_old  z-velocity on z-faces
 * \param[in]  source source term on cell centers
 *
 * \param[out] S_new  state on cell centers
 * \param[out] U_new  x-velocity on x-faces
 * \param[out] V_new  y-velocity on y-faces
 * \param[out] W_new  z-velocity on z-faces
 *
 */
void erf_advance(int level,
                 MultiFab& cons_old, MultiFab& cons_new,
                 MultiFab& U_old, MultiFab& V_old, MultiFab& W_old,
                 MultiFab& U_new, MultiFab& V_new, MultiFab& W_new,
                 MultiFab& xmom_crse, MultiFab& ymom_crse, MultiFab& zmom_crse,
                 MultiFab& source,
                 std::array<MultiFab, AMREX_SPACEDIM>& flux,
                 const amrex::Geometry crse_geom,
                 const amrex::Geometry fine_geom,
                 const amrex::IntVect ref_ratio,
                 const amrex::Real dt,
                 const amrex::Real time,
                       amrex::InterpFaceRegister* ifr,
                 const SolverChoice& solverChoice,
                 const amrex::Real* dptr_dens_hse,
                 const amrex::Real* dptr_pres_hse,
                 const amrex::Real* dptr_rayleigh_tau,
                 const amrex::Real* dptr_rayleigh_ubar,
                 const amrex::Real* dptr_rayleigh_vbar,
                 const amrex::Real* dptr_rayleigh_thetabar);

void erf_rhs(int level,
             Vector<MultiFab>& S_rhs,
             const Vector<MultiFab>& S_data,
             MultiFab& source,
             std::array< MultiFab, AMREX_SPACEDIM>&  advflux,
             std::array< MultiFab, AMREX_SPACEDIM>& diffflux,
             const amrex::Geometry geom, const amrex::Real dt,
                   amrex::InterpFaceRegister* ifr,
             const SolverChoice& solverChoice,
             const bool lo_z_is_no_slip,
             const bool hi_z_is_no_slip,
             const amrex::Real* dptr_dens_hse,
             const amrex::Real* dptr_pres_hse,
             const amrex::Real* dptr_rayleigh_tau,
             const amrex::Real* dptr_rayleigh_ubar,
             const amrex::Real* dptr_rayleigh_vbar,
             const amrex::Real* dptr_rayleigh_thetabar);

AMREX_GPU_DEVICE
Real InterpolateFromCellOrFace(const int &i, const int &j, const int &k,
                           const Array4<const Real>& qty, const int & qty_index,
                           const enum Coord& coordDir,
                           const int &spatial_order);

AMREX_GPU_DEVICE
Real InterpolatePertFromCell(const int &i, const int &j, const int &k,
                           const Array4<const Real>& qty, const int & qty_index,
                           const enum Coord& coordDir,
                           const int &spatial_order,
                           const amrex::Real* dptr_hse);

AMREX_GPU_DEVICE
Real InterpolateDensityPertFromCellToFace(const int &i, const int &j, const int &k,
                               const Array4<const Real>& cons_in,
                               const enum Coord& coordDir,
                               const int &spatial_order,
                               const amrex::Real* dptr_hse);

/** Meant for {x, y, z}- momentum equations */
AMREX_GPU_DEVICE
Real ComputeAdvectedQuantityForMom(const int &i, const int &j, const int &k,
                                   const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                                   const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                                   const enum AdvectedQuantity &advectedQuantity,
                                   const enum AdvectingQuantity &advectingQuantity,
                                   const int &spatial_order);

/** Meant for {x, y, z}- momentum equations */
AMREX_GPU_DEVICE
Real AdvectionContributionForMom(const int &i, const int &j, const int &k,
                                 const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                                 const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                                 const enum MomentumEqn &momentumEqn,
                                 const GpuArray<Real, AMREX_SPACEDIM>& cellSize,
                                 const SolverChoice &solverChoice);

/** Meant for the continuity, energy, and scalar equations */
AMREX_GPU_DEVICE
Real ComputeAdvectedQuantityForState(const int &i, const int &j, const int &k,
                             const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                             const Array4<Real>& cell_data,
                             const enum AdvectedQuantity &advectedQuantity,
                             const enum AdvectingQuantity &advectingQuantity,
                             const int &spatial_order);

/** Meant for the continuity, energy, and scalar equations */
AMREX_GPU_DEVICE
Real AdvectionContributionForState(const int &i, const int &j, const int &k,
                                   const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                                   const Array4<const Real>& cell_data, const int &qty_index,
                                   const Array4<Real>& flux_x, const Array4<Real>& flux_y,
                                   const Array4<Real>& flux_z,
                                   const GpuArray<Real, AMREX_SPACEDIM>& cellSize,
                                   const int &spatial_order);

/** Compute K (i-1/2, j+1/2, k) etc given K(i, j, k) or nut (i, j, k) is known */
AMREX_GPU_DEVICE
Real InterpolateTurbulentViscosity(const int &i, const int &j, const int &k,
                            const enum MomentumEqn &momentumEqn,
                            const enum DiffusionDir &diffDir,
                            const Array4<Real>& Ksmag);

/** Compute tau_ij (m + 1/2), tau_ij (m - 1/2) where m = {i, j, k} for DNS or Smagorinsky */
AMREX_GPU_DEVICE
Real ComputeStressTerm (const int &i, const int &j, const int &k,
                        const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                        const enum MomentumEqn &momentumEqn,
                        const enum DiffusionDir &diffDir,
                        const GpuArray<Real, AMREX_SPACEDIM>& cellSize,
                        const Array4<Real>& Ksmag,
                        const SolverChoice &solverChoice);

AMREX_GPU_DEVICE
Real DiffusionContributionForMom(const int &i, const int &j, const int &k,
                                 const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                                 const enum MomentumEqn &momentumEqn,
                                 const GpuArray<Real, AMREX_SPACEDIM>& cellSize,
                                 const Array4<Real>& Ksmag,
                                 const SolverChoice &solverChoice,
                                 bool use_no_slip_stencil_at_lo_k,
                                 bool use_no_slip_stencil_at_hi_k);

AMREX_GPU_DEVICE
Real ComputeDiffusionFluxForState(const int &i, const int &j, const int &k,
                                  const Array4<const Real>& cell_data, const int & qty_index,
                                  const amrex::Real invCellWidth,
                                  const Array4<Real>& Ksmag,
                                  const SolverChoice &solverChoice,
                                  const enum Coord& coordDir);

AMREX_GPU_DEVICE
Real DiffusionContributionForState(const int &i, const int &j, const int &k,
                                   const Array4<const Real>& cell_data, const int & qty_index,
                                   const Array4<Real>& flux_x, const Array4<Real>& flux_y,
                                   const Array4<Real>& flux_z,
                                   const GpuArray<Real, AMREX_SPACEDIM>& cellSize,
                                   const Array4<Real>& Ksmag,
                                   const SolverChoice &solverChoice);

int ComputeGhostCells(const int &spatial_order);

#endif
