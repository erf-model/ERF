#ifndef _STENCILS_H_
#define _STENCILS_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_InterpFaceRegister.H>
#include "DataStruct.H"
#include "IndexDefines.H"

using namespace amrex;

void MomentumToVelocity(MultiFab& xvel_out, MultiFab& yvel_out, MultiFab& zvel_out,
                        const MultiFab& cons_in,
                        const MultiFab& xmom_in, const MultiFab& ymom_in, const MultiFab& zmom_in,
                        const IntVect& ngrow);

void VelocityToMomentum(const MultiFab& xvel_in , const MultiFab& yvel_in, const MultiFab& zvel_in,
                        const MultiFab& cons_in,
                        MultiFab& xmom_out, MultiFab& ymom_out, MultiFab& zmom_out,
                        const IntVect& ngrow);


AMREX_GPU_DEVICE
Real InterpolateFromCellOrFace(const int &i, const int &j, const int &k,
                           const Array4<const Real>& qty, const int & qty_index,
                           const Real& uadv,
                           const enum Coord& coordDir,
                           const int& spatial_order);

AMREX_GPU_DEVICE
Real InterpolatePertFromCell(const int &i, const int &j, const int &k,
                           const Array4<const Real>& qty, const int & qty_index,
                           const Real& uadv,
                           const enum Coord& coordDir,
                           const int& spatial_order,
                           const amrex::Real* dptr_hse);

AMREX_GPU_DEVICE
Real InterpolateDensityPertFromCellToFace(const int &i, const int &j, const int &k,
                               const Array4<const Real>& cons_in,
                               const Real& uadv,
                               const enum Coord& coordDir,
                               const int& spatial_order,
                               const amrex::Real* dptr_hse);

/** Meant for {x, y, z}- momentum equations */
AMREX_GPU_DEVICE
Real AdvectionContributionForXMom(const int &i, const int &j, const int &k,
                                  const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                                  const Array4<const Real>& u,
                                  const GpuArray<Real, AMREX_SPACEDIM>& cellSize,
                                  const int& spatial_order);
AMREX_GPU_DEVICE
Real AdvectionContributionForYMom(const int &i, const int &j, const int &k,
                                  const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                                  const Array4<const Real>& v,
                                  const GpuArray<Real, AMREX_SPACEDIM>& cellSize,
                                  const int& spatial_order);
AMREX_GPU_DEVICE
Real AdvectionContributionForZMom(const int &i, const int &j, const int &k,
                                  const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                                  const Array4<const Real>& w,
                                  const GpuArray<Real, AMREX_SPACEDIM>& cellSize,
                                  const int& spatial_order);

/** Meant for the continuity, energy, and scalar equations */
AMREX_GPU_DEVICE
Real AdvectionContributionForState(const int &i, const int &j, const int &k,
                                   const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                                   const Array4<const Real>& cell_prim, const int &qty_index,
                                   const Array4<Real>& flux_x, const Array4<Real>& flux_y,
                                   const Array4<Real>& flux_z,
                                   const GpuArray<Real, AMREX_SPACEDIM>& cellSize,
                                   const int& spatial_order);

/** Compute K (i-1/2, j+1/2, k) etc given K(i, j, k) or nut (i, j, k) is known */
AMREX_GPU_DEVICE
Real InterpolateTurbulentViscosity(const int &i, const int &j, const int &k,
                            const enum MomentumEqn &momentumEqn,
                            const enum DiffusionDir &diffDir,
                            const Array4<Real>& Ksmag);

/** Compute tau_ij (m + 1/2), tau_ij (m - 1/2) where m = {i, j, k} for DNS or Smagorinsky */
AMREX_GPU_DEVICE
Real ComputeStressTerm (const int &i, const int &j, const int &k,
                        const Array4<const Real>& u, const Array4<const Real>& v, const Array4<const Real>& w,
                        const enum MomentumEqn &momentumEqn,
                        const enum DiffusionDir &diffDir,
                        const GpuArray<Real, AMREX_SPACEDIM>& cellSize,
                        const Array4<Real>& Ksmag,
                        const SolverChoice &solverChoice);

AMREX_GPU_DEVICE
Real DiffusionContributionForMom(const int &i, const int &j, const int &k,
                                 const Array4<const Real>& u, const Array4<const Real>& v, const Array4<const Real>& w,
                                 const Array4<const Real>& cons,
                                 const enum MomentumEqn &momentumEqn,
                                 const GpuArray<Real, AMREX_SPACEDIM>& cellSize,
                                 const Array4<Real>& Ksmag,
                                 const SolverChoice &solverChoice,
                                 bool dirichlet_at_lo_k,
                                 bool dirichlet_at_hi_k);

AMREX_GPU_DEVICE
Real ComputeDiffusionFluxForState(const int &i, const int &j, const int &k,
                                  const Array4<const Real>& cell_data,
                                  const Array4<const Real>& cell_prim, const int & prim_index,
                                  const amrex::Real invCellWidth,
                                  const Array4<Real>& Ksmag,
                                  const SolverChoice &solverChoice,
                                  const enum Coord& coordDir);

AMREX_GPU_DEVICE
Real DiffusionContributionForState(const int &i, const int &j, const int &k,
                                   const Array4<const Real>& cell_data,
                                   const Array4<const Real>& cell_prim, const int & qty_index,
                                   const Array4<Real>& flux_x, const Array4<Real>& flux_y,
                                   const Array4<Real>& flux_z,
                                   const GpuArray<Real, AMREX_SPACEDIM>& cellSize,
                                   const Array4<Real>& Ksmag,
                                   const SolverChoice &solverChoice);
#endif
