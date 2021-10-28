#ifndef _RK3_H_
#define _RK3_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <ERF.H>
#include "IndexDefines.H"

using namespace amrex;

// cons: cell-centered state
// U,V,W: face-based velocity
// source: cell-centered source term

void MomentumToVelocity(MultiFab& xvel_out, MultiFab& yvel_out, MultiFab& zvel_out,
                        MultiFab& cons_in,
                        MultiFab& xmom_in, MultiFab& ymom_in, MultiFab& zmom_in,
                        const SolverChoice& solverChoice);

void VelocityToMomentum(MultiFab& xvel_in , MultiFab& yvel_in, MultiFab& zvel_in,
                        MultiFab& cons_in,
                        MultiFab& xmom_out, MultiFab& ymom_out, MultiFab& zmom_out,
                        const SolverChoice& solverChoice);

void RK3_advance(MultiFab& cons_old, MultiFab& cons_new,
                 MultiFab& U_old, MultiFab& V_old, MultiFab& W_old,
                 MultiFab& U_new, MultiFab& V_new, MultiFab& W_new,
                 MultiFab& source,
                 MultiFab& eta, MultiFab& zeta, MultiFab& kappa,
                 std::array<MultiFab, AMREX_SPACEDIM>& faceflux,
                 std::array< MultiFab, 2 >& edgeflux_x,
                 std::array< MultiFab, 2 >& edgeflux_y,
                 std::array< MultiFab, 2 >& edgeflux_z,
                 std::array< MultiFab, AMREX_SPACEDIM>& cenflux,
                 const amrex::Geometry geom, const amrex::Real* dxp, const amrex::Real dt,
                 const SolverChoice& solverChoice);

void RK3_stage  (MultiFab& cons_old, MultiFab& cons_new,
                 MultiFab& xmom_old, MultiFab& ymom_old, MultiFab& zmom_old,
                 MultiFab& xmom_new, MultiFab& ymom_new, MultiFab& zmom_new,
                 MultiFab& xvel    , MultiFab& yvel    , MultiFab& zvel    ,
                 MultiFab& prim    , MultiFab& source,
                 MultiFab& eta, MultiFab& zeta, MultiFab& kappa,
                 std::array<MultiFab, AMREX_SPACEDIM>& faceflux,
                 std::array< MultiFab, 2 >& edgeflux_x,
                 std::array< MultiFab, 2 >& edgeflux_y,
                 std::array< MultiFab, 2 >& edgeflux_z,
                 std::array< MultiFab, AMREX_SPACEDIM>& cenflux,
                 const amrex::Geometry geom, const amrex::Real* dxp, const amrex::Real dt,
                 const SolverChoice& solverChoice);

AMREX_GPU_DEVICE
Real InterpolateFromCellOrFace(const int &i, const int &j, const int &k,
                           const Array4<Real>& qty, const int & qty_index,
                           const enum NextOrPrev &nextOrPrev,
                           const enum Coord& coordDir,
                           const int &spatial_order);

AMREX_GPU_DEVICE
Real InterpolateDensityFromCellToFace(const int &i, const int &j, const int &k,
                               const Array4<Real>& cons_in,
                               const enum NextOrPrev &nextOrPrev,
                               const enum Coord& coordDir,
                               const int &spatial_order);

AMREX_GPU_DEVICE
Real InterpolateRhoThetaFromCellToFace(const int &i, const int &j, const int &k,
                                      const Array4<Real>& cons_in,
                                      const enum NextOrPrev &nextOrPrev,
                                      const enum Coord& coordDir,
                                      const int &spatial_order);

AMREX_GPU_DEVICE
Real InterpolateRhoScalarFromCellToFace(const int &i, const int &j, const int &k,
                                      const Array4<Real>& cons_in,
                                      const enum NextOrPrev &nextOrPrev,
                                      const enum Coord& coordDir,
                                      const int &spatial_order);

// Meant for {x, y, z}- momentum equations
AMREX_GPU_DEVICE
Real ComputeAdvectedQuantityForMom(const int &i, const int &j, const int &k,
                                   const Array4<Real>& rho_u, const Array4<Real>& rho_v, const Array4<Real>& rho_w,
                                   const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                                   const enum NextOrPrev &nextOrPrev,
                                   const enum AdvectedQuantity &advectedQuantity,
                                   const enum AdvectingQuantity &advectingQuantity,
                                   const int &spatial_order);

// Meant for {x, y, z}- momentum equations
AMREX_GPU_DEVICE
Real AdvectionContributionForMom(const int &i, const int &j, const int &k,
                                 const Array4<Real>& rho_u, const Array4<Real>& rho_v, const Array4<Real>& rho_w,
                                 const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                                 const enum MomentumEqn &momentumEqn,
                                 const GpuArray<Real, AMREX_SPACEDIM>& cellSize,
                                 const SolverChoice &solverChoice);

// Meant for the continuity, energy, and scalar equations
AMREX_GPU_DEVICE
Real ComputeAdvectedQuantityForState(const int &i, const int &j, const int &k,
                             const Array4<Real>& rho_u, const Array4<Real>& rho_v, const Array4<Real>& rho_w,
                             const Array4<Real>& cell_data,
                             const enum NextOrPrev &nextOrPrev,
                             const enum AdvectedQuantity &advectedQuantity,
                             const enum AdvectingQuantity &advectingQuantity,
                             const int &spatial_order);

// Meant for the continuity, energy, and scalar equations
AMREX_GPU_DEVICE
Real AdvectionContributionForState(const int &i, const int &j, const int &k,
                                   const Array4<Real>& rho_u, const Array4<Real>& rho_v, const Array4<Real>& rho_w,
                                   const Array4<Real>& cell_data, const int &qty_index,
                                   const GpuArray<Real, AMREX_SPACEDIM>& cellSize,
                                   const int &spatial_order);


// Compute Sij
AMREX_GPU_DEVICE
Real ComputeStrainRate(const int &i, const int &j, const int &k,
                        const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                        const enum NextOrPrev &nextOrPrev,
                        const enum MomentumEqn &momentumEqn,
                        const enum DiffusionDir &diffDir,
                        const GpuArray<Real, AMREX_SPACEDIM>& cellSize);

// Given a velocity field and a geometry, compute nut or K(i, j, k) at cell-centers
//AMREX_GPU_DEVICE
void ComputeTurbulentViscosity(MultiFab& xvel, MultiFab& yvel, MultiFab& zvel,
                               MultiFab& cons_in, MultiFab& eddyViscosity,
                               const GpuArray<Real, AMREX_SPACEDIM>& cellSize,
                               const SolverChoice& solverChoice);

// Compute K (i-1/2, j+1/2, k) etc given K(i, j, k) or nut (i, j, k) is known
AMREX_GPU_DEVICE
Real InterpolateTurbulentViscosity(const int &i, const int &j, const int &k,
                            const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                            const enum NextOrPrev &nextOrPrev,
                            const enum MomentumEqn &momentumEqn,
                            const enum DiffusionDir &diffDir,
                            const GpuArray<Real, AMREX_SPACEDIM>& cellSize,
                            const Array4<Real>& nut);

// Compute tau_ij (m + 1/2), tau_ij (m - 1/2) where m = {i, j, k} for DNS or Smagorinsky
AMREX_GPU_DEVICE
Real ComputeStressTerm (const int &i, const int &j, const int &k,
                        const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                        const enum NextOrPrev &nextOrPrev,
                        const enum MomentumEqn &momentumEqn,
                        const enum DiffusionDir &diffDir,
                        const GpuArray<Real, AMREX_SPACEDIM>& cellSize,
                        const Array4<Real>& nut,
                        const SolverChoice &solverChoice);

AMREX_GPU_DEVICE
Real DiffusionContributionForMom(const int &i, const int &j, const int &k,
                                 const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                                 const enum MomentumEqn &momentumEqn,
                                 const GpuArray<Real, AMREX_SPACEDIM>& cellSize,
                                 const Array4<Real>& nut,
                                 const SolverChoice &solverChoice);

AMREX_GPU_DEVICE
Real ComputeDiffusionTermForState(const int &i, const int &j, const int &k,
                                  const Array4<Real>& cell_data, const int & qty_index,
                                  const enum Coord& coordDir);

AMREX_GPU_DEVICE
Real DiffusionContributionForState(const int &i, const int &j, const int &k,
                                   const Array4<Real>& cell_data, const int & qty_index,
                                   const GpuArray<Real, AMREX_SPACEDIM>& cellSize,
                                   const SolverChoice &solverChoice);

int ComputeGhostCells(const int &spatial_order);

#endif
