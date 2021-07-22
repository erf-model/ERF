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

void CalcAdvFlux (const MultiFab& cons_in, 
                  const MultiFab& cons_x, const MultiFab& cons_y, const MultiFab& cons_z, 
                  const MultiFab& u_x, const MultiFab& v_y, const MultiFab& w_z, 
                  std::array<MultiFab, AMREX_SPACEDIM>& faceflux_in,
                  std::array< MultiFab, 2 >& edgeflux_x_in,
                  std::array< MultiFab, 2 >& edgeflux_y_in,
                  std::array< MultiFab, 2 >& edgeflux_z_in,
                  std::array< MultiFab, AMREX_SPACEDIM>& cenflux_in,
                  const amrex::Geometry geom,
                  const amrex::Real* dx, const amrex::Real dt,
                  const SolverChoice& solverChoice);

void CalcDiffFlux(const MultiFab& cons_in, 
                  const MultiFab& cons_x, const MultiFab& cons_y, const MultiFab& cons_z, 
                  const MultiFab& u_x, const MultiFab& v_y, const MultiFab& w_z, 
                  const MultiFab& eta_in, const MultiFab& zeta_in, const MultiFab& kappa_in,
                  std::array<MultiFab, AMREX_SPACEDIM>& faceflux_in,
                  std::array< MultiFab, 2 >& edgeflux_x_in,
                  std::array< MultiFab, 2 >& edgeflux_y_in,
                  std::array< MultiFab, 2 >& edgeflux_z_in,
                  std::array< MultiFab, AMREX_SPACEDIM>& cenflux_in,
                  const amrex::Geometry geom,
                  const amrex::Real* dx, const amrex::Real dt,
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
                                 const amrex::Geometry &geom,
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
                                   const amrex::Geometry &geom,
                                   const int &spatial_order);


// Compute Sij
AMREX_GPU_DEVICE
Real ComputeStrainRate(const int &i, const int &j, const int &k,
                        const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                        const enum NextOrPrev &nextOrPrev,
                        const enum MomentumEqn &momentumEqn,
                        const enum DiffusionDir &diffDir,
                        const Geometry &geom);

// Compute turbVisc, i.e., 'nut' or Kij at cell centres. Initialize it to 1.0
//TODO: Check if it is better to  fill all Kij and use it or compute on the fly
AMREX_GPU_DEVICE
Real ComputeTurbulentViscosity(const int &i, const int &j, const int &k,
                           const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                           const enum NextOrPrev &nextOrPrev,
                           const enum MomentumEqn &momentumEqn,
                           const enum DiffusionDir &diffDir,
                           const Geometry &geom,
                           Array4<Real>& nut);

// Given a velocity field and a geometry, compute nut or K(i, j, k) are cell-centers
//AMREX_GPU_DEVICE
void ComputeTurbulentViscosity(MultiFab& xvel, MultiFab& yvel, MultiFab& zvel,
                               MultiFab& eddyViscosity,
                               const Geometry &geom, const SolverChoice& solverChoice);

// Compute K (i-1/2, j+1/2, k) etc given K(i, j, k) or nut (i, j, k) is known
AMREX_GPU_DEVICE
Real InterpolateTurbulentViscosity(const int &i, const int &j, const int &k,
                            const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                            const enum NextOrPrev &nextOrPrev,
                            const enum MomentumEqn &momentumEqn,
                            const enum DiffusionDir &diffDir,
                            const Geometry &geom,
                            const Array4<Real>& nut);

// Compute tau_ij (m + 1/2), tau_ij (m - 1/2) where m = {i, j, k} for DNS or Smagorinsky
AMREX_GPU_DEVICE
Real ComputeStressTerm (const int &i, const int &j, const int &k,
                        const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                        const enum NextOrPrev &nextOrPrev,
                        const enum MomentumEqn &momentumEqn,
                        const enum DiffusionDir &diffDir,
                        const Geometry &geom,
                        const Array4<Real>& nut,
                        const SolverChoice &solverChoice);

AMREX_GPU_DEVICE
Real DiffusionContributionForMom(const int &i, const int &j, const int &k,
                                 const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                                 const enum MomentumEqn &momentumEqn,
                                 const amrex::Geometry &geom,
                                 const Array4<Real>& nut,
                                 const SolverChoice &solverChoice);

AMREX_GPU_DEVICE
Real ComputeDiffusionTermForState(const int &i, const int &j, const int &k,
                                  const Array4<Real>& cell_data, const int & qty_index,
                                  const enum Coord& coordDir);

AMREX_GPU_DEVICE
Real DiffusionContributionForState(const int &i, const int &j, const int &k,
                                   const Array4<Real>& cell_data, const int & qty_index,
                                   const amrex::Geometry &geom,
                                   const SolverChoice &solverChoice);

int ComputeGhostCells(const int &spatial_order);

#endif
