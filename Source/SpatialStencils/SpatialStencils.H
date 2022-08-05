#ifndef _STENCILS_H_
#define _STENCILS_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_InterpFaceRegister.H>
#include <DataStruct.H>
#include <IndexDefines.H>
#include <ABLMost.H>
#include <PBLModels.H>

void MomentumToVelocity(amrex::MultiFab& xvel_out,
                        const amrex::IntVect& xvel_ngrow,
                        amrex::MultiFab& yvel_out,
                        const amrex::IntVect& yvel_ngrow,
                        amrex::MultiFab& zvel_out,
                        const amrex::IntVect& zvel_ngrow,
                        const amrex::MultiFab& cons_in,
                        const amrex::MultiFab& xmom_in,
                        const amrex::MultiFab& ymom_in,
                        const amrex::MultiFab& zmom_in);

void VelocityToMomentum(const amrex::MultiFab& xvel_in,
                        const amrex::IntVect& xvel_ngrow,
                        const amrex::MultiFab& yvel_in,
                        const amrex::IntVect& yvel_ngrow,
                        const amrex::MultiFab& zvel_in,
                        const amrex::IntVect& zvel_ngrow,
                        const amrex::MultiFab& cons_in,
                        amrex::MultiFab& xmom_out,
                        amrex::MultiFab& ymom_out,
                        amrex::MultiFab& zmom_out);

/** Meant for {x, y, z}- momentum equations */
AMREX_GPU_DEVICE
amrex::Real AdvectionSrcForXMom(
                       const int &i, const int &j, const int &k,
                       const amrex::Array4<const amrex::Real>& rho_u,
                       const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& rho_w,
                       const amrex::Array4<const amrex::Real>& u,
                       const amrex::Array4<const amrex::Real>& z_nd,
                       const amrex::Array4<const amrex::Real>& detJ,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSize,
                       const int& spatial_order, const int& use_terrain);
AMREX_GPU_DEVICE
amrex::Real AdvectionSrcForYMom(
                       const int &i, const int &j, const int &k,
                       const amrex::Array4<const amrex::Real>& rho_u,
                       const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& rho_w,
                       const amrex::Array4<const amrex::Real>& v,
                       const amrex::Array4<const amrex::Real>& z_nd,
                       const amrex::Array4<const amrex::Real>& detJ,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSize,
                       const int& spatial_order, const int& use_terrain);
AMREX_GPU_DEVICE
amrex::Real AdvectionSrcForZMom(
                       const int &i, const int &j, const int &k,
                       const amrex::Array4<const amrex::Real>& rho_u,
                       const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& rho_w,
                       const amrex::Array4<const amrex::Real>& w,
                       const amrex::Array4<const amrex::Real>& z_nd,
                       const amrex::Array4<const amrex::Real>& detJ,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSize,
                       const int& spatial_order, const int& use_terrain,
                       const int& domhi_z);

/** Compute advection source for the continuity, energy, and scalar equations */
void AdvectionSrcForState(
                       const amrex::Box& bx,
                       const int &start_comp, const int &num_comp,
                       const amrex::Array4<const amrex::Real>& rho_u,
                       const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& rho_w,
                       const amrex::Array4<const amrex::Real>& cell_prim,
                       const amrex::Array4<amrex::Real>& src,
                       const amrex::Array4<amrex::Real>& flux_x,
                       const amrex::Array4<amrex::Real>& flux_y,
                       const amrex::Array4<amrex::Real>& flux_z,
                       const amrex::Array4<const amrex::Real>& z_nd,
                       const amrex::Array4<const amrex::Real>& detJ,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSize,
                       const int& spatial_order, const int& use_terrain,
                       const int& use_deardorff, const int& use_QKE);

/** Compute K (i-1/2, j+1/2, k) etc given K(i, j, k) or nut (i, j, k) is known */
AMREX_GPU_DEVICE
amrex::Real InterpolateTurbulentViscosity(const int &i, const int &j, const int &k,
                            const enum MomentumEqn &momentumEqn,
                            const enum DiffusionDir &diffDir,
                            const amrex::Array4<amrex::Real>& Ksmag);

/** Compute tau_ij (m + 1/2), tau_ij (m - 1/2) where m = {i, j, k} for DNS or Smagorinsky */
AMREX_GPU_DEVICE
amrex::Real ComputeStressTerm (
                       const int &i, const int &j, const int &k,
                       const amrex::Array4<const amrex::Real>& u,
                       const amrex::Array4<const amrex::Real>& v,
                       const amrex::Array4<const amrex::Real>& w,
                       const enum MomentumEqn &momentumEqn,
                       const enum DiffusionDir &diffDir,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSize,
                       const amrex::Array4<amrex::Real>& Ksmag,
                       const SolverChoice &solverChoice,
                       const amrex::Box& domain, const amrex::BCRec* bc_ptr,
                       const amrex::Real er_at_face);

AMREX_GPU_DEVICE
amrex::Real DiffusionSrcForMomWithTerrain(
                       const int &i, const int &j, const int &k,
                       const amrex::Array4<const amrex::Real>& u,
                       const amrex::Array4<const amrex::Real>& v,
                       const amrex::Array4<const amrex::Real>& w,
                       const amrex::Array4<const amrex::Real>& cons,
                       const enum MomentumEqn &momentumEqn,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSize,
                       const amrex::Array4<amrex::Real>& Ksmag,
                       const SolverChoice &solverChoice,
                       const amrex::Array4<const amrex::Real>& z_nd,
                       const amrex::Array4<const amrex::Real>& detJ,
                       const amrex::Box& domain, const amrex::BCRec* bc_ptr);

AMREX_GPU_DEVICE
amrex::Real DiffusionSrcForMom(
                       const int &i, const int &j, const int &k,
                       const amrex::Array4<const amrex::Real>& u,
                       const amrex::Array4<const amrex::Real>& v,
                       const amrex::Array4<const amrex::Real>& w,
                       const amrex::Array4<const amrex::Real>& cons,
                       const enum MomentumEqn &momentumEqn,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSize,
                       const amrex::Array4<amrex::Real>& Ksmag,
                       const SolverChoice &solverChoice,
                       const amrex::Box& domain, const amrex::BCRec* bc_ptr,
                       const amrex::Array4<const amrex::Real>& er_arr);

AMREX_GPU_DEVICE
amrex::Real ComputeDiffusionFluxForState(
                       const int &i, const int &j, const int &k,
                       const amrex::Array4<const amrex::Real>& cell_data,
                       const amrex::Array4<const amrex::Real>& cell_prim,
                       const int & prim_index,
                       const amrex::Real invCellWidth,
                       const amrex::Array4<amrex::Real>& Ksmag,
                       const SolverChoice &solverChoice,
                       const enum Coord& coordDir);

AMREX_GPU_DEVICE
void DiffusionSrcForState(const amrex::Box& bx, const amrex::Box& domain, int qty_index,
                          const amrex::Array4<const amrex::Real>& u,
                          const amrex::Array4<const amrex::Real>& v,
                          const amrex::Array4<const amrex::Real>& w,
                          const amrex::Array4<const amrex::Real>& cell_data,
                          const amrex::Array4<const amrex::Real>& cell_prim,
                          const amrex::Array4<const amrex::Real>& source_fab,
                          const amrex::Array4<amrex::Real>& cell_rhs,
                          const amrex::Array4<amrex::Real>& xflux,
                          const amrex::Array4<amrex::Real>& yflux,
                          const amrex::Array4<amrex::Real>& zflux,
                          const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                          const amrex::Array4<amrex::Real>& K_turb,
                          const SolverChoice &solverChoice,
                          std::unique_ptr<ABLMost>& most,
                          const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> grav_gpu,
                          const amrex::BCRec* bc_ptr,
                          const amrex::Real*  dptr_rayleigh_tau,
                          const amrex::Real*  dptr_rayleigh_thetabar);
#endif
