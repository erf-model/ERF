#ifndef _STENCILS_H_
#define _STENCILS_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_InterpFaceRegister.H>
#include <DataStruct.H>
#include <IndexDefines.H>
#include <ABLMost.H>

void make_metrics (const amrex::Geometry& geom,
                   amrex::MultiFab& z_phys_nd,
                   amrex::MultiFab& z_phys_cc,
                   amrex::MultiFab& detJ_cc);

void MomentumToVelocity (amrex::MultiFab& xvel_out,
                         const amrex::IntVect& xvel_ngrow,
                         amrex::MultiFab& yvel_out,
                         const amrex::IntVect& yvel_ngrow,
                         amrex::MultiFab& zvel_out,
                         const amrex::IntVect& zvel_ngrow,
                         const amrex::MultiFab& cons_in,
                         const amrex::MultiFab& xmom_in,
                         const amrex::MultiFab& ymom_in,
                         const amrex::MultiFab& zmom_in);

void VelocityToMomentum (const amrex::MultiFab& xvel_in,
                         const amrex::IntVect& xvel_ngrow,
                         const amrex::MultiFab& yvel_in,
                         const amrex::IntVect& yvel_ngrow,
                         const amrex::MultiFab& zvel_in,
                         const amrex::IntVect& zvel_ngrow,
                         const amrex::MultiFab& cons_in,
                         amrex::MultiFab& xmom_out,
                         amrex::MultiFab& ymom_out,
                         amrex::MultiFab& zmom_out);

/** Compute advection source for the continuity, energy, and scalar equations */
void AdvectionSrcForRhoAndTheta (const amrex::Box& bx, const amrex::Box& valid_bx,
                                 const amrex::Array4<amrex::Real>& src,
                                 const amrex::Array4<const amrex::Real>& rho_u,    // These are being used
                                 const amrex::Array4<const amrex::Real>& rho_v,    //  to define the fluxes
                                 const amrex::Array4<const amrex::Real>& rho_w,
                                 amrex::Real fac,
                                 const amrex::Array4<      amrex::Real>& avg_xmom, // These are being defined
                                 const amrex::Array4<      amrex::Real>& avg_ymom, //  from the rho fluxes
                                 const amrex::Array4<      amrex::Real>& avg_zmom,
                                 const amrex::Array4<const amrex::Real>& z_t,
                                 const amrex::Array4<const amrex::Real>& cell_prim,
                                 const amrex::Array4<const amrex::Real>& z_nd,
                                 const amrex::Array4<const amrex::Real>& detJ,
                                 const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSize,
                                 const int& spatial_order, const int& use_terrain);

void AdvectionSrcForScalars (const amrex::Box& bx,
                             const int &start_comp, const int &num_comp,
                             const amrex::Array4<const amrex::Real>& rho_u,
                             const amrex::Array4<const amrex::Real>& rho_v,
                             const amrex::Array4<const amrex::Real>& rho_w,
                             const amrex::Array4<const amrex::Real>& cell_prim,
                             const amrex::Array4<amrex::Real>& src,
                             const amrex::Array4<const amrex::Real>& detJ,
                             const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSize,
                             const int& spatial_order, const int& use_terrain,
                             const int& use_deardorff, const int& use_QKE);

void AdvectionSrcForMom (int level, const amrex::Box& bx, const amrex::Box& valid_bx,
                         const amrex::Array4<const int>& mlo_x            , const amrex::Array4<const int>& mlo_y,
                         const amrex::Array4<const int>& mhi_x            , const amrex::Array4<const int>& mhi_y,
                         const amrex::Array4<      amrex::Real>& rho_u_rhs, const amrex::Array4<      amrex::Real>& rho_v_rhs,
                         const amrex::Array4<      amrex::Real>& rho_w_rhs,
                         const amrex::Array4<const amrex::Real>& u        , const amrex::Array4<const amrex::Real>& v,
                         const amrex::Array4<const amrex::Real>& w        ,
                         const amrex::Array4<const amrex::Real>& rho_u    , const amrex::Array4<const amrex::Real>& rho_v,
                         const amrex::Array4<const amrex::Real>& rho_w    , const amrex::Array4<const amrex::Real>& z_t,
                         const amrex::Array4<const amrex::Real>& z_nd     , const amrex::Array4<const amrex::Real>& detJ,
                         const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                         const int spatial_order, const int use_terrain, const int domhi_z);

void DiffusionSrcForMom_N (int level, const amrex::Box& bx, const amrex::Box& valid_bx, const amrex::Box& domain,
                           const amrex::Array4<const int>& mlo_x            , const amrex::Array4<const int>& mlo_y,
                           const amrex::Array4<const int>& mhi_x            , const amrex::Array4<const int>& mhi_y,
                           const amrex::Array4<      amrex::Real>& rho_u_rhs, const amrex::Array4<      amrex::Real>& rho_v_rhs,
                           const amrex::Array4<      amrex::Real>& rho_w_rhs,
                           const amrex::Array4<const amrex::Real>& u        , const amrex::Array4<const amrex::Real>& v,
                           const amrex::Array4<const amrex::Real>& w        , const amrex::Array4<const amrex::Real>& K_turb,
                           const amrex::Array4<const amrex::Real>& cell_data, const amrex::Array4<const amrex::Real>& er_arr,
                           const SolverChoice& solverChoice                 , const amrex::BCRec* bc_ptr,
                           const amrex::Array4<const amrex::Real>& z_nd     , const amrex::Array4<const amrex::Real>& detJ,
                           const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv);

void DiffusionSrcForMom_T (int level, const amrex::Box& bx, const amrex::Box& valid_bx, const amrex::Box& domain,
                           const amrex::Array4<const int>& mlo_x            , const amrex::Array4<const int>& mlo_y,
                           const amrex::Array4<const int>& mhi_x            , const amrex::Array4<const int>& mhi_y,
                           const amrex::Array4<      amrex::Real>& rho_u_rhs, const amrex::Array4<      amrex::Real>& rho_v_rhs,
                           const amrex::Array4<      amrex::Real>& rho_w_rhs,
                           const amrex::Array4<const amrex::Real>& u        , const amrex::Array4<const amrex::Real>& v,
                           const amrex::Array4<const amrex::Real>& w        , const amrex::Array4<const amrex::Real>& K_turb,
                           const amrex::Array4<const amrex::Real>& cell_data, const amrex::Array4<const amrex::Real>& er_arr,
                           const SolverChoice& solverChoice                 , const amrex::BCRec* bc_ptr,
                           const amrex::Array4<const amrex::Real>& z_nd     , const amrex::Array4<const amrex::Real>& detJ,
                           const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv);

void DiffusionSrcForState (const amrex::Box& bx, const amrex::Box& domain, int n_start, int n_end,
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
                           const amrex::Array4<const amrex::Real>& K_turb,
                           const SolverChoice &solverChoice,
                           const amrex::Real& theta_mean,
                           const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> grav_gpu,
                           const amrex::BCRec* bc_ptr);
#endif
