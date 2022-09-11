#ifndef _DIFFUSION_H_
#define _DIFFUSION_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_InterpFaceRegister.H>
#include <DataStruct.H>
#include <IndexDefines.H>
#include <ABLMost.H>

void DiffusionSrcForMom_N (int level, const amrex::Box& bx, const amrex::Box& valid_bx, const amrex::Box& domain,
                           const amrex::Array4<const int>& mlo_x            , const amrex::Array4<const int>& mlo_y,
                           const amrex::Array4<const int>& mhi_x            , const amrex::Array4<const int>& mhi_y,
                           const amrex::Array4<      amrex::Real>& rho_u_rhs, const amrex::Array4<      amrex::Real>& rho_v_rhs,
                           const amrex::Array4<      amrex::Real>& rho_w_rhs,
                           const amrex::Array4<const amrex::Real>& u        , const amrex::Array4<const amrex::Real>& v,
                           const amrex::Array4<const amrex::Real>& w        , const amrex::Array4<const amrex::Real>& K_turb,
                           const amrex::Array4<const amrex::Real>& cell_data, const amrex::Array4<const amrex::Real>& er_arr,
                           const SolverChoice& solverChoice                 , const amrex::BCRec* bc_ptr,
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
