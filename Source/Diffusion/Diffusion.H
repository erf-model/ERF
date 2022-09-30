#ifndef _DIFFUSION_H_
#define _DIFFUSION_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_InterpFaceRegister.H>
#include <DataStruct.H>
#include <IndexDefines.H>
#include <ABLMost.H>

void DiffusionSrcForMom_N (const amrex::Box& bxx, const amrex::Box& bxy, const amrex::Box& bxz,
                           const amrex::Array4<      amrex::Real>& rho_u_rhs, const amrex::Array4<      amrex::Real>& rho_v_rhs, const amrex::Array4<      amrex::Real>& rho_w_rhs,
                           const amrex::Array4<const amrex::Real>& tau11    , const amrex::Array4<const amrex::Real>& tau22    , const amrex::Array4<const amrex::Real>& tau33    ,
                           const amrex::Array4<const amrex::Real>& tau12    , const amrex::Array4<const amrex::Real>& tau13    , const amrex::Array4<const amrex::Real>& tau23    ,
                           const amrex::Array4<const amrex::Real>& cell_data, const SolverChoice& solverChoice,  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv);

void DiffusionSrcForMom_T (const amrex::Box& bxx, const amrex::Box& bxy, const amrex::Box& bxz,
                           const amrex::Array4<      amrex::Real>& rho_u_rhs, const amrex::Array4<      amrex::Real>& rho_v_rhs, const amrex::Array4<      amrex::Real>& rho_w_rhs,
                           const amrex::Array4<const amrex::Real>& tau11    , const amrex::Array4<const amrex::Real>& tau22    , const amrex::Array4<const amrex::Real>& tau33    ,
                           const amrex::Array4<const amrex::Real>& tau12    , const amrex::Array4<const amrex::Real>& tau13,
                           const amrex::Array4<const amrex::Real>& tau21    , const amrex::Array4<const amrex::Real>& tau23,
                           const amrex::Array4<const amrex::Real>& tau31    , const amrex::Array4<const amrex::Real>& tau32,
                           const amrex::Array4<const amrex::Real>& cell_data, const amrex::Array4<const amrex::Real>& detJ,
                           const SolverChoice& solverChoice                 , const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv);

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

void DiffusionSrcForState_N (const amrex::Box& bx, const amrex::Box& domain, int n_start, int n_end,
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

void DiffusionSrcForState_T (const amrex::Box& bx, const amrex::Box& domain, int n_start, int n_end,
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
                             const amrex::Array4<const amrex::Real>& z_nd,
                             const amrex::Array4<const amrex::Real>& detJ,
                             const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxInv,
                             const amrex::Array4<const amrex::Real>& K_turb,
                             const SolverChoice &solverChoice,
                             const amrex::Real& theta_mean,
                             const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> grav_gpu,
                             const amrex::BCRec* bc_ptr);

void ComputeStressConsVisc_N(amrex::Box& bxcc, amrex::Box& tbxxy, amrex::Box& tbxxz, amrex::Box& tbxyz, amrex::Real mu_eff,
                             amrex::Array4<amrex::Real>& tau11, amrex::Array4<amrex::Real>& tau22, amrex::Array4<amrex::Real>& tau33,
                             amrex::Array4<amrex::Real>& tau12, amrex::Array4<amrex::Real>& tau13, amrex::Array4<amrex::Real>& tau23,
                             const amrex::Array4<const amrex::Real>& er_arr,
                             const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv);

void ComputeStressConsVisc_T(amrex::Box& bxcc, amrex::Box& tbxxy, amrex::Box& tbxxz, amrex::Box& tbxyz, amrex::Real mu_eff,
                             amrex::Array4<amrex::Real>& tau11, amrex::Array4<amrex::Real>& tau22, amrex::Array4<amrex::Real>& tau33,
                             amrex::Array4<amrex::Real>& tau12, amrex::Array4<amrex::Real>& tau13,
                             amrex::Array4<amrex::Real>& tau21, amrex::Array4<amrex::Real>& tau23,
                             amrex::Array4<amrex::Real>& tau31, amrex::Array4<amrex::Real>& tau32,
                             const amrex::Array4<const amrex::Real>& er_arr,
                             const amrex::Array4<const amrex::Real>& z_nd  ,
                             const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv);

void ComputeStressVarVisc_N(amrex::Box& bxcc, amrex::Box& tbxxy, amrex::Box& tbxxz, amrex::Box& tbxyz, amrex::Real mu_eff,
                            const amrex::Array4<const amrex::Real>& K_turb,
                            amrex::Array4<amrex::Real>& tau11, amrex::Array4<amrex::Real>& tau22, amrex::Array4<amrex::Real>& tau33,
                            amrex::Array4<amrex::Real>& tau12, amrex::Array4<amrex::Real>& tau13, amrex::Array4<amrex::Real>& tau23,
                            const amrex::Array4<const amrex::Real>& er_arr,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv);

void ComputeStressVarVisc_T(amrex::Box& bxcc, amrex::Box& tbxxy, amrex::Box& tbxxz, amrex::Box& tbxyz, amrex::Real mu_eff,
                            const amrex::Array4<const amrex::Real>& K_turb,
                            amrex::Array4<amrex::Real>& tau11, amrex::Array4<amrex::Real>& tau22, amrex::Array4<amrex::Real>& tau33,
                            amrex::Array4<amrex::Real>& tau12, amrex::Array4<amrex::Real>& tau13,
                            amrex::Array4<amrex::Real>& tau21, amrex::Array4<amrex::Real>& tau23,
                            amrex::Array4<amrex::Real>& tau31, amrex::Array4<amrex::Real>& tau32,
                            const amrex::Array4<const amrex::Real>& er_arr,
                            const amrex::Array4<const amrex::Real>& z_nd  ,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv);

void ComputeStrain_N(amrex::Box& bxcc, amrex::Box& tbxxy, amrex::Box& tbxxz, amrex::Box& tbxyz,
                     const amrex::Array4<const amrex::Real>& u,
                     const amrex::Array4<const amrex::Real>& v,
                     const amrex::Array4<const amrex::Real>& w,
                     amrex::Array4<amrex::Real>& tau11, amrex::Array4<amrex::Real>& tau22, amrex::Array4<amrex::Real>& tau33,
                     amrex::Array4<amrex::Real>& tau12, amrex::Array4<amrex::Real>& tau13, amrex::Array4<amrex::Real>& tau23,
                     const amrex::BCRec* bc_ptr, const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv);

void ComputeStrain_T(amrex::Box& bxcc, amrex::Box& tbxxy, amrex::Box& tbxxz, amrex::Box& tbxyz,
                     const amrex::Array4<const amrex::Real>& u,
                     const amrex::Array4<const amrex::Real>& v,
                     const amrex::Array4<const amrex::Real>& w,
                     amrex::Array4<amrex::Real>& tau11, amrex::Array4<amrex::Real>& tau22, amrex::Array4<amrex::Real>& tau33,
                     amrex::Array4<amrex::Real>& tau12, amrex::Array4<amrex::Real>& tau13,
                     amrex::Array4<amrex::Real>& tau21, amrex::Array4<amrex::Real>& tau23,
                     amrex::Array4<amrex::Real>& tau31, amrex::Array4<amrex::Real>& tau32,
                     const amrex::Array4<const amrex::Real>& z_nd  ,
                     const amrex::BCRec* bc_ptr, const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv);

void ComputeTurbVisc_SMAG(amrex::Box& bxcc,
                          const amrex::Array4<amrex::Real>& K_turb,
                          const amrex::Array4<const amrex::Real>& cell_data,
                          const amrex::Array4<const amrex::Real>& tau11,
                          const amrex::Array4<const amrex::Real>& tau22,
                          const amrex::Array4<const amrex::Real>& tau33,
                          const amrex::Array4<const amrex::Real>& tau12,
                          const amrex::Array4<const amrex::Real>& tau13,
                          const amrex::Array4<const amrex::Real>& tau23,
                          const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxInv,
                          const SolverChoice& solverChoice);
#endif
