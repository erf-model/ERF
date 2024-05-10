#ifndef _DIFFUSION_H_
#define _DIFFUSION_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <DataStruct.H>
#include <IndexDefines.H>
#include <ABLMost.H>

void DiffusionSrcForMom_N (const amrex::Box& bxx, const amrex::Box& bxy, const amrex::Box& bxz,
                           const amrex::Array4<      amrex::Real>& rho_u_rhs,
                           const amrex::Array4<      amrex::Real>& rho_v_rhs,
                           const amrex::Array4<      amrex::Real>& rho_w_rhs,
                           const amrex::Array4<const amrex::Real>& tau11    ,
                           const amrex::Array4<const amrex::Real>& tau22    ,
                           const amrex::Array4<const amrex::Real>& tau33    ,
                           const amrex::Array4<const amrex::Real>& tau12    ,
                           const amrex::Array4<const amrex::Real>& tau13    ,
                           const amrex::Array4<const amrex::Real>& tau23    ,
                           const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxInv,
                           const amrex::Array4<const amrex::Real>& mf_m      ,
                           const amrex::Array4<const amrex::Real>& mf_u      ,
                           const amrex::Array4<const amrex::Real>& mf_v      );

void DiffusionSrcForMom_T (const amrex::Box& bxx, const amrex::Box& bxy, const amrex::Box& bxz,
                           const amrex::Array4<      amrex::Real>& rho_u_rhs,
                           const amrex::Array4<      amrex::Real>& rho_v_rhs,
                           const amrex::Array4<      amrex::Real>& rho_w_rhs,
                           const amrex::Array4<const amrex::Real>& tau11    ,
                           const amrex::Array4<const amrex::Real>& tau22    , const amrex::Array4<const amrex::Real>& tau33    ,
                           const amrex::Array4<const amrex::Real>& tau12    , const amrex::Array4<const amrex::Real>& tau13,
                           const amrex::Array4<const amrex::Real>& tau21    , const amrex::Array4<const amrex::Real>& tau23,
                           const amrex::Array4<const amrex::Real>& tau31    , const amrex::Array4<const amrex::Real>& tau32,
                           const amrex::Array4<const amrex::Real>& detJ,
                           const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxInv,
                           const amrex::Array4<const amrex::Real>& mf_m      ,
                           const amrex::Array4<const amrex::Real>& mf_u      ,
                           const amrex::Array4<const amrex::Real>& mf_v      );



void DiffusionSrcForState_N (const amrex::Box& bx, const amrex::Box& domain,
                             int start_comp, int num_comp,
                             const bool& exp_most,
                             const amrex::Array4<const amrex::Real>& u,
                             const amrex::Array4<const amrex::Real>& v,
                             const amrex::Array4<const amrex::Real>& cell_data,
                             const amrex::Array4<const amrex::Real>& cell_prim,
                             const amrex::Array4<amrex::Real>& cell_rhs,
                             const amrex::Array4<amrex::Real>& xflux,
                             const amrex::Array4<amrex::Real>& yflux,
                             const amrex::Array4<amrex::Real>& zflux,
                             const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                             const amrex::Array4<const amrex::Real>& SmnSmn_a,
                             const amrex::Array4<const amrex::Real>& mf_m,
                             const amrex::Array4<const amrex::Real>& mf_u,
                             const amrex::Array4<const amrex::Real>& mf_v ,
                                   amrex::Array4<      amrex::Real>& hfx_z,
                                   amrex::Array4<      amrex::Real>& diss,
                             const amrex::Array4<const amrex::Real>& mu_turb,
                             const DiffChoice& diffChoice,
                             const TurbChoice& turbChoice,
                             const amrex::Array4<const amrex::Real>& tm_arr,
                             const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> grav_gpu,
                             const amrex::BCRec* bc_ptr,
                             const bool use_most);

void DiffusionSrcForState_T (const amrex::Box& bx, const amrex::Box& domain,
                             int start_comp, int num_comp,
                             const bool& exp_most,
                             const amrex::Array4<const amrex::Real>& u,
                             const amrex::Array4<const amrex::Real>& v,
                             const amrex::Array4<const amrex::Real>& cell_data,
                             const amrex::Array4<const amrex::Real>& cell_prim,
                             const amrex::Array4<amrex::Real>& cell_rhs,
                             const amrex::Array4<amrex::Real>& xflux,
                             const amrex::Array4<amrex::Real>& yflux,
                             const amrex::Array4<amrex::Real>& zflux,
                             const amrex::Array4<const amrex::Real>& z_nd,
                             const amrex::Array4<const amrex::Real>& ax,
                             const amrex::Array4<const amrex::Real>& ay,
                             const amrex::Array4<const amrex::Real>& az,
                             const amrex::Array4<const amrex::Real>& detJ,
                             const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxInv,
                             const amrex::Array4<const amrex::Real>& SmnSmn_a,
                             const amrex::Array4<const amrex::Real>& mf_m,
                             const amrex::Array4<const amrex::Real>& mf_u,
                             const amrex::Array4<const amrex::Real>& mf_v ,
                                   amrex::Array4<      amrex::Real>& hfx_z,
                                   amrex::Array4<      amrex::Real>& diss,
                             const amrex::Array4<const amrex::Real>& mu_turb,
                             const DiffChoice& diffChoice,
                             const TurbChoice& turbChoice,
                             const amrex::Array4<const amrex::Real>& tm_arr,
                             const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> grav_gpu,
                             const amrex::BCRec* bc_ptr,
                             const bool use_most);



void ComputeStressConsVisc_N (amrex::Box bxcc, amrex::Box tbxxy, amrex::Box tbxxz, amrex::Box tbxyz, amrex::Real mu_eff,
                              const amrex::Array4<const amrex::Real>& cell_data,
                              amrex::Array4<amrex::Real>& tau11, amrex::Array4<amrex::Real>& tau22, amrex::Array4<amrex::Real>& tau33,
                              amrex::Array4<amrex::Real>& tau12, amrex::Array4<amrex::Real>& tau13, amrex::Array4<amrex::Real>& tau23,
                              const amrex::Array4<const amrex::Real>& er_arr);

void ComputeStressConsVisc_T (amrex::Box bxcc, amrex::Box tbxxy, amrex::Box tbxxz, amrex::Box tbxyz, amrex::Real mu_eff,
                              const amrex::Array4<const amrex::Real>& cell_data,
                              amrex::Array4<amrex::Real>& tau11, amrex::Array4<amrex::Real>& tau22, amrex::Array4<amrex::Real>& tau33,
                              amrex::Array4<amrex::Real>& tau12, amrex::Array4<amrex::Real>& tau13,
                              amrex::Array4<amrex::Real>& tau21, amrex::Array4<amrex::Real>& tau23,
                              amrex::Array4<amrex::Real>& tau31, amrex::Array4<amrex::Real>& tau32,
                              const amrex::Array4<const amrex::Real>& er_arr,
                              const amrex::Array4<const amrex::Real>& z_nd,
                              const amrex::Array4<const amrex::Real>& detJ,
                              const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxInv);



void ComputeStressVarVisc_N (amrex::Box bxcc, amrex::Box tbxxy, amrex::Box tbxxz, amrex::Box tbxyz, amrex::Real mu_eff,
                             const amrex::Array4<const amrex::Real>& mu_turb,
                             const amrex::Array4<const amrex::Real>& cell_data,
                             amrex::Array4<amrex::Real>& tau11, amrex::Array4<amrex::Real>& tau22, amrex::Array4<amrex::Real>& tau33,
                             amrex::Array4<amrex::Real>& tau12, amrex::Array4<amrex::Real>& tau13, amrex::Array4<amrex::Real>& tau23,
                             const amrex::Array4<const amrex::Real>& er_arr);

void ComputeStressVarVisc_T (amrex::Box bxcc, amrex::Box tbxxy, amrex::Box tbxxz, amrex::Box tbxyz, amrex::Real mu_eff,
                             const amrex::Array4<const amrex::Real>& mu_turb,
                             const amrex::Array4<const amrex::Real>& cell_data,
                             amrex::Array4<amrex::Real>& tau11, amrex::Array4<amrex::Real>& tau22, amrex::Array4<amrex::Real>& tau33,
                             amrex::Array4<amrex::Real>& tau12, amrex::Array4<amrex::Real>& tau13,
                             amrex::Array4<amrex::Real>& tau21, amrex::Array4<amrex::Real>& tau23,
                             amrex::Array4<amrex::Real>& tau31, amrex::Array4<amrex::Real>& tau32,
                             const amrex::Array4<const amrex::Real>& er_arr,
                             const amrex::Array4<const amrex::Real>& z_nd,
                             const amrex::Array4<const amrex::Real>& detJ,
                             const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxInv);



void ComputeStrain_N (amrex::Box bxcc, amrex::Box tbxxy, amrex::Box tbxxz, amrex::Box tbxyz, amrex::Box domain,
                     const amrex::Array4<const amrex::Real>& u,
                     const amrex::Array4<const amrex::Real>& v,
                     const amrex::Array4<const amrex::Real>& w,
                     amrex::Array4<amrex::Real>& tau11, amrex::Array4<amrex::Real>& tau22, amrex::Array4<amrex::Real>& tau33,
                     amrex::Array4<amrex::Real>& tau12, amrex::Array4<amrex::Real>& tau13, amrex::Array4<amrex::Real>& tau23,
                     const amrex::BCRec* bc_ptr, const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxInv,
                     const amrex::Array4<const amrex::Real>& mf_m, const amrex::Array4<const amrex::Real>& mf_u, const amrex::Array4<const amrex::Real>& mf_v);

void ComputeStrain_T (amrex::Box bxcc, amrex::Box tbxxy, amrex::Box tbxxz, amrex::Box tbxyz, amrex::Box domain,
                     const amrex::Array4<const amrex::Real>& u,
                     const amrex::Array4<const amrex::Real>& v,
                     const amrex::Array4<const amrex::Real>& w,
                     amrex::Array4<amrex::Real>& tau11, amrex::Array4<amrex::Real>& tau22, amrex::Array4<amrex::Real>& tau33,
                     amrex::Array4<amrex::Real>& tau12, amrex::Array4<amrex::Real>& tau13,
                     amrex::Array4<amrex::Real>& tau21, amrex::Array4<amrex::Real>& tau23,
                     amrex::Array4<amrex::Real>& tau31, amrex::Array4<amrex::Real>& tau32,
                     const amrex::Array4<const amrex::Real>& z_nd,
                     const amrex::Array4<const amrex::Real>& detJ,
                     const amrex::BCRec* bc_ptr, const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxInv,
                     const amrex::Array4<const amrex::Real>& mf_m, const amrex::Array4<const amrex::Real>& mf_u, const amrex::Array4<const amrex::Real>& mf_v);
#endif
