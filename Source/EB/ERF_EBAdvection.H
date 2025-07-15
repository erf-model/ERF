#ifndef ERF_EBADVECTION_H_
#define ERF_EBADVECTION_H_

#include "AMReX.H"
#include "AMReX_MultiFab.H"
#include "AMReX_BCRec.H"
#include "AMReX_YAFluxRegister.H"

#include "ERF_DataStruct.H"
#include "ERF_IndexDefines.H"
#include "ERF_SurfaceLayer.H"

/** Compute advection tendency for density and potential temperature */
void EBAdvectionSrcForRho (const amrex::Box& bx,
                           const amrex::Array4<amrex::Real>& src,
                           const amrex::Array4<const amrex::Real>& rho_u,    // These are being used
                           const amrex::Array4<const amrex::Real>& rho_v,    //  to define the fluxes
                           const amrex::Array4<const amrex::Real>& omega,
                           const amrex::Array4<      amrex::Real>& avg_xmom, // These are being defined
                           const amrex::Array4<      amrex::Real>& avg_ymom, //  from the rho fluxes
                           const amrex::Array4<      amrex::Real>& avg_zmom,
                           const amrex::Array4<const int>& mask_arr,
                           const amrex::Array4<const amrex::EBCellFlag>& cfg_arr,
                           const amrex::Array4<const amrex::Real>& ax_arr,
                           const amrex::Array4<const amrex::Real>& ay_arr,
                           const amrex::Array4<const amrex::Real>& az_arr,
                           const amrex::Array4<const amrex::Real>& fcx_arr,
                           const amrex::Array4<const amrex::Real>& fcy_arr,
                           const amrex::Array4<const amrex::Real>& fcz_arr,
                           const amrex::Array4<const amrex::Real>& detJ,
                           const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                           const amrex::Array4<const amrex::Real>& mf_mx,
                           const amrex::Array4<const amrex::Real>& mf_my,
                           const amrex::Array4<const amrex::Real>& mf_uy,
                           const amrex::Array4<const amrex::Real>& mf_vx,
                           const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_arr,
                           const bool fixed_rho,
                           bool already_on_centroids);

/** Compute advection tendency for all scalars other than density and potential temperature */
void EBAdvectionSrcForScalars (const amrex::Box& bx,
                               const int icomp, const int ncomp,
                               const amrex::Array4<const amrex::Real>& avg_xmom,
                               const amrex::Array4<const amrex::Real>& avg_ymom,
                               const amrex::Array4<const amrex::Real>& avg_zmom,
                               const amrex::Array4<const amrex::Real>& cell_prim,
                               const amrex::Array4<amrex::Real>& src,
                               const amrex::Array4<const int>& mask_arr,
                               const amrex::Array4<const amrex::EBCellFlag>& cfg_arr,
                               const amrex::Array4<const amrex::Real>& ax_arr,
                               const amrex::Array4<const amrex::Real>& ay_arr,
                               const amrex::Array4<const amrex::Real>& az_arr,
                               const amrex::Array4<const amrex::Real>& fcx_arr,
                               const amrex::Array4<const amrex::Real>& fcy_arr,
                               const amrex::Array4<const amrex::Real>& fcz_arr,
                               const amrex::Array4<const amrex::Real>& vf_arr,
                               const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                               const amrex::Array4<const amrex::Real>& mf_mx,
                               const amrex::Array4<const amrex::Real>& mf_my,
                               const AdvType horiz_adv_type, const AdvType vert_adv_type,
                               const amrex::Real horiz_upw_frac, const amrex::Real vert_upw_frac,
                               const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_arr,
                               const amrex::Box& domain,
                               const amrex::BCRec* bc_ptr_h,
                               bool already_on_centroids);

#endif
