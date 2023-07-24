#ifndef _ADVECTION_H_
#define _ADVECTION_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_InterpFaceRegister.H>
#include <DataStruct.H>
#include <IndexDefines.H>
#include <ABLMost.H>


/** Compute advection tendency for density and potential temperature */
void AdvectionSrcForRhoAndTheta (const amrex::Box& bx, const amrex::Box& valid_bx,
                                 const amrex::Array4<amrex::Real>& src,
                                 const amrex::Array4<const amrex::Real>& rho_u,    // These are being used
                                 const amrex::Array4<const amrex::Real>& rho_v,    //  to define the fluxes
                                 const amrex::Array4<const amrex::Real>& omega,
                                 amrex::Real fac,
                                 const amrex::Array4<      amrex::Real>& avg_xmom, // These are being defined
                                 const amrex::Array4<      amrex::Real>& avg_ymom, //  from the rho fluxes
                                 const amrex::Array4<      amrex::Real>& avg_zmom,
                                 const amrex::Array4<const amrex::Real>& cell_prim,
                                 const amrex::Array4<const amrex::Real>& z_nd,
                                 const amrex::Array4<const amrex::Real>& detJ,
                                 const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSize,
                                 const amrex::Array4<const amrex::Real>& mf_m,
                                 const amrex::Array4<const amrex::Real>& mf_u,
                                 const amrex::Array4<const amrex::Real>& mf_v,
                                 const int horiz_adv_type, const int vert_adv_type,
                                 const int use_terrain);

/** Compute advection tendency for all scalars other than density and potential temperature */
void AdvectionSrcForScalars (const amrex::Box& bx,
                             const int icomp, const int ncomp,
                             const amrex::Array4<const amrex::Real>& avg_xmom,
                             const amrex::Array4<const amrex::Real>& avg_ymom,
                             const amrex::Array4<const amrex::Real>& avg_zmom,
                             const amrex::Array4<const amrex::Real>& cell_prim,
                             const amrex::Array4<amrex::Real>& src,
                             const amrex::Array4<const amrex::Real>& detJ,
                             const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSize,
                             const amrex::Array4<const amrex::Real>& mf_m,
                             const int horiz_adv_type, const int vert_adv_type,
                             const int use_terrain);

/** Compute advection tendencies for all components of momentum */
void AdvectionSrcForMom (const amrex::Box& bxx, const amrex::Box& bxy, const amrex::Box& bxz,
                         const amrex::Array4<      amrex::Real>& rho_u_rhs, const amrex::Array4<      amrex::Real>& rho_v_rhs,
                         const amrex::Array4<      amrex::Real>& rho_w_rhs,
                         const amrex::Array4<const amrex::Real>& u        , const amrex::Array4<const amrex::Real>& v,
                         const amrex::Array4<const amrex::Real>& w        ,
                         const amrex::Array4<const amrex::Real>& rho_u    , const amrex::Array4<const amrex::Real>& rho_v,
                         const amrex::Array4<const amrex::Real>& Omega    ,
                         const amrex::Array4<const amrex::Real>& z_nd     , const amrex::Array4<const amrex::Real>& detJ,
                         const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                         const amrex::Array4<const amrex::Real>& mf_m,
                         const amrex::Array4<const amrex::Real>& mf_u,
                         const amrex::Array4<const amrex::Real>& mf_v,
                         const int horiz_adv_type, const int vert_adv_type,
                         const int use_terrain, const int domhi_z);

#endif
