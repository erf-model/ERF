#ifndef _ADVECTION_H_
#define _ADVECTION_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_YAFluxRegister.H>
#include <DataStruct.H>
#include <IndexDefines.H>
#include <ABLMost.H>


/** Compute advection tendency for density and potential temperature */
void AdvectionSrcForRho (const amrex::Box& bx,
                         const amrex::Array4<amrex::Real>& src,
                         const amrex::Array4<const amrex::Real>& rho_u,    // These are being used
                         const amrex::Array4<const amrex::Real>& rho_v,    //  to define the fluxes
                         const amrex::Array4<const amrex::Real>& omega,
                         const amrex::Array4<      amrex::Real>& avg_xmom, // These are being defined
                         const amrex::Array4<      amrex::Real>& avg_ymom, //  from the rho fluxes
                         const amrex::Array4<      amrex::Real>& avg_zmom,
                         const amrex::Array4<const amrex::Real>& z_nd,
                         const amrex::Array4<const amrex::Real>& detJ,
                         const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                         const amrex::Array4<const amrex::Real>& mf_m,
                         const amrex::Array4<const amrex::Real>& mf_u,
                         const amrex::Array4<const amrex::Real>& mf_v,
                         const bool use_terrain,
                         const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_arr
#ifdef ERF_USE_POISSON_SOLVE
                         ,const bool const_rho = false
#endif
                        );

/** Compute advection tendency for all scalars other than density and potential temperature */
void AdvectionSrcForScalars (const amrex::Box& bx,
                             const int icomp, const int ncomp,
                             const amrex::Array4<const amrex::Real>& avg_xmom,
                             const amrex::Array4<const amrex::Real>& avg_ymom,
                             const amrex::Array4<const amrex::Real>& avg_zmom,
                             const amrex::Array4<const amrex::Real>& cell_prim,
                             const amrex::Array4<amrex::Real>& src,
                             const amrex::Array4<const amrex::Real>& detJ,
                             const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                             const amrex::Array4<const amrex::Real>& mf_m,
                             const AdvType horiz_adv_type, const AdvType vert_adv_type,
                             const amrex::Real horiz_upw_frac, const amrex::Real vert_upw_frac,
                             const bool use_terrain,
                             const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_arr);

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
                         const AdvType horiz_adv_type, const AdvType vert_adv_type,
                         const amrex::Real horiz_upw_frac, const amrex::Real vert_upw_frac,
                         const bool use_terrain, const int lo_z_face, const int hi_z_face);

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
AdvType
EfficientAdvType (int nrk, AdvType adv_type)
{
    AdvType eff_adv_type;

    if (nrk == 0) {

        eff_adv_type = AdvType::Centered_2nd;

    } else if (nrk == 1) {

        if ( (adv_type == AdvType::Centered_2nd) ||
             (adv_type == AdvType::Upwind_3rd)   ||
             (adv_type == AdvType::Centered_4th) ||
             (adv_type == AdvType::Weno_3)       ||
             (adv_type == AdvType::Weno_3Z)      ||
             (adv_type == AdvType::Weno_3MZQ)    )
            {
                eff_adv_type = AdvType::Centered_2nd;

            } else if ( (adv_type == AdvType::Upwind_5th) ||
                        (adv_type == AdvType::Weno_5)     ||
                        (adv_type == AdvType::Weno_5Z)    )
            {
                eff_adv_type = AdvType::Upwind_3rd;

            } else { // (adv_type == AdvType::Centered_6th)

            eff_adv_type = AdvType::Centered_4th;
        }

    } else {

        eff_adv_type = adv_type;
    }

    return(eff_adv_type);
}
#endif
