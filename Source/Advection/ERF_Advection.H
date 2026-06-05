#ifndef ERF_ADVECTION_H_
#define ERF_ADVECTION_H_

#include "AMReX.H"
#include "AMReX_MultiFab.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_BCRec.H"
#include "AMReX_YAFluxRegister.H"

#include "ERF_DataStruct.H"
#include "ERF_IndexDefines.H"
#include "ERF_SurfaceLayer.H"
#include "ERF_EB.H"

/** Compute advection tendency for density and potential temperature */
void AdvectionSrcForRho (const amrex::Box& bx,
                         const amrex::Array4<amrex::Real>& src,
                         const amrex::Array4<const amrex::Real>& rho_u,    // These are being used
                         const amrex::Array4<const amrex::Real>& rho_v,    //  to define the fluxes
                         const amrex::Array4<const amrex::Real>& omega,
                         const amrex::Array4<      amrex::Real>& avg_xmom, // These are being defined
                         const amrex::Array4<      amrex::Real>& avg_ymom, //  from the rho fluxes
                         const amrex::Array4<      amrex::Real>& avg_zmom,
                         const amrex::Array4<const amrex::Real>& ax_arr,
                         const amrex::Array4<const amrex::Real>& ay_arr,
                         const amrex::Array4<const amrex::Real>& az_arr,
                         const amrex::Array4<const amrex::Real>& detJ,
                         const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                         const amrex::Array4<const amrex::Real>& mf_mx,
                         const amrex::Array4<const amrex::Real>& mf_my,
                         const amrex::Array4<const amrex::Real>& mf_uy,
                         const amrex::Array4<const amrex::Real>& mf_vx,
                         const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_arr,
                         const bool fixed_rho);

/** Compute advection tendency for all scalars other than density and potential temperature */
void AdvectionSrcForScalars (const amrex::Box& bx,
                             const int icomp, const int ncomp,
                             const amrex::Array4<const amrex::Real>& avg_xmom,
                             const amrex::Array4<const amrex::Real>& avg_ymom,
                             const amrex::Array4<const amrex::Real>& avg_zmom,
                             const amrex::Array4<const amrex::Real>& cell_prim,
                             const amrex::Array4<amrex::Real>& src,
                             const amrex::Array4<const amrex::Real>& vf_arr,
                             const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                             const amrex::Array4<const amrex::Real>& mf_mx,
                             const amrex::Array4<const amrex::Real>& mf_my,
                             const AdvType horiz_adv_type, const AdvType vert_adv_type,
                             const amrex::Real horiz_upw_frac, const amrex::Real vert_upw_frac,
                             const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_arr,
                             const amrex::Box& domain,
                             const amrex::BCRec* bc_ptr_h);

/** Compute advection tendencies for all components of momentum */
void AdvectionSrcForMom (const amrex::MFIter& mfi,
                         const amrex::Box& bx,
                         const amrex::Box& bxx, const amrex::Box& bxy, const amrex::Box& bxz,
                         const amrex::Vector<amrex::Box>& bxx_grown,
                         const amrex::Vector<amrex::Box>& bxy_grown,
                         const amrex::Vector<amrex::Box>& bxz_grown,
                         const amrex::Array4<      amrex::Real>& rho_u_rhs, const amrex::Array4<      amrex::Real>& rho_v_rhs,
                         const amrex::Array4<      amrex::Real>& rho_w_rhs,
                         const amrex::Array4<const amrex::Real>& rho,
                         const amrex::Array4<const amrex::Real>& u        , const amrex::Array4<const amrex::Real>& v,
                         const amrex::Array4<const amrex::Real>& w        ,
                         const amrex::Array4<const amrex::Real>& rho_u    , const amrex::Array4<const amrex::Real>& rho_v,
                         const amrex::Array4<const amrex::Real>& Omega    ,
                         const amrex::Array4<const amrex::Real>& z_nd,
                         const amrex::Array4<const amrex::Real>& ax,
                         const amrex::Array4<const amrex::Real>& ay,
                         const amrex::Array4<const amrex::Real>& az,
                         const amrex::Array4<const amrex::Real>& detJ,
                               amrex::Gpu::DeviceVector<amrex::Real>& stretched_dz_d,
                         const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                         const amrex::Array4<const amrex::Real>& mf_mx,
                         const amrex::Array4<const amrex::Real>& mf_ux,
                         const amrex::Array4<const amrex::Real>& mf_vx,
                         const amrex::Array4<const amrex::Real>& mf_my,
                         const amrex::Array4<const amrex::Real>& mf_uy,
                         const amrex::Array4<const amrex::Real>& mf_vy,
                         const AdvType horiz_adv_type, const AdvType vert_adv_type,
                         const amrex::Real horiz_upw_frac, const amrex::Real vert_upw_frac,
                         MeshType& mesh_type, TerrainType& terrain_type,
                         const eb_& ebfact,
                               amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_u_arr,
                               amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_v_arr,
                               amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_w_arr,
                         const amrex::Vector<amrex::iMultiFab>& physbnd_mask,
                         const bool already_on_centroids,
                         const int lo_z_face, const int hi_z_face,
                         const amrex::Box& domain,
                         const amrex::BCRec* bc_ptr_h);

/** Compute advection tendencies for all components of momentum when dz is constant */
void AdvectionSrcForMom_ConstantDz (
                            const amrex::Box& bxx, const amrex::Box& bxy, const amrex::Box& bxz,
                            const amrex::Array4<      amrex::Real>& rho_u_rhs, const amrex::Array4<      amrex::Real>& rho_v_rhs,
                            const amrex::Array4<      amrex::Real>& rho_w_rhs,
                            const amrex::Array4<const amrex::Real>& u        , const amrex::Array4<const amrex::Real>& v,
                            const amrex::Array4<const amrex::Real>& w        ,
                            const amrex::Array4<const amrex::Real>& rho_u    , const amrex::Array4<const amrex::Real>& rho_v,
                            const amrex::Array4<const amrex::Real>& Omega    ,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                            const amrex::Gpu::DeviceVector<amrex::Real>& stretched_dz_d,
                            const amrex::Array4<const amrex::Real>& mf_mx,
                            const amrex::Array4<const amrex::Real>& mf_ux,
                            const amrex::Array4<const amrex::Real>& mf_vx,
                            const amrex::Array4<const amrex::Real>& mf_my,
                            const amrex::Array4<const amrex::Real>& mf_uy,
                            const amrex::Array4<const amrex::Real>& mf_vy,
                            const AdvType horiz_adv_type, const AdvType vert_adv_type,
                            const amrex::Real horiz_upw_frac, const amrex::Real vert_upw_frac,
                            TerrainType& terrain_type,
                            const int lo_z_face, const int hi_z_face);

/** Compute advection tendencies for all components of momentum when dz is stretched */
void AdvectionSrcForMom_StretchedDz (
                            const amrex::Box& bxx, const amrex::Box& bxy, const amrex::Box& bxz,
                            const amrex::Array4<      amrex::Real>& rho_u_rhs, const amrex::Array4<      amrex::Real>& rho_v_rhs,
                            const amrex::Array4<      amrex::Real>& rho_w_rhs,
                            const amrex::Array4<const amrex::Real>& u        , const amrex::Array4<const amrex::Real>& v,
                            const amrex::Array4<const amrex::Real>& w        ,
                            const amrex::Array4<const amrex::Real>& rho_u    , const amrex::Array4<const amrex::Real>& rho_v,
                            const amrex::Array4<const amrex::Real>& Omega    ,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                            const amrex::Gpu::DeviceVector<amrex::Real>& stretched_dz_d,
                            const amrex::Array4<const amrex::Real>& mf_mx,
                            const amrex::Array4<const amrex::Real>& mf_ux,
                            const amrex::Array4<const amrex::Real>& mf_vx,
                            const amrex::Array4<const amrex::Real>& mf_my,
                            const amrex::Array4<const amrex::Real>& mf_uy,
                            const amrex::Array4<const amrex::Real>& mf_vy,
                            const AdvType horiz_adv_type, const AdvType vert_adv_type,
                            const amrex::Real horiz_upw_frac, const amrex::Real vert_upw_frac,
                            const int lo_z_face, const int hi_z_face);

/** Compute advection tendencies for all components of momentum for TERRAIN-FITTED COORDS only */
void AdvectionSrcForMom_TF (
                            const amrex::Box& bxx, const amrex::Box& bxy, const amrex::Box& bxz,
                            const amrex::Array4<      amrex::Real>& rho_u_rhs, const amrex::Array4<      amrex::Real>& rho_v_rhs,
                            const amrex::Array4<      amrex::Real>& rho_w_rhs,
                            const amrex::Array4<const amrex::Real>& u        , const amrex::Array4<const amrex::Real>& v,
                            const amrex::Array4<const amrex::Real>& w        ,
                            const amrex::Array4<const amrex::Real>& rho_u    , const amrex::Array4<const amrex::Real>& rho_v,
                            const amrex::Array4<const amrex::Real>& Omega    ,
                            const amrex::Array4<const amrex::Real>& z_nd,
                            const amrex::Array4<const amrex::Real>& ax,
                            const amrex::Array4<const amrex::Real>& ay,
                            const amrex::Array4<const amrex::Real>& az,
                            const amrex::Array4<const amrex::Real>& detJ,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                            const amrex::Array4<const amrex::Real>& mf_mx,
                            const amrex::Array4<const amrex::Real>& mf_ux,
                            const amrex::Array4<const amrex::Real>& mf_vx,
                            const amrex::Array4<const amrex::Real>& mf_my,
                            const amrex::Array4<const amrex::Real>& mf_uy,
                            const amrex::Array4<const amrex::Real>& mf_vy,
                            const AdvType horiz_adv_type, const AdvType vert_adv_type,
                            const amrex::Real horiz_upw_frac, const amrex::Real vert_upw_frac,
                            const int lo_z_face, const int hi_z_face);

/** Compute advection tendencies for all components of momentum when using EB for terrain */
void AdvectionSrcForMom_EB (const amrex::MFIter& mfi,
                            const amrex::Box& bxx, const amrex::Box& bxy, const amrex::Box& bxz,
                            const amrex::Vector<amrex::Box>& bxx_grown,
                            const amrex::Vector<amrex::Box>& bxy_grown,
                            const amrex::Vector<amrex::Box>& bxz_grown,
                            const amrex::Array4<      amrex::Real>& rho_u_rhs, const amrex::Array4<      amrex::Real>& rho_v_rhs,
                            const amrex::Array4<      amrex::Real>& rho_w_rhs,
                            const amrex::Array4<const amrex::Real>& u        , const amrex::Array4<const amrex::Real>& v,
                            const amrex::Array4<const amrex::Real>& w        ,
                            const amrex::Array4<const amrex::Real>& rho_u    , const amrex::Array4<const amrex::Real>& rho_v,
                            const amrex::Array4<const amrex::Real>& Omega    ,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                            const amrex::Array4<const amrex::Real>& mf_mx,
                            const amrex::Array4<const amrex::Real>& mf_ux,
                            const amrex::Array4<const amrex::Real>& mf_vx,
                            const amrex::Array4<const amrex::Real>& mf_my,
                            const amrex::Array4<const amrex::Real>& mf_uy,
                            const amrex::Array4<const amrex::Real>& mf_vy,
                            const AdvType horiz_adv_type, const AdvType vert_adv_type,
                            const amrex::Real horiz_upw_frac, const amrex::Real vert_upw_frac,
                            const eb_& ebfact,
                                  amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_u_arr,
                                  amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_v_arr,
                                  amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_w_arr,
                            const amrex::Vector<amrex::iMultiFab>& physbnd_mask,
                            const bool already_on_centroids,
                            const int lo_z_face, const int hi_z_face,
                            const amrex::Box& domain);

/** Compute advection tendencies for all normal components of momentum */
void
AdvectionSrcForOpenBC_Normal (const amrex::Box& bx,
                              const int& dir,
                              const amrex::Array4<      amrex::Real>& rhs_arr,
                              const amrex::Array4<const amrex::Real>& vel_norm_arr,
                              const amrex::Array4<const amrex::Real>& cell_data_arr,
                              const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxInv,
                              const bool do_lo=false);

/** Compute advection tendencies for u momentum tangent to open bc*/
void
AdvectionSrcForOpenBC_Tangent_Xmom (const amrex::Box& bxx,
                                    const int& dir,
                                    const amrex::Array4<      amrex::Real>& rho_u_rhs,
                                    const amrex::Array4<const amrex::Real>& u,
                                    const amrex::Array4<const amrex::Real>& rho_u,
                                    const amrex::Array4<const amrex::Real>& rho_v,
                                    const amrex::Array4<const amrex::Real>& Omega,
                                    const amrex::Array4<const amrex::Real>& ax,
                                    const amrex::Array4<const amrex::Real>& az,
                                    const amrex::Array4<const amrex::Real>& detJ,
                                    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                                    const bool do_lo=false);

/** Compute advection tendencies for v momentum tangent to open bc*/
void
AdvectionSrcForOpenBC_Tangent_Ymom (const amrex::Box& bxy,
                                    const int& dir,
                                    const amrex::Array4<      amrex::Real>& rho_v_rhs,
                                    const amrex::Array4<const amrex::Real>& v,
                                    const amrex::Array4<const amrex::Real>& rho_u,
                                    const amrex::Array4<const amrex::Real>& rho_v,
                                    const amrex::Array4<const amrex::Real>& Omega,
                                    const amrex::Array4<const amrex::Real>& ay,
                                    const amrex::Array4<const amrex::Real>& az,
                                    const amrex::Array4<const amrex::Real>& detJ,
                                    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                                    const bool do_lo=false);

/** Compute advection tendencies for w momentum tangent to open bc*/
void
AdvectionSrcForOpenBC_Tangent_Zmom (const amrex::Box& bxz,
                                    const int& dir,
                                    const amrex::Array4<      amrex::Real>& rho_w_rhs,
                                    const amrex::Array4<const amrex::Real>& w,
                                    const amrex::Array4<const amrex::Real>& rho_u,
                                    const amrex::Array4<const amrex::Real>& rho_v,
                                    const amrex::Array4<const amrex::Real>& Omega,
                                    const amrex::Array4<const amrex::Real>& ax,
                                    const amrex::Array4<const amrex::Real>& ay,
                                    const amrex::Array4<const amrex::Real>& az,
                                    const amrex::Array4<const amrex::Real>& detJ,
                                    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                                    const int domhi_z,
                                    const bool do_lo=false);

/** Compute advection tendencies for cons tangential to BC (2nd order)*/
void
AdvectionSrcForOpenBC_Tangent_Cons (const amrex::Box& bx,
                                    const int& dir,
                                    const int& icomp,
                                    const int& ncomp,
                                    const amrex::Array4<      amrex::Real>& cell_rhs,
                                    const amrex::Array4<const amrex::Real>& cell_prim,
                                    const amrex::Array4<const amrex::Real>& avg_xmom,
                                    const amrex::Array4<const amrex::Real>& avg_ymom,
                                    const amrex::Array4<const amrex::Real>& avg_zmom,
                                    const amrex::Array4<const amrex::Real>& detJ,
                                    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                                    const bool do_lo=false);

/** Compute advection tendencies in normal dir for vars tangent to open bc */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForOpenBC_Tangent (const int& i,
                               const int& j,
                               const int& k,
                               const int& nprim,
                               const int& dir,
                               const amrex::Array4<const amrex::Real>& prim_tang_arr,
                               const amrex::Array4<const amrex::Real>& mom_norm_arr,
                               const amrex::Real& dxInv,
                               const bool do_lo=false);

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
