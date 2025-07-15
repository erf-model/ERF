#include "AMReX_BCRec.H"

#include <ERF_Advection.H>
#include <ERF_AdvectionSrcForMom_N.H>
#include <ERF_AdvectionSrcForMom_T.H>

using namespace amrex;

/**
 * Function for computing the advective tendency for the momentum equations
 * ONLY when using terrain-fitted coordinates
 *
 * @param[in] bxx box over which the x-momentum is updated
 * @param[in] bxy box over which the y-momentum is updated
 * @param[in] bxz box over which the z-momentum is updated
 * @param[out] rho_u_rhs tendency for the x-momentum equation
 * @param[out] rho_v_rhs tendency for the y-momentum equation
 * @param[out] rho_w_rhs tendency for the z-momentum equation
 * @param[in] u x-component of the velocity
 * @param[in] v y-component of the velocity
 * @param[in] w z-component of the velocity
 * @param[in] rho_u x-component of the momentum
 * @param[in] rho_v y-component of the momentum
 * @param[in] Omega component of the momentum normal to the z-coordinate surface
 * @param[in] z_nd height coordinate at nodes
 * @param[in] ax   Area fraction of x-faces
 * @param[in] ay   Area fraction of y-faces
 * @param[in] az   Area fraction of z-faces
 * @param[in] detJ Jacobian of the metric transformation
 * @param[in] cellSizeInv inverse of the mesh spacing
 * @param[in] mf_m map factor at cell centers
 * @param[in] mf_u map factor at x-faces
 * @param[in] mf_v map factor at y-faces
 * @param[in] horiz_adv_type sets the spatial order to be used for lateral derivatives
 * @param[in] vert_adv_type  sets the spatial order to be used for vertical derivatives
 */
void
AdvectionSrcForMom_TF (const Box& bxx, const Box& bxy, const Box& bxz,
                       const Array4<      Real>& rho_u_rhs,
                       const Array4<      Real>& rho_v_rhs,
                       const Array4<      Real>& rho_w_rhs,
                       const Array4<const Real>& u,
                       const Array4<const Real>& v,
                       const Array4<const Real>& w,
                       const Array4<const Real>& rho_u,
                       const Array4<const Real>& rho_v,
                       const Array4<const Real>& Omega,
                       const Array4<const Real>& z_nd,
                       const Array4<const Real>& ax,
                       const Array4<const Real>& ay,
                       const Array4<const Real>& az,
                       const Array4<const Real>& detJ,
                       const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                       const Array4<const Real>& mf_mx,
                       const Array4<const Real>& mf_ux,
                       const Array4<const Real>& mf_vx,
                       const Array4<const Real>& mf_my,
                       const Array4<const Real>& mf_uy,
                       const Array4<const Real>& mf_vy,
                       const AdvType horiz_adv_type,
                       const AdvType vert_adv_type,
                       const Real horiz_upw_frac,
                       const Real vert_upw_frac,
                       const int lo_z_face, const int hi_z_face)
{
    BL_PROFILE_VAR("AdvectionSrcForMom_TF", AdvectionSrcForMom_TF);

    AMREX_ALWAYS_ASSERT(bxz.smallEnd(2) > 0);

    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    // Compute mapfactor inverses
    Box box2d_u(bxx);   box2d_u.setRange(2,0);   box2d_u.grow({3,3,0});
    Box box2d_v(bxy);   box2d_v.setRange(2,0);   box2d_v.grow({3,3,0});

    FArrayBox mf_ux_invFAB(box2d_u,1,The_Async_Arena());
    FArrayBox mf_uy_invFAB(box2d_u,1,The_Async_Arena());
    const Array4<Real>& mf_ux_inv = mf_ux_invFAB.array();
    const Array4<Real>& mf_uy_inv = mf_uy_invFAB.array();

    FArrayBox mf_vx_invFAB(box2d_v,1,The_Async_Arena());
    FArrayBox mf_vy_invFAB(box2d_v,1,The_Async_Arena());
    const Array4<Real>& mf_vx_inv = mf_vx_invFAB.array();
    const Array4<Real>& mf_vy_inv = mf_vy_invFAB.array();

    ParallelFor(box2d_u, box2d_v,
    [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        mf_ux_inv(i,j,0) = 1. / mf_ux(i,j,0);
        mf_uy_inv(i,j,0) = 1. / mf_uy(i,j,0);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        mf_vx_inv(i,j,0) = 1. / mf_vx(i,j,0);
        mf_vy_inv(i,j,0) = 1. / mf_vy(i,j,0);
    });

    // Inline with 2nd order for efficiency
    if (horiz_adv_type == AdvType::Centered_2nd && vert_adv_type == AdvType::Centered_2nd)
    {
            ParallelFor(bxx, bxy, bxz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real xflux_hi = 0.25 * (rho_u(i,j,k) * mf_uy_inv(i,j,0) + rho_u(i+1,j,k) * mf_uy_inv(i+1,j,0)) *
                                       (u(i+1,j,k) + u(i,j,k)) * 0.5 * (ax(i,j,k) + ax(i+1,j,k));

                Real xflux_lo = 0.25 * (rho_u(i,j,k) * mf_uy_inv(i,j,0) + rho_u(i-1,j,k) * mf_uy_inv(i-1,j,0)) *
                                       (u(i-1,j,k) + u(i,j,k)) * 0.5 * (ax(i,j,k) + ax(i-1,j,k));

                Real met_h_zeta_yhi = Compute_h_zeta_AtEdgeCenterK(i,j+1,k,cellSizeInv,z_nd);
                Real yflux_hi = 0.25 * (rho_v(i,j+1,k)*mf_vx_inv(i,j+1,0) + rho_v(i-1,j+1,k)*mf_vx_inv(i-1,j+1,0)) *
                                       (u(i,j+1,k) + u(i,j,k)) * met_h_zeta_yhi;

                Real met_h_zeta_ylo = Compute_h_zeta_AtEdgeCenterK(i,j  ,k,cellSizeInv,z_nd);
                Real yflux_lo = 0.25 * (rho_v(i,j  ,k)*mf_vx_inv(i,j  ,0) + rho_v(i-1,j  ,k)*mf_vx_inv(i-1,j  ,0)) *
                                       (u(i,j-1,k) + u(i,j,k)) * met_h_zeta_ylo;

                Real zflux_hi = 0.25 * (Omega(i,j,k+1) + Omega(i-1,j,k+1)) * (u(i,j,k+1) + u(i,j,k)) *
                                        0.5 * (az(i,j,k+1) + az(i-1,j,k+1));
                Real zflux_lo = 0.25 * (Omega(i,j,k  ) + Omega(i-1,j,k  )) * (u(i,j,k-1) + u(i,j,k)) *
                                        0.5 * (az(i,j,k  ) + az(i-1,j,k  ));

                Real mfsq = mf_ux(i,j,0) * mf_uy(i,j,0);

                Real advectionSrc = (xflux_hi - xflux_lo) * dxInv * mfsq
                                  + (yflux_hi - yflux_lo) * dyInv * mfsq
                                  + (zflux_hi - zflux_lo) * dzInv;
                rho_u_rhs(i, j, k) = -advectionSrc / (0.5 * (detJ(i,j,k) + detJ(i-1,j,k)));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

                Real met_h_zeta_xhi = Compute_h_zeta_AtEdgeCenterK(i+1,j,k,cellSizeInv,z_nd);
                Real xflux_hi = 0.25 * (rho_u(i+1,j,k)*mf_uy_inv(i+1,j,0) + rho_u(i+1,j-1,k)*mf_uy_inv(i+1,j-1,0)) *
                                       (v(i+1,j,k) + v(i,j,k)) * met_h_zeta_xhi;

                Real met_h_zeta_xlo = Compute_h_zeta_AtEdgeCenterK(i  ,j,k,cellSizeInv,z_nd);
                Real xflux_lo = 0.25 * (rho_u(i, j, k)*mf_uy_inv(i  ,j,0) + rho_u(i  ,j-1,k)*mf_uy_inv(i-1,j  ,0)) *
                                       (v(i-1,j,k) + v(i,j,k)) * met_h_zeta_xlo;

                Real yflux_hi = 0.25 * (rho_v(i,j+1,k)*mf_vx_inv(i,j+1,0) + rho_v(i,j  ,k) * mf_vx_inv(i,j  ,0)) *
                                       (v(i,j+1,k) + v(i,j,k)) * 0.5 * (ay(i,j,k) + ay(i,j+1,k));

                Real yflux_lo = 0.25 * (rho_v(i,j  ,k)*mf_vx_inv(i,j  ,0) + rho_v(i,j-1,k) * mf_vx_inv(i,j-1,0)) *
                                       (v(i,j-1,k) + v(i,j,k)) * 0.5 * (ay(i,j,k) + ay(i,j-1,k));

                Real zflux_hi = 0.25 * (Omega(i,j,k+1) + Omega(i, j-1, k+1)) * (v(i,j,k+1) + v(i,j,k)) *
                                        0.5 * (az(i,j,k+1) + az(i,j-1,k+1));
                Real zflux_lo = 0.25 * (Omega(i,j,k  ) + Omega(i, j-1, k  )) * (v(i,j,k-1) + v(i,j,k)) *
                                        0.5 * (az(i,j,k  ) + az(i,j-1,k  ));

                Real mfsq = mf_vx(i,j,0) * mf_vy(i,j,0);

                Real advectionSrc = (xflux_hi - xflux_lo) * dxInv * mfsq
                                  + (yflux_hi - yflux_lo) * dyInv * mfsq
                                  + (zflux_hi - zflux_lo) * dzInv;
                rho_v_rhs(i, j, k) = -advectionSrc / (0.5 * (detJ(i,j,k) + detJ(i,j-1,k)));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
                Real met_h_zeta_xhi = Compute_h_zeta_AtEdgeCenterJ(i+1,j  ,k  ,cellSizeInv,z_nd);
                Real xflux_hi = 0.25*(rho_u(i+1,j  ,k) + rho_u(i+1,j,k-1)) * mf_uy_inv(i+1,j,0) *
                                     (w(i+1,j,k) + w(i,j,k)) * met_h_zeta_xhi;

                Real met_h_zeta_xlo = Compute_h_zeta_AtEdgeCenterJ(i  ,j  ,k  ,cellSizeInv,z_nd);
                Real xflux_lo = 0.25*(rho_u(i  ,j  ,k) + rho_u(i  ,j,k-1)) * mf_uy_inv(i  ,j,0) *
                                     (w(i-1,j,k) + w(i,j,k)) * met_h_zeta_xlo;

                Real met_h_zeta_yhi = Compute_h_zeta_AtEdgeCenterI(i  ,j+1,k  ,cellSizeInv,z_nd);
                Real yflux_hi = 0.25*(rho_v(i,j+1,k) + rho_v(i,j+1,k-1)) * mf_vx_inv(i,j+1,0) *
                                     (w(i,j+1,k) + w(i,j,k)) * met_h_zeta_yhi;

                Real met_h_zeta_ylo = Compute_h_zeta_AtEdgeCenterI(i  ,j  ,k  ,cellSizeInv,z_nd);
                Real yflux_lo = 0.25*(rho_v(i,j  ,k) + rho_v(i,j  ,k-1)) * mf_vx_inv(i,j  ,0) *
                                     (w(i,j-1,k) + w(i,j,k)) * met_h_zeta_ylo;

                Real zflux_lo = 0.25 * (Omega(i,j,k) + Omega(i,j,k-1)) * (w(i,j,k) + w(i,j,k-1));

                Real zflux_hi = (k == hi_z_face) ? Omega(i,j,k) * w(i,j,k)  * az(i,j,k):
                    0.25 * (Omega(i,j,k) + Omega(i,j,k+1)) * (w(i,j,k) + w(i,j,k+1)) *
                    0.5  * (az(i,j,k) + az(i,j,k+1));

                Real mfsq = mf_mx(i,j,0) * mf_my(i,j,0);

                Real advectionSrc = (xflux_hi - xflux_lo) * dxInv * mfsq
                                  + (yflux_hi - yflux_lo) * dyInv * mfsq
                                  + (zflux_hi - zflux_lo) * dzInv;
                rho_w_rhs(i, j, k) = -advectionSrc / (0.5*(detJ(i,j,k) + detJ(i,j,k-1)));
            });

    // Template higher order methods
    } else {
        if (horiz_adv_type == AdvType::Centered_2nd) {
                AdvectionSrcForMomVert<CENTERED2>(bxx, bxy, bxz,
                                                  rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                  rho_u, rho_v, Omega, u, v, w,
                                                  z_nd, ax, ay, az, detJ, cellSizeInv,
                                                  mf_mx, mf_ux_inv, mf_vx_inv,
                                                  mf_my, mf_uy_inv, mf_vy_inv,
                                                  horiz_upw_frac, vert_upw_frac,
                                                  vert_adv_type, lo_z_face, hi_z_face);
        } else if (horiz_adv_type == AdvType::Upwind_3rd) {
                AdvectionSrcForMomVert<UPWIND3>(bxx, bxy, bxz,
                                                rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                rho_u, rho_v, Omega, u, v, w,
                                                z_nd, ax, ay, az, detJ, cellSizeInv,
                                                mf_mx, mf_ux_inv, mf_vx_inv,
                                                mf_my, mf_uy_inv, mf_vy_inv,
                                                horiz_upw_frac, vert_upw_frac,
                                                vert_adv_type, lo_z_face, hi_z_face);
        } else if (horiz_adv_type == AdvType::Centered_4th) {
                AdvectionSrcForMomVert<CENTERED4>(bxx, bxy, bxz,
                                                  rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                  rho_u, rho_v, Omega, u, v, w,
                                                  z_nd, ax, ay, az, detJ, cellSizeInv,
                                                  mf_mx, mf_ux_inv, mf_vx_inv,
                                                  mf_my, mf_uy_inv, mf_vy_inv,
                                                  horiz_upw_frac, vert_upw_frac,
                                                  vert_adv_type, lo_z_face, hi_z_face);
        } else if (horiz_adv_type == AdvType::Upwind_5th) {
                AdvectionSrcForMomVert<UPWIND5>(bxx, bxy, bxz,
                                                rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                rho_u, rho_v, Omega, u, v, w,
                                                z_nd, ax, ay, az, detJ, cellSizeInv,
                                                mf_mx, mf_ux_inv, mf_vx_inv,
                                                mf_my, mf_uy_inv, mf_vy_inv,
                                                horiz_upw_frac, vert_upw_frac,
                                                vert_adv_type, lo_z_face, hi_z_face);
        } else if (horiz_adv_type == AdvType::Centered_6th) {
                AdvectionSrcForMomVert<CENTERED6>(bxx, bxy, bxz,
                                                  rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                  rho_u, rho_v, Omega, u, v, w,
                                                  z_nd, ax, ay, az, detJ, cellSizeInv,
                                                  mf_mx, mf_ux_inv, mf_vx_inv,
                                                  mf_my, mf_uy_inv, mf_vy_inv,
                                                  horiz_upw_frac, vert_upw_frac,
                                                  vert_adv_type, lo_z_face, hi_z_face);
        } else if (horiz_adv_type == AdvType::Weno_3) {
                AdvectionSrcForMomVert<WENO3>(bxx, bxy, bxz,
                                            rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                            rho_u, rho_v, Omega, u, v, w,
                                            z_nd, ax, ay, az, detJ, cellSizeInv,
                                            mf_mx, mf_ux_inv, mf_vx_inv,
                                            mf_my, mf_uy_inv, mf_vy_inv,
                                            horiz_upw_frac, vert_upw_frac,
                                            vert_adv_type, lo_z_face, hi_z_face);
        } else if (horiz_adv_type == AdvType::Weno_3Z) {
                AdvectionSrcForMomVert<WENO_Z3>(bxx, bxy, bxz,
                                                rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                rho_u, rho_v, Omega, u, v, w,
                                                z_nd, ax, ay, az, detJ, cellSizeInv,
                                                mf_mx, mf_ux_inv, mf_vx_inv,
                                                mf_my, mf_uy_inv, mf_vy_inv,
                                                horiz_upw_frac, vert_upw_frac,
                                                vert_adv_type, lo_z_face, hi_z_face);
        } else if (horiz_adv_type == AdvType::Weno_3MZQ) {
                AdvectionSrcForMomVert<WENO_MZQ3>(bxx, bxy, bxz,
                                                rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                rho_u, rho_v, Omega, u, v, w,
                                                z_nd, ax, ay, az, detJ, cellSizeInv,
                                                mf_mx, mf_ux_inv, mf_vx_inv,
                                                mf_my, mf_uy_inv, mf_vy_inv,
                                                horiz_upw_frac, vert_upw_frac,
                                                vert_adv_type, lo_z_face, hi_z_face);
        } else if (horiz_adv_type == AdvType::Weno_5) {
                AdvectionSrcForMomVert<WENO5>(bxx, bxy, bxz,
                                            rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                            rho_u, rho_v, Omega, u, v, w,
                                            z_nd, ax, ay, az, detJ, cellSizeInv,
                                            mf_mx, mf_ux_inv, mf_vx_inv,
                                            mf_my, mf_uy_inv, mf_vy_inv,
                                            horiz_upw_frac, vert_upw_frac,
                                            vert_adv_type, lo_z_face, hi_z_face);
        } else if (horiz_adv_type == AdvType::Weno_5Z) {
                AdvectionSrcForMomVert<WENO_Z5>(bxx, bxy, bxz,
                                                rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                rho_u, rho_v, Omega, u, v, w,
                                                z_nd, ax, ay, az, detJ, cellSizeInv,
                                                mf_mx, mf_ux_inv, mf_vx_inv,
                                                mf_my, mf_uy_inv, mf_vy_inv,
                                                horiz_upw_frac, vert_upw_frac,
                                                vert_adv_type, lo_z_face, hi_z_face);
        } else if (horiz_adv_type == AdvType::Weno_7) {
                AdvectionSrcForMomVert<WENO7>(bxx, bxy, bxz,
                                            rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                            rho_u, rho_v, Omega, u, v, w,
                                            z_nd, ax, ay, az, detJ, cellSizeInv,
                                            mf_mx, mf_ux_inv, mf_vx_inv,
                                            mf_my, mf_uy_inv, mf_vy_inv,
                                            horiz_upw_frac, vert_upw_frac,
                                            vert_adv_type, lo_z_face, hi_z_face);
        } else if (horiz_adv_type == AdvType::Weno_7Z) {
                AdvectionSrcForMomVert<WENO_Z7>(bxx, bxy, bxz,
                                                rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                rho_u, rho_v, Omega, u, v, w,
                                                z_nd, ax, ay, az, detJ, cellSizeInv,
                                                mf_mx, mf_ux_inv, mf_vx_inv,
                                                mf_my, mf_uy_inv, mf_vy_inv,
                                                horiz_upw_frac, vert_upw_frac,
                                                vert_adv_type, lo_z_face, hi_z_face);
        } else {
            AMREX_ASSERT_WITH_MESSAGE(false, "Unknown advection scheme!");
        }
    }
}
