#include "AMReX_BCRec.H"

#include <ERF_Advection.H>
#include <ERF_AdvectionSrcForMom_N.H>
#include <ERF_AdvectionSrcForMom_T.H>

using namespace amrex;

/**
 * Function for computing the advective tendency for the momentum equations
 * when using EB
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
 * @param[in] vf   Volume fraction
 * @param[in] cellSizeInv inverse of the mesh spacing
 * @param[in] mf_m map factor at cell centers
 * @param[in] mf_u map factor at x-faces
 * @param[in] mf_v map factor at y-faces
 * @param[in] horiz_adv_type sets the spatial order to be used for lateral derivatives
 * @param[in] vert_adv_type  sets the spatial order to be used for vertical derivatives
 */
void
AdvectionSrcForMom_EB (const Box& bxx, const Box& bxy, const Box& bxz,
                       const Array4<      Real>& rho_u_rhs,
                       const Array4<      Real>& rho_v_rhs,
                       const Array4<      Real>& rho_w_rhs,
                       const Array4<const Real>& /*u*/,
                       const Array4<const Real>& /*v*/,
                       const Array4<const Real>& /*w*/,
                       const Array4<const Real>& /*rho_u*/,
                       const Array4<const Real>& /*rho_v*/,
                       const Array4<const Real>& /*omega*/,
                       const Array4<const Real>& /*ax*/,
                       const Array4<const Real>& /*ay*/,
                       const Array4<const Real>& /*az*/,
                       const Array4<const Real>& /*vf*/,
                       const GpuArray<Real, AMREX_SPACEDIM>& /*cellSizeInv*/,
                       const Array4<const Real>& /*mf_m*/,
                       const Array4<const Real>& /*mf_u*/,
                       const Array4<const Real>& /*mf_v*/,
                       const AdvType /*horiz_adv_type*/,
                       const AdvType /*vert_adv_type*/,
                       const Real /*horiz_upw_frac*/,
                       const Real /*vert_upw_frac*/,
                       const int /*lo_z_face*/, const int /*hi_z_face*/,
                       const Box& /*domain*/)
{
    BL_PROFILE_VAR("AdvectionSrcForMom_EB", AdvectionSrcForMom_EB);

//     auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    AMREX_ALWAYS_ASSERT(bxz.smallEnd(2) > 0);

#if 0
    // compute mapfactor inverses
    Box box2d_u(bxx);   box2d_u.setRange(2,0);   box2d_u.grow({3,3,0});
    Box box2d_v(bxy);   box2d_v.setRange(2,0);   box2d_v.grow({3,3,0});

    FArrayBox mf_u_invFAB(box2d_u,1,The_Async_Arena());
    FArrayBox mf_v_invFAB(box2d_v,1,The_Async_Arena());

    const Array4<Real>& mf_u_inv = mf_u_invFAB.array();
    const Array4<Real>& mf_v_inv = mf_v_invFAB.array();

    ParallelFor(box2d_u, box2d_v,
    [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        mf_u_inv(i,j,0) = 1. / mf_u(i,j,0);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        mf_v_inv(i,j,0) = 1. / mf_v(i,j,0);
    });
#endif

    ParallelFor(bxx, bxy, bxz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        rho_u_rhs(i, j, k) = 0.0;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        rho_v_rhs(i, j, k) = 0.0;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        rho_w_rhs(i, j, k) = 0.0;
    });
}
