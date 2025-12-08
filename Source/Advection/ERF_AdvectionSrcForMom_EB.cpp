#include "AMReX_BCRec.H"

#include <ERF_Advection.H>
#include <ERF_AdvectionSrcForMom_N.H>
#include <ERF_AdvectionSrcForMom_T.H>
#include <ERF_EBAdvectionSrcForMom.H>

using namespace amrex;

/**
 * Function for computing the advective tendency for the momentum equations
 * when using EB
 *
 * @param[in] bxx box over which the x-momentum is updated
 * @param[in] bxy box over which the y-momentum is updated
 * @param[in] bxz box over which the z-momentum is updated
 * @param[in] bxx_grown grown boxes of bxx to loop over the nodal grids of bxx
 * @param[in] bxy_grown grown boxes of bxy to loop over the nodal grids of bxy
 * @param[in] bxz_grown grown boxes of bxz to loop over the nodal grids of bxz
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
 * @param[in] cellSizeInv inverse of the mesh spacing
 * @param[in] mf_m map factor at cell centers
 * @param[in] mf_u map factor at x-faces
 * @param[in] mf_v map factor at y-faces
 * @param[in] ebfact EB factories for cell- and face-centered variables
 * @param[in] flx_u_arr fluxes for x-momentum
 * @param[in] flx_v_arr fluxes for y-momentum
 * @param[in] flx_w_arr fluxes for z-momentum
 * @param[in] physbnd_mask Vector of masks for flux interpolation (=1 otherwise, =0 if physbnd)
 * @param[in] already_on_centroids flag whether flux interpolation is unnecessary
 * @param[in] horiz_adv_type sets the spatial order to be used for lateral derivatives
 * @param[in] vert_adv_type  sets the spatial order to be used for vertical derivatives
 */
void
AdvectionSrcForMom_EB ( const MFIter& mfi,
                        const Box& bxx, const Box& bxy, const Box& bxz,
                        const Vector<Box>& bxx_grown,
                        const Vector<Box>& bxy_grown,
                        const Vector<Box>& bxz_grown,
                        const Array4<      Real>& rho_u_rhs,
                        const Array4<      Real>& rho_v_rhs,
                        const Array4<      Real>& rho_w_rhs,
                        const Array4<const Real>& u,
                        const Array4<const Real>& v,
                        const Array4<const Real>& w,
                        const Array4<const Real>& rho_u,
                        const Array4<const Real>& rho_v,
                        const Array4<const Real>& omega,
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
                        const eb_& ebfact,
                              GpuArray<Array4<Real>, AMREX_SPACEDIM>& flx_u_arr,
                              GpuArray<Array4<Real>, AMREX_SPACEDIM>& flx_v_arr,
                              GpuArray<Array4<Real>, AMREX_SPACEDIM>& flx_w_arr,
                        const Vector<iMultiFab>& physbnd_mask,
                        const bool already_on_centroids,
                        const int lo_z_face, const int hi_z_face,
                        const Box& /*domain*/)
{
    BL_PROFILE_VAR("AdvectionSrcForMom_EB", AdvectionSrcForMom_EB);

    AMREX_ALWAYS_ASSERT(bxz.smallEnd(2) > 0);

    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    // compute mapfactor inverses
    Box box2d_u(bxx);   box2d_u.setRange(2,0);   box2d_u.grow({3,3,0});
    Box box2d_v(bxy);   box2d_v.setRange(2,0);   box2d_v.grow({3,3,0});

    FArrayBox mf_ux_invFAB(box2d_u,1,The_Async_Arena());
    FArrayBox mf_uy_invFAB(box2d_u,1,The_Async_Arena());
    FArrayBox mf_vx_invFAB(box2d_v,1,The_Async_Arena());
    FArrayBox mf_vy_invFAB(box2d_v,1,The_Async_Arena());
    const Array4<Real>& mf_ux_inv = mf_ux_invFAB.array();
    const Array4<Real>& mf_uy_inv = mf_uy_invFAB.array();
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

    // EB u-factory
    Array4<const EBCellFlag> u_cflag   = (ebfact.get_u_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();
    Array4<const Real      > u_vfrac   = (ebfact.get_u_const_factory())->getVolFrac().const_array(mfi);
    Array4<const Real      > u_afrac_x = (ebfact.get_u_const_factory())->getAreaFrac()[0]->const_array(mfi);
    Array4<const Real      > u_afrac_y = (ebfact.get_u_const_factory())->getAreaFrac()[1]->const_array(mfi);
    Array4<const Real      > u_afrac_z = (ebfact.get_u_const_factory())->getAreaFrac()[2]->const_array(mfi);
    Array4<const Real      > u_fcx     = (ebfact.get_u_const_factory())->getFaceCent()[0]->const_array(mfi);
    Array4<const Real      > u_fcy     = (ebfact.get_u_const_factory())->getFaceCent()[1]->const_array(mfi);
    Array4<const Real      > u_fcz     = (ebfact.get_u_const_factory())->getFaceCent()[2]->const_array(mfi);

    // EB v-factory
    Array4<const EBCellFlag> v_cflag   = (ebfact.get_v_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();
    Array4<const Real      > v_vfrac   = (ebfact.get_v_const_factory())->getVolFrac().const_array(mfi);
    Array4<const Real      > v_afrac_x = (ebfact.get_v_const_factory())->getAreaFrac()[0]->const_array(mfi);
    Array4<const Real      > v_afrac_y = (ebfact.get_v_const_factory())->getAreaFrac()[1]->const_array(mfi);
    Array4<const Real      > v_afrac_z = (ebfact.get_v_const_factory())->getAreaFrac()[2]->const_array(mfi);
    Array4<const Real      > v_fcx     = (ebfact.get_v_const_factory())->getFaceCent()[0]->const_array(mfi);
    Array4<const Real      > v_fcy     = (ebfact.get_v_const_factory())->getFaceCent()[1]->const_array(mfi);
    Array4<const Real      > v_fcz     = (ebfact.get_v_const_factory())->getFaceCent()[2]->const_array(mfi);

    // EB w-factory
    Array4<const EBCellFlag> w_cflag   = (ebfact.get_w_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();
    Array4<const Real      > w_vfrac   = (ebfact.get_w_const_factory())->getVolFrac().const_array(mfi);
    Array4<const Real      > w_afrac_x = (ebfact.get_w_const_factory())->getAreaFrac()[0]->const_array(mfi);
    Array4<const Real      > w_afrac_y = (ebfact.get_w_const_factory())->getAreaFrac()[1]->const_array(mfi);
    Array4<const Real      > w_afrac_z = (ebfact.get_w_const_factory())->getAreaFrac()[2]->const_array(mfi);
    Array4<const Real      > w_fcx     = (ebfact.get_w_const_factory())->getFaceCent()[0]->const_array(mfi);
    Array4<const Real      > w_fcy     = (ebfact.get_w_const_factory())->getFaceCent()[1]->const_array(mfi);
    Array4<const Real      > w_fcz     = (ebfact.get_w_const_factory())->getFaceCent()[2]->const_array(mfi);

    // Inline with 2nd order for efficiency
    if (horiz_adv_type == AdvType::Centered_2nd && vert_adv_type == AdvType::Centered_2nd)
    {
        // Fluxes for x-momentum
        ParallelFor(bxx_grown[0], bxx_grown[1], bxx_grown[2],
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if ( u_afrac_x(i,j,k)>0.){
                flx_u_arr[0](i,j,k) = 0.25 * u_afrac_x(i,j,k)
                                    * (rho_u(i,j,k) * mf_ux_inv(i,j,0) + rho_u(i-1,j,k) * mf_ux_inv(i-1,j,0))
                                    * (u(i-1,j,k) + u(i,j,k));
            } else {
                flx_u_arr[0](i,j,k) = 0.;
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if ( u_afrac_y(i,j,k)>0.){
                flx_u_arr[1](i,j,k) = 0.25 * u_afrac_y(i,j,k)
                                    * (rho_v(i,j,k) * mf_vy_inv(i,j,0) + rho_v(i-1,j,k) * mf_vy_inv(i-1,j,0))
                                    * (u(i,j-1,k) + u(i,j,k));
            } else {
                flx_u_arr[1](i,j,k) = 0.;
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if ( u_afrac_z(i,j,k)>0.){
                flx_u_arr[2](i,j,k) = 0.25 * (omega(i,j,k) + omega(i-1,j,k)) * (u(i,j,k-1) + u(i,j,k));
            } else {
                flx_u_arr[2](i,j,k) = 0.;
            }
        });
        // Fluxes for y-momentum
        ParallelFor(bxy_grown[0], bxy_grown[1], bxy_grown[2],
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if ( v_afrac_x(i,j,k)>0.){
                flx_v_arr[0](i,j,k) = 0.25 * v_afrac_x(i,j,k)
                                    * (rho_u(i,j,k) * mf_uy_inv(i,j,0) + rho_u(i,j-1,k) * mf_uy_inv(i,j-1,0))
                                    * (v(i-1,j,k) + v(i,j,k));
            } else {
                flx_v_arr[0](i,j,k) = 0.;
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if ( v_afrac_y(i,j,k)>0.){
                flx_v_arr[1](i,j,k) = 0.25 * v_afrac_y(i,j,k)
                                    * (rho_v(i,j,k) * mf_vy_inv(i,j,0) + rho_v(i,j-1,k) * mf_vy_inv(i,j-1,0))
                                    * (v(i,j-1,k) + v(i,j,k));
            } else {
                flx_v_arr[1](i,j,k) = 0.;
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if ( v_afrac_z(i,j,k)>0.){
                flx_v_arr[2](i,j,k) = 0.25 * (omega(i,j,k) + omega(i,j-1,k)) * (v(i,j,k-1) + v(i,j,k));
            } else {
                flx_v_arr[2](i,j,k) = 0.;
            }
        });
        // Fluxes for z-momentum
        ParallelFor(bxz_grown[0], bxz_grown[1], bxz_grown[2],
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if ( w_afrac_x(i,j,k)>0.){
                flx_w_arr[0](i,j,k) = 0.25 * w_afrac_x(i,j,k)
                                    * (rho_u(i,j,k) + rho_u(i,j, k-1)) * mf_ux_inv(i,j,0)
                                    * (w(i-1,j,k) + w(i,j,k));
            } else {
                flx_w_arr[0](i,j,k) = 0.;
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if ( w_afrac_y(i,j,k)>0.){
                flx_w_arr[1](i,j,k) = 0.25 * w_afrac_y(i,j,k)
                                    * (rho_v(i,j,k) + rho_v(i,j,k-1)) * mf_vy_inv(i,j,0)
                                    * (w(i,j-1,k) + w(i,j,k));
            } else {
                flx_w_arr[1](i,j,k) = 0.;
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if ( w_afrac_z(i,j,k)>0.){
                flx_w_arr[2](i,j,k) = (k==hi_z_face+1) ? omega(i,j,k) * w(i,j,k) : // Not sure for this line
                                    0.25 * (omega(i,j,k) + omega(i,j,k-1)) * (w(i,j,k) + w(i,j,k-1));
            } else {
                flx_w_arr[2](i,j,k) = 0.;
            }
        });

    // Template higher order methods
    } else {

        if (horiz_adv_type == AdvType::Centered_2nd) {
            EBAdvectionSrcForMomVert<CENTERED2>(bxx_grown, bxy_grown, bxz_grown,
                                                rho_u, rho_v, omega, u, v, w,
                                                u_cflag, u_afrac_x, u_afrac_y, u_afrac_z,
                                                v_cflag, v_afrac_x, v_afrac_y, v_afrac_z,
                                                w_cflag, w_afrac_x, w_afrac_y, w_afrac_z,
                                                mf_ux_inv, mf_vx_inv,
                                                mf_uy_inv, mf_vy_inv,
                                                horiz_upw_frac, vert_upw_frac, vert_adv_type,
                                                flx_u_arr, flx_v_arr, flx_w_arr,
                                                lo_z_face, hi_z_face);
        } else if (horiz_adv_type == AdvType::Upwind_3rd) {
            EBAdvectionSrcForMomVert<UPWIND3>(  bxx_grown, bxy_grown, bxz_grown,
                                                rho_u, rho_v, omega, u, v, w,
                                                u_cflag, u_afrac_x, u_afrac_y, u_afrac_z,
                                                v_cflag, v_afrac_x, v_afrac_y, v_afrac_z,
                                                w_cflag, w_afrac_x, w_afrac_y, w_afrac_z,
                                                mf_ux_inv, mf_vx_inv,
                                                mf_uy_inv, mf_vy_inv,
                                                horiz_upw_frac, vert_upw_frac, vert_adv_type,
                                                flx_u_arr, flx_v_arr, flx_w_arr,
                                                lo_z_face, hi_z_face);
        } else if (horiz_adv_type == AdvType::Centered_4th) {
            EBAdvectionSrcForMomVert<CENTERED4>(bxx_grown, bxy_grown, bxz_grown,
                                                rho_u, rho_v, omega, u, v, w,
                                                u_cflag, u_afrac_x, u_afrac_y, u_afrac_z,
                                                v_cflag, v_afrac_x, v_afrac_y, v_afrac_z,
                                                w_cflag, w_afrac_x, w_afrac_y, w_afrac_z,
                                                mf_ux_inv, mf_vx_inv,
                                                mf_uy_inv, mf_vy_inv,
                                                horiz_upw_frac, vert_upw_frac, vert_adv_type,
                                                flx_u_arr, flx_v_arr, flx_w_arr,
                                                lo_z_face, hi_z_face);
        } else if (horiz_adv_type == AdvType::Upwind_5th) {
            EBAdvectionSrcForMomVert<UPWIND5>(  bxx_grown, bxy_grown, bxz_grown,
                                                rho_u, rho_v, omega, u, v, w,
                                                u_cflag, u_afrac_x, u_afrac_y, u_afrac_z,
                                                v_cflag, v_afrac_x, v_afrac_y, v_afrac_z,
                                                w_cflag, w_afrac_x, w_afrac_y, w_afrac_z,
                                                mf_ux_inv, mf_vx_inv,
                                                mf_uy_inv, mf_vy_inv,
                                                horiz_upw_frac, vert_upw_frac, vert_adv_type,
                                                flx_u_arr, flx_v_arr, flx_w_arr,
                                                lo_z_face, hi_z_face);
        } else if (horiz_adv_type == AdvType::Centered_6th) {
            EBAdvectionSrcForMomVert<CENTERED6>(bxx_grown, bxy_grown, bxz_grown,
                                                rho_u, rho_v, omega, u, v, w,
                                                u_cflag, u_afrac_x, u_afrac_y, u_afrac_z,
                                                v_cflag, v_afrac_x, v_afrac_y, v_afrac_z,
                                                w_cflag, w_afrac_x, w_afrac_y, w_afrac_z,
                                                mf_ux_inv, mf_vx_inv,
                                                mf_uy_inv, mf_vy_inv,
                                                horiz_upw_frac, vert_upw_frac, vert_adv_type,
                                                flx_u_arr, flx_v_arr, flx_w_arr,
                                                lo_z_face, hi_z_face);
        } else {
            AMREX_ASSERT_WITH_MESSAGE(false, "Unknown advection scheme!");
        }
    } // horiz_adv_type

    // Update momentum RHS using the fluxes
    if (already_on_centroids) {

        ParallelFor(bxx, bxy, bxz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (u_vfrac(i,j,k)>0.) {
                Real mfsq = mf_ux(i,j,0) * mf_uy(i,j,0);

                Real advectionSrc = ( (flx_u_arr[0](i+1, j  , k  ) - flx_u_arr[0](i, j, k)) * dxInv * mfsq
                                    + (flx_u_arr[1](i  , j+1, k  ) - flx_u_arr[1](i, j, k)) * dyInv * mfsq
                                    + (flx_u_arr[2](i  , j  , k+1) - flx_u_arr[2](i, j, k)) * dzInv ) / u_vfrac(i,j,k);
                rho_u_rhs(i, j, k) = -advectionSrc;
            } else {
                rho_u_rhs(i, j, k) = 0.;
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (v_vfrac(i,j,k)>0.) {
                Real mfsq = mf_vx(i,j,0) * mf_vy(i,j,0);

                Real advectionSrc = ( (flx_v_arr[0](i+1, j  , k  ) - flx_v_arr[0](i, j, k)) * dxInv * mfsq
                                    + (flx_v_arr[1](i  , j+1, k  ) - flx_v_arr[1](i, j, k)) * dyInv * mfsq
                                    + (flx_v_arr[2](i  , j  , k+1) - flx_v_arr[2](i, j, k)) * dzInv ) / v_vfrac(i,j,k);
                rho_v_rhs(i, j, k) = -advectionSrc;
            } else {
                rho_v_rhs(i, j, k) = 0.;
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (w_vfrac(i,j,k)>0.) {
                Real mfsq = mf_mx(i,j,0) * mf_my(i,j,0);

                rho_w_rhs(i, j, k) = - ( (w_afrac_x(i+1, j  , k  ) * flx_w_arr[0](i+1, j  , k  ) - w_afrac_x(i, j, k) * flx_w_arr[0](i, j, k)) * dxInv * mfsq
                                       + (w_afrac_y(i  , j+1, k  ) * flx_w_arr[1](i  , j+1, k  ) - w_afrac_y(i, j, k) * flx_w_arr[1](i, j, k)) * dyInv * mfsq
                                       + (w_afrac_z(i  , j  , k+1) * flx_w_arr[2](i  , j  , k+1) - w_afrac_z(i, j, k) * flx_w_arr[2](i, j, k)) * dzInv ) / w_vfrac(i,j,k);
            } else {
                rho_w_rhs(i, j, k) = 0;
            }
        });

    } else {
        // !already_on_centroids

        Array4<const int> u_mask = physbnd_mask[IntVars::xmom].const_array(mfi);
        Array4<const int> v_mask = physbnd_mask[IntVars::ymom].const_array(mfi);
        Array4<const int> w_mask = physbnd_mask[IntVars::zmom].const_array(mfi);

        ParallelFor(bxx, bxy, bxz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (u_vfrac(i,j,k)>0.) {
                Real mfsq = mf_ux(i,j,0) * mf_uy(i,j,0);

                if (u_cflag(i,j,k).isCovered())
                {
                    rho_u_rhs(i, j, k) = 0.;
                }
                else if (u_cflag(i,j,k).isRegular())
                {
                    rho_u_rhs(i, j, k) = - ( (u_afrac_x(i+1, j  , k  ) * flx_u_arr[0](i+1, j  , k  ) - u_afrac_x(i, j, k) * flx_u_arr[0](i, j, k)) * dxInv * mfsq
                                           + (u_afrac_y(i  , j+1, k  ) * flx_u_arr[1](i  , j+1, k  ) - u_afrac_y(i, j, k) * flx_u_arr[1](i, j, k)) * dyInv * mfsq
                                           + (u_afrac_z(i  , j  , k+1) * flx_u_arr[2](i  , j  , k+1) - u_afrac_z(i, j, k) * flx_u_arr[2](i, j, k)) * dzInv ) / u_vfrac(i,j,k);
                }
                else
                {
                    // Bilinear interpolation
                    Real fxm = flx_u_arr[0](i,j,k);
                    if (u_afrac_x(i,j,k) != Real(0.0) && u_afrac_x(i,j,k) != Real(1.0)) {
                        int jj = j + static_cast<int>(std::copysign(Real(1.0), u_fcx(i,j,k,0)));
                        int kk = k + static_cast<int>(std::copysign(Real(1.0), u_fcx(i,j,k,1)));
                        Real fracy = (u_mask(i-1,jj,k) || u_mask(i,jj,k)) ? std::abs(u_fcx(i,j,k,0)) : Real(0.0);
                        Real fracz = (u_mask(i-1,j,kk) || u_mask(i,j,kk)) ? std::abs(u_fcx(i,j,k,1)) : Real(0.0);
                        fxm = (Real(1.0)-fracy)*(Real(1.0)-fracz)*fxm
                            +      fracy *(Real(1.0)-fracz)*flx_u_arr[0](i,jj,k )
                            +      fracz *(Real(1.0)-fracy)*flx_u_arr[0](i,j ,kk)
                            +      fracy *           fracz *flx_u_arr[0](i,jj,kk);
                    }

                    Real fxp = flx_u_arr[0](i+1,j,k);
                    if (u_afrac_x(i+1,j,k) != Real(0.0) && u_afrac_x(i+1,j,k) != Real(1.0)) {
                        int jj = j + static_cast<int>(std::copysign(Real(1.0),u_fcx(i+1,j,k,0)));
                        int kk = k + static_cast<int>(std::copysign(Real(1.0),u_fcx(i+1,j,k,1)));
                        Real fracy = (u_mask(i,jj,k) || u_mask(i+1,jj,k)) ? std::abs(u_fcx(i+1,j,k,0)) : Real(0.0);
                        Real fracz = (u_mask(i,j,kk) || u_mask(i+1,j,kk)) ? std::abs(u_fcx(i+1,j,k,1)) : Real(0.0);
                        fxp = (Real(1.0)-fracy)*(Real(1.0)-fracz)*fxp
                            +      fracy *(Real(1.0)-fracz)*flx_u_arr[0](i+1,jj,k )
                            +      fracz *(Real(1.0)-fracy)*flx_u_arr[0](i+1,j ,kk)
                            +      fracy *     fracz *flx_u_arr[0](i+1,jj,kk);
                    }

                    Real fym = flx_u_arr[1](i,j,k);
                    if (u_afrac_y(i,j,k) != Real(0.0) && u_afrac_y(i,j,k) != Real(1.0)) {
                        int ii = i + static_cast<int>(std::copysign(Real(1.0),u_fcy(i,j,k,0)));
                        int kk = k + static_cast<int>(std::copysign(Real(1.0),u_fcy(i,j,k,1)));
                        Real fracx = (u_mask(ii,j-1,k) || u_mask(ii,j,k)) ? std::abs(u_fcy(i,j,k,0)) : Real(0.0);
                        Real fracz = (u_mask(i,j-1,kk) || u_mask(i,j,kk)) ? std::abs(u_fcy(i,j,k,1)) : Real(0.0);
                        fym = (Real(1.0)-fracx)*(Real(1.0)-fracz)*fym
                            +      fracx *(Real(1.0)-fracz)*flx_u_arr[1](ii,j,k )
                            +      fracz *(Real(1.0)-fracx)*flx_u_arr[1](i ,j,kk)
                            +      fracx *     fracz *flx_u_arr[1](ii,j,kk);
                    }

                    Real fyp = flx_u_arr[1](i,j+1,k);
                    if (u_afrac_y(i,j+1,k) != Real(0.0) && u_afrac_y(i,j+1,k) != Real(1.0)) {
                        int ii = i + static_cast<int>(std::copysign(Real(1.0),u_fcy(i,j+1,k,0)));
                        int kk = k + static_cast<int>(std::copysign(Real(1.0),u_fcy(i,j+1,k,1)));
                        Real fracx = (u_mask(ii,j,k) || u_mask(ii,j+1,k)) ? std::abs(u_fcy(i,j+1,k,0)) : Real(0.0);
                        Real fracz = (u_mask(i,j,kk) || u_mask(i,j+1,kk)) ? std::abs(u_fcy(i,j+1,k,1)) : Real(0.0);
                        fyp = (Real(1.0)-fracx)*(Real(1.0)-fracz)*fyp
                            +      fracx *(Real(1.0)-fracz)*flx_u_arr[1](ii,j+1,k )
                            +      fracz *(Real(1.0)-fracx)*flx_u_arr[1](i ,j+1,kk)
                            +      fracx *     fracz *flx_u_arr[1](ii,j+1,kk);
                    }

                    Real fzm = flx_u_arr[2](i,j,k);
                    if (u_afrac_z(i,j,k) != Real(0.0) && u_afrac_z(i,j,k) != Real(1.0)) {
                        int ii = i + static_cast<int>(std::copysign(Real(1.0),u_fcz(i,j,k,0)));
                        int jj = j + static_cast<int>(std::copysign(Real(1.0),u_fcz(i,j,k,1)));
                        Real fracx = (u_mask(ii,j,k-1) || u_mask(ii,j,k)) ? std::abs(u_fcz(i,j,k,0)) : Real(0.0);
                        Real fracy = (u_mask(i,jj,k-1) || u_mask(i,jj,k)) ? std::abs(u_fcz(i,j,k,1)) : Real(0.0);
                        fzm = (Real(1.0)-fracx)*(Real(1.0)-fracy)*fzm
                            +      fracx *(Real(1.0)-fracy)*flx_u_arr[2](ii,j ,k)
                            +      fracy *(Real(1.0)-fracx)*flx_u_arr[2](i ,jj,k)
                            +      fracx *     fracy *flx_u_arr[2](ii,jj,k);
                    }

                    Real fzp = flx_u_arr[2](i,j,k+1);
                    if (u_afrac_z(i,j,k+1) != Real(0.0) && u_afrac_z(i,j,k+1) != Real(1.0)) {
                        int ii = i + static_cast<int>(std::copysign(Real(1.0),u_fcz(i,j,k+1,0)));
                        int jj = j + static_cast<int>(std::copysign(Real(1.0),u_fcz(i,j,k+1,1)));
                        Real fracx = (u_mask(ii,j,k) || u_mask(ii,j,k+1)) ? std::abs(u_fcz(i,j,k+1,0)) : Real(0.0);
                        Real fracy = (u_mask(i,jj,k) || u_mask(i,jj,k+1)) ? std::abs(u_fcz(i,j,k+1,1)) : Real(0.0);
                        fzp = (Real(1.0)-fracx)*(Real(1.0)-fracy)*fzp
                            +      fracx *(Real(1.0)-fracy)*flx_u_arr[2](ii,j ,k+1)
                            +      fracy *(Real(1.0)-fracx)*flx_u_arr[2](i ,jj,k+1)
                            +      fracx *     fracy *flx_u_arr[2](ii,jj,k+1);
                    }

                    rho_u_rhs(i, j, k) = - ( (u_afrac_x(i+1, j  , k  ) * fxp - u_afrac_x(i, j, k) * fxm) * dxInv * mfsq
                                           + (u_afrac_y(i  , j+1, k  ) * fyp - u_afrac_y(i, j, k) * fym) * dyInv * mfsq
                                           + (u_afrac_z(i  , j  , k+1) * fzp - u_afrac_z(i, j, k) * fzm) * dzInv ) / u_vfrac(i,j,k);
                }

            } else {
                rho_u_rhs(i, j, k) = 0.;
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (v_vfrac(i,j,k)>0.) {
                Real mfsq = mf_vx(i,j,0) * mf_vy(i,j,0);

                if (v_cflag(i,j,k).isCovered())
                {
                    rho_v_rhs(i, j, k) = 0.;
                }
                else if (v_cflag(i,j,k).isRegular())
                {
                    rho_v_rhs(i, j, k) = - ( (v_afrac_x(i+1, j  , k  ) * flx_v_arr[0](i+1, j  , k  ) - v_afrac_x(i, j, k) * flx_v_arr[0](i, j, k)) * dxInv * mfsq
                                           + (v_afrac_y(i  , j+1, k  ) * flx_v_arr[1](i  , j+1, k  ) - v_afrac_y(i, j, k) * flx_v_arr[1](i, j, k)) * dyInv * mfsq
                                           + (v_afrac_z(i  , j  , k+1) * flx_v_arr[2](i  , j  , k+1) - v_afrac_z(i, j, k) * flx_v_arr[2](i, j, k)) * dzInv ) / v_vfrac(i,j,k);
                }
                else
                {
                    // Bilinear interpolation
                    Real fxm = flx_v_arr[0](i,j,k);
                    if (v_afrac_x(i,j,k) != Real(0.0) && v_afrac_x(i,j,k) != Real(1.0)) {
                        int jj = j + static_cast<int>(std::copysign(Real(1.0), v_fcx(i,j,k,0)));
                        int kk = k + static_cast<int>(std::copysign(Real(1.0), v_fcx(i,j,k,1)));
                        Real fracy = (v_mask(i-1,jj,k) || v_mask(i,jj,k)) ? std::abs(v_fcx(i,j,k,0)) : Real(0.0);
                        Real fracz = (v_mask(i-1,j,kk) || v_mask(i,j,kk)) ? std::abs(v_fcx(i,j,k,1)) : Real(0.0);
                        fxm = (Real(1.0)-fracy)*(Real(1.0)-fracz)*fxm
                            +      fracy *(Real(1.0)-fracz)*flx_v_arr[0](i,jj,k )
                            +      fracz *(Real(1.0)-fracy)*flx_v_arr[0](i,j ,kk)
                            +      fracy *     fracz *flx_v_arr[0](i,jj,kk);
                    }

                    Real fxp = flx_v_arr[0](i+1,j,k);
                    if (v_afrac_x(i+1,j,k) != Real(0.0) && v_afrac_x(i+1,j,k) != Real(1.0)) {
                        int jj = j + static_cast<int>(std::copysign(Real(1.0),v_fcx(i+1,j,k,0)));
                        int kk = k + static_cast<int>(std::copysign(Real(1.0),v_fcx(i+1,j,k,1)));
                        Real fracy = (v_mask(i,jj,k) || v_mask(i+1,jj,k)) ? std::abs(v_fcx(i+1,j,k,0)) : Real(0.0);
                        Real fracz = (v_mask(i,j,kk) || v_mask(i+1,j,kk)) ? std::abs(v_fcx(i+1,j,k,1)) : Real(0.0);
                        fxp = (Real(1.0)-fracy)*(Real(1.0)-fracz)*fxp
                            +      fracy *(Real(1.0)-fracz)*flx_v_arr[0](i+1,jj,k )
                            +      fracz *(Real(1.0)-fracy)*flx_v_arr[0](i+1,j ,kk)
                            +      fracy *     fracz *flx_v_arr[0](i+1,jj,kk);
                    }

                    Real fym = flx_v_arr[1](i,j,k);
                    if (v_afrac_y(i,j,k) != Real(0.0) && v_afrac_y(i,j,k) != Real(1.0)) {
                        int ii = i + static_cast<int>(std::copysign(Real(1.0),v_fcy(i,j,k,0)));
                        int kk = k + static_cast<int>(std::copysign(Real(1.0),v_fcy(i,j,k,1)));
                        Real fracx = (v_mask(ii,j-1,k) || v_mask(ii,j,k)) ? std::abs(v_fcy(i,j,k,0)) : Real(0.0);
                        Real fracz = (v_mask(i,j-1,kk) || v_mask(i,j,kk)) ? std::abs(v_fcy(i,j,k,1)) : Real(0.0);
                        fym = (Real(1.0)-fracx)*(Real(1.0)-fracz)*fym
                            +      fracx *(Real(1.0)-fracz)*flx_v_arr[1](ii,j,k )
                            +      fracz *(Real(1.0)-fracx)*flx_v_arr[1](i ,j,kk)
                            +      fracx *           fracz *flx_v_arr[1](ii,j,kk);
                    }

                    Real fyp = flx_v_arr[1](i,j+1,k);
                    if (v_afrac_y(i,j+1,k) != Real(0.0) && v_afrac_y(i,j+1,k) != Real(1.0)) {
                        int ii = i + static_cast<int>(std::copysign(Real(1.0),v_fcy(i,j+1,k,0)));
                        int kk = k + static_cast<int>(std::copysign(Real(1.0),v_fcy(i,j+1,k,1)));
                        Real fracx = (v_mask(ii,j,k) || v_mask(ii,j+1,k)) ? std::abs(v_fcy(i,j+1,k,0)) : Real(0.0);
                        Real fracz = (v_mask(i,j,kk) || v_mask(i,j+1,kk)) ? std::abs(v_fcy(i,j+1,k,1)) : Real(0.0);
                        fyp = (Real(1.0)-fracx)*(Real(1.0)-fracz)*fyp
                            +      fracx *(Real(1.0)-fracz)*flx_v_arr[1](ii,j+1,k )
                            +      fracz *(Real(1.0)-fracx)*flx_v_arr[1](i ,j+1,kk)
                            +      fracx *     fracz *flx_v_arr[1](ii,j+1,kk);
                    }

                    Real fzm = flx_v_arr[2](i,j,k);
                    if (v_afrac_z(i,j,k) != Real(0.0) && v_afrac_z(i,j,k) != Real(1.0)) {
                        int ii = i + static_cast<int>(std::copysign(Real(1.0),v_fcz(i,j,k,0)));
                        int jj = j + static_cast<int>(std::copysign(Real(1.0),v_fcz(i,j,k,1)));
                        Real fracx = (v_mask(ii,j,k-1) || v_mask(ii,j,k)) ? std::abs(v_fcz(i,j,k,0)) : Real(0.0);
                        Real fracy = (v_mask(i,jj,k-1) || v_mask(i,jj,k)) ? std::abs(v_fcz(i,j,k,1)) : Real(0.0);
                        fzm = (Real(1.0)-fracx)*(Real(1.0)-fracy)*fzm
                            +      fracx *(Real(1.0)-fracy)*flx_v_arr[2](ii,j ,k)
                            +      fracy *(Real(1.0)-fracx)*flx_v_arr[2](i ,jj,k)
                            +      fracx *     fracy *flx_v_arr[2](ii,jj,k);
                    }

                    Real fzp = flx_v_arr[2](i,j,k+1);
                    if (v_afrac_z(i,j,k+1) != Real(0.0) && v_afrac_z(i,j,k+1) != Real(1.0)) {
                        int ii = i + static_cast<int>(std::copysign(Real(1.0),v_fcz(i,j,k+1,0)));
                        int jj = j + static_cast<int>(std::copysign(Real(1.0),v_fcz(i,j,k+1,1)));
                        Real fracx = (v_mask(ii,j,k) || v_mask(ii,j,k+1)) ? std::abs(v_fcz(i,j,k+1,0)) : Real(0.0);
                        Real fracy = (v_mask(i,jj,k) || v_mask(i,jj,k+1)) ? std::abs(v_fcz(i,j,k+1,1)) : Real(0.0);
                        fzp = (Real(1.0)-fracx)*(Real(1.0)-fracy)*fzp
                            +      fracx *(Real(1.0)-fracy)*flx_v_arr[2](ii,j ,k+1)
                            +      fracy *(Real(1.0)-fracx)*flx_v_arr[2](i ,jj,k+1)
                            +      fracx *     fracy *flx_v_arr[2](ii,jj,k+1);
                    }

                    rho_v_rhs(i, j, k) = - ( (v_afrac_x(i+1, j  , k  ) * fxp - v_afrac_x(i, j, k) * fxm) * dxInv * mfsq
                                           + (v_afrac_y(i  , j+1, k  ) * fyp - v_afrac_y(i, j, k) * fym) * dyInv * mfsq
                                           + (v_afrac_z(i  , j  , k+1) * fzp - v_afrac_z(i, j, k) * fzm) * dzInv ) / v_vfrac(i,j,k);
                }

            } else {
                rho_v_rhs(i, j, k) = 0.;
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (w_vfrac(i,j,k)>0.) {
                Real mfsq = mf_mx(i,j,0) * mf_my(i,j,0);

                if (w_cflag(i,j,k).isCovered())
                {
                    rho_w_rhs(i, j, k) = 0.;
                }
                else if (w_cflag(i,j,k).isRegular())
                {
                    rho_w_rhs(i, j, k) = - ( (w_afrac_x(i+1, j  , k  ) * flx_w_arr[0](i+1, j  , k  ) - w_afrac_x(i, j, k) * flx_w_arr[0](i, j, k)) * dxInv * mfsq
                                           + (w_afrac_y(i  , j+1, k  ) * flx_w_arr[1](i  , j+1, k  ) - w_afrac_y(i, j, k) * flx_w_arr[1](i, j, k)) * dyInv * mfsq
                                           + (w_afrac_z(i  , j  , k+1) * flx_w_arr[2](i  , j  , k+1) - w_afrac_z(i, j, k) * flx_w_arr[2](i, j, k)) * dzInv ) / w_vfrac(i,j,k);
                }
                else
                {
                    // Bilinear interpolation
                    Real fxm = flx_w_arr[0](i,j,k);
                    if (w_afrac_x(i,j,k) != Real(0.0) && w_afrac_x(i,j,k) != Real(1.0)) {
                        int jj = j + static_cast<int>(std::copysign(Real(1.0), w_fcx(i,j,k,0)));
                        int kk = k + static_cast<int>(std::copysign(Real(1.0), w_fcx(i,j,k,1)));
                        Real fracy = (w_mask(i-1,jj,k) || w_mask(i,jj,k)) ? std::abs(w_fcx(i,j,k,0)) : Real(0.0);
                        Real fracz = (w_mask(i-1,j,kk) || w_mask(i,j,kk)) ? std::abs(w_fcx(i,j,k,1)) : Real(0.0);
                        fxm = (Real(1.0)-fracy)*(Real(1.0)-fracz)*fxm
                            +      fracy *(Real(1.0)-fracz)*flx_w_arr[0](i,jj,k )
                            +      fracz *(Real(1.0)-fracy)*flx_w_arr[0](i,j ,kk)
                            +      fracy *     fracz *flx_w_arr[0](i,jj,kk);
                    }

                    Real fxp = flx_w_arr[0](i+1,j,k);
                    if (w_afrac_x(i+1,j,k) != Real(0.0) && w_afrac_x(i+1,j,k) != Real(1.0)) {
                        int jj = j + static_cast<int>(std::copysign(Real(1.0),w_fcx(i+1,j,k,0)));
                        int kk = k + static_cast<int>(std::copysign(Real(1.0),w_fcx(i+1,j,k,1)));
                        Real fracy = (w_mask(i,jj,k) || w_mask(i+1,jj,k)) ? std::abs(w_fcx(i+1,j,k,0)) : Real(0.0);
                        Real fracz = (w_mask(i,j,kk) || w_mask(i+1,j,kk)) ? std::abs(w_fcx(i+1,j,k,1)) : Real(0.0);
                        fxp = (Real(1.0)-fracy)*(Real(1.0)-fracz)*fxp
                            +      fracy *(Real(1.0)-fracz)*flx_w_arr[0](i+1,jj,k )
                            +      fracz *(Real(1.0)-fracy)*flx_w_arr[0](i+1,j ,kk)
                            +      fracy *     fracz *flx_w_arr[0](i+1,jj,kk);
                    }

                    Real fym = flx_w_arr[1](i,j,k);
                    if (w_afrac_y(i,j,k) != Real(0.0) && w_afrac_y(i,j,k) != Real(1.0)) {
                        int ii = i + static_cast<int>(std::copysign(Real(1.0),w_fcy(i,j,k,0)));
                        int kk = k + static_cast<int>(std::copysign(Real(1.0),w_fcy(i,j,k,1)));
                        Real fracx = (w_mask(ii,j-1,k) || w_mask(ii,j,k)) ? std::abs(w_fcy(i,j,k,0)) : Real(0.0);
                        Real fracz = (w_mask(i,j-1,kk) || w_mask(i,j,kk)) ? std::abs(w_fcy(i,j,k,1)) : Real(0.0);
                        fym = (Real(1.0)-fracx)*(Real(1.0)-fracz)*fym
                            +      fracx *(Real(1.0)-fracz)*flx_w_arr[1](ii,j,k )
                            +      fracz *(Real(1.0)-fracx)*flx_w_arr[1](i ,j,kk)
                            +      fracx *     fracz *flx_w_arr[1](ii,j,kk);
                    }

                    Real fyp = flx_w_arr[1](i,j+1,k);
                    if (w_afrac_y(i,j+1,k) != Real(0.0) && w_afrac_y(i,j+1,k) != Real(1.0)) {
                        int ii = i + static_cast<int>(std::copysign(Real(1.0),w_fcy(i,j+1,k,0)));
                        int kk = k + static_cast<int>(std::copysign(Real(1.0),w_fcy(i,j+1,k,1)));
                        Real fracx = (w_mask(ii,j,k) || w_mask(ii,j+1,k)) ? std::abs(w_fcy(i,j+1,k,0)) : Real(0.0);
                        Real fracz = (w_mask(i,j,kk) || w_mask(i,j+1,kk)) ? std::abs(w_fcy(i,j+1,k,1)) : Real(0.0);
                        fyp = (Real(1.0)-fracx)*(Real(1.0)-fracz)*fyp
                            +      fracx *(Real(1.0)-fracz)*flx_w_arr[1](ii,j+1,k )
                            +      fracz *(Real(1.0)-fracx)*flx_w_arr[1](i ,j+1,kk)
                            +      fracx *     fracz *flx_w_arr[1](ii,j+1,kk);
                    }

                    Real fzm = flx_w_arr[2](i,j,k);
                    if (w_afrac_z(i,j,k) != Real(0.0) && w_afrac_z(i,j,k) != Real(1.0)) {
                        int ii = i + static_cast<int>(std::copysign(Real(1.0),w_fcz(i,j,k,0)));
                        int jj = j + static_cast<int>(std::copysign(Real(1.0),w_fcz(i,j,k,1)));
                        Real fracx = (w_mask(ii,j,k-1) || w_mask(ii,j,k)) ? std::abs(w_fcz(i,j,k,0)) : Real(0.0);
                        Real fracy = (w_mask(i,jj,k-1) || w_mask(i,jj,k)) ? std::abs(w_fcz(i,j,k,1)) : Real(0.0);
                        fzm = (Real(1.0)-fracx)*(Real(1.0)-fracy)*fzm
                            +      fracx *(Real(1.0)-fracy)*flx_w_arr[2](ii,j ,k)
                            +      fracy *(Real(1.0)-fracx)*flx_w_arr[2](i ,jj,k)
                            +      fracx *     fracy *flx_w_arr[2](ii,jj,k);
                    }

                    Real fzp = flx_w_arr[2](i,j,k+1);
                    if (w_afrac_z(i,j,k+1) != Real(0.0) && w_afrac_z(i,j,k+1) != Real(1.0)) {
                        int ii = i + static_cast<int>(std::copysign(Real(1.0),w_fcz(i,j,k+1,0)));
                        int jj = j + static_cast<int>(std::copysign(Real(1.0),w_fcz(i,j,k+1,1)));
                        Real fracx = (w_mask(ii,j,k) || w_mask(ii,j,k+1)) ? std::abs(w_fcz(i,j,k+1,0)) : Real(0.0);
                        Real fracy = (w_mask(i,jj,k) || w_mask(i,jj,k+1)) ? std::abs(w_fcz(i,j,k+1,1)) : Real(0.0);
                        fzp = (Real(1.0)-fracx)*(Real(1.0)-fracy)*fzp
                            +      fracx *(Real(1.0)-fracy)*flx_w_arr[2](ii,j ,k+1)
                            +      fracy *(Real(1.0)-fracx)*flx_w_arr[2](i ,jj,k+1)
                            +      fracx *     fracy *flx_w_arr[2](ii,jj,k+1);
                    }

                    rho_w_rhs(i, j, k) = - ( (w_afrac_x(i+1, j  , k  ) * fxp - w_afrac_x(i, j, k) * fxm) * dxInv * mfsq
                                           + (w_afrac_y(i  , j+1, k  ) * fyp - w_afrac_y(i, j, k) * fym) * dyInv * mfsq
                                           + (w_afrac_z(i  , j  , k+1) * fzp - w_afrac_z(i, j, k) * fzm) * dzInv ) / w_vfrac(i,j,k);
                }

            } else {
                rho_w_rhs(i, j, k) = 0;
            }
        });

    }

}
