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
                        const Array4<const Real>& mf_m,
                        const Array4<const Real>& mf_u,
                        const Array4<const Real>& mf_v,
                        const AdvType horiz_adv_type,
                        const AdvType vert_adv_type,
                        const Real horiz_upw_frac,
                        const Real vert_upw_frac,
                        const eb_& ebfact,
                        const GpuArray<const Array4<Real>, AMREX_SPACEDIM>& flx_u_arr,
                        const GpuArray<const Array4<Real>, AMREX_SPACEDIM>& flx_v_arr,
                        const GpuArray<const Array4<Real>, AMREX_SPACEDIM>& flx_w_arr,
                        const int lo_z_face, const int hi_z_face,
                        const Box& /*domain*/)
{
    BL_PROFILE_VAR("AdvectionSrcForMom_EB", AdvectionSrcForMom_EB);

    AMREX_ALWAYS_ASSERT(bxz.smallEnd(2) > 0);

    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

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

    // EB u-factory
    Array4<const EBCellFlag> u_cflag   = (ebfact.get_u_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();
    Array4<const Real      > u_vfrac   = (ebfact.get_u_const_factory())->getVolFrac().const_array(mfi);
    Array4<const Real      > u_afrac_x = (ebfact.get_u_const_factory())->getAreaFrac()[0]->const_array(mfi);
    Array4<const Real      > u_afrac_y = (ebfact.get_u_const_factory())->getAreaFrac()[1]->const_array(mfi);
    Array4<const Real      > u_afrac_z = (ebfact.get_u_const_factory())->getAreaFrac()[2]->const_array(mfi);

    // EB v-factory
    Array4<const EBCellFlag> v_cflag   = (ebfact.get_v_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();
    Array4<const Real      > v_vfrac   = (ebfact.get_v_const_factory())->getVolFrac().const_array(mfi);
    Array4<const Real      > v_afrac_x = (ebfact.get_v_const_factory())->getAreaFrac()[0]->const_array(mfi);
    Array4<const Real      > v_afrac_y = (ebfact.get_v_const_factory())->getAreaFrac()[1]->const_array(mfi);
    Array4<const Real      > v_afrac_z = (ebfact.get_v_const_factory())->getAreaFrac()[2]->const_array(mfi);

    // EB w-factory
    Array4<const EBCellFlag> w_cflag   = (ebfact.get_w_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();
    Array4<const Real      > w_vfrac   = (ebfact.get_w_const_factory())->getVolFrac().const_array(mfi);
    Array4<const Real      > w_afrac_x = (ebfact.get_w_const_factory())->getAreaFrac()[0]->const_array(mfi);
    Array4<const Real      > w_afrac_y = (ebfact.get_w_const_factory())->getAreaFrac()[1]->const_array(mfi);
    Array4<const Real      > w_afrac_z = (ebfact.get_w_const_factory())->getAreaFrac()[2]->const_array(mfi);

    // Inline with 2nd order for efficiency
    if (horiz_adv_type == AdvType::Centered_2nd && vert_adv_type == AdvType::Centered_2nd)
    {
        // Fluxes for x-momentum
        ParallelFor(bxx_grown[0], bxx_grown[1], bxx_grown[2],
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if ( u_afrac_x(i,j,k)>0.){
                flx_u_arr[0](i,j,k) = 0.25 * u_afrac_x(i,j,k)
                                    * (rho_u(i,j,k) * mf_u_inv(i,j,0) + rho_u(i-1,j,k) * mf_u_inv(i-1,j,0))
                                    * (u(i-1,j,k) + u(i,j,k));
            } else {
                flx_u_arr[0](i,j,k) = 0.;
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if ( u_afrac_y(i,j,k)>0.){
                flx_u_arr[1](i,j,k) = 0.25 * u_afrac_y(i,j,k)
                                    * (rho_v(i,j,k) * mf_v_inv(i,j,0) + rho_v(i-1,j,k) * mf_v_inv(i-1,j,0))
                                    * (u(i,j-1,k) + u(i,j,k));
            } else {
                flx_u_arr[1](i,j,k) = 0.;
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if ( u_afrac_z(i,j,k)>0.){
                flx_u_arr[2](i,j,k) = 0.25 * u_afrac_z(i,j,k)
                                    * (omega(i,j,k) + omega(i-1,j,k)) * (u(i,j,k-1) + u(i,j,k));
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
                                    * (rho_u(i,j,k) * mf_u_inv(i,j,0) + rho_u(i,j-1,k) * mf_u_inv(i,j-1,0))
                                    * (v(i-1,j,k) + v(i,j,k));
            } else {
                flx_v_arr[0](i,j,k) = 0.;
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if ( v_afrac_y(i,j,k)>0.){
                flx_v_arr[1](i,j,k) = 0.25 * v_afrac_y(i,j,k)
                                    * (rho_v(i,j,k) * mf_v_inv(i,j,0) + rho_v(i,j-1,k) * mf_v_inv(i,j-1,0))
                                    * (v(i,j-1,k) + v(i,j,k));
            } else {
                flx_v_arr[1](i,j,k) = 0.;
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if ( v_afrac_z(i,j,k)>0.){
                flx_v_arr[2](i,j,k) = 0.25 * v_afrac_z(i,j,k  )
                                    * (omega(i,j,k) + omega(i,j-1,k)) * (v(i,j,k-1) + v(i,j,k));
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
                                    * (rho_u(i,j,k) + rho_u(i,j, k-1)) * mf_u_inv(i,j,0)
                                    * (w(i-1,j,k) + w(i,j,k));
            } else {
                flx_w_arr[0](i,j,k) = 0.;
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if ( w_afrac_y(i,j,k)>0.){
                flx_w_arr[1](i,j,k) = 0.25 * w_afrac_y(i,j,k)
                                    * (rho_v(i,j,k) + rho_v(i,j,k-1)) * mf_v_inv(i,j,0)
                                    * (w(i,j-1,k) + w(i,j,k));
            } else {
                flx_w_arr[1](i,j,k) = 0.;
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if ( w_afrac_z(i,j,k)>0.){
                flx_w_arr[2](i,j,k) = (k==hi_z_face+1) ? omega(i,j,k) * w(i,j,k) : // Not sure for this line
                                    0.25 * w_afrac_z(i,j,k)
                                    * (omega(i,j,k) + omega(i,j,k-1)) * (w(i,j,k) + w(i,j,k-1));
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
                                                mf_u_inv, mf_v_inv,
                                                horiz_upw_frac, vert_upw_frac, vert_adv_type,
                                                flx_u_arr, flx_v_arr, flx_w_arr,
                                                lo_z_face, hi_z_face);
        } else if (horiz_adv_type == AdvType::Upwind_3rd) {
            EBAdvectionSrcForMomVert<UPWIND3>(  bxx_grown, bxy_grown, bxz_grown,
                                                rho_u, rho_v, omega, u, v, w,
                                                u_cflag, u_afrac_x, u_afrac_y, u_afrac_z,
                                                v_cflag, v_afrac_x, v_afrac_y, v_afrac_z,
                                                w_cflag, w_afrac_x, w_afrac_y, w_afrac_z,
                                                mf_u_inv, mf_v_inv,
                                                horiz_upw_frac, vert_upw_frac, vert_adv_type,
                                                flx_u_arr, flx_v_arr, flx_w_arr,
                                                lo_z_face, hi_z_face);
        } else if (horiz_adv_type == AdvType::Centered_4th) {
            EBAdvectionSrcForMomVert<CENTERED4>(bxx_grown, bxy_grown, bxz_grown,
                                                rho_u, rho_v, omega, u, v, w,
                                                u_cflag, u_afrac_x, u_afrac_y, u_afrac_z,
                                                v_cflag, v_afrac_x, v_afrac_y, v_afrac_z,
                                                w_cflag, w_afrac_x, w_afrac_y, w_afrac_z,
                                                mf_u_inv, mf_v_inv,
                                                horiz_upw_frac, vert_upw_frac, vert_adv_type,
                                                flx_u_arr, flx_v_arr, flx_w_arr,
                                                lo_z_face, hi_z_face);
        } else if (horiz_adv_type == AdvType::Upwind_5th) {
            EBAdvectionSrcForMomVert<UPWIND5>(  bxx_grown, bxy_grown, bxz_grown,
                                                rho_u, rho_v, omega, u, v, w,
                                                u_cflag, u_afrac_x, u_afrac_y, u_afrac_z,
                                                v_cflag, v_afrac_x, v_afrac_y, v_afrac_z,
                                                w_cflag, w_afrac_x, w_afrac_y, w_afrac_z,
                                                mf_u_inv, mf_v_inv,
                                                horiz_upw_frac, vert_upw_frac, vert_adv_type,
                                                flx_u_arr, flx_v_arr, flx_w_arr,
                                                lo_z_face, hi_z_face);
        } else if (horiz_adv_type == AdvType::Centered_6th) {
            EBAdvectionSrcForMomVert<CENTERED6>(bxx_grown, bxy_grown, bxz_grown,
                                                rho_u, rho_v, omega, u, v, w,
                                                u_cflag, u_afrac_x, u_afrac_y, u_afrac_z,
                                                v_cflag, v_afrac_x, v_afrac_y, v_afrac_z,
                                                w_cflag, w_afrac_x, w_afrac_y, w_afrac_z,
                                                mf_u_inv, mf_v_inv,
                                                horiz_upw_frac, vert_upw_frac, vert_adv_type,
                                                flx_u_arr, flx_v_arr, flx_w_arr,
                                                lo_z_face, hi_z_face);
        } else {
            AMREX_ASSERT_WITH_MESSAGE(false, "Unknown advection scheme!");
        }
    } // horiz_adv_type

    // Update momentum RHS using the fluxes

    ParallelFor(bxx, bxy, bxz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (u_vfrac(i,j,k)>0.) {
            Real mfsq = mf_u(i,j,0) * mf_u(i,j,0);

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
            Real mfsq = mf_v(i,j,0) * mf_v(i,j,0);

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
            Real mfsq = mf_m(i,j,0) * mf_m(i,j,0);

            Real advectionSrc = ( (flx_w_arr[0](i+1, j  , k  ) - flx_w_arr[0](i, j, k)) * dxInv * mfsq
                                + (flx_w_arr[1](i  , j+1, k  ) - flx_w_arr[1](i, j, k)) * dyInv * mfsq
                                + (flx_w_arr[2](i  , j  , k+1) - flx_w_arr[2](i, j, k)) * dzInv ) / w_vfrac(i,j,k);
            rho_w_rhs(i, j, k) = -advectionSrc;
        } else {
            rho_w_rhs(i, j, k) = 0;
        }
    });

}