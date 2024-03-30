#include <IndexDefines.H>
#include <TerrainMetrics.H>
#include <Advection.H>
#include <AdvectionSrcForScalars.H>

using namespace amrex;

/**
 * Function for computing the advective tendency for the update equations for rho and (rho theta)
 * This routine has explicit expressions for all cases (terrain or not) when
 * the horizontal and vertial spatial orders are <= 2, and calls more specialized
 * functions when either (or both) spatial order(s) is greater than 2.
 *
 * @param[in] bx box over which the scalars are updated
 * @param[out] advectionSrc tendency for the scalar update equation
 * @param[in] rho_u x-component of momentum
 * @param[in] rho_v y-component of momentum
 * @param[in] Omega component of momentum normal to the z-coordinate surface
 * @param[out] avg_xmom x-component of time-averaged momentum defined in this routine
 * @param[out] avg_ymom y-component of time-averaged momentum defined in this routine
 * @param[out] avg_zmom z-component of time-averaged momentum defined in this routine
 * @param[in] z_nd height coordinate at nodes
 * @param[in] detJ Jacobian of the metric transformation (= 1 if use_terrain is false)
 * @param[in] cellSizeInv inverse of the mesh spacing
 * @param[in] mf_m map factor at cell centers
 * @param[in] mf_u map factor at x-faces
 * @param[in] mf_v map factor at y-faces
 * @param[in] horiz_adv_type advection scheme to be used in horiz. directions for dry scalars
 * @param[in] vert_adv_type advection scheme to be used in vert. directions for dry scalars
 * @param[in] horiz_upw_frac upwinding fraction to be used in horiz. directions for dry scalars (for Blended schemes only)
 * @param[in] vert_upw_frac upwinding fraction to be used in vert. directions for dry scalars (for Blended schemes only)
 * @param[in] use_terrain if true, use the terrain-aware derivatives (with metric terms)
 */

void
AdvectionSrcForRho (const Box& bx,
                    const Array4<Real>& advectionSrc,
                    const Array4<const Real>& rho_u,
                    const Array4<const Real>& rho_v,
                    const Array4<const Real>& Omega,
                    const Array4<      Real>& avg_xmom,
                    const Array4<      Real>& avg_ymom,
                    const Array4<      Real>& avg_zmom,
                    const Array4<const Real>& z_nd,
                    const Array4<const Real>& detJ,
                    const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                    const Array4<const Real>& mf_m,
                    const Array4<const Real>& mf_u,
                    const Array4<const Real>& mf_v,
#ifdef ERF_USE_EB
                    const Array4<const Real>& ax_arr,
                    const Array4<const Real>& ay_arr,
                    const Array4<const Real>& az_arr,
                    const Array4<const Real>& vf_arr,
#endif
                    const bool use_terrain,
                    const GpuArray<const Array4<Real>, AMREX_SPACEDIM>& flx_arr)
{
    BL_PROFILE_VAR("AdvectionSrcForRhoAndTheta", AdvectionSrcForRhoAndTheta);
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    const Box xbx = surroundingNodes(bx,0);
    const Box ybx = surroundingNodes(bx,1);
    const Box zbx = surroundingNodes(bx,2);

    ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        (flx_arr[0])(i,j,k,0) = rho_u(i,j,k) / mf_u(i,j,0);
        avg_xmom(i,j,k) = (flx_arr[0])(i,j,k,0);
    });
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        (flx_arr[1])(i,j,k,0) = rho_v(i,j,k) / mf_v(i,j,0);
        avg_ymom(i,j,k) = (flx_arr[1])(i,j,k,0);
    });
    ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfsq = mf_m(i,j,0) * mf_m(i,j,0);
        (flx_arr[2])(i,j,k,0) = Omega(i,j,k) / mfsq;
        avg_zmom(i,j,k  ) = (flx_arr[2])(i,j,k,0);
    });

#ifdef ERF_USE_EB
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        (flx_arr[0])(i,j,k,0) *= ax_arr(i,j,k);
        avg_xmom(i,j,k)       *= ax_arr(i,j,k);
    });
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        (flx_arr[1])(i,j,k,0) *= ay_arr(i,j,k);
        avg_ymom(i,j,k)       *= ay_arr(i,j,k);
    });
    ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        (flx_arr[2])(i,j,k,0) *= az_arr(i,j,k);
        avg_zmom(i,j,k)       *= az_arr(i,j,k);
    });
#else
    if (use_terrain) {
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real h_zeta = Compute_h_zeta_AtIface(i,j,k,cellSizeInv,z_nd);
            (flx_arr[0])(i,j,k,0) *= h_zeta;
            avg_xmom(i,j,k)       *= h_zeta;
        });
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real h_zeta =  Compute_h_zeta_AtJface(i,j,k,cellSizeInv,z_nd);
            (flx_arr[1])(i,j,k,0) *= h_zeta;
            avg_ymom(i,j,k)       *= h_zeta;
        });
    }
#endif

    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
#ifdef ERF_USE_EB
        Real invdetJ = (vf_arr(i,j,k) > 0.) ?  1. / vf_arr(i,j,k) : 1.;
#else
        Real invdetJ = (use_terrain) ?  1. / detJ(i,j,k) : 1.;
#endif

        Real mfsq = mf_m(i,j,0) * mf_m(i,j,0);

        advectionSrc(i,j,k,0) = - invdetJ * mfsq * (
          ( (flx_arr[0])(i+1,j,k,0) - (flx_arr[0])(i  ,j,k,0) ) * dxInv +
          ( (flx_arr[1])(i,j+1,k,0) - (flx_arr[1])(i,j  ,k,0) ) * dyInv +
          ( (flx_arr[2])(i,j,k+1,0) - (flx_arr[2])(i,j,k  ,0) ) * dzInv );
    });
}

/**
 * Function for computing the advective tendency for the update equations for all scalars other than rho and (rho theta)
 * This routine has explicit expressions for all cases (terrain or not) when
 * the horizontal and vertial spatial orders are <= 2, and calls more specialized
 * functions when either (or both) spatial order(s) is greater than 2.
 *
 * @param[in] bx box over which the scalars are updated if no external boundary conditions
 * @param[in] icomp component of first scalar to be updated
 * @param[in] ncomp number of components to be updated
 * @param[in] avg_xmom x-component of time-averaged momentum defined in this routine
 * @param[in] avg_ymom y-component of time-averaged momentum defined in this routine
 * @param[in] avg_zmom z-component of time-averaged momentum defined in this routine
 * @param[in] cell_prim primtive form of scalar variales, here only potential temperature theta
 * @param[out] advectionSrc tendency for the scalar update equation
 * @param[in] detJ Jacobian of the metric transformation (= 1 if use_terrain is false)
 * @param[in] cellSizeInv inverse of the mesh spacing
 * @param[in] mf_m map factor at cell centers
 * @param[in] horiz_adv_type advection scheme to be used in horiz. directions for dry scalars
 * @param[in] vert_adv_type advection scheme to be used in vert. directions for dry scalars
 * @param[in] horiz_upw_frac upwinding fraction to be used in horiz. directions for dry scalars (for Blended schemes only)
 * @param[in] vert_upw_frac upwinding fraction to be used in vert. directions for dry scalars (for Blended schemes only)
 * @param[in] use_terrain if true, use the terrain-aware derivatives (with metric terms)
 */

void
AdvectionSrcForScalars (const Box& bx, const int icomp, const int ncomp,
                        const Array4<const Real>& avg_xmom,
                        const Array4<const Real>& avg_ymom,
                        const Array4<const Real>& avg_zmom,
                        const Array4<const Real>& cell_prim,
                        const Array4<Real>& advectionSrc,
                        const Array4<const Real>& detJ,
                        const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                        const Array4<const Real>& mf_m,
                        const AdvType horiz_adv_type,
                        const AdvType vert_adv_type,
                        const Real horiz_upw_frac,
                        const Real vert_upw_frac,
#ifdef ERF_USE_EB
                        const Array4<const Real>& vf_arr,
#endif
                        const bool use_terrain,
                        const GpuArray<const Array4<Real>, AMREX_SPACEDIM>& flx_arr)
{
    BL_PROFILE_VAR("AdvectionSrcForScalars", AdvectionSrcForScalars);
    auto dxInv =     cellSizeInv[0], dyInv =     cellSizeInv[1], dzInv =     cellSizeInv[2];

    const Box xbx = surroundingNodes(bx,0);
    const Box ybx = surroundingNodes(bx,1);
    const Box zbx = surroundingNodes(bx,2);

    // Inline with 2nd order for efficiency
    // NOTE: we don't need to weight avg_xmom, avg_ymom, avg_zmom with terrain metrics
    //       (or with EB area fractions)
    //       because that was done when they were constructed in AdvectionSrcForRhoAndTheta
    if (horiz_adv_type == AdvType::Centered_2nd && vert_adv_type == AdvType::Centered_2nd)
    {
        ParallelFor(xbx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int cons_index = icomp + n;
            const int prim_index = cons_index - 1;
            const Real prim_on_face = 0.5 * (cell_prim(i,j,k,prim_index) + cell_prim(i-1,j,k,prim_index));
            (flx_arr[0])(i,j,k,cons_index) = avg_xmom(i,j,k) * prim_on_face;
        });
        ParallelFor(ybx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int cons_index = icomp + n;
            const int prim_index = cons_index - 1;
            const Real prim_on_face = 0.5 * (cell_prim(i,j,k,prim_index) + cell_prim(i,j-1,k,prim_index));
            (flx_arr[1])(i,j,k,cons_index) = avg_ymom(i,j,k) * prim_on_face;
        });
        ParallelFor(zbx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int cons_index = icomp + n;
            const int prim_index = cons_index - 1;
            const Real prim_on_face = 0.5 * (cell_prim(i,j,k,prim_index) + cell_prim(i,j,k-1,prim_index));
            (flx_arr[2])(i,j,k,cons_index) = avg_zmom(i,j,k) * prim_on_face;
        });

    // Template higher order methods (horizontal first)
    } else {
        switch(horiz_adv_type) {
        case AdvType::Centered_2nd:
            AdvectionSrcForScalarsVert<CENTERED2>(bx, ncomp, icomp, flx_arr, cell_prim,
                                                  avg_xmom, avg_ymom, avg_zmom,
                                                  horiz_upw_frac, vert_upw_frac, vert_adv_type);
            break;
        case AdvType::Upwind_3rd:
            AdvectionSrcForScalarsVert<UPWIND3>(bx, ncomp, icomp, flx_arr, cell_prim,
                                                avg_xmom, avg_ymom, avg_zmom,
                                                horiz_upw_frac, vert_upw_frac, vert_adv_type);
            break;
        case AdvType::Centered_4th:
            AdvectionSrcForScalarsVert<CENTERED4>(bx, ncomp, icomp, flx_arr, cell_prim,
                                                  avg_xmom, avg_ymom, avg_zmom,
                                                  horiz_upw_frac, vert_upw_frac, vert_adv_type);
            break;
        case AdvType::Upwind_5th:
            AdvectionSrcForScalarsVert<UPWIND5>(bx, ncomp, icomp, flx_arr, cell_prim,
                                                avg_xmom, avg_ymom, avg_zmom,
                                                horiz_upw_frac, vert_upw_frac, vert_adv_type);
            break;
        case AdvType::Centered_6th:
            AdvectionSrcForScalarsVert<CENTERED6>(bx, ncomp, icomp, flx_arr, cell_prim,
                                                  avg_xmom, avg_ymom, avg_zmom,
                                                  horiz_upw_frac, vert_upw_frac, vert_adv_type);
            break;
        case AdvType::Weno_3:
            AdvectionSrcForScalarsWrapper<WENO3,WENO3>(bx, ncomp, icomp, flx_arr, cell_prim,
                                                       avg_xmom, avg_ymom, avg_zmom,
                                                       horiz_upw_frac, vert_upw_frac);
            break;
        case AdvType::Weno_5:
            AdvectionSrcForScalarsWrapper<WENO5,WENO5>(bx, ncomp, icomp, flx_arr, cell_prim,
                                                       avg_xmom, avg_ymom, avg_zmom,
                                                       horiz_upw_frac, vert_upw_frac);
            break;
        case AdvType::Weno_3Z:
            AdvectionSrcForScalarsWrapper<WENO_Z3,WENO_Z3>(bx, ncomp, icomp, flx_arr, cell_prim,
                                                           avg_xmom, avg_ymom, avg_zmom,
                                                           horiz_upw_frac, vert_upw_frac);
            break;
        case AdvType::Weno_3MZQ:
            AdvectionSrcForScalarsWrapper<WENO_MZQ3,WENO_MZQ3>(bx, ncomp, icomp, flx_arr, cell_prim,
                                                               avg_xmom, avg_ymom, avg_zmom,
                                                               horiz_upw_frac, vert_upw_frac);
            break;
        case AdvType::Weno_5Z:
            AdvectionSrcForScalarsWrapper<WENO_Z5,WENO_Z5>(bx, ncomp, icomp, flx_arr, cell_prim,
                                                           avg_xmom, avg_ymom, avg_zmom,
                                                           horiz_upw_frac, vert_upw_frac);
            break;
        default:
            AMREX_ASSERT_WITH_MESSAGE(false, "Unknown advection scheme!");
        }
    }

    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
#ifdef ERF_USE_EB
        Real invdetJ = (vf_arr(i,j,k) > 0.) ?  1. / vf_arr(i,j,k) : 1.;
#else
        Real invdetJ = (use_terrain) ?  1. / detJ(i,j,k) : 1.;
#endif

        Real mfsq = mf_m(i,j,0) * mf_m(i,j,0);

        const int cons_index = icomp + n;
        advectionSrc(i,j,k,cons_index) = - invdetJ * mfsq * (
          ( (flx_arr[0])(i+1,j,k,cons_index) - (flx_arr[0])(i  ,j,k,cons_index) ) * dxInv +
          ( (flx_arr[1])(i,j+1,k,cons_index) - (flx_arr[1])(i,j  ,k,cons_index) ) * dyInv +
          ( (flx_arr[2])(i,j,k+1,cons_index) - (flx_arr[2])(i,j,k  ,cons_index) ) * dzInv );
    });
}
