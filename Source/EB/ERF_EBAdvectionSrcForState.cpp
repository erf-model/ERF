#include <ERF_EBAdvection.H>
#include <ERF_EBAdvectionSrcForScalars.H>
#include <ERF_IndexDefines.H>
#include <ERF_TerrainMetrics.H>
#include <ERF_Advection.H>
#include <ERF_AdvectionSrcForScalars.H>

using namespace amrex;

/**
 * Function for computing the advective tendency for the update equations for rho and (rho theta)
 * This routine has explicit expressions for all cases (terrain or not) when
 * the horizontal and vertical spatial orders are <= 2, and calls more specialized
 * functions when either (or both) spatial order(s) is greater than two
 *
 * @param[in] bx box over which the scalars are updated
 * @param[out] advectionSrc tendency for the scalar update equation
 * @param[in] rho_u x-component of momentum
 * @param[in] rho_v y-component of momentum
 * @param[in] Omega component of momentum normal to the z-coordinate surface
 * @param[out] avg_xmom x-component of time-averaged momentum defined in this routine
 * @param[out] avg_ymom y-component of time-averaged momentum defined in this routine
 * @param[out] avg_zmom z-component of time-averaged momentum defined in this routine
 * @param[in] detJ Jacobian of the metric transformation (= 1 if use_terrain is false)
 * @param[in] cellSizeInv inverse of the mesh spacing
 * @param[in] mf_m map factor at cell centers
 * @param[in] mf_u map factor at x-faces
 * @param[in] mf_v map factor at y-faces
 */

void
EBAdvectionSrcForRho (const Box& bx,
                      const Array4<Real>& advectionSrc,
                      const Array4<const Real>& rho_u,
                      const Array4<const Real>& rho_v,
                      const Array4<const Real>& Omega,
                      const Array4<      Real>& avg_xmom,
                      const Array4<      Real>& avg_ymom,
                      const Array4<      Real>& avg_zmom,
                      const Array4<const int>& mask_arr,
                      const Array4<const EBCellFlag>& cfg_arr,
                      const Array4<const Real>& ax_arr,
                      const Array4<const Real>& ay_arr,
                      const Array4<const Real>& az_arr,
                      const Array4<const Real>& fcx_arr,
                      const Array4<const Real>& fcy_arr,
                      const Array4<const Real>& fcz_arr,
                      const Array4<const Real>& detJ,
                      const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                      const Array4<const Real>& mf_mx,
                      const Array4<const Real>& mf_my,
                      const Array4<const Real>& mf_uy,
                      const Array4<const Real>& mf_vx,
                      const GpuArray<const Array4<Real>, AMREX_SPACEDIM>& flx_arr,
                      const bool fixed_rho,
                      bool already_on_centroids)
{
    BL_PROFILE_VAR("EBAdvectionSrcForRho", EBAdvectionSrcForRho);
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    const Box xbx = surroundingNodes(bx,0).grow(IntVect(0, 1, 1));
    const Box ybx = surroundingNodes(bx,1).grow(IntVect(1, 0, 1));
    const Box zbx = surroundingNodes(bx,2).grow(IntVect(1, 1, 0));

    ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        flx_arr[0](i,j,k,0) = rho_u(i,j,k) / mf_uy(i,j,0);
        avg_xmom(i,j,k) = flx_arr[0](i,j,k,0);
    });
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        flx_arr[1](i,j,k,0) = rho_v(i,j,k) / mf_vx(i,j,0);
        avg_ymom(i,j,k) = flx_arr[1](i,j,k,0);
    });
    ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfsq = mf_mx(i,j,0) * mf_my(i,j,0);
        flx_arr[2](i,j,k,0) = Omega(i,j,k) / mfsq;
        avg_zmom(i,j,k) = flx_arr[2](i,j,k,0);
    });

    if (fixed_rho) {
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            advectionSrc(i,j,k,0) = zero;
        });
    } else
    {
        if (already_on_centroids) {

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (detJ(i,j,k) > zero) {
                    Real mfsq = mf_mx(i,j,0) * mf_my(i,j,0);
                    advectionSrc(i,j,k,0) = - mfsq / detJ(i,j,k) * (
                                            ( ax_arr(i+1,j,k) * flx_arr[0](i+1,j,k,0) - ax_arr(i,j,k) * flx_arr[0](i,j,k,0) ) * dxInv +
                                            ( ay_arr(i,j+1,k) * flx_arr[1](i,j+1,k,0) - ay_arr(i,j,k) * flx_arr[1](i,j,k,0) ) * dyInv +
                                            ( az_arr(i,j,k+1) * flx_arr[2](i,j,k+1,0) - az_arr(i,j,k) * flx_arr[2](i,j,k,0) ) * dzInv );
                } else {
                    advectionSrc(i,j,k,0) = zero;
                }
            });

        } else {

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (detJ(i,j,k) > zero) {
                    Real mfsq = mf_mx(i,j,0) * mf_my(i,j,0);
                    if (cfg_arr(i,j,k).isCovered())
                    {
                        advectionSrc(i,j,k,0) = zero;
                    }
                    else if (cfg_arr(i,j,k).isRegular())
                    {
                        advectionSrc(i,j,k,0) = - mfsq / detJ(i,j,k) * (
                                                ( ax_arr(i+1,j,k) * flx_arr[0](i+1,j,k,0) - ax_arr(i,j,k) * flx_arr[0](i,j,k,0) ) * dxInv +
                                                ( ay_arr(i,j+1,k) * flx_arr[1](i,j+1,k,0) - ay_arr(i,j,k) * flx_arr[1](i,j,k,0) ) * dyInv +
                                                ( az_arr(i,j,k+1) * flx_arr[2](i,j,k+1,0) - az_arr(i,j,k) * flx_arr[2](i,j,k,0) ) * dzInv );
                    }
                    else
                    {
                    // Bilinear interpolation
                    Real fxm = flx_arr[0](i,j,k,0);
                    if (ax_arr(i,j,k) != zero && ax_arr(i,j,k) != one) {
                        int jj = j + static_cast<int>(std::copysign(one, fcx_arr(i,j,k,0)));
                        int kk = k + static_cast<int>(std::copysign(one, fcx_arr(i,j,k,1)));
                        Real fracy = (mask_arr(i-1,jj,k) || mask_arr(i,jj,k)) ? std::abs(fcx_arr(i,j,k,0)) : zero;
                        Real fracz = (mask_arr(i-1,j,kk) || mask_arr(i,j,kk)) ? std::abs(fcx_arr(i,j,k,1)) : zero;
                        fxm = (one-fracy)*(one-fracz)*fxm
                            +      fracy *(one-fracz)*flx_arr[0](i,jj,k ,0)
                            +      fracz *(one-fracy)*flx_arr[0](i,j ,kk,0)
                            +      fracy *     fracz *flx_arr[0](i,jj,kk,0);
                    }

                    Real fxp = flx_arr[0](i+1,j,k,0);
                    if (ax_arr(i+1,j,k) != zero && ax_arr(i+1,j,k) != one) {
                        int jj = j + static_cast<int>(std::copysign(one,fcx_arr(i+1,j,k,0)));
                        int kk = k + static_cast<int>(std::copysign(one,fcx_arr(i+1,j,k,1)));
                        Real fracy = (mask_arr(i,jj,k) || mask_arr(i+1,jj,k)) ? std::abs(fcx_arr(i+1,j,k,0)) : zero;
                        Real fracz = (mask_arr(i,j,kk) || mask_arr(i+1,j,kk)) ? std::abs(fcx_arr(i+1,j,k,1)) : zero;
                        fxp = (one-fracy)*(one-fracz)*fxp
                            +      fracy *(one-fracz)*flx_arr[0](i+1,jj,k ,0)
                            +      fracz *(one-fracy)*flx_arr[0](i+1,j ,kk,0)
                            +      fracy *     fracz *flx_arr[0](i+1,jj,kk,0);
                    }

                    Real fym = flx_arr[1](i,j,k,0);
                    if (ay_arr(i,j,k) != zero && ay_arr(i,j,k) != one) {
                        int ii = i + static_cast<int>(std::copysign(one,fcy_arr(i,j,k,0)));
                        int kk = k + static_cast<int>(std::copysign(one,fcy_arr(i,j,k,1)));
                        Real fracx = (mask_arr(ii,j-1,k) || mask_arr(ii,j,k)) ? std::abs(fcy_arr(i,j,k,0)) : zero;
                        Real fracz = (mask_arr(i,j-1,kk) || mask_arr(i,j,kk)) ? std::abs(fcy_arr(i,j,k,1)) : zero;
                        fym = (one-fracx)*(one-fracz)*fym
                            +      fracx *(one-fracz)*flx_arr[1](ii,j,k ,0)
                            +      fracz *(one-fracx)*flx_arr[1](i ,j,kk,0)
                            +      fracx *     fracz *flx_arr[1](ii,j,kk,0);
                    }

                    Real fyp = flx_arr[1](i,j+1,k,0);
                    if (ay_arr(i,j+1,k) != zero && ay_arr(i,j+1,k) != one) {
                        int ii = i + static_cast<int>(std::copysign(one,fcy_arr(i,j+1,k,0)));
                        int kk = k + static_cast<int>(std::copysign(one,fcy_arr(i,j+1,k,1)));
                        Real fracx = (mask_arr(ii,j,k) || mask_arr(ii,j+1,k)) ? std::abs(fcy_arr(i,j+1,k,0)) : zero;
                        Real fracz = (mask_arr(i,j,kk) || mask_arr(i,j+1,kk)) ? std::abs(fcy_arr(i,j+1,k,1)) : zero;
                        fyp = (one-fracx)*(one-fracz)*fyp
                            +      fracx *(one-fracz)*flx_arr[1](ii,j+1,k ,0)
                            +      fracz *(one-fracx)*flx_arr[1](i ,j+1,kk,0)
                            +      fracx *     fracz *flx_arr[1](ii,j+1,kk,0);
                    }

                    Real fzm = flx_arr[2](i,j,k,0);
                    if (az_arr(i,j,k) != zero && az_arr(i,j,k) != one) {
                        int ii = i + static_cast<int>(std::copysign(one,fcz_arr(i,j,k,0)));
                        int jj = j + static_cast<int>(std::copysign(one,fcz_arr(i,j,k,1)));
                        Real fracx = (mask_arr(ii,j,k-1) || mask_arr(ii,j,k)) ? std::abs(fcz_arr(i,j,k,0)) : zero;
                        Real fracy = (mask_arr(i,jj,k-1) || mask_arr(i,jj,k)) ? std::abs(fcz_arr(i,j,k,1)) : zero;
                        fzm = (one-fracx)*(one-fracy)*fzm
                            +      fracx *(one-fracy)*flx_arr[2](ii,j ,k,0)
                            +      fracy *(one-fracx)*flx_arr[2](i ,jj,k,0)
                            +      fracx *     fracy *flx_arr[2](ii,jj,k,0);
                    }

                    Real fzp = flx_arr[2](i,j,k+1,0);
                    if (az_arr(i,j,k+1) != zero && az_arr(i,j,k+1) != one) {
                        int ii = i + static_cast<int>(std::copysign(one,fcz_arr(i,j,k+1,0)));
                        int jj = j + static_cast<int>(std::copysign(one,fcz_arr(i,j,k+1,1)));
                        Real fracx = (mask_arr(ii,j,k) || mask_arr(ii,j,k+1)) ? std::abs(fcz_arr(i,j,k+1,0)) : zero;
                        Real fracy = (mask_arr(i,jj,k) || mask_arr(i,jj,k+1)) ? std::abs(fcz_arr(i,j,k+1,1)) : zero;
                        fzp = (one-fracx)*(one-fracy)*fzp
                            +      fracx *(one-fracy)*flx_arr[2](ii,j ,k+1,0)
                            +      fracy *(one-fracx)*flx_arr[2](i ,jj,k+1,0)
                            +      fracx *     fracy *flx_arr[2](ii,jj,k+1,0);
                    }

                    advectionSrc(i,j,k,0) = -  mfsq / detJ(i,j,k) * (
                                ( ax_arr(i+1,j,k) * fxp - ax_arr(i,j,k) * fxm ) * dxInv +
                                ( ay_arr(i,j+1,k) * fyp - ay_arr(i,j,k) * fym ) * dyInv +
                                ( az_arr(i,j,k+1) * fzp - az_arr(i,j,k) * fzm ) * dzInv );
                    }
                } else {
                    advectionSrc(i,j,k,0) = zero;
                }
            });

        } // already_on_centroids
    } // fixed_rho
}

/**
 * Function for computing the advective tendency for scalars when terrain_type is EB.
 *
 * @param[in] bx box over which the scalars are updated if no external boundary conditions
 * @param[in] icomp component of first scalar to be updated
 * @param[in] ncomp number of components to be updated
 * @param[in] avg_xmom x-component of time-averaged momentum defined in this routine
 * @param[in] avg_ymom y-component of time-averaged momentum defined in this routine
 * @param[in] avg_zmom z-component of time-averaged momentum defined in this routine
 * @param[in] cell_prim primitive form of scalar variables, here only potential temperature theta
 * @param[out] advectionSrc tendency for the scalar update equation
 * @param[in] mask_arr Cell-centered masks (=1 otherwise, =0 if physbnd)
 * @param[in] cfg_arr Cell-centered flags
 * @param[in] ax_arr Area fraction of x-face
 * @param[in] ay_arr Area fraction of y-face
 * @param[in] az_arr Area fraction of z-face
 * @param[in] fcx_arr Face centroid of x-face
 * @param[in] fcy_arr Face centroid of y-face
 * @param[in] fcz_arr Face centroid of z-face
 * @param[in] detJ Jacobian of the metric transformation (= 1 if use_terrain is false)
 * @param[in] cellSizeInv inverse of the mesh spacing
 * @param[in] mf_mx map factor at cell centers
 * @param[in] mf_my map factor at cell centers
 * @param[in] horiz_adv_type advection scheme to be used in horiz. directions for dry scalars
 * @param[in] vert_adv_type advection scheme to be used in vert. directions for dry scalars
 * @param[in] horiz_upw_frac upwinding fraction to be used in horiz. directions for dry scalars (for Blended schemes only)
 * @param[in] vert_upw_frac upwinding fraction to be used in vert. directions for dry scalars (for Blended schemes only)
 */

void
EBAdvectionSrcForScalars (const Box& bx,
                        const int icomp,
                        const int ncomp,
                        const Array4<const Real>& avg_xmom,
                        const Array4<const Real>& avg_ymom,
                        const Array4<const Real>& avg_zmom,
                        const Array4<const Real>& cell_prim,
                        const Array4<Real>& advectionSrc,
                        const Array4<const int>& mask_arr,
                        const Array4<const EBCellFlag>& cfg_arr,
                        const Array4<const Real>& ax_arr,
                        const Array4<const Real>& ay_arr,
                        const Array4<const Real>& az_arr,
                        const Array4<const Real>& fcx_arr,
                        const Array4<const Real>& fcy_arr,
                        const Array4<const Real>& fcz_arr,
                        const Array4<const Real>& detJ,
                        const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                        const Array4<const Real>& mf_mx,
                        const Array4<const Real>& mf_my,
                        const AdvType horiz_adv_type,
                        const AdvType vert_adv_type,
                        const Real horiz_upw_frac,
                        const Real vert_upw_frac,
                        const GpuArray<const Array4<Real>, AMREX_SPACEDIM>& flx_arr,
                        const Box& domain,
                        const BCRec* bc_ptr_h,
                        bool already_on_centroids)
{
    BL_PROFILE_VAR("EBAdvectionSrcForScalars", EBAdvectionSrcForScalars);
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    const Box xbx = surroundingNodes(bx,0);
    const Box ybx = surroundingNodes(bx,1);
    const Box zbx = surroundingNodes(bx,2);

    // Open bc will be imposed upon all vars (we only access cons here for simplicity)
    const bool xlo_open = (bc_ptr_h[BCVars::cons_bc].lo(0) == ERFBCType::open);
    const bool xhi_open = (bc_ptr_h[BCVars::cons_bc].hi(0) == ERFBCType::open);
    const bool ylo_open = (bc_ptr_h[BCVars::cons_bc].lo(1) == ERFBCType::open);
    const bool yhi_open = (bc_ptr_h[BCVars::cons_bc].hi(1) == ERFBCType::open);

    // Only advection operations in bndry normal direction with OPEN BC
    Box  bx_xlo,  bx_xhi,  bx_ylo,  bx_yhi;
    if (xlo_open) {
        if ( bx.smallEnd(0) == domain.smallEnd(0)) {  bx_xlo = makeSlab( bx,0,domain.smallEnd(0));}
    }
    if (xhi_open) {
        if ( bx.bigEnd(0) == domain.bigEnd(0))     {  bx_xhi = makeSlab( bx,0,domain.bigEnd(0)  );}
    }
    if (ylo_open) {
        if ( bx.smallEnd(1) == domain.smallEnd(1)) {  bx_ylo = makeSlab( bx,1,domain.smallEnd(1));}
    }
    if (yhi_open) {
        if ( bx.bigEnd(1) == domain.bigEnd(1))     {  bx_yhi = makeSlab( bx,1,domain.bigEnd(1)  );}
    }

    // Inline with 2nd order for efficiency
    // NOTE: For EB, avg_xmom, avg_ymom, avg_zmom were are weighted by area fractions in AdvectionSrcForRho
    //       The flux is weighted by area fraction after its interpolation.
    if (horiz_adv_type == AdvType::Centered_2nd && vert_adv_type == AdvType::Centered_2nd)
    {
        ParallelFor(xbx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int cons_index = icomp + n;
            const int prim_index = cons_index - 1;
            const Real prim_on_face = myhalf * (cell_prim(i,j,k,prim_index) + cell_prim(i-1,j,k,prim_index));
            flx_arr[0](i,j,k,cons_index) = avg_xmom(i,j,k) * prim_on_face;
        });
        ParallelFor(ybx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int cons_index = icomp + n;
            const int prim_index = cons_index - 1;
            const Real prim_on_face = myhalf * (cell_prim(i,j,k,prim_index) + cell_prim(i,j-1,k,prim_index));
            flx_arr[1](i,j,k,cons_index) = avg_ymom(i,j,k) * prim_on_face;
        });
        ParallelFor(zbx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int cons_index = icomp + n;
            const int prim_index = cons_index - 1;
            const Real prim_on_face = myhalf * (cell_prim(i,j,k,prim_index) + cell_prim(i,j,k-1,prim_index));
            flx_arr[2](i,j,k,cons_index) = avg_zmom(i,j,k) * prim_on_face;
        });

    // Template higher order methods (horizontal first)
    } else {
        switch(horiz_adv_type) {
        case AdvType::Centered_2nd:
            EBAdvectionSrcForScalarsVert<CENTERED2>(bx, ncomp, icomp, flx_arr, cell_prim,
                                                avg_xmom, avg_ymom, avg_zmom, cfg_arr, ax_arr, ay_arr, az_arr,
                                                horiz_upw_frac, vert_upw_frac, vert_adv_type);
            break;
        case AdvType::Upwind_3rd:
            EBAdvectionSrcForScalarsVert<UPWIND3>(bx, ncomp, icomp, flx_arr, cell_prim,
                                                avg_xmom, avg_ymom, avg_zmom, cfg_arr, ax_arr, ay_arr, az_arr,
                                                horiz_upw_frac, vert_upw_frac, vert_adv_type);
            break;
        case AdvType::Centered_4th:
            EBAdvectionSrcForScalarsVert<CENTERED4>(bx, ncomp, icomp, flx_arr, cell_prim,
                                                avg_xmom, avg_ymom, avg_zmom, cfg_arr, ax_arr, ay_arr, az_arr,
                                                horiz_upw_frac, vert_upw_frac, vert_adv_type);
            break;
        case AdvType::Upwind_5th:
            EBAdvectionSrcForScalarsVert<UPWIND5>(bx, ncomp, icomp, flx_arr, cell_prim,
                                                avg_xmom, avg_ymom, avg_zmom, cfg_arr, ax_arr, ay_arr, az_arr,
                                                horiz_upw_frac, vert_upw_frac, vert_adv_type);
            break;
        case AdvType::Centered_6th:
            EBAdvectionSrcForScalarsVert<CENTERED6>(bx, ncomp, icomp, flx_arr, cell_prim,
                                                avg_xmom, avg_ymom, avg_zmom, cfg_arr, ax_arr, ay_arr, az_arr,
                                                horiz_upw_frac, vert_upw_frac, vert_adv_type);
            break;
        default:
            AMREX_ASSERT_WITH_MESSAGE(false, "Unknown advection scheme! WENO is currently not supported for EB.");
        }
    }

    // Compute divergence

    if (already_on_centroids) {

        ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int cons_index = icomp + n;
            if (detJ(i,j,k) > zero)
            {
                Real invdetJ = one / detJ(i,j,k);
                Real mfsq    = mf_mx(i,j,0) * mf_my(i,j,0);

                advectionSrc(i,j,k,cons_index) = - invdetJ * mfsq * (
                  ( ax_arr(i+1,j,k) * flx_arr[0](i+1,j,k,cons_index) - ax_arr(i,j,k) * flx_arr[0](i  ,j,k,cons_index) ) * dxInv +
                  ( ay_arr(i,j+1,k) * flx_arr[1](i,j+1,k,cons_index) - ay_arr(i,j,k) * flx_arr[1](i,j  ,k,cons_index) ) * dyInv +
                  ( az_arr(i,j,k+1) * flx_arr[2](i,j,k+1,cons_index) - az_arr(i,j,k) * flx_arr[2](i,j,k  ,cons_index) ) * dzInv );
            } else {
                advectionSrc(i,j,k,cons_index) = zero;
            }
        });

    } else {

        AMREX_HOST_DEVICE_FOR_4D(bx,ncomp,i,j,k,n,
        {
            const int cons_index = icomp + n;
            if (detJ(i,j,k) > zero)
            {
                Real invdetJ = one / detJ(i,j,k);
                Real mfsq    = mf_mx(i,j,0) * mf_my(i,j,0);
                if (cfg_arr(i,j,k).isCovered())
                {
                    advectionSrc(i,j,k,cons_index) = zero;
                }
                else if (cfg_arr(i,j,k).isRegular())
                {
                    advectionSrc(i,j,k,cons_index) = - invdetJ * mfsq * (
                        ( ax_arr(i+1,j,k) * flx_arr[0](i+1,j,k,cons_index) - ax_arr(i,j,k) * flx_arr[0](i  ,j,k,cons_index) ) * dxInv +
                        ( ay_arr(i,j+1,k) * flx_arr[1](i,j+1,k,cons_index) - ay_arr(i,j,k) * flx_arr[1](i,j  ,k,cons_index) ) * dyInv +
                        ( az_arr(i,j,k+1) * flx_arr[2](i,j,k+1,cons_index) - az_arr(i,j,k) * flx_arr[2](i,j,k  ,cons_index) ) * dzInv );
                }
                else
                {
                    // Bilinear interpolation
                    Real fxm = flx_arr[0](i,j,k,cons_index);
                    if (ax_arr(i,j,k) != zero && ax_arr(i,j,k) != one) {
                        int jj = j + static_cast<int>(std::copysign(one, fcx_arr(i,j,k,0)));
                        int kk = k + static_cast<int>(std::copysign(one, fcx_arr(i,j,k,1)));
                        Real fracy = (mask_arr(i-1,jj,k) || mask_arr(i,jj,k)) ? std::abs(fcx_arr(i,j,k,0)) : zero;
                        Real fracz = (mask_arr(i-1,j,kk) || mask_arr(i,j,kk)) ? std::abs(fcx_arr(i,j,k,1)) : zero;
                        fxm = (one-fracy)*(one-fracz)*fxm
                            +      fracy *(one-fracz)*flx_arr[0](i,jj,k ,cons_index)
                            +      fracz *(one-fracy)*flx_arr[0](i,j ,kk,cons_index)
                            +      fracy *     fracz *flx_arr[0](i,jj,kk,cons_index);
                    }

                    Real fxp = flx_arr[0](i+1,j,k,cons_index);
                    if (ax_arr(i+1,j,k) != zero && ax_arr(i+1,j,k) != one) {
                        int jj = j + static_cast<int>(std::copysign(one,fcx_arr(i+1,j,k,0)));
                        int kk = k + static_cast<int>(std::copysign(one,fcx_arr(i+1,j,k,1)));
                        Real fracy = (mask_arr(i,jj,k) || mask_arr(i+1,jj,k)) ? std::abs(fcx_arr(i+1,j,k,0)) : zero;
                        Real fracz = (mask_arr(i,j,kk) || mask_arr(i+1,j,kk)) ? std::abs(fcx_arr(i+1,j,k,1)) : zero;
                        fxp = (one-fracy)*(one-fracz)*fxp
                            +      fracy *(one-fracz)*flx_arr[0](i+1,jj,k ,cons_index)
                            +      fracz *(one-fracy)*flx_arr[0](i+1,j ,kk,cons_index)
                            +      fracy *     fracz *flx_arr[0](i+1,jj,kk,cons_index);
                    }

                    Real fym = flx_arr[1](i,j,k,cons_index);
                    if (ay_arr(i,j,k) != zero && ay_arr(i,j,k) != one) {
                        int ii = i + static_cast<int>(std::copysign(one,fcy_arr(i,j,k,0)));
                        int kk = k + static_cast<int>(std::copysign(one,fcy_arr(i,j,k,1)));
                        Real fracx = (mask_arr(ii,j-1,k) || mask_arr(ii,j,k)) ? std::abs(fcy_arr(i,j,k,0)) : zero;
                        Real fracz = (mask_arr(i,j-1,kk) || mask_arr(i,j,kk)) ? std::abs(fcy_arr(i,j,k,1)) : zero;
                        fym = (one-fracx)*(one-fracz)*fym
                            +      fracx *(one-fracz)*flx_arr[1](ii,j,k ,cons_index)
                            +      fracz *(one-fracx)*flx_arr[1](i ,j,kk,cons_index)
                            +      fracx *     fracz *flx_arr[1](ii,j,kk,cons_index);
                    }

                    Real fyp = flx_arr[1](i,j+1,k,cons_index);
                    if (ay_arr(i,j+1,k) != zero && ay_arr(i,j+1,k) != one) {
                        int ii = i + static_cast<int>(std::copysign(one,fcy_arr(i,j+1,k,0)));
                        int kk = k + static_cast<int>(std::copysign(one,fcy_arr(i,j+1,k,1)));
                        Real fracx = (mask_arr(ii,j,k) || mask_arr(ii,j+1,k)) ? std::abs(fcy_arr(i,j+1,k,0)) : zero;
                        Real fracz = (mask_arr(i,j,kk) || mask_arr(i,j+1,kk)) ? std::abs(fcy_arr(i,j+1,k,1)) : zero;
                        fyp = (one-fracx)*(one-fracz)*fyp
                            +      fracx *(one-fracz)*flx_arr[1](ii,j+1,k ,cons_index)
                            +      fracz *(one-fracx)*flx_arr[1](i ,j+1,kk,cons_index)
                            +      fracx *     fracz *flx_arr[1](ii,j+1,kk,cons_index);
                    }

                    Real fzm = flx_arr[2](i,j,k,cons_index);
                    if (az_arr(i,j,k) != zero && az_arr(i,j,k) != one) {
                        int ii = i + static_cast<int>(std::copysign(one,fcz_arr(i,j,k,0)));
                        int jj = j + static_cast<int>(std::copysign(one,fcz_arr(i,j,k,1)));
                        Real fracx = (mask_arr(ii,j,k-1) || mask_arr(ii,j,k)) ? std::abs(fcz_arr(i,j,k,0)) : zero;
                        Real fracy = (mask_arr(i,jj,k-1) || mask_arr(i,jj,k)) ? std::abs(fcz_arr(i,j,k,1)) : zero;
                        fzm = (one-fracx)*(one-fracy)*fzm
                            +      fracx *(one-fracy)*flx_arr[2](ii,j ,k,cons_index)
                            +      fracy *(one-fracx)*flx_arr[2](i ,jj,k,cons_index)
                            +      fracx *     fracy *flx_arr[2](ii,jj,k,cons_index);
                    }

                    Real fzp = flx_arr[2](i,j,k+1,cons_index);
                    if (az_arr(i,j,k+1) != zero && az_arr(i,j,k+1) != one) {
                        int ii = i + static_cast<int>(std::copysign(one,fcz_arr(i,j,k+1,0)));
                        int jj = j + static_cast<int>(std::copysign(one,fcz_arr(i,j,k+1,1)));
                        Real fracx = (mask_arr(ii,j,k) || mask_arr(ii,j,k+1)) ? std::abs(fcz_arr(i,j,k+1,0)) : zero;
                        Real fracy = (mask_arr(i,jj,k) || mask_arr(i,jj,k+1)) ? std::abs(fcz_arr(i,j,k+1,1)) : zero;
                        fzp = (one-fracx)*(one-fracy)*fzp
                            +      fracx *(one-fracy)*flx_arr[2](ii,j ,k+1,cons_index)
                            +      fracy *(one-fracx)*flx_arr[2](i ,jj,k+1,cons_index)
                            +      fracx *     fracy *flx_arr[2](ii,jj,k+1,cons_index);
                    }

                    advectionSrc(i,j,k,cons_index) = - invdetJ * mfsq * (
                          ( ax_arr(i+1,j,k) * fxp - ax_arr(i,j,k) * fxm ) * dxInv
                        + ( ay_arr(i,j+1,k) * fyp - ay_arr(i,j,k) * fym ) * dyInv
                        + ( az_arr(i,j,k+1) * fzp - az_arr(i,j,k) * fzm ) * dzInv );
                }

                // eb_compute_divergence(i,j,k,n,advectionSrc,AMREX_D_DECL(flx_arr[0],flx_arr[1],flx_arr[2]),
                //                     mask_arr, cfg_arr, detJ, AMREX_D_DECL(ax_arr,ay_arr,az_arr),
                //                     AMREX_D_DECL(fcx_arr,fcy_arr,fcz_arr), cellSizeInv, already_on_centroids);


            } else {
                advectionSrc(i,j,k,cons_index) = zero;
            }
        });

    }

    // Special advection operator for open BC (bndry tangent operations)
    if (xlo_open) {
        bool do_lo = true;
        AdvectionSrcForOpenBC_Tangent_Cons(bx_xlo, 0, icomp, ncomp, advectionSrc, cell_prim,
                                           avg_xmom, avg_ymom, avg_zmom,
                                           detJ, cellSizeInv, do_lo);
    }
    if (xhi_open) {
        AdvectionSrcForOpenBC_Tangent_Cons(bx_xhi, 0, icomp, ncomp, advectionSrc, cell_prim,
                                           avg_xmom, avg_ymom, avg_zmom,
                                           detJ, cellSizeInv);
    }
    if (ylo_open) {
        bool do_lo = true;
        AdvectionSrcForOpenBC_Tangent_Cons(bx_ylo, 1, icomp, ncomp, advectionSrc, cell_prim,
                                           avg_xmom, avg_ymom, avg_zmom,
                                           detJ, cellSizeInv, do_lo);
    }
    if (yhi_open) {
        AdvectionSrcForOpenBC_Tangent_Cons(bx_yhi, 1, icomp, ncomp, advectionSrc, cell_prim,
                                           avg_xmom, avg_ymom, avg_zmom,
                                           detJ, cellSizeInv);
    }
}
