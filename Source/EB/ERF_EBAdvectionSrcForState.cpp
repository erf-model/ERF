#include <ERF_EBAdvection.H>
#include <ERF_EBAdvectionSrcForScalars.H>
#include <ERF_IndexDefines.H>
#include <ERF_TerrainMetrics.H>
#include <ERF_Advection.H>
#include <ERF_AdvectionSrcForScalars.H>

using namespace amrex;

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
 * @param[in] detJ Jacobian of the metric transformation (= 1 if use_terrain is false)
 * @param[in] cellSizeInv inverse of the mesh spacing
 * @param[in] mf_m map factor at cell centers
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
                        const Array4<const EBCellFlag>& cfg_arr,
                        const Array4<const Real>& ax_arr,
                        const Array4<const Real>& ay_arr,
                        const Array4<const Real>& az_arr,
                        const Array4<const Real>& detJ,
                        const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                        const Array4<const Real>& mf_m,
                        const AdvType horiz_adv_type,
                        const AdvType vert_adv_type,
                        const Real horiz_upw_frac,
                        const Real vert_upw_frac,
                        const GpuArray<const Array4<Real>, AMREX_SPACEDIM>& flx_arr,
                        const Box& domain,
                        const BCRec* bc_ptr_h)
{
    BL_PROFILE_VAR("EBAdvectionSrcForScalars", EBAdvectionSrcForScalars);
    auto dxInv =     cellSizeInv[0], dyInv =     cellSizeInv[1], dzInv =     cellSizeInv[2];

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

    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const int cons_index = icomp + n;
        if (detJ(i,j,k) > 0.)
        {
            Real invdetJ = 1.0 / detJ(i,j,k);
            Real mfsq    = mf_m(i,j,0) * mf_m(i,j,0);

            advectionSrc(i,j,k,cons_index) = - invdetJ * mfsq * (
              ( (flx_arr[0])(i+1,j,k,cons_index) - (flx_arr[0])(i  ,j,k,cons_index) ) * dxInv +
              ( (flx_arr[1])(i,j+1,k,cons_index) - (flx_arr[1])(i,j  ,k,cons_index) ) * dyInv +
              ( (flx_arr[2])(i,j,k+1,cons_index) - (flx_arr[2])(i,j,k  ,cons_index) ) * dzInv );
        } else {
            advectionSrc(i,j,k,cons_index) = 0.;
        }
    });

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