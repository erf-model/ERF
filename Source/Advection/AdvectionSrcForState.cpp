#include <IndexDefines.H>
#include <TerrainMetrics.H>
#include <Advection.H>
#include <AdvectionSrcForState_N.H>
#include <AdvectionSrcForState_T.H>

using namespace amrex;

/**
 * Function for computing the advective tendency for the update equations for rho and (rho theta)
 * This routine has explicit expressions for all cases (terrain or not) when
 * the horizontal and vertial spatial orders are <= 2, and calls more specialized
 * functions when either (or both) spatial order(s) is greater than 2.
 *
 * @param[in] bx box over which the scalars are updated
 * @param[in] valid_bx box that contains only the cells not in the specified or relaxation zones
 * @param[out] advectionSrc tendency for the scalar update equation
 * @param[in] rho_u x-component of momentum
 * @param[in] rho_v y-component of momentum
 * @param[in] Omega component of momentum normal to the z-coordinate surface
 * @param[in] fac weighting factor for use in defining time-averaged momentum
 * @param[out] avg_xmom x-component of time-averaged momentum defined in this routine
 * @param[out] avg_ymom y-component of time-averaged momentum defined in this routine
 * @param[out] avg_zmom z-component of time-averaged momentum defined in this routine
 * @param[in] cell_prim primtive form of scalar variales, here only potential temperature theta
 * @param[in] z_nd height coordinate at nodes
 * @param[in] detJ Jacobian of the metric transformation (= 1 if use_terrain is false)
 * @param[in] cellSizeInv inverse of the mesh spacing
 * @param[in] mf_m map factor at cell centers
 * @param[in] mf_u map factor at x-faces
 * @param[in] mf_v map factor at y-faces
 * @param[in] horiz_adv_type advection scheme to be used in horiz. directions for dry scalars
 * @param[in] vert_adv_type advection scheme to be used in horiz. directions for dry scalars
 * @param[in] use_terrain if true, use the terrain-aware derivatives (with metric terms)
 */

void
AdvectionSrcForRhoAndTheta (const Box& bx, const Box& valid_bx,
                            const Array4<Real>& advectionSrc,
                            const Array4<const Real>& rho_u,
                            const Array4<const Real>& rho_v,
                            const Array4<const Real>& Omega, Real fac,
                            const Array4<      Real>& avg_xmom,
                            const Array4<      Real>& avg_ymom,
                            const Array4<      Real>& avg_zmom,
                            const Array4<const Real>& cell_prim,
                            const Array4<const Real>& z_nd, const Array4<const Real>& detJ,
                            const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                            const Array4<const Real>& mf_m,
                            const Array4<const Real>& mf_u,
                            const Array4<const Real>& mf_v,
                            const AdvType horiz_adv_type,
                            const AdvType vert_adv_type,
                            const int use_terrain)
{
    BL_PROFILE_VAR("AdvectionSrcForRhoAndTheta", AdvectionSrcForRhoAndTheta);
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    // We note that valid_bx is the actual grid, while bx may be a tile within that grid
    const auto& vbx_hi = ubound(valid_bx);

    if (!use_terrain) {
        // Inline with 2nd order for efficiency
        if (horiz_adv_type == AdvType::Centered_2nd && vert_adv_type == AdvType::Centered_2nd)
        {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real xflux_lo = rho_u(i  ,j,k) / mf_u(i  ,j  ,0);
                Real xflux_hi = rho_u(i+1,j,k) / mf_u(i+1,j  ,0);
                Real yflux_lo = rho_v(i,j  ,k) / mf_v(i  ,j  ,0);
                Real yflux_hi = rho_v(i,j+1,k) / mf_v(i  ,j+1,0);
                Real zflux_lo = Omega(i,j,k  );
                Real zflux_hi = Omega(i,j,k+1);

                avg_xmom(i  ,j,k) += fac*xflux_lo;
                if (i == vbx_hi.x)
                    avg_xmom(i+1,j,k) += fac*xflux_hi;
                avg_ymom(i,j  ,k) += fac*yflux_lo;
                if (j == vbx_hi.y)
                    avg_ymom(i,j+1,k) += fac*yflux_hi;
                avg_zmom(i,j,k  ) += fac*zflux_lo;
                if (k == vbx_hi.z)
                    avg_zmom(i,j,k+1) += fac*zflux_hi;

                Real mfsq = mf_m(i,j,0) * mf_m(i,j,0);

                advectionSrc(i,j,k,0) = -(
                                          ( xflux_hi - xflux_lo ) * dxInv * mfsq +
                                          ( yflux_hi - yflux_lo ) * dyInv * mfsq +
                                          ( zflux_hi - zflux_lo ) * dzInv );

                const int prim_index = 0;
                advectionSrc(i,j,k,1) = - 0.5 * (
              ( xflux_hi * (cell_prim(i+1,j,k,prim_index) + cell_prim(i,j,k,prim_index)) -
                xflux_lo * (cell_prim(i-1,j,k,prim_index) + cell_prim(i,j,k,prim_index)) ) * dxInv * mfsq +
              ( yflux_hi * (cell_prim(i,j+1,k,prim_index) + cell_prim(i,j,k,prim_index)) -
                yflux_lo * (cell_prim(i,j-1,k,prim_index) + cell_prim(i,j,k,prim_index)) ) * dyInv * mfsq +
              ( zflux_hi * (cell_prim(i,j,k+1,prim_index) + cell_prim(i,j,k,prim_index)) -
                zflux_lo * (cell_prim(i,j,k-1,prim_index) + cell_prim(i,j,k,prim_index)) ) * dzInv);
            });
        // Template higher order methods
        } else {
            if (horiz_adv_type == AdvType::Centered_2nd) {
                AdvectionSrcForRhoThetaVert_N<CENTERED2>(bx, vbx_hi, fac, advectionSrc,
                                                         cell_prim, rho_u, rho_v, Omega,
                                                         avg_xmom, avg_ymom, avg_zmom,
                                                         cellSizeInv, mf_m, mf_u, mf_v,
                                                         vert_adv_type);
            } else if (horiz_adv_type == AdvType::Upwind_3rd) {
                AdvectionSrcForRhoThetaVert_N<UPWIND3>(bx, vbx_hi, fac, advectionSrc,
                                                       cell_prim, rho_u, rho_v, Omega,
                                                       avg_xmom, avg_ymom, avg_zmom,
                                                       cellSizeInv, mf_m, mf_u, mf_v,
                                                       vert_adv_type);
            } else if (horiz_adv_type == AdvType::Centered_4th) {
                AdvectionSrcForRhoThetaVert_N<CENTERED4>(bx, vbx_hi, fac, advectionSrc,
                                                         cell_prim, rho_u, rho_v, Omega,
                                                         avg_xmom, avg_ymom, avg_zmom,
                                                         cellSizeInv, mf_m, mf_u, mf_v,
                                                         vert_adv_type);
            } else if (horiz_adv_type == AdvType::Upwind_5th) {
                AdvectionSrcForRhoThetaVert_N<UPWIND5>(bx, vbx_hi, fac, advectionSrc,
                                                       cell_prim, rho_u, rho_v, Omega,
                                                       avg_xmom, avg_ymom, avg_zmom,
                                                       cellSizeInv, mf_m, mf_u, mf_v,
                                                       vert_adv_type);
            } else if (horiz_adv_type == AdvType::Centered_6th) {
                AdvectionSrcForRhoThetaVert_N<CENTERED6>(bx, vbx_hi, fac, advectionSrc,
                                                         cell_prim, rho_u, rho_v, Omega,
                                                         avg_xmom, avg_ymom, avg_zmom,
                                                         cellSizeInv, mf_m, mf_u, mf_v,
                                                         vert_adv_type);
            } else {
                AMREX_ASSERT_WITH_MESSAGE(false, "Unknown advection scheme!");
            }
        }

    } else {
        // Inline with 2nd order for efficiency
        if (horiz_adv_type == AdvType::Centered_2nd && vert_adv_type == AdvType::Centered_2nd)
        {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real invdetJ = 1./ detJ(i,j,k);

                Real xflux_lo = rho_u(i  ,j,k) / mf_u(i  ,j  ,0);
                Real xflux_hi = rho_u(i+1,j,k) / mf_u(i+1,j  ,0);
                Real yflux_lo = rho_v(i,j  ,k) / mf_v(i  ,j  ,0);
                Real yflux_hi = rho_v(i,j+1,k) / mf_v(i  ,j+1,0);
                Real zflux_lo = Omega(i,j,k  );
                Real zflux_hi = Omega(i,j,k+1);

                Real met_h_zeta_xlo = Compute_h_zeta_AtIface(i  ,j  ,k,cellSizeInv,z_nd);
                xflux_lo *= met_h_zeta_xlo;
                Real met_h_zeta_xhi = Compute_h_zeta_AtIface(i+1,j  ,k,cellSizeInv,z_nd);
                xflux_hi *= met_h_zeta_xhi;

                Real met_h_zeta_ylo = Compute_h_zeta_AtJface(i  ,j  ,k,cellSizeInv,z_nd);
                yflux_lo *= met_h_zeta_ylo;
                Real met_h_zeta_yhi = Compute_h_zeta_AtJface(i  ,j+1,k,cellSizeInv,z_nd);
                yflux_hi *= met_h_zeta_yhi;

                avg_xmom(i  ,j,k) += fac*xflux_lo;
                if (i == vbx_hi.x)
                    avg_xmom(i+1,j,k) += fac*xflux_hi;
                avg_ymom(i,j  ,k) += fac*yflux_lo;
                if (j == vbx_hi.y)
                    avg_ymom(i,j+1,k) += fac*yflux_hi;
                avg_zmom(i,j,k  ) += fac*zflux_lo;
                if (k == vbx_hi.z)
                    avg_zmom(i,j,k+1) += fac*zflux_hi;

                Real mfsq = mf_m(i,j,0) * mf_m(i,j,0);

                advectionSrc(i,j,k,0) = - invdetJ * (
                                                     ( xflux_hi - xflux_lo ) * dxInv * mfsq +
                                                     ( yflux_hi - yflux_lo ) * dyInv * mfsq +
                                                     ( zflux_hi - zflux_lo ) * dzInv );

                const int prim_index = 0;
                advectionSrc(i,j,k,1) = - invdetJ * 0.5 * (
                ( xflux_hi * (cell_prim(i,j,k,prim_index) + cell_prim(i+1,j,k,prim_index)) -
                  xflux_lo * (cell_prim(i,j,k,prim_index) + cell_prim(i-1,j,k,prim_index)) ) * dxInv * mfsq +
                ( yflux_hi * (cell_prim(i,j,k,prim_index) + cell_prim(i,j+1,k,prim_index)) -
                  yflux_lo * (cell_prim(i,j,k,prim_index) + cell_prim(i,j-1,k,prim_index)) ) * dyInv * mfsq +
                ( zflux_hi * (cell_prim(i,j,k,prim_index) + cell_prim(i,j,k+1,prim_index)) -
                  zflux_lo * (cell_prim(i,j,k,prim_index) + cell_prim(i,j,k-1,prim_index)) ) * dzInv);
            });
        // Template higher order methods (horizontal first)
        } else {
            if (horiz_adv_type == AdvType::Centered_2nd) {
                AdvectionSrcForRhoThetaVert_T<CENTERED2>(bx, vbx_hi, fac, advectionSrc,
                                                         cell_prim, rho_u, rho_v, Omega,
                                                         avg_xmom, avg_ymom, avg_zmom,
                                                         z_nd, detJ, cellSizeInv, mf_m,
                                                         mf_u, mf_v, vert_adv_type);
            } else if (horiz_adv_type == AdvType::Upwind_3rd) {
                AdvectionSrcForRhoThetaVert_T<UPWIND3>(bx, vbx_hi, fac, advectionSrc,
                                                       cell_prim, rho_u, rho_v, Omega,
                                                       avg_xmom, avg_ymom, avg_zmom,
                                                       z_nd, detJ, cellSizeInv, mf_m,
                                                       mf_u, mf_v, vert_adv_type);
            } else if (horiz_adv_type == AdvType::Centered_4th) {
                AdvectionSrcForRhoThetaVert_T<CENTERED4>(bx, vbx_hi, fac, advectionSrc,
                                                         cell_prim, rho_u, rho_v, Omega,
                                                         avg_xmom, avg_ymom, avg_zmom,
                                                         z_nd, detJ, cellSizeInv, mf_m,
                                                         mf_u, mf_v, vert_adv_type);
            } else if (horiz_adv_type == AdvType::Upwind_5th) {
                AdvectionSrcForRhoThetaVert_T<UPWIND5>(bx, vbx_hi, fac, advectionSrc,
                                                       cell_prim, rho_u, rho_v, Omega,
                                                       avg_xmom, avg_ymom, avg_zmom,
                                                       z_nd, detJ, cellSizeInv, mf_m,
                                                       mf_u, mf_v, vert_adv_type);
            } else if (horiz_adv_type == AdvType::Centered_6th) {
                AdvectionSrcForRhoThetaVert_T<CENTERED6>(bx, vbx_hi, fac, advectionSrc,
                                                         cell_prim, rho_u, rho_v, Omega,
                                                         avg_xmom, avg_ymom, avg_zmom,
                                                         z_nd, detJ, cellSizeInv, mf_m,
                                                         mf_u, mf_v, vert_adv_type);
            } else {
                AMREX_ASSERT_WITH_MESSAGE(false, "Unknown advection scheme!");
            }
        }
    }
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
 * @param[in] vert_adv_type advection scheme to be used in horiz. directions for dry scalars
 * @param[in] use_terrain if true, use the terrain-aware derivatives (with metric terms)
 */

void
AdvectionSrcForScalars (const Box& bx, const int icomp, const int ncomp,
                        const Array4<const Real>& avg_xmom, const Array4<const Real>& avg_ymom,
                        const Array4<const Real>& avg_zmom,
                        const Array4<const Real>& cell_prim,
                        const Array4<Real>& advectionSrc,
                        const Array4<const Real>& detJ,
                        const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                        const Array4<const Real>& mf_m,
                        const AdvType horiz_adv_type,
                        const AdvType vert_adv_type,
                        const int use_terrain)
{
    BL_PROFILE_VAR("AdvectionSrcForScalars", AdvectionSrcForScalars);
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    // Inline with 2nd order for efficiency
    if (horiz_adv_type == AdvType::Centered_2nd && vert_adv_type == AdvType::Centered_2nd)
    {
        amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real invdetJ = (use_terrain) ?  1. / detJ(i,j,k) : 1.;

            // NOTE: we don't need to weight avg_xmom, avg_ymom, avg_zmom with terrain metrics
            //       because that was done when they were constructed in AdvectionSrcForRhoAndTheta

            const int cons_index = icomp + n;
            const int prim_index = cons_index - 1;

            Real mfsq = mf_m(i,j,0) * mf_m(i,j,0);

            advectionSrc(i,j,k,cons_index) = - 0.5 * invdetJ * (
              ( avg_xmom(i+1,j,k) * (cell_prim(i,j,k,prim_index) + cell_prim(i+1,j,k,prim_index)) -
                avg_xmom(i  ,j,k) * (cell_prim(i,j,k,prim_index) + cell_prim(i-1,j,k,prim_index)) ) * dxInv * mfsq +
              ( avg_ymom(i,j+1,k) * (cell_prim(i,j,k,prim_index) + cell_prim(i,j+1,k,prim_index)) -
                avg_ymom(i,j  ,k) * (cell_prim(i,j,k,prim_index) + cell_prim(i,j-1,k,prim_index)) ) * dyInv * mfsq +
              ( avg_zmom(i,j,k+1) * (cell_prim(i,j,k,prim_index) + cell_prim(i,j,k+1,prim_index)) -
                avg_zmom(i,j,k  ) * (cell_prim(i,j,k,prim_index) + cell_prim(i,j,k-1,prim_index)) ) * dzInv);
        });
    // Template higher order methods (horizontal first)
    } else {
        if (horiz_adv_type == AdvType::Centered_2nd) {
            AdvectionSrcForScalarsVert_N<CENTERED2>(bx, ncomp, icomp,
                                                    use_terrain, advectionSrc, cell_prim,
                                                    avg_xmom, avg_ymom, avg_zmom, detJ,
                                                    cellSizeInv, mf_m, vert_adv_type);
        } else if (horiz_adv_type == AdvType::Upwind_3rd) {
            AdvectionSrcForScalarsVert_N<UPWIND3>(bx, ncomp, icomp,
                                                  use_terrain, advectionSrc, cell_prim,
                                                  avg_xmom, avg_ymom, avg_zmom, detJ,
                                                  cellSizeInv, mf_m, vert_adv_type);
        } else if (horiz_adv_type == AdvType::Centered_4th) {
            AdvectionSrcForScalarsVert_N<CENTERED4>(bx, ncomp, icomp,
                                                    use_terrain, advectionSrc, cell_prim,
                                                    avg_xmom, avg_ymom, avg_zmom, detJ,
                                                    cellSizeInv, mf_m, vert_adv_type);
        } else if (horiz_adv_type == AdvType::Upwind_5th) {
            AdvectionSrcForScalarsVert_N<UPWIND5>(bx, ncomp, icomp,
                                                  use_terrain, advectionSrc, cell_prim,
                                                  avg_xmom, avg_ymom, avg_zmom, detJ,
                                                  cellSizeInv, mf_m, vert_adv_type);
        } else if (horiz_adv_type == AdvType::Centered_6th) {
            AdvectionSrcForScalarsVert_N<CENTERED6>(bx, ncomp, icomp,
                                                    use_terrain, advectionSrc, cell_prim,
                                                    avg_xmom, avg_ymom, avg_zmom, detJ,
                                                    cellSizeInv, mf_m, vert_adv_type);
        } else if (horiz_adv_type == AdvType::Weno_3) {
            AdvectionSrcForScalarsWrapper_N<WENO3,WENO3>(bx, ncomp, icomp,
                                                         use_terrain, advectionSrc, cell_prim,
                                                         avg_xmom, avg_ymom, avg_zmom, detJ,
                                                         cellSizeInv, mf_m);
        } else if (horiz_adv_type == AdvType::Weno_5) {
            AdvectionSrcForScalarsWrapper_N<WENO5,WENO5>(bx, ncomp, icomp,
                                                         use_terrain, advectionSrc, cell_prim,
                                                         avg_xmom, avg_ymom, avg_zmom, detJ,
                                                         cellSizeInv, mf_m);
        } else if (horiz_adv_type == AdvType::Weno_3Z) {
            AdvectionSrcForScalarsWrapper_N<WENO_Z3,WENO_Z3>(bx, ncomp, icomp,
                                                             use_terrain, advectionSrc, cell_prim,
                                                             avg_xmom, avg_ymom, avg_zmom, detJ,
                                                             cellSizeInv, mf_m);
        } else if (horiz_adv_type == AdvType::Weno_3MZQ) {
            AdvectionSrcForScalarsWrapper_N<WENO_MZQ3,WENO_MZQ3>(bx, ncomp, icomp,
                                                                 use_terrain, advectionSrc, cell_prim,
                                                                 avg_xmom, avg_ymom, avg_zmom, detJ,
                                                                 cellSizeInv, mf_m);
        } else if (horiz_adv_type == AdvType::Weno_5Z) {
            AdvectionSrcForScalarsWrapper_N<WENO_Z5,WENO_Z5>(bx, ncomp, icomp,
                                                             use_terrain, advectionSrc, cell_prim,
                                                             avg_xmom, avg_ymom, avg_zmom, detJ,
                                                             cellSizeInv, mf_m);
        } else {
            AMREX_ASSERT_WITH_MESSAGE(false, "Unknown advection scheme!");
        }
    }
}
