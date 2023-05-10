#include <AdvectionSrcForMom_N.H>
#include <AdvectionSrcForMom_T.H>

using namespace amrex;

/**
 * Function for computing the advective tendency for the momentum equations
 * This routine has explicit expressions for all cases (terrain or not) when
 * the horizontal and vertial spatial orders are <= 2, and calls more specialized
 * functions when either (or both) spatial order(s) is greater than 2.
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
 * @param[in] detJ Jacobian of the metric transformation (= 1 if use_terrain is false)
 * @param[in] cellSizeInv inverse of the mesh spacing
 * @param[in] mf_m map factor at cell centers
 * @param[in] mf_u map factor at x-faces
 * @param[in] mf_v map factor at y-faces
 * @param[in] all_use_WENO defines whether all variables (or just moisture variables) use WENO advection scheme
 * @param[in] spatial_order_WENO sets the spatial order if using WENO (3,5, or 7)
 * @param[in] horiz_spatial_order sets the spatial order to be used for lateral derivatives if not using WENO (2-6)
 * @param[in] vert_spatial_order sets the spatial order to be used for vertical derivatives if not using WENO (2-6)
 * @param[in] use_terrain if true, use the terrain-aware derivatives (with metric terms)
 * @param[in] domhi_z maximum k value in the domain
 */
void
AdvectionSrcForMom (const Box& bxx, const Box& bxy, const Box& bxz,
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
                    const Array4<const Real>& detJ,
                    const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                    const Array4<const Real>& mf_m,
                    const Array4<const Real>& mf_u,
                    const Array4<const Real>& mf_v,
                    const bool all_use_WENO,
                    const int spatial_order_WENO,
                    const int horiz_spatial_order,
                    const int vert_spatial_order,
                    const int use_terrain,
                    const int domhi_z)
{
    BL_PROFILE_VAR("AdvectionSrcForMom", AdvectionSrcForMom);

    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    AMREX_ALWAYS_ASSERT(bxz.smallEnd(2) > 0);

    // Not using terrain height coords
    if (!use_terrain) {
        // Directly compute RHS for second order for efficiency
        if (std::max(horiz_spatial_order,vert_spatial_order) == 2) {
            ParallelFor(bxx, bxy, bxz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real mf_u_inv_hi = 1. / mf_u(i+1,j  ,0); Real mf_u_inv_mid = 1. / mf_u(i  ,j  ,0);
                Real mf_u_inv_lo = 1. / mf_u(i-1,j  ,0);
                Real mf_v_inv_1  = 1. / mf_v(i  ,j+1,0); Real mf_v_inv_2   = 1. / mf_v(i-1,j+1,0);
                Real mf_v_inv_3  = 1. / mf_v(i  ,j  ,0); Real mf_v_inv_4   = 1. / mf_v(i-1,j  ,0);

                Real xflux_hi = 0.25 * (rho_u(i, j  , k) * mf_u_inv_mid + rho_u(i+1, j  , k) * mf_u_inv_hi) * (u(i+1,j,k) + u(i,j,k));
                Real xflux_lo = 0.25 * (rho_u(i, j  , k) * mf_u_inv_mid + rho_u(i-1, j  , k) * mf_u_inv_lo) * (u(i-1,j,k) + u(i,j,k));

                Real yflux_hi = 0.25 * (rho_v(i, j+1, k) * mf_v_inv_1 + rho_v(i-1, j+1, k) * mf_v_inv_2) * (u(i,j+1,k) + u(i,j,k));
                Real yflux_lo = 0.25 * (rho_v(i, j  , k) * mf_v_inv_3 + rho_v(i-1, j  , k) * mf_v_inv_4) * (u(i,j-1,k) + u(i,j,k));

                Real zflux_hi = 0.25 * (Omega(i, j, k+1) + Omega(i-1, j, k+1)) * (u(i,j,k+1) + u(i,j,k));
                Real zflux_lo = 0.25 * (Omega(i, j, k  ) + Omega(i-1, j, k  )) * (u(i,j,k-1) + u(i,j,k));

                Real mfsq = mf_u(i,j,0) * mf_u(i,j,0);

                Real advectionSrc = (xflux_hi - xflux_lo) * dxInv * mfsq
                                  + (yflux_hi - yflux_lo) * dyInv * mfsq
                                  + (zflux_hi - zflux_lo) * dzInv;
                rho_u_rhs(i, j, k) = -advectionSrc;
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real mf_v_inv_hi = 1. / mf_v(i  ,j+1,0); Real mf_v_inv_mid = 1. / mf_v(i  ,j  ,0);
                Real mf_v_inv_lo = 1. / mf_v(i  ,j-1,0);
                Real mf_u_inv_1  = 1. / mf_u(i+1,j  ,0); Real mf_u_inv_2   = 1. / mf_u(i+1,j-1,0);
                Real mf_u_inv_3  = 1. / mf_u(i  ,j  ,0); Real mf_u_inv_4   = 1. / mf_u(i  ,j-1,0);

                Real xflux_hi = 0.25 * (rho_u(i+1, j, k) * mf_u_inv_1 + rho_u(i+1, j-1, k) * mf_u_inv_2) * (v(i+1,j,k) + v(i,j,k));
                Real xflux_lo = 0.25 * (rho_u(i  , j, k) * mf_u_inv_3 + rho_u(i  , j-1, k) * mf_u_inv_4) * (v(i-1,j,k) + v(i,j,k));

                Real yflux_hi = 0.25 * (rho_v(i  ,j+1,k) * mf_v_inv_hi  + rho_v(i  ,j  ,k) * mf_v_inv_mid) * (v(i,j+1,k) + v(i,j,k));
                Real yflux_lo = 0.25 * (rho_v(i  ,j  ,k) * mf_v_inv_mid + rho_v(i  ,j-1,k) * mf_v_inv_lo ) * (v(i,j-1,k) + v(i,j,k));

                Real zflux_hi = 0.25 * (Omega(i, j, k+1) + Omega(i, j-1, k+1)) * (v(i,j,k+1) + v(i,j,k));
                Real zflux_lo = 0.25 * (Omega(i, j, k  ) + Omega(i, j-1, k  )) * (v(i,j,k-1) + v(i,j,k));

                Real mfsq = mf_v(i,j,0) * mf_v(i,j,0);

                Real advectionSrc = (xflux_hi - xflux_lo) * dxInv * mfsq
                                  + (yflux_hi - yflux_lo) * dyInv * mfsq
                                  + (zflux_hi - zflux_lo) * dzInv;
                rho_v_rhs(i, j, k) = -advectionSrc;
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real mf_u_inv_hi = 1. / mf_u(i+1,j  ,0); Real mf_u_inv_lo = 1. / mf_u(i  ,j  ,0);
                Real mf_v_inv_hi = 1. / mf_v(i  ,j+1,0); Real mf_v_inv_lo = 1. / mf_v(i  ,j  ,0);

                Real xflux_hi = 0.25*(rho_u(i+1,j  ,k) + rho_u(i+1, j, k-1)) * mf_u_inv_hi * (w(i+1,j,k) + w(i,j,k));
                Real xflux_lo = 0.25*(rho_u(i  ,j  ,k) + rho_u(i  , j, k-1)) * mf_u_inv_lo * (w(i-1,j,k) + w(i,j,k));

                Real yflux_hi = 0.25*(rho_v(i  ,j+1,k) + rho_v(i, j+1, k-1)) * mf_v_inv_hi * (w(i,j+1,k) + w(i,j,k));
                Real yflux_lo = 0.25*(rho_v(i  ,j  ,k) + rho_v(i, j  , k-1)) * mf_v_inv_lo * (w(i,j-1,k) + w(i,j,k));

                Real zflux_lo = 0.25 * (Omega(i,j,k) + Omega(i,j,k-1)) * (w(i,j,k) + w(i,j,k-1));

                Real zflux_hi = (k == domhi_z+1) ? Omega(i,j,k) * w(i,j,k) :
                    0.25 * (Omega(i,j,k) + Omega(i,j,k+1)) * (w(i,j,k) + w(i,j,k+1));

                Real mfsq = mf_m(i,j,0) * mf_m(i,j,0);

                Real advectionSrc = (xflux_hi - xflux_lo) * dxInv * mfsq
                                  + (yflux_hi - yflux_lo) * dyInv * mfsq
                                  + (zflux_hi - zflux_lo) * dzInv;
                rho_w_rhs(i, j, k) = -advectionSrc;
            });

        // Get struct for higher order interpolation
        } else {
            if (std::max(horiz_spatial_order,vert_spatial_order) == 3) {
                AdvectionSrcForMomWrapper_N<UPWIND3, UPWINDALL>(bxx, bxy, bxz,
                                                                rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                                rho_u, rho_v, Omega, u, v, w,
                                                                cellSizeInv, mf_m, mf_u, mf_v,
                                                                vert_spatial_order, domhi_z);

            } else if (std::max(horiz_spatial_order,vert_spatial_order) == 4) {
                AdvectionSrcForMomWrapper_N<UPWIND4, UPWINDALL>(bxx, bxy, bxz,
                                                                rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                                rho_u, rho_v, Omega, u, v, w,
                                                                cellSizeInv, mf_m, mf_u, mf_v,
                                                                vert_spatial_order, domhi_z);
            } else if (std::max(horiz_spatial_order,vert_spatial_order) == 5) {
                AdvectionSrcForMomWrapper_N<UPWIND5, UPWINDALL>(bxx, bxy, bxz,
                                                                rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                                rho_u, rho_v, Omega, u, v, w,
                                                                cellSizeInv, mf_m, mf_u, mf_v,
                                                                vert_spatial_order, domhi_z);
            } else if (std::max(horiz_spatial_order,vert_spatial_order) == 6) {
                AdvectionSrcForMomWrapper_N<UPWIND6, UPWINDALL>(bxx, bxy, bxz,
                                                                rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                                rho_u, rho_v, Omega, u, v, w,
                                                                cellSizeInv, mf_m, mf_u, mf_v,
                                                                vert_spatial_order, domhi_z);
            } else if (all_use_WENO && spatial_order_WENO==3) {
                AdvectionSrcForMomWrapper_N<WENO3, UPWINDALL>(bxx, bxy, bxz,
                                                              rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                              rho_u, rho_v, Omega, u, v, w,
                                                              cellSizeInv, mf_m, mf_u, mf_v,
                                                              vert_spatial_order, domhi_z);
            } else if (all_use_WENO && spatial_order_WENO==5) {
                AdvectionSrcForMomWrapper_N<WENO5, UPWINDALL>(bxx, bxy, bxz,
                                                              rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                              rho_u, rho_v, Omega, u, v, w,
                                                              cellSizeInv, mf_m, mf_u, mf_v,
                                                              vert_spatial_order, domhi_z);

            } else {
                AMREX_ASSERT_WITH_MESSAGE(false, "Unknown advection scheme!");
            }
        }

    // Using terrain height coordinates
    } else {
        // Directly compute RHS for second order for efficiency
        if (std::max(horiz_spatial_order,vert_spatial_order) == 2) {
            ParallelFor(bxx, bxy, bxz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real mf_u_inv_hi  = 1. / mf_u(i+1,j  ,0);
                Real mf_u_inv_mid = 1. / mf_u(i  ,j  ,0);
                Real mf_u_inv_lo  = 1. / mf_u(i-1,j  ,0);
                Real mf_v_inv_1   = 1. / mf_v(i  ,j+1,0); Real mf_v_inv_2  = 1. / mf_v(i-1,j+1,0);
                Real mf_v_inv_3   = 1. / mf_v(i  ,j  ,0); Real mf_v_inv_4  = 1. / mf_v(i-1,j  ,0);

                Real met_h_zeta_xhi = Compute_h_zeta_AtCellCenter(i  ,j  ,k  ,cellSizeInv,z_nd);
                Real xflux_hi = 0.25 * (rho_u(i, j  , k) * mf_u_inv_mid + rho_u(i+1, j  , k) * mf_u_inv_hi) * (u(i+1,j,k) + u(i,j,k)) * met_h_zeta_xhi;

                Real met_h_zeta_xlo = Compute_h_zeta_AtCellCenter(i-1,j  ,k  ,cellSizeInv,z_nd);
                Real xflux_lo = 0.25 * (rho_u(i, j  , k) * mf_u_inv_mid + rho_u(i-1, j  , k) * mf_u_inv_lo) * (u(i-1,j,k) + u(i,j,k)) * met_h_zeta_xlo;

                Real met_h_zeta_yhi = Compute_h_zeta_AtEdgeCenterK(i  ,j+1,k  ,cellSizeInv,z_nd);
                Real yflux_hi = 0.25 * (rho_v(i, j+1, k) * mf_v_inv_1 + rho_v(i-1, j+1, k) * mf_v_inv_2) * (u(i,j+1,k) + u(i,j,k)) * met_h_zeta_yhi;

                Real met_h_zeta_ylo = Compute_h_zeta_AtEdgeCenterK(i  ,j  ,k  ,cellSizeInv,z_nd);
                Real yflux_lo = 0.25 * (rho_v(i, j  , k) * mf_v_inv_3 + rho_v(i-1, j  , k) * mf_v_inv_4) * (u(i,j-1,k) + u(i,j,k)) * met_h_zeta_ylo;

                Real zflux_hi = 0.25 * (Omega(i, j, k+1) + Omega(i-1, j, k+1)) * (u(i,j,k+1) + u(i,j,k));
                Real zflux_lo = 0.25 * (Omega(i, j, k  ) + Omega(i-1, j, k  )) * (u(i,j,k-1) + u(i,j,k));

                Real mfsq = mf_u(i,j,0) * mf_u(i,j,0);

                Real advectionSrc = (xflux_hi - xflux_lo) * dxInv * mfsq
                                  + (yflux_hi - yflux_lo) * dyInv * mfsq
                                  + (zflux_hi - zflux_lo) * dzInv;

                rho_u_rhs(i, j, k) = -advectionSrc / (0.5 * (detJ(i,j,k) + detJ(i-1,j,k)));
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real mf_v_inv_hi = 1. / mf_v(i  ,j+1,0); Real mf_v_inv_mid = 1. / mf_v(i  ,j  ,0); Real mf_v_inv_lo = 1. / mf_v(i  ,j-1,0);
                Real mf_u_inv_1  = 1. / mf_u(i+1,j  ,0); Real mf_u_inv_2   = 1. / mf_u(i+1,j-1,0); Real mf_u_inv_3  = 1. / mf_u(i  ,j  ,0); Real mf_u_inv_4 = 1. / mf_u(i-1,j  ,0);

                Real met_h_zeta_xhi = Compute_h_zeta_AtEdgeCenterK(i+1,j  ,k  ,cellSizeInv,z_nd);
                Real xflux_hi = 0.25 * (rho_u(i+1,j  ,k) * mf_u_inv_1 + rho_u(i+1,j-1, k) * mf_u_inv_2) * (v(i+1,j,k) + v(i,j,k)) * met_h_zeta_xhi;

                Real met_h_zeta_xlo = Compute_h_zeta_AtEdgeCenterK(i  ,j  ,k  ,cellSizeInv,z_nd);
                Real xflux_lo = 0.25 * (rho_u(i, j  , k) * mf_u_inv_3 + rho_u(i  ,j-1, k) * mf_u_inv_4) * (v(i-1,j,k) + v(i,j,k)) * met_h_zeta_xlo;

                Real met_h_zeta_yhi = Compute_h_zeta_AtCellCenter(i  ,j  ,k  ,cellSizeInv,z_nd);
                Real yflux_hi = 0.25 * (rho_v(i  ,j+1, k) * mf_v_inv_hi + rho_v(i  ,j  ,k) * mf_v_inv_mid) * (v(i,j+1,k) + v(i,j,k)) * met_h_zeta_yhi;

                Real met_h_zeta_ylo = Compute_h_zeta_AtCellCenter(i  ,j-1,k  ,cellSizeInv,z_nd);
                Real yflux_lo = 0.25 * (rho_v(i  ,j  ,k) * mf_v_inv_mid + rho_v(i  , j-1, k) * mf_v_inv_lo) * (v(i,j-1,k) + v(i,j,k)) * met_h_zeta_ylo;

                Real zflux_hi = 0.25 * (Omega(i, j, k+1) + Omega(i, j-1, k+1)) * (v(i,j,k+1) + v(i,j,k));
                Real zflux_lo = 0.25 * (Omega(i, j, k  ) + Omega(i, j-1, k  )) * (v(i,j,k-1) + v(i,j,k));

                Real mfsq = mf_v(i,j,0) * mf_v(i,j,0);

                Real advectionSrc = (xflux_hi - xflux_lo) * dxInv * mfsq
                                  + (yflux_hi - yflux_lo) * dyInv * mfsq
                                  + (zflux_hi - zflux_lo) * dzInv;

                rho_v_rhs(i, j, k) = -advectionSrc / (0.5 * (detJ(i,j,k) + detJ(i,j-1,k)));
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real mf_u_inv_hi = 1. / mf_u(i+1,j  ,0); Real mf_u_inv_lo = 1. / mf_u(i  ,j  ,0);
                Real mf_v_inv_hi = 1. / mf_v(i  ,j+1,0); Real mf_v_inv_lo = 1. / mf_v(i  ,j  ,0);

                Real met_h_zeta_xhi = Compute_h_zeta_AtEdgeCenterJ(i+1,j  ,k  ,cellSizeInv,z_nd);
                Real xflux_hi = 0.25*(rho_u(i+1,j  ,k) + rho_u(i+1, j, k-1)) * mf_u_inv_hi * (w(i+1,j,k) + w(i,j,k)) * met_h_zeta_xhi;

                Real met_h_zeta_xlo = Compute_h_zeta_AtEdgeCenterJ(i  ,j  ,k  ,cellSizeInv,z_nd);
                Real xflux_lo = 0.25*(rho_u(i  ,j  ,k) + rho_u(i  , j, k-1)) * mf_u_inv_lo * (w(i-1,j,k) + w(i,j,k)) * met_h_zeta_xlo;

                Real met_h_zeta_yhi = Compute_h_zeta_AtEdgeCenterI(i  ,j+1,k  ,cellSizeInv,z_nd);
                Real yflux_hi = 0.25*(rho_v(i  ,j+1,k) + rho_v(i, j+1, k-1)) * mf_v_inv_hi * (w(i,j+1,k) + w(i,j,k)) * met_h_zeta_yhi;

                Real met_h_zeta_ylo = Compute_h_zeta_AtEdgeCenterI(i  ,j  ,k  ,cellSizeInv,z_nd);
                Real yflux_lo = 0.25*(rho_v(i  ,j  ,k) + rho_v(i, j  , k-1)) * mf_v_inv_lo * (w(i,j-1,k) + w(i,j,k)) * met_h_zeta_ylo;

                Real zflux_lo = 0.25 * (Omega(i,j,k) + Omega(i,j,k-1)) * (w(i,j,k) + w(i,j,k-1));

                Real zflux_hi = (k == domhi_z+1) ? Omega(i,j,k) * w(i,j,k) :
                    0.25 * (Omega(i,j,k) + Omega(i,j,k+1)) * (w(i,j,k) + w(i,j,k+1));

                Real mfsq = mf_m(i,j,0) * mf_m(i,j,0);

                Real advectionSrc = (xflux_hi - xflux_lo) * dxInv * mfsq
                                  + (yflux_hi - yflux_lo) * dyInv * mfsq
                                  + (zflux_hi - zflux_lo) * dzInv;

                rho_w_rhs(i, j, k) = -advectionSrc / (0.5*(detJ(i,j,k) + detJ(i,j,k-1)));
            });

        // Get struct for higher order interpolation
        } else {
            if (std::max(horiz_spatial_order,vert_spatial_order) == 3) {
                AdvectionSrcForMomWrapper_T<UPWIND3, UPWINDALL>(bxx, bxy, bxz,
                                                                rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                                rho_u, rho_v, Omega, u, v, w, z_nd, detJ,
                                                                cellSizeInv, mf_m, mf_u, mf_v,
                                                                vert_spatial_order, domhi_z);
            } else if (std::max(horiz_spatial_order,vert_spatial_order) == 4) {
                AdvectionSrcForMomWrapper_T<UPWIND4, UPWINDALL>(bxx, bxy, bxz,
                                                                rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                                rho_u, rho_v, Omega, u, v, w, z_nd, detJ,
                                                                cellSizeInv, mf_m, mf_u, mf_v,
                                                                vert_spatial_order, domhi_z);
            } else if (std::max(horiz_spatial_order,vert_spatial_order) == 5) {
                AdvectionSrcForMomWrapper_T<UPWIND5, UPWINDALL>(bxx, bxy, bxz,
                                                                rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                                rho_u, rho_v, Omega, u, v, w, z_nd, detJ,
                                                                cellSizeInv, mf_m, mf_u, mf_v,
                                                                vert_spatial_order, domhi_z);
            } else if (std::max(horiz_spatial_order,vert_spatial_order) == 6) {
                AdvectionSrcForMomWrapper_T<UPWIND6, UPWINDALL>(bxx, bxy, bxz,
                                                                rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                                rho_u, rho_v, Omega, u, v, w, z_nd, detJ,
                                                                cellSizeInv, mf_m, mf_u, mf_v,
                                                                vert_spatial_order, domhi_z);
            } else if (all_use_WENO && spatial_order_WENO==3) {
                AdvectionSrcForMomWrapper_T<WENO3, UPWINDALL>(bxx, bxy, bxz,
                                                              rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                              rho_u, rho_v, Omega, u, v, w, z_nd, detJ,
                                                              cellSizeInv, mf_m, mf_u, mf_v,
                                                              vert_spatial_order, domhi_z);
            } else if (all_use_WENO && spatial_order_WENO==5) {
                AdvectionSrcForMomWrapper_T<WENO5, UPWINDALL>(bxx, bxy, bxz,
                                                              rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                              rho_u, rho_v, Omega, u, v, w, z_nd, detJ,
                                                              cellSizeInv, mf_m, mf_u, mf_v,
                                                              vert_spatial_order, domhi_z);
            } else {
                AMREX_ASSERT_WITH_MESSAGE(false, "Unknown advection scheme!");
            }
        }
    }
}

