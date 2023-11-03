#include <IndexDefines.H>
#include <TerrainMetrics.H>
#include <Interpolation.H>

/**
 * Function for computing the advective tendency for the x-component of momentum
 * with metric terms and for spatial order > 2
 *
 * @param[in] i,j,k indices of x-face at which to create tendency
 * @param[in] rho_u x-component of momentum
 * @param[in] rho_v y-component of momentum
 * @param[in] Omega component of the momentum normal to the z-coordinate surface
 * @param[in] u     x-component of velocity
 * @param[in] z_nd height coordinate at nodes
 * @param[in] detJ Jacobian of the metric transformation (= 1 if use_terrain is false)
 * @param[in] cellSizeInv inverse of the mesh spacing
 * @param[in] mf_u map factor on x-faces
 * @param[in] mf_v map factor on y-faces
 */
template<typename InterpType_H, typename InterpType_V>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForXMom_T (int i, int j, int k,
                       const amrex::Array4<const amrex::Real>& rho_u,
                       const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& Omega,
                       const amrex::Array4<const amrex::Real>& z_nd,
                       const amrex::Array4<const amrex::Real>& detJ,
                       InterpType_H interp_u_h,
                       InterpType_V interp_u_v,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Array4<const amrex::Real>& mf_u_inv,
                       const amrex::Array4<const amrex::Real>& mf_v_inv)
{
    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    amrex::Real rho_u_avg_lo, rho_u_avg_hi;
    amrex::Real rho_v_avg_lo, rho_v_avg_hi;
    amrex::Real Omega_avg_lo, Omega_avg_hi;

    amrex::Real interp_hi(0.), interp_lo(0.);

    // ****************************************************************************************
    // X-fluxes (at cell centers)
    // ****************************************************************************************

    rho_u_avg_hi = 0.5 * (rho_u(i+1, j, k) * mf_u_inv(i+1,j  ,0) + rho_u(i, j, k) * mf_u_inv(i  ,j  ,0));
    interp_u_h.InterpolateInX(i,j,k+1,0,interp_hi,rho_u_avg_hi);

    rho_u_avg_lo = 0.5 * (rho_u(i-1, j, k) * mf_u_inv(i-1,j  ,0) + rho_u(i, j, k) * mf_u_inv(i  ,j  ,0));
    interp_u_h.InterpolateInX(i,j,k,0,interp_lo,rho_u_avg_lo);

    amrex::Real centFluxXXNext = rho_u_avg_hi * Compute_h_zeta_AtCellCenter(i  ,j  ,k  ,cellSizeInv,z_nd) * interp_hi;
    amrex::Real centFluxXXPrev = rho_u_avg_lo * Compute_h_zeta_AtCellCenter(i-1,j  ,k  ,cellSizeInv,z_nd) * interp_lo;

    // ****************************************************************************************
    // Y-fluxes (at edges in k-direction)
    // ****************************************************************************************
    rho_v_avg_hi = 0.5 * (rho_v(i, j+1, k) * mf_v_inv(i  ,j+1,0) + rho_v(i-1, j+1, k) * mf_v_inv(i-1,j+1,0));
    interp_u_h.InterpolateInY(i,j+1,k,0,interp_hi,rho_v_avg_hi);

    rho_v_avg_lo = 0.5 * (rho_v(i, j  , k) * mf_v_inv(i  ,j  ,0) + rho_v(i-1, j  , k) * mf_v_inv(i-1,j  ,0));
    interp_u_h.InterpolateInY(i,j,k,0,interp_lo,rho_v_avg_lo);

    amrex::Real edgeFluxXYNext = rho_v_avg_hi * Compute_h_zeta_AtEdgeCenterK(i  ,j+1,k  ,cellSizeInv,z_nd) * interp_hi;
    amrex::Real edgeFluxXYPrev = rho_v_avg_lo * Compute_h_zeta_AtEdgeCenterK(i  ,j  ,k  ,cellSizeInv,z_nd) * interp_lo;

    // ****************************************************************************************
    // Z-fluxes (at edges in j-direction)
    // ****************************************************************************************
    Omega_avg_hi = 0.5 * (Omega(i, j, k+1) + Omega(i-1, j, k+1));
    Omega_avg_lo = 0.5 * (Omega(i, j, k  ) + Omega(i-1, j, k  ));

    interp_u_v.InterpolateInZ(i,j,k+1,0,interp_hi,Omega_avg_hi);
    interp_u_v.InterpolateInZ(i,j,k  ,0,interp_lo,Omega_avg_lo);

    amrex::Real edgeFluxXZNext = Omega_avg_hi * interp_hi;
    amrex::Real edgeFluxXZPrev = Omega_avg_lo * interp_lo;

    // ****************************************************************************************

    amrex::Real mfsq = 1 / (mf_u_inv(i,j,0) * mf_u_inv(i,j,0));

    advectionSrc = (centFluxXXNext - centFluxXXPrev) * dxInv * mfsq
                 + (edgeFluxXYNext - edgeFluxXYPrev) * dyInv * mfsq
                 + (edgeFluxXZNext - edgeFluxXZPrev) * dzInv;
    advectionSrc /= 0.5*(detJ(i,j,k) + detJ(i-1,j,k));

    return advectionSrc;
}

/**
 * Function for computing the advective tendency for the y-component of momentum
 * with metric terms and for spatial order > 2
 *
 * @param[in] i,j,k indices of y-face at which to create tendency
 * @param[in] rho_u x-component of momentum
 * @param[in] rho_v y-component of momentum
 * @param[in] Omega component of the momentum normal to the z-coordinate surface
 * @param[in] v     y-component of velocity
 * @param[in] z_nd height coordinate at nodes
 * @param[in] detJ Jacobian of the metric transformation (= 1 if use_terrain is false)
 * @param[in] cellSizeInv inverse of the mesh spacing
 * @param[in] mf_u map factor on x-faces
 * @param[in] mf_v map factor on y-faces
 */
template<typename InterpType_H, typename InterpType_V>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForYMom_T (int i, int j, int k,
                       const amrex::Array4<const amrex::Real>& rho_u,
                       const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& Omega,
                       const amrex::Array4<const amrex::Real>& z_nd,
                       const amrex::Array4<const amrex::Real>& detJ,
                       InterpType_H interp_v_h,
                       InterpType_V interp_v_v,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Array4<const amrex::Real>& mf_u_inv,
                       const amrex::Array4<const amrex::Real>& mf_v_inv)
{
    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    amrex::Real rho_u_avg_lo, rho_u_avg_hi;
    amrex::Real rho_v_avg_lo, rho_v_avg_hi;
    amrex::Real Omega_avg_lo, Omega_avg_hi;

    amrex::Real interp_hi(0.), interp_lo(0.);

    // ****************************************************************************************
    // x-fluxes (at edges in k-direction)
    // ****************************************************************************************
    rho_u_avg_hi = 0.5 * (rho_u(i+1, j, k) * mf_u_inv(i+1,j  ,0) + rho_u(i+1, j-1, k) * mf_u_inv(i+1,j-1,0));
    interp_v_h.InterpolateInX(i+1,j,k,0,interp_hi,rho_u_avg_hi);

    rho_u_avg_lo = 0.5 * (rho_u(i  , j, k) * mf_u_inv(i  ,j  ,0) + rho_u(i  , j-1, k) * mf_u_inv(i  ,j-1,0));
    interp_v_h.InterpolateInX(i,j,k,0,interp_lo,rho_u_avg_lo);

    amrex::Real edgeFluxYXNext = rho_u_avg_hi * Compute_h_zeta_AtEdgeCenterK(i+1,j  ,k  ,cellSizeInv,z_nd) * interp_hi;
    amrex::Real edgeFluxYXPrev = rho_u_avg_lo * Compute_h_zeta_AtEdgeCenterK(i  ,j  ,k  ,cellSizeInv,z_nd) * interp_lo;

    // ****************************************************************************************
    // y-fluxes (at cell centers)
    // ****************************************************************************************
    rho_v_avg_hi = 0.5 * (rho_v(i, j, k) * mf_v_inv(i  ,j  ,0) + rho_v(i, j+1, k) * mf_v_inv(i  ,j+1,0));
    interp_v_h.InterpolateInY(i,j+1,k,0,interp_hi,rho_v_avg_hi);

    rho_v_avg_lo = 0.5 * (rho_v(i, j, k) * mf_v_inv(i  ,j  ,0) + rho_v(i, j-1, k) * mf_v_inv(i  ,j-1,0));
    interp_v_h.InterpolateInY(i,j,k,0,interp_lo,rho_v_avg_lo);

    amrex::Real centFluxYYNext = rho_v_avg_hi * Compute_h_zeta_AtCellCenter(i  ,j  ,k  ,cellSizeInv,z_nd) * interp_hi;
    amrex::Real centFluxYYPrev = rho_v_avg_lo * Compute_h_zeta_AtCellCenter(i  ,j-1,k  ,cellSizeInv,z_nd) * interp_lo;

    // ****************************************************************************************
    // Z-fluxes (at edges in j-direction)
    // ****************************************************************************************
    Omega_avg_hi = 0.5 * (Omega(i, j, k+1) + Omega(i, j-1, k+1));
    Omega_avg_lo = 0.5 * (Omega(i, j, k  ) + Omega(i, j-1, k  ));

    interp_v_v.InterpolateInZ(i,j,k+1,0,interp_hi,Omega_avg_hi);
    interp_v_v.InterpolateInZ(i,j,k  ,0,interp_lo,Omega_avg_lo);

    amrex::Real edgeFluxYZNext = Omega_avg_hi * interp_hi;
    amrex::Real edgeFluxYZPrev = Omega_avg_lo * interp_lo;

    // ****************************************************************************************

    amrex::Real mfsq = 1 / (mf_v_inv(i,j,0) * mf_v_inv(i,j,0));

    advectionSrc = (edgeFluxYXNext - edgeFluxYXPrev) * dxInv * mfsq
                 + (centFluxYYNext - centFluxYYPrev) * dyInv * mfsq
                 + (edgeFluxYZNext - edgeFluxYZPrev) * dzInv;
    advectionSrc /= 0.5*(detJ(i,j,k) + detJ(i,j-1,k));

    return advectionSrc;
}

/**
 * Function for computing the advective tendency for the z-component of momentum
 * with metric terms and for spatial order > 2
 *
 * @param[in] i,j,k indices of z-face at which to create tendency
 * @param[in] rho_u x-component of momentum
 * @param[in] rho_v y-component of momentum
 * @param[in] Omega component of the momentum normal to the z-coordinate surface
 * @param[in] w     z-component of velocity
 * @param[in] z_nd height coordinate at nodes
 * @param[in] detJ Jacobian of the metric transformation (= 1 if use_terrain is false)
 * @param[in] cellSizeInv inverse of the mesh spacing
 * @param[in] mf_m map factor on cell centers
 * @param[in] mf_u map factor on x-faces
 * @param[in] mf_v map factor on y-faces
 * @param[in] vert_adv_type int that defines advection stencil
 * @param[in] domhi_z maximum k value in the domain
 */
template<typename InterpType_H, typename InterpType_V, typename WallInterpType>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForZMom_T (int i, int j, int k,
                       const amrex::Array4<const amrex::Real>& rho_u,
                       const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& Omega,
                       const amrex::Array4<const amrex::Real>& w,
                       const amrex::Array4<const amrex::Real>& z_nd,
                       const amrex::Array4<const amrex::Real>& detJ,
                       InterpType_H interp_omega_h,
                       InterpType_V interp_omega_v,
                       WallInterpType interp_omega_wall,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Array4<const amrex::Real>& mf_m,
                       const amrex::Array4<const amrex::Real>& mf_u_inv,
                       const amrex::Array4<const amrex::Real>& mf_v_inv,
                       const AdvType vert_adv_type,
                       const int domhi_z)
{
    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    amrex::Real rho_u_avg_lo, rho_u_avg_hi;
    amrex::Real rho_v_avg_lo, rho_v_avg_hi;
    amrex::Real Omega_avg_lo, Omega_avg_hi;

    amrex::Real interp_hi(0.), interp_lo(0.);

    // ****************************************************************************************
    // x-fluxes (at edges in j-direction)
    // ****************************************************************************************
    rho_u_avg_hi = 0.5 * (rho_u(i+1, j, k) + rho_u(i+1, j, k-1)) * mf_u_inv(i+1,j  ,0);
    interp_omega_h.InterpolateInX(i+1,j,k,0,interp_hi,rho_u_avg_hi);

    rho_u_avg_lo = 0.5 * (rho_u(i  , j, k) + rho_u(i  , j, k-1)) * mf_u_inv(i  ,j  ,0);
    interp_omega_h.InterpolateInX(i,j,k,0,interp_lo,rho_u_avg_lo);

    amrex::Real edgeFluxZXNext = rho_u_avg_hi * Compute_h_zeta_AtEdgeCenterJ(i+1,j  ,k  ,cellSizeInv,z_nd) * interp_hi;
    amrex::Real edgeFluxZXPrev = rho_u_avg_lo * Compute_h_zeta_AtEdgeCenterJ(i  ,j  ,k  ,cellSizeInv,z_nd) * interp_lo;

    // ****************************************************************************************
    // y-fluxes (at edges in i-direction)
    // ****************************************************************************************
    rho_v_avg_hi = 0.5 * (rho_v(i, j+1, k) + rho_v(i, j+1, k-1)) * mf_v_inv(i  ,j+1,0);
    interp_omega_h.InterpolateInY(i,j+1,k,0,interp_hi,rho_v_avg_hi);

    rho_v_avg_lo = 0.5 * (rho_v(i, j  , k) + rho_v(i, j  , k-1)) * mf_v_inv(i  ,j  ,0);
    interp_omega_h.InterpolateInY(i,j,k,0,interp_lo,rho_v_avg_lo);

    amrex::Real edgeFluxZYNext = rho_v_avg_hi * Compute_h_zeta_AtEdgeCenterI(i  ,j+1,k  ,cellSizeInv,z_nd) * interp_hi;
    amrex::Real edgeFluxZYPrev = rho_v_avg_lo * Compute_h_zeta_AtEdgeCenterI(i  ,j  ,k  ,cellSizeInv,z_nd) * interp_lo;

    // ****************************************************************************************
    // z-fluxes (at cell centers)
    // ****************************************************************************************

    Omega_avg_hi = (k == domhi_z+1) ? Omega(i,j,k) : 0.5 * (Omega(i,j,k) + Omega(i,j,k+1));
    amrex::Real centFluxZZNext = Omega_avg_hi;

    // int l_spatial_order_hi = std::min(std::min(vert_spatial_order, 2*(domhi_z+1-k)), 2*(k+1));
    // If k == domhi_z   , l_spatial_order_hi = 2
    // If k == domhi_z-1 , l_spatial_order_hi = std::min(vert_spatial_order, 4);
    // If k ==         1 , l_spatial_order_hi = std::min(vert_spatial_order, 4);
    if (k == domhi_z+1) {
        centFluxZZNext *=  w(i,j,k);
    } else {
        if (k == domhi_z) {
            interp_omega_wall.InterpolateInZ(i,j,k+1,0,interp_hi,Omega_avg_hi,AdvType::Centered_2nd);
        } else if (k == domhi_z-1 || k == 1) {
            if (vert_adv_type != AdvType::Centered_2nd && vert_adv_type != AdvType::Upwind_3rd) {
                interp_omega_wall.InterpolateInZ(i,j,k+1,0,interp_hi,Omega_avg_hi,AdvType::Centered_4th);
            } else {
                interp_omega_wall.InterpolateInZ(i,j,k+1,0,interp_hi,Omega_avg_hi,vert_adv_type);
            }
        } else {
            interp_omega_v.InterpolateInZ(i,j,k+1,0,interp_hi,Omega_avg_hi);
        }
        centFluxZZNext *=  interp_hi;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    Omega_avg_lo = (k == 0) ? Omega(i,j,k) : 0.5 * (Omega(i,j,k) + Omega(i,j,k-1));
    amrex::Real centFluxZZPrev = Omega_avg_lo;

    // int l_spatial_order_lo = std::min(std::min(vert_spatial_order, 2*(domhi_z+2-k)), 2*k);
    // If k == domhi_z   , l_spatial_order_hi = 2
    // If k == domhi_z-1 , l_spatial_order_hi = std::min(vert_spatial_order, 4);
    // If k ==         1 , l_spatial_order_hi = std::min(vert_spatial_order, 4);
    if (k == 0) {
        centFluxZZPrev *=  w(i,j,k);
    } else {
        if (k == 1) {
            interp_omega_wall.InterpolateInZ(i,j,k,0,interp_lo,Omega_avg_lo,AdvType::Centered_2nd);
        } else if (k == 2 || k == domhi_z) {
            if (vert_adv_type != AdvType::Centered_2nd && vert_adv_type != AdvType::Upwind_3rd) {
                interp_omega_wall.InterpolateInZ(i,j,k,0,interp_lo,Omega_avg_lo,AdvType::Centered_4th);
            } else {
                interp_omega_wall.InterpolateInZ(i,j,k,0,interp_lo,Omega_avg_lo,vert_adv_type);
            }
        } else {
            interp_omega_v.InterpolateInZ(i,j,k,0,interp_lo,Omega_avg_lo);
        }
        centFluxZZPrev *=  interp_lo;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    amrex::Real mfsq = mf_m(i,j,0) * mf_m(i,j,0);

    advectionSrc = (edgeFluxZXNext - edgeFluxZXPrev) * dxInv * mfsq
                 + (edgeFluxZYNext - edgeFluxZYPrev) * dyInv * mfsq
                 + (centFluxZZNext - centFluxZZPrev) * dzInv;

    amrex::Real denom = 0.5*(detJ(i,j,k) + detJ(i,j,k-1));
    advectionSrc /= denom;

    return advectionSrc;
}

/**
 * Wrapper function for computing the advective tendency w/ spatial order > 2.
 */
template<typename InterpType_H, typename InterpType_V, typename WallInterpType>
void
AdvectionSrcForMomWrapper_T(const amrex::Box& bxx, const amrex::Box& bxy, const amrex::Box& bxz,
                            const amrex::Array4<amrex::Real>& rho_u_rhs,
                            const amrex::Array4<amrex::Real>& rho_v_rhs,
                            const amrex::Array4<amrex::Real>& rho_w_rhs,
                            const amrex::Array4<const amrex::Real>& rho_u,
                            const amrex::Array4<const amrex::Real>& rho_v,
                            const amrex::Array4<const amrex::Real>& Omega,
                            const amrex::Array4<const amrex::Real>& u,
                            const amrex::Array4<const amrex::Real>& v,
                            const amrex::Array4<const amrex::Real>& w,
                            const amrex::Array4<const amrex::Real>& z_nd,
                            const amrex::Array4<const amrex::Real>& detJ,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                            const amrex::Array4<const amrex::Real>& mf_m,
                            const amrex::Array4<const amrex::Real>& mf_u_inv,
                            const amrex::Array4<const amrex::Real>& mf_v_inv,
                            const AdvType vert_adv_type,
                            const int domhi_z)
{
    // Instantiate the appropriate structs
    InterpType_H interp_u_h(u); InterpType_V interp_u_v(u); // X-MOM
    InterpType_H interp_v_h(v); InterpType_V interp_v_v(v); // Y-MOM
    InterpType_H interp_w_h(w); InterpType_V interp_w_v(w); // Z-MOM
    WallInterpType interp_w_wall(w); // Z-MOM @ wall

    amrex::ParallelFor(bxx, bxy, bxz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        rho_u_rhs(i, j, k) = -AdvectionSrcForXMom_T(i, j, k, rho_u, rho_v, Omega, z_nd, detJ,
                                                    interp_u_h, interp_u_v, cellSizeInv, mf_u_inv, mf_v_inv);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        rho_v_rhs(i, j, k) = -AdvectionSrcForYMom_T(i, j, k, rho_u, rho_v, Omega, z_nd, detJ,
                                                    interp_v_h, interp_v_v, cellSizeInv, mf_u_inv, mf_v_inv);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        rho_w_rhs(i, j, k) = -AdvectionSrcForZMom_T(i, j, k, rho_u, rho_v, Omega, w, z_nd, detJ,
                                                    interp_w_h, interp_w_v, interp_w_wall,
                                                    cellSizeInv, mf_m, mf_u_inv, mf_v_inv,
                                                    vert_adv_type, domhi_z);
    });
}

/**
 * Wrapper function for computing the advective tendency w/ spatial order > 2.
 */
template<typename InterpType_H>
void
AdvectionSrcForMomVert_T(const amrex::Box& bxx, const amrex::Box& bxy, const amrex::Box& bxz,
                            const amrex::Array4<amrex::Real>& rho_u_rhs,
                            const amrex::Array4<amrex::Real>& rho_v_rhs,
                            const amrex::Array4<amrex::Real>& rho_w_rhs,
                            const amrex::Array4<const amrex::Real>& rho_u,
                            const amrex::Array4<const amrex::Real>& rho_v,
                            const amrex::Array4<const amrex::Real>& Omega,
                            const amrex::Array4<const amrex::Real>& u,
                            const amrex::Array4<const amrex::Real>& v,
                            const amrex::Array4<const amrex::Real>& w,
                            const amrex::Array4<const amrex::Real>& z_nd,
                            const amrex::Array4<const amrex::Real>& detJ,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                            const amrex::Array4<const amrex::Real>& mf_m,
                            const amrex::Array4<const amrex::Real>& mf_u_inv,
                            const amrex::Array4<const amrex::Real>& mf_v_inv,
                            const AdvType vert_adv_type, const int domhi_z)
{
    if (vert_adv_type == AdvType::Centered_2nd) {
        AdvectionSrcForMomWrapper_T<InterpType_H,CENTERED2,UPWINDALL>(bxx, bxy, bxz,
                                                                      rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                                      rho_u, rho_v, Omega, u, v, w, z_nd, detJ,
                                                                      cellSizeInv, mf_m, mf_u_inv, mf_v_inv,
                                                                      vert_adv_type, domhi_z);
    } else if (vert_adv_type == AdvType::Upwind_3rd) {
        AdvectionSrcForMomWrapper_T<InterpType_H,UPWIND3,UPWINDALL>(bxx, bxy, bxz,
                                                                    rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                                    rho_u, rho_v, Omega, u, v, w, z_nd, detJ,
                                                                    cellSizeInv, mf_m, mf_u_inv, mf_v_inv,
                                                                    vert_adv_type, domhi_z);
    } else if (vert_adv_type == AdvType::Centered_4th) {
        AdvectionSrcForMomWrapper_T<InterpType_H,CENTERED4,UPWINDALL>(bxx, bxy, bxz,
                                                                      rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                                      rho_u, rho_v, Omega, u, v, w, z_nd, detJ,
                                                                      cellSizeInv, mf_m, mf_u_inv, mf_v_inv,
                                                                      vert_adv_type, domhi_z);
    } else if (vert_adv_type == AdvType::Upwind_5th) {
        AdvectionSrcForMomWrapper_T<InterpType_H,UPWIND5,UPWINDALL>(bxx, bxy, bxz,
                                                                    rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                                    rho_u, rho_v, Omega, u, v, w, z_nd, detJ,
                                                                    cellSizeInv, mf_m, mf_u_inv, mf_v_inv,
                                                                    vert_adv_type, domhi_z);
    } else if (vert_adv_type == AdvType::Centered_6th) {
        AdvectionSrcForMomWrapper_T<InterpType_H,CENTERED6,UPWINDALL>(bxx, bxy, bxz,
                                                                      rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                                      rho_u, rho_v, Omega, u, v, w, z_nd, detJ,
                                                                      cellSizeInv, mf_m, mf_u_inv, mf_v_inv,
                                                                      vert_adv_type, domhi_z);
    } else {
        AMREX_ASSERT_WITH_MESSAGE(false, "Unknown advection scheme!");
    }
}
