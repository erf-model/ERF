#include <IndexDefines.H>
#include <TerrainMetrics.H>
#include <Interpolation.H>
#include <Interpolation_WENO.H>

/**
 * Function for computing the advective tendency for the x-component of momentum
 * with metric terms and for spatial order > 2 when not using WENO
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
 * @param[in] horiz_spatial_order spatial order to be used for lateral derivatives if not using WENO (2-6)
 * @param[in] vert_spatial_order spatial order to be used for vertical derivatives if not using WENO (2-6)
 */
template<typename InterpType>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForXMom_T (int i, int j, int k,
                       const amrex::Array4<const amrex::Real>& rho_u,
                       const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& Omega,
                       InterpType interp_u,
                       const amrex::Array4<const amrex::Real>& z_nd,
                       const amrex::Array4<const amrex::Real>& detJ,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Array4<const amrex::Real>& mf_u, const
                       amrex::Array4<const amrex::Real>& mf_v,
                       int horiz_spatial_order, int vert_spatial_order)
{
    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    amrex::Real rho_u_avg_lo, rho_u_avg_hi;
    amrex::Real rho_v_avg_lo, rho_v_avg_hi;
    amrex::Real Omega_avg_lo, Omega_avg_hi;

    amrex::Real mf_u_inv_hi = 1. / mf_u(i+1,j  ,0); amrex::Real mf_u_inv_mid = 1. / mf_u(i  ,j  ,0);
    amrex::Real mf_u_inv_lo = 1. / mf_u(i-1,j  ,0);
    amrex::Real mf_v_inv_1  = 1. / mf_v(i  ,j+1,0); amrex::Real mf_v_inv_2   = 1. / mf_v(i-1,j+1,0);
    amrex::Real mf_v_inv_3  = 1. / mf_v(i  ,j  ,0); amrex::Real mf_v_inv_4   = 1. / mf_v(i-1,j  ,0);

    amrex::Real interp_hi(0.), interp_lo(0.);

    // ****************************************************************************************
    // X-fluxes (at cell centers)
    // ****************************************************************************************

    rho_u_avg_hi = 0.5 * (rho_u(i+1, j, k) * mf_u_inv_hi + rho_u(i, j, k) * mf_u_inv_mid);
    rho_u_avg_lo = 0.5 * (rho_u(i-1, j, k) * mf_u_inv_lo + rho_u(i, j, k) * mf_u_inv_mid);

    interp_u.InterpolateInX(i,j,k,0,interp_hi,interp_lo,rho_u_avg_hi,rho_u_avg_lo);

    amrex::Real centFluxXXNext = rho_u_avg_hi * Compute_h_zeta_AtCellCenter(i  ,j  ,k  ,cellSizeInv,z_nd) * interp_hi;
    amrex::Real centFluxXXPrev = rho_u_avg_lo * Compute_h_zeta_AtCellCenter(i-1,j  ,k  ,cellSizeInv,z_nd) * interp_lo;

    // ****************************************************************************************
    // Y-fluxes (at edges in k-direction)
    // ****************************************************************************************
    rho_v_avg_hi = 0.5 * (rho_v(i, j+1, k) * mf_v_inv_1 + rho_v(i-1, j+1, k) * mf_v_inv_2);
    rho_v_avg_lo = 0.5 * (rho_v(i, j  , k) * mf_v_inv_3 + rho_v(i-1, j  , k) * mf_v_inv_4);

    interp_u.InterpolateInY(i,j,k,0,interp_hi,interp_lo,rho_v_avg_hi,rho_v_avg_lo);

    amrex::Real edgeFluxXYNext = rho_v_avg_hi * Compute_h_zeta_AtEdgeCenterK(i  ,j+1,k  ,cellSizeInv,z_nd) * interp_hi;
    amrex::Real edgeFluxXYPrev = rho_v_avg_lo * Compute_h_zeta_AtEdgeCenterK(i  ,j  ,k  ,cellSizeInv,z_nd) * interp_lo;

    // ****************************************************************************************
    // Z-fluxes (at edges in j-direction)
    // ****************************************************************************************
    Omega_avg_hi = 0.5 * (Omega(i, j, k+1) + Omega(i-1, j, k+1));
    Omega_avg_lo = 0.5 * (Omega(i, j, k  ) + Omega(i-1, j, k  ));

    interp_u.InterpolateInZ(i,j,k,0,interp_hi,interp_lo,Omega_avg_hi,Omega_avg_lo);

    amrex::Real edgeFluxXZNext = Omega_avg_hi * interp_hi;
    amrex::Real edgeFluxXZPrev = Omega_avg_lo * interp_lo;

    // ****************************************************************************************

    amrex::Real mfsq = mf_u(i,j,0) * mf_u(i,j,0);

    advectionSrc = (centFluxXXNext - centFluxXXPrev) * dxInv * mfsq
                 + (edgeFluxXYNext - edgeFluxXYPrev) * dyInv * mfsq
                 + (edgeFluxXZNext - edgeFluxXZPrev) * dzInv;
    advectionSrc /= 0.5*(detJ(i,j,k) + detJ(i-1,j,k));

    return advectionSrc;
}

/**
 * Function for computing the advective tendency for the y-component of momentum
 * with metric terms and for spatial order > 2 when not using WENO
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
 * @param[in] horiz_spatial_order spatial order to be used for lateral derivatives if not using WENO (2-6)
 * @param[in] vert_spatial_order spatial order to be used for vertical derivatives if not using WENO (2-6)
 */
template<typename InterpType>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForYMom_T (int i, int j, int k,
                       const amrex::Array4<const amrex::Real>& rho_u,
                       const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& Omega,
                       InterpType interp_v,
                       const amrex::Array4<const amrex::Real>& z_nd,
                       const amrex::Array4<const amrex::Real>& detJ,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Array4<const amrex::Real>& mf_u,
                       const amrex::Array4<const amrex::Real>& mf_v,
                       int horiz_spatial_order, int vert_spatial_order)
{
    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    amrex::Real rho_u_avg_lo, rho_u_avg_hi;
    amrex::Real rho_v_avg_lo, rho_v_avg_hi;
    amrex::Real Omega_avg_lo, Omega_avg_hi;

    amrex::Real mf_v_inv_hi = 1. / mf_v(i  ,j+1,0); amrex::Real mf_v_inv_mid = 1. / mf_v(i  ,j  ,0);
    amrex::Real mf_v_inv_lo = 1. / mf_v(i  ,j-1,0);
    amrex::Real mf_u_inv_1  = 1. / mf_u(i+1,j  ,0); amrex::Real mf_u_inv_2   = 1. / mf_u(i+1,j-1,0);
    amrex::Real mf_u_inv_3  = 1. / mf_u(i  ,j  ,0); amrex::Real mf_u_inv_4   = 1. / mf_u(i  ,j-1,0);

    amrex::Real interp_hi(0.), interp_lo(0.);

    // ****************************************************************************************
    // x-fluxes (at edges in k-direction)
    // ****************************************************************************************
    rho_u_avg_hi = 0.5 * (rho_u(i+1, j, k) * mf_u_inv_1 + rho_u(i+1, j-1, k) * mf_u_inv_2);
    rho_u_avg_lo = 0.5 * (rho_u(i  , j, k) * mf_u_inv_3 + rho_u(i  , j-1, k) * mf_u_inv_4);

    interp_v.InterpolateInX(i,j,k,0,interp_hi,interp_lo,rho_u_avg_hi,rho_u_avg_lo);

    amrex::Real edgeFluxYXNext = rho_u_avg_hi * Compute_h_zeta_AtEdgeCenterK(i+1,j  ,k  ,cellSizeInv,z_nd) * interp_hi;
    amrex::Real edgeFluxYXPrev = rho_u_avg_lo * Compute_h_zeta_AtEdgeCenterK(i  ,j  ,k  ,cellSizeInv,z_nd) * interp_lo;

    // ****************************************************************************************
    // y-fluxes (at cell centers)
    // ****************************************************************************************
    rho_v_avg_hi = 0.5 * (rho_v(i, j, k) * mf_v_inv_mid + rho_v(i, j+1, k) * mf_v_inv_hi);
    rho_v_avg_lo = 0.5 * (rho_v(i, j, k) * mf_v_inv_mid + rho_v(i, j-1, k) * mf_v_inv_lo);

    interp_v.InterpolateInY(i,j,k,0,interp_hi,interp_lo,rho_v_avg_hi,rho_v_avg_lo);

    amrex::Real centFluxYYNext = rho_v_avg_hi * Compute_h_zeta_AtCellCenter(i  ,j  ,k  ,cellSizeInv,z_nd) * interp_hi;
    amrex::Real centFluxYYPrev = rho_v_avg_lo * Compute_h_zeta_AtCellCenter(i  ,j-1,k  ,cellSizeInv,z_nd) * interp_lo;

    // ****************************************************************************************
    // Z-fluxes (at edges in j-direction)
    // ****************************************************************************************
    Omega_avg_hi = 0.5 * (Omega(i, j, k+1) + Omega(i, j-1, k+1));
    Omega_avg_lo = 0.5 * (Omega(i, j, k  ) + Omega(i, j-1, k  ));

    interp_v.InterpolateInZ(i,j,k,0,interp_hi,interp_lo,Omega_avg_hi,Omega_avg_lo);

    amrex::Real edgeFluxYZNext = Omega_avg_hi * interp_hi;
    amrex::Real edgeFluxYZPrev = Omega_avg_lo * interp_lo;

    // ****************************************************************************************

    amrex::Real mfsq = mf_v(i,j,0) * mf_v(i,j,0);

    advectionSrc = (edgeFluxYXNext - edgeFluxYXPrev) * dxInv * mfsq
                 + (centFluxYYNext - centFluxYYPrev) * dyInv * mfsq
                 + (edgeFluxYZNext - edgeFluxYZPrev) * dzInv;
    advectionSrc /= 0.5*(detJ(i,j,k) + detJ(i,j-1,k));

    return advectionSrc;
}

/**
 * Function for computing the advective tendency for the z-component of momentum
 * with metric terms and for spatial order > 2 when not using WENO
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
 * @param[in] horiz_spatial_order spatial order to be used for lateral derivatives if not using WENO (2-6)
 * @param[in] vert_spatial_order spatial order to be used for vertical derivatives if not using WENO (2-6)
 * @param[in] domhi_z maximum k value in the domain
 */
template<typename InterpType>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForZMom_T (int i, int j, int k,
                       const amrex::Array4<const amrex::Real>& rho_u,
                       const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& Omega,
                       InterpType interp_omega,
                       const amrex::Array4<const amrex::Real>& z_nd,
                       const amrex::Array4<const amrex::Real>& detJ,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Array4<const amrex::Real>& mf_m,
                       const amrex::Array4<const amrex::Real>& mf_u,
                       const amrex::Array4<const amrex::Real>& mf_v,
                       int domhi_z)
{
    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    amrex::Real rho_u_avg_lo, rho_u_avg_hi;
    amrex::Real rho_v_avg_lo, rho_v_avg_hi;
    amrex::Real Omega_avg_lo, Omega_avg_hi;

    amrex::Real mf_u_inv_hi = 1. / mf_u(i+1,j  ,0); amrex::Real mf_u_inv_lo = 1. / mf_u(i  ,j  ,0);
    amrex::Real mf_v_inv_hi = 1. / mf_v(i  ,j+1,0); amrex::Real mf_v_inv_lo = 1. / mf_v(i  ,j  ,0);

    amrex::Real interp_hi(0.), interp_lo(0.);

    // ****************************************************************************************
    // x-fluxes (at edges in j-direction)
    // ****************************************************************************************
    rho_u_avg_hi = 0.5 * (rho_u(i+1, j, k) + rho_u(i+1, j, k-1)) * mf_u_inv_hi;
    rho_u_avg_lo = 0.5 * (rho_u(i  , j, k) + rho_u(i  , j, k-1)) * mf_u_inv_lo;

    interp_omega.InterpolateInX(i,j,k,0,interp_hi,interp_lo,rho_u_avg_hi,rho_u_avg_lo);

    amrex::Real edgeFluxZXNext = rho_u_avg_hi * Compute_h_zeta_AtEdgeCenterJ(i+1,j  ,k  ,cellSizeInv,z_nd) * interp_hi;
    amrex::Real edgeFluxZXPrev = rho_u_avg_lo * Compute_h_zeta_AtEdgeCenterJ(i  ,j  ,k  ,cellSizeInv,z_nd) * interp_lo;

    // ****************************************************************************************
    // y-fluxes (at edges in i-direction)
    // ****************************************************************************************
    rho_v_avg_hi = 0.5 * (rho_v(i, j+1, k) + rho_v(i, j+1, k-1)) * mf_v_inv_hi;
    rho_v_avg_lo = 0.5 * (rho_v(i, j  , k) + rho_v(i, j  , k-1)) * mf_v_inv_lo;

    interp_omega.InterpolateInY(i,j,k,0,interp_hi,interp_lo,rho_v_avg_hi,rho_v_avg_lo);

    amrex::Real edgeFluxZYNext = rho_v_avg_hi * Compute_h_zeta_AtEdgeCenterI(i  ,j+1,k  ,cellSizeInv,z_nd) * interp_hi;
    amrex::Real edgeFluxZYPrev = rho_v_avg_lo * Compute_h_zeta_AtEdgeCenterI(i  ,j  ,k  ,cellSizeInv,z_nd) * interp_lo;

    // ****************************************************************************************
    // z-fluxes (at cell centers)
    // ****************************************************************************************
    Omega_avg_hi = (k == domhi_z+1) ? Omega(i,j,k) : 0.5 * (Omega(i,j,k) + Omega(i,j,k+1));
    Omega_avg_lo = (k == 0)         ? Omega(i,j,k) : 0.5 * (Omega(i,j,k) + Omega(i,j,k-1));

    interp_omega.InterpolateInZ(i,j,k,0,interp_hi,interp_lo,Omega_avg_hi,Omega_avg_lo);

    amrex::Real centFluxZZNext = Omega_avg_hi * interp_hi;
    amrex::Real centFluxZZPrev = Omega_avg_lo * interp_lo;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    amrex::Real mfsq = mf_m(i,j,0) * mf_m(i,j,0);

    advectionSrc = (edgeFluxZXNext - edgeFluxZXPrev) * dxInv * mfsq
                 + (edgeFluxZYNext - edgeFluxZYPrev) * dyInv * mfsq
                 + (centFluxZZNext - centFluxZZPrev) * dzInv;

    amrex::Real denom = 0.5*(detJ(i,j,k) + detJ(i,j,k-1));
    advectionSrc /= denom;

    return advectionSrc;
}
