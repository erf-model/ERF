#include <IndexDefines.H>
#include <TerrainMetrics.H>
#include <Interpolation.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForXMom_T (int i, int j, int k,
                       const amrex::Array4<const amrex::Real>& rho_u, const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& Omega, const amrex::Array4<const amrex::Real>& u,
                       const amrex::Array4<const amrex::Real>& z_nd,  const amrex::Array4<const amrex::Real>& detJ,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Array4<const amrex::Real>& mf_m,  const amrex::Array4<const amrex::Real>& mf_u,
                       const amrex::Array4<const amrex::Real>& mf_v,
                       int spatial_order)
{
    //used when spatial order > 2

    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    amrex::Real rho_u_avg, rho_v_avg, Omega_avg_lo, Omega_avg_hi;

    amrex::Real met_h_zeta;

    amrex::Real mf_u_inv_hi = 1. / mf_u(i+1,j  ,0); amrex::Real mf_u_inv_mid = 1. / mf_u(i  ,j  ,0); amrex::Real mf_u_inv_lo = 1. / mf_u(i-1,j  ,0);
    amrex::Real mf_v_inv_1  = 1. / mf_v(i  ,j+1,0); amrex::Real mf_v_inv_2   = 1. / mf_v(i-1,j+1,0); amrex::Real mf_v_inv_3  = 1. / mf_v(i  ,j  ,0); amrex::Real mf_v_inv_4 = 1. / mf_v(i-1,j  ,0);
    amrex::Real mf_m_inv_hi = 1. / mf_m(i  ,j  ,0); amrex::Real mf_m_inv_lo  = 1. / mf_m(i-1,j  ,0);

    // ****************************************************************************************
    // X-fluxes (at cell centers)
    // ****************************************************************************************

    met_h_zeta = Compute_h_zeta_AtCellCenter(i  ,j  ,k  ,cellSizeInv,z_nd);
    rho_u_avg = 0.5 * (rho_u(i+1, j, k) * mf_u_inv_hi + rho_u(i, j, k) * mf_u_inv_mid);
    amrex::Real centFluxXXNext = rho_u_avg * met_h_zeta *
                          InterpolateInX(i+1, j, k, u, 0, rho_u_avg, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    met_h_zeta = Compute_h_zeta_AtCellCenter(i-1,j  ,k  ,cellSizeInv,z_nd);
    rho_u_avg = 0.5 * (rho_u(i-1, j, k) * mf_u_inv_lo + rho_u(i, j, k) * mf_u_inv_mid);
    amrex::Real centFluxXXPrev = rho_u_avg * met_h_zeta *
                          InterpolateInX(i  , j, k, u, 0, rho_u_avg, spatial_order);

    // ****************************************************************************************
    // Y-fluxes (at edges in k-direction)
    // ****************************************************************************************

    // Metric is at edge and center Z (red pentagon)
    met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i  ,j+1,k  ,cellSizeInv,z_nd);
    rho_v_avg = 0.5 * (rho_v(i, j+1, k) * mf_v_inv_1 + rho_v(i-1, j+1, k) * mf_v_inv_2);
    amrex::Real edgeFluxXYNext = rho_v_avg * met_h_zeta *
                          InterpolateInY(i, j+1, k, u, 0, rho_v_avg, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // Metric is at edge and center Z (red pentagon)
    met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i  ,j  ,k  ,cellSizeInv,z_nd);
    rho_v_avg = 0.5 * (rho_v(i, j  , k) * mf_v_inv_3 + rho_v(i-1, j  , k) * mf_v_inv_4);
    amrex::Real edgeFluxXYPrev = rho_v_avg * met_h_zeta *
                          InterpolateInY(i, j  , k, u, 0, rho_v_avg, spatial_order);

    // ****************************************************************************************
    // Z-fluxes (at edges in j-direction)
    // ****************************************************************************************

    Omega_avg_hi = 0.5 * (Omega(i, j, k+1) * mf_m_inv_hi + Omega(i-1, j, k+1) * mf_m_inv_lo);
    amrex::Real edgeFluxXZNext = Omega_avg_hi * InterpolateInZ(i,j,k+1,u,0,Omega_avg_hi,spatial_order);
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    Omega_avg_lo = 0.5 * (Omega(i, j, k) * mf_m_inv_hi + Omega(i-1, j, k) * mf_m_inv_lo);
    amrex::Real edgeFluxXZPrev = Omega_avg_lo * InterpolateInZ(i,j,k  ,u,0,Omega_avg_lo,spatial_order);

    // ****************************************************************************************

    advectionSrc = (centFluxXXNext - centFluxXXPrev) * dxInv * mf_u(i,j,0)
                 + (edgeFluxXYNext - edgeFluxXYPrev) * dyInv * mf_u(i,j,0)
                 + (edgeFluxXZNext - edgeFluxXZPrev) * dzInv;
    advectionSrc /= 0.5*(detJ(i,j,k) + detJ(i-1,j,k));

    return advectionSrc;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForYMom_T (int i, int j, int k,
                       const amrex::Array4<const amrex::Real>& rho_u, const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& Omega, const amrex::Array4<const amrex::Real>& v,
                       const amrex::Array4<const amrex::Real>& z_nd, const amrex::Array4<const amrex::Real>& detJ,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Array4<const amrex::Real>& mf_m,  const amrex::Array4<const amrex::Real>& mf_u,
                       const amrex::Array4<const amrex::Real>& mf_v,
                       int spatial_order)
{
    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    amrex::Real rho_u_avg, rho_v_avg, Omega_avg_lo, Omega_avg_hi;

    amrex::Real met_h_zeta;

    amrex::Real mf_v_inv_hi = 1. / mf_v(i  ,j+1,0); amrex::Real mf_v_inv_mid = 1. / mf_v(i  ,j  ,0); amrex::Real mf_v_inv_lo = 1. / mf_v(i  ,j-1,0);
    amrex::Real mf_u_inv_1  = 1. / mf_u(i+1,j  ,0); amrex::Real mf_u_inv_2   = 1. / mf_u(i+1,j-1,0); amrex::Real mf_u_inv_3  = 1. / mf_u(i  ,j  ,0); amrex::Real mf_u_inv_4 = 1. / mf_u(i  ,j-1,0);
    amrex::Real mf_m_inv_hi = 1. / mf_m(i  ,j  ,0); amrex::Real mf_m_inv_lo  = 1. / mf_m(i  ,j-1,0);

    // ****************************************************************************************
    // x-fluxes (at edges in k-direction)
    // ****************************************************************************************

    // Metric is at edge and center Z (red pentagon)
    met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i+1,j  ,k  ,cellSizeInv,z_nd);
    rho_u_avg = 0.5 * (rho_u(i+1, j, k) * mf_u_inv_1 + rho_u(i+1, j-1, k) * mf_u_inv_2);
    amrex::Real edgeFluxYXNext = rho_u_avg * met_h_zeta *
                          InterpolateInX(i+1, j, k, v, 0, rho_u_avg, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // Metric is at edge and center Z (red pentagon)
    met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i  ,j  ,k  ,cellSizeInv,z_nd);
    rho_u_avg = 0.5 * (rho_u(i  , j, k) * mf_u_inv_3 + rho_u(i  , j-1, k) * mf_u_inv_4);
    amrex::Real edgeFluxYXPrev = rho_u_avg * met_h_zeta *
                          InterpolateInX(i  , j, k, v, 0, rho_u_avg, spatial_order);

    // ****************************************************************************************
    // y-fluxes (at cell centers)
    // ****************************************************************************************

    met_h_zeta = Compute_h_zeta_AtCellCenter(i  ,j  ,k  ,cellSizeInv,z_nd);
    rho_v_avg = 0.5 * (rho_v(i, j, k) * mf_v_inv_mid + rho_v(i, j+1, k) * mf_v_inv_hi);
    amrex::Real centFluxYYNext = rho_v_avg * met_h_zeta *
                          InterpolateInY(i, j+1, k, v, 0, rho_v_avg, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    met_h_zeta = Compute_h_zeta_AtCellCenter(i  ,j-1,k  ,cellSizeInv,z_nd);
    rho_v_avg = 0.5 * (rho_v(i, j, k) * mf_v_inv_mid + rho_v(i, j-1, k) * mf_v_inv_lo);
    amrex::Real centFluxYYPrev = rho_v_avg * met_h_zeta *
                          InterpolateInY(i  , j, k, v, 0, rho_v_avg, spatial_order);


    // ****************************************************************************************
    // Z-fluxes (at edges in j-direction)
    // ****************************************************************************************

    Omega_avg_hi = 0.5 * (Omega(i, j, k+1) * mf_m_inv_hi + Omega(i, j-1, k+1) * mf_m_inv_lo);
    amrex::Real edgeFluxYZNext = Omega_avg_hi *
                          InterpolateInZ(i, j, k+1, v, 0, Omega_avg_hi, spatial_order);
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    Omega_avg_lo = 0.5 * (Omega(i, j, k) * mf_m_inv_hi + Omega(i, j-1, k) * mf_m_inv_lo);
    amrex::Real edgeFluxYZPrev = Omega_avg_lo*
                          InterpolateInZ(i, j, k  , v, 0, Omega_avg_lo, spatial_order);

    // ****************************************************************************************

    advectionSrc = (edgeFluxYXNext - edgeFluxYXPrev) * dxInv * mf_v(i,j,0)
                 + (centFluxYYNext - centFluxYYPrev) * dyInv * mf_v(i,j,0)
                 + (edgeFluxYZNext - edgeFluxYZPrev) * dzInv;
    advectionSrc /= 0.5*(detJ(i,j,k) + detJ(i,j-1,k));

    return advectionSrc;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForZMom_T (int i, int j, int k,
                       const amrex::Array4<const amrex::Real>& rho_u, const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& Omega, const amrex::Array4<const amrex::Real>& w,
                       const amrex::Array4<const amrex::Real>& z_nd, const amrex::Array4<const amrex::Real>& detJ,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Array4<const amrex::Real>& mf_m,  const amrex::Array4<const amrex::Real>& mf_u,
                       const amrex::Array4<const amrex::Real>& mf_v,
                       int spatial_order, int domhi_z)
{
    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    amrex::Real rho_u_avg, rho_v_avg, Omega_avg;

    amrex::Real met_h_zeta;

    amrex::Real mf_u_inv_hi = 1. / mf_u(i+1,j  ,0); amrex::Real mf_u_inv_lo = 1. / mf_u(i  ,j  ,0);
    amrex::Real mf_v_inv_hi = 1. / mf_v(i  ,j+1,0); amrex::Real mf_v_inv_lo = 1. / mf_v(i  ,j  ,0);
    amrex::Real mf_m_inv    = 1. / mf_m(i  ,j  ,0);

    // ****************************************************************************************
    // x-fluxes (at edges in j-direction)
    // ****************************************************************************************

    // Metric is at edge and center Y (magenta cross)
    met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i+1,j  ,k  ,cellSizeInv,z_nd);
    rho_u_avg = 0.5 * (rho_u(i+1, j, k) + rho_u(i+1, j, k-1)) * mf_u_inv_hi;
    amrex::Real edgeFluxZXNext = rho_u_avg * met_h_zeta *
                          InterpolateInX(i+1, j, k, w, 0, rho_u_avg, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // Metric is at edge and center Y (magenta cross)
    met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i  ,j  ,k  ,cellSizeInv,z_nd);
    rho_u_avg = 0.5 * (rho_u(i  , j, k) + rho_u(i  , j, k-1)) * mf_u_inv_lo;
    amrex::Real edgeFluxZXPrev = rho_u_avg * met_h_zeta *
                          InterpolateInX(i  , j, k, w, 0, rho_u_avg, spatial_order);

    // ****************************************************************************************
    // y-fluxes (at edges in i-direction)
    // ****************************************************************************************

    // Metric is at edge and center I (purple hexagon)
    met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i  ,j+1,k  ,cellSizeInv,z_nd);
    rho_v_avg = 0.5 * (rho_v(i, j+1, k) + rho_v(i, j+1, k-1)) * mf_v_inv_hi;
    amrex::Real edgeFluxZYNext = rho_v_avg * met_h_zeta *
                          InterpolateInY(i, j+1, k, w, 0, rho_v_avg, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // Metric is at edge and center I (purple hexagon)
    met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i  ,j  ,k  ,cellSizeInv,z_nd);
    rho_v_avg = 0.5 * (rho_v(i, j  , k) + rho_v(i, j  , k-1)) * mf_v_inv_lo;
    amrex::Real edgeFluxZYPrev = rho_v_avg * met_h_zeta *
                          InterpolateInY(i, j  , k, w, 0, rho_v_avg, spatial_order);

    // ****************************************************************************************
    // z-fluxes (at cell centers)
    // ****************************************************************************************

    Omega_avg = (k == domhi_z+1) ? Omega(i,j,k) : 0.5 * (Omega(i,j,k) + Omega(i,j,k+1));
    amrex::Real centFluxZZNext = Omega_avg * mf_m_inv;

    // If k == domhi_z   and spatial_order >= 3 we would reach to k = domhi_z+2 so we set to spatial_order = min(spatial_order,2)
    // If k == domhi_z-1 and spatial_order >= 5 we would reach to k = domhi_z+2 so we set to spatial_order = min(spatial_order,4)
    int l_spatial_order_hi = std::min(std::min(spatial_order, 2*(domhi_z+1-k)), 2*(k+1));
    centFluxZZNext *= (k == domhi_z+1) ? w(i,j,k) :
        InterpolateInZ(i, j, k+1, w, 0, Omega_avg, l_spatial_order_hi);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    Omega_avg = (k == 0) ? Omega(i,j,k) : 0.5 * (Omega(i,j,k) + Omega(i,j,k-1));
    amrex::Real centFluxZZPrev = Omega_avg * mf_m_inv;

    // If k == 1 and spatial_order >= 3 we would reach to k = -1 so we set to spatial_order = min(spatial_order,2)
    // If k == 2 and spatial_order >= 5 we would reach to k = -1 so we set to spatial_order = min(spatial_order,4)
    int l_spatial_order_lo = std::min(std::min(spatial_order, 2*(domhi_z+2-k)), 2*k);
    centFluxZZPrev *= (k == 0) ? w(i,j,k) : InterpolateInZ(i, j, k  , w, 0, Omega_avg, l_spatial_order_lo);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    advectionSrc = (edgeFluxZXNext - edgeFluxZXPrev) * dxInv * mf_m(i,j,0)
                 + (edgeFluxZYNext - edgeFluxZYPrev) * dyInv * mf_m(i,j,0)
                 + (centFluxZZNext - centFluxZZPrev) * dzInv;

    amrex::Real denom = 0.5*(detJ(i,j,k) + detJ(i,j,k-1));
    advectionSrc /= denom;

    return advectionSrc;
}
