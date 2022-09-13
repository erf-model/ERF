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
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv, int spatial_order)
{
    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    amrex::Real rho_u_avg, rho_v_avg, Omega_avg_lo, Omega_avg_hi;

    amrex::Real met_h_xi, met_h_eta,  met_h_zeta;

    // ****************************************************************************************
    // X-fluxes (at cell centers)
    // ****************************************************************************************

    met_h_zeta = detJ(i,j,k);
    rho_u_avg = 0.5 * (rho_u(i+1, j, k) + rho_u(i, j, k));
    amrex::Real centFluxXXNext = rho_u_avg * met_h_zeta *
                          InterpolateInX(i+1, j, k, u, 0, rho_u_avg, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    met_h_zeta = detJ(i-1,j,k);
    rho_u_avg = 0.5 * (rho_u(i-1, j, k) + rho_u(i, j, k));
    amrex::Real centFluxXXPrev = rho_u_avg * met_h_zeta *
                          InterpolateInX(i  , j, k, u, 0, rho_u_avg, spatial_order);

    // ****************************************************************************************
    // Y-fluxes (at edges in k-direction)
    // ****************************************************************************************

    // Metric is at edge and center Z (red pentagon)
    ComputeMetricAtEdgeCenterK(i  ,j+1,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);

    rho_v_avg = 0.5 * (rho_v(i, j+1, k) + rho_v(i-1, j+1, k));
    amrex::Real edgeFluxXYNext = rho_v_avg * met_h_zeta *
                          InterpolateInY(i, j+1, k, u, 0, rho_v_avg, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // Metric is at edge and center Z (red pentagon)
    ComputeMetricAtEdgeCenterK(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);

    rho_v_avg = 0.5 * (rho_v(i, j  , k) + rho_v(i-1, j  , k));
    amrex::Real edgeFluxXYPrev = rho_v_avg * met_h_zeta *
                          InterpolateInY(i, j  , k, u, 0, rho_v_avg, spatial_order);

    // ****************************************************************************************
    // Z-fluxes (at edges in j-direction)
    // ****************************************************************************************

    Omega_avg_hi = 0.5 * (Omega(i, j, k+1) + Omega(i-1, j, k+1));
    amrex::Real edgeFluxXZNext = Omega_avg_hi * InterpolateInZ(i,j,k+1,u,0,Omega_avg_hi,spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    Omega_avg_lo = 0.5 * (Omega(i, j, k) + Omega(i-1, j, k));
    amrex::Real edgeFluxXZPrev = Omega_avg_lo * InterpolateInZ(i,j,k  ,u,0,Omega_avg_lo,spatial_order);

    // ****************************************************************************************

    advectionSrc = (centFluxXXNext - centFluxXXPrev) * dxInv
                 + (edgeFluxXYNext - edgeFluxXYPrev) * dyInv
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
                       int spatial_order)
{
    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    amrex::Real rho_u_avg, rho_v_avg, Omega_avg_lo, Omega_avg_hi;

    amrex::Real met_h_xi, met_h_eta,  met_h_zeta;

    // ****************************************************************************************
    // x-fluxes (at edges in k-direction)
    // ****************************************************************************************

    // Metric is at edge and center Z (red pentagon)
    ComputeMetricAtEdgeCenterK(i+1,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);

    rho_u_avg = 0.5 * (rho_u(i+1, j, k) + rho_u(i+1, j-1, k));
    amrex::Real edgeFluxYXNext = rho_u_avg * met_h_zeta *
                          InterpolateInX(i+1, j, k, v, 0, rho_u_avg, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // Metric is at edge and center Z (red pentagon)
    ComputeMetricAtEdgeCenterK(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);

    rho_u_avg = 0.5 * (rho_u(i, j, k) + rho_u(i, j-1, k));
    amrex::Real edgeFluxYXPrev = rho_u_avg * met_h_zeta *
                          InterpolateInX(i  , j, k, v, 0, rho_u_avg, spatial_order);

    // ****************************************************************************************
    // y-fluxes (at cell centers)
    // ****************************************************************************************

    met_h_zeta = detJ(i,j,k),
    rho_v_avg = 0.5 * (rho_v(i, j+1, k) + rho_v(i, j, k));
    amrex::Real centFluxYYNext = rho_v_avg * met_h_zeta *
                          InterpolateInY(i, j+1, k, v, 0, rho_v_avg, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    met_h_zeta = detJ(i,j-1,k),
    rho_v_avg = 0.5 * (rho_v(i, j-1, k) + rho_v(i, j, k));
    amrex::Real centFluxYYPrev = rho_v_avg * met_h_zeta *
                          InterpolateInY(i  , j, k, v, 0, rho_v_avg, spatial_order);


    // ****************************************************************************************
    // Z-fluxes (at edges in j-direction)
    // ****************************************************************************************

    Omega_avg_hi = 0.5*(Omega(i, j, k+1) + Omega(i, j-1, k+1));
    amrex::Real edgeFluxYZNext = Omega_avg_hi *
        InterpolateInZ(i, j, k+1, v, 0, Omega_avg_hi, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    Omega_avg_lo = 0.5 * (Omega(i, j, k) + Omega(i, j-1, k));
    amrex::Real edgeFluxYZPrev = Omega_avg_lo*
                InterpolateInZ(i, j, k  , v, 0, Omega_avg_lo, spatial_order);

    // ****************************************************************************************

    advectionSrc = (edgeFluxYXNext - edgeFluxYXPrev) * dxInv
                 + (centFluxYYNext - centFluxYYPrev) * dyInv
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
                       int spatial_order, int domhi_z)
{
    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    amrex::Real rho_u_avg, rho_v_avg, Omega_avg, vec;

    amrex::Real met_h_xi, met_h_eta,  met_h_zeta;

    // ****************************************************************************************
    // x-fluxes (at edges in j-direction)
    // ****************************************************************************************

    // Metric is at edge and center Y (magenta cross)
    ComputeMetricAtEdgeCenterJ(i+1,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);

    rho_u_avg = 0.5*(rho_u(i+1,j,k) + rho_u(i+1,j,k-1));
    amrex::Real edgeFluxZXNext = rho_u_avg * met_h_zeta *
                          InterpolateInX(i+1, j, k, w, 0, rho_u_avg, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // Metric is at edge and center Y (magenta cross)
    ComputeMetricAtEdgeCenterJ(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);

    rho_u_avg = 0.5*(rho_u(i,j,k) + rho_u(i,j,k-1));
    amrex::Real edgeFluxZXPrev = rho_u_avg * met_h_zeta *
                          InterpolateInX(i  , j, k, w, 0, rho_u_avg, spatial_order);

    // ****************************************************************************************
    // y-fluxes (at edges in i-direction)
    // ****************************************************************************************

    // Metric is at edge and center I (purple hexagon)
    ComputeMetricAtEdgeCenterI(i  ,j+1,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);

    rho_v_avg = 0.5*(rho_v(i,j+1,k) + rho_v(i,j+1,k-1));
    amrex::Real edgeFluxZYNext = rho_v_avg * met_h_zeta *
                          InterpolateInY(i, j+1, k, w, 0, rho_v_avg, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // Metric is at edge and center I (purple hexagon)
    ComputeMetricAtEdgeCenterI(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);

    rho_v_avg = 0.5*(rho_v(i,j,k) + rho_v(i,j,k-1));
    amrex::Real edgeFluxZYPrev = rho_v_avg * met_h_zeta *
                          InterpolateInY(i, j  , k, w, 0, rho_v_avg, spatial_order);

    // ****************************************************************************************
    // z-fluxes (at cell centers)
    // ****************************************************************************************

    Omega_avg = (k == domhi_z+1) ? Omega(i,j,k) : 0.5 * (Omega(i,j,k) + Omega(i,j,k+1));
    amrex::Real centFluxZZNext = Omega_avg;

    centFluxZZNext *= (k == domhi_z+1) ? w(i,j,k) :
        InterpolateInZ(i, j, k+1, w, 0, Omega_avg, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    Omega_avg = (k == 0) ? Omega(i,j,k) : 0.5 * (Omega(i,j,k) + Omega(i,j,k-1));
    amrex::Real centFluxZZPrev = Omega_avg;

    centFluxZZPrev *= (k == 0) ? w(i,j,k) : InterpolateInZ(i, j, k  , w, 0, Omega_avg, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    advectionSrc = (edgeFluxZXNext - edgeFluxZXPrev) * dxInv
                 + (edgeFluxZYNext - edgeFluxZYPrev) * dyInv
                 + (centFluxZZNext - centFluxZZPrev) * dzInv;

    amrex::Real denom = (k == 0) ? detJ(i,j,k) : 0.5*(detJ(i,j,k) + detJ(i,j,k-1));
    advectionSrc /= denom;

    return advectionSrc;
}
