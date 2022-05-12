#include <SpatialStencils.H>

#ifdef ERF_USE_TERRAIN
#include "TerrainMetrics.H"
#include "IndexDefines.H"
#endif

using namespace amrex;

#ifdef ERF_USE_TERRAIN
AMREX_GPU_DEVICE
Real
AdvectionSrcForXMom(const int &i, const int &j, const int &k,
                    const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                    const Array4<const Real>& u,
                    const Array4<const Real>& z_nd, const Array4<const Real>& detJ,
                    const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                    const int& spatial_order)
{
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    Real rho_u_avg, rho_v_avg, rho_w_avg, vec;

    Real met_h_xi, met_h_eta,  met_h_zeta;
     
    // ****************************************************************************************
    // X-fluxes (at cell centers)
    // ****************************************************************************************
    
    // Cell-Center is staggered
    ComputeMetricAtCellCenter(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);

    rho_u_avg = 0.5 * (rho_u(i+1, j, k) + rho_u(i, j, k));
    Real centFluxXXNext = rho_u_avg * met_h_zeta *
                          InterpolateFromCellOrFace(i+1, j, k, u, 0, rho_u_avg, Coord::x, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
    // Cell-Center is staggered
    ComputeMetricAtCellCenter(i-1,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);

    rho_u_avg = 0.5 * (rho_u(i-1, j, k) + rho_u(i, j, k));
    Real centFluxXXPrev = rho_u_avg * met_h_zeta *
                          InterpolateFromCellOrFace(i  , j, k, u, 0, rho_u_avg, Coord::x, spatial_order);

    // ****************************************************************************************
    // Y-fluxes (at edges in k-direction)
    // ****************************************************************************************

    // Metric is at edge and center Z (red pentagon)
    ComputeMetricAtEdgeCenterK(i  ,j+1,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);

    rho_v_avg = 0.5 * (rho_v(i, j+1, k) + rho_v(i-1, j+1, k));
    Real edgeFluxXYNext = rho_v_avg * met_h_zeta *
                          InterpolateFromCellOrFace(i, j+1, k, u, 0, rho_v_avg, Coord::y, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
    // Metric is at edge and center Z (red pentagon)
    ComputeMetricAtEdgeCenterK(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);

    rho_v_avg = 0.5 * (rho_v(i, j  , k) + rho_v(i-1, j  , k));
    Real edgeFluxXYPrev = rho_v_avg * met_h_zeta *
                          InterpolateFromCellOrFace(i, j  , k, u, 0, rho_v_avg, Coord::y, spatial_order);

    // ****************************************************************************************
    // Z-fluxes (at edges in j-direction)
    // ****************************************************************************************

    // Metric is at edge and center Y (magenta cross)
    ComputeMetricAtEdgeCenterJ(i  ,j  ,k+1,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_xi_eta);

    vec = -met_h_xi  * 0.5   * (rho_u(i,j,k) + rho_u(i,j,k+1))
          -met_h_eta * 0.125 * (
             rho_v(i,j,k  ) + rho_v(i-1,j,k  ) + rho_v(i,j+1,k  ) + rho_v(i-1,j+1,k  ) +
             rho_v(i,j,k+1) + rho_v(i-1,j,k+1) + rho_v(i,j+1,k+1) + rho_v(i-1,j+1,k+1) )
          +0.5 * (rho_w(i,j,k+1) + rho_w(i-1,j,k+1));
    Real edgeFluxXZNext = vec *
                          InterpolateFromCellOrFace(i, j, k+1, u, 0, rho_w_avg, Coord::z, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // Metric is at edge and center Y (magenta cross)
    ComputeMetricAtEdgeCenterJ(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_xi_eta);

    vec = -met_h_xi  * 0.5 * (rho_u(i,j,k) + rho_u(i,j,k-1))
          -met_h_eta * 0.125 * (
             rho_v(i,j,k  ) + rho_v(i-1,j,k  ) + rho_v(i,j+1,k  ) + rho_v(i-1,j+1,k  ) +
             rho_v(i,j,k-1) + rho_v(i-1,j,k-1) + rho_v(i,j+1,k-1) + rho_v(i-1,j+1,k-1) )
          +0.5 * (rho_w(i,j,k) + rho_w(i-1,j,k));
    Real edgeFluxXZPrev = vec *
                          InterpolateFromCellOrFace(i, j, k  , u, 0, rho_w_avg, Coord::z, spatial_order);

    // ****************************************************************************************
    // ****************************************************************************************

    Real advectionSrc = (centFluxXXNext - centFluxXXPrev) * dxInv
                      + (edgeFluxXYNext - edgeFluxXYPrev) * dyInv
                      + (edgeFluxXZNext - edgeFluxXZPrev) * dzInv;
    advectionSrc /= 0.5*(detJ(i,j,k) + detJ(i-1,j,k));


    return advectionSrc;
}
#else
AMREX_GPU_DEVICE
Real
AdvectionSrcForXMom(const int &i, const int &j, const int &k,
                    const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                    const Array4<const Real>& u,
                    const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                    const int& spatial_order)
{
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    Real rho_u_avg, rho_v_avg, rho_w_avg;

    rho_u_avg = 0.5 * (rho_u(i+1, j, k) + rho_u(i, j, k));
    Real centFluxXXNext = rho_u_avg *
                          InterpolateFromCellOrFace(i+1, j, k, u, 0, rho_u_avg, Coord::x, spatial_order);

    rho_u_avg = 0.5 * (rho_u(i-1, j, k) + rho_u(i, j, k));
    Real centFluxXXPrev = rho_u_avg *
                          InterpolateFromCellOrFace(i  , j, k, u, 0, rho_u_avg, Coord::x, spatial_order);

    rho_v_avg = 0.5 * (rho_v(i, j+1, k) + rho_v(i-1, j+1, k));
    Real edgeFluxXYNext = rho_v_avg *
                          InterpolateFromCellOrFace(i, j+1, k, u, 0, rho_v_avg, Coord::y, spatial_order);

    rho_v_avg = 0.5 * (rho_v(i, j  , k) + rho_v(i-1, j  , k));
    Real edgeFluxXYPrev = rho_v_avg *
                          InterpolateFromCellOrFace(i, j  , k, u, 0, rho_v_avg, Coord::y, spatial_order);

    rho_w_avg = 0.5 * (rho_w(i, j, k+1) + rho_w(i-1, j, k+1));
    Real edgeFluxXZNext = rho_w_avg *
                          InterpolateFromCellOrFace(i, j, k+1, u, 0, rho_w_avg, Coord::z, spatial_order);

    rho_w_avg = 0.5 * (rho_w(i, j, k) + rho_w(i-1, j, k));
    Real edgeFluxXZPrev = rho_w_avg *
                          InterpolateFromCellOrFace(i, j, k  , u, 0, rho_w_avg, Coord::z, spatial_order);

    Real advectionSrc = (centFluxXXNext - centFluxXXPrev) * dxInv
                               + (edgeFluxXYNext - edgeFluxXYPrev) * dyInv
                               + (edgeFluxXZNext - edgeFluxXZPrev) * dzInv;

    return advectionSrc;
}
#endif

#ifdef ERF_USE_TERRAIN
AMREX_GPU_DEVICE
Real
AdvectionSrcForYMom(const int &i, const int &j, const int &k,
                    const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                    const Array4<const Real>& v,
                    const Array4<const Real>& z_nd, const Array4<const Real>& detJ,
                    const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                    const int& spatial_order)
{
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    Real rho_u_avg, rho_v_avg, rho_w_avg, vec;

    Real met_h_xi, met_h_eta,  met_h_zeta;
    
    // ****************************************************************************************
    // x-fluxes (at edges in k-direction)
    // ****************************************************************************************

    // Metric is at edge and center Z (red pentagon)
    ComputeMetricAtEdgeCenterK(i+1,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);

    rho_u_avg = 0.5 * (rho_u(i+1, j, k) + rho_u(i+1, j-1, k));
    Real edgeFluxYXNext = rho_u_avg * met_h_zeta *
                          InterpolateFromCellOrFace(i+1, j, k, v, 0, rho_u_avg, Coord::x, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // Metric is at edge and center Z (red pentagon)
    ComputeMetricAtEdgeCenterK(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);
    
    rho_u_avg = 0.5 * (rho_u(i, j, k) + rho_u(i, j-1, k));
    Real edgeFluxYXPrev = rho_u_avg * met_h_zeta *
                          InterpolateFromCellOrFace(i  , j, k, v, 0, rho_u_avg, Coord::x, spatial_order);

    // ****************************************************************************************
    // y-fluxes (at cell centers)
    // ****************************************************************************************

    // Cell-Center is staggered
    ComputeMetricAtCellCenter(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);
    
    rho_v_avg = 0.5 * (rho_v(i, j+1, k) + rho_v(i, j, k));
    Real centFluxYYNext = rho_v_avg * met_h_zeta *
                          InterpolateFromCellOrFace(i, j+1, k, v, 0, rho_v_avg, Coord::y, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // Cell-Center is staggered
    ComputeMetricAtCellCenter(i  ,j-1,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);

    rho_v_avg = 0.5 * (rho_v(i, j-1, k) + rho_v(i, j, k));
    Real centFluxYYPrev = rho_v_avg * met_h_zeta *
                          InterpolateFromCellOrFace(i  , j, k, v, 0, rho_v_avg, Coord::y, spatial_order);


    // ****************************************************************************************
    // Z-fluxes (at edges in i-direction)
    // ****************************************************************************************

    // Metric is at edge and center I (purple hexagon)
    ComputeMetricAtEdgeCenterI(i  ,j  ,k+1,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_xi_eta);

    vec = -met_h_eta * 0.5 * (rho_v(i,j,k) + rho_v(i,j,k+1))
          -met_h_xi  * 0.125 * (
             rho_u(i,j,k  ) + rho_u(i,j-1,k  ) + rho_u(i+1,j,k  ) + rho_u(i+1,j-1,k  ) +
             rho_u(i,j,k+1) + rho_u(i,j-1,k+1) + rho_u(i+1,j,k+1) + rho_u(i+1,j-1,k+1) )
          +0.5 * (rho_w(i,j,k+1) + rho_w(i,j-1,k+1));
    Real edgeFluxYZNext = vec *
                          InterpolateFromCellOrFace(i, j, k+1, v, 0, rho_w_avg, Coord::z, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // Metric is at edge and center I (purple hexagon)
    ComputeMetricAtEdgeCenterI(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_xi_eta);
    
    vec = -met_h_eta * 0.5 * (rho_v(i,j,k) + rho_v(i,j,k-1))
          -met_h_xi  * 0.125 * (
             rho_u(i,j,k  ) + rho_u(i,j-1,k  ) + rho_u(i+1,j,k  ) + rho_u(i+1,j-1,k  ) +
             rho_u(i,j,k-1) + rho_u(i,j-1,k-1) + rho_u(i+1,j,k-1) + rho_u(i+1,j-1,k-1) )
          +0.5 * (rho_w(i,j,k) + rho_w(i,j-1,k));
    Real edgeFluxYZPrev = vec *
                          InterpolateFromCellOrFace(i, j, k  , v, 0, rho_w_avg, Coord::z, spatial_order);

    // ****************************************************************************************
    // ****************************************************************************************

    Real advectionSrc = (edgeFluxYXNext - edgeFluxYXPrev) * dxInv
                      + (centFluxYYNext - centFluxYYPrev) * dyInv
                      + (edgeFluxYZNext - edgeFluxYZPrev) * dzInv;
    advectionSrc /= 0.5*(detJ(i,j,k) + detJ(i,j-1,k));

    return advectionSrc;
}
#else
AMREX_GPU_DEVICE
Real
AdvectionSrcForYMom(const int &i, const int &j, const int &k,
                    const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                    const Array4<const Real>& v,
                    const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                    const int& spatial_order)
{
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    Real rho_u_avg, rho_v_avg, rho_w_avg;

    rho_u_avg = 0.5*(rho_u(i+1, j, k) + rho_u(i+1, j-1, k));
    Real edgeFluxYXNext = rho_u_avg *
                          InterpolateFromCellOrFace(i+1, j, k, v, 0, rho_u_avg, Coord::x, spatial_order);

    rho_u_avg = 0.5*(rho_u(i  , j, k) + rho_u(i  , j-1, k));
    Real edgeFluxYXPrev = rho_u_avg *
                          InterpolateFromCellOrFace(i  , j, k, v, 0, rho_u_avg, Coord::x, spatial_order);

    rho_v_avg = 0.5*(rho_v(i, j, k) + rho_v(i, j+1, k));
    Real centFluxYYNext = rho_v_avg *
                          InterpolateFromCellOrFace(i, j+1, k, v, 0, rho_v_avg, Coord::y, spatial_order);

    rho_v_avg = 0.5*(rho_v(i, j, k) + rho_v(i, j-1, k));
    Real centFluxYYPrev = rho_v_avg *
                          InterpolateFromCellOrFace(i, j  , k, v, 0, rho_v_avg, Coord::y, spatial_order);

    rho_w_avg = 0.5*(rho_w(i, j, k+1) + rho_w(i, j-1, k+1));
    Real edgeFluxYZNext = rho_w_avg *
                          InterpolateFromCellOrFace(i, j, k+1, v, 0, rho_w_avg, Coord::z, spatial_order);

    rho_w_avg = 0.5*(rho_w(i, j, k) + rho_w(i, j-1, k));
    Real edgeFluxYZPrev = rho_w_avg *
                          InterpolateFromCellOrFace(i, j, k  , v, 0, rho_w_avg, Coord::z, spatial_order);

    Real advectionSrc = (edgeFluxYXNext - edgeFluxYXPrev) * dxInv
                      + (centFluxYYNext - centFluxYYPrev) * dyInv
                      + (edgeFluxYZNext - edgeFluxYZPrev) * dzInv;

    return advectionSrc;
}
#endif

#ifdef ERF_USE_TERRAIN
AMREX_GPU_DEVICE
Real
AdvectionSrcForZMom(const int &i, const int &j, const int &k,
                    const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                    const Array4<const Real>& w,
                    const Array4<const Real>& z_nd, const Array4<const Real>& detJ,
                    const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                    const int& spatial_order)
{
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    Real rho_u_avg, rho_v_avg, rho_w_avg, vec;

    Real met_h_xi, met_h_eta,  met_h_zeta;

    // ****************************************************************************************
    // x-fluxes (at edges in j-direction)
    // ****************************************************************************************

    // Metric is at edge and center Y (magenta cross)
    ComputeMetricAtEdgeCenterJ(i+1,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);

    rho_u_avg = 0.5*(rho_u(i+1,j,k) + rho_u(i+1,j,k-1));
    Real edgeFluxZXNext = rho_u_avg * met_h_zeta *
                          InterpolateFromCellOrFace(i+1, j, k, w, 0, rho_u_avg, Coord::x, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // Metric is at edge and center Y (magenta cross)
    ComputeMetricAtEdgeCenterJ(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);

    rho_u_avg = 0.5*(rho_u(i,j,k) + rho_u(i,j,k-1));
    Real edgeFluxZXPrev = rho_u_avg * met_h_zeta *
                          InterpolateFromCellOrFace(i  , j, k, w, 0, rho_u_avg, Coord::x, spatial_order);

    // ****************************************************************************************
    // y-fluxes (at edges in i-direction)
    // ****************************************************************************************

    // Metric is at edge and center I (purple hexagon)
    ComputeMetricAtEdgeCenterI(i  ,j+1,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);

    rho_v_avg = 0.5*(rho_v(i,j+1,k) + rho_v(i,j+1,k-1));
    Real edgeFluxZYNext = rho_v_avg * met_h_zeta *
                          InterpolateFromCellOrFace(i, j+1, k, w, 0, rho_v_avg, Coord::y, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // Metric is at edge and center I (purple hexagon)
    ComputeMetricAtEdgeCenterI(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);

    rho_v_avg = 0.5*(rho_v(i,j,k) + rho_v(i,j,k-1));
    Real edgeFluxZYPrev = rho_v_avg * met_h_zeta *
                          InterpolateFromCellOrFace(i, j  , k, w, 0, rho_v_avg, Coord::y, spatial_order);

    // ****************************************************************************************
    // z-fluxes (at cell centers)
    // ****************************************************************************************

    // Cell-Center is staggered
    ComputeMetricAtCellCenter(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_xi_eta);

    // This is at the cell (i,j,k)
    vec  = -met_h_xi  * 0.5 * ( rho_u(i,j,k) + rho_u(i+1,j  ,k))
           -met_h_eta * 0.5 * ( rho_v(i,j,k) + rho_v(i  ,j+1,k))
           +            0.5 * ( rho_w(i,j,k) + rho_w(i  ,j  ,k+1));

    Real centFluxZZNext = vec *
                          InterpolateFromCellOrFace(i, j, k+1, w, 0, rho_w_avg, Coord::z, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
    // Cell-Center is staggered
    ComputeMetricAtCellCenter(i  ,j  ,k-1,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_xi_eta);

    // This is at the cell (i,j,k-1)
    vec  = -met_h_xi  * 0.5 * ( rho_u(i,j,k-1) + rho_u(i+1,j  ,k-1))
           -met_h_eta * 0.5 * ( rho_v(i,j,k-1) + rho_v(i  ,j+1,k-1))
           +            0.5 * ( rho_w(i,j,k-1) + rho_w(i  ,j  ,k  ));

    Real centFluxZZPrev = vec *
                          InterpolateFromCellOrFace(i, j, k  , w, 0, rho_w_avg, Coord::z, spatial_order);

    Real advectionSrc = (edgeFluxZXNext - edgeFluxZXPrev) * dxInv
                      + (edgeFluxZYNext - edgeFluxZYPrev) * dyInv
                      + (centFluxZZNext - centFluxZZPrev) * dzInv;
    Real denom = (k == 0) ? detJ(i,j,k) : 0.5*(detJ(i,j,k) + detJ(i,j,k-1));
    advectionSrc /= denom;

    return advectionSrc;
}
#else
AMREX_GPU_DEVICE
Real
AdvectionSrcForZMom(const int &i, const int &j, const int &k,
                    const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                    const Array4<const Real>& w,
                    const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                    const int& spatial_order)
{
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    Real rho_u_avg, rho_v_avg, rho_w_avg;

    rho_u_avg = 0.5*(rho_u(i+1, j, k) + rho_u(i+1, j, k-1));
    Real edgeFluxZXNext = rho_u_avg *
                          InterpolateFromCellOrFace(i+1, j, k, w, 0, rho_u_avg, Coord::x, spatial_order);

    rho_u_avg = 0.5*(rho_u(i  , j, k) + rho_u(i  , j, k-1));
    Real edgeFluxZXPrev = rho_u_avg *
                          InterpolateFromCellOrFace(i  , j, k, w, 0, rho_u_avg, Coord::x, spatial_order);

    rho_v_avg = 0.5*(rho_v(i, j+1, k) + rho_v(i, j+1, k-1));
    Real edgeFluxZYNext = rho_v_avg *
                          InterpolateFromCellOrFace(i, j+1, k, w, 0, rho_v_avg, Coord::y, spatial_order);

    rho_v_avg = 0.5*(rho_v(i, j  , k) + rho_v(i, j  , k-1));
    Real edgeFluxZYPrev = rho_v_avg *
                          InterpolateFromCellOrFace(i, j  , k, w, 0, rho_v_avg, Coord::y, spatial_order);

    rho_w_avg = 0.5*(rho_w(i, j  , k+1) + rho_w(i, j, k));
    Real centFluxZZNext = rho_w_avg *
                          InterpolateFromCellOrFace(i, j, k+1, w, 0, rho_w_avg, Coord::z, spatial_order);

    rho_w_avg = 0.5*(rho_w(i, j  , k-1) + rho_w(i, j, k));
    Real centFluxZZPrev = rho_w_avg *
                          InterpolateFromCellOrFace(i, j, k  , w, 0, rho_w_avg, Coord::z, spatial_order);

    Real advectionSrc = (edgeFluxZXNext - edgeFluxZXPrev) * dxInv
                      + (edgeFluxZYNext - edgeFluxZYPrev) * dyInv
                      + (centFluxZZNext - centFluxZZPrev) * dzInv;

    return advectionSrc;
}
#endif

AMREX_GPU_DEVICE
Real
AdvectionSrcForState(const int &i, const int &j, const int &k,
                     const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                     const Array4<const Real>& cell_prim, const int &qty_index,
                     const Array4<Real>& xflux, const Array4<Real>& yflux, const Array4<Real>& zflux,
#ifdef ERF_USE_TERRAIN
                     const Array4<const Real>& z_nd, const Array4<const Real>& detJ,
#endif
                     const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                     const int &spatial_order)
{
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

#ifdef ERF_USE_TERRAIN
     Real met_h_xi, met_h_eta,  met_h_zeta;
     
    // ****************************************************************************************
    // X-faces
    // ****************************************************************************************

    // Metric at U location
    ComputeMetricAtIface(i+1,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);
    Real xflux_hi = rho_u(i+1,j,k) * met_h_zeta;

    // Metric at U location
    ComputeMetricAtIface(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);
    Real xflux_lo = rho_u(i  ,j,k) * met_h_zeta;

    // ****************************************************************************************
    // Y-faces
    // ****************************************************************************************
    
    // Metric at V location
    ComputeMetricAtJface(i  ,j+1,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);
    Real yflux_hi = rho_v(i,j+1,k) * met_h_zeta;
    // Metric at V location
    ComputeMetricAtJface(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);
    Real yflux_lo = rho_v(i,j  ,k) * met_h_zeta;

    // ****************************************************************************************
    // Z-faces
    // ****************************************************************************************

    Real zflux_hi = OmegaFromW(i,j,k+1,rho_w(i,j,k+1),rho_u,rho_v,z_nd,cellSizeInv);
    Real zflux_lo = OmegaFromW(i,j,k  ,rho_w(i,j,k  ),rho_u,rho_v,z_nd,cellSizeInv);
#else
    Real xflux_hi = rho_u(i+1,j,k);
    Real xflux_lo = rho_u(i  ,j,k);
    Real yflux_hi = rho_v(i,j+1,k);
    Real yflux_lo = rho_v(i,j  ,k);
    Real zflux_hi = rho_w(i,j,k+1);
    Real zflux_lo = rho_w(i,j,k  );
#endif

    // ****************************************************************************************
    // Now that we have the correctly weighted vector components, we can multiply by the
    //     scalar (such as theta or C) on the respective faces
    // ****************************************************************************************

    if (qty_index != Rho_comp)
    {
        const int prim_index = qty_index - RhoTheta_comp;

        // These are only used to construct the sign to be used in upwinding
        Real uadv_hi = rho_u(i+1,j,k);
        Real uadv_lo = rho_u(i  ,j,k);
        Real vadv_hi = rho_v(i,j+1,k);
        Real vadv_lo = rho_v(i,j  ,k);
        Real wadv_hi = rho_w(i,j,k+1);
        Real wadv_lo = rho_w(i,j,k  );

        xflux_hi *= InterpolateFromCellOrFace(i+1, j, k, cell_prim, prim_index, uadv_hi, Coord::x, spatial_order);
        xflux_lo *= InterpolateFromCellOrFace(i  , j, k, cell_prim, prim_index, uadv_lo, Coord::x, spatial_order);

        yflux_hi *= InterpolateFromCellOrFace(i, j+1, k, cell_prim, prim_index, vadv_hi, Coord::y, spatial_order);
        yflux_lo *= InterpolateFromCellOrFace(i, j  , k, cell_prim, prim_index, vadv_lo, Coord::y, spatial_order);

        zflux_hi *= InterpolateFromCellOrFace(i, j, k+1, cell_prim, prim_index, wadv_hi, Coord::z, spatial_order);
        zflux_lo *= InterpolateFromCellOrFace(i, j, k  , cell_prim, prim_index, wadv_lo, Coord::z, spatial_order);
    }

    xflux(i+1,j,k,qty_index) = xflux_hi;
    xflux(i  ,j,k,qty_index) = xflux_lo;
    yflux(i,j+1,k,qty_index) = yflux_hi;
    yflux(i,j  ,k,qty_index) = yflux_lo;
    zflux(i,j,k+1,qty_index) = zflux_hi;
    zflux(i,j,k  ,qty_index) = zflux_lo;

    Real advectionSrc = ( (xflux_hi - xflux_lo) * dxInv
                         +(yflux_hi - yflux_lo) * dyInv
                         +(zflux_hi - zflux_lo) * dzInv );
#ifdef ERF_USE_TERRAIN
    advectionSrc /= detJ(i,j,k);
#endif

    return advectionSrc;
}
#if 0
AMREX_GPU_DEVICE
Real
AdvectionSrcForState(const int &i, const int &j, const int &k,
                     const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                     const Array4<const Real>& cell_prim, const int &qty_index,
                     const Array4<Real>& xflux, const Array4<Real>& yflux, const Array4<Real>& zflux,
                     const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                     const int &spatial_order) {

    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    Real xflux_hi = rho_u(i+1,j,k);
    Real xflux_lo = rho_u(i  ,j,k);
    Real yflux_hi = rho_v(i,j+1,k);
    Real yflux_lo = rho_v(i,j  ,k);
    Real zflux_hi = rho_w(i,j,k+1);
    Real zflux_lo = rho_w(i,j,k  );

    if (qty_index != Rho_comp)
    {
        const int prim_index = qty_index - RhoTheta_comp;

        // These are only used to construct the sign to be used in upwinding
        Real uadv_hi = rho_u(i+1,j,k);
        Real uadv_lo = rho_u(i  ,j,k);
        Real vadv_hi = rho_v(i,j+1,k);
        Real vadv_lo = rho_v(i,j  ,k);
        Real wadv_hi = rho_w(i,j,k+1);
        Real wadv_lo = rho_w(i,j,k  );

        xflux_hi *= InterpolateFromCellOrFace(i+1, j, k, cell_prim, prim_index, uadv_hi, Coord::x, spatial_order);
        xflux_lo *= InterpolateFromCellOrFace(i  , j, k, cell_prim, prim_index, uadv_lo, Coord::x, spatial_order);

        yflux_hi *= InterpolateFromCellOrFace(i, j+1, k, cell_prim, prim_index, vadv_hi, Coord::y, spatial_order);
        yflux_lo *= InterpolateFromCellOrFace(i, j  , k, cell_prim, prim_index, vadv_lo, Coord::y, spatial_order);

        zflux_hi *= InterpolateFromCellOrFace(i, j, k+1, cell_prim, prim_index, wadv_hi, Coord::z, spatial_order);
        zflux_lo *= InterpolateFromCellOrFace(i, j, k  , cell_prim, prim_index, wadv_lo, Coord::z, spatial_order);
    }

    xflux(i+1,j,k,qty_index) = xflux_hi;
    xflux(i  ,j,k,qty_index) = xflux_lo;
    yflux(i,j+1,k,qty_index) = yflux_hi;
    yflux(i,j  ,k,qty_index) = yflux_lo;
    zflux(i,j,k+1,qty_index) = zflux_hi;
    zflux(i,j,k  ,qty_index) = zflux_lo;

    Real advectionSrc = (xflux_hi - xflux_lo) * dxInv
                      + (yflux_hi - yflux_lo) * dyInv
                      + (zflux_hi - zflux_lo) * dzInv;

    return advectionSrc;
}
#endif
