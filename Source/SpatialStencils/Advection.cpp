#include <IndexDefines.H>
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
    Real rho_u_avg, rho_v_avg, rho_w_avg;

    Real met_h_xi, met_h_eta,  met_h_zeta;

    // ****************************************************************************************
    // X-fluxes (at cell centers)
    // ****************************************************************************************

    rho_u_avg = 0.5 * (rho_u(i+1, j, k) + rho_u(i, j, k));
    Real centFluxXXNext = rho_u_avg * Compute_h_zeta_AtCellCenter(i,j,k,cellSizeInv,z_nd) *
                          InterpolateFromCellOrFace(i+1, j, k, u, 0, rho_u_avg, Coord::x, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    rho_u_avg = 0.5 * (rho_u(i-1, j, k) + rho_u(i, j, k));
    Real centFluxXXPrev = rho_u_avg * Compute_h_zeta_AtCellCenter(i-1,j,k,cellSizeInv,z_nd) *
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

    Real edgeFluxXZNext = 0.5 *
        ( OmegaFromW(i  ,j,k+1,rho_w(i  ,j,k+1),rho_u,rho_v,z_nd,cellSizeInv)
         +OmegaFromW(i-1,j,k+1,rho_w(i-1,j,k+1),rho_u,rho_v,z_nd,cellSizeInv) ) *
        InterpolateFromCellOrFace(i, j, k+1, u, 0, rho_w_avg, Coord::z, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    Real edgeFluxXZPrev = 0.5 *
        ( OmegaFromW(i  ,j,k,rho_w(i  ,j,k),rho_u,rho_v,z_nd,cellSizeInv)
         +OmegaFromW(i-1,j,k,rho_w(i-1,j,k),rho_u,rho_v,z_nd,cellSizeInv) ) *
        InterpolateFromCellOrFace(i, j, k  , u, 0, rho_w_avg, Coord::z, spatial_order);

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
    Real rho_u_avg, rho_v_avg, rho_w_avg;

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

    rho_v_avg = 0.5 * (rho_v(i, j+1, k) + rho_v(i, j, k));
    Real centFluxYYNext = rho_v_avg * Compute_h_zeta_AtCellCenter(i,j,k,cellSizeInv,z_nd) *
                          InterpolateFromCellOrFace(i, j+1, k, v, 0, rho_v_avg, Coord::y, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    rho_v_avg = 0.5 * (rho_v(i, j-1, k) + rho_v(i, j, k));
    Real centFluxYYPrev = rho_v_avg * Compute_h_zeta_AtCellCenter(i,j-1,k,cellSizeInv,z_nd) *
                          InterpolateFromCellOrFace(i  , j, k, v, 0, rho_v_avg, Coord::y, spatial_order);


    // ****************************************************************************************
    // Z-fluxes (at edges in j-direction)
    // ****************************************************************************************

    Real edgeFluxYZNext = 0.5 *
        ( OmegaFromW(i,j  ,k,rho_w(i,j  ,k+1),rho_u,rho_v,z_nd,cellSizeInv)
         +OmegaFromW(i,j-1,k,rho_w(i,j-1,k+1),rho_u,rho_v,z_nd,cellSizeInv) ) *
        InterpolateFromCellOrFace(i, j, k+1, v, 0, rho_w_avg, Coord::z, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    Real edgeFluxYZPrev = 0.5 *
        ( OmegaFromW(i,j  ,k,rho_w(i,j  ,k),rho_u,rho_v,z_nd,cellSizeInv)
         +OmegaFromW(i,j-1,k,rho_w(i,j-1,k),rho_u,rho_v,z_nd,cellSizeInv) ) *
        InterpolateFromCellOrFace(i, j, k  , v, 0, rho_w_avg, Coord::z, spatial_order);

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
                    const int& spatial_order, const int& domhi_z)
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
           -met_h_eta * 0.5 * ( rho_v(i,j,k) + rho_v(i  ,j+1,k));
    vec += (k == domhi_z+1) ? rho_w(i,j,k) : 0.5 * ( rho_w(i,j,k) + rho_w(i  ,j  ,k+1));

    Real centFluxZZNext = vec;
    centFluxZZNext *= (k == domhi_z+1) ? rho_w(i,j,k) :
        InterpolateFromCellOrFace(i, j, k+1, w, 0, rho_w_avg, Coord::z, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // Cell-Center is staggered
    ComputeMetricAtCellCenter(i  ,j  ,k-1,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_xi_eta);

    // This is at the cell (i,j,k-1)
    vec  = -met_h_xi  * 0.5 * ( rho_u(i,j,k-1) + rho_u(i+1,j  ,k-1))
           -met_h_eta * 0.5 * ( rho_v(i,j,k-1) + rho_v(i  ,j+1,k-1));
    vec += (k == 0) ? rho_w(i,j,k) : 0.5 * (rho_w(i,j,k) + rho_w(i,j,k-1));

    Real centFluxZZPrev = vec;
    centFluxZZPrev *= (k == 0) ? rho_w(i,j,k) :
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
                    const int& spatial_order, const int& domhi_z)
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

    Real centFluxZZPrev;
    Real centFluxZZNext;

    int local_spatial_order = spatial_order;
    if (k <= 1 || k >= domhi_z) {
        local_spatial_order = std::min(local_spatial_order,2);
    } else if (k == 2 || k == domhi_z-1) {
        local_spatial_order = std::min(local_spatial_order,4);
    }

    if (k == 0) {
        centFluxZZPrev = rho_w(i,j,k) * rho_w(i,j,k);
    } else {
        centFluxZZPrev = 0.5 * (rho_w(i,j,k) + rho_w(i,j,k-1)) *
            InterpolateFromCellOrFace(i, j, k  , w, 0, rho_w_avg, Coord::z, local_spatial_order);
    }

    if (k == domhi_z+1) {
        centFluxZZNext =  rho_w(i,j,k) * rho_w(i,j,k);
    } else {
        centFluxZZNext = 0.5 * (rho_w(i,j,k) + rho_w(i,j,k+1)) *
            InterpolateFromCellOrFace(i, j, k+1, w, 0, rho_w_avg, Coord::z, local_spatial_order);
    }

    Real advectionSrc = (edgeFluxZXNext - edgeFluxZXPrev) * dxInv
                      + (edgeFluxZYNext - edgeFluxZYPrev) * dyInv
                      + (centFluxZZNext - centFluxZZPrev) * dzInv;

    return advectionSrc;
}
#endif

void
AdvectionSrcForState(const Box& bx, const int &icomp, const int &ncomp,
                     const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                     const Array4<const Real>& cell_prim,
                     const Array4<Real>& advectionSrc,
                     const Array4<Real>& xflux, const Array4<Real>& yflux, const Array4<Real>& zflux,
#ifdef ERF_USE_TERRAIN
                     const Array4<const Real>& z_nd, const Array4<const Real>& detJ,
#endif
                     const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                     const int &spatial_order, const int &use_deardorff, const int &use_QKE)
{
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {

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
    Real zflux_lo = (k == 0) ? 0.0_rt : OmegaFromW(i,j,k  ,rho_w(i,j,k  ),rho_u,rho_v,z_nd,cellSizeInv);
#else
    Real xflux_hi = rho_u(i+1,j,k);
    Real xflux_lo = rho_u(i  ,j,k);
    Real yflux_hi = rho_v(i,j+1,k);
    Real yflux_lo = rho_v(i,j  ,k);
    Real zflux_hi = rho_w(i,j,k+1);
    Real zflux_lo = rho_w(i,j,k  );
#endif

    // These are only used to construct the sign to be used in upwinding
    Real uadv_hi = rho_u(i+1,j,k);
    Real uadv_lo = rho_u(i  ,j,k);
    Real vadv_hi = rho_v(i,j+1,k);
    Real vadv_lo = rho_v(i,j  ,k);
    Real wadv_hi = rho_w(i,j,k+1);
    Real wadv_lo = rho_w(i,j,k  );

    // ****************************************************************************************
    // Now that we have the correctly weighted vector components, we can multiply by the
    //     scalar (such as theta or C) on the respective faces
    // ****************************************************************************************

#ifdef ERF_USE_TERRAIN
    Real invdetJ = 1.0 / detJ(i,j,k);
#endif

    for (int n = icomp; n < icomp+ncomp; n++)
    {
        if ((n != RhoKE_comp && n != RhoQKE_comp) ||
            (  use_deardorff && n == RhoKE_comp) ||
            (  use_QKE       && n == RhoQKE_comp) )
        {
            Real xflux_hi_n, xflux_lo_n, yflux_hi_n, yflux_lo_n, zflux_hi_n, zflux_lo_n;
            if (n != Rho_comp)
            {
                const int prim_index = n - RhoTheta_comp;

                xflux_hi_n = xflux_hi * InterpolateFromCellOrFace(i+1, j, k, cell_prim, prim_index, uadv_hi, Coord::x, spatial_order);
                xflux_lo_n = xflux_lo * InterpolateFromCellOrFace(i  , j, k, cell_prim, prim_index, uadv_lo, Coord::x, spatial_order);

                yflux_hi_n = yflux_hi * InterpolateFromCellOrFace(i, j+1, k, cell_prim, prim_index, vadv_hi, Coord::y, spatial_order);
                yflux_lo_n = yflux_lo * InterpolateFromCellOrFace(i, j  , k, cell_prim, prim_index, vadv_lo, Coord::y, spatial_order);

                zflux_hi_n = zflux_hi * InterpolateFromCellOrFace(i, j, k+1, cell_prim, prim_index, wadv_hi, Coord::z, spatial_order);
                zflux_lo_n = zflux_lo * InterpolateFromCellOrFace(i, j, k  , cell_prim, prim_index, wadv_lo, Coord::z, spatial_order);

            } else {

                xflux_hi_n = xflux_hi; xflux_lo_n = xflux_lo;
                yflux_hi_n = yflux_hi; yflux_lo_n = yflux_lo;
                zflux_hi_n = zflux_hi; zflux_lo_n = zflux_lo;
            }

            xflux(i+1,j,k,n) = xflux_hi_n;
            xflux(i  ,j,k,n) = xflux_lo_n;
            yflux(i,j+1,k,n) = yflux_hi_n;
            yflux(i,j  ,k,n) = yflux_lo_n;
            zflux(i,j,k+1,n) = zflux_hi_n;
            zflux(i,j,k  ,n) = zflux_lo_n;

            advectionSrc(i,j,k,n) = -( (xflux_hi_n - xflux_lo_n) * dxInv
                                      +(yflux_hi_n - yflux_lo_n) * dyInv
                                      +(zflux_hi_n - zflux_lo_n) * dzInv );
#ifdef ERF_USE_TERRAIN
            advectionSrc(i,j,k,n) *= invdetJ;
#endif
        } else {

            xflux(i+1,j,k,n) = 0.;
            xflux(i  ,j,k,n) = 0.;
            yflux(i,j+1,k,n) = 0.;
            yflux(i,j  ,k,n) = 0.;
            zflux(i,j,k+1,n) = 0.;
            zflux(i,j,k  ,n) = 0.;

            advectionSrc(i,j,k,n) = 0.;
        }
    } // n
    });
}
