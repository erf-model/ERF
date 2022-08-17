#include <IndexDefines.H>
#include <SpatialStencils.H>
#include <TerrainMetrics.H>
#include <Interpolation.H>

AMREX_FORCE_INLINE
AMREX_GPU_DEVICE
amrex::Real
AdvectionSrcForXMom(int i, int j, int k,
                    const amrex::Array4<const amrex::Real>& rho_u, const amrex::Array4<const amrex::Real>& rho_v,
                    const amrex::Array4<const amrex::Real>& rho_w, const amrex::Array4<const amrex::Real>& z_t,
                    const amrex::Array4<const amrex::Real>& u,
                    const amrex::Array4<const amrex::Real>& z_nd, const amrex::Array4<const amrex::Real>& detJ,
                    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                    int spatial_order, int use_terrain)
{
    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    amrex::Real rho_u_avg, rho_v_avg, rho_w_avg;

    if (use_terrain) {

    amrex::Real met_h_xi, met_h_eta,  met_h_zeta;

    // ****************************************************************************************
    // X-fluxes (at cell centers)
    // ****************************************************************************************

    rho_u_avg = 0.5 * (rho_u(i+1, j, k) + rho_u(i, j, k));
    amrex::Real centFluxXXNext = rho_u_avg * Compute_h_zeta_AtCellCenter(i,j,k,cellSizeInv,z_nd) *
                          InterpolateFromCellOrFace(i+1, j, k, u, 0, rho_u_avg, Coord::x, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    rho_u_avg = 0.5 * (rho_u(i-1, j, k) + rho_u(i, j, k));
    amrex::Real centFluxXXPrev = rho_u_avg * Compute_h_zeta_AtCellCenter(i-1,j,k,cellSizeInv,z_nd) *
                          InterpolateFromCellOrFace(i  , j, k, u, 0, rho_u_avg, Coord::x, spatial_order);

    // ****************************************************************************************
    // Y-fluxes (at edges in k-direction)
    // ****************************************************************************************

    // Metric is at edge and center Z (red pentagon)
    ComputeMetricAtEdgeCenterK(i  ,j+1,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);

    rho_v_avg = 0.5 * (rho_v(i, j+1, k) + rho_v(i-1, j+1, k));
    amrex::Real edgeFluxXYNext = rho_v_avg * met_h_zeta *
                          InterpolateFromCellOrFace(i, j+1, k, u, 0, rho_v_avg, Coord::y, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // Metric is at edge and center Z (red pentagon)
    ComputeMetricAtEdgeCenterK(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);

    rho_v_avg = 0.5 * (rho_v(i, j  , k) + rho_v(i-1, j  , k));
    amrex::Real edgeFluxXYPrev = rho_v_avg * met_h_zeta *
                          InterpolateFromCellOrFace(i, j  , k, u, 0, rho_v_avg, Coord::y, spatial_order);

    // ****************************************************************************************
    // Z-fluxes (at edges in j-direction)
    // ****************************************************************************************

    amrex::Real z_t_ihi_zhi = (z_t) ? z_t(i  ,j,k+1) : 0.;
    amrex::Real z_t_ilo_zhi = (z_t) ? z_t(i-1,j,k+1) : 0.;

    rho_w_avg = 0.5 * (rho_w(i, j, k+1) + rho_w(i-1, j, k+1));

    amrex::Real edgeFluxXZNext = 0.5 *
        ( (OmegaFromW(i  ,j,k+1,rho_w(i  ,j,k+1),rho_u,rho_v,z_nd,cellSizeInv) - z_t_ihi_zhi)
         +(OmegaFromW(i-1,j,k+1,rho_w(i-1,j,k+1),rho_u,rho_v,z_nd,cellSizeInv) - z_t_ilo_zhi) ) *
        InterpolateFromCellOrFace(i, j, k+1, u, 0, rho_w_avg, Coord::z, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    amrex::Real z_t_ihi_zlo = (z_t) ? z_t(i  ,j,k) : 0.;
    amrex::Real z_t_ilo_zlo = (z_t) ? z_t(i-1,j,k) : 0.;

    rho_w_avg = 0.5 * (rho_w(i, j, k) + rho_w(i-1, j, k));

    amrex::Real edgeFluxXZPrev = 0.5 *
        ( (OmegaFromW(i  ,j,k,rho_w(i  ,j,k),rho_u,rho_v,z_nd,cellSizeInv) - z_t_ihi_zlo)
         +(OmegaFromW(i-1,j,k,rho_w(i-1,j,k),rho_u,rho_v,z_nd,cellSizeInv) - z_t_ilo_zlo) ) *
        InterpolateFromCellOrFace(i, j, k  , u, 0, rho_w_avg, Coord::z, spatial_order);

    // ****************************************************************************************

    advectionSrc = (centFluxXXNext - centFluxXXPrev) * dxInv
                 + (edgeFluxXYNext - edgeFluxXYPrev) * dyInv
                 + (edgeFluxXZNext - edgeFluxXZPrev) * dzInv;
    advectionSrc /= 0.5*(detJ(i,j,k) + detJ(i-1,j,k));

    } else {

        rho_u_avg = 0.5 * (rho_u(i+1, j, k) + rho_u(i, j, k));
        amrex::Real centFluxXXNext = rho_u_avg *
                              InterpolateFromCellOrFace(i+1, j, k, u, 0, rho_u_avg, Coord::x, spatial_order);

        rho_u_avg = 0.5 * (rho_u(i-1, j, k) + rho_u(i, j, k));
        amrex::Real centFluxXXPrev = rho_u_avg *
                              InterpolateFromCellOrFace(i  , j, k, u, 0, rho_u_avg, Coord::x, spatial_order);

        rho_v_avg = 0.5 * (rho_v(i, j+1, k) + rho_v(i-1, j+1, k));
        amrex::Real edgeFluxXYNext = rho_v_avg *
                              InterpolateFromCellOrFace(i, j+1, k, u, 0, rho_v_avg, Coord::y, spatial_order);

        rho_v_avg = 0.5 * (rho_v(i, j  , k) + rho_v(i-1, j  , k));
        amrex::Real edgeFluxXYPrev = rho_v_avg *
                              InterpolateFromCellOrFace(i, j  , k, u, 0, rho_v_avg, Coord::y, spatial_order);

        rho_w_avg = 0.5 * (rho_w(i, j, k+1) + rho_w(i-1, j, k+1));
        amrex::Real edgeFluxXZNext = rho_w_avg *
                              InterpolateFromCellOrFace(i, j, k+1, u, 0, rho_w_avg, Coord::z, spatial_order);

        rho_w_avg = 0.5 * (rho_w(i, j, k) + rho_w(i-1, j, k));
        amrex::Real edgeFluxXZPrev = rho_w_avg *
                              InterpolateFromCellOrFace(i, j, k  , u, 0, rho_w_avg, Coord::z, spatial_order);

        advectionSrc = (centFluxXXNext - centFluxXXPrev) * dxInv
                     + (edgeFluxXYNext - edgeFluxXYPrev) * dyInv
                     + (edgeFluxXZNext - edgeFluxXZPrev) * dzInv;
    }

    return advectionSrc;
}

AMREX_FORCE_INLINE
AMREX_GPU_DEVICE
amrex::Real
AdvectionSrcForYMom(int i, int j, int k,
                    const amrex::Array4<const amrex::Real>& rho_u, const amrex::Array4<const amrex::Real>& rho_v,
                    const amrex::Array4<const amrex::Real>& rho_w, const amrex::Array4<const amrex::Real>& z_t,
                    const amrex::Array4<const amrex::Real>& v,
                    const amrex::Array4<const amrex::Real>& z_nd, const amrex::Array4<const amrex::Real>& detJ,
                    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                    int spatial_order, int use_terrain)
{
    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    amrex::Real rho_u_avg, rho_v_avg, rho_w_avg;

    if (use_terrain) {

    amrex::Real met_h_xi, met_h_eta,  met_h_zeta;

    // ****************************************************************************************
    // x-fluxes (at edges in k-direction)
    // ****************************************************************************************

    // Metric is at edge and center Z (red pentagon)
    ComputeMetricAtEdgeCenterK(i+1,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);

    rho_u_avg = 0.5 * (rho_u(i+1, j, k) + rho_u(i+1, j-1, k));
    amrex::Real edgeFluxYXNext = rho_u_avg * met_h_zeta *
                          InterpolateFromCellOrFace(i+1, j, k, v, 0, rho_u_avg, Coord::x, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // Metric is at edge and center Z (red pentagon)
    ComputeMetricAtEdgeCenterK(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);

    rho_u_avg = 0.5 * (rho_u(i, j, k) + rho_u(i, j-1, k));
    amrex::Real edgeFluxYXPrev = rho_u_avg * met_h_zeta *
                          InterpolateFromCellOrFace(i  , j, k, v, 0, rho_u_avg, Coord::x, spatial_order);

    // ****************************************************************************************
    // y-fluxes (at cell centers)
    // ****************************************************************************************

    rho_v_avg = 0.5 * (rho_v(i, j+1, k) + rho_v(i, j, k));
    amrex::Real centFluxYYNext = rho_v_avg * Compute_h_zeta_AtCellCenter(i,j,k,cellSizeInv,z_nd) *
                          InterpolateFromCellOrFace(i, j+1, k, v, 0, rho_v_avg, Coord::y, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    rho_v_avg = 0.5 * (rho_v(i, j-1, k) + rho_v(i, j, k));
    amrex::Real centFluxYYPrev = rho_v_avg * Compute_h_zeta_AtCellCenter(i,j-1,k,cellSizeInv,z_nd) *
                          InterpolateFromCellOrFace(i  , j, k, v, 0, rho_v_avg, Coord::y, spatial_order);


    // ****************************************************************************************
    // Z-fluxes (at edges in j-direction)
    // ****************************************************************************************

    amrex::Real z_t_jhi_zhi = (z_t) ? z_t(i,j  ,k+1) : 0.;
    amrex::Real z_t_jlo_zhi = (z_t) ? z_t(i,j-1,k+1) : 0.;

    rho_w_avg = 0.5*(rho_w(i, j, k+1) + rho_w(i, j-1, k+1));

    amrex::Real edgeFluxYZNext = 0.5 *
        ( (OmegaFromW(i,j  ,k,rho_w(i,j  ,k+1),rho_u,rho_v,z_nd,cellSizeInv) - z_t_jhi_zhi)
         +(OmegaFromW(i,j-1,k,rho_w(i,j-1,k+1),rho_u,rho_v,z_nd,cellSizeInv) - z_t_jlo_zhi) ) *
        InterpolateFromCellOrFace(i, j, k+1, v, 0, rho_w_avg, Coord::z, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    amrex::Real z_t_jhi_zlo = (z_t) ? z_t(i,j  ,k) : 0.;
    amrex::Real z_t_jlo_zlo = (z_t) ? z_t(i,j-1,k) : 0.;

    rho_w_avg = 0.5*(rho_w(i, j, k) + rho_w(i, j-1, k));

    amrex::Real edgeFluxYZPrev = 0.5 *
        ( (OmegaFromW(i,j  ,k,rho_w(i,j  ,k),rho_u,rho_v,z_nd,cellSizeInv) - z_t_jhi_zlo)
         +(OmegaFromW(i,j-1,k,rho_w(i,j-1,k),rho_u,rho_v,z_nd,cellSizeInv) - z_t_jlo_zlo) ) *
        InterpolateFromCellOrFace(i, j, k  , v, 0, rho_w_avg, Coord::z, spatial_order);

    // ****************************************************************************************

    advectionSrc = (edgeFluxYXNext - edgeFluxYXPrev) * dxInv
                 + (centFluxYYNext - centFluxYYPrev) * dyInv
                 + (edgeFluxYZNext - edgeFluxYZPrev) * dzInv;
    advectionSrc /= 0.5*(detJ(i,j,k) + detJ(i,j-1,k));

    } else {

        rho_u_avg = 0.5*(rho_u(i+1, j, k) + rho_u(i+1, j-1, k));
        amrex::Real edgeFluxYXNext = rho_u_avg *
                              InterpolateFromCellOrFace(i+1, j, k, v, 0, rho_u_avg, Coord::x, spatial_order);

        rho_u_avg = 0.5*(rho_u(i  , j, k) + rho_u(i  , j-1, k));
        amrex::Real edgeFluxYXPrev = rho_u_avg *
                              InterpolateFromCellOrFace(i  , j, k, v, 0, rho_u_avg, Coord::x, spatial_order);

        rho_v_avg = 0.5*(rho_v(i, j, k) + rho_v(i, j+1, k));
        amrex::Real centFluxYYNext = rho_v_avg *
                              InterpolateFromCellOrFace(i, j+1, k, v, 0, rho_v_avg, Coord::y, spatial_order);

        rho_v_avg = 0.5*(rho_v(i, j, k) + rho_v(i, j-1, k));
        amrex::Real centFluxYYPrev = rho_v_avg *
                              InterpolateFromCellOrFace(i, j  , k, v, 0, rho_v_avg, Coord::y, spatial_order);

        rho_w_avg = 0.5*(rho_w(i, j, k+1) + rho_w(i, j-1, k+1));
        amrex::Real edgeFluxYZNext = rho_w_avg *
                              InterpolateFromCellOrFace(i, j, k+1, v, 0, rho_w_avg, Coord::z, spatial_order);

        rho_w_avg = 0.5*(rho_w(i, j, k) + rho_w(i, j-1, k));
        amrex::Real edgeFluxYZPrev = rho_w_avg *
                              InterpolateFromCellOrFace(i, j, k  , v, 0, rho_w_avg, Coord::z, spatial_order);

        advectionSrc = (edgeFluxYXNext - edgeFluxYXPrev) * dxInv
                     + (centFluxYYNext - centFluxYYPrev) * dyInv
                     + (edgeFluxYZNext - edgeFluxYZPrev) * dzInv;
    }
    return advectionSrc;
}

AMREX_FORCE_INLINE
AMREX_GPU_DEVICE
amrex::Real
AdvectionSrcForZMom(int i, int j, int k,
                    const amrex::Array4<const amrex::Real>& rho_u, const amrex::Array4<const amrex::Real>& rho_v,
                    const amrex::Array4<const amrex::Real>& rho_w, const amrex::Array4<const amrex::Real>& z_t,
                    const amrex::Array4<const amrex::Real>& w,
                    const amrex::Array4<const amrex::Real>& z_nd, const amrex::Array4<const amrex::Real>& detJ,
                    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                    int spatial_order, int use_terrain, int domhi_z)
{
    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    amrex::Real rho_u_avg, rho_v_avg, rho_w_avg, vec;

    if (use_terrain) {

    amrex::Real met_h_xi, met_h_eta,  met_h_zeta;

    // ****************************************************************************************
    // x-fluxes (at edges in j-direction)
    // ****************************************************************************************

    // Metric is at edge and center Y (magenta cross)
    ComputeMetricAtEdgeCenterJ(i+1,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);

    rho_u_avg = 0.5*(rho_u(i+1,j,k) + rho_u(i+1,j,k-1));
    amrex::Real edgeFluxZXNext = rho_u_avg * met_h_zeta *
                          InterpolateFromCellOrFace(i+1, j, k, w, 0, rho_u_avg, Coord::x, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // Metric is at edge and center Y (magenta cross)
    ComputeMetricAtEdgeCenterJ(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);

    rho_u_avg = 0.5*(rho_u(i,j,k) + rho_u(i,j,k-1));
    amrex::Real edgeFluxZXPrev = rho_u_avg * met_h_zeta *
                          InterpolateFromCellOrFace(i  , j, k, w, 0, rho_u_avg, Coord::x, spatial_order);

    // ****************************************************************************************
    // y-fluxes (at edges in i-direction)
    // ****************************************************************************************

    // Metric is at edge and center I (purple hexagon)
    ComputeMetricAtEdgeCenterI(i  ,j+1,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);

    rho_v_avg = 0.5*(rho_v(i,j+1,k) + rho_v(i,j+1,k-1));
    amrex::Real edgeFluxZYNext = rho_v_avg * met_h_zeta *
                          InterpolateFromCellOrFace(i, j+1, k, w, 0, rho_v_avg, Coord::y, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // Metric is at edge and center I (purple hexagon)
    ComputeMetricAtEdgeCenterI(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);

    rho_v_avg = 0.5*(rho_v(i,j,k) + rho_v(i,j,k-1));
    amrex::Real edgeFluxZYPrev = rho_v_avg * met_h_zeta *
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

    amrex::Real z_t_zhi = 0.;
    if (z_t)
        z_t_zhi = (k == domhi_z+1) ? z_t(i,j,k) : 0.5 * (z_t(i,j,k) + z_t(i,j,k+1));

    amrex::Real centFluxZZNext = vec - z_t_zhi;

    centFluxZZNext *= (k == domhi_z+1) ? rho_w(i,j,k) :
        InterpolateFromCellOrFace(i, j, k+1, w, 0, rho_w_avg, Coord::z, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // Cell-Center is staggered
    ComputeMetricAtCellCenter(i  ,j  ,k-1,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_xi_eta);

    // This is at the cell (i,j,k-1)
    vec  = -met_h_xi  * 0.5 * ( rho_u(i,j,k-1) + rho_u(i+1,j  ,k-1))
           -met_h_eta * 0.5 * ( rho_v(i,j,k-1) + rho_v(i  ,j+1,k-1));
    vec += (k == 0) ? rho_w(i,j,k) : 0.5 * (rho_w(i,j,k) + rho_w(i,j,k-1));

    amrex::Real z_t_zlo = 0.;
    if (z_t)
        z_t_zlo = (k == 0) ? z_t(i,j,k) : 0.5 * (z_t(i,j,k) + z_t(i,j,k-1));

    amrex::Real centFluxZZPrev = vec - z_t_zlo;

    centFluxZZPrev *= (k == 0) ? rho_w(i,j,k) :
                          InterpolateFromCellOrFace(i, j, k  , w, 0, rho_w_avg, Coord::z, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    advectionSrc = (edgeFluxZXNext - edgeFluxZXPrev) * dxInv
                 + (edgeFluxZYNext - edgeFluxZYPrev) * dyInv
                 + (centFluxZZNext - centFluxZZPrev) * dzInv;

    amrex::Real denom = (k == 0) ? detJ(i,j,k) : 0.5*(detJ(i,j,k) + detJ(i,j,k-1));
    advectionSrc /= denom;

    } else {

        rho_u_avg = 0.5*(rho_u(i+1, j, k) + rho_u(i+1, j, k-1));
        amrex::Real edgeFluxZXNext = rho_u_avg *
                              InterpolateFromCellOrFace(i+1, j, k, w, 0, rho_u_avg, Coord::x, spatial_order);

        rho_u_avg = 0.5*(rho_u(i  , j, k) + rho_u(i  , j, k-1));
        amrex::Real edgeFluxZXPrev = rho_u_avg *
                              InterpolateFromCellOrFace(i  , j, k, w, 0, rho_u_avg, Coord::x, spatial_order);

        rho_v_avg = 0.5*(rho_v(i, j+1, k) + rho_v(i, j+1, k-1));
        amrex::Real edgeFluxZYNext = rho_v_avg *
                              InterpolateFromCellOrFace(i, j+1, k, w, 0, rho_v_avg, Coord::y, spatial_order);

        rho_v_avg = 0.5*(rho_v(i, j  , k) + rho_v(i, j  , k-1));
        amrex::Real edgeFluxZYPrev = rho_v_avg *
                              InterpolateFromCellOrFace(i, j  , k, w, 0, rho_v_avg, Coord::y, spatial_order);

        amrex::Real centFluxZZPrev;
        amrex::Real centFluxZZNext;

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

        advectionSrc = (edgeFluxZXNext - edgeFluxZXPrev) * dxInv
                     + (edgeFluxZYNext - edgeFluxZYPrev) * dyInv
                     + (centFluxZZNext - centFluxZZPrev) * dzInv;
    }
    return advectionSrc;
}
