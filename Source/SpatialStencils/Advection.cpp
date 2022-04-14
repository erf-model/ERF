#include <SpatialStencils.H>

using namespace amrex;

#ifdef ERF_USE_TERRAIN
AMREX_GPU_DEVICE
Real
AdvectionContributionForXMom(const int &i, const int &j, const int &k,
                             const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                             const Array4<const Real>& u,
                             const Array4<const Real>& z_nd, const Array4<const Real>& detJ,
                             const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                             const int& spatial_order)
{
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    Real rho_u_avg, rho_v_avg, rho_w_avg, vec;

    // ****************************************************************************************
    // X-fluxes (at cell centers)
    // ****************************************************************************************

    // This is dh/dzeta at cell center (i,j,k)
    Real met_xhi = 0.25 * dzInv * // diff in k at cell(i,j)
        ( z_nd(i+1,j  ,k+1) + z_nd(i+1,j+1,k+1) + z_nd(i,j,k+1) + z_nd(i,j+1,k+1)    // avg in i and j, hi k
         -z_nd(i+1,j  ,k  ) - z_nd(i+1,j+1,k  ) - z_nd(i,j,k  ) - z_nd(i,j+1,k  ) ); // avg in i and j, lo k

    rho_u_avg = 0.5 * (rho_u(i+1, j, k) + rho_u(i, j, k));
    Real centFluxXXNext = rho_u_avg * met_xhi *
                          InterpolateFromCellOrFace(i+1, j, k, u, 0, rho_u_avg, Coord::x, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // This is dh/dzeta at cell center (i-1,j,k)
    Real met_xlo = 0.25 * dzInv * // diff in k at cell (i-1,j)
        ( z_nd(i-1,j  ,k+1) + z_nd(i-1,j+1,k+1) + z_nd(i,j,k+1) + z_nd(i,j+1,k+1)    // avg in i and j, hi k
         -z_nd(i-1,j  ,k  ) - z_nd(i-1,j+1,k  ) - z_nd(i,j,k  ) - z_nd(i,j+1,k  ) ); // avg in i and j, lo k

    rho_u_avg = 0.5 * (rho_u(i-1, j, k) + rho_u(i, j, k));
    Real centFluxXXPrev = rho_u_avg * met_xlo *
                          InterpolateFromCellOrFace(i  , j, k, u, 0, rho_u_avg, Coord::x, spatial_order);

    // ****************************************************************************************
    // Y-fluxes (at edges in k-direction)
    // ****************************************************************************************

    // This is dh/dzeta at the edge (i-1/2, j+1/2,k)
    Real met_yhi = dzInv * (z_nd(i,j+1,k+1) - z_nd(i,j+1,k)); // diff in k at node (i-1/2,j+1/2)

    rho_v_avg = 0.5 * (rho_v(i, j+1, k) + rho_v(i-1, j+1, k));
    Real edgeFluxXYNext = rho_v_avg * met_yhi *
                          InterpolateFromCellOrFace(i, j+1, k, u, 0, rho_v_avg, Coord::y, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // This is dh/dzeta at the edge (i-1/2, j-1/2,k)
    Real met_ylo = dzInv * (z_nd(i,j,k+1) - z_nd(i,j,k)); // diff in k at node (i-1/2,j-1/2)

    rho_v_avg = 0.5 * (rho_v(i, j  , k) + rho_v(i-1, j  , k));
    Real edgeFluxXYPrev = rho_v_avg * met_ylo *
                          InterpolateFromCellOrFace(i, j  , k, u, 0, rho_v_avg, Coord::y, spatial_order);

    // ****************************************************************************************
    // Z-fluxes (at edges in j-direction)
    // ****************************************************************************************

    // This is dh/dxi at the edge (i-1/2,j,k+1/2)
    Real met_zhi_xi   = 0.5 * dxInv *
        ( z_nd(i+1,j+1,k+1) + z_nd(i+1,j  ,k+1)    // hi i, avg in j, hi k
         -z_nd(i-1,j+1,k+1) - z_nd(i-1,j  ,k+1) ); // lo i, avg in j, hi k

    // This is dh/deta at the edge (i-1/2,j,k+1/2)
    Real met_zhi_eta  =  dyInv * (z_nd(i,j+1,k+1) - z_nd(i,j,k+1)); // diff in j at node (i-1/2,k+1/2)

    vec = -met_zhi_xi  * 0.5   * (u(i,j,k) + u(i,j,k+1))
          -met_zhi_eta * 0.125 * (
             rho_v(i,j,k  ) + rho_v(i-1,j,k  ) + rho_v(i,j+1,k  ) + rho_v(i-1,j+1,k  ) +
             rho_v(i,j,k+1) + rho_v(i-1,j,k+1) + rho_v(i,j+1,k+1) + rho_v(i-1,j+1,k+1) )
          +0.5 * (rho_w(i,j,k+1) + rho_w(i-1,j,k+1));
    Real edgeFluxXZNext = vec *
                          InterpolateFromCellOrFace(i, j, k+1, u, 0, rho_w_avg, Coord::z, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // This is dh/dxi at the edge (i-1/2,j,k-1/2)
    Real met_zlo_xi   = 0.5 * dxInv *
        ( z_nd(i+1,j+1,k) + z_nd(i+1,j,k)    // hi i, avg in j, lo k
         -z_nd(i-1,j+1,k) - z_nd(i-1,j,k) ); // lo i, avg in j, lo k

    // This is dh/deta at the edge (i-1/2,j,k-1/2)
    Real met_zlo_eta  =  dyInv * (z_nd(i,j+1,k) - z_nd(i,j,k)); // diff in j at node (i-1/2,k-1/2)

    vec = -met_zlo_xi  * 0.5 * (u(i,j,k) + u(i,j,k-1))
          -met_zlo_eta * 0.125 * (
             rho_v(i,j,k  ) + rho_v(i-1,j,k  ) + rho_v(i,j+1,k  ) + rho_v(i-1,j+1,k  ) +
             rho_v(i,j,k-1) + rho_v(i-1,j,k-1) + rho_v(i,j+1,k-1) + rho_v(i-1,j+1,k-1) )
          +0.5 * (rho_w(i,j,k) + rho_w(i-1,j,k));
    Real edgeFluxXZPrev = vec *
                          InterpolateFromCellOrFace(i, j, k  , u, 0, rho_w_avg, Coord::z, spatial_order);

    // ****************************************************************************************
    // ****************************************************************************************

    Real advectionContribution = (centFluxXXNext - centFluxXXPrev) * dxInv
                               + (edgeFluxXYNext - edgeFluxXYPrev) * dyInv
                               + (edgeFluxXZNext - edgeFluxXZPrev) * dzInv;
    advectionContribution /= 0.5*(detJ(i,j,k) + detJ(i-1,j,k));

    return advectionContribution;
}
#else
AMREX_GPU_DEVICE
Real
AdvectionContributionForXMom(const int &i, const int &j, const int &k,
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

    Real advectionContribution = (centFluxXXNext - centFluxXXPrev) * dxInv
                               + (edgeFluxXYNext - edgeFluxXYPrev) * dyInv
                               + (edgeFluxXZNext - edgeFluxXZPrev) * dzInv;

    return advectionContribution;
}
#endif

#ifdef ERF_USE_TERRAIN
AMREX_GPU_DEVICE
Real
AdvectionContributionForYMom(const int &i, const int &j, const int &k,
                             const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                             const Array4<const Real>& v,
                             const Array4<const Real>& z_nd, const Array4<const Real>& detJ,
                             const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                             const int& spatial_order)
{
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    Real rho_u_avg, rho_v_avg, rho_w_avg, vec;

    // ****************************************************************************************
    // x-fluxes (at edges in k-direction)
    // ****************************************************************************************

    // This is dh/dzeta at the edge (i+1/2, j-1/2,k)
    Real met_xhi = dzInv * (z_nd(i+1,j,k+1) - z_nd(i+1,j,k)); // diff in k at node (i+1/2,j-1/2)

    rho_u_avg = 0.5 * (rho_u(i+1, j, k) + rho_u(i+1, j-1, k));
    Real edgeFluxYXNext = rho_u_avg * met_xhi *
                          InterpolateFromCellOrFace(i+1, j, k, v, 0, rho_u_avg, Coord::x, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // This is dh/dzeta at the edge (i-1/2, j-1/2,k)
    Real met_xlo = dzInv * (z_nd(i,j,k+1) - z_nd(i,j,k)); // diff in k at node (i-1/2,j-1/2)

    rho_u_avg = 0.5 * (rho_u(i, j, k) + rho_u(i, j-1, k));
    Real edgeFluxYXPrev = rho_u_avg * met_xlo *
                          InterpolateFromCellOrFace(i  , j, k, v, 0, rho_u_avg, Coord::x, spatial_order);

    // ****************************************************************************************
    // y-fluxes (at cell centers)
    // ****************************************************************************************

    // This is dh/dzeta at cell center (i,j,k)
    Real met_yhi = 0.25 * dzInv * // diff in k at cell (i,j)
        ( z_nd(i+1,j  ,k+1) + z_nd(i+1,j+1,k+1) + z_nd(i,j,k+1) + z_nd(i,j+1,k+1)    // avg in i and j, hi k
         -z_nd(i+1,j  ,k  ) - z_nd(i+1,j+1,k  ) - z_nd(i,j,k  ) - z_nd(i,j+1,k  ) ); // avg in i and j, lo k
    rho_v_avg = 0.5 * (rho_v(i, j+1, k) + rho_v(i, j, k));
    Real centFluxYYNext = rho_v_avg * met_yhi *
                          InterpolateFromCellOrFace(i, j+1, k, v, 0, rho_v_avg, Coord::y, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // This is dh/dzeta at cell center (i,j-1,k)
    Real met_ylo = 0.25 * dzInv * // diff in k at cell (i,j-1)
        ( z_nd(i+1,j  ,k+1) + z_nd(i+1,j-1,k+1) + z_nd(i,j,k+1) + z_nd(i,j-1,k+1)    // avg in i and j, hi k
         -z_nd(i+1,j  ,k  ) - z_nd(i+1,j-1,k  ) - z_nd(i,j,k  ) - z_nd(i,j-1,k  ) ); // avg in i and j, lo k

    rho_v_avg = 0.5 * (rho_v(i, j-1, k) + rho_v(i, j, k));
    Real centFluxYYPrev = rho_v_avg * met_ylo *
                          InterpolateFromCellOrFace(i  , j, k, v, 0, rho_v_avg, Coord::y, spatial_order);


    // ****************************************************************************************
    // Z-fluxes (at edges in i-direction)
    // ****************************************************************************************

    // This is dh/deta at the edge (i,j-1/2,k+1/2)
    Real met_zhi_eta = 0.5 * dyInv *
        ( z_nd(i+1,j+1,k+1) + z_nd(i,j+1,k+1)    // average in i, hi j, hi k
         -z_nd(i+1,j-1,k+1) - z_nd(i,j-1,k+1) ); // average in i, lo j, hi k

    // This is dh/dxi at the edge (i,j-1/2,k+1/2)
    Real met_zhi_xi  =  dxInv * (z_nd(i+1,j,k+1) - z_nd(i,j,k+1)); // diff in i-direction

    vec = -met_zhi_eta * 0.5 * (v(i,j,k) + v(i,j,k+1))
          -met_zhi_xi  * 0.125 * (
             rho_u(i,j,k  ) + rho_u(i,j-1,k  ) + rho_u(i+1,j,k  ) + rho_u(i+1,j-1,k  ) +
             rho_u(i,j,k+1) + rho_u(i,j-1,k+1) + rho_u(i+1,j,k+1) + rho_u(i+1,j-1,k+1) )
          +0.5 * (rho_w(i,j,k+1) + rho_w(i,j-1,k+1));
    Real edgeFluxYZNext = vec *
                          InterpolateFromCellOrFace(i, j, k+1, v, 0, rho_w_avg, Coord::z, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // This is dh/deta at the edge (i,j-1/2,k-1/2)
    Real met_zlo_eta = 0.5 * dyInv *
        ( z_nd(i+1,j+1,k) + z_nd(i,j+1,k)    // average in i, hi j, lo k
         -z_nd(i+1,j-1,k) - z_nd(i,j-1,k) ); // average in i, lo j, lo k

    // This is dh/dxi at the edge (i,j-1/2,k-1/2)
    Real met_zlo_xi  =  dxInv * (z_nd(i+1,j,k) - z_nd(i,j,k)); // diff in i-direction

    vec = -met_zlo_eta * 0.5 * (v(i,j,k) + v(i,j,k-1))
          -met_zlo_xi  * 0.125 * (
             rho_u(i,j,k  ) + rho_u(i,j-1,k  ) + rho_u(i+1,j,k  ) + rho_u(i+1,j-1,k  ) +
             rho_u(i,j,k-1) + rho_u(i,j-1,k-1) + rho_u(i+1,j,k-1) + rho_u(i+1,j-1,k-1) )
          +0.5 * (rho_w(i,j,k) + rho_w(i,j-1,k));
    Real edgeFluxYZPrev = vec *
                          InterpolateFromCellOrFace(i, j, k  , v, 0, rho_w_avg, Coord::z, spatial_order);

    // ****************************************************************************************
    // ****************************************************************************************

    Real advectionContribution = (edgeFluxYXNext - edgeFluxYXPrev) * dxInv
                               + (centFluxYYNext - centFluxYYPrev) * dyInv
                               + (edgeFluxYZNext - edgeFluxYZPrev) * dzInv;
    advectionContribution /= 0.5*(detJ(i,j,k) + detJ(i,j-1,k));

    return advectionContribution;
}
#else
AMREX_GPU_DEVICE
Real
AdvectionContributionForYMom(const int &i, const int &j, const int &k,
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

    Real advectionContribution = (edgeFluxYXNext - edgeFluxYXPrev) * dxInv
                               + (centFluxYYNext - centFluxYYPrev) * dyInv
                               + (edgeFluxYZNext - edgeFluxYZPrev) * dzInv;

    return advectionContribution;
}
#endif

#ifdef ERF_USE_TERRAIN
AMREX_GPU_DEVICE
Real
AdvectionContributionForZMom(const int &i, const int &j, const int &k,
                             const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                             const Array4<const Real>& w,
                             const Array4<const Real>& z_nd, const Array4<const Real>& detJ,
                             const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                             const int& spatial_order)
{
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    Real rho_u_avg, rho_v_avg, rho_w_avg, vec;

    // ****************************************************************************************
    // x-fluxes (at edges in j-direction)
    // ****************************************************************************************

    // This is dh/dzeta at the edge (i+1/2, j, k-1/2)
    Real met_xhi = 0.25 * dzInv * ( z_nd(i+1,j,k+1) + z_nd(i+1,j+1,k+1)   // hi i, avg in j, hi k
                                   -z_nd(i+1,j,k-1) - z_nd(i+1,j+1,k-1)); // hi i, avg in j, lo k

    rho_u_avg = 0.5*(rho_u(i+1,j,k) + rho_u(i+1,j,k-1));
    Real edgeFluxZXNext = rho_u_avg * met_xhi *
                          InterpolateFromCellOrFace(i+1, j, k, w, 0, rho_u_avg, Coord::x, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // This is dh/dzeta at the edge (i-1/2, j, k-1/2)
    Real met_xlo = 0.25 * dzInv * ( z_nd(i,j,k+1) + z_nd(i,j+1,k+1)   // lo i, avg in j, hi k
                                   -z_nd(i,j,k-1) - z_nd(i,j+1,k-1)); // lo i, avg in j, lo k

    rho_u_avg = 0.5*(rho_u(i,j,k) + rho_u(i,j,k-1));
    Real edgeFluxZXPrev = rho_u_avg * met_xlo *
                          InterpolateFromCellOrFace(i  , j, k, w, 0, rho_u_avg, Coord::x, spatial_order);

    // ****************************************************************************************
    // y-fluxes (at edges in i-direction)
    // ****************************************************************************************

    // This is dh/dzeta at the edge (i, j+1/2, k-1/2)
    Real met_yhi = 0.25 * dzInv * ( z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k+1)   // avg in i, hi j, hi k
                                   -z_nd(i,j+1,k-1) - z_nd(i+1,j+1,k-1)); // avg in i, hi j, lo k

    rho_v_avg = 0.5*(rho_v(i,j+1,k) + rho_v(i,j+1,k-1));
    Real edgeFluxZYNext = rho_v_avg * met_yhi *
                          InterpolateFromCellOrFace(i, j+1, k, w, 0, rho_v_avg, Coord::y, spatial_order);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    // This is dh/dzeta at the edge (i, j-1/2, k-1/2)
    Real met_ylo = 0.25 * dzInv * ( z_nd(i,j,k+1) + z_nd(i+1,j,k+1)   // avg in i, lo j, hi k
                                   -z_nd(i,j,k-1) - z_nd(i+1,j,k-1)); // avg in i, lo j, lo k

    rho_v_avg = 0.5*(rho_v(i,j,k) + rho_v(i,j,k-1));
    Real edgeFluxZYPrev = rho_v_avg * met_ylo *
                          InterpolateFromCellOrFace(i, j  , k, w, 0, rho_v_avg, Coord::y, spatial_order);

    // ****************************************************************************************
    // z-fluxes (at cell centers)
    // ****************************************************************************************

    // This is dh/dxi at the cell (i,j,k)
    Real met_zhi_xi  =  0.25 * dxInv *
        ( z_nd(i+1,j,k+1) + z_nd(i+1,j+1,k+1) + z_nd(i+1,j,k) + z_nd(i+1,j+1,k)
         -z_nd(i  ,j,k+1) - z_nd(i  ,j+1,k+1) - z_nd(i  ,j,k) - z_nd(i  ,j+1,k) );

    // This is dh/deta at the cell (i,j,k)
    Real met_zhi_eta  =  0.25 * dyInv *
        ( z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k+1) + z_nd(i,j+1,k) + z_nd(i+1,j+1,k)
         -z_nd(i,j  ,k+1) - z_nd(i+1,j  ,k+1) - z_nd(i,j  ,k) - z_nd(i+1,j  ,k) );

    // This is at the cell (i,j,k)
    vec  = -met_zhi_xi  * 0.5 * ( rho_u(i,j,k) + rho_u(i+1,j,k))
           -met_zhi_eta * 0.5 * ( rho_v(i,j,k) + rho_v(i,j+1,k))
           +              0.5 * ( rho_w(i,j,k) + rho_w(i,j,k+1));

    Real centFluxZZNext = vec *
                          InterpolateFromCellOrFace(i, j, k+1, w, 0, rho_w_avg, Coord::z, spatial_order);

    // This is dh/dxi at the cell (i,j,k-1)
    Real met_zlo_xi  =  0.25 * dxInv *
        ( z_nd(i+1,j,k-1) + z_nd(i+1,j+1,k-1) + z_nd(i+1,j,k) + z_nd(i+1,j+1,k)
         -z_nd(i  ,j,k-1) - z_nd(i  ,j+1,k-1) - z_nd(i  ,j,k) - z_nd(i  ,j+1,k) );

    // This is dh/deta at the cell (i,j,k-1)
    Real met_zlo_eta  =  0.25 * dyInv *
        ( z_nd(i,j+1,k-1) + z_nd(i+1,j+1,k-1) + z_nd(i,j+1,k) + z_nd(i+1,j+1,k)
         -z_nd(i,j  ,k-1) - z_nd(i+1,j  ,k-1) - z_nd(i,j  ,k) - z_nd(i+1,j  ,k) );

    // This is at the cell (i,j,k-1)
    vec  = -met_zlo_xi  * 0.5 * ( rho_u(i,j,k-1) + rho_u(i+1,j  ,k-1))
           -met_zlo_eta * 0.5 * ( rho_v(i,j,k-1) + rho_v(i  ,j+1,k-1))
           +              0.5 * ( rho_w(i,j,k-1) + rho_w(i  ,j  ,k  ));

    Real centFluxZZPrev = vec *
                          InterpolateFromCellOrFace(i, j, k  , w, 0, rho_w_avg, Coord::z, spatial_order);

    Real advectionContribution = (edgeFluxZXNext - edgeFluxZXPrev) * dxInv
                               + (edgeFluxZYNext - edgeFluxZYPrev) * dyInv
                               + (centFluxZZNext - centFluxZZPrev) * dzInv;
    advectionContribution /= 0.5*(detJ(i,j,k) + detJ(i,j,k-1));

    return advectionContribution;
}
#else
AMREX_GPU_DEVICE
Real
AdvectionContributionForZMom(const int &i, const int &j, const int &k,
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

    Real advectionContribution = (edgeFluxZXNext - edgeFluxZXPrev) * dxInv
                               + (edgeFluxZYNext - edgeFluxZYPrev) * dyInv
                               + (centFluxZZNext - centFluxZZPrev) * dzInv;

    return advectionContribution;
}
#endif

#ifdef ERF_USE_TERRAIN
AMREX_GPU_DEVICE
Real
AdvectionContributionForState(const int &i, const int &j, const int &k,
                              const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                              const Array4<const Real>& cell_prim, const int &qty_index,
                              const Array4<Real>& xflux, const Array4<Real>& yflux, const Array4<Real>& zflux,
                              const Array4<const Real>& z_nd, const Array4<const Real>& detJ,
                              const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                              const int &spatial_order) {

    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    Real advectionContribution;

    // ****************************************************************************************
    // Y-faces
    // ****************************************************************************************

    // This is dh/dzeta at x-face (i+1/2,j,k)
    Real met_xhi = 0.5 * dzInv *
                         ( z_nd(i+1,j  ,k+1) + z_nd(i+1,j+1,k+1)    // hi i, hi k
                          -z_nd(i+1,j  ,k  ) - z_nd(i+1,j+1,k  ) ); // hi i, lo k

    // This is dh/dzeta at x-face (i-1/2,j,k)
    Real met_xlo = 0.5 * dzInv *
                         ( z_nd(i  ,j  ,k+1) + z_nd(i  ,j+1,k+1)    // lo i, hi k
                          -z_nd(i  ,j  ,k  ) - z_nd(i  ,j+1,k  ) ); // lo i, lo k

    xflux(i+1,j,k,qty_index) = rho_u(i+1,j,k) * met_xhi;
    xflux(i  ,j,k,qty_index) = rho_u(i  ,j,k) * met_xlo;

    // ****************************************************************************************
    // Y-faces
    // ****************************************************************************************

    // This is dh/dzeta at y-face (i,j+1/2,k)
    Real met_yhi = 0.5 * dzInv *
                         ( z_nd(i  ,j+1,k+1) + z_nd(i+1,j+1,k+1)    // hi j, hi k
                          -z_nd(i  ,j+1,k  ) - z_nd(i+1,j+1,k  ) ); // hi j, lo k

    // This is dh/dzeta at y-face (i,j-1/2,k)
    Real met_ylo = 0.5 * dzInv *
                         ( z_nd(i  ,j  ,k+1) + z_nd(i+1,j  ,k+1)    // lo j, hi k
                          -z_nd(i  ,j  ,k  ) - z_nd(i+1,j  ,k  ) ); // lo j, lo k

    yflux(i,j+1,k,qty_index) = rho_v(i,j+1,k) * met_yhi;
    yflux(i,j  ,k,qty_index) = rho_v(i,j  ,k) * met_ylo;

    // ****************************************************************************************
    // Z-faces
    // ****************************************************************************************

    // This is dh/dxi at z-face (i,j,k+1/2)
    Real met_zhi_xi   = 0.5 * dxInv *
                              ( z_nd(i+1,j+1,k+1) + z_nd(i+1,j  ,k+1)    // hi i, hi k
                               -z_nd(i  ,j+1,k+1) - z_nd(i  ,j  ,k+1) ); // lo i, hi k

    // This is dh/deta at z-face (i,j,k+1/2)
    Real met_zhi_eta  = 0.5 * dyInv *
                              ( z_nd(i+1,j+1,k+1) + z_nd(i  ,j+1,k+1)    // hi j, hi k
                               -z_nd(i+1,j  ,k+1) - z_nd(i  ,j  ,k+1) ); // lo j, hi k

    // This is dh/dxi at z-face (i,j,k-1/2)
    Real met_zlo_xi   = 0.5 * dxInv *
                              ( z_nd(i+1,j+1,k  ) + z_nd(i+1,j  ,k  )    // hi i, lo k
                               -z_nd(i  ,j+1,k  ) - z_nd(i  ,j  ,k  ) ); // lo i, lo k

    // This is dh/deta at z-face (i,j,k-1/2)
    Real met_zlo_eta  = 0.5 * dyInv *
                              ( z_nd(i+1,j+1,k  ) + z_nd(i  ,j+1,k  )    // hi j, lo k
                               -z_nd(i+1,j  ,k  ) - z_nd(i  ,j  ,k  ) ); // lo j, lo k

    Real vec_zhi_xi   = 0.25 * ( rho_u(i,j,k+1) + rho_u(i+1,j,k+1) + rho_u(i,j,k) + rho_u(i+1,j,k));
    Real vec_zhi_eta  = 0.25 * ( rho_v(i,j,k+1) + rho_v(i,j+1,k+1) + rho_v(i,j,k) + rho_v(i,j+1,k));

    Real vec_zlo_xi   = 0.25 * ( rho_u(i,j,k-1) + rho_u(i+1,j,k-1) + rho_u(i,j,k) + rho_u(i+1,j,k));
    Real vec_zlo_eta  = 0.25 * ( rho_v(i,j,k-1) + rho_v(i,j+1,k-1) + rho_v(i,j,k) + rho_v(i,j+1,k));

    zflux(i,j,k+1,qty_index) = -met_zhi_xi * vec_zhi_xi - met_zhi_eta * vec_zhi_eta + rho_w(i,j,k+1);
    zflux(i,j,k  ,qty_index) = -met_zlo_xi * vec_zlo_xi - met_zlo_eta * vec_zlo_eta + rho_w(i,j,k  );

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

        xflux(i+1,j,k,qty_index) *=
            InterpolateFromCellOrFace(i+1, j, k, cell_prim, prim_index, uadv_hi, Coord::x, spatial_order);

        xflux(i  ,j,k,qty_index) *=
            InterpolateFromCellOrFace(i  , j, k, cell_prim, prim_index, uadv_lo, Coord::x, spatial_order);

        yflux(i,j+1,k,qty_index) *=
            InterpolateFromCellOrFace(i, j+1, k, cell_prim, prim_index, vadv_hi, Coord::y, spatial_order);

        yflux(i,j  ,k,qty_index) *=
            InterpolateFromCellOrFace(i, j  , k, cell_prim, prim_index, vadv_lo, Coord::y, spatial_order);

        zflux(i,j,k+1,qty_index) *=
            InterpolateFromCellOrFace(i, j, k+1, cell_prim, prim_index, wadv_hi, Coord::z, spatial_order);

        zflux(i,j,k  ,qty_index) *=
            InterpolateFromCellOrFace(i, j, k  , cell_prim, prim_index, wadv_lo, Coord::z, spatial_order);
    }

    advectionContribution = (xflux(i+1,j,k,qty_index) - xflux(i  ,j,k,qty_index)) * dxInv
                          + (yflux(i,j+1,k,qty_index) - yflux(i,j  ,k,qty_index)) * dyInv
                          + (zflux(i,j,k+1,qty_index) - zflux(i,j,k  ,qty_index)) * dzInv;

    advectionContribution /= detJ(i,j,k);

    return advectionContribution;
}
#else
AMREX_GPU_DEVICE
Real
AdvectionContributionForState(const int &i, const int &j, const int &k,
                              const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                              const Array4<const Real>& cell_prim, const int &qty_index,
                              const Array4<Real>& xflux, const Array4<Real>& yflux, const Array4<Real>& zflux,
                              const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                              const int &spatial_order) {

    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    Real advectionContribution;

    if (qty_index == Rho_comp)
    {
        xflux(i+1,j,k,qty_index) = rho_u(i+1,j,k);
        xflux(i  ,j,k,qty_index) = rho_u(i  ,j,k);
        yflux(i,j+1,k,qty_index) = rho_v(i,j+1,k);
        yflux(i,j  ,k,qty_index) = rho_v(i,j  ,k);
        zflux(i,j,k+1,qty_index) = rho_w(i,j,k+1);
        zflux(i,j,k  ,qty_index) = rho_w(i,j,k  );

    } else {
        const int prim_index = qty_index - RhoTheta_comp;

        Real uadv_hi = rho_u(i+1,j,k);
        Real uadv_lo = rho_u(i  ,j,k);
        Real vadv_hi = rho_v(i,j+1,k);
        Real vadv_lo = rho_v(i,j  ,k);
        Real wadv_hi = rho_w(i,j,k+1);
        Real wadv_lo = rho_w(i,j,k  );
        xflux(i+1,j,k,qty_index) = rho_u(i+1,j,k) *
            InterpolateFromCellOrFace(i+1, j, k, cell_prim, prim_index, uadv_hi, Coord::x, spatial_order);

        xflux(i  ,j,k,qty_index) = rho_u(i  ,j,k) *
            InterpolateFromCellOrFace(i  , j, k, cell_prim, prim_index, uadv_lo, Coord::x, spatial_order);

        yflux(i,j+1,k,qty_index) = rho_v(i,j+1,k) *
            InterpolateFromCellOrFace(i, j+1, k, cell_prim, prim_index, vadv_hi, Coord::y, spatial_order);

        yflux(i,j  ,k,qty_index) = rho_v(i,j  ,k) *
            InterpolateFromCellOrFace(i, j  , k, cell_prim, prim_index, vadv_lo, Coord::y, spatial_order);

        zflux(i,j,k+1,qty_index) = rho_w(i,j,k+1)*
            InterpolateFromCellOrFace(i, j, k+1, cell_prim, prim_index, wadv_hi, Coord::z, spatial_order);

        zflux(i,j,k  ,qty_index) = rho_w(i,j,k  )*
            InterpolateFromCellOrFace(i, j, k  , cell_prim, prim_index, wadv_lo, Coord::z, spatial_order);
    }

    advectionContribution = (xflux(i+1,j,k,qty_index) - xflux(i  ,j,k,qty_index)) * dxInv
                          + (yflux(i,j+1,k,qty_index) - yflux(i,j  ,k,qty_index)) * dyInv
                          + (zflux(i,j,k+1,qty_index) - zflux(i,j,k  ,qty_index)) * dzInv;

    return advectionContribution;
}
#endif
