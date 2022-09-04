#include <IndexDefines.H>
#include <SpatialStencils.H>
#include <TerrainMetrics.H>
#include <Interpolation.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForXMom_N (int i, int j, int k,
                       const amrex::Array4<const amrex::Real>& rho_u, const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& rho_w, const amrex::Array4<const amrex::Real>& z_t,
                       const amrex::Array4<const amrex::Real>& u,
                       const amrex::Array4<const amrex::Real>& z_nd, const amrex::Array4<const amrex::Real>& detJ,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       int spatial_order)
{
    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    amrex::Real rho_u_avg, rho_v_avg, rho_w_avg;


    rho_u_avg = 0.5 * (rho_u(i+1, j, k) + rho_u(i, j, k));
    amrex::Real centFluxXXNext = rho_u_avg *
                          InterpolateInX(i+1, j, k, u, 0, rho_u_avg, spatial_order);

    rho_u_avg = 0.5 * (rho_u(i-1, j, k) + rho_u(i, j, k));
    amrex::Real centFluxXXPrev = rho_u_avg *
                          InterpolateInX(i  , j, k, u, 0, rho_u_avg, spatial_order);

    rho_v_avg = 0.5 * (rho_v(i, j+1, k) + rho_v(i-1, j+1, k));
    amrex::Real edgeFluxXYNext = rho_v_avg *
                          InterpolateInY(i, j+1, k, u, 0, rho_v_avg, spatial_order);

    rho_v_avg = 0.5 * (rho_v(i, j  , k) + rho_v(i-1, j  , k));
    amrex::Real edgeFluxXYPrev = rho_v_avg *
                          InterpolateInY(i, j  , k, u, 0, rho_v_avg, spatial_order);

    rho_w_avg = 0.5 * (rho_w(i, j, k+1) + rho_w(i-1, j, k+1));
    amrex::Real edgeFluxXZNext = rho_w_avg *
                          InterpolateInZ(i, j, k+1, u, 0, rho_w_avg, spatial_order);

    rho_w_avg = 0.5 * (rho_w(i, j, k) + rho_w(i-1, j, k));
    amrex::Real edgeFluxXZPrev = rho_w_avg *
                              InterpolateInZ(i, j, k  , u, 0, rho_w_avg, spatial_order);

    advectionSrc = (centFluxXXNext - centFluxXXPrev) * dxInv
                 + (edgeFluxXYNext - edgeFluxXYPrev) * dyInv
                 + (edgeFluxXZNext - edgeFluxXZPrev) * dzInv;

    return advectionSrc;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForYMom_N (int i, int j, int k,
                       const amrex::Array4<const amrex::Real>& rho_u, const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& rho_w, const amrex::Array4<const amrex::Real>& z_t,
                       const amrex::Array4<const amrex::Real>& v,
                       const amrex::Array4<const amrex::Real>& z_nd, const amrex::Array4<const amrex::Real>& detJ,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       int spatial_order)
{
    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    amrex::Real rho_u_avg, rho_v_avg, rho_w_avg;

    rho_u_avg = 0.5*(rho_u(i+1, j, k) + rho_u(i+1, j-1, k));
    amrex::Real edgeFluxYXNext = rho_u_avg *
                          InterpolateInX(i+1, j, k, v, 0, rho_u_avg, spatial_order);

    rho_u_avg = 0.5*(rho_u(i  , j, k) + rho_u(i  , j-1, k));
    amrex::Real edgeFluxYXPrev = rho_u_avg *
                          InterpolateInX(i  , j, k, v, 0, rho_u_avg, spatial_order);

    rho_v_avg = 0.5*(rho_v(i, j, k) + rho_v(i, j+1, k));
    amrex::Real centFluxYYNext = rho_v_avg *
                          InterpolateInY(i, j+1, k, v, 0, rho_v_avg, spatial_order);

    rho_v_avg = 0.5*(rho_v(i, j, k) + rho_v(i, j-1, k));
    amrex::Real centFluxYYPrev = rho_v_avg *
                          InterpolateInY(i, j  , k, v, 0, rho_v_avg, spatial_order);

    rho_w_avg = 0.5*(rho_w(i, j, k+1) + rho_w(i, j-1, k+1));
    amrex::Real edgeFluxYZNext = rho_w_avg *
                          InterpolateInZ(i, j, k+1, v, 0, rho_w_avg, spatial_order);

    rho_w_avg = 0.5*(rho_w(i, j, k) + rho_w(i, j-1, k));
    amrex::Real edgeFluxYZPrev = rho_w_avg *
                          InterpolateInZ(i, j, k  , v, 0, rho_w_avg, spatial_order);

    advectionSrc = (edgeFluxYXNext - edgeFluxYXPrev) * dxInv
                 + (centFluxYYNext - centFluxYYPrev) * dyInv
                 + (edgeFluxYZNext - edgeFluxYZPrev) * dzInv;

    return advectionSrc;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForZMom_N (int i, int j, int k,
                       const amrex::Array4<const amrex::Real>& rho_u, const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& rho_w, const amrex::Array4<const amrex::Real>& z_t,
                       const amrex::Array4<const amrex::Real>& w,
                       const amrex::Array4<const amrex::Real>& z_nd, const amrex::Array4<const amrex::Real>& detJ,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       int spatial_order, int domhi_z)
{
    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    amrex::Real rho_u_avg, rho_v_avg, rho_w_avg, vec;

    rho_u_avg = 0.5*(rho_u(i+1, j, k) + rho_u(i+1, j, k-1));
    amrex::Real edgeFluxZXNext = rho_u_avg *
                          InterpolateInX(i+1, j, k, w, 0, rho_u_avg, spatial_order);

    rho_u_avg = 0.5*(rho_u(i  , j, k) + rho_u(i  , j, k-1));
    amrex::Real edgeFluxZXPrev = rho_u_avg *
                          InterpolateInX(i  , j, k, w, 0, rho_u_avg, spatial_order);

    rho_v_avg = 0.5*(rho_v(i, j+1, k) + rho_v(i, j+1, k-1));
    amrex::Real edgeFluxZYNext = rho_v_avg *
                          InterpolateInY(i, j+1, k, w, 0, rho_v_avg, spatial_order);

    rho_v_avg = 0.5*(rho_v(i, j  , k) + rho_v(i, j  , k-1));
    amrex::Real edgeFluxZYPrev = rho_v_avg *
                          InterpolateInY(i, j  , k, w, 0, rho_v_avg, spatial_order);

    amrex::Real centFluxZZPrev;
    amrex::Real centFluxZZNext;

    int local_spatial_order = spatial_order;
    if (k <= 1 || k >= domhi_z) {
            local_spatial_order = std::min(local_spatial_order,2);
    } else if (k == 2 || k == domhi_z-1) {
        local_spatial_order = std::min(local_spatial_order,4);
    }

    if (k == 0) {
        centFluxZZPrev = rho_w(i,j,k) * w(i,j,k);
    } else {
        rho_w_avg = 0.5 * (rho_w(i,j,k) + rho_w(i,j,k-1));
        centFluxZZPrev = rho_w_avg * InterpolateInZ(i, j, k  , w, 0, rho_w_avg, local_spatial_order);
    }

    if (k == domhi_z+1) {
        centFluxZZNext =  rho_w(i,j,k) * w(i,j,k);
    } else {
        rho_w_avg = 0.5 * (rho_w(i,j,k) + rho_w(i,j,k+1));
        centFluxZZNext = rho_w_avg * InterpolateInZ(i, j, k+1, w, 0, rho_w_avg, local_spatial_order);
    }

    advectionSrc = (edgeFluxZXNext - edgeFluxZXPrev) * dxInv
                 + (edgeFluxZYNext - edgeFluxZYPrev) * dyInv
                 + (centFluxZZNext - centFluxZZPrev) * dzInv;

    return advectionSrc;
}
