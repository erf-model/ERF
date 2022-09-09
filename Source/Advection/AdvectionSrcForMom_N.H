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

    amrex::Real xflux_hi; amrex::Real xflux_lo;
    amrex::Real yflux_hi; amrex::Real yflux_lo;
    amrex::Real zflux_hi; amrex::Real zflux_lo;

    if (spatial_order == 2) {
        xflux_hi = 0.25 * (rho_u(i, j  , k) + rho_u(i+1, j  , k)) * (u(i+1,j,k) + u(i,j,k));
        xflux_lo = 0.25 * (rho_u(i, j  , k) + rho_u(i-1, j  , k)) * (u(i-1,j,k) + u(i,j,k));

        yflux_hi = 0.25 * (rho_v(i, j+1, k) + rho_v(i-1, j+1, k)) * (u(i,j+1,k) + u(i,j,k));
        yflux_lo = 0.25 * (rho_v(i, j  , k) + rho_v(i-1, j  , k)) * (u(i,j-1,k) + u(i,j,k));

        zflux_hi = 0.25 * (rho_w(i, j, k+1) + rho_w(i-1, j, k+1)) * (u(i,j,k+1) + u(i,j,k));
        zflux_lo = 0.25 * (rho_w(i, j, k  ) + rho_w(i-1, j, k  )) * (u(i,j,k-1) + u(i,j,k));

    } else {
        rho_u_avg = 0.5 * (rho_u(i+1, j, k) + rho_u(i, j, k));
        xflux_hi = rho_u_avg * InterpolateInX(i+1, j, k, u, 0, rho_u_avg, spatial_order);

        rho_u_avg = 0.5 * (rho_u(i-1, j, k) + rho_u(i, j, k));
        xflux_lo = rho_u_avg * InterpolateInX(i  , j, k, u, 0, rho_u_avg, spatial_order);

        rho_v_avg = 0.5 * (rho_v(i, j+1, k) + rho_v(i-1, j+1, k));
        yflux_hi = rho_v_avg * InterpolateInY(i, j+1, k, u, 0, rho_v_avg, spatial_order);

        rho_v_avg = 0.5 * (rho_v(i, j  , k) + rho_v(i-1, j  , k));
        yflux_lo = rho_v_avg * InterpolateInY(i, j  , k, u, 0, rho_v_avg, spatial_order);

        rho_w_avg = 0.5 * (rho_w(i, j, k+1) + rho_w(i-1, j, k+1));
        zflux_hi = rho_w_avg * InterpolateInZ(i, j, k+1, u, 0, rho_w_avg, spatial_order);

        rho_w_avg = 0.5 * (rho_w(i, j, k) + rho_w(i-1, j, k));
        zflux_lo = rho_w_avg * InterpolateInZ(i, j, k  , u, 0, rho_w_avg, spatial_order);
    }

    advectionSrc = (xflux_hi - xflux_lo) * dxInv
                 + (yflux_hi - yflux_lo) * dyInv
                 + (zflux_hi - zflux_lo) * dzInv;

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

    amrex::Real xflux_hi; amrex::Real xflux_lo;
    amrex::Real yflux_hi; amrex::Real yflux_lo;
    amrex::Real zflux_hi; amrex::Real zflux_lo;

    if (spatial_order == 2) {

        xflux_hi = 0.25 * (rho_u(i+1, j, k) + rho_u(i+1, j-1, k)) * (v(i+1,j,k) + v(i,j,k));
        xflux_lo = 0.25 * (rho_u(i  , j, k) + rho_u(i  , j-1, k)) * (v(i-1,j,k) + v(i,j,k));

        yflux_hi = 0.25 * (rho_v(i, j, k  ) + rho_v(i  , j+1, k)) * (v(i,j+1,k) + v(i,j,k));
        yflux_lo = 0.25 * (rho_v(i, j, k  ) + rho_v(i  , j-1, k)) * (v(i,j-1,k) + v(i,j,k));

        zflux_hi = 0.25 * (rho_w(i, j, k+1) + rho_w(i, j-1, k+1)) * (v(i,j,k+1) + v(i,j,k));
        zflux_lo = 0.25 * (rho_w(i, j, k  ) + rho_w(i, j-1, k  )) * (v(i,j,k-1) + v(i,j,k));

    } else {

        rho_u_avg = 0.5*(rho_u(i+1, j, k) + rho_u(i+1, j-1, k));
        xflux_hi = rho_u_avg * InterpolateInX(i+1, j, k, v, 0, rho_u_avg, spatial_order);

        rho_u_avg = 0.5*(rho_u(i  , j, k) + rho_u(i  , j-1, k));
        xflux_lo = rho_u_avg * InterpolateInX(i  , j, k, v, 0, rho_u_avg, spatial_order);

        rho_v_avg = 0.5*(rho_v(i, j, k) + rho_v(i, j+1, k));
        yflux_hi = rho_v_avg * InterpolateInY(i, j+1, k, v, 0, rho_v_avg, spatial_order);

        rho_v_avg = 0.5*(rho_v(i, j, k) + rho_v(i, j-1, k));
        yflux_lo = rho_v_avg * InterpolateInY(i, j  , k, v, 0, rho_v_avg, spatial_order);

        rho_w_avg = 0.5*(rho_w(i, j, k+1) + rho_w(i, j-1, k+1));
        zflux_hi = rho_w_avg * InterpolateInZ(i, j, k+1, v, 0, rho_w_avg, spatial_order);

        rho_w_avg = 0.5*(rho_w(i, j, k) + rho_w(i, j-1, k));
        zflux_lo = rho_w_avg * InterpolateInZ(i, j, k  , v, 0, rho_w_avg, spatial_order);
    }

    advectionSrc = (xflux_hi - xflux_lo) * dxInv
                 + (yflux_hi - yflux_lo) * dyInv
                 + (zflux_hi - zflux_lo) * dzInv;

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

    amrex::Real xflux_hi; amrex::Real xflux_lo;
    amrex::Real yflux_hi; amrex::Real yflux_lo;
    amrex::Real zflux_hi; amrex::Real zflux_lo;

    if (spatial_order == 2) {

        xflux_hi = 0.25*(rho_u(i+1, j, k) + rho_u(i+1, j, k-1)) * (w(i+1,j,k) + w(i,j,k));
        xflux_lo = 0.25*(rho_u(i  , j, k) + rho_u(i  , j, k-1)) * (w(i-1,j,k) + w(i,j,k));

        yflux_hi = 0.25*(rho_v(i, j+1, k) + rho_v(i, j+1, k-1)) * (w(i,j+1,k) + w(i,j,k));
        yflux_lo = 0.25*(rho_v(i, j  , k) + rho_v(i, j  , k-1)) * (w(i,j-1,k) + w(i,j,k));

        zflux_lo = (k == 0) ? rho_w(i,j,k) * w(i,j,k) :
            0.25 * (rho_w(i,j,k) + rho_w(i,j,k-1)) * (w(i,j,k) + w(i,j,k-1));

        zflux_hi = (k == domhi_z+1) ? rho_w(i,j,k) * w(i,j,k) :
            0.25 * (rho_w(i,j,k) + rho_w(i,j,k+1)) * (w(i,j,k) + w(i,j,k+1));

    } else {

        rho_u_avg = 0.5*(rho_u(i+1, j, k) + rho_u(i+1, j, k-1));
        xflux_hi = rho_u_avg * InterpolateInX(i+1, j, k, w, 0, rho_u_avg, spatial_order);

        rho_u_avg = 0.5*(rho_u(i  , j, k) + rho_u(i  , j, k-1));
        xflux_lo = rho_u_avg * InterpolateInX(i  , j, k, w, 0, rho_u_avg, spatial_order);

        rho_v_avg = 0.5*(rho_v(i, j+1, k) + rho_v(i, j+1, k-1));
        yflux_hi = rho_v_avg * InterpolateInY(i, j+1, k, w, 0, rho_v_avg, spatial_order);

        rho_v_avg = 0.5*(rho_v(i, j  , k) + rho_v(i, j  , k-1));
        yflux_lo = rho_v_avg * InterpolateInY(i, j  , k, w, 0, rho_v_avg, spatial_order);

        int local_spatial_order = spatial_order;
        if (k <= 1 || k >= domhi_z) {
                local_spatial_order = std::min(local_spatial_order,2);
        } else if (k == 2 || k == domhi_z-1) {
            local_spatial_order = std::min(local_spatial_order,4);
        }

        if (k == 0) {
            zflux_lo = rho_w(i,j,k) * w(i,j,k);
        } else {
            rho_w_avg = 0.5 * (rho_w(i,j,k) + rho_w(i,j,k-1));
            zflux_lo = rho_w_avg * InterpolateInZ(i, j, k  , w, 0, rho_w_avg, local_spatial_order);
        }

        if (k == domhi_z+1) {
            zflux_hi =  rho_w(i,j,k) * w(i,j,k);
        } else {
            rho_w_avg = 0.5 * (rho_w(i,j,k) + rho_w(i,j,k+1));
            zflux_hi = rho_w_avg * InterpolateInZ(i, j, k+1, w, 0, rho_w_avg, local_spatial_order);
        }
    }

    advectionSrc = (xflux_hi - xflux_lo) * dxInv
                 + (yflux_hi - yflux_lo) * dyInv
                 + (zflux_hi - zflux_lo) * dzInv;

    return advectionSrc;
}
