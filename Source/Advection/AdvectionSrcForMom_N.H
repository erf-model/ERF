#include <IndexDefines.H>
#include <TerrainMetrics.H>
#include <Interpolation.H>
#include <Interpolation_WENO.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForXMom_N (int i, int j, int k,
                       const amrex::Array4<const amrex::Real>& rho_u, const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& rho_w, const amrex::Array4<const amrex::Real>& u,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Array4<const amrex::Real>& mf_u,
                       const amrex::Array4<const amrex::Real>& mf_v,
                       int horiz_spatial_order, int vert_spatial_order)
{
    //used when spatial order > 2

    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    amrex::Real rho_u_avg, rho_v_avg, rho_w_avg;

    amrex::Real xflux_hi; amrex::Real xflux_lo;
    amrex::Real yflux_hi; amrex::Real yflux_lo;
    amrex::Real zflux_hi; amrex::Real zflux_lo;

    amrex::Real mf_u_inv_hi = 1. / mf_u(i+1,j  ,0); amrex::Real mf_u_inv_mid = 1. / mf_u(i  ,j  ,0); amrex::Real mf_u_inv_lo = 1. / mf_u(i-1,j  ,0);
    amrex::Real mf_v_inv_1  = 1. / mf_v(i  ,j+1,0); amrex::Real mf_v_inv_2   = 1. / mf_v(i-1,j+1,0); amrex::Real mf_v_inv_3  = 1. / mf_v(i  ,j  ,0); amrex::Real mf_v_inv_4 = 1. / mf_v(i-1,j  ,0);

    rho_u_avg = 0.5 * (rho_u(i+1, j, k) * mf_u_inv_hi + rho_u(i, j, k) * mf_u_inv_mid);
    xflux_hi = rho_u_avg * InterpolateInX(i+1, j, k, u, 0, rho_u_avg, horiz_spatial_order);

    rho_u_avg = 0.5 * (rho_u(i-1, j, k) * mf_u_inv_lo + rho_u(i, j, k) * mf_u_inv_mid);
    xflux_lo = rho_u_avg * InterpolateInX(i  , j, k, u, 0, rho_u_avg, horiz_spatial_order);

    rho_v_avg = 0.5 * (rho_v(i, j+1, k) * mf_v_inv_1 + rho_v(i-1, j+1, k) * mf_v_inv_2);
    yflux_hi = rho_v_avg * InterpolateInY(i, j+1, k, u, 0, rho_v_avg, horiz_spatial_order);

    rho_v_avg = 0.5 * (rho_v(i, j  , k) * mf_v_inv_3 + rho_v(i-1, j  , k) * mf_v_inv_4);
    yflux_lo = rho_v_avg * InterpolateInY(i, j  , k, u, 0, rho_v_avg, horiz_spatial_order);

    rho_w_avg = 0.5 * (rho_w(i, j, k+1) + rho_w(i-1, j, k+1));
    zflux_hi = rho_w_avg * InterpolateInZ(i, j, k+1, u, 0, rho_w_avg, vert_spatial_order);

    rho_w_avg = 0.5 * (rho_w(i, j, k) + rho_w(i-1, j, k));
    zflux_lo = rho_w_avg * InterpolateInZ(i, j, k  , u, 0, rho_w_avg, vert_spatial_order);

    amrex::Real mfsq = mf_u(i,j,0) * mf_u(i,j,0);

    advectionSrc = (xflux_hi - xflux_lo) * dxInv * mfsq
                 + (yflux_hi - yflux_lo) * dyInv * mfsq
                 + (zflux_hi - zflux_lo) * dzInv;

    return advectionSrc;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForYMom_N (int i, int j, int k,
                       const amrex::Array4<const amrex::Real>& rho_u, const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& rho_w, const amrex::Array4<const amrex::Real>& v,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Array4<const amrex::Real>& mf_u,
                       const amrex::Array4<const amrex::Real>& mf_v,
                       int horiz_spatial_order, int vert_spatial_order)
{
    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    amrex::Real rho_u_avg, rho_v_avg, rho_w_avg;

    amrex::Real xflux_hi; amrex::Real xflux_lo;
    amrex::Real yflux_hi; amrex::Real yflux_lo;
    amrex::Real zflux_hi; amrex::Real zflux_lo;

    amrex::Real mf_v_inv_hi = 1. / mf_v(i  ,j+1,0); amrex::Real mf_v_inv_mid = 1. / mf_v(i  ,j  ,0); amrex::Real mf_v_inv_lo = 1. / mf_v(i  ,j-1,0);
    amrex::Real mf_u_inv_1  = 1. / mf_u(i+1,j  ,0); amrex::Real mf_u_inv_2   = 1. / mf_u(i+1,j-1,0); amrex::Real mf_u_inv_3  = 1. / mf_u(i  ,j  ,0); amrex::Real mf_u_inv_4 = 1. / mf_u(i  ,j-1,0);

    //used when spatial order > 2

    rho_u_avg = 0.5*(rho_u(i+1, j, k) * mf_u_inv_1 + rho_u(i+1, j-1, k) * mf_u_inv_2);
    xflux_hi = rho_u_avg * InterpolateInX(i+1, j, k, v, 0, rho_u_avg, horiz_spatial_order);

    rho_u_avg = 0.5*(rho_u(i  , j, k) * mf_u_inv_3 + rho_u(i  , j-1, k) * mf_u_inv_4);
    xflux_lo = rho_u_avg * InterpolateInX(i  , j, k, v, 0, rho_u_avg, horiz_spatial_order);

    rho_v_avg = 0.5*(rho_v(i, j, k) * mf_v_inv_mid + rho_v(i, j+1, k) * mf_v_inv_hi);
    yflux_hi = rho_v_avg * InterpolateInY(i, j+1, k, v, 0, rho_v_avg, horiz_spatial_order);

    rho_v_avg = 0.5*(rho_v(i, j, k) * mf_v_inv_mid + rho_v(i, j-1, k) * mf_v_inv_lo);
    yflux_lo = rho_v_avg * InterpolateInY(i, j  , k, v, 0, rho_v_avg, horiz_spatial_order);

    rho_w_avg = 0.5*(rho_w(i, j, k+1) + rho_w(i, j-1, k+1));
    zflux_hi = rho_w_avg * InterpolateInZ(i, j, k+1, v, 0, rho_w_avg, vert_spatial_order);

    rho_w_avg = 0.5*(rho_w(i, j, k) + rho_w(i, j-1, k));
    zflux_lo = rho_w_avg * InterpolateInZ(i, j, k  , v, 0, rho_w_avg, vert_spatial_order);

    amrex::Real mfsq = mf_v(i,j,0) * mf_v(i,j,0);

    advectionSrc = (xflux_hi - xflux_lo) * dxInv * mfsq
                 + (yflux_hi - yflux_lo) * dyInv * mfsq
                 + (zflux_hi - zflux_lo) * dzInv;

   return advectionSrc;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForZMom_N (int i, int j, int k,
                       const amrex::Array4<const amrex::Real>& rho_u, const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& rho_w, const amrex::Array4<const amrex::Real>& w,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Array4<const amrex::Real>& mf_m,
                       const amrex::Array4<const amrex::Real>& mf_u,
                       const amrex::Array4<const amrex::Real>& mf_v,
                       int horiz_spatial_order, int vert_spatial_order, int domhi_z)
{

    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    amrex::Real rho_u_avg, rho_v_avg, rho_w_avg;

    amrex::Real xflux_hi; amrex::Real xflux_lo;
    amrex::Real yflux_hi; amrex::Real yflux_lo;
    amrex::Real zflux_hi; amrex::Real zflux_lo;

    amrex::Real mf_u_inv_hi = 1. / mf_u(i+1,j  ,0); amrex::Real mf_u_inv_lo = 1. / mf_u(i  ,j  ,0);
    amrex::Real mf_v_inv_hi = 1. / mf_v(i  ,j+1,0); amrex::Real mf_v_inv_lo = 1. / mf_v(i  ,j  ,0);

    //used when spatial order > 2

    rho_u_avg = 0.5*(rho_u(i+1, j, k) + rho_u(i+1, j, k-1)) * mf_u_inv_hi;
    xflux_hi = rho_u_avg * InterpolateInX(i+1, j, k, w, 0, rho_u_avg, horiz_spatial_order);

    rho_u_avg = 0.5*(rho_u(i  , j, k) + rho_u(i  , j, k-1)) * mf_u_inv_lo;
    xflux_lo = rho_u_avg * InterpolateInX(i  , j, k, w, 0, rho_u_avg, horiz_spatial_order);

    rho_v_avg = 0.5*(rho_v(i, j+1, k) + rho_v(i, j+1, k-1)) * mf_v_inv_hi;
    yflux_hi = rho_v_avg * InterpolateInY(i, j+1, k, w, 0, rho_v_avg, horiz_spatial_order);

    rho_v_avg = 0.5*(rho_v(i, j  , k) + rho_v(i, j  , k-1)) * mf_v_inv_lo;
    yflux_lo = rho_v_avg * InterpolateInY(i, j  , k, w, 0, rho_v_avg, horiz_spatial_order);

    // If k == 1 and spatial_order >= 3 we would reach to k = -1 so we set to spatial_order = min(spatial_order,2)
    // If k == 2 and spatial_order >= 5 we would reach to k = -1 so we set to spatial_order = min(spatial_order,4)
    int l_spatial_order_lo = std::min(std::min(vert_spatial_order, 2*(domhi_z+2-k)), 2*k);

    if (k == 0) {
        zflux_lo = rho_w(i,j,k) * w(i,j,k);
    } else {
        rho_w_avg = 0.5 * (rho_w(i,j,k) + rho_w(i,j,k-1));
        zflux_lo = rho_w_avg * InterpolateInZ(i, j, k  , w, 0, rho_w_avg, l_spatial_order_lo);
    }

    // If k+1 == domhi_z   and spatial_order >= 3 we would reach to k = domhi_z+2 so we set to spatial_order = min(spatial_order,2)
    // If k+1 == domhi_z-1 and spatial_order >= 5 we would reach to k = domhi_z+2 so we set to spatial_order = min(spatial_order,4)
    int l_spatial_order_hi = std::min(std::min(vert_spatial_order, 2*(domhi_z+1-k)), 2*(k+1));

    if (k == domhi_z+1) {
        zflux_hi =  rho_w(i,j,k) * w(i,j,k);
    } else {
        rho_w_avg = 0.5 * (rho_w(i,j,k) + rho_w(i,j,k+1));
        zflux_hi = rho_w_avg * InterpolateInZ(i, j, k+1, w, 0, rho_w_avg, l_spatial_order_hi);
    }

    amrex::Real mfsq = mf_m(i,j,0) * mf_m(i,j,0);

    advectionSrc = (xflux_hi - xflux_lo) * dxInv * mfsq
                 + (yflux_hi - yflux_lo) * dyInv * mfsq
                 + (zflux_hi - zflux_lo) * dzInv;

   return advectionSrc;
}


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForXMom_WENO_N (int i, int j, int k,
                            const amrex::Array4<const amrex::Real>& rho_u, const amrex::Array4<const amrex::Real>& rho_v,
                            const amrex::Array4<const amrex::Real>& rho_w, const amrex::Array4<const amrex::Real>& u,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                            const amrex::Array4<const amrex::Real>& mf_u,
                            const amrex::Array4<const amrex::Real>& mf_v,
                            int spatial_order_WENO)
{
    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    amrex::Real rho_u_avg, rho_v_avg, rho_w_avg;

    amrex::Real xflux_hi; amrex::Real xflux_lo;
    amrex::Real yflux_hi; amrex::Real yflux_lo;
    amrex::Real zflux_hi; amrex::Real zflux_lo;

    amrex::Real mf_u_inv_hi = 1. / mf_u(i+1,j  ,0); amrex::Real mf_u_inv_mid = 1. / mf_u(i  ,j  ,0); amrex::Real mf_u_inv_lo = 1. / mf_u(i-1,j  ,0);
    amrex::Real mf_v_inv_1  = 1. / mf_v(i  ,j+1,0); amrex::Real mf_v_inv_2   = 1. / mf_v(i-1,j+1,0); amrex::Real mf_v_inv_3  = 1. / mf_v(i  ,j  ,0); amrex::Real mf_v_inv_4 = 1. / mf_v(i-1,j  ,0);

    rho_u_avg = 0.5 * (rho_u(i+1, j, k) * mf_u_inv_hi + rho_u(i, j, k) * mf_u_inv_mid);
    xflux_hi = rho_u_avg * InterpolateInX_WENO(i+1, j, k, u, 0, spatial_order_WENO);

    rho_u_avg = 0.5 * (rho_u(i-1, j, k) * mf_u_inv_lo + rho_u(i, j, k) * mf_u_inv_mid);
    xflux_lo = rho_u_avg * InterpolateInX_WENO(i  , j, k, u, 0, spatial_order_WENO);

    rho_v_avg = 0.5 * (rho_v(i, j+1, k) * mf_v_inv_1 + rho_v(i-1, j+1, k) * mf_v_inv_2);
    yflux_hi = rho_v_avg * InterpolateInY_WENO(i, j+1, k, u, 0, spatial_order_WENO);

    rho_v_avg = 0.5 * (rho_v(i, j  , k) * mf_v_inv_3 + rho_v(i-1, j  , k) * mf_v_inv_4);
    yflux_lo = rho_v_avg * InterpolateInY_WENO(i, j  , k, u, 0, spatial_order_WENO);

    rho_w_avg = 0.5 * (rho_w(i, j, k+1) + rho_w(i-1, j, k+1));
    zflux_hi = rho_w_avg * InterpolateInZ_WENO(i, j, k+1, u, 0, spatial_order_WENO);

    rho_w_avg = 0.5 * (rho_w(i, j, k) + rho_w(i-1, j, k));
    zflux_lo = rho_w_avg * InterpolateInZ_WENO(i, j, k  , u, 0, spatial_order_WENO);

    amrex::Real mfsq = mf_u(i,j,0) * mf_u(i,j,0);

    advectionSrc = (xflux_hi - xflux_lo) * dxInv * mfsq
                 + (yflux_hi - yflux_lo) * dyInv * mfsq
                 + (zflux_hi - zflux_lo) * dzInv;

    return advectionSrc;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForYMom_WENO_N (int i, int j, int k,
                            const amrex::Array4<const amrex::Real>& rho_u, const amrex::Array4<const amrex::Real>& rho_v,
                            const amrex::Array4<const amrex::Real>& rho_w, const amrex::Array4<const amrex::Real>& v,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                            const amrex::Array4<const amrex::Real>& mf_u,
                            const amrex::Array4<const amrex::Real>& mf_v,
                            int spatial_order_WENO)
{
    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    amrex::Real rho_u_avg, rho_v_avg, rho_w_avg;

    amrex::Real xflux_hi; amrex::Real xflux_lo;
    amrex::Real yflux_hi; amrex::Real yflux_lo;
    amrex::Real zflux_hi; amrex::Real zflux_lo;

    amrex::Real mf_v_inv_hi = 1. / mf_v(i  ,j+1,0); amrex::Real mf_v_inv_mid = 1. / mf_v(i  ,j  ,0); amrex::Real mf_v_inv_lo = 1. / mf_v(i  ,j-1,0);
    amrex::Real mf_u_inv_1  = 1. / mf_u(i+1,j  ,0); amrex::Real mf_u_inv_2   = 1. / mf_u(i+1,j-1,0); amrex::Real mf_u_inv_3  = 1. / mf_u(i  ,j  ,0); amrex::Real mf_u_inv_4 = 1. / mf_u(i  ,j-1,0);

    rho_u_avg = 0.5*(rho_u(i+1, j, k) * mf_u_inv_1 + rho_u(i+1, j-1, k) * mf_u_inv_2);
    xflux_hi = rho_u_avg * InterpolateInX_WENO(i+1, j, k, v, 0, spatial_order_WENO);

    rho_u_avg = 0.5*(rho_u(i  , j, k) * mf_u_inv_3 + rho_u(i  , j-1, k) * mf_u_inv_4);
    xflux_lo = rho_u_avg * InterpolateInX_WENO(i  , j, k, v, 0, spatial_order_WENO);

    rho_v_avg = 0.5*(rho_v(i, j, k) * mf_v_inv_mid + rho_v(i, j+1, k) * mf_v_inv_hi);
    yflux_hi = rho_v_avg * InterpolateInY_WENO(i, j+1, k, v, 0, spatial_order_WENO);

    rho_v_avg = 0.5*(rho_v(i, j, k) * mf_v_inv_mid + rho_v(i, j-1, k) * mf_v_inv_lo);
    yflux_lo = rho_v_avg * InterpolateInY_WENO(i, j  , k, v, 0, spatial_order_WENO);

    rho_w_avg = 0.5*(rho_w(i, j, k+1) + rho_w(i, j-1, k+1));
    zflux_hi = rho_w_avg * InterpolateInZ_WENO(i, j, k+1, v, 0, spatial_order_WENO);

    rho_w_avg = 0.5*(rho_w(i, j, k) + rho_w(i, j-1, k));
    zflux_lo = rho_w_avg * InterpolateInZ_WENO(i, j, k  , v, 0, spatial_order_WENO);

    amrex::Real mfsq = mf_v(i,j,0) * mf_v(i,j,0);

    advectionSrc = (xflux_hi - xflux_lo) * dxInv * mfsq
                 + (yflux_hi - yflux_lo) * dyInv * mfsq
                 + (zflux_hi - zflux_lo) * dzInv;

   return advectionSrc;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForZMom_WENO_N (int i, int j, int k,
                            const amrex::Array4<const amrex::Real>& rho_u, const amrex::Array4<const amrex::Real>& rho_v,
                            const amrex::Array4<const amrex::Real>& rho_w, const amrex::Array4<const amrex::Real>& w,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                            const amrex::Array4<const amrex::Real>& mf_m,
                            const amrex::Array4<const amrex::Real>& mf_u,
                            const amrex::Array4<const amrex::Real>& mf_v,
                            int spatial_order_WENO, int domhi_z)
{

    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    amrex::Real rho_u_avg, rho_v_avg, rho_w_avg;

    amrex::Real xflux_hi; amrex::Real xflux_lo;
    amrex::Real yflux_hi; amrex::Real yflux_lo;
    amrex::Real zflux_hi; amrex::Real zflux_lo;

    amrex::Real mf_u_inv_hi = 1. / mf_u(i+1,j  ,0); amrex::Real mf_u_inv_lo = 1. / mf_u(i  ,j  ,0);
    amrex::Real mf_v_inv_hi = 1. / mf_v(i  ,j+1,0); amrex::Real mf_v_inv_lo = 1. / mf_v(i  ,j  ,0);

    rho_u_avg = 0.5*(rho_u(i+1, j, k) + rho_u(i+1, j, k-1)) * mf_u_inv_hi;
    xflux_hi = rho_u_avg * InterpolateInX_WENO(i+1, j, k, w, 0, spatial_order_WENO);

    rho_u_avg = 0.5*(rho_u(i  , j, k) + rho_u(i  , j, k-1)) * mf_u_inv_lo;
    xflux_lo = rho_u_avg * InterpolateInX_WENO(i  , j, k, w, 0, spatial_order_WENO);

    rho_v_avg = 0.5*(rho_v(i, j+1, k) + rho_v(i, j+1, k-1)) * mf_v_inv_hi;
    yflux_hi = rho_v_avg * InterpolateInY_WENO(i, j+1, k, w, 0, spatial_order_WENO);

    rho_v_avg = 0.5*(rho_v(i, j  , k) + rho_v(i, j  , k-1)) * mf_v_inv_lo;
    yflux_lo = rho_v_avg * InterpolateInY_WENO(i, j  , k, w, 0, spatial_order_WENO);

    // Constrain to 2nd order near top and bottom
    if (k == 0) {
        zflux_lo = rho_w(i,j,k) * w(i,j,k);
    } else if ((k > domhi_z-2) || (k < 3)) {
        rho_w_avg = 0.5 * (rho_w(i,j,k) + rho_w(i,j,k-1));
        zflux_lo = rho_w_avg * InterpolateInZ(i, j, k  , w, 0, rho_w_avg, 2);
    } else {
        rho_w_avg = 0.5 * (rho_w(i,j,k) + rho_w(i,j,k-1));
        zflux_lo = rho_w_avg * InterpolateInZ_WENO(i, j, k  , w, 0, spatial_order_WENO);
    }

    if (k == domhi_z+1) {
        zflux_hi =  rho_w(i,j,k) * w(i,j,k);
    } else if ((k > domhi_z-2) || (k < 3)) {
        rho_w_avg = 0.5 * (rho_w(i,j,k) + rho_w(i,j,k+1));
        zflux_hi = rho_w_avg * InterpolateInZ(i, j, k+1, w, 0, rho_w_avg, 2);
    } else {
        rho_w_avg = 0.5 * (rho_w(i,j,k) + rho_w(i,j,k+1));
        zflux_hi = rho_w_avg * InterpolateInZ_WENO(i, j, k+1, w, 0, spatial_order_WENO);
    }

    amrex::Real mfsq = mf_m(i,j,0) * mf_m(i,j,0);

    advectionSrc = (xflux_hi - xflux_lo) * dxInv * mfsq
                 + (yflux_hi - yflux_lo) * dyInv * mfsq
                 + (zflux_hi - zflux_lo) * dzInv;

   return advectionSrc;
}
