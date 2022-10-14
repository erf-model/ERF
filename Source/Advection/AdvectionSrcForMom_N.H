#include <IndexDefines.H>
#include <TerrainMetrics.H>
#include <Interpolation.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForXMom_N (int i, int j, int k,
                       const amrex::Array4<const amrex::Real>& rho_u, const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& rho_w, const amrex::Array4<const amrex::Real>& u,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Array4<const amrex::Real>& mf_m,
                       const amrex::Array4<const amrex::Real>& mf_u,
                       const amrex::Array4<const amrex::Real>& mf_v,
                       int spatial_order)
{
    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    amrex::Real rho_u_avg, rho_v_avg, rho_w_avg;

    amrex::Real xflux_hi; amrex::Real xflux_lo;
    amrex::Real yflux_hi; amrex::Real yflux_lo;
    amrex::Real zflux_hi; amrex::Real zflux_lo;

    if (spatial_order == 2) {
        xflux_hi = 0.25 * (rho_u(i, j  , k) + rho_u(i+1, j  , k)) * (u(i+1,j,k) + u(i,j,k)) * mf_m(i  ,j,0);
        xflux_lo = 0.25 * (rho_u(i, j  , k) + rho_u(i-1, j  , k)) * (u(i-1,j,k) + u(i,j,k)) * mf_m(i-1,j,0);

        yflux_hi = 0.25 * (rho_v(i, j+1, k) + rho_v(i-1, j+1, k)) * (u(i,j+1,k) + u(i,j,k)) *
                   0.5 *  (mf_u(i,j,0) + mf_u(i,j+1,0));
        yflux_lo = 0.25 * (rho_v(i, j  , k) + rho_v(i-1, j  , k)) * (u(i,j-1,k) + u(i,j,k)) *
                   0.5 *  (mf_u(i,j,0) + mf_u(i,j-1,0));

        zflux_hi = 0.25 * (rho_w(i, j, k+1) + rho_w(i-1, j, k+1)) * (u(i,j,k+1) + u(i,j,k)); //mapfacs are moved from here to advSource term
        zflux_lo = 0.25 * (rho_w(i, j, k  ) + rho_w(i-1, j, k  )) * (u(i,j,k-1) + u(i,j,k)); //to reduce number of operations

    } else {
        rho_u_avg = 0.5 * (rho_u(i+1, j, k) + rho_u(i, j, k));
        xflux_hi = rho_u_avg * InterpolateInX(i+1, j, k, u, 0, rho_u_avg, spatial_order) * mf_m(i  ,j,0);

        rho_u_avg = 0.5 * (rho_u(i-1, j, k) + rho_u(i, j, k));
        xflux_lo = rho_u_avg * InterpolateInX(i  , j, k, u, 0, rho_u_avg, spatial_order) * mf_m(i-1,j,0);

        rho_v_avg = 0.5 * (rho_v(i, j+1, k) + rho_v(i-1, j+1, k));
        yflux_hi = rho_v_avg * InterpolateInY(i, j+1, k, u, 0, rho_v_avg, spatial_order) *
                   0.5 *  (mf_u(i,j,0) + mf_u(i,j+1,0));

        rho_v_avg = 0.5 * (rho_v(i, j  , k) + rho_v(i-1, j  , k));
        yflux_lo = rho_v_avg * InterpolateInY(i, j  , k, u, 0, rho_v_avg, spatial_order) *
                   0.5 *  (mf_u(i,j,0) + mf_u(i,j-1,0));

        rho_w_avg = 0.5 * (rho_w(i, j, k+1) + rho_w(i-1, j, k+1));
        zflux_hi = rho_w_avg * InterpolateInZ(i, j, k+1, u, 0, rho_w_avg, spatial_order);

        rho_w_avg = 0.5 * (rho_w(i, j, k) + rho_w(i-1, j, k));
        zflux_lo = rho_w_avg * InterpolateInZ(i, j, k  , u, 0, rho_w_avg, spatial_order);
    }

    advectionSrc = (xflux_hi - xflux_lo) * dxInv 
                 + (yflux_hi - yflux_lo) * dyInv 
                 + (zflux_hi - zflux_lo) * dzInv;
    advectionSrc *= mf_u(i,j,0);

    return advectionSrc;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForYMom_N (int i, int j, int k,
                       const amrex::Array4<const amrex::Real>& rho_u, const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& rho_w, const amrex::Array4<const amrex::Real>& v,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Array4<const amrex::Real>& mf_m,
                       const amrex::Array4<const amrex::Real>& mf_u,
                       const amrex::Array4<const amrex::Real>& mf_v,
                       int spatial_order)
{
    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    amrex::Real rho_u_avg, rho_v_avg, rho_w_avg;

    amrex::Real xflux_hi; amrex::Real xflux_lo;
    amrex::Real yflux_hi; amrex::Real yflux_lo;
    amrex::Real zflux_hi; amrex::Real zflux_lo;

    if (spatial_order == 2) {

        xflux_hi = 0.25 * (rho_u(i+1, j, k) + rho_u(i+1, j-1, k)) * (v(i+1,j,k) + v(i,j,k)) *
                   0.5 *  (mf_v(i,j,0) + mf_v(i+1,j,0));
        xflux_lo = 0.25 * (rho_u(i  , j, k) + rho_u(i  , j-1, k)) * (v(i-1,j,k) + v(i,j,k)) *
                   0.5 *  (mf_v(i,j,0) + mf_v(i-1,j,0));

        yflux_hi = 0.25 * (rho_v(i, j, k  ) + rho_v(i  , j+1, k)) * (v(i,j+1,k) + v(i,j,k)) * mf_m(i,j  ,0);
        yflux_lo = 0.25 * (rho_v(i, j, k  ) + rho_v(i  , j-1, k)) * (v(i,j-1,k) + v(i,j,k)) * mf_m(i,j-1,0);

        zflux_hi = 0.25 * (rho_w(i, j, k+1) + rho_w(i, j-1, k+1)) * (v(i,j,k+1) + v(i,j,k)); //mapfacs are moved from here to advSource term
        zflux_lo = 0.25 * (rho_w(i, j, k  ) + rho_w(i, j-1, k  )) * (v(i,j,k-1) + v(i,j,k)); //to reduce number of operations

    } else {

        rho_u_avg = 0.5*(rho_u(i+1, j, k) + rho_u(i+1, j-1, k));
        xflux_hi = rho_u_avg * InterpolateInX(i+1, j, k, v, 0, rho_u_avg, spatial_order) *
                   0.5 *  (mf_v(i,j,0) + mf_v(i+1,j,0));

        rho_u_avg = 0.5*(rho_u(i  , j, k) + rho_u(i  , j-1, k));
        xflux_lo = rho_u_avg * InterpolateInX(i  , j, k, v, 0, rho_u_avg, spatial_order) *
                   0.5 *  (mf_v(i,j,0) + mf_v(i-1,j,0));

        rho_v_avg = 0.5*(rho_v(i, j, k) + rho_v(i, j+1, k));
        yflux_hi = rho_v_avg * InterpolateInY(i, j+1, k, v, 0, rho_v_avg, spatial_order) * mf_m(i,j  ,0);

        rho_v_avg = 0.5*(rho_v(i, j, k) + rho_v(i, j-1, k));
        yflux_lo = rho_v_avg * InterpolateInY(i, j  , k, v, 0, rho_v_avg, spatial_order) * mf_m(i,j-1,0);

        rho_w_avg = 0.5*(rho_w(i, j, k+1) + rho_w(i, j-1, k+1));
        zflux_hi = rho_w_avg * InterpolateInZ(i, j, k+1, v, 0, rho_w_avg, spatial_order);

        rho_w_avg = 0.5*(rho_w(i, j, k) + rho_w(i, j-1, k));
        zflux_lo = rho_w_avg * InterpolateInZ(i, j, k  , v, 0, rho_w_avg, spatial_order);
    }

    advectionSrc = (xflux_hi - xflux_lo) * dxInv
                 + (yflux_hi - yflux_lo) * dyInv
                 + (zflux_hi - zflux_lo) * dzInv;
    advectionSrc *= mf_v(i,j,0);

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
                       int spatial_order, int domhi_z)
{
    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    amrex::Real rho_u_avg, rho_v_avg, rho_w_avg;

    amrex::Real xflux_hi; amrex::Real xflux_lo;
    amrex::Real yflux_hi; amrex::Real yflux_lo;
    amrex::Real zflux_hi; amrex::Real zflux_lo;

    if (spatial_order == 2) {

        xflux_hi = 0.25*(rho_u(i+1, j, k) + rho_u(i+1, j, k-1)) * (w(i+1,j,k) + w(i,j,k)) * mf_u(i+1,j,0);
        xflux_lo = 0.25*(rho_u(i  , j, k) + rho_u(i  , j, k-1)) * (w(i-1,j,k) + w(i,j,k)) * mf_u(i  ,j,0);

        yflux_hi = 0.25*(rho_v(i, j+1, k) + rho_v(i, j+1, k-1)) * (w(i,j+1,k) + w(i,j,k)) * mf_v(i,j+1,0);
        yflux_lo = 0.25*(rho_v(i, j  , k) + rho_v(i, j  , k-1)) * (w(i,j-1,k) + w(i,j,k)) * mf_v(i,j  ,0);

        zflux_lo = (k == 0) ? rho_w(i,j,k) * w(i,j,k) :
            0.25 * (rho_w(i,j,k) + rho_w(i,j,k-1)) * (w(i,j,k) + w(i,j,k-1)); //mapfacs are moved from here to advSource term
                                                                              //to reduce number of operations
        zflux_hi = (k == domhi_z+1) ? rho_w(i,j,k) * w(i,j,k) :
            0.25 * (rho_w(i,j,k) + rho_w(i,j,k+1)) * (w(i,j,k) + w(i,j,k+1));

    } else {

        rho_u_avg = 0.5*(rho_u(i+1, j, k) + rho_u(i+1, j, k-1));
        xflux_hi = rho_u_avg * InterpolateInX(i+1, j, k, w, 0, rho_u_avg, spatial_order) * mf_u(i+1,j,0);

        rho_u_avg = 0.5*(rho_u(i  , j, k) + rho_u(i  , j, k-1));
        xflux_lo = rho_u_avg * InterpolateInX(i  , j, k, w, 0, rho_u_avg, spatial_order) * mf_u(i  ,j,0);

        rho_v_avg = 0.5*(rho_v(i, j+1, k) + rho_v(i, j+1, k-1));
        yflux_hi = rho_v_avg * InterpolateInY(i, j+1, k, w, 0, rho_v_avg, spatial_order) * mf_v(i,j+1,0);

        rho_v_avg = 0.5*(rho_v(i, j  , k) + rho_v(i, j  , k-1));
        yflux_lo = rho_v_avg * InterpolateInY(i, j  , k, w, 0, rho_v_avg, spatial_order) * mf_v(i,j  ,0);

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
    advectionSrc *= mf_m(i,j,0);

    return advectionSrc;
}
