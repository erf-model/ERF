#include <IndexDefines.H>
#include <TerrainMetrics.H>
#include <Interpolation.H>

/**
 * Function for computing the advective tendency for the x-component of momentum
 * without metric terms and for spatial order > 2
 *
 * @param[in] i,j,k indices of x-face at which to create tendency
 * @param[in] rho_u x-component of momentum
 * @param[in] rho_v y-component of momentum
 * @param[in] rho_w z-component of momentum
 * @param[in] u     x-component of velocity
 * @param[in] cellSizeInv inverse of the mesh spacing
 * @param[in] mf_u map factor on x-faces
 * @param[in] mf_v map factor on y-faces
 * @param[in] horiz_spatial_order spatial order to be used for lateral derivatives if not using WENO (2-6)
 * @param[in] vert_spatial_order spatial order to be used for vertical derivatives if not using WENO (2-6)
 */
template<typename InterpType>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForXMom_N (int i, int j, int k,
                       const amrex::Array4<const amrex::Real>& rho_u,
                       const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& rho_w,
                       InterpType interp_u,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Array4<const amrex::Real>& mf_u,
                       const amrex::Array4<const amrex::Real>& mf_v)
{
    //used when spatial order > 2

    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    amrex::Real rho_u_avg_lo, rho_u_avg_hi;
    amrex::Real rho_v_avg_lo, rho_v_avg_hi;
    amrex::Real rho_w_avg_lo, rho_w_avg_hi;

    amrex::Real xflux_hi; amrex::Real xflux_lo;
    amrex::Real yflux_hi; amrex::Real yflux_lo;
    amrex::Real zflux_hi; amrex::Real zflux_lo;

    amrex::Real mf_u_inv_hi = 1. / mf_u(i+1,j  ,0); amrex::Real mf_u_inv_mid = 1. / mf_u(i  ,j  ,0); amrex::Real mf_u_inv_lo = 1. / mf_u(i-1,j  ,0);
    amrex::Real mf_v_inv_1  = 1. / mf_v(i  ,j+1,0); amrex::Real mf_v_inv_2   = 1. / mf_v(i-1,j+1,0); amrex::Real mf_v_inv_3  = 1. / mf_v(i  ,j  ,0); amrex::Real mf_v_inv_4 = 1. / mf_v(i-1,j  ,0);

    amrex::Real interp_hi(0.), interp_lo(0.);

    rho_u_avg_hi = 0.5 * (rho_u(i+1, j, k) * mf_u_inv_hi + rho_u(i, j, k) * mf_u_inv_mid);
    rho_u_avg_lo = 0.5 * (rho_u(i-1, j, k) * mf_u_inv_lo + rho_u(i, j, k) * mf_u_inv_mid);
    interp_u.InterpolateInX(i,j,k,0,interp_hi,interp_lo,rho_u_avg_hi,rho_u_avg_lo);
    xflux_hi = rho_u_avg_hi * interp_hi;
    xflux_lo = rho_u_avg_lo * interp_lo;

    rho_v_avg_hi = 0.5 * (rho_v(i, j+1, k) * mf_v_inv_1 + rho_v(i-1, j+1, k) * mf_v_inv_2);
    rho_v_avg_lo = 0.5 * (rho_v(i, j  , k) * mf_v_inv_3 + rho_v(i-1, j  , k) * mf_v_inv_4);
    interp_u.InterpolateInY(i,j,k,0,interp_hi,interp_lo,rho_v_avg_hi,rho_v_avg_lo);
    yflux_hi = rho_v_avg_hi * interp_hi;
    yflux_lo = rho_v_avg_lo * interp_lo;

    rho_w_avg_hi = 0.5 * (rho_w(i, j, k+1) + rho_w(i-1, j, k+1));
    rho_w_avg_lo = 0.5 * (rho_w(i, j, k) + rho_w(i-1, j, k));
    interp_u.InterpolateInZ_lo(i,j,k,0,interp_lo,rho_u_avg_lo);
    interp_u.InterpolateInZ_hi(i,j,k,0,interp_hi,rho_u_avg_hi);
    zflux_hi = rho_w_avg_hi * interp_hi;
    zflux_lo = rho_w_avg_lo * interp_lo;

    amrex::Real mfsq = mf_u(i,j,0) * mf_u(i,j,0);

    advectionSrc = (xflux_hi - xflux_lo) * dxInv * mfsq
                 + (yflux_hi - yflux_lo) * dyInv * mfsq
                 + (zflux_hi - zflux_lo) * dzInv;

    return advectionSrc;
}

/**
 * Function for computing the advective tendency for the y-component of momentum
 * without metric terms and for spatial order > 2
 *
 * @param[in] i,j,k indices of y-face at which to create tendency
 * @param[in] rho_u x-component of momentum
 * @param[in] rho_v y-component of momentum
 * @param[in] rho_w z-component of momentum
 * @param[in] v     y-component of velocity
 * @param[in] cellSizeInv inverse of the mesh spacing
 * @param[in] mf_u map factor on x-faces
 * @param[in] mf_v map factor on y-faces
 * @param[in] horiz_spatial_order spatial order to be used for lateral derivatives if not using WENO (2-6)
 * @param[in] vert_spatial_order spatial order to be used for vertical derivatives if not using WENO (2-6)
 */
template<typename InterpType>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForYMom_N (int i, int j, int k,
                       const amrex::Array4<const amrex::Real>& rho_u,
                       const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& rho_w,
                       InterpType interp_v,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Array4<const amrex::Real>& mf_u,
                       const amrex::Array4<const amrex::Real>& mf_v)
{
    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    amrex::Real rho_u_avg_lo, rho_u_avg_hi;
    amrex::Real rho_v_avg_lo, rho_v_avg_hi;
    amrex::Real rho_w_avg_lo, rho_w_avg_hi;

    amrex::Real xflux_hi; amrex::Real xflux_lo;
    amrex::Real yflux_hi; amrex::Real yflux_lo;
    amrex::Real zflux_hi; amrex::Real zflux_lo;

    amrex::Real mf_v_inv_hi = 1. / mf_v(i  ,j+1,0); amrex::Real mf_v_inv_mid = 1. / mf_v(i  ,j  ,0); amrex::Real mf_v_inv_lo = 1. / mf_v(i  ,j-1,0);
    amrex::Real mf_u_inv_1  = 1. / mf_u(i+1,j  ,0); amrex::Real mf_u_inv_2   = 1. / mf_u(i+1,j-1,0); amrex::Real mf_u_inv_3  = 1. / mf_u(i  ,j  ,0); amrex::Real mf_u_inv_4 = 1. / mf_u(i  ,j-1,0);

    amrex::Real interp_hi(0.), interp_lo(0.);

    rho_u_avg_hi = 0.5 * (rho_u(i+1, j, k) * mf_u_inv_1 + rho_u(i+1, j-1, k) * mf_u_inv_2);
    rho_u_avg_lo = 0.5 * (rho_u(i  , j, k) * mf_u_inv_3 + rho_u(i  , j-1, k) * mf_u_inv_4);
    interp_v.InterpolateInX(i,j,k,0,interp_hi,interp_lo,rho_u_avg_hi,rho_u_avg_lo);
    xflux_hi = rho_u_avg_hi * interp_hi;
    xflux_lo = rho_u_avg_lo * interp_lo;

    rho_v_avg_hi = 0.5 * (rho_v(i, j, k) * mf_v_inv_mid + rho_v(i, j+1, k) * mf_v_inv_hi);
    rho_v_avg_lo = 0.5 * (rho_v(i, j, k) * mf_v_inv_mid + rho_v(i, j-1, k) * mf_v_inv_lo);
    interp_v.InterpolateInY(i,j,k,0,interp_hi,interp_lo,rho_v_avg_hi,rho_v_avg_lo);
    yflux_hi = rho_v_avg_hi * interp_hi;
    yflux_lo = rho_v_avg_lo * interp_lo;

    rho_w_avg_hi = 0.5 * (rho_w(i, j, k+1) + rho_w(i, j-1, k+1));
    rho_w_avg_lo = 0.5 * (rho_w(i, j, k  ) + rho_w(i, j-1, k  ));
    interp_v.InterpolateInZ_lo(i,j,k,0,interp_lo,rho_v_avg_lo);
    interp_v.InterpolateInZ_hi(i,j,k,0,interp_hi,rho_v_avg_hi);
    zflux_hi = rho_w_avg_hi * interp_hi;
    zflux_lo = rho_w_avg_lo * interp_lo;

    amrex::Real mfsq = mf_v(i,j,0) * mf_v(i,j,0);

    advectionSrc = (xflux_hi - xflux_lo) * dxInv * mfsq
                 + (yflux_hi - yflux_lo) * dyInv * mfsq
                 + (zflux_hi - zflux_lo) * dzInv;

   return advectionSrc;
}

/**
 * Function for computing the advective tendency for the z-component of momentum
 * without metric terms and for spatial order > 2
 *
 * @param[in] i,j,k indices of z-face at which to create tendency
 * @param[in] rho_u x-component of momentum
 * @param[in] rho_v y-component of momentum
 * @param[in] rho_w z-component of momentum
 * @param[in] w     z-component of velocity
 * @param[in] cellSizeInv inverse of the mesh spacing
 * @param[in] mf_m map factor on cell centers
 * @param[in] mf_u map factor on x-faces
 * @param[in] mf_v map factor on y-faces
 * @param[in] horiz_spatial_order spatial order to be used for lateral derivatives if not using WENO (2-6)
 * @param[in] vert_spatial_order spatial order to be used for vertical derivatives if not using WENO (2-6)
 * @param[in] domhi_z maximum k value in the domain
 */
template<typename InterpType, typename WallInterpType>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForZMom_N (int i, int j, int k,
                       const amrex::Array4<const amrex::Real>& rho_u,
                       const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& rho_w,
                       const amrex::Array4<const amrex::Real>& w,
                       InterpType interp_w,
                       WallInterpType interp_w_wall,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Array4<const amrex::Real>& mf_m,
                       const amrex::Array4<const amrex::Real>& mf_u,
                       const amrex::Array4<const amrex::Real>& mf_v,
                       const int vert_spatial_order,
                       const int domhi_z)
{

    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    amrex::Real rho_u_avg_lo, rho_u_avg_hi;
    amrex::Real rho_v_avg_lo, rho_v_avg_hi;
    amrex::Real rho_w_avg_lo, rho_w_avg_hi;

    amrex::Real xflux_hi; amrex::Real xflux_lo;
    amrex::Real yflux_hi; amrex::Real yflux_lo;
    amrex::Real zflux_hi; amrex::Real zflux_lo;

    amrex::Real mf_u_inv_hi = 1. / mf_u(i+1,j  ,0); amrex::Real mf_u_inv_lo = 1. / mf_u(i  ,j  ,0);
    amrex::Real mf_v_inv_hi = 1. / mf_v(i  ,j+1,0); amrex::Real mf_v_inv_lo = 1. / mf_v(i  ,j  ,0);

    amrex::Real interp_hi(0.), interp_lo(0.);

    rho_u_avg_hi = 0.5 * (rho_u(i+1, j, k) + rho_u(i+1, j, k-1)) * mf_u_inv_hi;
    rho_u_avg_lo = 0.5 * (rho_u(i  , j, k) + rho_u(i  , j, k-1)) * mf_u_inv_lo;
    interp_w.InterpolateInX(i,j,k,0,interp_hi,interp_lo,rho_u_avg_hi,rho_u_avg_lo);
    xflux_hi = rho_u_avg_hi * interp_hi;
    xflux_lo = rho_u_avg_lo * interp_lo;

    rho_v_avg_hi = 0.5 * (rho_v(i, j+1, k) + rho_v(i, j+1, k-1)) * mf_v_inv_hi;
    rho_v_avg_lo = 0.5 * (rho_v(i, j  , k) + rho_v(i, j  , k-1)) * mf_v_inv_lo;
    interp_w.InterpolateInY(i,j,k,0,interp_hi,interp_lo,rho_v_avg_hi,rho_v_avg_lo);
    yflux_hi = rho_v_avg_hi * interp_hi;
    yflux_lo = rho_v_avg_lo * interp_lo;

    // If k+1 == domhi_z   and spatial_order >= 3 we would reach to k = domhi_z+2 so we set to spatial_order = min(spatial_order,2)
    // If k+1 == domhi_z-1 and spatial_order >= 5 we would reach to k = domhi_z+2 so we set to spatial_order = min(spatial_order,4)
    int l_spatial_order_hi = std::min(std::min(vert_spatial_order, 2*(domhi_z+1-k)), 2*(k+1));

    if (k == domhi_z+1) {
        zflux_hi =  rho_w(i,j,k) * w(i,j,k);
    } else if (l_spatial_order_hi < vert_spatial_order) {
        rho_w_avg_hi = 0.5 * (rho_w(i,j,k) + rho_w(i,j,k+1));
        interp_w_wall.InterpolateInZ_hi(i,j,k,0,interp_hi,rho_w_avg_hi,l_spatial_order_hi);
        zflux_hi = rho_w_avg_hi * interp_hi;
    } else {
        rho_w_avg_hi = 0.5 * (rho_w(i,j,k) + rho_w(i,j,k+1));
        interp_w.InterpolateInZ_hi(i,j,k,0,interp_hi,rho_w_avg_hi);
        zflux_hi = rho_w_avg_hi * interp_hi;
    }

    // If k == 1 and spatial_order >= 3 we would reach to k = -1 so we set to spatial_order = min(spatial_order,2)
    // If k == 2 and spatial_order >= 5 we would reach to k = -1 so we set to spatial_order = min(spatial_order,4)
    int l_spatial_order_lo = std::min(std::min(vert_spatial_order, 2*(domhi_z+2-k)), 2*k);

    if (k == 0) {
        zflux_lo = rho_w(i,j,k) * w(i,j,k);
    } else if (l_spatial_order_lo < vert_spatial_order) {
        rho_w_avg_lo = 0.5 * (rho_w(i,j,k) + rho_w(i,j,k-1));
        interp_w_wall.InterpolateInZ_lo(i,j,k,0,interp_lo,rho_w_avg_lo,l_spatial_order_lo);
        zflux_lo = rho_w_avg_lo * interp_lo;
    } else {
        rho_w_avg_lo = 0.5 * (rho_w(i,j,k) + rho_w(i,j,k-1));
        interp_w.InterpolateInZ_lo(i,j,k,0,interp_lo,rho_w_avg_lo);
        zflux_lo = rho_w_avg_lo * interp_lo;
    }

    amrex::Real mfsq = mf_m(i,j,0) * mf_m(i,j,0);

    advectionSrc = (xflux_hi - xflux_lo) * dxInv * mfsq
                 + (yflux_hi - yflux_lo) * dyInv * mfsq
                 + (zflux_hi - zflux_lo) * dzInv;

   return advectionSrc;
}

/**
 * Wrapper function for computing the advective tendency w/ spatial order > 2.
 */
template<typename InterpType, typename WallInterpType>
void
AdvectionSrcForMomWrapper_N(const amrex::Box& bxx, const amrex::Box& bxy, const amrex::Box& bxz,
                            const amrex::Array4<amrex::Real>& rho_u_rhs,
                            const amrex::Array4<amrex::Real>& rho_v_rhs,
                            const amrex::Array4<amrex::Real>& rho_w_rhs,
                            InterpType interp_u,
                            InterpType interp_v,
                            InterpType interp_w,
                            WallInterpType interp_w_wall,
                            const amrex::Array4<const amrex::Real>& rho_u,
                            const amrex::Array4<const amrex::Real>& rho_v,
                            const amrex::Array4<const amrex::Real>& rho_w,
                            const amrex::Array4<const amrex::Real>& w,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                            const amrex::Array4<const amrex::Real>& mf_m,
                            const amrex::Array4<const amrex::Real>& mf_u,
                            const amrex::Array4<const amrex::Real>& mf_v,
                            const int vert_spatial_order,
                            const int domhi_z)
{
    amrex::ParallelFor(bxx, bxy, bxz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        rho_u_rhs(i, j, k) = -AdvectionSrcForXMom_N(i, j, k, rho_u, rho_v, rho_w,
                                                    interp_u, cellSizeInv, mf_u, mf_v);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        rho_v_rhs(i, j, k) = -AdvectionSrcForYMom_N(i, j, k, rho_u, rho_v, rho_w,
                                                    interp_v, cellSizeInv, mf_u, mf_v);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        rho_w_rhs(i, j, k) = -AdvectionSrcForZMom_N(i, j, k, rho_u, rho_v, rho_w, w,
                                                    interp_w, interp_w_wall,
                                                    cellSizeInv, mf_m, mf_u, mf_v,
                                                    vert_spatial_order, domhi_z);
    });
}
