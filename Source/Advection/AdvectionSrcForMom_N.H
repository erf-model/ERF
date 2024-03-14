#include <IndexDefines.H>
#include <Interpolation.H>

/**
 * Function for computing the advective tendency for the x-component of momentum
 * without metric terms and for higher-order stencils
 *
 * @param[in] i,j,k indices of x-face at which to create tendency
 * @param[in] rho_u x-component of momentum
 * @param[in] rho_v y-component of momentum
 * @param[in] rho_w z-component of momentum
 * @param[in] u     x-component of velocity
 * @param[in] cellSizeInv inverse of the mesh spacing
 * @param[in] mf_u map factor on x-faces
 * @param[in] mf_v map factor on y-faces
 */
template<typename InterpType_H, typename InterpType_V>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForXMom_N (int i, int j, int k,
                       const amrex::Array4<const amrex::Real>& rho_u,
                       const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& rho_w,
                       InterpType_H interp_u_h,
                       InterpType_V interp_u_v,
                       const amrex::Real upw_frac_h,
                       const amrex::Real upw_frac_v,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Array4<const amrex::Real>& mf_u_inv,
                       const amrex::Array4<const amrex::Real>& mf_v_inv)
{
    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    amrex::Real rho_u_avg_lo, rho_u_avg_hi;
    amrex::Real rho_v_avg_lo, rho_v_avg_hi;
    amrex::Real rho_w_avg_lo, rho_w_avg_hi;

    amrex::Real xflux_hi; amrex::Real xflux_lo;
    amrex::Real yflux_hi; amrex::Real yflux_lo;
    amrex::Real zflux_hi; amrex::Real zflux_lo;

    amrex::Real interp_hi(0.), interp_lo(0.);

    rho_u_avg_hi = 0.5 * (rho_u(i+1, j, k) * mf_u_inv(i+1,j,0) + rho_u(i, j, k) * mf_u_inv(i,j,0));
    interp_u_h.InterpolateInX(i+1,j,k,0,interp_hi,rho_u_avg_hi,upw_frac_h);
    xflux_hi = rho_u_avg_hi * interp_hi;

    rho_u_avg_lo = 0.5 * (rho_u(i-1, j, k) * mf_u_inv(i-1,j,0) + rho_u(i, j, k) * mf_u_inv(i,j,0));
    interp_u_h.InterpolateInX(i,j,k,0,interp_lo,rho_u_avg_lo,upw_frac_h);
    xflux_lo = rho_u_avg_lo * interp_lo;

    rho_v_avg_hi = 0.5 * (rho_v(i, j+1, k) * mf_v_inv(i,j+1,0) + rho_v(i-1, j+1, k) * mf_v_inv(i-1,j+1,0));
    interp_u_h.InterpolateInY(i,j+1,k,0,interp_hi,rho_v_avg_hi,upw_frac_h);
    yflux_hi = rho_v_avg_hi * interp_hi;

    rho_v_avg_lo = 0.5 * (rho_v(i, j  , k) * mf_v_inv(i,j  ,0) + rho_v(i-1, j  , k) * mf_v_inv(i-1,j  ,0));
    interp_u_h.InterpolateInY(i,j,k,0,interp_lo,rho_v_avg_lo,upw_frac_h);
    yflux_lo = rho_v_avg_lo * interp_lo;

    rho_w_avg_hi = 0.5 * (rho_w(i, j, k+1) + rho_w(i-1, j, k+1));
    interp_u_v.InterpolateInZ(i,j,k+1,0,interp_hi,rho_w_avg_hi,upw_frac_v);
    zflux_hi = rho_w_avg_hi * interp_hi;

    rho_w_avg_lo = 0.5 * (rho_w(i, j, k  ) + rho_w(i-1, j, k  ));
    interp_u_v.InterpolateInZ(i,j,k,0,interp_lo,rho_w_avg_lo,upw_frac_v);
    zflux_lo = rho_w_avg_lo * interp_lo;

    amrex::Real mfsq = 1 / (mf_u_inv(i,j,0) * mf_u_inv(i,j,0));

    advectionSrc = (xflux_hi - xflux_lo) * dxInv * mfsq
                 + (yflux_hi - yflux_lo) * dyInv * mfsq
                 + (zflux_hi - zflux_lo) * dzInv;

    return advectionSrc;
}

/**
 * Function for computing the advective tendency for the y-component of momentum
 * without metric terms and for higher-order stencils
 *
 * @param[in] i,j,k indices of y-face at which to create tendency
 * @param[in] rho_u x-component of momentum
 * @param[in] rho_v y-component of momentum
 * @param[in] rho_w z-component of momentum
 * @param[in] v     y-component of velocity
 * @param[in] cellSizeInv inverse of the mesh spacing
 * @param[in] mf_u map factor on x-faces
 * @param[in] mf_v map factor on y-faces
 */
template<typename InterpType_H, typename InterpType_V>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForYMom_N (int i, int j, int k,
                       const amrex::Array4<const amrex::Real>& rho_u,
                       const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& rho_w,
                       InterpType_H interp_v_h,
                       InterpType_V interp_v_v,
                       const amrex::Real upw_frac_h,
                       const amrex::Real upw_frac_v,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Array4<const amrex::Real>& mf_u_inv,
                       const amrex::Array4<const amrex::Real>& mf_v_inv)
{
    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    amrex::Real rho_u_avg_lo, rho_u_avg_hi;
    amrex::Real rho_v_avg_lo, rho_v_avg_hi;
    amrex::Real rho_w_avg_lo, rho_w_avg_hi;

    amrex::Real xflux_hi; amrex::Real xflux_lo;
    amrex::Real yflux_hi; amrex::Real yflux_lo;
    amrex::Real zflux_hi; amrex::Real zflux_lo;

    amrex::Real interp_hi(0.), interp_lo(0.);

    rho_u_avg_hi = 0.5 * (rho_u(i+1, j, k) * mf_u_inv(i+1,j  ,0) + rho_u(i+1, j-1, k) * mf_u_inv(i+1,j-1,0));
    interp_v_h.InterpolateInX(i+1,j,k,0,interp_hi,rho_u_avg_hi,upw_frac_h);
    xflux_hi = rho_u_avg_hi * interp_hi;

    rho_u_avg_lo = 0.5 * (rho_u(i  , j, k) * mf_u_inv(i  ,j  ,0) + rho_u(i  , j-1, k) * mf_u_inv(i  ,j-1,0));
    interp_v_h.InterpolateInX(i,j,k,0,interp_lo,rho_u_avg_lo,upw_frac_h);
    xflux_lo = rho_u_avg_lo * interp_lo;

    rho_v_avg_hi = 0.5 * (rho_v(i, j, k) * mf_v_inv(i  ,j  ,0) + rho_v(i, j+1, k) * mf_v_inv(i  ,j+1,0));
    interp_v_h.InterpolateInY(i,j+1,k,0,interp_hi,rho_v_avg_hi,upw_frac_h);
    yflux_hi = rho_v_avg_hi * interp_hi;

    rho_v_avg_lo = 0.5 * (rho_v(i, j, k) * mf_v_inv(i  ,j  ,0) + rho_v(i, j-1, k) * mf_v_inv(i  ,j-1,0));
    interp_v_h.InterpolateInY(i,j,k,0,interp_lo,rho_v_avg_lo,upw_frac_h);
    yflux_lo = rho_v_avg_lo * interp_lo;

    rho_w_avg_hi = 0.5 * (rho_w(i, j, k+1) + rho_w(i, j-1, k+1));
    interp_v_v.InterpolateInZ(i,j,k+1,0,interp_hi,rho_w_avg_hi,upw_frac_v);
    zflux_hi = rho_w_avg_hi * interp_hi;

    rho_w_avg_lo = 0.5 * (rho_w(i, j, k  ) + rho_w(i, j-1, k  ));
    interp_v_v.InterpolateInZ(i,j,k  ,0,interp_lo,rho_w_avg_lo,upw_frac_v);
    zflux_lo = rho_w_avg_lo * interp_lo;

    amrex::Real mfsq = 1 / (mf_v_inv(i,j,0) * mf_v_inv(i,j,0));

    advectionSrc = (xflux_hi - xflux_lo) * dxInv * mfsq
                 + (yflux_hi - yflux_lo) * dyInv * mfsq
                 + (zflux_hi - zflux_lo) * dzInv;

   return advectionSrc;
}

/**
 * Function for computing the advective tendency for the z-component of momentum
 * without metric terms and for higher-order stencils
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
 * @param[in] domhi_z maximum k value in the domain
 */
template<typename InterpType_H, typename InterpType_V, typename WallInterpType>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForZMom_N (int i, int j, int k,
                       const amrex::Array4<const amrex::Real>& rho_u,
                       const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& rho_w,
                       const amrex::Array4<const amrex::Real>& w,
                       InterpType_H   interp_w_h,
                       InterpType_V   interp_w_v,
                       WallInterpType interp_w_wall,
                       const amrex::Real upw_frac_h,
                       const amrex::Real upw_frac_v,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Array4<const amrex::Real>& mf_m,
                       const amrex::Array4<const amrex::Real>& mf_u_inv,
                       const amrex::Array4<const amrex::Real>& mf_v_inv,
                       const AdvType vert_adv_type,
                       const int lo_z_face, const int hi_z_face)
{

    amrex::Real advectionSrc;
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    amrex::Real rho_u_avg_lo, rho_u_avg_hi;
    amrex::Real rho_v_avg_lo, rho_v_avg_hi;
    amrex::Real rho_w_avg_lo, rho_w_avg_hi;

    amrex::Real xflux_hi; amrex::Real xflux_lo;
    amrex::Real yflux_hi; amrex::Real yflux_lo;
    amrex::Real zflux_hi; amrex::Real zflux_lo;

    amrex::Real interp_hi(0.), interp_lo(0.);

    rho_u_avg_hi = 0.5 * (rho_u(i+1, j, k) + rho_u(i+1, j, k-1)) * mf_u_inv(i+1,j  ,0);
    interp_w_h.InterpolateInX(i+1,j,k,0,interp_hi,rho_u_avg_hi,upw_frac_h);
    xflux_hi = rho_u_avg_hi * interp_hi;

    rho_u_avg_lo = 0.5 * (rho_u(i  , j, k) + rho_u(i  , j, k-1)) * mf_u_inv(i  ,j  ,0);
    interp_w_h.InterpolateInX(i,j,k,0,interp_lo,rho_u_avg_lo,upw_frac_h);
    xflux_lo = rho_u_avg_lo * interp_lo;

    rho_v_avg_hi = 0.5 * (rho_v(i, j+1, k) + rho_v(i, j+1, k-1)) * mf_v_inv(i  ,j+1,0);
    interp_w_h.InterpolateInY(i,j+1,k,0,interp_hi,rho_v_avg_hi,upw_frac_h);
    yflux_hi = rho_v_avg_hi * interp_hi;

    rho_v_avg_lo = 0.5 * (rho_v(i, j  , k) + rho_v(i, j  , k-1)) * mf_v_inv(i  ,j  ,0);
    interp_w_h.InterpolateInY(i,j,k,0,interp_lo,rho_v_avg_lo,upw_frac_h);
    yflux_lo = rho_v_avg_lo * interp_lo;

    // int l_spatial_order_hi = std::min(std::min(vert_spatial_order, 2*(hi_z_face-k)), 2*(k+1));
    // If k == hi_z_face-1, l_spatial_order_hi = 2
    // If k == hi_z_face-2, l_spatial_order_hi = std::min(vert_spatial_order, 4);
    // If k == lo_z+face+1, l_spatial_order_hi = std::min(vert_spatial_order, 4);

    if (k == hi_z_face) {
        zflux_hi =  rho_w(i,j,k) * w(i,j,k);
    } else {
        rho_w_avg_hi = 0.5 * (rho_w(i,j,k) + rho_w(i,j,k+1));
        if (k == hi_z_face-1)
        {
            interp_w_wall.InterpolateInZ(i,j,k+1,0,interp_hi,rho_w_avg_hi,0,AdvType::Centered_2nd);
        } else if (k == hi_z_face-2 || k == lo_z_face+1) {
            if (vert_adv_type != AdvType::Centered_2nd && vert_adv_type != AdvType::Upwind_3rd) {
               interp_w_wall.InterpolateInZ(i,j,k+1,0,interp_hi,rho_w_avg_hi,0,AdvType::Centered_4th);
            } else {
               interp_w_wall.InterpolateInZ(i,j,k+1,0,interp_hi,rho_w_avg_hi,upw_frac_v,vert_adv_type);
            }
        } else {
            interp_w_v.InterpolateInZ(i,j,k+1,0,interp_hi,rho_w_avg_hi,upw_frac_v);
        }
        zflux_hi = rho_w_avg_hi * interp_hi;
    }

    // int l_spatial_order_lo = std::min(std::min(vert_spatial_order, 2*(hi_z_face+1-k)), 2*k);
    // If k == lo_z_face+1, l_spatial_order_hi = 2
    // If k == lo_z_face+2, l_spatial_order_hi = std::min(vert_spatial_order, 4);
    // If k == hi_z_face-1, l_spatial_order_hi = std::min(vert_spatial_order, 4);

    if (k == lo_z_face) {
        zflux_lo = rho_w(i,j,k) * w(i,j,k);
    } else {
        rho_w_avg_lo = 0.5 * (rho_w(i,j,k) + rho_w(i,j,k-1));
        if (k == lo_z_face+1) {
            interp_w_wall.InterpolateInZ(i,j,k,0,interp_lo,rho_w_avg_lo,0,AdvType::Centered_2nd);
        } else if (k == lo_z_face+2 || k == hi_z_face-1) {
            if (vert_adv_type != AdvType::Centered_2nd && vert_adv_type != AdvType::Upwind_3rd) {
                interp_w_wall.InterpolateInZ(i,j,k,0,interp_lo,rho_w_avg_lo,0,AdvType::Centered_4th);
            } else {
                interp_w_wall.InterpolateInZ(i,j,k,0,interp_lo,rho_w_avg_lo,upw_frac_v,vert_adv_type);
            }
        } else {
            interp_w_v.InterpolateInZ(i,j,k,0,interp_lo,rho_w_avg_lo,upw_frac_v);
        }
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
template<typename InterpType_H, typename InterpType_V, typename WallInterpType>
void
AdvectionSrcForMomWrapper_N (const amrex::Box& bxx, const amrex::Box& bxy, const amrex::Box& bxz,
                             const amrex::Array4<amrex::Real>& rho_u_rhs,
                             const amrex::Array4<amrex::Real>& rho_v_rhs,
                             const amrex::Array4<amrex::Real>& rho_w_rhs,
                             const amrex::Array4<const amrex::Real>& rho_u,
                             const amrex::Array4<const amrex::Real>& rho_v,
                             const amrex::Array4<const amrex::Real>& rho_w,
                             const amrex::Array4<const amrex::Real>& u,
                             const amrex::Array4<const amrex::Real>& v,
                             const amrex::Array4<const amrex::Real>& w,
                             const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                             const amrex::Array4<const amrex::Real>& mf_m,
                             const amrex::Array4<const amrex::Real>& mf_u_inv,
                             const amrex::Array4<const amrex::Real>& mf_v_inv,
                             const amrex::Real upw_frac_h,
                             const amrex::Real upw_frac_v,
                             const AdvType vert_adv_type,
                             const int lo_z_face, const int hi_z_face)
{
    // Instantiate the appropriate structs
    InterpType_H interp_u_h(u); InterpType_V interp_u_v(u); // X-MOM
    InterpType_H interp_v_h(v); InterpType_V interp_v_v(v); // Y-MOM
    InterpType_H interp_w_h(w); InterpType_V interp_w_v(w); // Z-MOM
    WallInterpType interp_w_wall(w); // Z-MOM @ wall

    amrex::ParallelFor(bxx, bxy, bxz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        rho_u_rhs(i, j, k) = -AdvectionSrcForXMom_N(i, j, k, rho_u, rho_v, rho_w,
                                                    interp_u_h, interp_u_v,
                                                    upw_frac_h, upw_frac_v,
                                                    cellSizeInv, mf_u_inv, mf_v_inv);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        rho_v_rhs(i, j, k) = -AdvectionSrcForYMom_N(i, j, k, rho_u, rho_v, rho_w,
                                                    interp_v_h, interp_v_v,
                                                    upw_frac_h, upw_frac_v,
                                                    cellSizeInv, mf_u_inv, mf_v_inv);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        rho_w_rhs(i, j, k) = -AdvectionSrcForZMom_N(i, j, k, rho_u, rho_v, rho_w, w,
                                                    interp_w_h, interp_w_v, interp_w_wall,
                                                    upw_frac_h, upw_frac_v,
                                                    cellSizeInv, mf_m, mf_u_inv, mf_v_inv,
                                                    vert_adv_type, lo_z_face, hi_z_face);
    });
}

/**
 * Wrapper function for computing the advective tendency w/ spatial order > 2.
 */
template<typename InterpType_H>
void
AdvectionSrcForMomVert_N (const amrex::Box& bxx, const amrex::Box& bxy, const amrex::Box& bxz,
                          const amrex::Array4<amrex::Real>& rho_u_rhs,
                          const amrex::Array4<amrex::Real>& rho_v_rhs,
                          const amrex::Array4<amrex::Real>& rho_w_rhs,
                          const amrex::Array4<const amrex::Real>& rho_u,
                          const amrex::Array4<const amrex::Real>& rho_v,
                          const amrex::Array4<const amrex::Real>& rho_w,
                          const amrex::Array4<const amrex::Real>& u,
                          const amrex::Array4<const amrex::Real>& v,
                          const amrex::Array4<const amrex::Real>& w,
                          const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                          const amrex::Array4<const amrex::Real>& mf_m,
                          const amrex::Array4<const amrex::Real>& mf_u_inv,
                          const amrex::Array4<const amrex::Real>& mf_v_inv,
                          const amrex::Real upw_frac_h,
                          const amrex::Real upw_frac_v,
                          const AdvType vert_adv_type,
                          const int lo_z_face, const int hi_z_face)
{
    if (vert_adv_type == AdvType::Centered_2nd) {
        AdvectionSrcForMomWrapper_N<InterpType_H,CENTERED2,UPWINDALL>(bxx, bxy, bxz,
                                                                      rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                                      rho_u, rho_v, rho_w, u, v, w,
                                                                      cellSizeInv, mf_m,
                                                                      mf_u_inv, mf_v_inv,
                                                                      upw_frac_h, upw_frac_v,
                                                                      vert_adv_type,
                                                                      lo_z_face, hi_z_face);
    } else if (vert_adv_type == AdvType::Upwind_3rd) {
        AdvectionSrcForMomWrapper_N<InterpType_H,UPWIND3,UPWINDALL>(bxx, bxy, bxz,
                                                                    rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                                    rho_u, rho_v, rho_w, u, v, w,
                                                                    cellSizeInv, mf_m,
                                                                    mf_u_inv, mf_v_inv,
                                                                    upw_frac_h, upw_frac_v,
                                                                    vert_adv_type,
                                                                    lo_z_face, hi_z_face);
    } else if (vert_adv_type == AdvType::Centered_4th) {
        AdvectionSrcForMomWrapper_N<InterpType_H,CENTERED4,UPWINDALL>(bxx, bxy, bxz,
                                                                      rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                                      rho_u, rho_v, rho_w, u, v, w,
                                                                      cellSizeInv, mf_m,
                                                                      mf_u_inv, mf_v_inv,
                                                                      upw_frac_h, upw_frac_v,
                                                                      vert_adv_type,
                                                                      lo_z_face, hi_z_face);
    } else if (vert_adv_type == AdvType::Upwind_5th) {
        AdvectionSrcForMomWrapper_N<InterpType_H,UPWIND5,UPWINDALL>(bxx, bxy, bxz,
                                                                    rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                                    rho_u, rho_v, rho_w, u, v, w,
                                                                    cellSizeInv, mf_m,
                                                                    mf_u_inv, mf_v_inv,
                                                                    upw_frac_h, upw_frac_v,
                                                                    vert_adv_type,
                                                                    lo_z_face, hi_z_face);
    } else if (vert_adv_type == AdvType::Centered_6th) {
        AdvectionSrcForMomWrapper_N<InterpType_H,CENTERED6,UPWINDALL>(bxx, bxy, bxz,
                                                                      rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                                      rho_u, rho_v, rho_w, u, v, w,
                                                                      cellSizeInv, mf_m,
                                                                      mf_u_inv, mf_v_inv,
                                                                      upw_frac_h, upw_frac_v,
                                                                      vert_adv_type,
                                                                      lo_z_face, hi_z_face);
    } else {
        AMREX_ASSERT_WITH_MESSAGE(false, "Unknown advection scheme!");
    }
}
