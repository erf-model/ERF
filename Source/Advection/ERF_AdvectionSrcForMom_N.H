#include "ERF_IndexDefines.H"
#include "ERF_Interpolation.H"

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
template<typename InterpType_H>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
AdvectionFluxXForXMom_N (int i, int j, int k,
                         const amrex::Array4<const amrex::Real>& rho_u,
                         const amrex::Array4<      amrex::Real>& xm_x_f_arr,
                         InterpType_H interp_u_h,
                         const amrex::Array4<const amrex::Real>& mf_u_inv)
{
    amrex::Real interp(0.), rho_u_avg(0.);
    rho_u_avg = 0.5 * (rho_u(i-1, j, k) * mf_u_inv(i-1,j,0) + rho_u(i, j, k) * mf_u_inv(i,j,0));
    interp_u_h.InterpolateInX(i,j,k,0,interp,rho_u_avg);
    xm_x_f_arr(i,j,k) = rho_u_avg * interp;
}
template<typename InterpType_H>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
AdvectionFluxYForXMom_N (int i, int j, int k,
                         const amrex::Array4<const amrex::Real>& rho_v,
                         const amrex::Array4<      amrex::Real>& xm_y_f_arr,
                         InterpType_H interp_u_h,
                         const amrex::Array4<const amrex::Real>& mf_v_inv)
{
    amrex::Real interp(0.), rho_v_avg(0.);
    rho_v_avg = 0.5 * (rho_v(i, j  , k) * mf_v_inv(i,j  ,0) + rho_v(i-1, j  , k) * mf_v_inv(i-1,j  ,0));
    interp_u_h.InterpolateInY(i,j,k,0,interp,rho_v_avg);
    xm_y_f_arr(i,j,k) = rho_v_avg * interp;
}
template<typename InterpType_V>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
AdvectionFluxZForXMom_N (int i, int j, int k,
                         const amrex::Array4<const amrex::Real>& rho_w,
                         const amrex::Array4<      amrex::Real>& xm_z_f_arr,
                         InterpType_V interp_u_v)
{
    amrex::Real interp(0.), rho_w_avg(0.);
    rho_w_avg = 0.5 * (rho_w(i, j, k  ) + rho_w(i-1, j, k  ));
    interp_u_v.InterpolateInZ(i,j,k,0,interp,rho_w_avg);
    xm_z_f_arr(i,j,k) = rho_w_avg * interp;
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
template<typename InterpType_H>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
AdvectionFluxXForYMom_N (int i, int j, int k,
                         const amrex::Array4<const amrex::Real>& rho_u,
                         const amrex::Array4<      amrex::Real>& ym_x_f_arr,
                         InterpType_H interp_v_h,
                         const amrex::Array4<const amrex::Real>& mf_u_inv)
{
    amrex::Real interp(0.), rho_u_avg(0.);
    rho_u_avg = 0.5 * (rho_u(i  , j, k) * mf_u_inv(i  ,j  ,0) + rho_u(i  , j-1, k) * mf_u_inv(i  ,j-1,0));
    interp_v_h.InterpolateInX(i,j,k,0,interp,rho_u_avg);
    ym_x_f_arr(i,j,k) = rho_u_avg * interp;
}
template<typename InterpType_H>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
AdvectionFluxYForYMom_N (int i, int j, int k,
                         const amrex::Array4<const amrex::Real>& rho_v,
                         const amrex::Array4<      amrex::Real>& ym_y_f_arr,
                         InterpType_H interp_v_h,
                         const amrex::Array4<const amrex::Real>& mf_v_inv)
{
    amrex::Real interp(0.), rho_v_avg(0.);
    rho_v_avg = 0.5 * (rho_v(i, j, k) * mf_v_inv(i  ,j  ,0) + rho_v(i, j-1, k) * mf_v_inv(i  ,j-1,0));
    interp_v_h.InterpolateInY(i,j,k,0,interp,rho_v_avg);
    ym_y_f_arr(i,j,k) = rho_v_avg * interp;
}
template<typename InterpType_V>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
AdvectionFluxZForYMom_N (int i, int j, int k,
                         const amrex::Array4<const amrex::Real>& rho_w,
                         const amrex::Array4<      amrex::Real>& ym_z_f_arr,
                         InterpType_V interp_v_v)
{
    amrex::Real interp(0.), rho_w_avg(0.);
    rho_w_avg = 0.5 * (rho_w(i, j, k  ) + rho_w(i, j-1, k  ));
    interp_v_v.InterpolateInZ(i,j,k  ,0,interp,rho_w_avg);
    ym_z_f_arr(i,j,k) = rho_w_avg * interp;
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
template<typename InterpType_H>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
AdvectionFluxXForZMom_N (int i, int j, int k,
                         const amrex::Array4<const amrex::Real>& rho_u,
                         const amrex::Array4<      amrex::Real>& zm_x_f_arr,
                         InterpType_H   interp_w_h,
                         const amrex::Array4<const amrex::Real>& mf_u_inv)
{
    amrex::Real interp(0.), rho_u_avg(0.);
    rho_u_avg = 0.5 * (rho_u(i  , j, k) + rho_u(i  , j, k-1)) * mf_u_inv(i  ,j  ,0);
    interp_w_h.InterpolateInX(i,j,k,0,interp,rho_u_avg);
    zm_x_f_arr(i,j,k) = rho_u_avg * interp;
}
template<typename InterpType_H>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
AdvectionFluxYForZMom_N (int i, int j, int k,
                         const amrex::Array4<const amrex::Real>& rho_v,
                         const amrex::Array4<      amrex::Real>& zm_y_f_arr,
                         InterpType_H   interp_w_h,
                         const amrex::Array4<const amrex::Real>& mf_v_inv)
{
    amrex::Real interp(0.), rho_v_avg(0.);
    rho_v_avg = 0.5 * (rho_v(i, j  , k) + rho_v(i, j  , k-1)) * mf_v_inv(i  ,j  ,0);
    interp_w_h.InterpolateInY(i,j,k,0,interp,rho_v_avg);
    zm_y_f_arr(i,j,k) = rho_v_avg * interp;
}
template<typename InterpType_V, typename WallInterpType>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
AdvectionFluxZForZMom_N (int i, int j, int k,
                         const amrex::Array4<const amrex::Real>& rho_w,
                         const amrex::Array4<      amrex::Real>& zm_z_f_arr,
                         const amrex::Array4<const amrex::Real>& w,
                         InterpType_V   interp_w_v,
                         WallInterpType interp_w_wall,
                         const AdvType vert_adv_type,
                         const int lo_z_face, const int hi_z_face)
{
    amrex::Real interp(0.), rho_w_avg(0.);
    if ( (k==lo_z_face) || (k==hi_z_face) ) {
        zm_z_f_arr(i,j,k) = rho_w(i,j,k) * w(i,j,k);
    } else {
        rho_w_avg = 0.5 * (rho_w(i,j,k) + rho_w(i,j,k-1));
        if ( (k==lo_z_face+1) || (k==hi_z_face-1) ) {
            interp_w_wall.InterpolateInZ(i,j,k,0,interp,rho_w_avg,AdvType::Centered_2nd);
        } else if ( (k==lo_z_face+2 || k==hi_z_face-1) ||
                    (k==hi_z_face-2 || k==lo_z_face+1) ) {
            if (vert_adv_type != AdvType::Centered_2nd && vert_adv_type != AdvType::Upwind_3rd) {
                interp_w_wall.InterpolateInZ(i,j,k,0,interp,rho_w_avg,AdvType::Centered_4th);
            } else {
                interp_w_wall.InterpolateInZ(i,j,k,0,interp,rho_w_avg,vert_adv_type);
            }
        } else {
            interp_w_v.InterpolateInZ(i,j,k,0,interp,rho_w_avg);
        }
        zm_z_f_arr(i,j,k) = rho_w_avg * interp;
    }
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
                             const amrex::Gpu::DeviceVector<amrex::Real>& stretched_dz_d,
                             const amrex::Array4<const amrex::Real>& mf_m,
                             const amrex::Array4<const amrex::Real>& mf_u_inv,
                             const amrex::Array4<const amrex::Real>& mf_v_inv,
                             const amrex::Real upw_frac_h,
                             const amrex::Real upw_frac_v,
                             const AdvType vert_adv_type,
                             const int lo_z_face, const int hi_z_face)
{
    // Instantiate the appropriate structs
    InterpType_H interp_u_h(u, upw_frac_h); InterpType_V interp_u_v(u, upw_frac_v); // X-MOM
    InterpType_H interp_v_h(v, upw_frac_h); InterpType_V interp_v_v(v, upw_frac_v); // Y-MOM
    InterpType_H interp_w_h(w, upw_frac_h); InterpType_V interp_w_v(w, upw_frac_v); // Z-MOM
    WallInterpType interp_w_wall(w, upw_frac_v); // Z-MOM @ wall

    auto dxInv  = cellSizeInv[0];
    auto dyInv  = cellSizeInv[1];
    auto dz_ptr = stretched_dz_d.data();

    // X momentum fluxes and advection source
    amrex::FArrayBox flux_x(amrex::grow(amrex::convert(bxx,amrex::IntVect(0)),0,1),1,amrex::The_Async_Arena());
    amrex::FArrayBox flux_y(amrex::surroundingNodes(bxx,1),1,amrex::The_Async_Arena());
    amrex::FArrayBox flux_z(amrex::surroundingNodes(bxx,2),1,amrex::The_Async_Arena());
    {
    amrex::Array4<amrex::Real> x_f_arr = flux_x.array();
    amrex::Array4<amrex::Real> y_f_arr = flux_y.array();
    amrex::Array4<amrex::Real> z_f_arr = flux_z.array();
    amrex::ParallelFor(amrex::grow(amrex::convert(bxx,amrex::IntVect(0)),0,1),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        AdvectionFluxXForXMom_N(i, j, k, rho_u, x_f_arr, interp_u_h, mf_u_inv);
    });
    amrex::ParallelFor(amrex::surroundingNodes(bxx,1),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        AdvectionFluxYForXMom_N(i, j, k, rho_v, y_f_arr, interp_u_h, mf_v_inv);
    });
    amrex::ParallelFor(amrex::surroundingNodes(bxx,2),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        AdvectionFluxZForXMom_N(i, j, k, rho_w, z_f_arr, interp_u_v);
    });
    amrex::ParallelFor(bxx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        auto dzInv = 1.0/dz_ptr[k];
        amrex::Real mfsq = 1 / (mf_u_inv(i,j,0) * mf_u_inv(i,j,0));
        rho_u_rhs(i, j, k) = -( (x_f_arr(i  ,j  ,k  ) - x_f_arr(i-1,j  ,k  )) * dxInv * mfsq
                              + (y_f_arr(i  ,j+1,k  ) - y_f_arr(i  ,j  ,k  )) * dyInv * mfsq
                              + (z_f_arr(i  ,j  ,k+1) - z_f_arr(i  ,j  ,k  )) * dzInv );
    });
    }


    // Y momentum fluxes and advection source
    flux_x.resize(amrex::surroundingNodes(bxy,0));
    flux_y.resize(amrex::grow(amrex::convert(bxy,amrex::IntVect(0)),1,1));
    flux_z.resize(amrex::surroundingNodes(bxy,2));
    {
    amrex::Array4<amrex::Real> x_f_arr = flux_x.array();
    amrex::Array4<amrex::Real> y_f_arr = flux_y.array();
    amrex::Array4<amrex::Real> z_f_arr = flux_z.array();
    amrex::ParallelFor(amrex::surroundingNodes(bxy,0),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        AdvectionFluxXForYMom_N(i, j, k, rho_u, x_f_arr, interp_v_h, mf_u_inv);
    });
    amrex::ParallelFor(amrex::grow(amrex::convert(bxy,amrex::IntVect(0)),1,1),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        AdvectionFluxYForYMom_N(i, j, k, rho_v, y_f_arr, interp_v_h, mf_v_inv);
    });
    amrex::ParallelFor(amrex::surroundingNodes(bxy,2),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        AdvectionFluxZForYMom_N(i, j, k, rho_w, z_f_arr, interp_v_v);
    });
    amrex::ParallelFor(bxy,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        auto dzInv = 1.0/dz_ptr[k];
        amrex::Real mfsq = 1 / (mf_v_inv(i,j,0) * mf_v_inv(i,j,0));
        rho_v_rhs(i, j, k) = -( (x_f_arr(i+1,j  ,k  ) - x_f_arr(i  ,j  ,k  )) * dxInv * mfsq
                              + (y_f_arr(i  ,j  ,k  ) - y_f_arr(i  ,j-1,k  )) * dyInv * mfsq
                              + (z_f_arr(i  ,j  ,k+1) - z_f_arr(i  ,j  ,k  )) * dzInv );
    });
    }


    // Z momentum fluxes and advection source
    flux_x.resize(amrex::surroundingNodes(bxz,0));
    flux_z.resize(amrex::surroundingNodes(bxz,1));
    flux_y.resize(amrex::grow(amrex::convert(bxz,amrex::IntVect(0)),2,1));
    amrex::Print() << "BX: " << bxz << ' '
                   << amrex::grow(amrex::convert(bxz,amrex::IntVect(0)),2,1) << ' '
                   << lo_z_face << ' ' << hi_z_face << "\n";
    {
    amrex::Array4<amrex::Real> x_f_arr = flux_x.array();
    amrex::Array4<amrex::Real> y_f_arr = flux_y.array();
    amrex::Array4<amrex::Real> z_f_arr = flux_z.array();
    amrex::ParallelFor(amrex::surroundingNodes(bxz,0),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        AdvectionFluxXForZMom_N(i, j, k, rho_u, x_f_arr, interp_w_h, mf_u_inv);
    });
    amrex::ParallelFor(amrex::surroundingNodes(bxz,1),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        AdvectionFluxYForZMom_N(i, j, k, rho_v, y_f_arr, interp_w_h, mf_v_inv);
    });
    amrex::ParallelFor(amrex::grow(amrex::convert(bxz,amrex::IntVect(0)),2,1),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        AdvectionFluxZForZMom_N(i, j, k, rho_w, z_f_arr, w, interp_w_v, interp_w_wall,
                                vert_adv_type, lo_z_face, hi_z_face);
    });
    amrex::ParallelFor(bxz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        auto dzInv = 2.0/(dz_ptr[k] + dz_ptr[k-1]);
        amrex::Real mfsq = mf_m(i,j,0) * mf_m(i,j,0);
        rho_w_rhs(i, j, k) = -( (x_f_arr(i+1,j  ,k  ) - x_f_arr(i  ,j  ,k  )) * dxInv * mfsq
                              + (y_f_arr(i  ,j+1,k  ) - y_f_arr(i  ,j  ,k  )) * dyInv * mfsq
                              + (z_f_arr(i  ,j  ,k  ) - z_f_arr(i  ,j  ,k-1)) * dzInv );
    });
    }
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
                          const amrex::Gpu::DeviceVector<amrex::Real>& stretched_dz_d,
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
                                                                      cellSizeInv, stretched_dz_d, mf_m,
                                                                      mf_u_inv, mf_v_inv,
                                                                      upw_frac_h, upw_frac_v,
                                                                      vert_adv_type,
                                                                      lo_z_face, hi_z_face);
    } else if (vert_adv_type == AdvType::Upwind_3rd) {
        AdvectionSrcForMomWrapper_N<InterpType_H,UPWIND3,UPWINDALL>(bxx, bxy, bxz,
                                                                    rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                                    rho_u, rho_v, rho_w, u, v, w,
                                                                    cellSizeInv, stretched_dz_d, mf_m,
                                                                    mf_u_inv, mf_v_inv,
                                                                    upw_frac_h, upw_frac_v,
                                                                    vert_adv_type,
                                                                    lo_z_face, hi_z_face);
    } else if (vert_adv_type == AdvType::Centered_4th) {
        AdvectionSrcForMomWrapper_N<InterpType_H,CENTERED4,UPWINDALL>(bxx, bxy, bxz,
                                                                      rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                                      rho_u, rho_v, rho_w, u, v, w,
                                                                      cellSizeInv, stretched_dz_d, mf_m,
                                                                      mf_u_inv, mf_v_inv,
                                                                      upw_frac_h, upw_frac_v,
                                                                      vert_adv_type,
                                                                      lo_z_face, hi_z_face);
    } else if (vert_adv_type == AdvType::Upwind_5th) {
        AdvectionSrcForMomWrapper_N<InterpType_H,UPWIND5,UPWINDALL>(bxx, bxy, bxz,
                                                                    rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                                    rho_u, rho_v, rho_w, u, v, w,
                                                                    cellSizeInv, stretched_dz_d, mf_m,
                                                                    mf_u_inv, mf_v_inv,
                                                                    upw_frac_h, upw_frac_v,
                                                                    vert_adv_type,
                                                                    lo_z_face, hi_z_face);
    } else if (vert_adv_type == AdvType::Centered_6th) {
        AdvectionSrcForMomWrapper_N<InterpType_H,CENTERED6,UPWINDALL>(bxx, bxy, bxz,
                                                                      rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                                                      rho_u, rho_v, rho_w, u, v, w,
                                                                      cellSizeInv, stretched_dz_d, mf_m,
                                                                      mf_u_inv, mf_v_inv,
                                                                      upw_frac_h, upw_frac_v,
                                                                      vert_adv_type,
                                                                      lo_z_face, hi_z_face);
    } else {
        AMREX_ASSERT_WITH_MESSAGE(false, "Unknown advection scheme!");
    }
}
