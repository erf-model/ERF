#include <IndexDefines.H>
#include <Interpolation.H>

/**
 * Wrapper function for computing the advective tendency w/ spatial order > 2.
 */
template<typename InterpType>
void
AdvectionSrcForRhoThetaWrapper_N(const amrex::Box& bx,
                                 const amrex::Dim3& vbx_hi,
                                 const amrex::Real fac,
                                 InterpType interp_prim,
                                 const amrex::Array4<      amrex::Real>& advectionSrc,
                                 const amrex::Array4<const amrex::Real>& rho_u,
                                 const amrex::Array4<const amrex::Real>& rho_v,
                                 const amrex::Array4<const amrex::Real>& rho_w,
                                 const amrex::Array4<      amrex::Real>& avg_xmom,
                                 const amrex::Array4<      amrex::Real>& avg_ymom,
                                 const amrex::Array4<      amrex::Real>& avg_zmom,
                                 const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                                 const amrex::Array4<const amrex::Real>& mf_m,
                                 const amrex::Array4<const amrex::Real>& mf_u,
                                 const amrex::Array4<const amrex::Real>& mf_v)
{
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      amrex::Real xflux_lo = rho_u(i  ,j,k) / mf_u(i  ,j  ,0);
      amrex::Real xflux_hi = rho_u(i+1,j,k) / mf_u(i+1,j  ,0);
      amrex::Real yflux_lo = rho_v(i,j  ,k) / mf_v(i  ,j  ,0);
      amrex::Real yflux_hi = rho_v(i,j+1,k) / mf_v(i  ,j+1,0);
      amrex::Real zflux_lo = rho_w(i,j,k  );
      amrex::Real zflux_hi = rho_w(i,j,k+1);
      
      avg_xmom(i  ,j,k) += fac*xflux_lo;
      if (i == vbx_hi.x) avg_xmom(i+1,j,k) += fac*xflux_hi;
        
      avg_ymom(i,j  ,k) += fac*yflux_lo;
      if (j == vbx_hi.y) avg_ymom(i,j+1,k) += fac*yflux_hi;
        
      avg_zmom(i,j,k  ) += fac*zflux_lo;
      if (k == vbx_hi.z) avg_zmom(i,j,k+1) += fac*zflux_hi;
      
      amrex::Real mf   = mf_m(i,j,0);
      amrex::Real mfsq = mf*mf;
      
      advectionSrc(i,j,k,0) = -(
                                ( xflux_hi - xflux_lo ) * dxInv * mfsq +
                                ( yflux_hi - yflux_lo ) * dyInv * mfsq +
                                ( zflux_hi - zflux_lo ) * dzInv);

      const int prim_index = 0;
      amrex::Real interpx_hi(0.), interpx_lo(0.);
      amrex::Real interpy_hi(0.), interpy_lo(0.);
      amrex::Real interpz_hi(0.), interpz_lo(0.);
      interp_prim.InterpolateInX(i,j,k,0,interpx_hi,interpx_lo,1.0,1.0);
      //interp_prim.InterpolateInX(i,j,k,0,interpy_hi,interpy_lo,rho_u(i  ,j+1,k  ),rho_v(i  ,j  ,k  ));
      //interp_prim.InterpolateInX(i,j,k,0,interpz_hi,interpz_lo,rho_w(i  ,j  ,k+1),rho_w(i  ,j  ,k  ));
      advectionSrc(i,j,k,1) = -(
                                ( xflux_hi * interpx_hi - xflux_lo * interpx_lo ) * dxInv * mfsq +
                                ( yflux_hi * interpy_hi - yflux_lo * interpy_lo ) * dyInv * mfsq +
                                ( zflux_hi * interpz_hi - zflux_lo * interpz_lo ) * dzInv);
    });
}
