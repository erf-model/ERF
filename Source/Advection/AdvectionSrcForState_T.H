#include <IndexDefines.H>
#include <TerrainMetrics.H>
#include <Interpolation.H>

/**
 * Wrapper function for computing the advective tendency w/ spatial order > 2.
 */
template<typename InterpType_H, typename InterpType_V>
void
AdvectionSrcForRhoThetaWrapper_T(const amrex::Box& bx,
                                 const amrex::Dim3& vbx_hi,
                                 const amrex::Real fac,
                                 const amrex::Array4<      amrex::Real>& advectionSrc,
                                 const amrex::Array4<const amrex::Real>& cell_prim,
                                 const amrex::Array4<const amrex::Real>& rho_u,
                                 const amrex::Array4<const amrex::Real>& rho_v,
                                 const amrex::Array4<const amrex::Real>& Omega,
                                 const amrex::Array4<      amrex::Real>& avg_xmom,
                                 const amrex::Array4<      amrex::Real>& avg_ymom,
                                 const amrex::Array4<      amrex::Real>& avg_zmom,
                                 const amrex::Array4<const amrex::Real>& z_nd,
                                 const amrex::Array4<const amrex::Real>& detJ,
                                 const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                                 const amrex::Array4<const amrex::Real>& mf_m,
                                 const amrex::Array4<const amrex::Real>& mf_u,
                                 const amrex::Array4<const amrex::Real>& mf_v)
{
    // Instantiate struct
    InterpType_H interp_prim_h(cell_prim);
    InterpType_V interp_prim_v(cell_prim);

    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      amrex::Real invdetJ = 1./ detJ(i,j,k);

      amrex::Real xflux_lo = rho_u(i  ,j,k) / mf_u(i  ,j  ,0);
      amrex::Real xflux_hi = rho_u(i+1,j,k) / mf_u(i+1,j  ,0);
      amrex::Real yflux_lo = rho_v(i,j  ,k) / mf_v(i  ,j  ,0);
      amrex::Real yflux_hi = rho_v(i,j+1,k) / mf_v(i  ,j+1,0);
      amrex::Real zflux_lo = Omega(i,j,k  );
      amrex::Real zflux_hi = Omega(i,j,k+1);

      amrex::Real met_h_zeta_xlo = Compute_h_zeta_AtIface(i  ,j  ,k,cellSizeInv,z_nd);
      xflux_lo *= met_h_zeta_xlo;
      amrex::Real met_h_zeta_xhi = Compute_h_zeta_AtIface(i+1,j  ,k,cellSizeInv,z_nd);
      xflux_hi *= met_h_zeta_xhi;

      amrex::Real met_h_zeta_ylo = Compute_h_zeta_AtJface(i  ,j  ,k,cellSizeInv,z_nd);
      yflux_lo *= met_h_zeta_ylo;
      amrex::Real met_h_zeta_yhi = Compute_h_zeta_AtJface(i  ,j+1,k,cellSizeInv,z_nd);
      yflux_hi *= met_h_zeta_yhi;

      avg_xmom(i  ,j,k) += fac*xflux_lo;
      if (i == vbx_hi.x)
        avg_xmom(i+1,j,k) += fac*xflux_hi;
      avg_ymom(i,j  ,k) += fac*yflux_lo;
      if (j == vbx_hi.y)
        avg_ymom(i,j+1,k) += fac*yflux_hi;
      avg_zmom(i,j,k  ) += fac*zflux_lo;
      if (k == vbx_hi.z)
        avg_zmom(i,j,k+1) += fac*zflux_hi;

      amrex::Real mfsq = mf_m(i,j,0) * mf_m(i,j,0);

      advectionSrc(i,j,k,0) = - invdetJ * (
                                           ( xflux_hi - xflux_lo ) * dxInv * mfsq +
                                           ( yflux_hi - yflux_lo ) * dyInv * mfsq +
                                           ( zflux_hi - zflux_lo ) * dzInv );

      const int prim_index = 0;
      amrex::Real interpx_hi(0.), interpx_lo(0.);
      amrex::Real interpy_hi(0.), interpy_lo(0.);
      amrex::Real interpz_hi(0.), interpz_lo(0.);

      interp_prim_h.InterpolateInX(i+1,j,k,prim_index,interpx_hi,rho_u(i+1,j  ,k  ));
      interp_prim_h.InterpolateInX(i  ,j,k,prim_index,interpx_lo,rho_u(i  ,j  ,k  ));

      interp_prim_h.InterpolateInY(i,j+1,k,prim_index,interpy_hi,rho_v(i  ,j+1,k  ));
      interp_prim_h.InterpolateInY(i,j  ,k,prim_index,interpy_lo,rho_v(i  ,j  ,k  ));

      interp_prim_v.InterpolateInZ(i,j,k+1,prim_index,interpz_hi,Omega(i  ,j  ,k+1));
      interp_prim_v.InterpolateInZ(i,j,k  ,prim_index,interpz_lo,Omega(i  ,j  ,k  ));

      advectionSrc(i,j,k,1) = - invdetJ * (
                                           ( xflux_hi * interpx_hi - xflux_lo * interpx_lo ) * dxInv * mfsq +
                                           ( yflux_hi * interpy_hi - yflux_lo * interpy_lo ) * dyInv * mfsq +
                                           ( zflux_hi * interpz_hi - zflux_lo * interpz_lo ) * dzInv);
    });
}

/**
 * Wrapper function for computing the advective tendency w/ spatial order > 2.
 */
template<typename InterpType_H>
void
AdvectionSrcForRhoThetaVert_T(const amrex::Box& bx,
                              const amrex::Dim3& vbx_hi,
                              const amrex::Real fac,
                              const amrex::Array4<      amrex::Real>& advectionSrc,
                              const amrex::Array4<const amrex::Real>& cell_prim,
                              const amrex::Array4<const amrex::Real>& rho_u,
                              const amrex::Array4<const amrex::Real>& rho_v,
                              const amrex::Array4<const amrex::Real>& Omega,
                              const amrex::Array4<      amrex::Real>& avg_xmom,
                              const amrex::Array4<      amrex::Real>& avg_ymom,
                              const amrex::Array4<      amrex::Real>& avg_zmom,
                              const amrex::Array4<const amrex::Real>& z_nd,
                              const amrex::Array4<const amrex::Real>& detJ,
                              const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                              const amrex::Array4<const amrex::Real>& mf_m,
                              const amrex::Array4<const amrex::Real>& mf_u,
                              const amrex::Array4<const amrex::Real>& mf_v,
                              const AdvType vert_adv_type)
{
    if (vert_adv_type == AdvType::Centered_2nd) {
        AdvectionSrcForRhoThetaWrapper_T<InterpType_H,CENTERED2>(bx, vbx_hi, fac, advectionSrc,
                                                                 cell_prim, rho_u, rho_v, Omega,
                                                                 avg_xmom, avg_ymom, avg_zmom,
                                                                 z_nd, detJ, cellSizeInv, mf_m,
                                                                 mf_u, mf_v);
    } else if (vert_adv_type == AdvType::Upwind_3rd) {
        AdvectionSrcForRhoThetaWrapper_T<InterpType_H,UPWIND3>(bx, vbx_hi, fac, advectionSrc,
                                                               cell_prim, rho_u, rho_v, Omega,
                                                               avg_xmom, avg_ymom, avg_zmom,
                                                               z_nd, detJ, cellSizeInv, mf_m,
                                                               mf_u, mf_v);
    } else if (vert_adv_type == AdvType::Centered_4th) {
        AdvectionSrcForRhoThetaWrapper_T<InterpType_H,CENTERED4>(bx, vbx_hi, fac, advectionSrc,
                                                                 cell_prim, rho_u, rho_v, Omega,
                                                                 avg_xmom, avg_ymom, avg_zmom,
                                                                 z_nd, detJ, cellSizeInv, mf_m,
                                                                 mf_u, mf_v);
    } else if (vert_adv_type == AdvType::Upwind_5th) {
        AdvectionSrcForRhoThetaWrapper_T<InterpType_H,UPWIND5>(bx, vbx_hi, fac, advectionSrc,
                                                               cell_prim, rho_u, rho_v, Omega,
                                                               avg_xmom, avg_ymom, avg_zmom,
                                                               z_nd, detJ, cellSizeInv, mf_m,
                                                               mf_u, mf_v);
    } else if (vert_adv_type == AdvType::Centered_6th) {
        AdvectionSrcForRhoThetaWrapper_T<InterpType_H,CENTERED6>(bx, vbx_hi, fac, advectionSrc,
                                                                 cell_prim, rho_u, rho_v, Omega,
                                                                 avg_xmom, avg_ymom, avg_zmom,
                                                                 z_nd, detJ, cellSizeInv, mf_m,
                                                                 mf_u, mf_v);
    } else {
        AMREX_ASSERT_WITH_MESSAGE(false, "Unknown advection scheme!");
    }
}
