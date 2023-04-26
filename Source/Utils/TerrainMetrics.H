#ifndef _TERRAIN_METRIC_H_
#define _TERRAIN_METRIC_H_

#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <IndexDefines.H>

/**
 * Utility routines for constructing terrain metric terms
 */

// Declare function for ERF.cpp
void init_terrain_grid( const amrex::Geometry& geom, amrex::MultiFab& z_phys_nd);

//*****************************************************************************************
// Compute terrain metric terms at cell-center
//*****************************************************************************************
// Metric is at cell center
AMREX_FORCE_INLINE
AMREX_GPU_DEVICE
amrex::Real
Compute_h_zeta_AtCellCenter (const int &i, const int &j, const int &k,
                             const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                             const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dzInv = cellSizeInv[2];
    amrex::Real met_h_zeta = 0.25 * dzInv *
      ( z_nd(i+1,j,k+1) + z_nd(i+1,j+1,k+1) + z_nd(i,j,k+1) + z_nd(i,j+1,k+1)
       -z_nd(i+1,j,k  ) - z_nd(i+1,j+1,k  ) - z_nd(i,j,k  ) - z_nd(i,j+1,k  ) );
    return met_h_zeta;
}

// Metric is at cell center
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_xi_AtCellCenter (const int &i, const int &j, const int &k,
                           const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                           const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dxInv = cellSizeInv[0];
    amrex::Real met_h_xi   = 0.25 * dxInv *
      ( z_nd(i+1,j,k) + z_nd(i+1,j+1,k) + z_nd(i+1,j,k+1) + z_nd(i+1,j+1,k+1)
       -z_nd(i  ,j,k) - z_nd(i  ,j+1,k) - z_nd(i  ,j,k+1) - z_nd(i  ,j+1,k+1) );
    return met_h_xi;
}

// Metric is at cell center
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_eta_AtCellCenter (const int &i, const int &j, const int &k,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                            const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dyInv = cellSizeInv[1];
    amrex::Real met_h_eta  = 0.25 * dyInv *
      ( z_nd(i,j+1,k) + z_nd(i+1,j+1,k) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k+1)
       -z_nd(i,j  ,k) - z_nd(i+1,j  ,k) - z_nd(i,j  ,k+1) - z_nd(i+1,j  ,k+1) );
    return met_h_eta;
}

//*****************************************************************************************
// Compute terrain metric terms at face-centers
//*****************************************************************************************
// Metric coincides with U location
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_zeta_AtIface (const int &i, const int &j, const int &k,
                        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                        const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dzInv = cellSizeInv[2];
    amrex::Real met_h_zeta = 0.5 * dzInv * ( z_nd(i,j,k+1) + z_nd(i,j+1,k+1)
                                           - z_nd(i,j,k  ) - z_nd(i,j+1,k  ) );
    return met_h_zeta;
}

// Metric coincides with U location
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_xi_AtIface (const int &i, const int &j, const int &k,
                      const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                      const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dxInv = cellSizeInv[0];
    amrex::Real met_h_xi   = 0.125 * dxInv *
      ( z_nd(i+1,j,k) + z_nd(i+1,j,k+1) + z_nd(i+1,j+1,k) + z_nd(i+1,j+1,k+1)
       -z_nd(i-1,j,k) - z_nd(i-1,j,k+1) - z_nd(i-1,j+1,k) - z_nd(i-1,j+1,k+1) );
    return met_h_xi;
}

// Metric coincides with U location
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_eta_AtIface (const int &i, const int &j, const int &k,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dyInv = cellSizeInv[1];
    amrex::Real met_h_eta  = 0.5 * dyInv * ( z_nd(i,j+1,k  ) + z_nd(i,j+1,k+1)
                                           - z_nd(i,j  ,k  ) - z_nd(i,j  ,k+1) );
    return met_h_eta;
}

// Metric coincides with V location
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_zeta_AtJface (const int &i, const int &j, const int &k,
                        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                        const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dzInv = cellSizeInv[2];
    amrex::Real met_h_zeta = 0.5 * dzInv * ( z_nd(i,j,k+1) + z_nd(i+1,j,k+1)
                                           - z_nd(i,j,k  ) - z_nd(i+1,j,k  ) );
    return met_h_zeta;
}

// Metric coincides with V location
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_xi_AtJface (const int &i, const int &j, const int &k,
                      const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                      const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dxInv = cellSizeInv[0];
    amrex::Real met_h_xi   = 0.5 * dxInv * ( z_nd(i+1,j,k) + z_nd(i+1,j,k+1)
                                            -z_nd(i  ,j,k) - z_nd(i  ,j,k+1) );
    return met_h_xi;
}

// Metric coincides with V location
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_eta_AtJface (const int &i, const int &j, const int &k,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dyInv = cellSizeInv[1];
    amrex::Real met_h_eta  = 0.125 * dyInv *
      ( z_nd(i,j+1,k) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k) + z_nd(i+1,j+1,k+1)
       -z_nd(i,j-1,k) - z_nd(i,j-1,k+1) - z_nd(i+1,j-1,k) - z_nd(i+1,j-1,k+1) );
    return met_h_eta;
}

// Metric coincides with K location
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_zeta_AtKface (const int &i, const int &j, const int &k,
                        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                        const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dzInv = cellSizeInv[2];
    amrex::Real met_h_zeta = 0.125 * dzInv *
      ( z_nd(i,j,k+1) + z_nd(i+1,j,k+1) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k+1)
       -z_nd(i,j,k-1) - z_nd(i+1,j,k-1) - z_nd(i,j+1,k-1) - z_nd(i+1,j+1,k-1) );
    return met_h_zeta;
}

// Metric coincides with K location
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_xi_AtKface (const int &i, const int &j, const int &k,
                      const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                      const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dxInv = cellSizeInv[0];
    amrex::Real met_h_xi   = 0.5 * dxInv * ( z_nd(i+1,j,k) + z_nd(i+1,j+1,k)
                                            -z_nd(i  ,j,k) - z_nd(i  ,j+1,k) );
    return met_h_xi;
}

// Metric coincides with K location
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_eta_AtKface (const int &i, const int &j, const int &k,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dyInv = cellSizeInv[1];
    amrex::Real met_h_eta  = 0.5 * dyInv * ( z_nd(i,j+1,k) + z_nd(i+1,j+1,k)
                                            -z_nd(i,j  ,k) - z_nd(i+1,j  ,k) );
    return met_h_eta;
}

//*****************************************************************************************
// Compute terrain metric terms at edge-centers
//*****************************************************************************************
// -- EdgeCenterK --

// Metric is at edge and center Z (red pentagon)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_zeta_AtEdgeCenterK (const int &i, const int &j, const int &k,
                              const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                              const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dzInv = cellSizeInv[2];
    amrex::Real met_h_zeta = dzInv * (z_nd(i,j,k+1) - z_nd(i,j,k));
    return met_h_zeta;
}

// Metric is at edge and center Z (red pentagon)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_xi_AtEdgeCenterK (const int &i, const int &j, const int &k,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                            const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dxInv = cellSizeInv[0];
    amrex::Real met_h_xi  = 0.25 * dxInv *
      ( z_nd(i+1,j,k) + z_nd(i+1,j,k+1)
       -z_nd(i-1,j,k) - z_nd(i-1,j,k+1) );
    return met_h_xi;
}

// Metric is at edge and center Z (red pentagon)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_eta_AtEdgeCenterK (const int &i, const int &j, const int &k,
                             const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                             const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dyInv = cellSizeInv[1];
    amrex::Real met_h_eta = 0.25 * dyInv *
      ( z_nd(i,j+1,k) + z_nd(i,j+1,k+1)
       -z_nd(i,j-1,k) - z_nd(i,j-1,k+1) );
    return met_h_eta;
}

// -- EdgeCenterJ --

// Metric is at edge and center Y (magenta cross)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_zeta_AtEdgeCenterJ (const int &i, const int &j, const int &k,
                              const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                              const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dzInv = cellSizeInv[2];
    amrex::Real met_h_zeta = 0.25 * dzInv * ( z_nd(i,j,k+1) + z_nd(i,j+1,k+1)
                                             -z_nd(i,j,k-1) - z_nd(i,j+1,k-1) );
    return met_h_zeta;
}

// Metric is at edge and center Y (magenta cross)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_xi_AtEdgeCenterJ (const int &i, const int &j, const int &k,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                            const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dxInv = cellSizeInv[0];
    amrex::Real met_h_xi = 0.25 * dxInv *
      ( z_nd(i+1,j+1,k) + z_nd(i+1,j  ,k)
       -z_nd(i-1,j+1,k) - z_nd(i-1,j  ,k) );
    return met_h_xi;
}

// Metric is at edge and center Y (magenta cross)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_eta_AtEdgeCenterJ (const int &i, const int &j, const int &k,
                             const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                             const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dyInv = cellSizeInv[1];
    amrex::Real met_h_eta = dyInv * ( z_nd(i,j+1,k) - z_nd(i,j,k) );
    return met_h_eta;
}

// -- EdgeCenterI --

// Metric is at edge and center Y (magenta cross)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_zeta_AtEdgeCenterI (const int &i, const int &j, const int &k,
                              const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                              const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dzInv = cellSizeInv[2];
    amrex::Real met_h_zeta = 0.25 * dzInv * ( z_nd(i,j,k+1) + z_nd(i+1,j,k+1)
                                             -z_nd(i,j,k-1) - z_nd(i+1,j,k-1) );
    return met_h_zeta;
}

// Metric is at edge and center Y (magenta cross)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_xi_AtEdgeCenterI (const int &i, const int &j, const int &k,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                            const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dxInv = cellSizeInv[0];
    amrex::Real met_h_xi  = dxInv * ( z_nd(i+1,j,k) - z_nd(i,j,k) );
    return met_h_xi;
}

// Metric is at edge and center Y (magenta cross)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_eta_AtEdgeCenterI (const int &i, const int &j, const int &k,
                             const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                             const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dyInv = cellSizeInv[1];
    amrex::Real met_h_eta = 0.25 * dyInv *
      ( z_nd(i+1,j+1,k) + z_nd(i,j+1,k)
       -z_nd(i+1,j-1,k) - z_nd(i,j-1,k) );
    return met_h_eta;
}

//*****************************************************************************************
// Map between W <--> Omega
//*****************************************************************************************
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
OmegaFromW (int i, int j, int k, amrex::Real w,
            const amrex::Array4<const amrex::Real> u,
            const amrex::Array4<const amrex::Real> v,
            const amrex::Array4<const amrex::Real> z_nd,
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxInv)
{
  // This is dh/dxi at z-face (i,j,k-1/2)
  amrex::Real met_zlo_xi   = 0.5 * dxInv[0] *
    ( z_nd(i+1,j+1,k  ) + z_nd(i+1,j  ,k  )    // hi i, lo k
     -z_nd(i  ,j+1,k  ) - z_nd(i  ,j  ,k  ) ); // lo i, lo k

  // This is dh/deta at z-face (i,j,k-1/2)
  amrex::Real met_zlo_eta  = 0.5 * dxInv[1] *
    ( z_nd(i+1,j+1,k  ) + z_nd(i  ,j+1,k  )    // hi j, lo k
     -z_nd(i+1,j  ,k  ) - z_nd(i  ,j  ,k  ) ); // lo j, lo k

  // No-slip BC
  /*
  amrex::Real omega;
  if (k == 0) {
      omega = 0.;
  } else {
      amrex::Real vec_zlo_xi   = 0.25 * ( u(i,j,k-1) + u(i+1,j,k-1) + u(i,j,k) + u(i+1,j,k) );
      amrex::Real vec_zlo_eta  = 0.25 * ( v(i,j,k-1) + v(i,j+1,k-1) + v(i,j,k) + v(i,j+1,k) );
      omega =  -met_zlo_xi * vec_zlo_xi - met_zlo_eta * vec_zlo_eta + w;
  }
  */

  // Slip BC or moving terrain
  // Use extrapolation instead of interpolation if at the bottom boundary
  amrex::Real vec_zlo_xi   = (k == 0) ? 1.5 * (0.5*(u(i,j,k)+u(i+1,j,k))) - 0.5*(0.5*(u(i,j,k+1)+u(i+1,j,k+1))) :
    0.25 * ( u(i,j,k-1) + u(i+1,j,k-1) + u(i,j,k) + u(i+1,j,k) );
  amrex::Real vec_zlo_eta  = (k == 0) ? 1.5 * (0.5*(v(i,j,k)+v(i,j+1,k))) - 0.5*(0.5*(v(i,j,k+1)+v(i,j+1,k+1))) :
    0.25 * ( v(i,j,k-1) + v(i,j+1,k-1) + v(i,j,k) + v(i,j+1,k) );

  amrex::Real omega =  -met_zlo_xi * vec_zlo_xi - met_zlo_eta * vec_zlo_eta + w;
  return omega;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
WFromOmega (int i, int j, int k, amrex::Real omega,
            const amrex::Array4<const amrex::Real>& z_nd,
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxInv)
{
  // This is dh/dxi at z-face (i,j,k-1/2)
  amrex::Real met_xi   = 0.5 * dxInv[0] *
    ( z_nd(i+1,j+1,k  ) + z_nd(i+1,j  ,k  )    // hi i, lo k
     -z_nd(i  ,j+1,k  ) - z_nd(i  ,j  ,k  ) ); // lo i, lo k

  // This is dh/deta at z-face (i,j,k-1/2)
  amrex::Real met_eta  = 0.5 * dxInv[1] *
    ( z_nd(i+1,j+1,k  ) + z_nd(i  ,j+1,k  )    // hi j, lo k
     -z_nd(i+1,j  ,k  ) - z_nd(i  ,j  ,k  ) ); // lo j, lo k

  // This dots the Omega vector (normal to the face) with the z-direction
  amrex::Real w = omega / std::sqrt(1.0 + met_xi*met_xi + met_eta*met_eta);

  return w;
}
#endif
