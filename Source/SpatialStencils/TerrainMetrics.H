#ifndef _TERRAIN_METRIC_H_
#define _TERRAIN_METRIC_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include "IndexDefines.H"

//*****************************************************************************************
// Compute terrain metric terms at cell-center
//*****************************************************************************************
// Metric is at cell center
AMREX_GPU_DEVICE
inline void
ComputeMetricAtCellCenter(const int &i, const int &j, const int &k,
                          amrex::Real& met_h_xi,
                          amrex::Real& met_h_eta,
                          amrex::Real& met_h_zeta,
                          const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                          const amrex::Array4<const amrex::Real>& z_nd, const int flag)
{
  met_h_xi   = 0.0;
  met_h_eta  = 0.0;
  met_h_zeta = 0.0;
  amrex::Real dxInv = cellSizeInv[0];
  amrex::Real dyInv = cellSizeInv[1];
  amrex::Real dzInv = cellSizeInv[2];

  switch(flag) {
  case TerrainMet::h_xi:
    met_h_xi   = 0.25 * dxInv *
      ( z_nd(i+1,j,k) + z_nd(i+1,j+1,k) + z_nd(i+1,j,k+1) + z_nd(i+1,j+1,k+1)
       -z_nd(i  ,j,k) - z_nd(i  ,j+1,k) - z_nd(i  ,j,k+1) - z_nd(i  ,j+1,k+1) );
    break;
  case TerrainMet::h_eta:
    met_h_eta  = 0.25 * dyInv *
      ( z_nd(i,j+1,k) + z_nd(i+1,j+1,k) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k+1)
       -z_nd(i,j  ,k) - z_nd(i+1,j  ,k) - z_nd(i,j  ,k+1) - z_nd(i+1,j  ,k+1) );
    break;
  case TerrainMet::h_zeta:
    met_h_zeta = 0.25 * dzInv *
      ( z_nd(i+1,j,k+1) + z_nd(i+1,j+1,k+1) + z_nd(i,j,k+1) + z_nd(i,j+1,k+1)
       -z_nd(i+1,j,k  ) - z_nd(i+1,j+1,k  ) - z_nd(i,j,k  ) - z_nd(i,j+1,k  ) );
    break;
  case TerrainMet::h_xi_eta:
    met_h_xi   = 0.25 * dxInv *
      ( z_nd(i+1,j,k) + z_nd(i+1,j+1,k) + z_nd(i+1,j,k+1) + z_nd(i+1,j+1,k+1)
       -z_nd(i  ,j,k) - z_nd(i  ,j+1,k) - z_nd(i  ,j,k+1) - z_nd(i  ,j+1,k+1) );
    met_h_eta  = 0.25 * dyInv *
      ( z_nd(i,j+1,k) + z_nd(i+1,j+1,k) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k+1)
       -z_nd(i,j  ,k) - z_nd(i+1,j  ,k) - z_nd(i,j  ,k+1) - z_nd(i+1,j  ,k+1) );
    break;
  case TerrainMet::h_xi_zeta:
    met_h_xi   = 0.25 * dxInv *
      ( z_nd(i+1,j,k) + z_nd(i+1,j+1,k) + z_nd(i+1,j,k+1) + z_nd(i+1,j+1,k+1)
       -z_nd(i  ,j,k) - z_nd(i  ,j+1,k) - z_nd(i  ,j,k+1) - z_nd(i  ,j+1,k+1) );
    met_h_zeta = 0.25 * dzInv *
      ( z_nd(i+1,j,k+1) + z_nd(i+1,j+1,k+1) + z_nd(i,j,k+1) + z_nd(i,j+1,k+1)
       -z_nd(i+1,j,k  ) - z_nd(i+1,j+1,k  ) - z_nd(i,j,k  ) - z_nd(i,j+1,k  ) );
     break;
  case TerrainMet::h_eta_zeta:
    met_h_eta  = 0.25 * dyInv *
      ( z_nd(i,j+1,k) + z_nd(i+1,j+1,k) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k+1)
       -z_nd(i,j  ,k) - z_nd(i+1,j  ,k) - z_nd(i,j  ,k+1) - z_nd(i+1,j  ,k+1) );
    met_h_zeta = 0.25 * dzInv *
      ( z_nd(i+1,j,k+1) + z_nd(i+1,j+1,k+1) + z_nd(i,j,k+1) + z_nd(i,j+1,k+1)
       -z_nd(i+1,j,k  ) - z_nd(i+1,j+1,k  ) - z_nd(i,j,k  ) - z_nd(i,j+1,k  ) );
    break;
  default:
    met_h_xi   = 0.25 * dxInv *
      ( z_nd(i+1,j,k) + z_nd(i+1,j+1,k) + z_nd(i+1,j,k+1) + z_nd(i+1,j+1,k+1)
       -z_nd(i  ,j,k) - z_nd(i  ,j+1,k) - z_nd(i  ,j,k+1) - z_nd(i  ,j+1,k+1) );
    met_h_eta  = 0.25 * dyInv *
      ( z_nd(i,j+1,k) + z_nd(i+1,j+1,k) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k+1)
       -z_nd(i,j  ,k) - z_nd(i+1,j  ,k) - z_nd(i,j  ,k+1) - z_nd(i+1,j  ,k+1) );
    met_h_zeta = 0.25 * dzInv *
      ( z_nd(i+1,j,k+1) + z_nd(i+1,j+1,k+1) + z_nd(i,j,k+1) + z_nd(i,j+1,k+1)
       -z_nd(i+1,j,k  ) - z_nd(i+1,j+1,k  ) - z_nd(i,j,k  ) - z_nd(i,j+1,k  ) );
  }
}

//*****************************************************************************************
// Compute terrain metric terms at face-centers
//*****************************************************************************************
// Metric coincides with U location
AMREX_GPU_DEVICE
inline void
ComputeMetricAtIface(const int &i, const int &j, const int &k,
                   amrex::Real& met_h_xi,
                   amrex::Real& met_h_eta,
                   amrex::Real& met_h_zeta,
                   const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                   const amrex::Array4<const amrex::Real>& z_nd, const int flag)
{
  met_h_xi   = 0.0;
  met_h_eta  = 0.0;
  met_h_zeta = 0.0;
  amrex::Real dxInv = cellSizeInv[0];
  amrex::Real dyInv = cellSizeInv[1];
  amrex::Real dzInv = cellSizeInv[2];

  switch(flag) {
  case TerrainMet::h_xi:
    met_h_xi   = 0.125 * dxInv *
      ( z_nd(i+1,j,k) + z_nd(i+1,j,k+1) + z_nd(i+1,j+1,k) + z_nd(i+1,j+1,k+1)
       -z_nd(i-1,j,k) - z_nd(i-1,j,k+1) - z_nd(i-1,j+1,k) - z_nd(i-1,j+1,k+1) );
    break;
  case TerrainMet::h_eta:
    met_h_eta  = 0.5 * dyInv * ( z_nd(i,j,k  ) + z_nd(i,j  ,k+1)
                               - z_nd(i,j+1,k) - z_nd(i,j+1,k+1) );
    break;
  case TerrainMet::h_zeta:
    met_h_zeta = 0.5 * dzInv * ( z_nd(i,j,k+1) + z_nd(i,j+1,k+1)
                               - z_nd(i,j,k  ) - z_nd(i,j+1,k  ) );
    break;
  case TerrainMet::h_xi_eta:
    met_h_xi   = 0.125 * dxInv *
      ( z_nd(i+1,j,k) + z_nd(i+1,j,k+1) + z_nd(i+1,j+1,k) + z_nd(i+1,j+1,k+1)
       -z_nd(i-1,j,k) - z_nd(i-1,j,k+1) - z_nd(i-1,j+1,k) - z_nd(i-1,j+1,k+1) );
    met_h_eta  = 0.5 * dyInv * ( z_nd(i,j,k  ) + z_nd(i,j  ,k+1)
                               - z_nd(i,j+1,k) - z_nd(i,j+1,k+1) );
    break;
  case TerrainMet::h_xi_zeta:
    met_h_xi   = 0.125 * dxInv *
      ( z_nd(i+1,j,k) + z_nd(i+1,j,k+1) + z_nd(i+1,j+1,k) + z_nd(i+1,j+1,k+1)
       -z_nd(i-1,j,k) - z_nd(i-1,j,k+1) - z_nd(i-1,j+1,k) - z_nd(i-1,j+1,k+1) );
    met_h_zeta = 0.5 * dzInv * ( z_nd(i,j,k+1) + z_nd(i,j+1,k+1)
                               - z_nd(i,j,k  ) - z_nd(i,j+1,k  ) );
     break;
  case TerrainMet::h_eta_zeta:
    met_h_eta  = 0.5 * dyInv * ( z_nd(i,j,k  ) + z_nd(i,j  ,k+1)
                               - z_nd(i,j+1,k) - z_nd(i,j+1,k+1) );
    met_h_zeta = 0.5 * dzInv * ( z_nd(i,j,k+1) + z_nd(i,j+1,k+1)
                               - z_nd(i,j,k  ) - z_nd(i,j+1,k  ) );
    break;
  default:
    met_h_xi   = 0.125 * dxInv *
      ( z_nd(i+1,j,k) + z_nd(i+1,j,k+1) + z_nd(i+1,j+1,k) + z_nd(i+1,j+1,k+1)
       -z_nd(i-1,j,k) - z_nd(i-1,j,k+1) - z_nd(i-1,j+1,k) - z_nd(i-1,j+1,k+1) );
    met_h_eta  = 0.5 * dyInv * ( z_nd(i,j,k  ) + z_nd(i,j  ,k+1)
                               - z_nd(i,j+1,k) - z_nd(i,j+1,k+1) );
    met_h_zeta = 0.5 * dzInv * ( z_nd(i,j,k+1) + z_nd(i,j+1,k+1)
                               - z_nd(i,j,k  ) - z_nd(i,j+1,k  ) );
  }
}

// Metric coincides with V location
AMREX_GPU_DEVICE
inline void
ComputeMetricAtJface(const int &i, const int &j, const int &k,
                   amrex::Real& met_h_xi,
                   amrex::Real& met_h_eta,
                   amrex::Real& met_h_zeta,
                   const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                   const amrex::Array4<const amrex::Real>& z_nd, const int flag)
{
  met_h_xi   = 0.0;
  met_h_eta  = 0.0;
  met_h_zeta = 0.0;
  amrex::Real dxInv = cellSizeInv[0];
  amrex::Real dyInv = cellSizeInv[1];
  amrex::Real dzInv = cellSizeInv[2];

  switch(flag) {
  case TerrainMet::h_xi:
    met_h_xi   = 0.5 * dxInv * ( z_nd(i+1,j,k) + z_nd(i+1,j,k+1)
                                -z_nd(i  ,j,k) - z_nd(i  ,j,k+1) );
    break;
  case TerrainMet::h_eta:
    met_h_eta  = 0.125 * dyInv *
      ( z_nd(i,j+1,k) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k) + z_nd(i+1,j+1,k+1)
       -z_nd(i,j-1,k) - z_nd(i,j-1,k+1) - z_nd(i+1,j-1,k) - z_nd(i+1,j-1,k+1) );
    break;
  case TerrainMet::h_zeta:
    met_h_zeta = 0.5 * dzInv * ( z_nd(i,j,k+1) + z_nd(i+1,j,k+1)
                               - z_nd(i,j,k  ) - z_nd(i+1,j,k  ) );
    break;
  case TerrainMet::h_xi_eta:
    met_h_xi   = 0.5 * dxInv * ( z_nd(i+1,j,k) + z_nd(i+1,j,k+1)
                                -z_nd(i  ,j,k) - z_nd(i  ,j,k+1) );
    met_h_eta  = 0.125 * dyInv *
      ( z_nd(i,j+1,k) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k) + z_nd(i+1,j+1,k+1)
       -z_nd(i,j-1,k) - z_nd(i,j-1,k+1) - z_nd(i+1,j-1,k) - z_nd(i+1,j-1,k+1) );
    break;
  case TerrainMet::h_xi_zeta:
    met_h_xi   = 0.5 * dxInv * ( z_nd(i+1,j,k) + z_nd(i+1,j,k+1)
                                -z_nd(i  ,j,k) - z_nd(i  ,j,k+1) );
    met_h_zeta = 0.5 * dzInv * ( z_nd(i,j,k+1) + z_nd(i+1,j,k+1)
                               - z_nd(i,j,k  ) - z_nd(i+1,j,k  ) );
     break;
  case TerrainMet::h_eta_zeta:
    met_h_eta  = 0.125 * dyInv *
      ( z_nd(i,j+1,k) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k) + z_nd(i+1,j+1,k+1)
       -z_nd(i,j-1,k) - z_nd(i,j-1,k+1) - z_nd(i+1,j-1,k) - z_nd(i+1,j-1,k+1) );
    met_h_zeta = 0.5 * dzInv * ( z_nd(i,j,k+1) + z_nd(i+1,j,k+1)
                               - z_nd(i,j,k  ) - z_nd(i+1,j,k  ) );
    break;
  default:
    met_h_xi   = 0.5 * dxInv * ( z_nd(i+1,j,k) + z_nd(i+1,j,k+1)
                                -z_nd(i  ,j,k) - z_nd(i  ,j,k+1) );
    met_h_eta  = 0.125 * dyInv *
      ( z_nd(i,j+1,k) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k) + z_nd(i+1,j+1,k+1)
       -z_nd(i,j-1,k) - z_nd(i,j-1,k+1) - z_nd(i+1,j-1,k) - z_nd(i+1,j-1,k+1) );
    met_h_zeta = 0.5 * dzInv * ( z_nd(i,j,k+1) + z_nd(i+1,j,k+1)
                               - z_nd(i,j,k  ) - z_nd(i+1,j,k  ) );
  }
}

// Metric coincides with W location
AMREX_GPU_DEVICE
inline void
ComputeMetricAtKface(const int &i, const int &j, const int &k,
                   amrex::Real& met_h_xi,
                   amrex::Real& met_h_eta,
                   amrex::Real& met_h_zeta,
                   const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                   const amrex::Array4<const amrex::Real>& z_nd, const int flag)
{
  met_h_xi   = 0.0;
  met_h_eta  = 0.0;
  met_h_zeta = 0.0;
  amrex::Real dxInv = cellSizeInv[0];
  amrex::Real dyInv = cellSizeInv[1];
  amrex::Real dzInv = cellSizeInv[2];

  switch(flag) {
  case TerrainMet::h_xi:
    met_h_xi   = 0.5 * dxInv * ( z_nd(i+1,j,k) + z_nd(i+1,j+1,k)
                                -z_nd(i  ,j,k) - z_nd(i  ,j+1,k) );
    break;
  case TerrainMet::h_eta:
    met_h_eta  = 0.5 * dyInv * ( z_nd(i,j+1,k) + z_nd(i+1,j+1,k)
                                -z_nd(i,j  ,k) - z_nd(i+1,j  ,k) );
    break;
  case TerrainMet::h_zeta:
    met_h_zeta = 0.125 * dzInv *
      ( z_nd(i,j,k+1) + z_nd(i+1,j,k+1) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k+1)
       -z_nd(i,j,k-1) - z_nd(i+1,j,k-1) - z_nd(i,j+1,k-1) - z_nd(i+1,j+1,k-1) );
    break;
  case TerrainMet::h_xi_eta:
    met_h_xi   = 0.5 * dxInv * ( z_nd(i+1,j,k) + z_nd(i+1,j+1,k)
                                -z_nd(i  ,j,k) - z_nd(i  ,j+1,k) );
    met_h_eta  = 0.5 * dyInv * ( z_nd(i,j+1,k) + z_nd(i+1,j+1,k)
                                -z_nd(i,j  ,k) - z_nd(i+1,j  ,k) );
    break;
  case TerrainMet::h_xi_zeta:
    met_h_xi   = 0.5 * dxInv * ( z_nd(i+1,j,k) + z_nd(i+1,j+1,k)
                                -z_nd(i  ,j,k) - z_nd(i  ,j+1,k) );
    met_h_zeta = 0.125 * dzInv *
      ( z_nd(i,j,k+1) + z_nd(i+1,j,k+1) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k+1)
       -z_nd(i,j,k-1) - z_nd(i+1,j,k-1) - z_nd(i,j+1,k-1) - z_nd(i+1,j+1,k-1) );
     break;
  case TerrainMet::h_eta_zeta:
    met_h_eta  = 0.5 * dyInv * ( z_nd(i,j+1,k) + z_nd(i+1,j+1,k)
                                -z_nd(i,j  ,k) - z_nd(i+1,j  ,k) );
    met_h_zeta = 0.125 * dzInv *
      ( z_nd(i,j,k+1) + z_nd(i+1,j,k+1) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k+1)
       -z_nd(i,j,k-1) - z_nd(i+1,j,k-1) - z_nd(i,j+1,k-1) - z_nd(i+1,j+1,k-1) );
    break;
  default:
    met_h_xi   = 0.5 * dxInv * ( z_nd(i+1,j,k) + z_nd(i+1,j+1,k)
                                -z_nd(i  ,j,k) - z_nd(i  ,j+1,k) );
    met_h_eta  = 0.5 * dyInv * ( z_nd(i,j+1,k) + z_nd(i+1,j+1,k)
                                -z_nd(i,j  ,k) - z_nd(i+1,j  ,k) );
    met_h_zeta = 0.125 * dzInv *
      ( z_nd(i,j,k+1) + z_nd(i+1,j,k+1) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k+1)
       -z_nd(i,j,k-1) - z_nd(i+1,j,k-1) - z_nd(i,j+1,k-1) - z_nd(i+1,j+1,k-1) );
  }
}

//*****************************************************************************************
// Compute terrain metric terms at edge-centers
//*****************************************************************************************
// Metric is at edge and center Z (red pentagon)
AMREX_GPU_DEVICE
inline void
ComputeMetricAtEdgeCenterK(const int &i, const int &j, const int &k,
                           amrex::Real& met_h_xi,
                           amrex::Real& met_h_eta,
                           amrex::Real& met_h_zeta,
                           const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                           const amrex::Array4<const amrex::Real>& z_nd, const int flag)
{
  met_h_xi   = 0.0;
  met_h_eta  = 0.0;
  met_h_zeta = 0.0;
  amrex::Real dxInv = cellSizeInv[0];
  amrex::Real dyInv = cellSizeInv[1];
  amrex::Real dzInv = cellSizeInv[2];

  switch(flag) {
  case TerrainMet::h_xi:
    met_h_xi  = 0.25 * dxInv *
      ( z_nd(i+1,j,k) + z_nd(i+1,j,k+1)
       -z_nd(i-1,j,k) - z_nd(i-1,j,k+1) );
    break;
  case TerrainMet::h_eta:
    met_h_eta = 0.25 * dyInv *
      ( z_nd(i,j+1,k) + z_nd(i,j+1,k+1)
       -z_nd(i,j-1,k) - z_nd(i,j-1,k+1) );
    break;
  case TerrainMet::h_zeta:
    met_h_zeta = dzInv * (z_nd(i,j,k+1) - z_nd(i,j,k));
    break;
  case TerrainMet::h_xi_eta:
    met_h_xi  = 0.25 * dxInv *
      ( z_nd(i+1,j,k) + z_nd(i+1,j,k+1)
       -z_nd(i-1,j,k) - z_nd(i-1,j,k+1) );
    met_h_eta = 0.25 * dyInv *
      ( z_nd(i,j+1,k) + z_nd(i,j+1,k+1)
       -z_nd(i,j-1,k) - z_nd(i,j-1,k+1) );
    break;
  case TerrainMet::h_xi_zeta:
    met_h_xi  = 0.25 * dxInv *
      ( z_nd(i+1,j,k) + z_nd(i+1,j,k+1)
       -z_nd(i-1,j,k) - z_nd(i-1,j,k+1) );
    met_h_zeta = dzInv * (z_nd(i,j,k+1) - z_nd(i,j,k));
     break;
  case TerrainMet::h_eta_zeta:
    met_h_eta = 0.25 * dyInv *
      ( z_nd(i,j+1,k) + z_nd(i,j+1,k+1)
       -z_nd(i,j-1,k) - z_nd(i,j-1,k+1) );
    met_h_zeta = dzInv * (z_nd(i,j,k+1) - z_nd(i,j,k));
    break;
  default:
    met_h_xi  = 0.25 * dxInv *
      ( z_nd(i+1,j,k) + z_nd(i+1,j,k+1)
       -z_nd(i-1,j,k) - z_nd(i-1,j,k+1) );
    met_h_eta = 0.25 * dyInv *
      ( z_nd(i,j+1,k) + z_nd(i,j+1,k+1)
       -z_nd(i,j-1,k) - z_nd(i,j-1,k+1) );
    met_h_zeta = dzInv * (z_nd(i,j,k+1) - z_nd(i,j,k));
  }
}

// Metric is at edge and center Y (magenta cross)
AMREX_GPU_DEVICE
inline void
ComputeMetricAtEdgeCenterJ(const int &i, const int &j, const int &k,
                           amrex::Real& met_h_xi,
                           amrex::Real& met_h_eta,
                           amrex::Real& met_h_zeta,
                           const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                           const amrex::Array4<const amrex::Real>& z_nd, const int flag)
{
  met_h_xi   = 0.0;
  met_h_eta  = 0.0;
  met_h_zeta = 0.0;
  amrex::Real dxInv = cellSizeInv[0];
  amrex::Real dyInv = cellSizeInv[1];
  amrex::Real dzInv = cellSizeInv[2];

  switch(flag) {
  case TerrainMet::h_xi:
    met_h_xi = 0.25 * dxInv *
      ( z_nd(i+1,j+1,k) + z_nd(i+1,j  ,k)
       -z_nd(i-1,j+1,k) - z_nd(i-1,j  ,k) );
    break;
  case TerrainMet::h_eta:
    met_h_eta = dyInv * ( z_nd(i,j+1,k) - z_nd(i,j,k) );
    break;
  case TerrainMet::h_zeta:
    met_h_zeta = 0.25 * dzInv * ( z_nd(i,j,k+1) + z_nd(i,j+1,k+1)
                                 -z_nd(i,j,k-1) - z_nd(i,j+1,k-1) );
    break;
  case TerrainMet::h_xi_eta:
    met_h_xi = 0.25 * dxInv *
      ( z_nd(i+1,j+1,k) + z_nd(i+1,j  ,k)
       -z_nd(i-1,j+1,k) - z_nd(i-1,j  ,k) );
    met_h_eta = dyInv * ( z_nd(i,j+1,k) - z_nd(i,j,k) );
    break;
  case TerrainMet::h_xi_zeta:
    met_h_xi = 0.25 * dxInv *
      ( z_nd(i+1,j+1,k) + z_nd(i+1,j  ,k)
       -z_nd(i-1,j+1,k) - z_nd(i-1,j  ,k) );
    met_h_zeta = 0.25 * dzInv * ( z_nd(i,j,k+1) + z_nd(i,j+1,k+1)
                                 -z_nd(i,j,k-1) - z_nd(i,j+1,k-1) );
     break;
  case TerrainMet::h_eta_zeta:
    met_h_eta = dyInv * ( z_nd(i,j+1,k) - z_nd(i,j,k) );
    met_h_zeta = 0.25 * dzInv * ( z_nd(i,j,k+1) + z_nd(i,j+1,k+1)
                                 -z_nd(i,j,k-1) - z_nd(i,j+1,k-1) );
    break;
  default:
    met_h_xi = 0.25 * dxInv *
      ( z_nd(i+1,j+1,k) + z_nd(i+1,j  ,k)
       -z_nd(i-1,j+1,k) - z_nd(i-1,j  ,k) );
    met_h_eta = dyInv * ( z_nd(i,j+1,k) - z_nd(i,j,k) );
    met_h_zeta = 0.25 * dzInv * ( z_nd(i,j,k+1) + z_nd(i,j+1,k+1)
                                 -z_nd(i,j,k-1) - z_nd(i,j+1,k-1) );
  }
}

// Metric is at edge and center I (purple hexagon)
AMREX_GPU_DEVICE
inline void
ComputeMetricAtEdgeCenterI(const int &i, const int &j, const int &k,
                           amrex::Real& met_h_xi,
                           amrex::Real& met_h_eta,
                           amrex::Real& met_h_zeta,
                           const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                           const amrex::Array4<const amrex::Real>& z_nd, const int flag)
{
  met_h_xi   = 0.0;
  met_h_eta  = 0.0;
  met_h_zeta = 0.0;
  amrex::Real dxInv = cellSizeInv[0];
  amrex::Real dyInv = cellSizeInv[1];
  amrex::Real dzInv = cellSizeInv[2];

  switch(flag) {
  case TerrainMet::h_xi:
    met_h_xi  = dxInv * ( z_nd(i+1,j,k) - z_nd(i,j,k) );
    break;
  case TerrainMet::h_eta:
    met_h_eta = 0.25 * dyInv *
      ( z_nd(i+1,j+1,k) + z_nd(i,j+1,k)
       -z_nd(i+1,j-1,k) - z_nd(i,j-1,k) );
    break;
  case TerrainMet::h_zeta:
    met_h_zeta = 0.25 * dzInv * ( z_nd(i,j,k+1) + z_nd(i+1,j,k+1)
                                 -z_nd(i,j,k-1) - z_nd(i+1,j,k-1) );
    break;
  case TerrainMet::h_xi_eta:
    met_h_xi  = dxInv * ( z_nd(i+1,j,k) - z_nd(i,j,k) );
    met_h_eta = 0.25 * dyInv *
      ( z_nd(i+1,j+1,k) + z_nd(i,j+1,k)
       -z_nd(i+1,j-1,k) - z_nd(i,j-1,k) );
    break;
  case TerrainMet::h_xi_zeta:
    met_h_xi  = dxInv * ( z_nd(i+1,j,k) - z_nd(i,j,k) );
    met_h_zeta = 0.25 * dzInv * ( z_nd(i,j,k+1) + z_nd(i+1,j,k+1)
                                 -z_nd(i,j,k-1) - z_nd(i+1,j,k-1) );
    break;
  case TerrainMet::h_eta_zeta:
    met_h_eta = 0.25 * dyInv *
      ( z_nd(i+1,j+1,k) + z_nd(i,j+1,k)
       -z_nd(i+1,j-1,k) - z_nd(i,j-1,k) );
    met_h_zeta = 0.25 * dzInv * ( z_nd(i,j,k+1) + z_nd(i+1,j,k+1)
                                 -z_nd(i,j,k-1) - z_nd(i+1,j,k-1) );
    break;
  default:
    met_h_xi  = dxInv * ( z_nd(i+1,j,k) - z_nd(i,j,k) );
    met_h_eta = 0.25 * dyInv *
      ( z_nd(i+1,j+1,k) + z_nd(i,j+1,k)
       -z_nd(i+1,j-1,k) - z_nd(i,j-1,k) );
    met_h_zeta = 0.25 * dzInv * ( z_nd(i,j,k+1) + z_nd(i+1,j,k+1)
                                 -z_nd(i,j,k-1) - z_nd(i+1,j,k-1) );
  }
}

//*****************************************************************************************
// Map between W <--> Omega
//*****************************************************************************************
AMREX_GPU_DEVICE
inline
amrex::Real OmegaFromW(int i, int j, int k, amrex::Real w,
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

  // Use extrapolation instead of interpolation if at the bottom boundary
  
  amrex::Real vec_zlo_xi   = (k == 0) ? 1.5 * (0.5*(u(i,j,k)+u(i+1,j,k))) - 0.5*(0.5*(u(i,j,k+1)+u(i+1,j,k+1))) :
    0.25 * ( u(i,j,k-1) + u(i+1,j,k-1) + u(i,j,k) + u(i+1,j,k) );
  amrex::Real vec_zlo_eta  = (k == 0) ? 1.5 * (0.5*( v(i,j,k)+v(i,j+1,k))) - 0.5*(0.5*(v(i,j,k+1)+v(i,j+1,k+1))) :
    0.25 * ( v(i,j,k-1) + v(i,j+1,k-1) + v(i,j,k) + v(i,j+1,k) );  

  amrex::Real omega = -met_zlo_xi * vec_zlo_xi - met_zlo_eta * vec_zlo_eta + w;
  return omega;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real WFromOmega(int i, int j, int k, amrex::Real omega,
                       const amrex::Array4<const amrex::Real>& u,
                       const amrex::Array4<const amrex::Real>& v,
                       const amrex::Array4<const amrex::Real>& z_nd,
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

  // Use extrapolation instead of interpolation if at the bottom boundary
  
  amrex::Real vec_zlo_xi   = (k == 0) ? 1.5 * (0.5*(u(i,j,k)+u(i+1,j,k))) - 0.5*(0.5*(u(i,j,k+1)+u(i+1,j,k+1))) :
    0.25 * ( u(i,j,k-1) + u(i+1,j,k-1) + u(i,j,k) + u(i+1,j,k) );
  amrex::Real vec_zlo_eta  = (k == 0) ? 1.5 * (0.5*( v(i,j,k)+v(i,j+1,k))) - 0.5*(0.5*(v(i,j,k+1)+v(i,j+1,k+1))) :
    0.25 * ( v(i,j,k-1) + v(i,j+1,k-1) + v(i,j,k) + v(i,j+1,k) );  

  amrex::Real w = met_zlo_xi * vec_zlo_xi + met_zlo_eta * vec_zlo_eta + omega;
  return w;
}

// This routine follows WFromOmega with clipping at the boundary:
// Terrain metrics at the boundary are used for ghost cells (0 gradient).
// U & V are extrapolated to the top & bottom k faces.
AMREX_GPU_DEVICE
inline
amrex::Real WFromOmegaBC(int i, int j, int k, amrex::Real omega,
                         const amrex::Array4<const amrex::Real>& u,
                         const amrex::Array4<const amrex::Real>& v,
                         const amrex::Array4<const amrex::Real>& z_nd,
                         const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxInv)
{
  // Clip metric indices for ghost cell BCs
  amrex::Dim3 hi = ubound(z_nd);
  amrex::Dim3 lo = lbound(z_nd);
  int  im = (i  >= hi.x) ? (hi.x-1) : i;
       im = (im <= lo.x) ? (lo.x+1) : im;
  int  jm = (j  >= hi.y) ? (hi.y-1) : j;
       jm = (jm <= lo.y) ? (lo.y+1) : jm;
  int  km = (k  >= hi.z) ? (hi.z-1) : k;
       km = (km <= lo.z) ? (lo.z+1) : km;

  // This is dh/dxi at z-face (i,j,k-1/2)
  amrex::Real met_zlo_xi   = 0.5 * dxInv[0] *
    ( z_nd(im+1,jm+1,km  ) + z_nd(im+1,jm  ,km  )    // hi i, lo k
     -z_nd(im  ,jm+1,km  ) - z_nd(im  ,jm  ,km  ) ); // lo i, lo k

  // This is dh/deta at z-face (i,j,k-1/2)
  amrex::Real met_zlo_eta  = 0.5 * dxInv[1] *
    ( z_nd(im+1,jm+1,km  ) + z_nd(im  ,jm+1,km  )    // hi j, lo k
     -z_nd(im+1,jm  ,km  ) - z_nd(im  ,jm  ,km  ) ); // lo j, lo k

  // Clip velocity indices for ghost cell BCs
  amrex::Real vec_zlo_xi;
  amrex::Real vec_zlo_eta;
  hi = ubound(u);
  lo = lbound(u);
  hi.z = hi.z + 1; // nodal for W-velocity

  if( k == hi.z) {
    // foextrap to highest k-face
    vec_zlo_xi  = 0.75*( u(i,j,k-1)+u(i+1,j,k-1) ) - 0.25*( u(i,j,k-2)+u(i+1,j,k-2) );
    vec_zlo_eta = 0.75*( v(i,j,k-1)+v(i,j+1,k-1) ) - 0.25*( v(i,j,k-2)+v(i,j+1,k-2) );
  } else if (k == lo.z) {
    // foextrap to lowest k-face
    vec_zlo_xi  = 0.75*( u(i,j,k)+u(i+1,j,k) ) - 0.25*( u(i,j,k+1)+u(i+1,j,k+1) );
    vec_zlo_eta = 0.75*( v(i,j,k)+v(i,j+1,k) ) - 0.25*( v(i,j,k+1)+v(i,j+1,k+1) );
  } else {
    // Valid ghost cell for u & v
    vec_zlo_xi  = 0.25 * ( u(i,j,k-1) + u(i+1,j,k-1) + u(i,j,k) + u(i+1,j,k) );
    vec_zlo_eta = 0.25 * ( v(i,j,k-1) + v(i,j+1,k-1) + v(i,j,k) + v(i,j+1,k) );
  }

  amrex::Real w = met_zlo_xi * vec_zlo_xi + met_zlo_eta * vec_zlo_eta + omega;
  return w;
}
#endif
