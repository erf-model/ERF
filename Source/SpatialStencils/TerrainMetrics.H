#ifndef _TERRAIN_METRIC_H_
#define _TERRAIN_METRIC_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include "IndexDefines.H"

//*****************************************************************************************
// Compute terrain metric terms for face-centered quantities
//*****************************************************************************************
AMREX_GPU_DEVICE
inline void
ComputeMetricXmomXface(const int &i, const int &j, const int &k,
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
    met_h_xi   = 0.0;
    break;
  case TerrainMet::h_eta:
    met_h_eta  = 0.0;
    break;
  case TerrainMet::h_zeta:

    // dh/dzeta for X-momentum at X-face
    // I index backwards for consistency with DiffusionSrcForMom
    met_h_zeta = 0.25 * dzInv *
      ( z_nd(i,j,k+1) + z_nd(i,j+1,k+1) + z_nd(i-1,j,k+1) + z_nd(i-1,j+1,k+1)
       -z_nd(i,j,k  ) - z_nd(i,j+1,k  ) - z_nd(i-1,j,k  ) - z_nd(i-1,j+1,k  ) );
    
    break;
  default:
    met_h_xi   = 0.0;
    met_h_eta  = 0.0;
    met_h_zeta = 0.25 * dzInv *
      ( z_nd(i,j,k+1) + z_nd(i,j+1,k+1) + z_nd(i-1,j,k+1) + z_nd(i-1,j+1,k+1)
       -z_nd(i,j,k  ) - z_nd(i,j+1,k  ) - z_nd(i-1,j,k  ) - z_nd(i-1,j+1,k  ) );
  }
}

AMREX_GPU_DEVICE
inline void
ComputeMetricXmomYface(const int &i, const int &j, const int &k,
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
    met_h_xi   = 0.0;
    break;
  case TerrainMet::h_eta:
    met_h_eta  = 0.0;
    break;
  case TerrainMet::h_zeta:
    
    // dh/dzeta for X-momentum at Y-face
    // J indexing done in DiffusionSrcForMom
    met_h_zeta = dzInv * (z_nd(i,j,k+1) - z_nd(i,j,k));

    break;
  default:
    met_h_xi   = 0.0;
    met_h_eta  = 0.0;
    met_h_zeta = dzInv * (z_nd(i,j,k+1) - z_nd(i,j,k));
  }
}

AMREX_GPU_DEVICE
inline void
ComputeMetricXmomZface(const int &i, const int &j, const int &k,
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
    met_h_xi   = 0.0;
    break;
  case TerrainMet::h_eta:
    met_h_eta  = 0.0;
    break;
  case TerrainMet::h_zeta:

    // dh/dzeta for X-momentum at Z-face
    // K indexing done in DiffusionSrcForMom
    met_h_zeta = 0.25 * dzInv *
      ( z_nd(i,j,k+1) + z_nd(i,j+1,k+1)
       -z_nd(i,j,k-1) - z_nd(i,j+1,k-1) );

    break;
  default:
    met_h_xi   = 0.0;
    met_h_eta  = 0.0;
    met_h_zeta = 0.25 * dzInv *
      ( z_nd(i,j,k+1) + z_nd(i,j+1,k+1)
       -z_nd(i,j,k-1) - z_nd(i,j+1,k-1) );
  }
}

AMREX_GPU_DEVICE
inline void
ComputeMetricYmomXface(const int &i, const int &j, const int &k,
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
    met_h_xi   = 0.0;
    break;
  case TerrainMet::h_eta:
    met_h_eta  = 0.0;
    break;
  case TerrainMet::h_zeta:

    // dh/dzeta for Y-momentum at X-face
    // I indexing done in DiffusionSrcForMom
    met_h_zeta = dzInv * (z_nd(i,j,k+1) - z_nd(i,j,k));

    break;
  default:
    met_h_xi   = 0.0;
    met_h_eta  = 0.0;
    met_h_zeta = dzInv * (z_nd(i,j,k+1) - z_nd(i,j,k));
  }
}

AMREX_GPU_DEVICE
inline void
ComputeMetricYmomYface(const int &i, const int &j, const int &k,
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
    met_h_xi   = 0.0;
    break;
  case TerrainMet::h_eta:
    met_h_eta  = 0.0;
    break;
  case TerrainMet::h_zeta:

    // dh/dzeta for Y-momentum at Y-face
    // J index backwards for consistency with DiffusionSrcForMom
    met_h_zeta = 0.25 * dzInv *
      ( z_nd(i,j,k+1) + z_nd(i+1,j,k+1) + z_nd(i,j-1,k+1) + z_nd(i+1,j-1,k+1)
       -z_nd(i,j,k  ) - z_nd(i+1,j,k  ) - z_nd(i,j-1,k  ) - z_nd(i+1,j-1,k  ) );

    break;
  default:
    met_h_xi   = 0.0;
    met_h_eta  = 0.0;
    met_h_zeta = 0.25 * dzInv *
      ( z_nd(i,j,k+1) + z_nd(i+1,j,k+1) + z_nd(i,j-1,k+1) + z_nd(i+1,j-1,k+1)
       -z_nd(i,j,k  ) - z_nd(i+1,j,k  ) - z_nd(i,j-1,k  ) - z_nd(i+1,j-1,k  ) );
  }
}

AMREX_GPU_DEVICE
inline void
ComputeMetricYmomZface(const int &i, const int &j, const int &k,
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
    met_h_xi   = 0.0;
    break;
  case TerrainMet::h_eta:
    met_h_eta  = 0.0;
    break;
  case TerrainMet::h_zeta:

    // dh/dzeta for Y-momentum at Z-face
    // K indexing done in DiffusionSrcForMom
    met_h_zeta = 0.25 * dzInv *
      ( z_nd(i,j,k+1) + z_nd(i+1,j,k+1)
       -z_nd(i,j,k-1) - z_nd(i+1,j,k-1) );

    break;
  default:
    met_h_xi   = 0.0;
    met_h_eta  = 0.0;
    met_h_zeta = 0.25 * dzInv *
      ( z_nd(i,j,k+1) + z_nd(i+1,j,k+1)
       -z_nd(i,j,k-1) - z_nd(i+1,j,k-1) );
  }
}

AMREX_GPU_DEVICE
inline void
ComputeMetricZmomXface(const int &i, const int &j, const int &k,
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
    met_h_xi   = 0.0;
    break;
  case TerrainMet::h_eta:
    met_h_eta  = 0.0;
    break;
  case TerrainMet::h_zeta:

    // dh/dzeta for Z-momentum at X-face
    // I indexing done in DiffusionSrcForMom
    met_h_zeta = 0.25 * dzInv *
      ( z_nd(i,j,k+1) + z_nd(i,j+1,k+1)
       -z_nd(i,j,k-1) - z_nd(i,j+1,k-1));

    break;
  default:
    met_h_xi   = 0.0;
    met_h_eta  = 0.0;
    met_h_zeta = 0.25 * dzInv *
      ( z_nd(i,j,k+1) + z_nd(i,j+1,k+1)
       -z_nd(i,j,k-1) - z_nd(i,j+1,k-1));
  }
}

AMREX_GPU_DEVICE
inline void
ComputeMetricZmomYface(const int &i, const int &j, const int &k,
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
    met_h_xi   = 0.0;
    break;
  case TerrainMet::h_eta:
    met_h_eta  = 0.0;
    break;
  case TerrainMet::h_zeta:

    // dh/dzeta for Z-momentum at Y-face
    // J indexing done in DiffusionSrcForMom
    met_h_zeta = 0.25 * dzInv *
      ( z_nd(i,j,k+1) + z_nd(i+1,j,k+1)
       -z_nd(i,j,k-1) - z_nd(i+1,j,k-1));

    break;
  default:
    met_h_xi   = 0.0;
    met_h_eta  = 0.0;
    met_h_zeta = 0.25 * dzInv *
      ( z_nd(i,j,k+1) + z_nd(i+1,j,k+1)
       -z_nd(i,j,k-1) - z_nd(i+1,j,k-1));
  }
}

AMREX_GPU_DEVICE
inline void
ComputeMetricZmomZface(const int &i, const int &j, const int &k,
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
    met_h_xi   = 0.0;
    break;
  case TerrainMet::h_eta:
    met_h_eta  = 0.0;
    break;
  case TerrainMet::h_zeta:

    // dh/dzeta for Z-momentum at Z-face
    // K index indexing done in DiffusionSrcForMom
    met_h_zeta = 0.25 * dzInv *
      ( z_nd(i,j,k  ) + z_nd(i+1,j,k  ) + z_nd(i,j+1,k  ) + z_nd(i+1,j+1,k  )
       -z_nd(i,j,k-1) - z_nd(i+1,j,k-1) - z_nd(i,j+1,k-1) - z_nd(i+1,j+1,k-1) );

    break;
  default:
    met_h_xi   = 0.0;
    met_h_eta  = 0.0;
    met_h_zeta = 0.25 * dzInv *
      ( z_nd(i,j,k  ) + z_nd(i+1,j,k  ) + z_nd(i,j+1,k  ) + z_nd(i+1,j+1,k  )
       -z_nd(i,j,k-1) - z_nd(i+1,j,k-1) - z_nd(i,j+1,k-1) - z_nd(i+1,j+1,k-1) );
  }
}
#endif
