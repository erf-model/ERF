#ifndef _COMPUTE_STRAINRATE_H_
#define _COMPUTE_STRAINRATE_H_

#include <AMReX_BCRec.H>
#include <IndexDefines.H>
#include <TerrainMetrics.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
bool
is_on_face(int bc_in)
{
    if ( (bc_in == ERFBCType::ext_dir) || (bc_in == ERFBCType::ext_dir_ingested) ) {
        return true;
    } else {
        return false;
    }
}

// Note:  if the boundary condition type is ext_dir, it means that the velocity is fixed
//        on that face, so we use a special stencil that understands the value is
//        at the face, not at half a grid cell away from the face.   The stencils below also
//        allow for the Dirichlet value to be non-zero.  A wall boundary condition will generate
//        Dirichlet bc in the z-direction; inflow will generate Dirichlet bc in the lateral directions.

// AML: Stress tensor only implemented for the case h_xi = h_eta = 0!

//*****************************************************************************************
// Compute strain rate terms at cell faces
//*****************************************************************************************
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ComputeS11WithTerrain (int i, int j, int k,
                       const amrex::Array4<amrex::Real const>& u,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       amrex::Real met_h_xi, amrex::Real /*met_h_eta*/, amrex::Real met_h_zeta)
{
  amrex::Real GradUz;
  amrex::Real strainRate = 0;
  amrex::Real dx_inv = cellSizeInv[0];
  amrex::Real dz_inv = cellSizeInv[2];

  // At cell center
  GradUz = 0.25 * dz_inv *
    ( u(i  ,j  ,k+1) + u(i-1,j  ,k+1)
     -u(i  ,j  ,k-1) - u(i-1,j  ,k-1) );

  strainRate = ( u(i  ,j  ,k  ) - u(i-1,j  ,k  ) ) * dx_inv
             - (met_h_xi/met_h_zeta) * GradUz;

  return strainRate;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ComputeS22WithTerrain (int i, int j, int k,
                       const amrex::Array4<amrex::Real const>& v,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       amrex::Real /*met_h_xi*/, amrex::Real met_h_eta, amrex::Real met_h_zeta)
{
  amrex::Real GradVz;
  amrex::Real strainRate = 0;
  amrex::Real dy_inv = cellSizeInv[1];
  amrex::Real dz_inv = cellSizeInv[2];

  // At cell center
  GradVz = 0.25 * dz_inv *
    ( v(i  ,j  ,k+1) + v(i  ,j-1,k+1)
     -v(i  ,j  ,k-1) - v(i  ,j-1,k-1) );

  strainRate = (v(i, j, k) - v(i, j-1, k)) * dy_inv
             - (met_h_eta/met_h_zeta) * GradVz;

  return strainRate;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ComputeS33WithTerrain (int i, int j, int k,
                       const amrex::Array4<amrex::Real const>& w,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       amrex::Real /*met_h_xi*/, amrex::Real /*met_h_eta*/, amrex::Real met_h_zeta,
                       const amrex::Box& domain)
{
  amrex::Real strainRate = 0;
  amrex::Real dz_inv = cellSizeInv[2];

  if (k > domain.smallEnd(2) && k <= domain.bigEnd(2)+1)
      strainRate = ( w(i, j, k) - w(i, j, k-1) ) * dz_inv / met_h_zeta ;
  else
      strainRate = 0.;

  return strainRate;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ComputeS12WithTerrain (int i, int j, int k,
                       const amrex::Array4<amrex::Real const>& u,
                       const amrex::Array4<amrex::Real const>& v,
                       const amrex::Array4<amrex::Real const>& /*w*/,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       amrex::Real met_h_xi, amrex::Real met_h_eta, amrex::Real met_h_zeta,
                       const amrex::Box& domain, const amrex::BCRec* bc_ptr)
{
  amrex::Real strainRate = 0;
  amrex::Real GradUz,GradVz;
  amrex::Real dx_inv = cellSizeInv[0];
  amrex::Real dy_inv = cellSizeInv[1];
  amrex::Real dz_inv = cellSizeInv[2];

  if (j == domain.smallEnd(1) && is_on_face(bc_ptr[BCVars::xvel_bc].lo(1))) {
    strainRate  = (-(8./3.) * u(i,j-1,k) + 3. * u(i,j,k) -
                    (1./3.) * u(i,j+1,k))*dy_inv;
  } else if (j == domain.bigEnd(1)+1 && is_on_face(bc_ptr[BCVars::xvel_bc].hi(1))) {
    strainRate  = -(-(8./3.) * u(i,j,k) + 3. * u(i,j-1,k) -
                     (1./3.) * u(i,j-2,k))*dy_inv;
  } else {
    strainRate  =  (u(i, j, k) - u(i, j-1, k)) * dy_inv;
  }
  if (i == domain.smallEnd(0) && is_on_face(bc_ptr[BCVars::yvel_bc].lo(0))) {
    strainRate += (-(8./3.) * v(i-1,j,k) + 3. * v(i,j,k) -
                    (1./3.) * v(i+1,j,k))*dx_inv;
  } else if (i == domain.bigEnd(0)+1 && is_on_face(bc_ptr[BCVars::yvel_bc].hi(0))) {
    strainRate += -(-(8./3.) * v(i,j,k) + 3. * v(i-1,j,k) -
                     (1./3.) * v(i-2,j,k))*dx_inv;
  } else {
    strainRate +=  (v(i, j, k) - v(i-1, j, k)) * dx_inv;
  }

  // At EdgeCenterK
  GradUz = 0.25 * dz_inv *
    ( u(i  ,j  ,k+1) + u(i  ,j-1,k+1)
     -u(i  ,j  ,k-1) - u(i  ,j-1,k-1) );
  GradVz = 0.25 * dz_inv *
    ( v(i  ,j  ,k+1) + v(i-1,j  ,k+1)
     -v(i  ,j  ,k-1) - v(i-1,j  ,k-1) );

  strainRate += -(met_h_eta/met_h_zeta)*GradUz - (met_h_xi/met_h_zeta)*GradVz;

  strainRate *= 0.5;

  return strainRate;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ComputeS13WithTerrain (int i, int j, int k,
                       const amrex::Array4<amrex::Real const>& u,
                       const amrex::Array4<amrex::Real const>& /*v*/,
                       const amrex::Array4<amrex::Real const>& w,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       amrex::Real met_h_xi, amrex::Real /*met_h_eta*/, amrex::Real met_h_zeta,
                       const amrex::Box& domain, const amrex::BCRec* bc_ptr)
{
  amrex::Real GradWz;
  amrex::Real strainRate = 0;
  amrex::Real dx_inv = cellSizeInv[0];
  amrex::Real dz_inv = cellSizeInv[2];

  if (k == domain.smallEnd(2) && is_on_face(bc_ptr[BCVars::xvel_bc].lo(2))) {
    strainRate  =  (-(8./3.) * u(i,j,k-1) + 3. * u(i,j,k) -
                     (1./3.) * u(i,j,k+1)) * dz_inv / met_h_zeta;
  } else if (k == domain.bigEnd(2)+1 && is_on_face(bc_ptr[BCVars::xvel_bc].hi(2))) {
    strainRate  = -(-(8./3.) * u(i,j,k) + 3. * u(i,j,k-1) -
                     (1./3.) * u(i,j,k-2)) * dz_inv / met_h_zeta;
  } else {
    strainRate  =  (u(i, j, k) - u(i, j, k-1)) * dz_inv / met_h_zeta;
  }
  if (i == domain.smallEnd(0) && is_on_face(bc_ptr[BCVars::zvel_bc].lo(0))) {
    strainRate +=  (-(8./3.) * w(i-1,j,k) + 3. * w(i,j,k) - (1./3.) * w(i+1,j,k))*dx_inv;
  } else if (i == domain.bigEnd(0)+1 && is_on_face(bc_ptr[BCVars::zvel_bc].hi(0))) {
    strainRate += -(-(8./3.) * w(i,j,k) + 3. * w(i-1,j,k) - (1./3.) * w(i-2,j,k))*dx_inv;
  } else {
    strainRate +=  (w(i, j, k) - w(i-1, j, k)) * dx_inv;
  }

  // At EdgeCenterJ
  if (k == 0)
      GradWz = 0.5  * dz_inv * ( w(i  ,j  ,k+1) + w(i-1,j  ,k+1) -w(i  ,j  ,k  ) - w(i-1,j  ,k  ) );
  else
      GradWz = 0.25 * dz_inv * ( w(i  ,j  ,k+1) + w(i-1,j  ,k+1) -w(i  ,j  ,k-1) - w(i-1,j  ,k-1) );

  strainRate += -(met_h_xi/met_h_zeta)*GradWz;

  strainRate *= 0.5;

  return strainRate;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ComputeS23WithTerrain (int i, int j, int k,
                       const amrex::Array4<amrex::Real const>& /*u*/,
                       const amrex::Array4<amrex::Real const>& v,
                       const amrex::Array4<amrex::Real const>& w,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       amrex::Real /*met_h_xi*/, amrex::Real met_h_eta, amrex::Real met_h_zeta,
                       const amrex::Box& domain, const amrex::BCRec* bc_ptr)
{
  amrex::Real GradWz;
  amrex::Real strainRate = 0;
  amrex::Real dy_inv = cellSizeInv[1];
  amrex::Real dz_inv = cellSizeInv[2];

  if (k == domain.smallEnd(2) && is_on_face(bc_ptr[BCVars::yvel_bc].lo(2))) {
    strainRate  =  (-(8./3.) * v(i,j,k-1) + 3. * v(i,j,k  ) -
                     (1./3.) * v(i,j,k+1)) * dz_inv / met_h_zeta;
  } else if (k == domain.bigEnd(2)+1 && is_on_face(bc_ptr[BCVars::yvel_bc].hi(2))) {
    strainRate  = -(-(8./3.) * v(i,j,k  ) + 3. * v(i,j,k-1) -
                     (1./3.) * v(i,j,k-2)) * dz_inv / met_h_zeta;
  } else {
    strainRate  =  (v(i, j, k) - v(i, j, k-1)) * dz_inv / met_h_zeta;
  }
  if (j == domain.smallEnd(1) && is_on_face(bc_ptr[BCVars::zvel_bc].lo(1))) {
    strainRate +=  (-(8./3.) * w(i,j-1,k) + 3. * w(i,j  ,k) - (1./3.) * w(i,j+1,k))*dy_inv;
  } else if (j == domain.bigEnd(1)+1 && is_on_face(bc_ptr[BCVars::zvel_bc].hi(1))) {
    strainRate += -(-(8./3.) * w(i,j  ,k) + 3. * w(i,j-1,k) - (1./3.) * w(i,j-2,k))*dy_inv;
  } else {
    strainRate +=  (w(i, j, k) - w(i, j-1, k)) * dy_inv;
  }

  // At EdgeCenterI
  if (k == 0)
      GradWz = 0.5  * dz_inv * ( w(i  ,j  ,k+1) + w(i  ,j-1,k+1) -w(i  ,j  ,k  ) - w(i  ,j-1,k  ) );
  else
      GradWz = 0.25 * dz_inv * ( w(i  ,j  ,k+1) + w(i  ,j-1,k+1) -w(i  ,j  ,k-1) - w(i  ,j-1,k-1) );

  strainRate += -(met_h_eta/met_h_zeta)*GradWz;

  strainRate *= 0.5;

  return strainRate;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ComputeS21WithTerrain (int i, int j, int k,
                       const amrex::Array4<amrex::Real const>& u,
                       const amrex::Array4<amrex::Real const>& v,
                       const amrex::Array4<amrex::Real const>& w,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       amrex::Real met_h_xi, amrex::Real met_h_eta, amrex::Real met_h_zeta,
                       const amrex::Box& domain, const amrex::BCRec* bc_ptr)
{
    return ComputeS12WithTerrain(i,j,k,u,v,w,cellSizeInv,
                                 met_h_xi,met_h_eta,met_h_zeta,domain,bc_ptr);
}

// Metrics introduce asymmetry into strain-rate
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ComputeS31WithTerrain (int i, int j, int k,
                       const amrex::Array4<amrex::Real const>& u,
                       const amrex::Array4<amrex::Real const>& v,
                       const amrex::Array4<amrex::Real const>& w,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       amrex::Real met_h_xi, amrex::Real met_h_eta, amrex::Real met_h_zeta,
                       const amrex::Box& domain, const amrex::BCRec* bc_ptr)
{
   return ComputeS13WithTerrain(i,j,k,u,v,w,cellSizeInv,
                                met_h_xi,met_h_eta,met_h_zeta,domain,bc_ptr);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ComputeS32WithTerrain (int i, int j, int k,
                       const amrex::Array4<amrex::Real const>& u,
                       const amrex::Array4<amrex::Real const>& v,
                       const amrex::Array4<amrex::Real const>& w,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       amrex::Real met_h_xi, amrex::Real met_h_eta, amrex::Real met_h_zeta,
                       const amrex::Box& domain, const amrex::BCRec* bc_ptr)
{
    return ComputeS23WithTerrain(i,j,k,u,v,w,cellSizeInv,
                                 met_h_xi,met_h_eta,met_h_zeta,domain,bc_ptr);
}

//*****************************************************************************************
// Driver for above strain rate functions
//****************************************************************************************
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ComputeStrainRateWithTerrain (int i, int j, int k,
                              const amrex::Array4<amrex::Real const>& u,
                              const amrex::Array4<amrex::Real const>& v,
                              const amrex::Array4<amrex::Real const>& w,
                              enum MomentumEqn  momentumEqn,
                              enum DiffusionDir diffDir,
                              const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                              const amrex::Array4<const amrex::Real>& z_nd,
                              const amrex::Box& domain, const amrex::BCRec* bc_ptr)
{
  amrex::Real strainRate;
  amrex::Real met_h_xi,met_h_eta,met_h_zeta;
  switch (momentumEqn) {
  case MomentumEqn::x:
    switch (diffDir) {
    case DiffusionDir::x: // S11
      // Cell-Center is staggered
      ComputeMetricAtCellCenter(i-1,j,k,met_h_xi,met_h_eta,met_h_zeta,
                                cellSizeInv,z_nd,TerrainMet::h_xi_zeta);
      strainRate = ComputeS11WithTerrain(i,j,k,u,cellSizeInv,met_h_xi,met_h_eta,met_h_zeta);
      break;
    case DiffusionDir::y: // S12
      ComputeMetricAtEdgeCenterK(i,j,k,met_h_xi,met_h_eta,met_h_zeta,
                                 cellSizeInv,z_nd,TerrainMet::all);
      strainRate = ComputeS12WithTerrain(i,j,k,u,v,w,cellSizeInv,
                              met_h_xi,met_h_eta,met_h_zeta,
                              domain,bc_ptr);
      break;
    case DiffusionDir::z: // S13
      ComputeMetricAtEdgeCenterJ(i,j,k,met_h_xi,met_h_eta,met_h_zeta,
                                 cellSizeInv,z_nd,TerrainMet::h_xi_zeta);
      strainRate = ComputeS13WithTerrain(i,j,k,u,v,w,cellSizeInv,
                              met_h_xi,met_h_eta,met_h_zeta,
                              domain,bc_ptr);
      break;
    }
    break;
  case MomentumEqn::y:
    switch (diffDir) {
    case DiffusionDir::x: // S21
      ComputeMetricAtEdgeCenterK(i,j,k,met_h_xi,met_h_eta,met_h_zeta,
                                 cellSizeInv,z_nd,TerrainMet::all);
      strainRate = ComputeS21WithTerrain(i,j,k,u,v,w,cellSizeInv,
                              met_h_xi,met_h_eta,met_h_zeta,
                              domain,bc_ptr);
      break;
    case DiffusionDir::y: // S22
      // Cell-Center is staggered
      ComputeMetricAtCellCenter(i,j-1,k,met_h_xi,met_h_eta,met_h_zeta,
                                cellSizeInv,z_nd,TerrainMet::h_eta_zeta);
      strainRate = ComputeS22WithTerrain(i,j,k,v,cellSizeInv,
                              met_h_xi,met_h_eta,met_h_zeta);
      break;
    case DiffusionDir::z: // S23
      ComputeMetricAtEdgeCenterI(i,j,k,met_h_xi,met_h_eta,met_h_zeta,
                                 cellSizeInv,z_nd,TerrainMet::h_eta_zeta);
      strainRate = ComputeS23WithTerrain(i,j,k,u,v,w,cellSizeInv,
                              met_h_xi,met_h_eta,met_h_zeta,
                              domain,bc_ptr);
      break;
    }
    break;
  case MomentumEqn::z:
    switch (diffDir) {
    case DiffusionDir::x: // S31
      ComputeMetricAtEdgeCenterJ(i,j,k,met_h_xi,met_h_eta,met_h_zeta,
                                 cellSizeInv,z_nd,TerrainMet::h_xi_zeta);
      strainRate = ComputeS31WithTerrain(i,j,k,u,v,w,cellSizeInv,
                              met_h_xi,met_h_eta,met_h_zeta,
                              domain,bc_ptr);
      break;
    case DiffusionDir::y: // S32
      ComputeMetricAtEdgeCenterI(i,j,k,met_h_xi,met_h_eta,met_h_zeta,
                                 cellSizeInv,z_nd,TerrainMet::h_eta_zeta);
      strainRate = ComputeS32WithTerrain(i,j,k,u,v,w,cellSizeInv,
                              met_h_xi,met_h_eta,met_h_zeta,
                              domain,bc_ptr);
      break;
    case DiffusionDir::z: // S33
      // Cell-Center is staggered
      ComputeMetricAtCellCenter(i,j,k-1,met_h_xi,met_h_eta,met_h_zeta,
                                cellSizeInv,z_nd,TerrainMet::h_zeta);
      strainRate = ComputeS33WithTerrain(i,j,k,w,cellSizeInv,
                              met_h_xi,met_h_eta,met_h_zeta,domain);
      break;
    }
    break;
  }

  return strainRate;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ComputeS11 (int i, int j, int k,
            const amrex::Array4<amrex::Real const>& u,
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv)
{
    amrex::Real dx_inv = cellSizeInv[0];
    return (u(i, j, k) - u(i-1, j, k))*dx_inv;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ComputeS22 (int i, int j, int k,
            const amrex::Array4<amrex::Real const>& v,
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv)
{
    amrex::Real dy_inv = cellSizeInv[1];
    return (v(i, j, k) - v(i, j-1, k))*dy_inv;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ComputeS33 (int i, int j, int k,
            const amrex::Array4<amrex::Real const>& w,
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
            const amrex::Box& domain)
{
    amrex::Real dz_inv = cellSizeInv[2];
    if (k > domain.smallEnd(2) && k <= domain.bigEnd(2)+1)
        return (w(i, j, k) - w(i, j, k-1))*dz_inv;
    else
        return 0.;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ComputeS12 (int i, int j, int k,
            const amrex::Array4<amrex::Real const>& u,
            const amrex::Array4<amrex::Real const>& v,
            const amrex::Array4<amrex::Real const>& /*w*/,
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
            const amrex::Box& domain, const amrex::BCRec* bc_ptr)
{
    amrex::Real dx_inv = cellSizeInv[0];
    amrex::Real dy_inv = cellSizeInv[1];
    amrex::Real strainRate = 0;

    if (j == domain.smallEnd(1) && is_on_face(bc_ptr[BCVars::xvel_bc].lo(1))) {
        strainRate =  (-(8./3.) * u(i,j-1,k) + 3. * u(i,j,k) - (1./3.) * u(i,j+1,k))*dy_inv;
    } else if (j == domain.bigEnd(1)+1 && is_on_face(bc_ptr[BCVars::xvel_bc].hi(1))) {
        strainRate = -(-(8./3.) * u(i,j,k) + 3. * u(i,j-1,k) - (1./3.) * u(i,j-2,k))*dy_inv;
    } else {
        strainRate = (u(i, j, k) - u(i, j-1, k))*dy_inv;
    }
    if (i == domain.smallEnd(0) && is_on_face(bc_ptr[BCVars::yvel_bc].lo(0))) {
        strainRate +=  (-(8./3.) * v(i-1,j,k) + 3. * v(i,j,k) - (1./3.) * v(i+1,j,k))*dx_inv;
    } else if (i == domain.bigEnd(0)+1 && is_on_face(bc_ptr[BCVars::yvel_bc].hi(0))) {
        strainRate += -(-(8./3.) * v(i,j,k) + 3. * v(i-1,j,k) - (1./3.) * v(i-2,j,k))*dx_inv;
    } else {
        strainRate += (v(i, j, k) - v(i-1, j, k)) * dx_inv;
    }
    strainRate *= 0.5;

    return strainRate;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ComputeS13 (int i, int j, int k,
            const amrex::Array4<amrex::Real const>& u,
            const amrex::Array4<amrex::Real const>& /*v*/,
            const amrex::Array4<amrex::Real const>& w,
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
            const amrex::Box& domain, const amrex::BCRec* bc_ptr)
{
    amrex::Real dx_inv = cellSizeInv[0];
    amrex::Real dz_inv = cellSizeInv[2];

    amrex::Real strainRate;

    if (k == domain.smallEnd(2) && is_on_face(bc_ptr[BCVars::xvel_bc].lo(2))) {
        strainRate =  (-(8./3.) * u(i,j,k-1) + 3. * u(i,j,k) - (1./3.) * u(i,j,k+1))*dz_inv;
    } else if (k == domain.bigEnd(2)+1 && is_on_face(bc_ptr[BCVars::xvel_bc].hi(2))) {
        strainRate = -(-(8./3.) * u(i,j,k) + 3. * u(i,j,k-1) - (1./3.) * u(i,j,k-2))*dz_inv;
    } else {
        strainRate = (u(i, j, k) - u(i, j, k-1))*dz_inv;
    }
    if (i == domain.smallEnd(0) && is_on_face(bc_ptr[BCVars::zvel_bc].lo(0))) {
        strainRate +=  (-(8./3.) * w(i-1,j,k) + 3. * w(i,j,k) - (1./3.) * w(i+1,j,k))*dx_inv;
    } else if (i == domain.bigEnd(0)+1 && is_on_face(bc_ptr[BCVars::zvel_bc].hi(0))) {
        strainRate += -(-(8./3.) * w(i,j,k) + 3. * w(i-1,j,k) - (1./3.) * w(i-2,j,k))*dx_inv;
    } else {
        if (k > domain.smallEnd(2) && k <= domain.bigEnd(2)+1) {
            strainRate += (w(i, j, k) - w(i-1, j, k)) * dx_inv;
        }
    }
    strainRate *= 0.5;

    return strainRate;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ComputeS23 (int i, int j, int k,
            const amrex::Array4<amrex::Real const>& /*u*/,
            const amrex::Array4<amrex::Real const>& v,
            const amrex::Array4<amrex::Real const>& w,
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
            const amrex::Box& domain, const amrex::BCRec* bc_ptr)
{
    amrex::Real dy_inv = cellSizeInv[1];
    amrex::Real dz_inv = cellSizeInv[2];

    amrex::Real strainRate;

    if (k == domain.smallEnd(2) && is_on_face(bc_ptr[BCVars::yvel_bc].lo(2))) {
        strainRate =  (-(8./3.) * v(i,j,k-1) + 3. * v(i,j,k  ) - (1./3.) * v(i,j,k+1))*dz_inv;
    } else if (k == domain.bigEnd(2)+1 && is_on_face(bc_ptr[BCVars::yvel_bc].hi(2))) {
        strainRate = -(-(8./3.) * v(i,j,k  ) + 3. * v(i,j,k-1) - (1./3.) * v(i,j,k-2))*dz_inv;
    } else {
        strainRate = (v(i, j, k) - v(i, j, k-1))*dz_inv;
    }
    if (j == domain.smallEnd(1) && is_on_face(bc_ptr[BCVars::zvel_bc].lo(1))) {
        strainRate +=  (-(8./3.) * w(i,j-1,k) + 3. * w(i,j  ,k) - (1./3.) * w(i,j+1,k))*dy_inv;
    } else if (j == domain.bigEnd(1)+1 && is_on_face(bc_ptr[BCVars::zvel_bc].hi(1))) {
        strainRate += -(-(8./3.) * w(i,j  ,k) + 3. * w(i,j-1,k) - (1./3.) * w(i,j-2,k))*dy_inv;
    } else {
        if (k > domain.smallEnd(2) && k <= domain.bigEnd(2)+1)
            strainRate += (w(i, j, k) - w(i, j-1, k))*dy_inv;
    }
    strainRate *= 0.5;

    return strainRate;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ComputeS21 (int i, int j, int k,
            const amrex::Array4<amrex::Real const>& u,
            const amrex::Array4<amrex::Real const>& v,
            const amrex::Array4<amrex::Real const>& w,
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
            const amrex::Box& domain, const amrex::BCRec* bc_ptr)
{
    return ComputeS12(i,j,k,u,v,w,cellSizeInv,domain,bc_ptr);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ComputeS31 (int i, int j, int k,
            const amrex::Array4<amrex::Real const>& u,
            const amrex::Array4<amrex::Real const>& v,
            const amrex::Array4<amrex::Real const>& w,
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
            const amrex::Box& domain, const amrex::BCRec* bc_ptr)
{
    return ComputeS13(i,j,k,u,v,w,cellSizeInv,domain,bc_ptr);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ComputeS32 (int i, int j, int k,
            const amrex::Array4<amrex::Real const>& u,
            const amrex::Array4<amrex::Real const>& v,
            const amrex::Array4<amrex::Real const>& w,
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
            const amrex::Box& domain, const amrex::BCRec* bc_ptr)
{
    return ComputeS23(i,j,k,u,v,w,cellSizeInv,domain,bc_ptr);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ComputeStrainRate (int i, int j, int k,
                   const amrex::Array4<amrex::Real const>& u,
                   const amrex::Array4<amrex::Real const>& v,
                   const amrex::Array4<amrex::Real const>& w,
                   enum MomentumEqn  momentumEqn,
                   enum DiffusionDir diffDir,
                   const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                   const amrex::Box& domain, const amrex::BCRec* bc_ptr)
{
    amrex::Real strainRate;
    switch (momentumEqn) {
      case MomentumEqn::x:
        switch (diffDir) {
          case DiffusionDir::x: // S11
            strainRate = ComputeS11(i,j,k,u,cellSizeInv);
            break;
          case DiffusionDir::y: // S12
            strainRate = ComputeS12(i,j,k,u,v,w,cellSizeInv,domain,bc_ptr);
            break;
          case DiffusionDir::z: // S13
             strainRate = ComputeS13(i,j,k,u,v,w,cellSizeInv,domain,bc_ptr);
             break;
          }
          break;
      case MomentumEqn::y:
        switch (diffDir) {
          case DiffusionDir::x: // S21
            strainRate = ComputeS21(i,j,k,u,v,w,cellSizeInv,domain,bc_ptr);
            break;
          case DiffusionDir::y: // S22
            strainRate = ComputeS22(i,j,k,v,cellSizeInv);
            break;
          case DiffusionDir::z: // S23
            strainRate = ComputeS23(i,j,k,u,v,w,cellSizeInv,domain,bc_ptr);
            break;
        }
        break;
      case MomentumEqn::z:
        switch (diffDir) {
          case DiffusionDir::x: // S31
            strainRate = ComputeS31(i,j,k,u,v,w,cellSizeInv,domain,bc_ptr);
            break;
          case DiffusionDir::y: // S32
            strainRate = ComputeS32(i,j,k,u,v,w,cellSizeInv,domain,bc_ptr);
            break;
          case DiffusionDir::z: // S33
            strainRate = ComputeS33(i,j,k,w,cellSizeInv,domain);
            break;
        }
        break;
  }

  return strainRate;
}
#endif
