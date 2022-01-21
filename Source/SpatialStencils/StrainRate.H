#ifndef _COMPUTE_STRAINRATE_H_
#define _COMPUTE_STRAINRATE_H_

using namespace amrex;

// Note:  if the boundary condition type is ext_dir, it means that the velocity is fixed
//        on that face, so we use a special stencil that understands the value is
//        at the face, not at half a grid cell away from the face.   The stencils below also
//        allow for the Dirichlet value to be non-zero.  A wall boundary condition will generate
//        Dirichlet bc in the z-direction; inflow will generate Dirichlet bc in the lateral directions.

AMREX_GPU_DEVICE
inline Real
ComputeS11(const int &i, const int &j, const int &k,
           const Array4<Real const>& u, const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv)
{
  Real dx_inv = cellSizeInv[0];
  return (u(i, j, k) - u(i-1, j, k))*dx_inv;
}

AMREX_GPU_DEVICE
inline Real
ComputeS22(const int &i, const int &j, const int &k,
           const Array4<Real const>& v, const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv)
{
  Real dy_inv = cellSizeInv[1];
  return (v(i, j, k) - v(i, j-1, k))*dy_inv;
}

AMREX_GPU_DEVICE
inline Real
ComputeS33(const int &i, const int &j, const int &k,
           const Array4<Real const>& w, const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv)
{
  Real dz_inv = cellSizeInv[2];
  return (w(i, j, k) - w(i, j, k-1))*dz_inv;
}

AMREX_GPU_DEVICE
inline Real
ComputeS12(const int &i, const int &j, const int &k,
           const Array4<Real const>& u, const Array4<Real const>& v, const Array4<Real const>& /*w*/,
           const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
           const Box& domain, const amrex::BCRec* bc_ptr)
{
    Real dx_inv = cellSizeInv[0];
    Real dy_inv = cellSizeInv[1];
    Real strainRate = 0;

    if (j == domain.smallEnd(1) && bc_ptr[BCVars::xvel_bc].lo(1) == ERFBCType::ext_dir) {
        strainRate =  (-(8./3.) * u(i,j-1,k) + 3. * u(i,j,k) - (1./3.) * u(i,j+1,k))*dy_inv;
    } else if (j == domain.bigEnd(1)+1 && bc_ptr[BCVars::xvel_bc].hi(1) == ERFBCType::ext_dir) {
        strainRate = -(-(8./3.) * u(i,j,k) + 3. * u(i,j-1,k) - (1./3.) * u(i,j-2,k))*dy_inv;
    } else {
        strainRate  = (u(i, j, k) - u(i, j-1, k))*dy_inv;
    }
    if (i == domain.smallEnd(0) && bc_ptr[BCVars::yvel_bc].lo(0) == ERFBCType::ext_dir) {
        strainRate =  (-(8./3.) * v(i-1,j,k) + 3. * v(i,j,k) - (1./3.) * v(i+1,j,k))*dx_inv;
    } else if (i == domain.bigEnd(0)+1 && bc_ptr[BCVars::yvel_bc].hi(0) == ERFBCType::ext_dir) {
        strainRate = -(-(8./3.) * v(i,j,k) + 3. * v(i-1,j,k) - (1./3.) * v(i-2,j,k))*dx_inv;
    } else {
        strainRate += (v(i, j, k) - v(i-1, j, k)) * dx_inv;
    }
    strainRate *= 0.5;

    return strainRate;
}

AMREX_GPU_DEVICE
inline Real
ComputeS13(const int &i, const int &j, const int &k,
           const Array4<Real const>& u, const Array4<Real const>& /*v*/, const Array4<Real const>& w,
           const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
           const Box& domain, const amrex::BCRec* bc_ptr)
{
    Real dx_inv = cellSizeInv[0];
    Real dz_inv = cellSizeInv[2];

    Real strainRate;

    if (k == domain.smallEnd(2) && bc_ptr[BCVars::xvel_bc].lo(2) == ERFBCType::ext_dir) {
        strainRate =  (-(8./3.) * u(i,j,k-1) + 3. * u(i,j,k) - (1./3.) * u(i,j,k+1))*dz_inv;
    } else if (k == domain.bigEnd(2)+1 && bc_ptr[BCVars::xvel_bc].hi(2) == ERFBCType::ext_dir) {
        strainRate = -(-(8./3.) * u(i,j,k) + 3. * u(i,j,k-1) - (1./3.) * u(i,j,k-2))*dz_inv;
    } else {
        strainRate = (u(i, j, k) - u(i, j, k-1))*dz_inv;
    }
    if (i == domain.smallEnd(0) && bc_ptr[BCVars::zvel_bc].lo(0) == ERFBCType::ext_dir) {
        strainRate =  (-(8./3.) * w(i-1,j,k) + 3. * w(i,j,k) - (1./3.) * w(i+1,j,k))*dx_inv;
    } else if (i == domain.bigEnd(0)+1 && bc_ptr[BCVars::zvel_bc].hi(0) == ERFBCType::ext_dir) {
        strainRate = -(-(8./3.) * w(i,j,k) + 3. * w(i-1,j,k) - (1./3.) * w(i-2,j,k))*dx_inv;
    } else {
        strainRate += (w(i, j, k) - w(i-1, j, k)) * dx_inv;
    }
    strainRate *= 0.5;

    return strainRate;
}

AMREX_GPU_DEVICE
inline Real
ComputeS23(const int &i, const int &j, const int &k,
           const Array4<Real const>& /*u*/, const Array4<Real const>& v, const Array4<Real const>& w,
           const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
           const Box& domain, const amrex::BCRec* bc_ptr)
{
    Real dy_inv = cellSizeInv[1];
    Real dz_inv = cellSizeInv[2];

    Real strainRate;

    if (k == domain.smallEnd(2) && bc_ptr[BCVars::yvel_bc].lo(2) == ERFBCType::ext_dir) {
        strainRate =  (-(8./3.) * v(i,j,k-1) + 3. * v(i,j,k  ) - (1./3.) * v(i,j,k+1))*dz_inv;
    } else if (k == domain.bigEnd(2)+1 && bc_ptr[BCVars::yvel_bc].hi(2) == ERFBCType::ext_dir) {
        strainRate = -(-(8./3.) * v(i,j,k  ) + 3. * v(i,j,k-1) - (1./3.) * v(i,j,k-2))*dz_inv;
    } else {
        strainRate = (v(i, j, k) - v(i, j, k-1))*dz_inv;
    }
    if (j == domain.smallEnd(1) && bc_ptr[BCVars::zvel_bc].lo(1) == ERFBCType::ext_dir) {
        strainRate =  (-(8./3.) * w(i,j-1,k) + 3. * w(i,j  ,k) - (1./3.) * w(i,j+1,k))*dy_inv;
    } else if (j == domain.bigEnd(1)+1 && bc_ptr[BCVars::zvel_bc].hi(1) == ERFBCType::ext_dir) {
        strainRate = -(-(8./3.) * w(i,j  ,k) + 3. * w(i,j-1,k) - (1./3.) * w(i,j-2,k))*dy_inv;
    } else {
        strainRate += (w(i, j, k) - w(i, j-1, k))*dy_inv;
    }
    strainRate *= 0.5;

    return strainRate;
}

AMREX_GPU_DEVICE
inline Real
ComputeS21(const int &i, const int &j, const int &k,
           const Array4<Real const>& u, const Array4<Real const>& v, const Array4<Real const>& w,
           const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
           const Box& domain, const amrex::BCRec* bc_ptr)
{
    return ComputeS12(i,j,k,u,v,w,cellSizeInv,domain,bc_ptr);
}

AMREX_GPU_DEVICE
inline Real
ComputeS31(const int &i, const int &j, const int &k,
           const Array4<Real const>& u, const Array4<Real const>& v, const Array4<Real const>& w,
           const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
           const Box& domain, const amrex::BCRec* bc_ptr)
{
    return ComputeS13(i,j,k,u,v,w,cellSizeInv,domain,bc_ptr);
}

AMREX_GPU_DEVICE
inline Real
ComputeS32(const int &i, const int &j, const int &k,
           const Array4<Real const>& u, const Array4<Real const>& v, const Array4<Real const>& w,
           const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
           const Box& domain, const amrex::BCRec* bc_ptr)
{
    return ComputeS23(i,j,k,u,v,w,cellSizeInv,domain,bc_ptr);
}

AMREX_GPU_DEVICE
inline Real
ComputeStrainRate(const int &i, const int &j, const int &k,
                  const Array4<Real const>& u, const Array4<Real const>& v, const Array4<Real const>& w,
                  const enum MomentumEqn &momentumEqn,
                  const enum DiffusionDir &diffDir,
                  const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                  const Box& domain, const amrex::BCRec* bc_ptr)
{
    Real strainRate;
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
          default:
             amrex::Abort("Error: Diffusion direction is unrecognized");
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
          default:
            amrex::Abort("Error: Diffusion direction is unrecognized");
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
            strainRate = ComputeS33(i,j,k,w,cellSizeInv);
            break;
          default:
            amrex::Abort("Error: Diffusion direction is unrecognized");
        }
        break;
      default:
          amrex::Abort("Error: Momentum equation is unrecognized");
  }

  return strainRate;
}

#endif
