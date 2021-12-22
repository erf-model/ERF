#ifndef _COMPUTE_STRAINRATE_H_
#define _COMPUTE_STRAINRATE_H_

using namespace amrex;

AMREX_GPU_DEVICE
inline Real
ComputeStrainRate(const int &i, const int &j, const int &k,
                  const Array4<Real const>& u, const Array4<Real const>& v, const Array4<Real const>& w,
                  const enum MomentumEqn &momentumEqn,
                  const enum DiffusionDir &diffDir,
                  const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                  bool dirichlet_lo, bool dirichlet_hi)
{

  // Note:  if dirichlet_lo/dirichlet_hi is true, it means that the velocity is fixed
  //        on the z-face, so we use a special stencil that understands the value is
  //        at the z-face, not at dz/2 above/below the face.   The stencils below also
  //        allow for the Dirichlet value to be non-zero

  Real dx_inv = cellSizeInv[0];
  Real dy_inv = cellSizeInv[1];
  Real dz_inv = cellSizeInv[2];

  Real strainRate = 0;

  switch (momentumEqn) {
  case MomentumEqn::x:
    switch (diffDir) {
    case DiffusionDir::x: // S11
      strainRate = (u(i, j, k) - u(i-1, j, k))*dx_inv;
      break;
    case DiffusionDir::y: // S12
      strainRate  = (u(i, j, k) - u(i, j-1, k))*dy_inv + (v(i, j, k) - v(i-1, j, k)) * dx_inv;
      strainRate *= 0.5;
      break;
    case DiffusionDir::z: // S13
      if (dirichlet_lo) {
          strainRate =   (-(8./3.) * u(i,j,k-1) + 3. * u(i,j,k) - (1./3.) * u(i,j,k+1))*dz_inv
                       + (w(i, j, k) - w(i-1, j, k))*dx_inv;
      } else if (dirichlet_hi) {
          strainRate =  -(-(8./3.) * u(i,j,k) + 3. * u(i,j,k-1) - (1./3.) * u(i,j,k-2))*dz_inv
                       + (w(i, j, k-1) - w(i-1, j, k-1))*dx_inv;
      } else {
          strainRate = (u(i, j, k) - u(i, j, k-1))*dz_inv + (w(i, j, k) - w(i-1, j, k))*dx_inv;
      }
      strainRate *= 0.5;
      break;
    default:
      amrex::Abort("Error: Diffusion direction is unrecognized");
    }
    break;
  case MomentumEqn::y:
    switch (diffDir) {
    case DiffusionDir::x: // S21
      strainRate = (u(i, j, k) - u(i, j-1, k))*dy_inv + (v(i, j, k) - v(i-1, j, k)) * dx_inv;
      strainRate *= 0.5;
      break;
    case DiffusionDir::y: // S22
      strainRate = (v(i, j, k) - v(i, j-1, k))*dy_inv;
      break;
    case DiffusionDir::z: // S23
      if (dirichlet_lo) {
          strainRate =   (-(8./3.) * v(i,j,k-1) + 3. * v(i,j,k) - (1./3.) * v(i,j,k+1))*dz_inv
                       + (w(i, j, k) - w(i, j-1, k))*dy_inv;
      } else if (dirichlet_hi) {
          strainRate =  -(-(8./3.) * v(i,j,k) + 3. * v(i,j,k-1) - (1./3.) * v(i,j,k-2))*dz_inv
                       + (w(i, j, k) - w(i, j-1, k))*dy_inv;
      } else {
          strainRate = (v(i, j, k) - v(i, j, k-1))*dz_inv + (w(i, j, k) - w(i, j-1, k))*dy_inv;
      }
      strainRate *= 0.5;
      break;
    default:
      amrex::Abort("Error: Diffusion direction is unrecognized");
    }
    break;
  case MomentumEqn::z:
    switch (diffDir) {
    case DiffusionDir::x: // S31
      strainRate = (u(i, j, k) - u(i, j, k-1))*dz_inv + (w(i, j, k) - w(i-1, j, k)) * dx_inv;
      strainRate *= 0.5;
      break;
    case DiffusionDir::y: // S32
      strainRate = (v(i, j, k) - v(i, j, k-1))*dz_inv + (w(i, j, k) - w(i, j-1, k)) * dy_inv;
      strainRate *= 0.5;
      break;
    case DiffusionDir::z: // S33
      strainRate = (w(i, j, k) - w(i, j, k-1))*dz_inv;
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
