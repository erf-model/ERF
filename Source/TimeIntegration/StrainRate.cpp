#include <TimeIntegration.H>

using namespace amrex;

AMREX_GPU_DEVICE
Real
ComputeStrainRate(const int &i, const int &j, const int &k,
                  const Array4<Real const>& u, const Array4<Real const>& v, const Array4<Real const>& w,
                  const enum NextOrPrev &nextOrPrev,
                  const enum MomentumEqn &momentumEqn,
                  const enum DiffusionDir &diffDir,
                  const GpuArray<Real, AMREX_SPACEDIM>& cellSize,
                  bool use_no_slip_stencil) {
  Real dx_inv = 1.0/cellSize[0];
  Real dy_inv = 1.0/cellSize[1];
  Real dz_inv = 1.0/cellSize[2];

  Real strainRate = 0;

  switch (momentumEqn) {
  case MomentumEqn::x:
    switch (diffDir) {
    case DiffusionDir::x: // S11
      if (nextOrPrev == NextOrPrev::next)
        strainRate = (u(i+1, j, k) - u(i, j, k))*dx_inv; // S11 (i+1/2)
      else // nextOrPrev == NextOrPrev::prev
        strainRate = (u(i, j, k) - u(i-1, j, k))*dx_inv; // S11 (i-1/2)
      break;
    case DiffusionDir::y: // S12
      if (nextOrPrev == NextOrPrev::next)
        strainRate = (u(i, j+1, k) - u(i, j, k))*dy_inv + (v(i, j+1, k) - v(i-1, j+1, k))*dx_inv; // S12 (j+1/2)
      else // nextOrPrev == NextOrPrev::prev
        strainRate = (u(i, j, k) - u(i, j-1, k))*dy_inv + (v(i, j, k) - v(i-1, j, k))*dx_inv; // S12 (j-1/2)
      strainRate *= 0.5;
      break;
    case DiffusionDir::z: // S13
      if (nextOrPrev == NextOrPrev::next)
      {
        if (use_no_slip_stencil) {
            strainRate =  -(3. * u(i,j,k) - (1./3.) * u(i,j,k-1))*dz_inv
                         + (w(i, j, k+1) - w(i-1, j, k+1))*dx_inv; // S13 (k-1/2); // S13 (k+1/2)
        } else {
            strainRate = (u(i, j, k+1) - u(i, j, k))*dz_inv + (w(i, j, k+1) - w(i-1, j, k+1))*dx_inv; // S13 (k+1/2)
        }
      }
      else // nextOrPrev == NextOrPrev::prev
      {
        if (use_no_slip_stencil) {
            strainRate =  (3. * u(i,j,k) - (1./3.) * u(i,j,k+1))*dz_inv
                        + (w(i, j, k) - w(i-1, j, k))*dx_inv; // S13 (k-1/2); // S13 (k-1/2)
        } else {
            strainRate = (u(i, j, k) - u(i, j, k-1))*dz_inv + (w(i, j, k) - w(i-1, j, k))*dx_inv; // S13 (k-1/2)
        }
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
      if (nextOrPrev == NextOrPrev::next)
        strainRate = (u(i+1, j, k) - u(i+1, j-1, k))*dy_inv + (v(i+1, j, k) - v(i, j, k))*dx_inv; // S21 (i+1/2)
      else // nextOrPrev == NextOrPrev::prev
        strainRate = (u(i, j, k) - u(i, j-1, k))*dy_inv + (v(i, j, k) - v(i-1, j, k))*dx_inv; // S21 (i-1/2)
      strainRate *= 0.5;
      break;
    case DiffusionDir::y: // S22
      if (nextOrPrev == NextOrPrev::next)
        strainRate = (v(i, j+1, k) - v(i, j, k))*dy_inv; // S22 (j+1/2)
      else // nextOrPrev == NextOrPrev::prev
        strainRate = (v(i, j, k) - v(i, j-1, k))*dy_inv; // S22 (j-1/2)
      break;
    case DiffusionDir::z: // S23
      if (nextOrPrev == NextOrPrev::next)
      {
        if (use_no_slip_stencil) {
            strainRate =  -(3. * v(i,j,k) - (1./3.) * v(i,j,k-1))*dz_inv
                         + (w(i, j, k+1) - w(i, j-1, k+1))*dy_inv; // S23 (k+1/2
        } else {
            strainRate = (v(i, j, k+1) - v(i, j, k))*dz_inv
                       + (w(i, j, k+1) - w(i, j-1, k+1))*dy_inv; // S23 (k+1/2)
        }
      }
      else // nextOrPrev == NextOrPrev::prev
      {
        if (use_no_slip_stencil) {
            strainRate =  (3. * v(i,j,k) - (1./3.) * v(i,j,k+1))*dz_inv
                        + (w(i, j, k) - w(i, j-1, k))*dy_inv; // S23 (k-1/2)
        } else {
            strainRate = (v(i, j, k) - v(i, j, k-1))*dz_inv
                       + (w(i, j, k) - w(i, j-1, k))*dy_inv; // S23 (k-1/2)
        }
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
      if (nextOrPrev == NextOrPrev::next)
        strainRate = (u(i+1, j, k) - u(i+1, j, k-1))*dz_inv + (w(i+1, j, k) - w(i, j, k))*dx_inv; // S31 (i+1/2)
      else // nextOrPrev == NextOrPrev::prev
        strainRate = (u(i, j, k) - u(i, j, k-1))*dz_inv + (w(i, j, k) - w(i-1, j, k))*dx_inv; // S31 (i-1/2)
      strainRate *= 0.5;
      break;
    case DiffusionDir::y: // S32
      if (nextOrPrev == NextOrPrev::next)
        strainRate = (v(i, j+1, k) - v(i, j+1, k-1))*dz_inv + (w(i, j+1, k) - w(i, j, k))*dy_inv; // S32 (j+1/2)
      else // nextOrPrev == NextOrPrev::prev
        strainRate = (v(i, j, k) - v(i, j, k-1))*dz_inv + (w(i, j, k) - w(i, j-1, k))*dy_inv; // S32 (j-1/2)
      strainRate *= 0.5;
      break;
    case DiffusionDir::z: // S33
      if (nextOrPrev == NextOrPrev::next)
        strainRate = (w(i, j, k+1) - w(i, j, k))*dz_inv; // S33 (k+1/2)
      else // nextOrPrev == NextOrPrev::prev
        strainRate = (w(i, j, k) - w(i, j, k-1))*dz_inv; // S33 (k-1/2)
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

AMREX_GPU_DEVICE
Real
ComputeExpansionRate(const int &i, const int &j, const int &k,
                     const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                     const enum NextOrPrev &nextOrPrev,
                     const enum MomentumEqn &momentumEqn,
                     const enum DiffusionDir &diffDir,
                     const GpuArray<Real, AMREX_SPACEDIM>& cellSize) {
    Real dx_inv = 1.0/cellSize[0];
    Real dy_inv = 1.0/cellSize[1];
    Real dz_inv = 1.0/cellSize[2];

    Real expansionRate = 0;

    //TODO: Check if it is more efficient to store computed divergence at cell-centers rather than
    // computing on the fly

    switch (momentumEqn) {
        case MomentumEqn::x:
            switch (diffDir) {
                case DiffusionDir::x: // D11
                    if (nextOrPrev == NextOrPrev::next) // cell (i, j, k) is at x-face (i+1/2, j, k)
                        // D11 (i+1/2)
                        expansionRate = (u(i+1, j, k) - u(i, j, k))*dx_inv +
                                        (v(i, j+1, k) - v(i, j, k))*dy_inv +
                                        (w(i, j, k+1) - w(i, j, k))*dz_inv;
                    else // nextOrPrev == NextOrPrev::prev // D11 (i-1/2)
                        expansionRate = (u(i, j, k) - u(i-1, j, k))*dx_inv +
                                        (v(i-1, j+1, k) - v(i-1, j, k))*dy_inv +
                                        (w(i-1, j, k+1) - w(i-1, j, k))*dz_inv;
                    break;
                case DiffusionDir::y: // D12
                    expansionRate = 0.0;
                    break;
                case DiffusionDir::z: // D13
                    expansionRate = 0.0;
                    break;
                default:
                    amrex::Abort("Error: Diffusion direction is unrecognized");
            }
            break;
        case MomentumEqn::y:
            switch (diffDir) {
                case DiffusionDir::x: // D21
                    expansionRate = 0.0;
                    break;
                case DiffusionDir::y: // D22
                    if (nextOrPrev == NextOrPrev::next) // cell (i, j, k) is at y-face (i, j+1/2, k)
                        // D22 (j+1/2)
                        expansionRate = (u(i+1, j, k) - u(i, j, k))*dx_inv +
                                        (v(i, j+1, k) - v(i, j, k))*dy_inv +
                                        (w(i, j, k+1) - w(i, j, k))*dz_inv;
                    else // nextOrPrev == NextOrPrev::prev // D22 (j-1/2)
                        expansionRate = (u(i+1, j-1, k) - u(i, j-1, k))*dx_inv +
                                        (v(i, j, k) - v(i, j-1, k))*dy_inv +
                                        (w(i, j-1, k+1) - w(i, j-1, k))*dz_inv;
                    break;
                case DiffusionDir::z: // D23
                    expansionRate = 0.0;
                    break;
                default:
                    amrex::Abort("Error: Diffusion direction is unrecognized");
            }
            break;
        case MomentumEqn::z:
            switch (diffDir) {
                case DiffusionDir::x: // D31
                    expansionRate = 0.0;
                    break;
                case DiffusionDir::y: // D32
                    expansionRate = 0.0;
                    break;
                case DiffusionDir::z: // D33
                    if (nextOrPrev == NextOrPrev::next) // cell (i, j, k) is at z-face (i, j, k+1/2)
                        // D33 (k+1/2)
                        expansionRate = (u(i+1, j, k) - u(i, j, k))*dx_inv +
                                        (v(i, j+1, k) - v(i, j, k))*dy_inv +
                                        (w(i, j, k+1) - w(i, j, k))*dz_inv;
                    else // nextOrPrev == NextOrPrev::prev // D33 (k-1/2)
                        expansionRate = (u(i+1, j, k-1) - u(i, j, k-1))*dx_inv +
                                        (v(i, j+1, k-1) - v(i, j, k-1))*dy_inv +
                                        (w(i, j, k) - w(i, j, k-1))*dz_inv;
                    break;
                default:
                    amrex::Abort("Error: Diffusion direction is unrecognized");
            }
            break;
        default:
            amrex::Abort("Error: Momentum equation is unrecognized");
    }

    return (1.0/3.0) * expansionRate;
}
