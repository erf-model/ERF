/** \file EddyViscosity.H */

#ifndef _EDDY_VISCOSITY_H_
#define _EDDY_VISCOSITY_H_

#include <TimeIntegration.H>
#include "StrainRate.H"

using namespace amrex;

/** Compute Eddy Viscosity */
//AMREX_GPU_DEVICE
inline
void ComputeTurbulentViscosity(const MultiFab& xvel, const MultiFab& yvel, const MultiFab& zvel,
                               const MultiFab& cons_in, MultiFab& eddyViscosity,
                               const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                               const SolverChoice& solverChoice,
                               const bool& dirichlet_at_lo_k, const int& klo,
                               const bool& dirichlet_at_hi_k, const int& khi)
{
    const Real cellVol = 1.0 / (cellSizeInv[0] * cellSizeInv[1] * cellSizeInv[2]);
    const Real Delta = std::pow(cellVol,1.0/3.0);

    if (solverChoice.les_type == LESType::Smagorinsky)
    {
      Real Cs = solverChoice.Cs;
      Real CsDeltaSqr = Cs*Cs*Delta*Delta;

      for ( MFIter mfi(eddyViscosity,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box &bx = mfi.growntilebox(1);

        const Array4<Real const > &cell_data = cons_in.array(mfi);
        const Array4<Real> &K = eddyViscosity.array(mfi);

        const Array4<Real const> &u = xvel.array(mfi);
        const Array4<Real const> &v = yvel.array(mfi);
        const Array4<Real const> &w = zvel.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real S11 = ComputeStrainRate(i+1, j, k, u, v, w, MomentumEqn::x, DiffusionDir::x, cellSizeInv, false, false);
            Real S22 = ComputeStrainRate(i, j+1, k, u, v, w, MomentumEqn::y, DiffusionDir::y, cellSizeInv, false, false);
            Real S33 = ComputeStrainRate(i, j, k+1, u, v, w, MomentumEqn::z, DiffusionDir::z, cellSizeInv, false, false);

            Real S12 = 0.25* (
                      ComputeStrainRate(i  , j  , k, u, v, w, MomentumEqn::x, DiffusionDir::y, cellSizeInv, false, false)
                    + ComputeStrainRate(i  , j+1, k, u, v, w, MomentumEqn::x, DiffusionDir::y, cellSizeInv, false, false)
                    + ComputeStrainRate(i+1, j  , k, u, v, w, MomentumEqn::x, DiffusionDir::y, cellSizeInv, false, false)
                    + ComputeStrainRate(i+1, j+1, k, u, v, w, MomentumEqn::x, DiffusionDir::y, cellSizeInv, false, false)
                    );

            bool dirichlet_lo = (dirichlet_at_lo_k && (k == klo));
            bool dirichlet_hi = (dirichlet_at_hi_k && (k == khi));

            Real S13 = 0.25* (
                      ComputeStrainRate(i  , j, k  , u, v, w, MomentumEqn::x, DiffusionDir::z, cellSizeInv,
                                       dirichlet_lo, false)
                    + ComputeStrainRate(i  , j, k+1, u, v, w, MomentumEqn::x, DiffusionDir::z, cellSizeInv,
                                       false, dirichlet_hi)
                    + ComputeStrainRate(i+1, j, k  , u, v, w, MomentumEqn::x, DiffusionDir::z, cellSizeInv,
                                       dirichlet_lo, false)
                    + ComputeStrainRate(i+1, j, k+1, u, v, w, MomentumEqn::x, DiffusionDir::z, cellSizeInv,
                                       false, dirichlet_hi)
                    );

            Real S23 = 0.25* (
                      ComputeStrainRate(i, j  , k  , u, v, w, MomentumEqn::y, DiffusionDir::z, cellSizeInv,
                                       dirichlet_lo, false)
                    + ComputeStrainRate(i, j  , k+1, u, v, w, MomentumEqn::y, DiffusionDir::z, cellSizeInv,
                                       false, dirichlet_hi)
                    + ComputeStrainRate(i, j+1, k  , u, v, w, MomentumEqn::y, DiffusionDir::z, cellSizeInv,
                                       dirichlet_lo, false)
                    + ComputeStrainRate(i, j+1, k+1, u, v, w, MomentumEqn::y, DiffusionDir::z, cellSizeInv,
                                       false, dirichlet_hi)
                    );

            Real SmnSmn = S11*S11 + S22*S22 + S33*S33 + 2.0*S12*S12 + 2.0*S13*S13 + 2.0*S23*S23;
            // Note the positive sign, which aligns well with the positive sign in the diffusion term for momentum equation
            K(i, j, k, 0) = 2.0 * CsDeltaSqr * cell_data(i, j, k, Rho_comp) * std::sqrt(2.0*SmnSmn);
        });

        // This is a complete hack for now
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (dirichlet_at_lo_k && k < klo) K(i,j,k,0) = K(i,j,0,0);
            if (dirichlet_at_hi_k && k > khi) K(i,j,k,0) = K(i,j,khi,0);
        });

      } //mfi
    }
    else if (solverChoice.les_type == LESType::Deardorff)
    {

      Real l_C_k = solverChoice.Ck;

      for ( MFIter mfi(eddyViscosity,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box &bx = mfi.tilebox();

        const Array4<Real const > &cell_data = cons_in.array(mfi);
        const Array4<Real> &K = eddyViscosity.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            K(i,j,k,0) = l_C_k * Delta *  cell_data(i,j,k,Rho_comp) *
                                std::sqrt(cell_data(i,j,k,RhoKE_comp));
        });

      } //mfi
    }
} // function call

/// Compute K_LES (i-1/2, j+1/2, k) etc given K_LES (i, j, k) is known
// Note: This should be at edges for momEqnDir != diffDir, cell centers otherwise
AMREX_GPU_DEVICE
inline Real
InterpolateTurbulentViscosity(const int &i, const int &j, const int &k,
                              const enum MomentumEqn &momentumEqn,
                              const enum DiffusionDir &diffDir,
                              const Array4<Real>& K_LES) {
  // Assuming we already have 'K_LES' computed for all (i, j, k)
  Real turbViscInterpolated = 1.0;

  switch (momentumEqn) {
  case MomentumEqn::x: // Reference face is x-face index (i, j, k)
    switch (diffDir) {
    case DiffusionDir::x:
      turbViscInterpolated = K_LES(i-1, j, k);
      break;
    case DiffusionDir::y:
      turbViscInterpolated = 0.25*( K_LES(i-1, j, k) + K_LES(i, j, k) + K_LES(i-1, j-1, k) + K_LES(i, j-1, k) );
      break;
    case DiffusionDir::z:
      turbViscInterpolated = 0.25*( K_LES(i-1, j, k) + K_LES(i, j, k) + K_LES(i-1, j, k-1) + K_LES(i, j, k-1) );
      break;
    default:
      amrex::Abort("Error: Diffusion direction is unrecognized");
    }
    break;
  case MomentumEqn::y: // Reference face is y-face index (i, j, k)
    switch (diffDir) {
    case DiffusionDir::x:
      turbViscInterpolated = 0.25*( K_LES(i, j-1, k) + K_LES(i, j, k) + K_LES(i-1, j-1, k) + K_LES(i-1, j, k) );
      break;
    case DiffusionDir::y:
      turbViscInterpolated = K_LES(i, j-1, k);
      break;
    case DiffusionDir::z:
      turbViscInterpolated = 0.25*( K_LES(i, j-1, k) + K_LES(i, j, k) + K_LES(i, j-1, k-1) + K_LES(i, j, k-1) );
      break;
    default:
      amrex::Abort("Error: Diffusion direction is unrecognized");
    }
    break;
  case MomentumEqn::z: // Reference face is z-face index (i, j, k)
    switch (diffDir) {
    case DiffusionDir::x:
      turbViscInterpolated = 0.25*( K_LES(i, j, k-1) + K_LES(i, j, k) + K_LES(i-1, j, k-1) + K_LES(i-1, j, k) );
      break;
    case DiffusionDir::y:
      turbViscInterpolated = 0.25*( K_LES(i, j, k-1) + K_LES(i, j, k) + K_LES(i, j-1, k-1) + K_LES(i, j-1, k) );
      break;
    case DiffusionDir::z:
      turbViscInterpolated = K_LES(i, j, k-1);
      break;
    default:
      amrex::Abort("Error: Diffusion direction is unrecognized");
    }
    break;
  default:
    amrex::Abort("Error: Momentum equation is unrecognized");
  }

  return turbViscInterpolated;
}

#endif
