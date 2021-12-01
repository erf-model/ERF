/** \file EddyViscosity.cpp */
#include <TimeIntegration.H>

using namespace amrex;

/** Compute Eddy Viscosity */
//AMREX_GPU_DEVICE
void ComputeTurbulentViscosity(const MultiFab& xvel, const MultiFab& yvel, const MultiFab& zvel,
                               const MultiFab& cons_in, MultiFab& eddyViscosity,
                               const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                               const SolverChoice& solverChoice,
                               const bool& use_no_slip_stencil_at_lo_k, const int& klo,
                               const bool& use_no_slip_stencil_at_hi_k, const int& khi)
{
    const Real cellVol = 1.0 / (cellSizeInv[0] * cellSizeInv[1] * cellSizeInv[2]);
    Real Cs = solverChoice.Cs;
    Real CsDeltaSqr = Cs*Cs * std::pow(cellVol, 2.0/3.0);

    for ( MFIter mfi(eddyViscosity,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box &bx = mfi.tilebox();

        const Array4<Real const > &cell_data = cons_in.array(mfi);
        const Array4<Real> &K = eddyViscosity.array(mfi);

        const Array4<Real const> &u = xvel.array(mfi);
        const Array4<Real const> &v = yvel.array(mfi);
        const Array4<Real const> &w = zvel.array(mfi);

        amrex::ParallelFor(bx, eddyViscosity.nComp(),[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
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

            bool use_no_slip_stencil_lo = (use_no_slip_stencil_at_lo_k && (k == klo));
            bool use_no_slip_stencil_hi = (use_no_slip_stencil_at_hi_k && (k == khi));

            Real S13 = 0.25* (
                      ComputeStrainRate(i  , j, k  , u, v, w, MomentumEqn::x, DiffusionDir::z, cellSizeInv,
                                       use_no_slip_stencil_lo, use_no_slip_stencil_hi)
                    + ComputeStrainRate(i  , j, k+1, u, v, w, MomentumEqn::x, DiffusionDir::z, cellSizeInv,
                                       use_no_slip_stencil_lo, use_no_slip_stencil_hi)
                    + ComputeStrainRate(i+1, j, k  , u, v, w, MomentumEqn::x, DiffusionDir::z, cellSizeInv,
                                       use_no_slip_stencil_lo, use_no_slip_stencil_hi)
                    + ComputeStrainRate(i+1, j, k+1, u, v, w, MomentumEqn::x, DiffusionDir::z, cellSizeInv,
                                       use_no_slip_stencil_lo, use_no_slip_stencil_hi)
                    );

            Real S23 = 0.25* (
                      ComputeStrainRate(i, j  , k  , u, v, w, MomentumEqn::y, DiffusionDir::z, cellSizeInv,
                                       use_no_slip_stencil_lo, use_no_slip_stencil_hi)
                    + ComputeStrainRate(i, j  , k+1, u, v, w, MomentumEqn::y, DiffusionDir::z, cellSizeInv,
                                       use_no_slip_stencil_lo, use_no_slip_stencil_hi)
                    + ComputeStrainRate(i, j+1, k  , u, v, w, MomentumEqn::y, DiffusionDir::z, cellSizeInv,
                                       use_no_slip_stencil_lo, use_no_slip_stencil_hi)
                    + ComputeStrainRate(i, j+1, k+1, u, v, w, MomentumEqn::y, DiffusionDir::z, cellSizeInv,
                                       use_no_slip_stencil_lo, use_no_slip_stencil_hi)
                    );

            Real SmnSmn = S11*S11 + S22*S22 + S33*S33 + 2.0*S12*S12 + 2.0*S13*S13 + 2.0*S23*S23;
            // Note the positive sign, which aligns well with the positive sign in the diffusion term for momentum equation
            K(i, j, k, 0) = 2.0 * CsDeltaSqr * cell_data(i, j, k, Rho_comp) * std::sqrt(2.0*SmnSmn);
        });

    } //mfi
} // function call

/// Compute Ksmag (i-1/2, j+1/2, k) etc given Ksmag (i, j, k) is known
// Note: This should be at edges for momEqnDir != diffDir, cell centers otherwise
AMREX_GPU_DEVICE
Real
InterpolateTurbulentViscosity(const int &i, const int &j, const int &k,
                              const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                              const enum MomentumEqn &momentumEqn,
                              const enum DiffusionDir &diffDir,
                              const Array4<Real>& Ksmag) {
  // Assuming we already have 'Ksmag' computed for all (i, j, k)
  Real turbViscInterpolated = 1.0;

  switch (momentumEqn) {
  case MomentumEqn::x: // Reference face is x-face index (i, j, k)
    switch (diffDir) {
    case DiffusionDir::x:
      turbViscInterpolated = Ksmag(i-1, j, k);
      break;
    case DiffusionDir::y:
      turbViscInterpolated = 0.25*( Ksmag(i-1, j, k) + Ksmag(i, j, k) + Ksmag(i-1, j-1, k) + Ksmag(i, j-1, k) );
      break;
    case DiffusionDir::z:
      turbViscInterpolated = 0.25*( Ksmag(i-1, j, k) + Ksmag(i, j, k) + Ksmag(i-1, j, k-1) + Ksmag(i, j, k-1) );
      break;
    default:
      amrex::Abort("Error: Diffusion direction is unrecognized");
    }
    break;
  case MomentumEqn::y: // Reference face is y-face index (i, j, k)
    switch (diffDir) {
    case DiffusionDir::x:
      turbViscInterpolated = 0.25*( Ksmag(i, j-1, k) + Ksmag(i, j, k) + Ksmag(i-1, j-1, k) + Ksmag(i-1, j, k) );
      break;
    case DiffusionDir::y:
      turbViscInterpolated = Ksmag(i, j-1, k);
      break;
    case DiffusionDir::z:
      turbViscInterpolated = 0.25*( Ksmag(i, j-1, k) + Ksmag(i, j, k) + Ksmag(i, j-1, k-1) + Ksmag(i, j, k-1) );
      break;
    default:
      amrex::Abort("Error: Diffusion direction is unrecognized");
    }
    break;
  case MomentumEqn::z: // Reference face is z-face index (i, j, k)
    switch (diffDir) {
    case DiffusionDir::x:
      turbViscInterpolated = 0.25*( Ksmag(i, j, k-1) + Ksmag(i, j, k) + Ksmag(i-1, j, k-1) + Ksmag(i-1, j, k) );
      break;
    case DiffusionDir::y:
      turbViscInterpolated = 0.25*( Ksmag(i, j, k-1) + Ksmag(i, j, k) + Ksmag(i, j-1, k-1) + Ksmag(i, j-1, k) );
      break;
    case DiffusionDir::z:
      turbViscInterpolated = Ksmag(i, j, k-1);
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
