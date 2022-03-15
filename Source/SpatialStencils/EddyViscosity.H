/** \file EddyViscosity.H */

#ifndef _EDDY_VISCOSITY_H_
#define _EDDY_VISCOSITY_H_

#include <TimeIntegration.H>
#include "StrainRate.H"

using namespace amrex;

AMREX_GPU_DEVICE
inline
Real ComputeSmnSmn(const int& i,const int& j,const int& k,
                   const Array4<Real const>& u,
                   const Array4<Real const>& v,
                   const Array4<Real const>& w,
                   const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                   const Box& domain, const amrex::BCRec* bc_ptr)
{
    Real S11 = ComputeS11(i+1, j, k, u, cellSizeInv);
    Real S22 = ComputeS22(i, j+1, k, v, cellSizeInv);
    Real S33 = ComputeS33(i, j, k+1, w, cellSizeInv);

    Real S12 = 0.25* (
              ComputeS12(i  , j  , k, u, v, w, cellSizeInv, domain, bc_ptr)
            + ComputeS12(i  , j+1, k, u, v, w, cellSizeInv, domain, bc_ptr)
            + ComputeS12(i+1, j  , k, u, v, w, cellSizeInv, domain, bc_ptr)
            + ComputeS12(i+1, j+1, k, u, v, w, cellSizeInv, domain, bc_ptr)
            );

    Real S13 = 0.25* (
              ComputeS13(i  , j, k  , u, v, w,  cellSizeInv, domain, bc_ptr)
            + ComputeS13(i  , j, k+1, u, v, w,  cellSizeInv, domain, bc_ptr)
            + ComputeS13(i+1, j, k  , u, v, w,  cellSizeInv, domain, bc_ptr)
            + ComputeS13(i+1, j, k+1, u, v, w,  cellSizeInv, domain, bc_ptr)
            );

    Real S23 = 0.25* (
              ComputeS23(i, j  , k  , u, v, w,  cellSizeInv, domain, bc_ptr)
            + ComputeS23(i, j  , k+1, u, v, w,  cellSizeInv, domain, bc_ptr)
            + ComputeS23(i, j+1, k  , u, v, w,  cellSizeInv, domain, bc_ptr)
            + ComputeS23(i, j+1, k+1, u, v, w,  cellSizeInv, domain, bc_ptr)
            );

    Real SmnSmn = S11*S11 + S22*S22 + S33*S33 + 2.0*S12*S12 + 2.0*S13*S13 + 2.0*S23*S23;

    return SmnSmn;
}

/** Compute Eddy Viscosity */
//AMREX_GPU_DEVICE
inline
void ComputeTurbulentViscosity(const MultiFab& xvel, const MultiFab& yvel, const MultiFab& zvel,
                               const MultiFab& cons_in, MultiFab& eddyViscosity,
                               const Geometry& geom,
                               const SolverChoice& solverChoice,
                               const Gpu::DeviceVector<amrex::BCRec> domain_bcs_type_d)
{
    const GpuArray<Real, AMREX_SPACEDIM> cellSizeInv = geom.InvCellSizeArray();

    const Real cellVol = 1.0 / (cellSizeInv[0] * cellSizeInv[1] * cellSizeInv[2]);
    const Real Delta = std::pow(cellVol,1.0/3.0);

    const auto& domain = geom.Domain();
    const auto& dom_lo = amrex::lbound(domain);
    const auto& dom_hi = amrex::ubound(domain);

    // Compute the turbulent viscosity
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

        const amrex::BCRec* bc_ptr = domain_bcs_type_d.data();

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real SmnSmn = ComputeSmnSmn(i,j,k,u,v,w,cellSizeInv,domain,bc_ptr);

            // Note the positive sign, which aligns well with the positive sign in the diffusion term for momentum equation
            K(i, j, k, EddyDiff::Mom_h) = 2.0 * CsDeltaSqr * cell_data(i, j, k, Rho_comp) * std::sqrt(2.0*SmnSmn);
        });
      } //mfi
    }
    else if (solverChoice.les_type == LESType::Deardorff)
    {
      Real l_C_k = solverChoice.Ck;

      for ( MFIter mfi(eddyViscosity,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box &bx = mfi.growntilebox(1);

        const Array4<Real const > &cell_data = cons_in.array(mfi);
        const Array4<Real> &K = eddyViscosity.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            K(i,j,k,EddyDiff::Mom_h) = l_C_k * Delta *  cell_data(i,j,k,Rho_comp) *
                                std::sqrt(cell_data(i,j,k,RhoKE_comp));
        });
      } //mfi
    }

    // Fill in the rest of the turbulent diffusivities
    if (solverChoice.les_type == LESType::Smagorinsky || solverChoice.les_type == LESType::Deardorff) {
      Real inv_Pr_t = solverChoice.Pr_t_inv;
      Real inv_Sc_t = solverChoice.Sc_t_inv;
      Real inv_sigma_k = 1.0 / solverChoice.sigma_k;
      bool use_KE = (solverChoice.les_type == LESType::Deardorff);
      for ( MFIter mfi(eddyViscosity,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box &bx = mfi.growntilebox(1);
        const Array4<Real> &K = eddyViscosity.array(mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
      // Get eddy diffusivities from the eddy viscosity
      // Additional factor of 0.5 because K = 2 mu_t
      K(i,j,k,EddyDiff::Theta_v) = 0.5 * K(i,j,k,EddyDiff::Mom_h) * inv_Pr_t;
      K(i,j,k,EddyDiff::Scalar_v) = 0.5 * K(i,j,k,EddyDiff::Mom_h) * inv_Sc_t;
      K(i,j,k,EddyDiff::KE_v) = 0.5 * K(i,j,k,EddyDiff::Mom_h) * inv_sigma_k;

      // For LES: vertical and horizontal components are the same
      K(i,j,k,EddyDiff::Mom_v) = K(i,j,k,EddyDiff::Mom_h);
      K(i,j,k,EddyDiff::Theta_v) = K(i,j,k,EddyDiff::Theta_h);
      K(i,j,k,EddyDiff::KE_v) = K(i,j,k,EddyDiff::KE_h);
      K(i,j,k,EddyDiff::Scalar_v) = K(i,j,k,EddyDiff::Scalar_h);
        });
      }
    }

    for ( MFIter mfi(eddyViscosity,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box &bx = mfi.growntilebox(1);

        if (!(domain.contains(bx))) {
            // Fill values outside the domain by straightforward extrapolation.  Note this must be
            // done separately from the loop above so all the interior values are filled.  We also
            // do this in three separate loops so that we don't have any race conditions.
            const Array4<Real> &K = eddyViscosity.array(mfi);
            amrex::ParallelFor(bx, (int) EddyDiff::NumDiffs, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                if (i < dom_lo.x) {
                    if (j < dom_lo.y)
                        K(i,j,k,n) = K(dom_lo.x,dom_lo.y,k,n);
                    else if (j > dom_hi.y)
                        K(i,j,k,n) = K(dom_lo.x,dom_hi.y,k,n);
                    else
                       K(i,j,k,n) = K(dom_lo.x,j,k,n);
                } else if (i > dom_hi.x) {
                    if (j < dom_lo.y)
                        K(i,j,k,n) = K(dom_hi.x,dom_lo.y,k,n);
                    else if (j > dom_hi.y)
                        K(i,j,k,n) = K(dom_hi.x,dom_hi.y,k,n);
                    else
                        K(i,j,k,n) = K(dom_hi.x,j,k,n);
                } else if (j < dom_lo.y) {
                        K(i,j,k,n) = K(i,dom_lo.y,k,n);
                } else if (j > dom_hi.y) {
                        K(i,j,k,n) = K(i,dom_hi.y,k,n);
                }
           });
        }
    } //mfi

    // Fill interior ghost cells and any ghost cells outside a periodic domain
    eddyViscosity.FillBoundary(geom.periodicity());

    // Now extend to low and high in vertical (assuming not periodic in z)
    if (!geom.isPeriodic(2))
    {
        for ( MFIter mfi(eddyViscosity,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

            const Box &bx = mfi.growntilebox(1);

            if (!(domain.contains(bx))) {
                const Array4<Real> &K = eddyViscosity.array(mfi);
                amrex::ParallelFor(bx, (int) EddyDiff::NumDiffs, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (k < dom_lo.z) K(i,j,k,n) = K(i,j,dom_lo.z,n);
                    if (k > dom_hi.z) K(i,j,k,n) = K(i,j,dom_hi.z,n);
                });
            }
        } //mfi
    }

} // function call

AMREX_GPU_DEVICE
inline Real
ComputeTKEProduction (const int &i, const int &j, const int &k,
                      const Array4<const Real>& u, const Array4<const Real>& v, const Array4<const Real>& w,
                      const Array4<Real>& K_turb,
                      const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                      const Box& domain, const amrex::BCRec* bc_ptr)
{
    amrex::Real TKE_production = K_turb(i,j,k,EddyDiff::Mom_h) * ComputeSmnSmn(i,j,k,u,v,w,cellSizeInv,domain,bc_ptr);

    return TKE_production;
}

/// Compute K_turb (i-1/2, j+1/2, k) etc given K_turb (i, j, k) is known
// Note: This should be at edges for momEqnDir != diffDir, cell centers otherwise
AMREX_GPU_DEVICE
inline Real
InterpolateTurbulentViscosity(const int &i, const int &j, const int &k,
                              const enum MomentumEqn &momentumEqn,
                              const enum DiffusionDir &diffDir,
                              const Array4<Real>& K_turb) {
  // Assuming we already have 'K_turb' computed for all (i, j, k)
  Real turbViscInterpolated = 1.0;

  switch (momentumEqn) {
  case MomentumEqn::x: // Reference face is x-face index (i, j, k)
    switch (diffDir) {
    case DiffusionDir::x:
      turbViscInterpolated = K_turb(i-1, j, k, EddyDiff::Mom_h);
      break;
    case DiffusionDir::y:
      turbViscInterpolated = 0.25*( K_turb(i-1, j, k, EddyDiff::Mom_h) + K_turb(i, j, k, EddyDiff::Mom_h)
                    + K_turb(i-1, j-1, k, EddyDiff::Mom_h) + K_turb(i, j-1, k, EddyDiff::Mom_h) );
      break;
    case DiffusionDir::z:
      turbViscInterpolated = 0.25*( K_turb(i-1, j, k, EddyDiff::Mom_v) + K_turb(i, j, k, EddyDiff::Mom_v)
                    + K_turb(i-1, j, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );
      break;
    default:
      amrex::Abort("Error: Diffusion direction is unrecognized");
    }
    break;
  case MomentumEqn::y: // Reference face is y-face index (i, j, k)
    switch (diffDir) {
    case DiffusionDir::x:
      turbViscInterpolated = 0.25*( K_turb(i, j-1, k, EddyDiff::Mom_h) + K_turb(i, j, k, EddyDiff::Mom_h)
                    + K_turb(i-1, j-1, k, EddyDiff::Mom_h) + K_turb(i-1, j, k, EddyDiff::Mom_h) );
      break;
    case DiffusionDir::y:
      turbViscInterpolated = K_turb(i, j-1, k, EddyDiff::Mom_h);
      break;
    case DiffusionDir::z:
      turbViscInterpolated = 0.25*( K_turb(i, j-1, k, EddyDiff::Mom_v) + K_turb(i, j, k, EddyDiff::Mom_v)
                    + K_turb(i, j-1, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );
      break;
    default:
      amrex::Abort("Error: Diffusion direction is unrecognized");
    }
    break;
  case MomentumEqn::z: // Reference face is z-face index (i, j, k)
    switch (diffDir) {
    case DiffusionDir::x:
      turbViscInterpolated = 0.25*( K_turb(i, j, k-1, EddyDiff::Mom_v) + K_turb(i, j, k, EddyDiff::Mom_v)
                    + K_turb(i-1, j, k-1, EddyDiff::Mom_v) + K_turb(i-1, j, k, EddyDiff::Mom_v) );
      break;
    case DiffusionDir::y:
      turbViscInterpolated = 0.25*( K_turb(i, j, k-1, EddyDiff::Mom_v) + K_turb(i, j, k, EddyDiff::Mom_v)
                    + K_turb(i, j-1, k-1, EddyDiff::Mom_v) + K_turb(i, j-1, k, EddyDiff::Mom_v) );
      break;
    case DiffusionDir::z:
      turbViscInterpolated = K_turb(i, j, k-1, EddyDiff::Mom_v);
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
