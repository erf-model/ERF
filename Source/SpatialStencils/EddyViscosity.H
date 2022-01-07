/** \file EddyViscosity.H */

#ifndef _EDDY_VISCOSITY_H_
#define _EDDY_VISCOSITY_H_

#include <TimeIntegration.H>
#include "StrainRate.H"

using namespace amrex;

//AMREX_GPU_DEVICE
inline
Real ComputeSmnSmn(const int& i,const int& j,const int& k,
                   const Array4<Real const>& u,
                   const Array4<Real const>& v,
                   const Array4<Real const>& w,
                   const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                   const Box& domain,
                   const Gpu::DeviceVector<amrex::BCRec> domain_bcs_type_d)
{
    Real S11 = ComputeS11(i+1, j, k, u, cellSizeInv);
    Real S22 = ComputeS22(i, j+1, k, v, cellSizeInv);
    Real S33 = ComputeS33(i, j, k+1, w, cellSizeInv);

    Real S12 = 0.25* (
              ComputeS12(i  , j  , k, u, v, w, cellSizeInv, domain, domain_bcs_type_d)
            + ComputeS12(i  , j+1, k, u, v, w, cellSizeInv, domain, domain_bcs_type_d)
            + ComputeS12(i+1, j  , k, u, v, w, cellSizeInv, domain, domain_bcs_type_d)
            + ComputeS12(i+1, j+1, k, u, v, w, cellSizeInv, domain, domain_bcs_type_d)
            );

    Real S13 = 0.25* (
              ComputeS13(i  , j, k  , u, v, w,  cellSizeInv, domain, domain_bcs_type_d)
            + ComputeS13(i  , j, k+1, u, v, w,  cellSizeInv, domain, domain_bcs_type_d)
            + ComputeS13(i+1, j, k  , u, v, w,  cellSizeInv, domain, domain_bcs_type_d)
            + ComputeS13(i+1, j, k+1, u, v, w,  cellSizeInv, domain, domain_bcs_type_d)
            );

    Real S23 = 0.25* (
              ComputeS23(i, j  , k  , u, v, w,  cellSizeInv, domain, domain_bcs_type_d)
            + ComputeS23(i, j  , k+1, u, v, w,  cellSizeInv, domain, domain_bcs_type_d)
            + ComputeS23(i, j+1, k  , u, v, w,  cellSizeInv, domain, domain_bcs_type_d)
            + ComputeS23(i, j+1, k+1, u, v, w,  cellSizeInv, domain, domain_bcs_type_d)
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
            Real SmnSmn = ComputeSmnSmn(i,j,k,u,v,w,cellSizeInv,domain,domain_bcs_type_d);

            // Note the positive sign, which aligns well with the positive sign in the diffusion term for momentum equation
            K(i, j, k, 0) = 2.0 * CsDeltaSqr * cell_data(i, j, k, Rho_comp) * std::sqrt(2.0*SmnSmn);
        });

        // Fill values outside the domain by straightforwad extrapolation.  Note this must be
        // done separately from the loop above so all the interior values are filled
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (i < dom_lo.x) K(i,j,k,0) = K(dom_lo.x,j,k,0);
            if (j < dom_lo.y) K(i,j,k,0) = K(i,dom_lo.y,k,0);
            if (k < dom_lo.z) K(i,j,k,0) = K(i,j,dom_lo.z,0);
            if (i > dom_hi.x) K(i,j,k,0) = K(dom_hi.x,j,k,0);
            if (j > dom_hi.y) K(i,j,k,0) = K(i,dom_hi.y,k,0);
            if (k > dom_hi.z) K(i,j,k,0) = K(i,j,dom_hi.z,0);
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
            K(i,j,k,0) = l_C_k * Delta *  cell_data(i,j,k,Rho_comp) *
                                std::sqrt(cell_data(i,j,k,RhoKE_comp));
        });

        // Fill values outside the domain by straightforwad extrapolation.  Note this must be
        // done separately from the loop above so all the interior values are filled
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (i < dom_lo.x) K(i,j,k,0) = K(dom_lo.x,j,k,0);
            if (j < dom_lo.y) K(i,j,k,0) = K(i,dom_lo.y,k,0);
            if (k < dom_lo.z) K(i,j,k,0) = K(i,j,dom_lo.z,0);
            if (i > dom_hi.x) K(i,j,k,0) = K(dom_hi.x,j,k,0);
            if (j > dom_hi.y) K(i,j,k,0) = K(i,dom_hi.y,k,0);
            if (k > dom_hi.z) K(i,j,k,0) = K(i,j,dom_hi.z,0);
        });
      } //mfi
    }

    // Fill interior ghost cells and any ghost cells outside a periodic domain
    eddyViscosity.FillBoundary(geom.periodicity());

} // function call

AMREX_GPU_DEVICE
inline Real
ComputeTKEProduction (const int &i, const int &j, const int &k,
                      const Array4<const Real>& u, const Array4<const Real>& v, const Array4<const Real>& w,
                      const Array4<Real>& K_LES,
                      const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                      const SolverChoice &solverChoice,
                      const Box& domain,
                      const Gpu::DeviceVector<amrex::BCRec> domain_bcs_type_d)
{
    amrex::Real TKE_production = K_LES(i,j,k) * ComputeSmnSmn(i,j,k,u,v,w,cellSizeInv,domain,domain_bcs_type_d);

    return TKE_production;
}

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
