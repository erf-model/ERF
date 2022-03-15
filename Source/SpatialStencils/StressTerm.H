#include <TimeIntegration.H>
#include <EddyViscosity.H>
#include <ExpansionRate.H>
#include <StrainRate.H>

using namespace amrex;

#ifndef _COMPUTE_STRESSTERM_H
#define _COMPUTE_STRESSTERM_H

using namespace amrex;

// Compute tau_ij (m + 1/2), tau_ij (m - 1/2) where m = {i, j, k} for DNS or Smagorinsky
//    given the strain rates
AMREX_GPU_DEVICE
inline Real
ComputeStressGivenRates (const int &i, const int &j, const int &k,
                         const enum MomentumEqn &momentumEqn,
                         const enum DiffusionDir &diffDir,
                         const Array4<Real>& K_turb,
                         const SolverChoice &solverChoice,
                         amrex::Real strainRateDeviatoric)
{
    Real mu_effective = 0.0;
    //
    // NOTE: we do not allow LES with MolecDiffType::ConstantAlpha so we don't have to worry about the
    //       density-weighting of the LES term in this case.  When using MolecDiffType::ConstantAlpha,
    //       we multiply the entire stress term by (rho / rho0_trans) in DiffusionContributionForMom
    //       (which is calling this routine)
    //
    if ( (solverChoice.molec_diff_type == MolecDiffType::Constant) ||
         (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha) ) {
        mu_effective += 2.0 * solverChoice.dynamicViscosity; // 2*mu
    }

    if ( (solverChoice.les_type == LESType::Smagorinsky) ||
         (solverChoice.les_type == LESType::Deardorff  ) ) {
        mu_effective += InterpolateTurbulentViscosity(i, j, k, momentumEqn, diffDir, K_turb); // 2*mu_t
    } else if (solverChoice.les_type != LESType::None) {
        amrex::Abort("Error:  LES model is unrecognized");
    }

    Real stressTerm = mu_effective * strainRateDeviatoric; // tau_ij = mu_effective * sigma_ij
    return stressTerm;
}

// Compute tau_ij (m + 1/2), tau_ij (m - 1/2) where m = {i, j, k} for DNS or Smagorinsky
AMREX_GPU_DEVICE
inline Real
ComputeStressTerm (const int &i, const int &j, const int &k,
                   const Array4<const Real>& u, const Array4<const Real>& v, const Array4<const Real>& w,
                   const enum MomentumEqn &momentumEqn,
                   const enum DiffusionDir &diffDir,
                   const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                   const Array4<Real>& K_turb,
                   const SolverChoice &solverChoice,
                   const Box& domain, const amrex::BCRec* bc_ptr)
{

    // Here, we have computed strain rate on the fly.
    // TODO: It may be better to store S11, S12 etc. at all the (m+1/2) and (m-1/2) grid points (edges) and use them here.
    Real strainRate = ComputeStrainRate(i, j, k, u, v, w, momentumEqn, diffDir, cellSizeInv, domain, bc_ptr);

    // D_ij term
    Real expansionRate = ComputeExpansionRate(i, j, k, u, v, w, momentumEqn, diffDir, cellSizeInv);

    Real strainRateDeviatoric = strainRate - expansionRate; // sigma_ij = S_ij - D_ij

    Real stressTerm = ComputeStressGivenRates(i,j,k,momentumEqn,diffDir,K_turb,solverChoice, strainRateDeviatoric);

    return stressTerm;
}
#endif
