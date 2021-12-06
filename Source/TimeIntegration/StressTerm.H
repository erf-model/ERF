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
                         const Array4<Real>& K_LES,
                         const SolverChoice &solverChoice,
                         amrex::Real strainRateDeviatoric)
{
    Real mu_effective = 0.0;
    //TODO: dynamic viscosity, mu, is assumed to be constant in the current implementation.
    // Future implementations may account for mu = mu(T) computed at the coordinate of interest.
    // That could be done with a new MolecDiffType
    switch (solverChoice.molec_diff_type) {
        case MolecDiffType::Constant:
            mu_effective += 2.0 * solverChoice.dynamicViscosity; // 2*mu
            break;
        case MolecDiffType::None:
            break;
        default:
            amrex::Abort("Error: Molecular diffusion/viscosity model is unrecognized");
    }

    Real turbViscInterpolated = 0.0;
    // TODO: Add Deardorff model, perhaps take advantage of turbulence model indicator
    switch (solverChoice.les_type) {
        case LESType::Smagorinsky:
            // mu_effective = 2*mu + 2*mu_t if MolecDiffType::Constant else 2*mu_t
            mu_effective += InterpolateTurbulentViscosity(i, j, k, momentumEqn, diffDir, K_LES); // 2*mu_t
            break;
        case LESType::Deardorff:
            // mu_effective = 2*mu + 2*mu_t if MolecDiffType::Constant else 2*mu_t
            mu_effective += InterpolateTurbulentViscosity(i, j, k, momentumEqn, diffDir, K_LES); // 2*mu_t
        case LESType::None: // // mu_effective = 2*mu if MolecDiffType::Constant else 0
            break;
        default:
            amrex::Abort("Error:  LES model is unrecognized");
    }

    Real stressTerm = mu_effective * strainRateDeviatoric; // tau_ij = mu_effective * sigma_ij
    return stressTerm;
}

// Compute tau_ij (m + 1/2), tau_ij (m - 1/2) where m = {i, j, k} for DNS or Smagorinsky
AMREX_GPU_DEVICE
inline Real
ComputeStressTerm (const int &i, const int &j, const int &k,
                   const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                   const enum MomentumEqn &momentumEqn,
                   const enum DiffusionDir &diffDir,
                   const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                   const Array4<Real>& K_LES,
                   const SolverChoice &solverChoice,
                   bool use_no_slip_stencil_lo,
                   bool use_no_slip_stencil_hi)
{

    // Here, we have computed strain rate on the fly.
    // TODO: It may be better to store S11, S12 etc. at all the (m+1/2) and (m-1/2) grid points (edges) and use them here.
    Real strainRate = ComputeStrainRate(i, j, k, u, v, w, momentumEqn, diffDir, cellSizeInv,
                                       use_no_slip_stencil_lo, use_no_slip_stencil_hi);

    // D_ij term
    Real expansionRate = ComputeExpansionRate(i, j, k, u, v, w, momentumEqn, diffDir, cellSizeInv);

    Real strainRateDeviatoric = strainRate - expansionRate; // sigma_ij = S_ij - D_ij

    Real stressTerm = ComputeStressGivenRates(i,j,k,momentumEqn,diffDir,K_LES,solverChoice,
                                                strainRateDeviatoric);

    return stressTerm;
}
#endif
