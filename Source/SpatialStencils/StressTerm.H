#include <EddyViscosity.H>
#include <ExpansionRate.H>
#include <StrainRate.H>

#ifndef _COMPUTE_STRESSTERM_H
#define _COMPUTE_STRESSTERM_H

// Compute tau_ij (m + 1/2), tau_ij (m - 1/2) where m = {i, j, k} for DNS or Smagorinsky
//    given the strain rates
AMREX_GPU_DEVICE
inline amrex::Real
ComputeStressGivenRates (const int &i, const int &j, const int &k,
                         const enum MomentumEqn &momentumEqn,
                         const enum DiffusionDir &diffDir,
                         const amrex::Array4<amrex::Real>& K_turb,
                         const SolverChoice &solverChoice,
                         amrex::Real strainRateDeviatoric)
{
    amrex::Real mu_effective = 0.0;
    //
    // NOTE: we do not allow LES with MolecDiffType::ConstantAlpha so we don't have to worry about the
    //       density-weighting of the LES term in this case.  When using MolecDiffType::ConstantAlpha,
    //       we multiply the entire stress term by (rho / rho0_trans) in DiffusionSrcForMom
    //       (which is calling this routine)
    //
    if ( (solverChoice.molec_diff_type == MolecDiffType::Constant) ||
         (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha) ) {
        mu_effective += 2.0 * solverChoice.dynamicViscosity; // 2*mu
    }

    if ( (solverChoice.les_type == LESType::Smagorinsky) ||
         (solverChoice.les_type == LESType::Deardorff  ) ||
         (solverChoice.pbl_type == PBLType::MYNN25     )  ) {
        mu_effective += InterpolateTurbulentViscosity(i, j, k, momentumEqn, diffDir, K_turb); // 2*mu_t
    }

    amrex::Real stressTerm = mu_effective * strainRateDeviatoric; // tau_ij = mu_effective * sigma_ij
    return stressTerm;
}

#ifdef ERF_USE_TERRAIN
// Compute tau_ij (m + 1/2), tau_ij (m - 1/2) where m = {i, j, k} for DNS or Smagorinsky
AMREX_GPU_DEVICE
inline amrex::Real
ComputeStressTerm (const int &i, const int &j, const int &k,
                   const amrex::Array4<const amrex::Real>& u,
                   const amrex::Array4<const amrex::Real>& v,
                   const amrex::Array4<const amrex::Real>& w,
                   const enum MomentumEqn &momentumEqn,
                   const enum DiffusionDir &diffDir,
                   const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                   const amrex::Array4<amrex::Real>& K_turb,
                   const SolverChoice &solverChoice,
                   const amrex::Array4<const amrex::Real>& z_nd,
                   const amrex::Box& domain, const amrex::BCRec* bc_ptr)
{
    // S_ij term 
    amrex::Real strainRate = ComputeStrainRate(i, j, k, u, v, w, momentumEqn, diffDir, cellSizeInv, z_nd, domain, bc_ptr);

    // D_ij term
    amrex::Real expansionRate = ComputeExpansionRate(i, j, k, u, v, w, momentumEqn, diffDir, z_nd, cellSizeInv);

    // sigma_ij = S_ij - D_ij
    amrex::Real strainRateDeviatoric = strainRate - expansionRate; 

    amrex::Real stressTerm = ComputeStressGivenRates(i, j, k, momentumEqn, diffDir, K_turb, solverChoice, strainRateDeviatoric);

    return stressTerm;
}

#else
// Compute tau_ij (m + 1/2), tau_ij (m - 1/2) where m = {i, j, k} for DNS or Smagorinsky
AMREX_GPU_DEVICE
inline amrex::Real
ComputeStressTerm (const int &i, const int &j, const int &k,
                   const amrex::Array4<const amrex::Real>& u,
                   const amrex::Array4<const amrex::Real>& v,
                   const amrex::Array4<const amrex::Real>& w,
                   const enum MomentumEqn &momentumEqn,
                   const enum DiffusionDir &diffDir,
                   const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                   const amrex::Array4<amrex::Real>& K_turb,
                   const SolverChoice &solverChoice,
                   const amrex::Box& domain, const amrex::BCRec* bc_ptr)
{
    amrex::Real strainRate = ComputeStrainRate(i, j, k, u, v, w, momentumEqn, diffDir, cellSizeInv, domain, bc_ptr);

    // D_ij term
    amrex::Real expansionRate = ComputeExpansionRate(i, j, k, u, v, w, momentumEqn, diffDir, cellSizeInv);

    amrex::Real strainRateDeviatoric = strainRate - expansionRate; // sigma_ij = S_ij - D_ij

    amrex::Real stressTerm = ComputeStressGivenRates(i,j,k,momentumEqn,diffDir,K_turb,solverChoice, strainRateDeviatoric);

    return stressTerm;
}
#endif

#endif
