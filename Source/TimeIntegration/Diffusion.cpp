// Created by Pankaj Jha on 6/15/21.
#include <TimeIntegration.H>

using namespace amrex;

// Compute tau_ij (m + 1/2), tau_ij (m - 1/2) where m = {i, j, k} for DNS or Smagorinsky
AMREX_GPU_DEVICE
Real ComputeStressTerm (const int &i, const int &j, const int &k,
                        const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                        const enum NextOrPrev &nextOrPrev,
                        const enum MomentumEqn &momentumEqn,
                        const enum DiffusionDir &diffDir,
                        const GpuArray<Real, AMREX_SPACEDIM>& cellSize,
                        const Array4<Real>& Ksmag,
                        const SolverChoice &solverChoice) {
    Real stressTerm = 0.0;
    Real turbViscInterpolated = 0.0;

    // Here, we have computed strain rate on the fly.
    // TODO: It may be better to store S11, S12 etc. at all the (m+1/2) and (m-1/2) grid points (edges) and use them here.
    Real strainRate = ComputeStrainRate(i, j, k, u, v, w, nextOrPrev, momentumEqn, diffDir, cellSize);

    // TODO: Consider passing turbModel to this function instead of computing it here from SolverChoice
    enum TurbulenceModel turbModel;

    //TODO: Update this to account for other turbulence models in future. This would alo require update in SolverChoice
    if (solverChoice.use_smagorinsky)
        turbModel = TurbulenceModel::Smagorinsky;
    else
        turbModel = TurbulenceModel::DNS;

    switch (turbModel) {
        case TurbulenceModel::DNS:
            //TODO: dynamic viscosity, mu, is assumed to be constant in the current implementation.
            // Future implementations may account for mu = mu(T) computed at the coordinate of interest.
            stressTerm = 2.0 * solverChoice.dynamicViscosity * strainRate; // 2*mu*Sij(m+1/2) or 2*mu*Sij(m-1/2)
            break;
        case TurbulenceModel::Smagorinsky:
            turbViscInterpolated = InterpolateTurbulentViscosity(i, j, k, u, v, w, nextOrPrev, momentumEqn, diffDir, cellSize, Ksmag);
            stressTerm = turbViscInterpolated * strainRate; // // K_interp*Sij(m+1/2) or K_interp*Sij(m-1/2)
            break;
        default:
            amrex::Abort("Error: Turbulence model is unrecognized");
    }

    return stressTerm;
}

AMREX_GPU_DEVICE
Real
DiffusionContributionForMom(const int &i, const int &j, const int &k,
                            const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                            const enum MomentumEqn &momentumEqn,
                            const GpuArray<Real, AMREX_SPACEDIM>& cellSize,
                            const Array4<Real>& Ksmag,
                            const SolverChoice &solverChoice) {

    auto dx = cellSize[0], dy = cellSize[1], dz = cellSize[2];
    Real diffusionContribution = 0.0;

    switch (momentumEqn) {
        case MomentumEqn::x:
            Real tau11Next, tau11Prev, tau12Next, tau12Prev, tau13Next, tau13Prev;
            tau11Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, momentumEqn, DiffusionDir::x, cellSize, Ksmag, solverChoice);
            tau11Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, momentumEqn, DiffusionDir::x, cellSize, Ksmag, solverChoice);
            tau12Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, momentumEqn, DiffusionDir::y, cellSize, Ksmag, solverChoice);
            tau12Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, momentumEqn, DiffusionDir::y, cellSize, Ksmag, solverChoice);
            tau13Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, momentumEqn, DiffusionDir::z, cellSize, Ksmag, solverChoice);
            tau13Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, momentumEqn, DiffusionDir::z, cellSize, Ksmag, solverChoice);

            diffusionContribution = (tau11Next - tau11Prev) / dx  // Contribution to x-mom eqn from diffusive flux in x-dir
                                  + (tau12Next - tau12Prev) / dy  // Contribution to x-mom eqn from diffusive flux in y-dir
                                  + (tau13Next - tau13Prev) / dz; // Contribution to x-mom eqn from diffusive flux in z-dir
            break;
        case MomentumEqn::y:
            Real tau21Next, tau21Prev, tau22Next, tau22Prev, tau23Next, tau23Prev;
            tau21Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, momentumEqn, DiffusionDir::x, cellSize, Ksmag, solverChoice);
            tau21Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, momentumEqn, DiffusionDir::x, cellSize, Ksmag, solverChoice);
            tau22Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, momentumEqn, DiffusionDir::y, cellSize, Ksmag, solverChoice);
            tau22Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, momentumEqn, DiffusionDir::y, cellSize, Ksmag, solverChoice);
            tau23Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, momentumEqn, DiffusionDir::z, cellSize, Ksmag, solverChoice);
            tau23Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, momentumEqn, DiffusionDir::z, cellSize, Ksmag, solverChoice);

            diffusionContribution = (tau21Next - tau21Prev) / dx  // Contribution to y-mom eqn from diffusive flux in x-dir
                                  + (tau22Next - tau22Prev) / dy  // Contribution to y-mom eqn from diffusive flux in y-dir
                                  + (tau23Next - tau23Prev) / dz; // Contribution to y-mom eqn from diffusive flux in z-dir
            break;
        case MomentumEqn::z:
            Real tau31Next, tau31Prev, tau32Next, tau32Prev, tau33Next, tau33Prev;
            tau31Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, momentumEqn, DiffusionDir::x, cellSize, Ksmag, solverChoice);
            tau31Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, momentumEqn, DiffusionDir::x, cellSize, Ksmag, solverChoice);
            tau32Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, momentumEqn, DiffusionDir::y, cellSize, Ksmag, solverChoice);
            tau32Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, momentumEqn, DiffusionDir::y, cellSize, Ksmag, solverChoice);
            tau33Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, momentumEqn, DiffusionDir::z, cellSize, Ksmag, solverChoice);
            tau33Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, momentumEqn, DiffusionDir::z, cellSize, Ksmag, solverChoice);

            diffusionContribution = (tau31Next - tau31Prev) / dx  // Contribution to z-mom eqn from diffusive flux in x-dir
                                  + (tau32Next - tau32Prev) / dy  // Contribution to z-mom eqn from diffusive flux in y-dir
                                  + (tau33Next - tau33Prev) / dz; // Contribution to z-mom eqn from diffusive flux in z-dir
            break;
        default:
            amrex::Abort("Error: Momentum equation is unrecognized");
    }

    return diffusionContribution;
}

AMREX_GPU_DEVICE
amrex::Real ComputeDiffusionFluxForState(const int &i, const int &j, const int &k,
					 const Array4<Real>& cell_data, const int & qty_index,
					 const amrex::Real invCellWidth,
					 const Array4<Real>& Ksmag,
					 const SolverChoice &solverChoice,
					 const enum NextOrPrev &nextOrPrev,
					 const enum Coord& coordDir) {
  //TODO: Discuss about the implementation changes needed (if any). The current implemenation is based on a previous
  // version of documentation.

  // Get indices of states to left and right of the face on which we want the flux
  // TODO: Make a templated class that does  this automatically
  int ileft = i; int iright = i;
  int jleft = j; int jright = j;
  int kleft = k; int kright = k;
  switch (coordDir) {
  case Coord::x:
    if (nextOrPrev == NextOrPrev::next) {
      iright = i+1;
    } else {
      ileft = i-1;
    }
    break;
  case Coord::y:
    if (nextOrPrev == NextOrPrev::next) {
      jright = j+1;
    } else {
      jleft = j-1;
    }
    break;
  case Coord::z:
    if (nextOrPrev == NextOrPrev::next) {
      kright = k+1;
    } else {
      kleft = k-1;
    }
    break;
  default:
    amrex::Abort("Error: Coord direction is unrecognized");
  }

  // Get diffusion coefficients
  amrex::Real rhoAlpha = 0.0;
  amrex::Real rhoAlpha_molec;
  amrex::Real Pr_or_Sc_turb_inv;
  switch(qty_index) {
  case RhoTheta_comp: // Temperature
    rhoAlpha_molec = solverChoice.rhoAlpha_T;
    Pr_or_Sc_turb_inv = solverChoice.Pr_t_inv;
    break;
  case RhoScalar_comp: // Scalar
    rhoAlpha_molec = solverChoice.rhoAlpha_C;
    Pr_or_Sc_turb_inv = solverChoice.Sc_t_inv;
    break;
  default:
    amrex::Abort("Error: Diffusion term for the data index isn't implemented");
  }
  
  //TODO: Update this to account for other turbulence models in future. This would alo require update in SolverChoice
  enum TurbulenceModel turbModel;
  amrex::Real rhoAlpha_right, rhoAlpha_left;
  if (solverChoice.use_smagorinsky)
    turbModel = TurbulenceModel::Smagorinsky;
  else
    turbModel = TurbulenceModel::DNS;
  switch (turbModel) {
  case TurbulenceModel::DNS:
    //TODO: molecular transport coefficients are assumed to be constant in the current implementation.
    // note this assumption applies to the product, e.g. rhoAlpha_T = lambda/c_p
    // Future implementations may account for rhoAlpha = rhoAlpha(T) computed at the coordinate of interest.
    rhoAlpha = rhoAlpha_molec;
    break;
  case TurbulenceModel::Smagorinsky:
    // Ksmag = 2*mu_t
    rhoAlpha_right = Ksmag(ileft, jleft, kleft) * 0.5 * Pr_or_Sc_turb_inv;
    rhoAlpha_left = Ksmag(ileft, jleft, kleft) * 0.5 * Pr_or_Sc_turb_inv;
    rhoAlpha = 0.5*(rhoAlpha_left + rhoAlpha_right);
    break;
  default:
    amrex::Abort("Error: Turbulence model is unrecognized");
  }

  // Compute the flux
  amrex::Real diffusionFlux = rhoAlpha * (cell_data(iright, jright, kright, qty_index) /
					  cell_data(iright, jright, kright, Rho_comp)
					  - cell_data(ileft, jleft, kleft, qty_index) /
					  cell_data(ileft, jleft, kleft, Rho_comp));
  
  return diffusionFlux;
}

AMREX_GPU_DEVICE
Real
DiffusionContributionForState(const int &i, const int &j, const int &k,
                              const Array4<Real>& cell_data, const int & qty_index,
                              const GpuArray<Real, AMREX_SPACEDIM>& cellSize,
                              const Array4<Real>& Ksmag,
                              const SolverChoice &solverChoice) {

  const amrex::Real dx_inv = 1.0/cellSize[0];
  const amrex::Real dy_inv = 1.0/cellSize[1];
  const amrex::Real dz_inv = 1.0/cellSize[2];

    amrex::Real xDiffFluxNext = ComputeDiffusionFluxForState(i, j, k, cell_data, qty_index, dx_inv, Ksmag, solverChoice, NextOrPrev::next, Coord::x);
    amrex::Real yDiffFluxNext = ComputeDiffusionFluxForState(i, j, k, cell_data, qty_index, dy_inv, Ksmag, solverChoice, NextOrPrev::next, Coord::y);
    amrex::Real zDiffFluxNext = ComputeDiffusionFluxForState(i, j, k, cell_data, qty_index, dz_inv, Ksmag, solverChoice, NextOrPrev::next, Coord::z);
    
    amrex::Real xDiffFluxPrev = ComputeDiffusionFluxForState(i, j, k, cell_data, qty_index, dx_inv, Ksmag, solverChoice, NextOrPrev::prev, Coord::x);
    amrex::Real yDiffFluxPrev = ComputeDiffusionFluxForState(i, j, k, cell_data, qty_index, dy_inv, Ksmag, solverChoice, NextOrPrev::prev, Coord::y);
    amrex::Real zDiffFluxPrev = ComputeDiffusionFluxForState(i, j, k, cell_data, qty_index, dz_inv, Ksmag, solverChoice, NextOrPrev::prev, Coord::z);

    //TODO: Discuss about the implementation changes needed (if any). Diffusion coefficients are assumed to be constant.
    switch(qty_index) {
        case RhoTheta_comp: // Temperature
	  //diffCoeff = solverChoice.alpha_T;
            break;
        case RhoScalar_comp: // Scalar
	  //diffCoeff = solverChoice.alpha_C;
            break;
        default:
            amrex::Abort("Error: Diffusion term for the data index isn't implemented");
    }

    // Assemble diffusion contribution.
    Real diffusionContribution = 
      (xDiffFluxNext - xDiffFluxPrev) * dx_inv   // Diffusive flux in x-dir
      + (yDiffFluxNext - yDiffFluxPrev) * dy_inv   // Diffusive flux in y-dir
      + (zDiffFluxNext - zDiffFluxPrev) * dz_inv;  // Diffusive flux in z-dir

    return diffusionContribution;
}
