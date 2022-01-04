#include <SpatialStencils.H>
#include <EddyViscosity.H>
#include <ExpansionRate.H>
#include <StrainRate.H>
#include <StressTerm.H>

using namespace amrex;

AMREX_GPU_DEVICE
Real
DiffusionContributionForMom(const int &i, const int &j, const int &k,
                            const Array4<const Real>& u, const Array4<const Real>& v, const Array4<const Real>& w,
                            const enum MomentumEqn &momentumEqn,
                            const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                            const Array4<Real>& K_LES,
                            const SolverChoice &solverChoice,
                            const bool dirichlet_at_lo_k,
                            const bool dirichlet_at_hi_k)
{
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    Real diffusionContribution = 0.0;

    switch (momentumEqn) {
        case MomentumEqn::x:
            Real tau11Next, tau11Prev, tau12Next, tau12Prev, tau13Next, tau13Prev;
            tau11Next = ComputeStressTerm(i+1, j, k, u, v, w, momentumEqn,
                                          DiffusionDir::x, cellSizeInv, K_LES, solverChoice, false, false);
            tau11Prev = ComputeStressTerm(i  , j, k, u, v, w, momentumEqn,
                                          DiffusionDir::x, cellSizeInv, K_LES, solverChoice, false, false);
            tau12Next = ComputeStressTerm(i, j+1, k, u, v, w, momentumEqn,
                                          DiffusionDir::y, cellSizeInv, K_LES, solverChoice, false, false);
            tau12Prev = ComputeStressTerm(i, j  , k, u, v, w, momentumEqn,
                                          DiffusionDir::y, cellSizeInv, K_LES, solverChoice, false, false);
            tau13Next = ComputeStressTerm(i, j, k+1, u, v, w, momentumEqn,
                                          DiffusionDir::z, cellSizeInv, K_LES, solverChoice,
                                          false, dirichlet_at_hi_k);
            tau13Prev = ComputeStressTerm(i, j, k  , u, v, w, momentumEqn,
                                          DiffusionDir::z, cellSizeInv, K_LES, solverChoice,
                                          dirichlet_at_lo_k, false);

            diffusionContribution = (tau11Next - tau11Prev) * dxInv  // Contribution to x-mom eqn from diffusive flux in x-dir
                                  + (tau12Next - tau12Prev) * dyInv  // Contribution to x-mom eqn from diffusive flux in y-dir
                                  + (tau13Next - tau13Prev) * dzInv; // Contribution to x-mom eqn from diffusive flux in z-dir
            break;
        case MomentumEqn::y:
            Real tau21Next, tau21Prev, tau22Next, tau22Prev, tau23Next, tau23Prev;
            tau21Next = ComputeStressTerm(i+1, j, k, u, v, w, momentumEqn,
                                          DiffusionDir::x, cellSizeInv, K_LES, solverChoice, false, false);
            tau21Prev = ComputeStressTerm(i  , j, k, u, v, w, momentumEqn,
                                           DiffusionDir::x, cellSizeInv, K_LES, solverChoice, false, false);
            tau22Next = ComputeStressTerm(i, j+1, k, u, v, w, momentumEqn,
                                          DiffusionDir::y, cellSizeInv, K_LES, solverChoice, false, false);
            tau22Prev = ComputeStressTerm(i, j  , k, u, v, w, momentumEqn,
                                          DiffusionDir::y, cellSizeInv, K_LES, solverChoice, false, false);
            tau23Next = ComputeStressTerm(i, j, k+1, u, v, w, momentumEqn,
                                          DiffusionDir::z, cellSizeInv, K_LES, solverChoice,
                                          false, dirichlet_at_hi_k);
            tau23Prev = ComputeStressTerm(i, j, k  , u, v, w, momentumEqn,
                                          DiffusionDir::z, cellSizeInv, K_LES, solverChoice,
                                          dirichlet_at_lo_k, false);

            diffusionContribution = (tau21Next - tau21Prev) * dxInv  // Contribution to y-mom eqn from diffusive flux in x-dir
                                  + (tau22Next - tau22Prev) * dyInv  // Contribution to y-mom eqn from diffusive flux in y-dir
                                  + (tau23Next - tau23Prev) * dzInv; // Contribution to y-mom eqn from diffusive flux in z-dir
            break;
        case MomentumEqn::z:
            Real tau31Next, tau31Prev, tau32Next, tau32Prev, tau33Next, tau33Prev;
            tau31Next = ComputeStressTerm(i+1, j, k, u, v, w, momentumEqn,
                                          DiffusionDir::x, cellSizeInv, K_LES, solverChoice, false, false);
            tau31Prev = ComputeStressTerm(i  , j, k, u, v, w, momentumEqn,
                                          DiffusionDir::x, cellSizeInv, K_LES, solverChoice, false, false);
            tau32Next = ComputeStressTerm(i, j+1, k, u, v, w, momentumEqn,
                                          DiffusionDir::y, cellSizeInv, K_LES, solverChoice, false, false);
            tau32Prev = ComputeStressTerm(i, j  , k, u, v, w, momentumEqn,
                                          DiffusionDir::y, cellSizeInv, K_LES, solverChoice, false, false);
            tau33Next = ComputeStressTerm(i, j, k+1, u, v, w, momentumEqn,
                                          DiffusionDir::z, cellSizeInv, K_LES, solverChoice, false, false);
            tau33Prev = ComputeStressTerm(i, j, k  , u, v, w, momentumEqn,
                                          DiffusionDir::z, cellSizeInv, K_LES, solverChoice, false, false);

            diffusionContribution = (tau31Next - tau31Prev) * dxInv  // Contribution to z-mom eqn from diffusive flux in x-dir
                                  + (tau32Next - tau32Prev) * dyInv  // Contribution to z-mom eqn from diffusive flux in y-dir
                                  + (tau33Next - tau33Prev) * dzInv; // Contribution to z-mom eqn from diffusive flux in z-dir
            break;
        default:
            amrex::Abort("Error: Momentum equation is unrecognized");
    }

    return diffusionContribution;
}

AMREX_GPU_DEVICE
amrex::Real ComputeDiffusionFluxForState(const int &i, const int &j, const int &k,
                     const Array4<const Real>& cell_data,
                     const Array4<const Real>& cell_prim, const int & prim_index,
                     const amrex::Real invCellWidth,
                     const Array4<Real>& K_LES,
                     const SolverChoice &solverChoice,
                     const enum Coord& coordDir)
{
  // Get indices of states to left and right of the face on which we want the flux
  const int il = i - (coordDir == Coord::x);
  const int ir = i;
  const int jl = j - (coordDir == Coord::y);
  const int jr = j;
  const int kl = k - (coordDir == Coord::z);
  const int kr = k;

  // Get diffusion coefficients
  amrex::Real rhoAlpha_molec;
  amrex::Real Pr_or_Sc_turb_inv;

  amrex::Real rhoFace;
  if (solverChoice.molec_diff_type == MolecDiffType::ConstantDiffusivity) {
    rhoFace = (cell_data(il, jl, kl, Rho_comp) + cell_data(ir, jr, kr, Rho_comp)) * 0.5;
  }
  switch(prim_index) {
  case PrimTheta_comp: // Potential Temperature
    if (solverChoice.molec_diff_type == MolecDiffType::ConstantDiffusivity) {
        rhoAlpha_molec = rhoFace * solverChoice.alpha_T;
    } else {
        // rhoAlpha_T == solverChoice.rho0_trans * solverChoice.alpha_T
        rhoAlpha_molec = solverChoice.rhoAlpha_T;
    }
    Pr_or_Sc_turb_inv = solverChoice.Pr_t_inv;
    break;
  case PrimKE_comp: // Turbulent KE
    if (solverChoice.molec_diff_type == MolecDiffType::ConstantDiffusivity) {
        rhoAlpha_molec = rhoFace * solverChoice.alpha_T;
    } else {
        // rhoAlpha_T == solverChoice.rho0_trans * solverChoice.alpha_T
        rhoAlpha_molec = solverChoice.rhoAlpha_T;
    }
    Pr_or_Sc_turb_inv = solverChoice.Pr_t_inv;
    break;
  case PrimScalar_comp: // Scalar
    if (solverChoice.molec_diff_type == MolecDiffType::ConstantDiffusivity) {
        rhoAlpha_molec = rhoFace * solverChoice.alpha_C;
    } else {
        // rhoAlpha_C == solverChoice.rho0_trans * solverChoice.alpha_C
        rhoAlpha_molec = solverChoice.rhoAlpha_C;
    }
    Pr_or_Sc_turb_inv = solverChoice.Sc_t_inv;
    break;
  default:
    amrex::Abort("Error: Diffusion term for the data index isn't implemented");
  }

  amrex::Real rhoAlpha = 0.0;
  switch (solverChoice.molec_diff_type) {
  case MolecDiffType::Constant:
  case MolecDiffType::ConstantDiffusivity:
    rhoAlpha += rhoAlpha_molec;
    break;
  case MolecDiffType::None:
    break;
  default:
    amrex::Abort("Error: Molecular diffusion/viscosity model is unrecognized");
  }

  amrex::Real rhoAlpha_r = 0.0, rhoAlpha_l = 0.0;
  amrex::Real l_sigma_k = solverChoice.sigma_k;

  switch (solverChoice.les_type) {
  case LESType::Smagorinsky:
    // K_LES = 2*mu_t -> extra factor of 0.5 when computing rhoAlpha
    rhoAlpha_r = K_LES(ir, jr, kr) * Pr_or_Sc_turb_inv;
    rhoAlpha_l = K_LES(il, jl, kl) * Pr_or_Sc_turb_inv;
    rhoAlpha += 0.25*(rhoAlpha_l + rhoAlpha_r);
    break;
  case LESType::Deardorff:
    // K_LES = 2*mu_t -> extra factor of 0.5 when computing rhoAlpha
    rhoAlpha_r = K_LES(ir, jr, kr) * Pr_or_Sc_turb_inv;
    rhoAlpha_l = K_LES(il, jl, kl) * Pr_or_Sc_turb_inv;
    rhoAlpha += 0.25*(rhoAlpha_l + rhoAlpha_r) / l_sigma_k;
    break;
  case LESType::None:
    break;
  default:
    amrex::Abort("Error:  LES model is unrecognized");
  }

  // Compute the flux
  amrex::Real diffusionFlux = rhoAlpha * invCellWidth *
      (cell_prim(ir, jr, kr, prim_index) - cell_prim(il, jl, kl, prim_index));

  return diffusionFlux;
}

AMREX_GPU_DEVICE
Real
DiffusionContributionForState(const int &i, const int &j, const int &k,
                              const Array4<const Real>& cell_data,
                              const Array4<const Real>& cell_prim, const int & qty_index,
                              const Array4<Real>& xflux, const Array4<Real>& yflux, const Array4<Real>& zflux,
                              const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                              const Array4<Real>& K_LES,
                              const SolverChoice &solverChoice)
{
  const amrex::Real dx_inv = cellSizeInv[0];
  const amrex::Real dy_inv = cellSizeInv[1];
  const amrex::Real dz_inv = cellSizeInv[2];

  const int prim_index = qty_index - RhoTheta_comp;

  // TODO : could be more efficient to compute and save all fluxes before taking divergence (now all fluxes are computed 2x);
  xflux(i+1,j,k,qty_index) = ComputeDiffusionFluxForState(i+1, j, k, cell_data, cell_prim, prim_index, dx_inv, K_LES, solverChoice, Coord::x);
  xflux(i  ,j,k,qty_index) = ComputeDiffusionFluxForState(i  , j, k, cell_data, cell_prim, prim_index, dx_inv, K_LES, solverChoice, Coord::x);

  yflux(i,j+1,k,qty_index) = ComputeDiffusionFluxForState(i, j+1, k, cell_data, cell_prim, prim_index, dy_inv, K_LES, solverChoice, Coord::y);
  yflux(i,j  ,k,qty_index) = ComputeDiffusionFluxForState(i, j  , k, cell_data, cell_prim, prim_index, dy_inv, K_LES, solverChoice, Coord::y);

  zflux(i,j,k+1,qty_index) = ComputeDiffusionFluxForState(i, j, k+1, cell_data, cell_prim, prim_index, dz_inv, K_LES, solverChoice, Coord::z);
  zflux(i,j,k  ,qty_index) = ComputeDiffusionFluxForState(i, j, k  , cell_data, cell_prim, prim_index, dz_inv, K_LES, solverChoice, Coord::z);

  Real diffusionContribution =
      (xflux(i+1,j,k,qty_index) - xflux(i  ,j,k,qty_index)) * dx_inv   // Diffusive flux in x-dir
     +(yflux(i,j+1,k,qty_index) - yflux(i,j  ,k,qty_index)) * dy_inv   // Diffusive flux in y-dir
     +(zflux(i,j,k+1,qty_index) - zflux(i,j,k  ,qty_index)) * dz_inv;  // Diffusive flux in z-dir

  return diffusionContribution;
}
