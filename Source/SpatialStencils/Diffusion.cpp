#include <SpatialStencils.H>
#include <StressTerm.H>

using namespace amrex;

AMREX_GPU_DEVICE
Real
DiffusionContributionForMom(const int &i, const int &j, const int &k,
                            const Array4<const Real>& u, const Array4<const Real>& v, const Array4<const Real>& w,
                            const Array4<const Real>& cons,
                            const enum MomentumEqn &momentumEqn,
                            const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                            const Array4<Real>& K_turb,
                            const SolverChoice &solverChoice,
                            const Box& domain, const amrex::BCRec* bc_ptr)
{
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    Real diffContrib = 0.0;

    int l_spatial_order = solverChoice.spatial_order;

    switch (momentumEqn) {
        case MomentumEqn::x:
            Real tau11Next, tau11Prev, tau12Next, tau12Prev, tau13Next, tau13Prev;
            tau11Next = ComputeStressTerm(i+1, j, k, u, v, w, momentumEqn,
                                          DiffusionDir::x, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau11Prev = ComputeStressTerm(i  , j, k, u, v, w, momentumEqn,
                                          DiffusionDir::x, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau12Next = ComputeStressTerm(i, j+1, k, u, v, w, momentumEqn,
                                          DiffusionDir::y, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau12Prev = ComputeStressTerm(i, j  , k, u, v, w, momentumEqn,
                                          DiffusionDir::y, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau13Next = ComputeStressTerm(i, j, k+1, u, v, w, momentumEqn,
                                          DiffusionDir::z, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau13Prev = ComputeStressTerm(i, j, k  , u, v, w, momentumEqn,
                                          DiffusionDir::z, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);

            diffContrib = (tau11Next - tau11Prev) * dxInv  // Contribution to x-mom eqn from diffusive flux in x-dir
                        + (tau12Next - tau12Prev) * dyInv  // Contribution to x-mom eqn from diffusive flux in y-dir
                        + (tau13Next - tau13Prev) * dzInv; // Contribution to x-mom eqn from diffusive flux in z-dir
            if (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha)
            {
                diffContrib *= InterpolateFromCellOrFace(i, j, k, cons, Rho_comp, u(i,j,k), Coord::x, l_spatial_order) /
                               solverChoice.rho0_trans;
            }
            break;
        case MomentumEqn::y:
            Real tau21Next, tau21Prev, tau22Next, tau22Prev, tau23Next, tau23Prev;
            tau21Next = ComputeStressTerm(i+1, j, k, u, v, w, momentumEqn,
                                          DiffusionDir::x, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau21Prev = ComputeStressTerm(i  , j, k, u, v, w, momentumEqn,
                                          DiffusionDir::x, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau22Next = ComputeStressTerm(i, j+1, k, u, v, w, momentumEqn,
                                          DiffusionDir::y, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau22Prev = ComputeStressTerm(i, j  , k, u, v, w, momentumEqn,
                                          DiffusionDir::y, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau23Next = ComputeStressTerm(i, j, k+1, u, v, w, momentumEqn,
                                          DiffusionDir::z, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau23Prev = ComputeStressTerm(i, j, k  , u, v, w, momentumEqn,
                                          DiffusionDir::z, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);

            diffContrib = (tau21Next - tau21Prev) * dxInv  // Contribution to y-mom eqn from diffusive flux in x-dir
                        + (tau22Next - tau22Prev) * dyInv  // Contribution to y-mom eqn from diffusive flux in y-dir
                        + (tau23Next - tau23Prev) * dzInv; // Contribution to y-mom eqn from diffusive flux in z-dir
            if (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha)
            {
                diffContrib *= InterpolateFromCellOrFace(i, j, k, cons, Rho_comp, v(i,j,k), Coord::y, l_spatial_order) /
                               solverChoice.rho0_trans;
            }
            break;
        case MomentumEqn::z:
            Real tau31Next, tau31Prev, tau32Next, tau32Prev, tau33Next, tau33Prev;
            tau31Next = ComputeStressTerm(i+1, j, k, u, v, w, momentumEqn,
                                          DiffusionDir::x, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau31Prev = ComputeStressTerm(i  , j, k, u, v, w, momentumEqn,
                                          DiffusionDir::x, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau32Next = ComputeStressTerm(i, j+1, k, u, v, w, momentumEqn,
                                          DiffusionDir::y, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau32Prev = ComputeStressTerm(i, j  , k, u, v, w, momentumEqn,
                                          DiffusionDir::y, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau33Next = ComputeStressTerm(i, j, k+1, u, v, w, momentumEqn,
                                          DiffusionDir::z, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau33Prev = ComputeStressTerm(i, j, k  , u, v, w, momentumEqn,
                                          DiffusionDir::z, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);

            diffContrib = (tau31Next - tau31Prev) * dxInv  // Contribution to z-mom eqn from diffusive flux in x-dir
                        + (tau32Next - tau32Prev) * dyInv  // Contribution to z-mom eqn from diffusive flux in y-dir
                        + (tau33Next - tau33Prev) * dzInv; // Contribution to z-mom eqn from diffusive flux in z-dir
            if (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha)
            {
                diffContrib *= InterpolateFromCellOrFace(i, j, k, cons, Rho_comp, w(i,j,k), Coord::z, l_spatial_order) /
                               solverChoice.rho0_trans;
            }
            break;
        default:
            amrex::Abort("Error: Momentum equation is unrecognized");
    }

    return diffContrib;
}

AMREX_GPU_DEVICE
amrex::Real ComputeDiffusionFluxForState(const int &i, const int &j, const int &k,
                     const Array4<const Real>& cell_data,
                     const Array4<const Real>& cell_prim, const int & prim_index,
                     const amrex::Real invCellWidth,
                     const Array4<Real>& K_turb,
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
  if (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha) {
    rhoFace = (cell_data(il, jl, kl, Rho_comp) + cell_data(ir, jr, kr, Rho_comp)) * 0.5;
  }

  switch(prim_index) {
      case PrimTheta_comp: // Potential Temperature
          if (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha) {
              rhoAlpha_molec = rhoFace * solverChoice.alpha_T;
          } else {
              rhoAlpha_molec = solverChoice.rhoAlpha_T;
          }
          Pr_or_Sc_turb_inv = solverChoice.Pr_t_inv;
          break;

      case PrimKE_comp: // Turbulent KE
          rhoAlpha_molec = 0.;
          Pr_or_Sc_turb_inv = 1.0 / solverChoice.sigma_k;
          break;

      case PrimScalar_comp: // Scalar
          if (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha) {
              rhoAlpha_molec = rhoFace * solverChoice.alpha_C;
          } else {
              rhoAlpha_molec = solverChoice.rhoAlpha_C;
          }
          Pr_or_Sc_turb_inv = solverChoice.Sc_t_inv;
          break;
      default:
        amrex::Abort("Error: Diffusion term for the data index isn't implemented");
  }

  amrex::Real rhoAlpha = 0.0;

  if ( (solverChoice.molec_diff_type == MolecDiffType::Constant) ||
       (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha) ) {
    rhoAlpha += rhoAlpha_molec;
  }

  if ( (solverChoice.les_type == LESType::Smagorinsky) ||
       (solverChoice.les_type == LESType::Deardorff  ) ) {
      // K_turb = 2*mu_t -> extra factor of 0.5 when computing rhoAlpha
      rhoAlpha += 0.25*(K_turb(ir,jr,kr) + K_turb(il,jl,kl)) * Pr_or_Sc_turb_inv;
  } else if (solverChoice.les_type != LESType::None) {
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
                              const Array4<Real>& K_turb,
                              const SolverChoice &solverChoice)
{
  const amrex::Real dx_inv = cellSizeInv[0];
  const amrex::Real dy_inv = cellSizeInv[1];
  const amrex::Real dz_inv = cellSizeInv[2];

  const int prim_index = qty_index - RhoTheta_comp;

  // TODO : could be more efficient to compute and save all fluxes before taking divergence (now all fluxes are computed 2x);
  xflux(i+1,j,k,qty_index) = ComputeDiffusionFluxForState(i+1, j, k, cell_data, cell_prim, prim_index, dx_inv, K_turb, solverChoice, Coord::x);
  xflux(i  ,j,k,qty_index) = ComputeDiffusionFluxForState(i  , j, k, cell_data, cell_prim, prim_index, dx_inv, K_turb, solverChoice, Coord::x);

  yflux(i,j+1,k,qty_index) = ComputeDiffusionFluxForState(i, j+1, k, cell_data, cell_prim, prim_index, dy_inv, K_turb, solverChoice, Coord::y);
  yflux(i,j  ,k,qty_index) = ComputeDiffusionFluxForState(i, j  , k, cell_data, cell_prim, prim_index, dy_inv, K_turb, solverChoice, Coord::y);

  zflux(i,j,k+1,qty_index) = ComputeDiffusionFluxForState(i, j, k+1, cell_data, cell_prim, prim_index, dz_inv, K_turb, solverChoice, Coord::z);
  zflux(i,j,k  ,qty_index) = ComputeDiffusionFluxForState(i, j, k  , cell_data, cell_prim, prim_index, dz_inv, K_turb, solverChoice, Coord::z);

  Real diffusionContribution =
      (xflux(i+1,j,k,qty_index) - xflux(i  ,j,k,qty_index)) * dx_inv   // Diffusive flux in x-dir
     +(yflux(i,j+1,k,qty_index) - yflux(i,j  ,k,qty_index)) * dy_inv   // Diffusive flux in y-dir
     +(zflux(i,j,k+1,qty_index) - zflux(i,j,k  ,qty_index)) * dz_inv;  // Diffusive flux in z-dir

  return diffusionContribution;
}
