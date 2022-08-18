#include <SpatialStencils.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ComputeDiffusionFluxForState (int i, int j, int k,
                              const amrex::Array4<const amrex::Real>& cell_data,
                              const amrex::Array4<const amrex::Real>& cell_prim, const int & prim_index,
                              amrex::Real invCellWidth,
                              const amrex::Array4<amrex::Real>& K_turb,
                              const SolverChoice &solverChoice,
                              enum Coord coordDir)
{
  if ( (solverChoice.molec_diff_type == MolecDiffType::None) &&
       (solverChoice.les_type        ==       LESType::None) &&
       (solverChoice.pbl_type        ==       PBLType::None) ) {
          return 0.;
  } else {
  // Get indices of states to left and right of the face on which we want the flux
  const int il = i - (coordDir == Coord::x);
  const int ir = i;
  const int jl = j - (coordDir == Coord::y);
  const int jr = j;
  const int kl = k - (coordDir == Coord::z);
  const int kr = k;

  // Get diffusion coefficients
  amrex::Real rhoAlpha_molec;
  int eddy_diff_idx;

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
          if (coordDir == Coord::z) {
              eddy_diff_idx = EddyDiff::Theta_v;
          } else {
              eddy_diff_idx = EddyDiff::Theta_h;
          }
          break;

      case PrimKE_comp: // Turbulent KE
          rhoAlpha_molec = 0.;
          if (coordDir == Coord::z) {
              eddy_diff_idx = EddyDiff::KE_v;
          } else {
              eddy_diff_idx = EddyDiff::KE_h;
          }
          break;

      case PrimQKE_comp: // Turbulent QKE
          rhoAlpha_molec = 0.;
          if (coordDir == Coord::z) {
              eddy_diff_idx = EddyDiff::QKE_v;
          } else {
              eddy_diff_idx = EddyDiff::QKE_h;
          }
          break;

      case PrimScalar_comp: // Scalar
          if (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha) {
              rhoAlpha_molec = rhoFace * solverChoice.alpha_C;
          } else {
              rhoAlpha_molec = solverChoice.rhoAlpha_C;
          }
          if (coordDir == Coord::z) {
              eddy_diff_idx = EddyDiff::Scalar_v;
          } else {
              eddy_diff_idx = EddyDiff::Scalar_h;
          }
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
       (solverChoice.les_type == LESType::Deardorff  ) ||
       (solverChoice.pbl_type == PBLType::MYNN25     ) ) {
    rhoAlpha += 0.5*(K_turb(ir,jr,kr,eddy_diff_idx) + K_turb(il,jl,kl,eddy_diff_idx));
  }

  // Compute the flux
  amrex::Real diffusionFlux = rhoAlpha * invCellWidth *
      (cell_prim(ir, jr, kr, prim_index) - cell_prim(il, jl, kl, prim_index));

  return diffusionFlux;
    }
}
