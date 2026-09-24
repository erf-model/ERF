#ifndef ERF_NATIVE_SCALAR_DIFFUSION_H_
#define ERF_NATIVE_SCALAR_DIFFUSION_H_

#include "ERF_ScalarDiffusion.H"
#include "ERF_IndexDefines.H"

/** Native-only state/layout decisions resolved before entering a spatial
 * kernel. This policy belongs to native conserved-state adapters and is not
 * used by the reusable field-level numerical operations. */
struct NativeScalarDiffusionPolicy
{
  ScalarDiffusionCoefficients coefficients;
  int scalar_comp = 0;
  int rhs_comp = 0;
  int flux_comp = 0;
  int bc_comp = 0;
  bool is_theta = false;
  bool is_q1 = false;
  bool is_q2 = false;
  bool scale_raw_vertical_flux = false;
};

/** Resolve conserved-state identity to ERF's native coefficient category. */
inline NativeScalarDiffusionPolicy
ResolveNativeScalarDiffusionPolicy(const int qty_comp, const DiffChoice& choice)
{
  NativeScalarDiffusionPolicy policy;
  policy.scalar_comp = qty_comp - 1;
  policy.rhs_comp = qty_comp;
  policy.flux_comp = 0;
  policy.bc_comp =
    (qty_comp >= RhoScalar_comp && qty_comp < RhoScalar_comp + NSCALARS)
      ? BCVars::RhoScalar_bc_comp
      : qty_comp;
  if (policy.bc_comp > BCVars::RhoScalar_bc_comp) {
    policy.bc_comp -= (NSCALARS - 1);
  }
  policy.is_theta = qty_comp == RhoTheta_comp;
  policy.is_q1 = qty_comp == RhoQ1_comp;
  policy.is_q2 = qty_comp == RhoQ2_comp;
  policy.scale_raw_vertical_flux =
    policy.is_theta || qty_comp == RhoKE_comp || policy.is_q1;

  if (policy.is_theta) {
    policy.coefficients.eddy_h_comp = EddyDiff::Theta_h;
    policy.coefficients.eddy_v_comp = EddyDiff::Theta_v;
    if (choice.molec_diff_type == MolecDiffType::ConstantAlpha) {
      policy.coefficients.molecular_coeff = choice.alpha_T;
    } else if (choice.molec_diff_type == MolecDiffType::Constant) {
      policy.coefficients.molecular_coeff = choice.rhoAlpha_T;
    }
  } else if (qty_comp == RhoKE_comp) {
    policy.coefficients.eddy_h_comp = EddyDiff::KE_h;
    policy.coefficients.eddy_v_comp = EddyDiff::KE_v;
  } else if (
    qty_comp >= RhoScalar_comp && qty_comp < RhoScalar_comp + NSCALARS) {
    policy.coefficients.eddy_h_comp = EddyDiff::Scalar_h;
    policy.coefficients.eddy_v_comp = EddyDiff::Scalar_v;
    if (choice.molec_diff_type == MolecDiffType::ConstantAlpha) {
      policy.coefficients.molecular_coeff = choice.alpha_C;
    } else if (choice.molec_diff_type == MolecDiffType::Constant) {
      policy.coefficients.molecular_coeff = choice.rhoAlpha_C;
    }
  } else {
    // Every q-state, including future appended species, uses the Q eddy
    // category. Molecular diffusion is defined only through Q6.
    policy.coefficients.eddy_h_comp = EddyDiff::Q_h;
    policy.coefficients.eddy_v_comp = EddyDiff::Q_v;
    if (qty_comp >= RhoQ1_comp && qty_comp <= RhoQ6_comp) {
      if (choice.molec_diff_type == MolecDiffType::ConstantAlpha) {
        policy.coefficients.molecular_coeff = choice.alpha_C;
      } else if (choice.molec_diff_type == MolecDiffType::Constant) {
        policy.coefficients.molecular_coeff = choice.rhoAlpha_C;
      }
    }
  }
  return policy;
}

#endif // ERF_NATIVE_SCALAR_DIFFUSION_H_
