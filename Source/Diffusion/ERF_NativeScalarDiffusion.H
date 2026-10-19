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
    int rhs_comp    = 0;
    int flux_comp   = 0;
    int bc_comp     = 0;
    bool is_theta   = false;
    bool is_moist   = false;
    bool is_q1      = false;
    bool is_q2      = false;
    bool scale_raw_vertical_flux = false;
};

/** Resolve conserved-state identity to ERF's native coefficient category. */
inline NativeScalarDiffusionPolicy
ResolveNativeScalarDiffusionPolicy (const int qty_comp, const DiffChoice& choice)
{
    NativeScalarDiffusionPolicy policy;
    policy.scalar_comp = qty_comp - 1;
    policy.rhs_comp    = qty_comp;
    policy.flux_comp   = 0;
    policy.bc_comp     = cons_bc_comp(qty_comp);
    policy.is_theta = (qty_comp == RhoTheta_comp);
    // Only the water species are in the implicit moisture solve; non-water species
    // packed after them keep a fully explicit vertical flux
    policy.is_moist = (qty_comp >= RhoQ1_comp && qty_comp < (RhoQ1_comp + choice.n_qstate_moist));
    policy.is_q1    = (qty_comp == RhoQ1_comp);
    policy.is_q2    = (qty_comp == RhoQ2_comp);
    policy.scale_raw_vertical_flux =
        (policy.is_theta || qty_comp == RhoKE_comp || policy.is_moist);

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
    } else if (qty_comp >= RhoScalar_comp && qty_comp < RhoScalar_comp + NSCALARS) {
        policy.coefficients.eddy_h_comp = EddyDiff::Scalar_h;
        policy.coefficients.eddy_v_comp = EddyDiff::Scalar_v;
        if (choice.molec_diff_type == MolecDiffType::ConstantAlpha) {
            policy.coefficients.molecular_coeff = choice.alpha_C;
        } else if (choice.molec_diff_type == MolecDiffType::Constant) {
            policy.coefficients.molecular_coeff = choice.rhoAlpha_C;
        }
    } else {
        // Every q-state, including future appended species, uses the Q eddy
        // category. Molecular diffusion is defined for the moisture species
        // (number concentrations included), not for appended non-moist species.
        policy.coefficients.eddy_h_comp = EddyDiff::Q_h;
        policy.coefficients.eddy_v_comp = EddyDiff::Q_v;
        if (policy.is_moist) {
            if (choice.molec_diff_type == MolecDiffType::ConstantAlpha) {
                policy.coefficients.molecular_coeff = choice.alpha_C;
            } else if (choice.molec_diff_type == MolecDiffType::Constant) {
                policy.coefficients.molecular_coeff = choice.rhoAlpha_C;
            }
        }
    }
    return policy;
}

/** Molecular diffusivity by primitive component for the implicit vertical
 * solve and the EB diffusion (ERF_SetupVertDiff.H): theta, every passive scalar
 * (all NSCALARS of them) and the moisture species (number concentrations
 * included); everything else -- KE, and any non-moist species a scheme appends
 * after its moisture species -- is left at zero. With erf.molec_diff_type =
 * Constant or ConstantAlpha these are the coefficients
 * ResolveNativeScalarDiffusionPolicy() gives the explicit paths. */
inline amrex::Vector<amrex::Real>
VerticalMolecularAlphaByPrim (const DiffChoice& choice)
{
    amrex::Vector<amrex::Real> alpha_eff(NPRIMVAR_max, amrex::Real(0));
    const bool consA = (choice.molec_diff_type == MolecDiffType::ConstantAlpha);
    const amrex::Real alpha_T_eff = consA ? choice.alpha_T : choice.rhoAlpha_T;
    const amrex::Real alpha_C_eff = consA ? choice.alpha_C : choice.rhoAlpha_C;
    alpha_eff[PrimTheta_comp] = alpha_T_eff;
    for (int n = 0; n < NSCALARS; ++n) {
        alpha_eff[PrimScalar_comp + n] = alpha_C_eff;
    }
    for (int i = PrimQ1_comp; i < PrimQ1_comp + choice.n_qstate_moist; ++i) {
        alpha_eff[i] = alpha_C_eff;
    }
    return alpha_eff;
}

#endif // ERF_NATIVE_SCALAR_DIFFUSION_H_
