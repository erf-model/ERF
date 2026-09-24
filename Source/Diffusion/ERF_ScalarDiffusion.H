#ifndef ERF_SCALAR_DIFFUSION_H_
#define ERF_SCALAR_DIFFUSION_H_

#include <type_traits>

#include "AMReX_Box.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_Gpu.H"
#include "ERF_DiffStruct.H"
#include "ERF_IndexDefines.H"
#include "ERF_TerrainMetrics.H"

/**
 * Coefficients and state-layout mapping at the native ERF adapter boundary.
 *
 * The numerical routines below receive the scalar, density, flux, and RHS
 * components independently.  `molecular_coeff` is kinematic alpha for
 * ConstantAlpha and density-weighted rhoAlpha for Constant.  The density
 * multiplication is selected by the host-side template specialization.
 */
struct ScalarDiffusionCoefficients
{
  amrex::Real molecular_coeff = amrex::Real(0.0);
  int eddy_h_comp = 0;
  int eddy_v_comp = 0;
};

static_assert(std::is_trivially_copyable_v<ScalarDiffusionCoefficients>);

/** Native-only state/layout decisions resolved before entering a spatial
 * kernel. */
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

/**
 * Resolve conserved-state identity to the coefficient category used by ERF.
 * This routine belongs to the native adapter; reusable spatial primitives do
 * not inspect conserved-state indices.
 */
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
    // category.  Molecular diffusion is defined only through Q6.
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

/**
 * Density-weighted face transport coefficient.  The booleans are compile-time
 * choices for ConstantAlpha/Constant/None and turbulent/no-turbulence modes;
 * turbulent coefficients are already density weighted and are never
 * multiplied by rho here.
 *
 * `rho_comp` and the eddy component are independent of the scalar component.
 */
template <bool MultiplyMolecularByDensity, bool AddTurbulence>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
ScalarDiffusionFaceCoefficient(
  const amrex::Array4<const amrex::Real>& rho,
  const int rho_comp,
  const amrex::Array4<const amrex::Real>& mu_turb,
  const ScalarDiffusionCoefficients& coefficients,
  const int i,
  const int j,
  const int k,
  const int di,
  const int dj,
  const int dk,
  const int eddy_comp) noexcept
{
  amrex::Real coeff = coefficients.molecular_coeff;
  if constexpr (MultiplyMolecularByDensity) {
    coeff *= amrex::Real(0.5) *
             (rho(i, j, k, rho_comp) + rho(i - di, j - dj, k - dk, rho_comp));
  }
  if constexpr (AddTurbulence) {
    coeff += amrex::Real(0.5) * (mu_turb(i, j, k, eddy_comp) +
                                 mu_turb(i - di, j - dj, k - dk, eddy_comp));
  }
  return coeff;
}

/** Face-value one-sided derivative weights for an S-grid vertical boundary. */
struct StretchedDiffusionBoundaryWeights
{
  amrex::Real c1, c2, c3;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE StretchedDiffusionBoundaryWeights
StretchedDiffusionDirichletWeights(
  const amrex::Real dz_near, const amrex::Real dz_next) noexcept
{
  const amrex::Real f = dz_next / dz_near + amrex::Real(2.0);
  const amrex::Real f2 = f * f;
  const amrex::Real c3 = amrex::Real(2.0) / (f - f2);
  return {-(amrex::Real(1.0) - f2) * c3, -f2 * c3, c3};
}

/**
 * Compute the existing N-grid signed face flux, F^D=-K grad(chi). `scalar_comp`
 * selects the intensive scalar and `flux_comp` is applied by the caller to its
 * independent face-flux view.  `metric` is mf_x/mf_y for N horizontal faces,
 * or one vertically.  Physical Dirichlet flags select the legacy one-sided
 * stencil; ext_dir_upwind selection stays in the adapter.
 */
template <int Direction>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
ScalarDiffusionFlux_N(
  const amrex::Array4<const amrex::Real>& scalar,
  const int scalar_comp,
  const int i,
  const int j,
  const int k,
  const amrex::Real coeff,
  const amrex::Real inv_spacing,
  const amrex::Real metric,
  const bool ext_dir_lo,
  const bool ext_dir_hi) noexcept
{
  static_assert(Direction >= 0 && Direction < 3);
  amrex::Real grad;
  if constexpr (Direction == 0) {
    if (ext_dir_lo) {
      grad = -(amrex::Real(8.0) / amrex::Real(3.0)) *
               scalar(i - 1, j, k, scalar_comp) +
             amrex::Real(3.0) * scalar(i, j, k, scalar_comp) -
             (amrex::Real(1.0) / amrex::Real(3.0)) *
               scalar(i + 1, j, k, scalar_comp);
    } else if (ext_dir_hi) {
      grad =
        (amrex::Real(8.0) / amrex::Real(3.0)) * scalar(i, j, k, scalar_comp) -
        amrex::Real(3.0) * scalar(i - 1, j, k, scalar_comp) +
        (amrex::Real(1.0) / amrex::Real(3.0)) *
          scalar(i - 2, j, k, scalar_comp);
    } else {
      grad = scalar(i, j, k, scalar_comp) - scalar(i - 1, j, k, scalar_comp);
    }
  } else if constexpr (Direction == 1) {
    if (ext_dir_lo) {
      grad = -(amrex::Real(8.0) / amrex::Real(3.0)) *
               scalar(i, j - 1, k, scalar_comp) +
             amrex::Real(3.0) * scalar(i, j, k, scalar_comp) -
             (amrex::Real(1.0) / amrex::Real(3.0)) *
               scalar(i, j + 1, k, scalar_comp);
    } else if (ext_dir_hi) {
      grad =
        (amrex::Real(8.0) / amrex::Real(3.0)) * scalar(i, j, k, scalar_comp) -
        amrex::Real(3.0) * scalar(i, j - 1, k, scalar_comp) +
        (amrex::Real(1.0) / amrex::Real(3.0)) *
          scalar(i, j - 2, k, scalar_comp);
    } else {
      grad = scalar(i, j, k, scalar_comp) - scalar(i, j - 1, k, scalar_comp);
    }
  } else {
    if (ext_dir_lo) {
      grad = -(amrex::Real(8.0) / amrex::Real(3.0)) *
               scalar(i, j, k - 1, scalar_comp) +
             amrex::Real(3.0) * scalar(i, j, k, scalar_comp) -
             (amrex::Real(1.0) / amrex::Real(3.0)) *
               scalar(i, j, k + 1, scalar_comp);
    } else if (ext_dir_hi) {
      grad =
        (amrex::Real(8.0) / amrex::Real(3.0)) * scalar(i, j, k, scalar_comp) -
        amrex::Real(3.0) * scalar(i, j, k - 1, scalar_comp) +
        (amrex::Real(1.0) / amrex::Real(3.0)) *
          scalar(i, j, k - 2, scalar_comp);
    } else {
      grad = scalar(i, j, k, scalar_comp) - scalar(i, j, k - 1, scalar_comp);
    }
  }
  return -coeff * grad * inv_spacing * metric;
}

/** Existing S-grid vertical-center gradient, including boundary and tile-edge
 * spacing. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
StretchedScalarGradient(
  const amrex::Array4<const amrex::Real>& scalar,
  const int scalar_comp,
  const int i,
  const int j,
  const int k,
  const amrex::Real* dz,
  const int klo,
  const int khi,
  const bool ext_dir_lo,
  const bool ext_dir_hi) noexcept
{
  if (ext_dir_lo) {
    const auto w = StretchedDiffusionDirichletWeights(dz[k], dz[k + 1]);
    return (amrex::Real(1.0) / dz[k]) *
           (w.c1 * scalar(i, j, k - 1, scalar_comp) +
            w.c2 * scalar(i, j, k, scalar_comp) +
            w.c3 * scalar(i, j, k + 1, scalar_comp));
  }
  if (ext_dir_hi) {
    const auto w = StretchedDiffusionDirichletWeights(dz[k - 1], dz[k - 2]);
    return (amrex::Real(1.0) / dz[k - 1]) *
           (-(
             w.c1 * scalar(i, j, k, scalar_comp) +
             w.c2 * scalar(i, j, k - 1, scalar_comp) +
             w.c3 * scalar(i, j, k - 2, scalar_comp)));
  }
  const amrex::Real inv = (k == klo) ? amrex::Real(1.0) / dz[k]
                          : (k == khi + 1)
                            ? amrex::Real(1.0) / dz[k - 1]
                            : amrex::Real(2.0) / (dz[k] + dz[k - 1]);
  return inv *
         (scalar(i, j, k, scalar_comp) - scalar(i, j, k - 1, scalar_comp));
}

/** Existing S-grid horizontal raw flux before its separate map-factor
 * adjustment kernel. */
template <int Direction>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
ScalarDiffusionFlux_S_Horizontal(
  const amrex::Array4<const amrex::Real>& scalar,
  const int scalar_comp,
  const int i,
  const int j,
  const int k,
  const amrex::Real coeff,
  const amrex::Real inv_spacing,
  const amrex::Real face_map) noexcept
{
  static_assert(Direction == 0 || Direction == 1);
  const amrex::Real delta =
    Direction == 0
      ? scalar(i, j, k, scalar_comp) - scalar(i - 1, j, k, scalar_comp)
      : scalar(i, j, k, scalar_comp) - scalar(i, j - 1, k, scalar_comp);
  const amrex::Real gradient = inv_spacing * delta;
  return -coeff * face_map * gradient;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
ScalarDiffusionFlux_Tx(
  const amrex::Real coeff,
  const amrex::Real face_map,
  const amrex::Real grad_xi,
  const amrex::Real h_xi,
  const amrex::Real grad_z) noexcept
{
  return -coeff * face_map * (grad_xi - h_xi * grad_z);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
ScalarDiffusionFlux_Ty(
  const amrex::Real coeff,
  const amrex::Real face_map,
  const amrex::Real grad_eta,
  const amrex::Real h_eta,
  const amrex::Real grad_z) noexcept
{
  return -coeff * face_map * (grad_eta - h_eta * grad_z);
}

/**
 * Raw terrain-coordinate F_z and its conservative transformed transfer are
 * distinct quantities. G_z is the numerator
 * F_z - m_x h_xi Fbar_x - m_y h_eta Fbar_y; T_z=G_z/(m_x m_y).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
TerrainDiffusionGz(
  const amrex::Real raw_Fz,
  const amrex::Real mx,
  const amrex::Real h_xi,
  const amrex::Real Fbar_x,
  const amrex::Real my,
  const amrex::Real h_eta,
  const amrex::Real Fbar_y) noexcept
{
  return raw_Fz - h_xi * mx * Fbar_x - h_eta * my * Fbar_y;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
TerrainDiffusionMappedTz(
  const amrex::Real Gz, const amrex::Real mx, const amrex::Real my) noexcept
{
  return Gz / (mx * my);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
TerrainDiffusionMappedTx(
  const amrex::Real raw_Fx,
  const amrex::Real h_zeta,
  const amrex::Real my) noexcept
{
  return raw_Fx * h_zeta / my;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
TerrainDiffusionMappedTy(
  const amrex::Real raw_Fy,
  const amrex::Real h_zeta,
  const amrex::Real mx) noexcept
{
  return raw_Fy * h_zeta / mx;
}

/**
 * T-grid cell divergence numerator in ERF's economical form. `Gz` is the
 * terrain numerator (not raw F_z and not mapped T_z); `h_zeta_cell` is stored
 * in ERF's detJ view. The caller subtracts the returned flux divergence from
 * `rhs[rhs_comp]`.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
TerrainDiffusionDivergence_T(
  const amrex::Real fx_hi,
  const amrex::Real ax_hi,
  const amrex::Real my_hi,
  const amrex::Real fx_lo,
  const amrex::Real ax_lo,
  const amrex::Real my_lo,
  const amrex::Real fy_hi,
  const amrex::Real ay_hi,
  const amrex::Real mx_hi,
  const amrex::Real fy_lo,
  const amrex::Real ay_lo,
  const amrex::Real mx_lo,
  const amrex::Real Gz_hi,
  const amrex::Real Gz_lo,
  const amrex::Real dx_inv,
  const amrex::Real dy_inv,
  const amrex::Real dz_inv,
  const amrex::Real mx_cell,
  const amrex::Real my_cell,
  const amrex::Real h_zeta_cell) noexcept
{
  const amrex::Real mfsq = mx_cell * my_cell;
  amrex::Real divergence = (TerrainDiffusionMappedTx(fx_hi, ax_hi, my_hi) -
                            TerrainDiffusionMappedTx(fx_lo, ax_lo, my_lo)) *
                             dx_inv * mfsq +
                           (TerrainDiffusionMappedTy(fy_hi, ay_hi, mx_hi) -
                            TerrainDiffusionMappedTy(fy_lo, ay_lo, mx_lo)) *
                             dy_inv * mfsq +
                           (Gz_hi - Gz_lo) * dz_inv;
  divergence /= h_zeta_cell;
  return divergence;
}

/**
 * Apply the conservative N-grid divergence to a separately selected RHS
 * component. Face scratch is component-selected independently in each view.
 */
inline void
ApplyScalarDiffusionFluxDivergence_N(
  const amrex::Box& bx,
  const amrex::Array4<const amrex::Real>& xflux,
  const amrex::Array4<const amrex::Real>& yflux,
  const amrex::Array4<const amrex::Real>& raw_Fz,
  const int flux_comp,
  const amrex::Array4<amrex::Real>& rhs,
  const int rhs_comp,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
  const amrex::Array4<const amrex::Real>& mf_mx,
  const amrex::Array4<const amrex::Real>& mf_my)
{
  const amrex::Real dx_inv = cellSizeInv[0];
  const amrex::Real dy_inv = cellSizeInv[1];
  const amrex::Real dz_inv = cellSizeInv[2];
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real mfsq = mf_mx(i, j, 0) * mf_my(i, j, 0);
    rhs(i, j, k, rhs_comp) -=
      (xflux(i + 1, j, k, flux_comp) - xflux(i, j, k, flux_comp)) * dx_inv *
        mfsq +
      (yflux(i, j + 1, k, flux_comp) - yflux(i, j, k, flux_comp)) * dy_inv *
        mfsq +
      (raw_Fz(i, j, k + 1, flux_comp) - raw_Fz(i, j, k, flux_comp)) * dz_inv;
  });
}

/** S-grid conservative divergence; vertical cell width is the receiver's dz[k].
 */
inline void
ApplyScalarDiffusionFluxDivergence_S(
  const amrex::Box& bx,
  const amrex::Array4<const amrex::Real>& xflux,
  const amrex::Array4<const amrex::Real>& yflux,
  const amrex::Array4<const amrex::Real>& raw_Fz,
  const int flux_comp,
  const amrex::Array4<amrex::Real>& rhs,
  const int rhs_comp,
  const amrex::Array4<const amrex::Real>& mf_mx,
  const amrex::Array4<const amrex::Real>& mf_my,
  const amrex::Real dx_inv,
  const amrex::Real dy_inv,
  const amrex::Real* dz)
{
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real mfsq = mf_mx(i, j, 0) * mf_my(i, j, 0);
    const amrex::Real dzk_inv = amrex::Real(1.0) / dz[k];
    const amrex::Real stateContrib =
      (xflux(i + 1, j, k, flux_comp) - xflux(i, j, k, flux_comp)) * dx_inv *
        mfsq +
      (yflux(i, j + 1, k, flux_comp) - yflux(i, j, k, flux_comp)) * dy_inv *
        mfsq +
      (raw_Fz(i, j, k + 1, flux_comp) - raw_Fz(i, j, k, flux_comp)) * dzk_inv;
    rhs(i, j, k, rhs_comp) -= stateContrib;
  });
}

/** Scale only raw F_z for the explicit integration fraction. */
inline void
ScaleScalarDiffusionVerticalFlux(
  const amrex::Box& zbx,
  const amrex::Array4<amrex::Real>& raw_Fz,
  const int flux_comp,
  const amrex::Real explicit_fac)
{
  amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    raw_Fz(i, j, k, flux_comp) *= explicit_fac;
  });
}

#endif // ERF_SCALAR_DIFFUSION_H_
