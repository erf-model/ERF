#ifndef ERF_SCALAR_DIFFUSION_H_
#define ERF_SCALAR_DIFFUSION_H_

#include <type_traits>

#include "AMReX_Box.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_Gpu.H"
#include "ERF_DiffStruct.H"
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

/** Scalar, density, RHS and flux views with independent component indices. */
struct ScalarDiffusionFieldViews
{
  amrex::Array4<const amrex::Real> scalar;
  int scalar_comp = 0;
  amrex::Array4<const amrex::Real> density;
  int rho_comp = 0;
  amrex::Array4<const amrex::Real> mu_turb;
  amrex::Array4<amrex::Real> xflux;
  amrex::Array4<amrex::Real> yflux;
  amrex::Array4<amrex::Real> zflux;
  int flux_comp = 0;
  amrex::Array4<amrex::Real> rhs;
  int rhs_comp = 0;
};

/** Host-side coefficient form; dispatches to compile-time face kernels. */
struct ScalarDiffusionCoefficientMode
{
  bool molecular_is_kinematic = false;
  bool add_turbulence = false;
};

/** Explicit face fluxes supplied by an application-specific boundary model. */
struct ScalarDiffusionPrescribedFluxes
{
  bool use_x_on_side = false;
  bool use_y_on_side = false;
  bool use_z_on_side = false;
  bool use_x_on_bottom_z = false;
  bool use_y_on_bottom_z = false;
  amrex::Array4<const amrex::Real> x;
  amrex::Array4<const amrex::Real> y;
  amrex::Array4<const amrex::Real> z;
  int x_comp = 0;
  int y_comp = 0;
  int z_comp = 0;
};

/** Optional raw vertical-flux diagnostic, written by the z-face kernel. */
struct ScalarDiffusionRawFluxDiagnostic
{
  bool enabled = false;
  bool write_on_surface = false;
  amrex::Array4<amrex::Real> z;
  int z_comp = 0;
};

/** Surface-layer faces that may select prescribed fluxes or zero flux. */
struct ScalarDiffusionSurfaceFaces
{
  bool xlo = false;
  bool xhi = false;
  bool ylo = false;
  bool yhi = false;
  bool zlo = false;
  bool zhi = false;
};

/** Configuration shared by component-explicit face builders. */
struct ScalarDiffusionFluxPolicy
{
  ScalarDiffusionCoefficients coefficients;
  ScalarDiffusionCoefficientMode coefficient_mode;
  ScalarDiffusionPrescribedFluxes prescribed;
  ScalarDiffusionRawFluxDiagnostic diagnostic;
  ScalarDiffusionSurfaceFaces surface;
};

static_assert(std::is_trivially_copyable_v<ScalarDiffusionFieldViews>);
static_assert(std::is_trivially_copyable_v<ScalarDiffusionFluxPolicy>);

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
 * Component-explicit N-grid face construction. Boundary-condition selection,
 * map factors, prescribed surface fluxes and raw-flux diagnostics are explicit
 * inputs; this operation does not inspect native conserved-state identity.
 */
template <bool MultiplyMolecularByDensity, bool AddTurbulence>
inline void
BuildScalarDiffusionFluxes_N_Impl(
  const amrex::Box& bx,
  const amrex::Box& domain,
  const ScalarDiffusionFieldViews& field,
  const ScalarDiffusionFluxPolicy& policy,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
  const amrex::Array4<const amrex::Real>& u,
  const amrex::Array4<const amrex::Real>& v,
  const amrex::Array4<const amrex::Real>& mf_ux,
  const amrex::Array4<const amrex::Real>& mf_uy,
  const amrex::Array4<const amrex::Real>& mf_vx,
  const amrex::Array4<const amrex::Real>& mf_vy,
  const amrex::BCRec* bc_ptr,
  const int bc_comp)
{
  using namespace amrex;
  const auto dom_lo = lbound(domain);
  const auto dom_hi = ubound(domain);
  const Box xbx = surroundingNodes(bx, 0);
  const Box ybx = surroundingNodes(bx, 1);
  const Box zbx = surroundingNodes(bx, 2);
  const Real dx_inv = cellSizeInv[0];
  const Real dy_inv = cellSizeInv[1];
  const Real dz_inv = cellSizeInv[2];
  const auto coefficients = policy.coefficients;
  const auto prescribed = policy.prescribed;
  const auto diagnostic = policy.diagnostic;
  const auto surface = policy.surface;
  const int scalar_comp = field.scalar_comp;
  const int rho_comp = field.rho_comp;
  const int flux_comp = field.flux_comp;
  const auto scalar = field.scalar;
  const auto density = field.density;
  const auto mu_turb = field.mu_turb;
  const auto xflux = field.xflux;
  const auto yflux = field.yflux;
  const auto zflux = field.zflux;

  ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real coeff = ScalarDiffusionFaceCoefficient<
      MultiplyMolecularByDensity, AddTurbulence>(
      density, rho_comp, mu_turb, coefficients, i, j, k, 1, 0, 0,
      coefficients.eddy_h_comp);
    const bool ext_lo =
      ((bc_ptr[bc_comp].lo(0) == ERFBCType::ext_dir) ||
       (bc_ptr[bc_comp].lo(0) == ERFBCType::ext_dir_prim) ||
       (bc_ptr[bc_comp].lo(0) == ERFBCType::ext_dir_upwind &&
        u(dom_lo.x, j, k) >= Real(0.0))) &&
      i == dom_lo.x;
    const bool ext_hi =
      ((bc_ptr[bc_comp].hi(0) == ERFBCType::ext_dir) ||
       (bc_ptr[bc_comp].hi(0) == ERFBCType::ext_dir_prim) ||
       (bc_ptr[bc_comp].hi(0) == ERFBCType::ext_dir_upwind &&
        u(dom_hi.x + 1, j, k) <= Real(0.0))) &&
      i == dom_hi.x + 1;
    const bool on_surface =
      (surface.xlo && i == dom_lo.x) ||
      (surface.xhi && i == dom_hi.x + 1);

    if (ext_lo) {
      xflux(i, j, k, flux_comp) = ScalarDiffusionFlux_N<0>(
        scalar, scalar_comp, i, j, k, coeff, dx_inv,
        mf_ux(i, j, 0) / mf_uy(i, j, 0), true, false);
    } else if (ext_hi) {
      xflux(i, j, k, flux_comp) = ScalarDiffusionFlux_N<0>(
        scalar, scalar_comp, i, j, k, coeff, dx_inv,
        mf_ux(i, j, 0) / mf_uy(i, j, 0), false, true);
    } else if (on_surface && prescribed.use_x_on_side) {
      xflux(i, j, k, flux_comp) = prescribed.x(i, j, k, prescribed.x_comp);
    } else {
      xflux(i, j, k, flux_comp) = ScalarDiffusionFlux_N<0>(
        scalar, scalar_comp, i, j, k, coeff, dx_inv,
        mf_ux(i, j, 0) / mf_uy(i, j, 0), false, false);
    }
  });

  ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real coeff = ScalarDiffusionFaceCoefficient<
      MultiplyMolecularByDensity, AddTurbulence>(
      density, rho_comp, mu_turb, coefficients, i, j, k, 0, 1, 0,
      coefficients.eddy_h_comp);
    const bool ext_lo =
      ((bc_ptr[bc_comp].lo(1) == ERFBCType::ext_dir) ||
       (bc_ptr[bc_comp].lo(1) == ERFBCType::ext_dir_prim) ||
       (bc_ptr[bc_comp].lo(1) == ERFBCType::ext_dir_upwind &&
        v(i, dom_lo.y, k) >= Real(0.0))) &&
      j == dom_lo.y;
    const bool ext_hi =
      ((bc_ptr[bc_comp].hi(1) == ERFBCType::ext_dir) ||
       (bc_ptr[bc_comp].hi(1) == ERFBCType::ext_dir_prim) ||
       (bc_ptr[bc_comp].hi(1) == ERFBCType::ext_dir_upwind &&
        v(i, dom_hi.y + 1, k) <= Real(0.0))) &&
      j == dom_hi.y + 1;
    const bool on_surface =
      (surface.ylo && j == dom_lo.y) ||
      (surface.yhi && j == dom_hi.y + 1);

    if (ext_lo) {
      yflux(i, j, k, flux_comp) = ScalarDiffusionFlux_N<1>(
        scalar, scalar_comp, i, j, k, coeff, dy_inv,
        mf_vy(i, j, 0) / mf_vx(i, j, 0), true, false);
    } else if (ext_hi) {
      yflux(i, j, k, flux_comp) = ScalarDiffusionFlux_N<1>(
        scalar, scalar_comp, i, j, k, coeff, dy_inv,
        mf_vy(i, j, 0) / mf_vx(i, j, 0), false, true);
    } else if (on_surface && prescribed.use_y_on_side) {
      yflux(i, j, k, flux_comp) = prescribed.y(i, j, k, prescribed.y_comp);
    } else {
      yflux(i, j, k, flux_comp) = ScalarDiffusionFlux_N<1>(
        scalar, scalar_comp, i, j, k, coeff, dy_inv,
        mf_vy(i, j, 0) / mf_vx(i, j, 0), false, false);
    }
  });

  ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real coeff = ScalarDiffusionFaceCoefficient<
      MultiplyMolecularByDensity, AddTurbulence>(
      density, rho_comp, mu_turb, coefficients, i, j, k, 0, 0, 1,
      coefficients.eddy_v_comp);
    const bool ext_lo =
      ((bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir) ||
       (bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir_prim)) &&
      k == dom_lo.z;
    const bool ext_hi =
      ((bc_ptr[bc_comp].hi(2) == ERFBCType::ext_dir) ||
       (bc_ptr[bc_comp].hi(2) == ERFBCType::ext_dir_prim)) &&
      k == dom_hi.z + 1;
    const bool on_surface =
      (surface.zlo && k == dom_lo.z) ||
      (surface.zhi && k == dom_hi.z + 1);

    if (ext_lo) {
      zflux(i, j, k, flux_comp) = ScalarDiffusionFlux_N<2>(
        scalar, scalar_comp, i, j, k, coeff, dz_inv, Real(1.0), true, false);
    } else if (ext_hi) {
      zflux(i, j, k, flux_comp) = ScalarDiffusionFlux_N<2>(
        scalar, scalar_comp, i, j, k, coeff, dz_inv, Real(1.0), false, true);
    } else if (on_surface && prescribed.use_z_on_side) {
      zflux(i, j, k, flux_comp) = prescribed.z(i, j, k, prescribed.z_comp);
    } else if (on_surface) {
      zflux(i, j, k, flux_comp) = Real(0.0);
    } else {
      zflux(i, j, k, flux_comp) = ScalarDiffusionFlux_N<2>(
        scalar, scalar_comp, i, j, k, coeff, dz_inv, Real(1.0), false, false);
    }

    if (diagnostic.enabled && (diagnostic.write_on_surface || !on_surface)) {
      diagnostic.z(i, j, k, diagnostic.z_comp) = zflux(i, j, k, flux_comp);
    }
  });
}

/** Build all native N-grid face fluxes for any intensive scalar field. */
inline void
BuildScalarDiffusionFluxes_N(
  const amrex::Box& bx,
  const amrex::Box& domain,
  const ScalarDiffusionFieldViews& field,
  const ScalarDiffusionFluxPolicy& policy,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
  const amrex::Array4<const amrex::Real>& u,
  const amrex::Array4<const amrex::Real>& v,
  const amrex::Array4<const amrex::Real>& mf_ux,
  const amrex::Array4<const amrex::Real>& mf_uy,
  const amrex::Array4<const amrex::Real>& mf_vx,
  const amrex::Array4<const amrex::Real>& mf_vy,
  const amrex::BCRec* bc_ptr,
  const int bc_comp)
{
  const auto mode = policy.coefficient_mode;
  if (mode.molecular_is_kinematic) {
    if (mode.add_turbulence) {
      BuildScalarDiffusionFluxes_N_Impl<true, true>(
        bx, domain, field, policy, cellSizeInv, u, v, mf_ux, mf_uy, mf_vx,
        mf_vy, bc_ptr, bc_comp);
    } else {
      BuildScalarDiffusionFluxes_N_Impl<true, false>(
        bx, domain, field, policy, cellSizeInv, u, v, mf_ux, mf_uy, mf_vx,
        mf_vy, bc_ptr, bc_comp);
    }
  } else if (mode.add_turbulence) {
    BuildScalarDiffusionFluxes_N_Impl<false, true>(
      bx, domain, field, policy, cellSizeInv, u, v, mf_ux, mf_uy, mf_vx,
      mf_vy, bc_ptr, bc_comp);
  } else {
    BuildScalarDiffusionFluxes_N_Impl<false, false>(
      bx, domain, field, policy, cellSizeInv, u, v, mf_ux, mf_uy, mf_vx,
      mf_vy, bc_ptr, bc_comp);
  }
}

template <bool MultiplyMolecularByDensity, bool AddTurbulence>
inline void
BuildScalarDiffusionFluxes_S_Impl(
  const amrex::Box& bx,
  const amrex::Box& domain,
  const ScalarDiffusionFieldViews& field,
  const ScalarDiffusionFluxPolicy& policy,
  const amrex::Real dx_inv,
  const amrex::Real dy_inv,
  const amrex::Real* dz,
  const int klo,
  const int khi,
  const amrex::Array4<const amrex::Real>& mf_ux,
  const amrex::Array4<const amrex::Real>& mf_uy,
  const amrex::Array4<const amrex::Real>& mf_vx,
  const amrex::Array4<const amrex::Real>& mf_vy,
  const amrex::BCRec* bc_ptr,
  const int bc_comp)
{
  using namespace amrex;
  const auto dom_lo = lbound(domain);
  const auto dom_hi = ubound(domain);
  const Box xbx = surroundingNodes(bx, 0);
  const Box ybx = surroundingNodes(bx, 1);
  const Box zbx = surroundingNodes(bx, 2);
  const auto coefficients = policy.coefficients;
  const auto prescribed = policy.prescribed;
  const auto diagnostic = policy.diagnostic;
  const auto surface = policy.surface;
  const int scalar_comp = field.scalar_comp;
  const int rho_comp = field.rho_comp;
  const int flux_comp = field.flux_comp;
  const auto scalar = field.scalar;
  const auto density = field.density;
  const auto mu_turb = field.mu_turb;
  const auto xflux = field.xflux;
  const auto yflux = field.yflux;
  const auto zflux = field.zflux;

  ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real coeff = ScalarDiffusionFaceCoefficient<
      MultiplyMolecularByDensity, AddTurbulence>(
      density, rho_comp, mu_turb, coefficients, i, j, k, 1, 0, 0,
      coefficients.eddy_h_comp);
    const bool on_surface =
      (surface.xlo && i == dom_lo.x) ||
      (surface.xhi && i == dom_hi.x + 1);
    if (on_surface && prescribed.use_x_on_side) {
      xflux(i, j, k, flux_comp) = prescribed.x(i, j, k, prescribed.x_comp);
    } else {
      xflux(i, j, k, flux_comp) = ScalarDiffusionFlux_S_Horizontal<0>(
        scalar, scalar_comp, i, j, k, coeff, dx_inv, mf_ux(i, j, 0));
    }
  });

  ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real coeff = ScalarDiffusionFaceCoefficient<
      MultiplyMolecularByDensity, AddTurbulence>(
      density, rho_comp, mu_turb, coefficients, i, j, k, 0, 1, 0,
      coefficients.eddy_h_comp);
    const bool on_surface =
      (surface.ylo && j == dom_lo.y) ||
      (surface.yhi && j == dom_hi.y + 1);
    if (on_surface && prescribed.use_y_on_side) {
      yflux(i, j, k, flux_comp) = prescribed.y(i, j, k, prescribed.y_comp);
    } else {
      yflux(i, j, k, flux_comp) = ScalarDiffusionFlux_S_Horizontal<1>(
        scalar, scalar_comp, i, j, k, coeff, dy_inv, mf_vy(i, j, 0));
    }
  });

  ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real coeff = ScalarDiffusionFaceCoefficient<
      MultiplyMolecularByDensity, AddTurbulence>(
      density, rho_comp, mu_turb, coefficients, i, j, k, 0, 0, 1,
      coefficients.eddy_v_comp);
    const bool ext_lo =
      ((bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir) ||
       (bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir_prim)) &&
      k == dom_lo.z;
    const bool ext_hi =
      ((bc_ptr[bc_comp].hi(2) == ERFBCType::ext_dir) ||
       (bc_ptr[bc_comp].hi(2) == ERFBCType::ext_dir_prim)) &&
      k == dom_hi.z + 1;
    const bool on_surface =
      (surface.zlo && k == dom_lo.z) ||
      (surface.zhi && k == dom_hi.z + 1);
    const Real grad = StretchedScalarGradient(
      scalar, scalar_comp, i, j, k, dz, klo, khi, ext_lo, ext_hi);
    if (on_surface && prescribed.use_z_on_side) {
      zflux(i, j, k, flux_comp) = prescribed.z(i, j, k, prescribed.z_comp);
    } else if (on_surface) {
      zflux(i, j, k, flux_comp) = Real(0.0);
    } else {
      zflux(i, j, k, flux_comp) = -coeff * grad;
    }
    if (diagnostic.enabled && (diagnostic.write_on_surface || !on_surface)) {
      diagnostic.z(i, j, k, diagnostic.z_comp) = zflux(i, j, k, flux_comp);
    }
  });

  // S-grid raw horizontal fluxes are adjusted in the same separate kernels as
  // the native implementation before scaling or applying divergence.
  ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    xflux(i, j, k, flux_comp) /= mf_uy(i, j, 0);
  });
  ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    yflux(i, j, k, flux_comp) /= mf_vx(i, j, 0);
  });
}

/** Build and map-adjust all native S-grid face fluxes for an intensive field. */
inline void
BuildScalarDiffusionFluxes_S(
  const amrex::Box& bx,
  const amrex::Box& domain,
  const ScalarDiffusionFieldViews& field,
  const ScalarDiffusionFluxPolicy& policy,
  const amrex::Real dx_inv,
  const amrex::Real dy_inv,
  const amrex::Real* dz,
  const int klo,
  const int khi,
  const amrex::Array4<const amrex::Real>& mf_ux,
  const amrex::Array4<const amrex::Real>& mf_uy,
  const amrex::Array4<const amrex::Real>& mf_vx,
  const amrex::Array4<const amrex::Real>& mf_vy,
  const amrex::BCRec* bc_ptr,
  const int bc_comp)
{
  const auto mode = policy.coefficient_mode;
  if (mode.molecular_is_kinematic) {
    if (mode.add_turbulence) {
      BuildScalarDiffusionFluxes_S_Impl<true, true>(
        bx, domain, field, policy, dx_inv, dy_inv, dz, klo, khi, mf_ux,
        mf_uy, mf_vx, mf_vy, bc_ptr, bc_comp);
    } else {
      BuildScalarDiffusionFluxes_S_Impl<true, false>(
        bx, domain, field, policy, dx_inv, dy_inv, dz, klo, khi, mf_ux,
        mf_uy, mf_vx, mf_vy, bc_ptr, bc_comp);
    }
  } else if (mode.add_turbulence) {
    BuildScalarDiffusionFluxes_S_Impl<false, true>(
      bx, domain, field, policy, dx_inv, dy_inv, dz, klo, khi, mf_ux,
      mf_uy, mf_vx, mf_vy, bc_ptr, bc_comp);
  } else {
    BuildScalarDiffusionFluxes_S_Impl<false, false>(
      bx, domain, field, policy, dx_inv, dy_inv, dz, klo, khi, mf_ux,
      mf_uy, mf_vx, mf_vy, bc_ptr, bc_comp);
  }
}

template <bool MultiplyMolecularByDensity, bool AddTurbulence>
inline void
BuildScalarDiffusionFluxes_T_Impl(
  const amrex::Box& bx,
  const amrex::Box& domain,
  const ScalarDiffusionFieldViews& field,
  const ScalarDiffusionFluxPolicy& policy,
  const bool rotate,
  const amrex::Array4<const amrex::Real>& z_nd,
  const amrex::Array4<const amrex::Real>& z_cc,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
  const amrex::Array4<const amrex::Real>& mf_ux,
  const amrex::Array4<const amrex::Real>& mf_uy,
  const amrex::Array4<const amrex::Real>& mf_vx,
  const amrex::Array4<const amrex::Real>& mf_vy,
  const amrex::BCRec* bc_ptr,
  const int bc_comp)
{
  using namespace amrex;
  const auto dom_lo = lbound(domain);
  const auto dom_hi = ubound(domain);
  Box xbx = surroundingNodes(bx, 0);
  Box ybx = surroundingNodes(bx, 1);
  const Box zbx = surroundingNodes(bx, 2);
  if (xbx.smallEnd(2) != dom_lo.z) xbx.growLo(2, 1);
  if (ybx.smallEnd(2) != dom_lo.z) ybx.growLo(2, 1);
  if (xbx.bigEnd(2) != dom_hi.z) xbx.growHi(2, 1);
  if (ybx.bigEnd(2) != dom_hi.z) ybx.growHi(2, 1);

  const Real dx_inv = cellSizeInv[0];
  const Real dy_inv = cellSizeInv[1];
  const Real dz_inv = cellSizeInv[2];
  const auto coefficients = policy.coefficients;
  const auto prescribed = policy.prescribed;
  const auto diagnostic = policy.diagnostic;
  const auto surface = policy.surface;
  const int scalar_comp = field.scalar_comp;
  const int rho_comp = field.rho_comp;
  const int flux_comp = field.flux_comp;
  const auto scalar = field.scalar;
  const auto density = field.density;
  const auto mu_turb = field.mu_turb;
  const auto xflux = field.xflux;
  const auto yflux = field.yflux;
  const auto zflux = field.zflux;

  ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real coeff = ScalarDiffusionFaceCoefficient<
      MultiplyMolecularByDensity, AddTurbulence>(
      density, rho_comp, mu_turb, coefficients, i, j, k, 1, 0, 0,
      coefficients.eddy_h_comp);
    const Real h_xi = Compute_h_xi_AtIface(i, j, k, cellSizeInv, z_nd);
    const bool on_side =
      (surface.xlo && i == dom_lo.x) ||
      (surface.xhi && i == dom_hi.x + 1);
    const bool on_bottom_z =
      surface.zlo && rotate && k == dom_lo.z;

    const Real idz_hi = Real(1.0) / (z_cc(i, j, k + 1) - z_cc(i, j, k - 1));
    const Real idz_lo =
      Real(1.0) / (z_cc(i - 1, j, k + 1) - z_cc(i - 1, j, k - 1));
    const Real grad_z = Real(0.5) *
      (scalar(i, j, k + 1, scalar_comp) * idz_hi +
       scalar(i - 1, j, k + 1, scalar_comp) * idz_lo -
       scalar(i, j, k - 1, scalar_comp) * idz_hi -
       scalar(i - 1, j, k - 1, scalar_comp) * idz_lo);
    const Real grad_xi = dx_inv *
      (scalar(i, j, k, scalar_comp) - scalar(i - 1, j, k, scalar_comp));

    if (on_side && prescribed.use_x_on_side) {
      xflux(i, j, k, flux_comp) = prescribed.x(i, j, k, prescribed.x_comp);
    } else if (on_bottom_z && prescribed.use_x_on_bottom_z) {
      xflux(i, j, k, flux_comp) = prescribed.x(i, j, 0, prescribed.x_comp);
    } else {
      xflux(i, j, k, flux_comp) = ScalarDiffusionFlux_Tx(
        coeff, mf_ux(i, j, 0), grad_xi, h_xi, grad_z);
    }
  });

  ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real coeff = ScalarDiffusionFaceCoefficient<
      MultiplyMolecularByDensity, AddTurbulence>(
      density, rho_comp, mu_turb, coefficients, i, j, k, 0, 1, 0,
      coefficients.eddy_h_comp);
    const Real h_eta = Compute_h_eta_AtJface(i, j, k, cellSizeInv, z_nd);
    const bool on_side =
      (surface.ylo && j == dom_lo.y) ||
      (surface.yhi && j == dom_hi.y + 1);
    const bool on_bottom_z =
      surface.zlo && rotate && k == dom_lo.z;

    const Real idz_hi = Real(1.0) / (z_cc(i, j, k + 1) - z_cc(i, j, k - 1));
    const Real idz_lo =
      Real(1.0) / (z_cc(i, j - 1, k + 1) - z_cc(i, j - 1, k - 1));
    const Real grad_z = Real(0.5) *
      (scalar(i, j, k + 1, scalar_comp) * idz_hi +
       scalar(i, j - 1, k + 1, scalar_comp) * idz_lo -
       scalar(i, j, k - 1, scalar_comp) * idz_hi -
       scalar(i, j - 1, k - 1, scalar_comp) * idz_lo);
    const Real grad_eta = dy_inv *
      (scalar(i, j, k, scalar_comp) - scalar(i, j - 1, k, scalar_comp));

    if (on_side && prescribed.use_y_on_side) {
      yflux(i, j, k, flux_comp) = prescribed.y(i, j, k, prescribed.y_comp);
    } else if (on_bottom_z && prescribed.use_y_on_bottom_z) {
      yflux(i, j, k, flux_comp) = prescribed.y(i, j, 0, prescribed.y_comp);
    } else {
      yflux(i, j, k, flux_comp) = ScalarDiffusionFlux_Ty(
        coeff, mf_vy(i, j, 0), grad_eta, h_eta, grad_z);
    }
  });

  ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real coeff = ScalarDiffusionFaceCoefficient<
      MultiplyMolecularByDensity, AddTurbulence>(
      density, rho_comp, mu_turb, coefficients, i, j, k, 0, 0, 1,
      coefficients.eddy_v_comp);
    const bool ext_lo =
      ((bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir) ||
       (bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir_prim)) &&
      k == dom_lo.z;
    const bool ext_hi =
      ((bc_ptr[bc_comp].hi(2) == ERFBCType::ext_dir) ||
       (bc_ptr[bc_comp].hi(2) == ERFBCType::ext_dir_prim)) &&
      k == dom_hi.z + 1;
    const bool on_surface =
      (surface.zlo && k == dom_lo.z) ||
      (surface.zhi && k == dom_hi.z + 1);

    Real grad_z;
    if (ext_lo) {
      const Real zm = Compute_Z_AtWFace(i, j, k + 1, z_nd);
      const Real dz0 = zm - Compute_Z_AtWFace(i, j, k, z_nd);
      const Real dz1 = Compute_Z_AtWFace(i, j, k + 2, z_nd) - zm;
      const Real idz0 = Real(1.0) / dz0;
      const Real f = dz1 / dz0 + Real(2.0);
      const Real f2 = f * f;
      const Real c3 = Real(2.0) / (f - f2);
      const Real c2 = -f2 * c3;
      const Real c1 = -(Real(1.0) - f2) * c3;
      grad_z = idz0 * (c1 * scalar(i, j, k - 1, scalar_comp) +
                       c2 * scalar(i, j, k, scalar_comp) +
                       c3 * scalar(i, j, k + 1, scalar_comp));
    } else if (ext_hi) {
      const Real zm = Compute_Z_AtWFace(i, j, k - 1, z_nd);
      const Real dz0 = Compute_Z_AtWFace(i, j, k, z_nd) - zm;
      const Real dz1 = zm - Compute_Z_AtWFace(i, j, k - 2, z_nd);
      const Real idz0 = Real(1.0) / dz0;
      const Real f = dz1 / dz0 + Real(2.0);
      const Real f2 = f * f;
      const Real c3 = Real(2.0) / (f - f2);
      const Real c2 = -f2 * c3;
      const Real c1 = -(Real(1.0) - f2) * c3;
      grad_z = idz0 * (-(c1 * scalar(i, j, k, scalar_comp) +
                          c2 * scalar(i, j, k - 1, scalar_comp) +
                          c3 * scalar(i, j, k - 2, scalar_comp)));
    } else {
      const Real h_zeta = Compute_h_zeta_AtKface(i, j, k, cellSizeInv, z_nd);
      grad_z = (dz_inv / h_zeta) *
        (scalar(i, j, k, scalar_comp) - scalar(i, j, k - 1, scalar_comp));
    }

    if (on_surface && prescribed.use_z_on_side) {
      zflux(i, j, k, flux_comp) = prescribed.z(i, j, k, prescribed.z_comp);
    } else if (on_surface) {
      zflux(i, j, k, flux_comp) = Real(0.0);
    } else {
      zflux(i, j, k, flux_comp) = -coeff * grad_z;
    }
    if (diagnostic.enabled && (diagnostic.write_on_surface || !on_surface)) {
      diagnostic.z(i, j, k, diagnostic.z_comp) = zflux(i, j, k, flux_comp);
    }
  });
}

/** Build all native terrain-grid face fluxes, including the grown horizontal
 * face regions needed by the mapped z-face transfer. */
inline void
BuildScalarDiffusionFluxes_T(
  const amrex::Box& bx,
  const amrex::Box& domain,
  const ScalarDiffusionFieldViews& field,
  const ScalarDiffusionFluxPolicy& policy,
  const bool rotate,
  const amrex::Array4<const amrex::Real>& z_nd,
  const amrex::Array4<const amrex::Real>& z_cc,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
  const amrex::Array4<const amrex::Real>& mf_ux,
  const amrex::Array4<const amrex::Real>& mf_uy,
  const amrex::Array4<const amrex::Real>& mf_vx,
  const amrex::Array4<const amrex::Real>& mf_vy,
  const amrex::BCRec* bc_ptr,
  const int bc_comp)
{
  const auto mode = policy.coefficient_mode;
  if (mode.molecular_is_kinematic) {
    if (mode.add_turbulence) {
      BuildScalarDiffusionFluxes_T_Impl<true, true>(
        bx, domain, field, policy, rotate, z_nd, z_cc, cellSizeInv, mf_ux,
        mf_uy, mf_vx, mf_vy, bc_ptr, bc_comp);
    } else {
      BuildScalarDiffusionFluxes_T_Impl<true, false>(
        bx, domain, field, policy, rotate, z_nd, z_cc, cellSizeInv, mf_ux,
        mf_uy, mf_vx, mf_vy, bc_ptr, bc_comp);
    }
  } else if (mode.add_turbulence) {
    BuildScalarDiffusionFluxes_T_Impl<false, true>(
      bx, domain, field, policy, rotate, z_nd, z_cc, cellSizeInv, mf_ux,
      mf_uy, mf_vx, mf_vy, bc_ptr, bc_comp);
  } else {
    BuildScalarDiffusionFluxes_T_Impl<false, false>(
      bx, domain, field, policy, rotate, z_nd, z_cc, cellSizeInv, mf_ux,
      mf_uy, mf_vx, mf_vy, bc_ptr, bc_comp);
  }
}

/** Apply the complete mapped terrain divergence, including x/y-to-z-face
 * transfer and one-sided physical bottom/top extrapolation. */
inline void
ApplyScalarDiffusionMappedDivergence_T(
  const amrex::Box& bx,
  const amrex::Box& domain,
  const ScalarDiffusionFieldViews& field,
  const amrex::Array4<const amrex::Real>& z_nd,
  const amrex::Array4<const amrex::Real>& ax,
  const amrex::Array4<const amrex::Real>& ay,
  const amrex::Array4<const amrex::Real>& detJ,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
  const amrex::Array4<const amrex::Real>& mf_mx,
  const amrex::Array4<const amrex::Real>& mf_uy,
  const amrex::Array4<const amrex::Real>& mf_my,
  const amrex::Array4<const amrex::Real>& mf_vx,
  const bool suppress_lower_surface_z_flux)
{
  using namespace amrex;
  const auto dom_lo = lbound(domain);
  const auto dom_hi = ubound(domain);
  const Real dx_inv = cellSizeInv[0];
  const Real dy_inv = cellSizeInv[1];
  const Real dz_inv = cellSizeInv[2];
  const Real half = Real(0.5);
  const Real fourth = Real(0.25);
  const auto xflux = field.xflux;
  const auto yflux = field.yflux;
  const auto zflux = field.zflux;
  const auto rhs = field.rhs;
  const int flux_comp = field.flux_comp;
  const int rhs_comp = field.rhs_comp;

  ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real h_xi_lo = Compute_h_xi_AtKface(i, j, k, cellSizeInv, z_nd);
    const Real h_eta_lo = Compute_h_eta_AtKface(i, j, k, cellSizeInv, z_nd);
    Real xbar_lo;
    Real ybar_lo;
    if (k == dom_lo.z) {
      const Real xlo = half * (xflux(i, j, k, flux_comp) +
                                xflux(i + 1, j, k, flux_comp));
      const Real xhi = half * (xflux(i, j, k + 1, flux_comp) +
                                xflux(i + 1, j, k + 1, flux_comp));
      xbar_lo = Real(1.5) * xlo - half * xhi;
      const Real ylo = half * (yflux(i, j, k, flux_comp) +
                                yflux(i, j + 1, k, flux_comp));
      const Real yhi = half * (yflux(i, j, k + 1, flux_comp) +
                                yflux(i, j + 1, k + 1, flux_comp));
      ybar_lo = Real(1.5) * ylo - half * yhi;
    } else {
      xbar_lo = fourth *
        (xflux(i, j, k, flux_comp) + xflux(i + 1, j, k, flux_comp) +
         xflux(i, j, k - 1, flux_comp) + xflux(i + 1, j, k - 1, flux_comp));
      ybar_lo = fourth *
        (yflux(i, j, k, flux_comp) + yflux(i, j + 1, k, flux_comp) +
         yflux(i, j, k - 1, flux_comp) + yflux(i, j + 1, k - 1, flux_comp));
    }

    const Real h_xi_hi = Compute_h_xi_AtKface(i, j, k + 1, cellSizeInv, z_nd);
    const Real h_eta_hi = Compute_h_eta_AtKface(i, j, k + 1, cellSizeInv, z_nd);
    Real xbar_hi;
    Real ybar_hi;
    if (k == dom_hi.z) {
      const Real xlo = half * (xflux(i, j, k - 1, flux_comp) +
                                xflux(i + 1, j, k - 1, flux_comp));
      const Real xhi = half * (xflux(i, j, k, flux_comp) +
                                xflux(i + 1, j, k, flux_comp));
      xbar_hi = Real(1.5) * xhi - half * xlo;
      const Real ylo = half * (yflux(i, j, k - 1, flux_comp) +
                                yflux(i, j + 1, k - 1, flux_comp));
      const Real yhi = half * (yflux(i, j, k, flux_comp) +
                                yflux(i, j + 1, k, flux_comp));
      ybar_hi = Real(1.5) * yhi - half * ylo;
    } else {
      xbar_hi = fourth *
        (xflux(i, j, k + 1, flux_comp) + xflux(i + 1, j, k + 1, flux_comp) +
         xflux(i, j, k, flux_comp) + xflux(i + 1, j, k, flux_comp));
      ybar_hi = fourth *
        (yflux(i, j, k + 1, flux_comp) + yflux(i, j + 1, k + 1, flux_comp) +
         yflux(i, j, k, flux_comp) + yflux(i, j + 1, k, flux_comp));
    }

    const Real Gz_lo =
      (suppress_lower_surface_z_flux && k == dom_lo.z)
        ? Real(0.0)
        : TerrainDiffusionGz(zflux(i, j, k, flux_comp), mf_mx(i, j, 0),
                             h_xi_lo, xbar_lo, mf_my(i, j, 0), h_eta_lo,
                             ybar_lo);
    const Real Gz_hi = TerrainDiffusionGz(
      zflux(i, j, k + 1, flux_comp), mf_mx(i, j, 0), h_xi_hi, xbar_hi,
      mf_my(i, j, 0), h_eta_hi, ybar_hi);

    const Real contribution = TerrainDiffusionDivergence_T(
      xflux(i + 1, j, k, flux_comp), ax(i + 1, j, k), mf_uy(i + 1, j, 0),
      xflux(i, j, k, flux_comp), ax(i, j, k), mf_uy(i, j, 0),
      yflux(i, j + 1, k, flux_comp), ay(i, j + 1, k), mf_vx(i, j + 1, 0),
      yflux(i, j, k, flux_comp), ay(i, j, k), mf_vx(i, j, 0), Gz_hi, Gz_lo,
      dx_inv, dy_inv, dz_inv, mf_mx(i, j, 0), mf_my(i, j, 0), detJ(i, j, k));
    rhs(i, j, k, rhs_comp) -= contribution;
  });
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
