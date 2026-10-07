/** \file ERF_TerrainDiffusionLimits.H */

#ifndef ERF_TERRAIN_DIFFUSION_LIMITS_H_
#define ERF_TERRAIN_DIFFUSION_LIMITS_H_

#include <AMReX_Array4.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

#include "ERF_NumericalConstants.H"

#include <cmath>

/**
 * Terrain drop across one cell of a terrain-fitted mesh.
 *
 * dx_drop and dy_drop are the height differences across the cell along the
 * coordinate surfaces, each the mean of the absolute differences on the four
 * cell edges in that direction; dz is the physical thickness of the cell, the
 * mean of its four vertical edges (detJ * dz at the cell center).  All three
 * are in metres and read only the eight nodes of the cell.
 */
struct TerrainCellDrops
{
    amrex::Real dx_drop;
    amrex::Real dy_drop;
    amrex::Real dz;
};

/**
 * Compute the terrain drops across cell (i,j,k) from the nodal heights.
 *
 * @param[in] i,j,k cell-centered indices
 * @param[in] z_nd  nodal physical heights
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
TerrainCellDrops
ComputeTerrainCellDrops (int i, int j, int k,
                         const amrex::Array4<const amrex::Real>& z_nd) noexcept
{
    using amrex::Math::abs;
    TerrainCellDrops d;
    d.dx_drop = fourth * ( abs(z_nd(i+1,j  ,k  ) - z_nd(i,j  ,k  ))
                         + abs(z_nd(i+1,j+1,k  ) - z_nd(i,j+1,k  ))
                         + abs(z_nd(i+1,j  ,k+1) - z_nd(i,j  ,k+1))
                         + abs(z_nd(i+1,j+1,k+1) - z_nd(i,j+1,k+1)) );
    d.dy_drop = fourth * ( abs(z_nd(i  ,j+1,k  ) - z_nd(i  ,j,k  ))
                         + abs(z_nd(i+1,j+1,k  ) - z_nd(i+1,j,k  ))
                         + abs(z_nd(i  ,j+1,k+1) - z_nd(i  ,j,k+1))
                         + abs(z_nd(i+1,j+1,k+1) - z_nd(i+1,j,k+1)) );
    d.dz      = fourth * ( (z_nd(i  ,j  ,k+1) - z_nd(i  ,j  ,k))
                         + (z_nd(i+1,j  ,k+1) - z_nd(i+1,j  ,k))
                         + (z_nd(i  ,j+1,k+1) - z_nd(i  ,j+1,k))
                         + (z_nd(i+1,j+1,k+1) - z_nd(i+1,j+1,k)) );
    return d;
}

/**
 * Slope factor alpha = max( sqrt(dx_drop^2 + dy_drop^2) / dz, 1 ) of a cell:
 * the terrain drop across one cell measured in cell thicknesses, which is
 * h * Delta_x / Delta_z with h the physical coordinate slope.  This is the
 * factor of WRF's Smagorinsky-2D slope limiter (dyn_em/module_diffusion_em.F,
 * smag2d_km, diff_opt = 2).  WRF averages |z_x| at the two x-faces and the two
 * w-levels of the cell and multiplies by dx / msftx; here the drop is taken on
 * the four x-edges of the cell (ERF heights are nodal) and carries no map
 * factor, so that alpha is exactly h * Delta_x / Delta_z with the physical
 * spacing.  The two agree when the map factor is 1.
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
TerrainSlopeFactor (const TerrainCellDrops& d) noexcept
{
    const amrex::Real a = std::sqrt(d.dx_drop*d.dx_drop + d.dy_drop*d.dy_drop) / d.dz;
    return amrex::max(a, amrex::Real(1.0));
}

/** Deformation threshold speed of WRF's slope limiter, def_limit = max(10/Delta_h, 1e-3) (m/s). */
constexpr amrex::Real smag2d_def_limit_speed = amrex::Real(10.0);
/** Floor of WRF's deformation threshold (1/s). */
constexpr amrex::Real smag2d_def_limit_min   = amrex::Real(1.0e-3);

/**
 * Apply the opt-in limits of WRF's smag2d_km to a Smagorinsky-2D horizontal
 * kinematic eddy viscosity, in WRF's order:
 *   1. the cap  nu <= kh_cap * Delta_h          (only if kh_cap > 0; WRF uses 10 m/s);
 *   2. the slope limiter (only if use_slope):
 *        nu /= alpha^2  if |D_h| > max(10/Delta_h, 1e-3),
 *        nu /= alpha    otherwise.
 * |D_h| is the horizontal deformation sqrt((u_x - v_y)^2 + (u_y + v_x)^2),
 * which equals ERF's sqrt(2 S_mn S_mn) of ComputeSmnSmn2D and WRF's
 * sqrt(def2).  Delta_h = sqrt(dx/m_x * dy/m_y), WRF's mlen_h.
 *
 * @param[in] nu        Cs^2 Delta_h^2 |D_h| (m^2/s)
 * @param[in] strain    |D_h| (1/s)
 * @param[in] DeltaH    horizontal filter width (m)
 * @param[in] alpha     slope factor from TerrainSlopeFactor (>= 1)
 * @param[in] kh_cap    cap speed (m/s); <= 0 leaves nu uncapped
 * @param[in] use_slope apply the slope limiter
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Smag2DLimitedViscosity (amrex::Real nu, amrex::Real strain, amrex::Real DeltaH,
                        amrex::Real alpha, amrex::Real kh_cap, bool use_slope) noexcept
{
    if (kh_cap > zero) {
        const amrex::Real nu_max = kh_cap * DeltaH;
        nu = (nu < nu_max) ? nu : nu_max;
    }
    if (use_slope) {
        const amrex::Real def_limit = amrex::max(amrex::Real(smag2d_def_limit_speed) / DeltaH,
                                                 amrex::Real(smag2d_def_limit_min));
        nu = (strain > def_limit) ? nu / (alpha*alpha) : nu / alpha;
    }
    return nu;
}

/**
 * Explicit diffusion rates of one cell, for the diffusive time-step estimate
 * dt * rate <= Fourier limit (1/2 for forward Euler in 3-D; the three-stage
 * Runge-Kutta scheme's real-axis bound is about 2.5 / 4 = 0.63).
 *
 * Momentum: mu_h * max(2 a + b, a + 2 b) for the horizontal stresses
 * (tau11 = 2 mu S11), 2 mu_p h^2 / dz^2 for the terrain-metric part, and
 * max(e_uv, 2 e_w) mu_v / dz^2 for the vertical part (tau33 = 2 mu S33).  The
 * metric term is an upper bound: the projected stresses of the terrain stress
 * (ComputeStressVarVisc_T) carry mu_v today (ERF issue #4214) and would carry
 * mu_h once that is changed, hence mu_p = max(mu_h, mu_v); with mu_v there, the
 * K_h part of that coupling enters only through the mixed xi-zeta terms, so
 * the bound can overstate that coupling by about alpha / 2 on steep slopes.
 *
 * Scalar: K_h (a + b) + K_h h^2 / dz^2 + e K_v / dz^2; the K_h h^2 part
 * (TerrainDiffusionGz) is explicit.
 *
 * a = (m_x / dx)^2, b = (m_y / dy)^2 and h^2 = (m_x dx_drop / dx)^2 +
 * (m_y dy_drop / dy)^2 is the squared physical slope; on a mesh without
 * terrain pass h2 = 0.  e is the explicit fraction of the vertical diffusion
 * (ExplicitVerticalFraction): 0 when the implicit solve makes it A-stable,
 * 1 otherwise.  All diffusivities are dynamic (rho * nu); rates are in 1/s.
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
MomentumDiffusiveRate (amrex::Real rho, amrex::Real mu_h, amrex::Real mu_v,
                       amrex::Real a, amrex::Real b, amrex::Real h2, amrex::Real dzinv2,
                       amrex::Real e_uv, amrex::Real e_w) noexcept
{
    const amrex::Real horiz = (a > b) ? (two*a + b) : (a + two*b);
    const amrex::Real mu_p  = (mu_h > mu_v) ? mu_h : mu_v;
    const amrex::Real vert  = (e_uv > two*e_w) ? e_uv : two*e_w;
    return ( mu_h * horiz + two * mu_p * h2 * dzinv2 + vert * mu_v * dzinv2 ) / rho;
}

/** Scalar counterpart of MomentumDiffusiveRate (see there). */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ScalarDiffusiveRate (amrex::Real rho, amrex::Real K_h, amrex::Real K_v,
                     amrex::Real a, amrex::Real b, amrex::Real h2, amrex::Real dzinv2,
                     amrex::Real e_vert) noexcept
{
    return ( K_h * (a + b) + K_h * h2 * dzinv2 + e_vert * K_v * dzinv2 ) / rho;
}

/**
 * Explicit fraction of the vertical diffusion for the diffusive time-step check: 0 if the
 * partly implicit Runge-Kutta scheme is A-stable for the vertical diffusion, 1 otherwise.
 *
 * Each stage k starts from S^n with stage step c_k dt; its explicit part (1 - f_k) acts on
 * the previous stage's state and its implicit part f_k on its own, so for dS/dt = lambda S,
 * z = lambda dt,  S_k = (1 + (1 - f_k) c_k z S_{k-1}) / (1 - f_k c_k z),  S_0 = 1, and the
 * amplification factor is the last S_k.  The vertical diffusion is treated as implicit
 * (e = 0) when |S_k| <= 1 at every stage for z in [-1e8, -1e-4]: each stage's state feeds
 * the next stage's explicit terms, so an intermediate stage that amplifies a mode is not
 * acceptable even if the last one damps it.  With ERF's three-stage scheme (c = 1/3, 1/2, 1)
 * the default erf.vert_implicit_fac = 1 1 0 gives G = (1 + z/2) / (1 - z/2) (Crank-Nicolson)
 * and 1 1 1 backward Euler, both with bounded stages; the anelastic MidPoint scheme is
 * c = 1/2, 1 with factors 1 0, also Crank-Nicolson.
 *
 * @param[in] stage_frac c_k of each stage
 * @param[in] fac        implicit factor f_k of each stage (erf.vert_implicit_fac)
 * @param[in] implicit   whether the implicit solve treats this component at all
 */
template <typename V>
amrex::Real
ExplicitVerticalFraction (const V& stage_frac, const V& fac, bool implicit)
{
    const int nstages = static_cast<int>(stage_frac.size());
    if (!implicit || static_cast<int>(fac.size()) < nstages) { return one; }
    for (int n = 0; n <= 120; ++n) {
        const double z = -std::pow(10.0, -4.0 + 0.1*n);
        double S = 1.0;
        for (int k = 0; k < nstages; ++k) {
            const double f = amrex::min(amrex::max(static_cast<double>(fac[k]), 0.0), 1.0);
            const double c = static_cast<double>(stage_frac[k]);
            S = (1.0 + (1.0 - f) * c * z * S) / (1.0 - f * c * z);
            // every stage, not only the last: each stage's state feeds the explicit terms
            // (advection, horizontal and metric diffusion, the strain) of the next one
            if (std::abs(S) > 1.0 + 1.0e-12) { return one; }
        }
    }
    return zero;
}

#endif
