/** \file ERF_TerrainDiffusionLimits.H */

#ifndef ERF_TERRAIN_DIFFUSION_LIMITS_H_
#define ERF_TERRAIN_DIFFUSION_LIMITS_H_

#include <AMReX_Array4.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Array.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

#include "ERF_NumericalConstants.H"
#include "ERF_IndexDefines.H"

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
 * Momentum: the largest of the three component rows,
 *   u: mu_h (2 a + b) + e_uv mu_v / dz^2 + mu_v sqrt(a c)
 *   v: mu_h (a + 2 b) + e_uv mu_v / dz^2 + mu_v sqrt(b c)
 *   w: mu_v (a + b)   + 2 e_w mu_v / dz^2 + mu_v (sqrt(a c) + sqrt(b c)),
 * with c = 1/dz^2, plus 2 mu_p h^2 / dz^2 for the terrain-metric part.  The horizontal
 * diffusion of w uses mu_v: tau31 and tau32 carry the same edge-averaged mu_v as tau13 and
 * tau23.  The sqrt terms are the u-w and v-w cross derivatives (dw/dx in tau13, du/dz in
 * tau31, ...), which stay explicit whether or not the implicit solve takes the second
 * vertical derivative; they are always counted, since with dz > dx/m they exceed the
 * vertical term (dx = 100 m, dz = 500 m: five times).  (The expansion
 * rate in tau11 also couples u to dw/dz, with mu_h; that is covered by "an estimate when
 * mu_h != mu_v" below.)  The metric term is an upper bound: the
 * projected stresses of the terrain stress (ComputeStressVarVisc_T) carry mu_v today (ERF
 * issue #4214) and would carry mu_h once that is changed, hence mu_p = max(mu_h, mu_v);
 * with mu_v there, the K_h part of that coupling enters only through the mixed xi-zeta
 * terms, so the bound can overstate that coupling by about alpha / 2 on steep slopes.
 *
 * Scalar: K_h (a + b) + K_h h^2 / dz^2 + e K_v / dz^2; the K_h h^2 part
 * (TerrainDiffusionGz) is explicit.
 *
 * a = (m_x / dx)^2, b = (m_y / dy)^2 and h^2 = (m_x dx_drop / dx)^2 +
 * (m_y dy_drop / dy)^2 is the squared physical slope; on a mesh without
 * terrain pass h2 = 0.  e is the explicit fraction of the vertical diffusion
 * (ExplicitVerticalFraction): 0 when the implicit solve keeps every stage bounded,
 * 1 otherwise.  All diffusivities are dynamic (rho * nu); rates are in 1/s.
 * The caller passes the largest diffusivities and the smallest density around the
 * cell (NeighbourhoodMax, NeighbourhoodMinPositive), which bounds the row whatever
 * face or edge averaging the operator uses.
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
MomentumDiffusiveRate (amrex::Real rho, amrex::Real mu_h, amrex::Real mu_v,
                       amrex::Real a, amrex::Real b, amrex::Real h2, amrex::Real dzinv2,
                       amrex::Real e_uv, amrex::Real e_w) noexcept
{
    const amrex::Real xa   = std::sqrt(a * dzinv2);
    const amrex::Real xb   = std::sqrt(b * dzinv2);
    const amrex::Real r_u  = mu_h * (two*a + b) + e_uv * mu_v * dzinv2 + mu_v * xa;
    const amrex::Real r_v  = mu_h * (a + two*b) + e_uv * mu_v * dzinv2 + mu_v * xb;
    const amrex::Real r_w  = mu_v * (a + b)     + two * e_w * mu_v * dzinv2 + mu_v * (xa + xb);
    const amrex::Real r_uv = (r_u > r_v) ? r_u : r_v;
    const amrex::Real row  = (r_uv > r_w) ? r_uv : r_w;
    const amrex::Real mu_p = (mu_h > mu_v) ? mu_h : mu_v;
    return ( row + two * mu_p * h2 * dzinv2 ) / rho;
}

/**
 * Largest value of component n over cell (i,j,k) and its 26 neighbours.  The diffusion
 * operators average the cell diffusivities to faces and edges, so a row of the operator
 * never carries more than this.  With NeighbourhoodMinPositive for the density and the
 * coefficients frozen, the scalar rate is a Gershgorin bound on the centred interior stencil
 * of the Cartesian operator (the higher-order one-sided gradients next to Dirichlet
 * boundaries, ERF_ScalarDiffusion.H, whose row is about 4/3 larger, are not covered); for
 * momentum it is a bound for
 * uniform coefficients when mu_h = mu_v (tight on
 * cubic cells) and an estimate otherwise, since the operator is then not symmetric.
 * Needs one ghost cell.
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
NeighbourhoodMax (int i, int j, int k, int n, const amrex::Array4<const amrex::Real>& a) noexcept
{
    amrex::Real m = a(i,j,k,n);
    for (int kk = -1; kk <= 1; ++kk) {
    for (int jj = -1; jj <= 1; ++jj) {
    for (int ii = -1; ii <= 1; ++ii) {
        const amrex::Real v = a(i+ii,j+jj,k+kk,n);
        m = (v > m) ? v : m;
    }}}
    return m;
}

/**
 * Smallest positive value of component n over cell (i,j,k) and its 26 neighbours; ghost
 * cells that hold no density (<= 0) are skipped.  Returns `fallback` if none is positive.
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
NeighbourhoodMinPositive (int i, int j, int k, int n, const amrex::Array4<const amrex::Real>& a,
                          amrex::Real fallback) noexcept
{
    amrex::Real m = fallback;
    bool found = false;
    for (int kk = -1; kk <= 1; ++kk) {
    for (int jj = -1; jj <= 1; ++jj) {
    for (int ii = -1; ii <= 1; ++ii) {
        const amrex::Real v = a(i+ii,j+jj,k+kk,n);
        if (v > zero && (!found || v < m)) { m = v; found = true; }
    }}}
    return m;
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

/** Which components carry diffusion, and their explicit vertical fractions, for CellDiffusiveRates. */
struct DiffusiveRateSettings
{
    bool variable_dz = false;  // terrain-fitted or stretched mesh: dz and slope from z_nd
    bool has_q       = false;  // moisture carried
    bool has_ke      = false;  // turbulent kinetic energy carried
    bool has_scalar  = false;  // advected scalar carried
    amrex::Real e_uv = one, e_w = one, e_th = one, e_q = one, e_ke = one;
};

/**
 * Explicit eddy-diffusion rates of cell (i,j,k), momentum and scalars (1/s), as the
 * diffusive time-step check reduces them: MomentumDiffusiveRate and ScalarDiffusiveRate with
 * the largest diffusivities and the smallest positive density over the cell and its 26
 * neighbours.  Both are zero in a cell without density (covered).
 *
 * @param[in]  s_arr  conserved state (Rho_comp), one ghost cell
 * @param[in]  k_arr  eddy diffusivities (EddyDiff components), one ghost cell
 * @param[in]  z_nd   nodal heights (read only when set.variable_dz)
 * @param[in]  mf_mx, mf_my map factors of the cell
 * @param[in]  dxinv  inverse cell sizes
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
CellDiffusiveRates (int i, int j, int k,
                    const amrex::Array4<const amrex::Real>& s_arr,
                    const amrex::Array4<const amrex::Real>& k_arr,
                    const amrex::Array4<const amrex::Real>& z_nd,
                    amrex::Real mf_mx, amrex::Real mf_my,
                    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxinv,
                    const DiffusiveRateSettings& set,
                    amrex::Real& rate_mom, amrex::Real& rate_scal) noexcept
{
    const amrex::Real ax = mf_mx * dxinv[0];
    const amrex::Real ay = mf_my * dxinv[1];
    amrex::Real dzinv = dxinv[2];
    amrex::Real h2    = zero;
    if (set.variable_dz) {
        const TerrainCellDrops d = ComputeTerrainCellDrops(i,j,k,z_nd);
        dzinv = one / d.dz;
        const amrex::Real hx = ax * d.dx_drop;
        const amrex::Real hy = ay * d.dy_drop;
        h2 = hx*hx + hy*hy;
    }
    const amrex::Real rho = s_arr(i,j,k,Rho_comp);
    if (!(rho > zero)) { rate_mom = zero; rate_scal = zero; return; }
    const amrex::Real rdiv   = NeighbourhoodMinPositive(i,j,k,Rho_comp,s_arr,rho);
    const amrex::Real a2     = ax*ax;
    const amrex::Real b2     = ay*ay;
    const amrex::Real dzinv2 = dzinv*dzinv;
    rate_mom = MomentumDiffusiveRate(rdiv, NeighbourhoodMax(i,j,k,EddyDiff::Mom_h,k_arr),
                                     NeighbourhoodMax(i,j,k,EddyDiff::Mom_v,k_arr),
                                     a2, b2, h2, dzinv2, set.e_uv, set.e_w);
    amrex::Real rs = ScalarDiffusiveRate(rdiv, NeighbourhoodMax(i,j,k,EddyDiff::Theta_h,k_arr),
                                         NeighbourhoodMax(i,j,k,EddyDiff::Theta_v,k_arr),
                                         a2, b2, h2, dzinv2, set.e_th);
    if (set.has_q) {
        rs = amrex::max(rs, ScalarDiffusiveRate(rdiv, NeighbourhoodMax(i,j,k,EddyDiff::Q_h,k_arr),
                                                NeighbourhoodMax(i,j,k,EddyDiff::Q_v,k_arr),
                                                a2, b2, h2, dzinv2, set.e_q));
    }
    if (set.has_ke) {
        rs = amrex::max(rs, ScalarDiffusiveRate(rdiv, NeighbourhoodMax(i,j,k,EddyDiff::KE_h,k_arr),
                                                NeighbourhoodMax(i,j,k,EddyDiff::KE_v,k_arr),
                                                a2, b2, h2, dzinv2, set.e_ke));
    }
    if (set.has_scalar) {
        rs = amrex::max(rs, ScalarDiffusiveRate(rdiv, NeighbourhoodMax(i,j,k,EddyDiff::Scalar_h,k_arr),
                                                NeighbourhoodMax(i,j,k,EddyDiff::Scalar_v,k_arr),
                                                a2, b2, h2, dzinv2, one));
    }
    rate_scal = rs;
}

/**
 * Whether the terrain slope-factor report is evaluated this coarse step: always on a new set
 * of grids, and every step on a moving terrain (MovingFittedMesh), where alpha changes with time.
 */
inline bool
SlopeReportEvaluate (bool new_grids, bool moving_terrain)
{
    return new_grids || moving_terrain;
}

/**
 * Whether an evaluated slope-factor report is printed: always when forced (new grids),
 * otherwise only if the largest alpha moved by more than 1 % from the last printed value.
 * alpha >= 1 by construction, so the relative test never fires on changes below 0.01.
 */
inline bool
SlopeReportPrint (bool force, amrex::Real alpha_max, amrex::Real alpha_last)
{
    return force || std::abs(alpha_max - alpha_last) > amrex::Real(0.01) * alpha_last;
}

/**
 * Explicit fraction of the vertical diffusion for the diffusive time-step check: 0 if every
 * stage of the partly implicit Runge-Kutta scheme stays bounded for the vertical diffusion,
 * 1 otherwise (including every pattern for which that has not been shown).
 *
 * Each stage k starts from S^n with stage step c_k dt; its explicit part (1 - f_k) acts on
 * the previous stage's state and its implicit part f_k on its own, so for dS/dt = lambda S
 * and y = -lambda dt >= 0,  S_k = (1 - (1 - f_k) c_k y S_{k-1}) / (1 + f_k c_k y),  S_0 = 1.
 * Each stage's state feeds the next stage's explicit terms, so every |S_k| must stay <= 1,
 * not only the last.  The result is decided analytically:
 *
 *  - three stages (c = 1/3, 1/2, 1; compressible): bounded for f_1 >= 1/2, f_2 = 1 and any
 *    f_3 in [0,1].  |S_1| <= 1 iff f_1 >= 1/2: |1 - (1-f) x| <= 1 + (1-f) x <= 1 + f x for
 *    f >= 1/2, and for f < 1/2, S_1 -> -(1-f)/f < -1 as x -> infinity;
 *    S_2 = 1 / (1 + y/2); and S_3 = (1 + (f_3 - 1/2) y) / ((1 + y/2)(1 + f_3 y)), whose
 *    numerator is at most 1 + y/2 in size.  The default 1 1 0 is Crank-Nicolson, 1 1 1
 *    backward Euler.
 *  - two stages (c = 1/2, 1; anelastic MidPoint): bounded for f_1 = 1 and any f_2 in [0,1]
 *    (the same S_3 form); with f_1 < 1 and f_2 = 0, as MidPoint forces, the second stage
 *    grows without bound.
 *
 * Every other pattern returns 1: some are not bounded (0 0 1, 1 0 0, and 0.995 0.1 0.0666,
 * whose last stage reaches 1.0008 near y = 23.5), and the rest are not shown to be.
 *
 * @param[in] nstages  3 for the compressible scheme, 2 for anelastic MidPoint
 * @param[in] fac      implicit factor f_k of each stage (erf.vert_implicit_fac)
 * @param[in] implicit whether the implicit solve treats this component at all
 */
template <typename V>
amrex::Real
ExplicitVerticalFraction (int nstages, const V& fac, bool implicit)
{
    if (!implicit || static_cast<int>(fac.size()) < nstages) { return one; }
    // The proof assumes factors in [0,1]; anything else (the inputs are not range-checked,
    // and 1 1.5 0 amplifies by 4/3) counts as explicit rather than being clamped.
    for (int k = 0; k < nstages; ++k) {
        const double fk = static_cast<double>(fac[k]);
        if (!(fk >= 0.0 && fk <= 1.0)) { return one; }
    }
    auto f = [&] (int k) { return static_cast<double>(fac[k]); };
    if (nstages == 3) {
        return (f(0) >= 0.5 && f(1) >= 1.0) ? zero : one;
    }
    if (nstages == 2) {
        return (f(0) >= 1.0) ? zero : one;
    }
    return one;
}

#endif
