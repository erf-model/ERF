#ifndef ERF_IMMERSED_WALL_CELL_H_
#define ERF_IMMERSED_WALL_CELL_H_

#include <cmath>
#include <limits>

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Algorithm.H>

#include "ERF_Constants.H"
#include "ERF_MOSTUtils.H"

/**
 * The wall cell of immersed-forcing terrain and the fraction-stress wall law
 * (erf.if_wall_form = fraction_stress), shared by the momentum sources, the wall-face
 * masks of the diffusion and the RANS wall distance so they cannot disagree.
 *
 * A column of solid fractions beta (terrain_blanking) meets the air in its WALL CELL: a
 * partial cell with fluid above it, or a fluid cell sitting on a solid cell. The surface is
 * at h = z_bot + beta dz of that cell (beta = 0 for the fluid cell), the wall cell's fluid
 * centroid at d1 = (1 - beta) dz / 2 above it and the centre of the cell above at
 * d2 = (3/2 - beta) dz. The fraction-stress law takes the friction velocity
 *     u*^2 = (1 - beta) u*_own^2 + beta u*_ref^2,
 *     u*_own = kappa |u_t| / ln(d1/z0),   u*_ref = kappa |u_t,ref| / ln(d2/z0),
 * so a nearly full fluid cell uses its own velocity (on a face-aligned surface, beta = 0,
 * that is the plain surface-layer law of flat ground at dz/2) and a thin fluid sliver the
 * better resolved cell above, and applies it as a wall stress spread over the cell height,
 * integrated exactly over the step:
 *     F = -rho C u_t,  C = (1 - exp(-r dt)) / dt,  r = u*^2 / (dz |u_t|).
 * The resolved diffusive flux through the wall face (the bottom face of the wall cell) is
 * switched off, so the column loses exactly u*^2 through the wall. This is form 7 of the
 * kynema-sgf immersed terrain, written for the vertical columns of ERF's immersed forcing
 * (horizontal velocity, vertical distances) and neutral stratification.
 */
namespace ib_wall {

/** A cell whose solid fraction is below small counts as fluid */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool is_fluid (amrex::Real beta, amrex::Real small) noexcept
{
    return beta < small;
}

/** A cell whose solid fraction is within small of one counts as solid */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool is_solid (amrex::Real beta, amrex::Real small) noexcept
{
    return beta >= amrex::Real(1.0) - small;
}

/**
 * Cell k, with solid fractions b_km1, b_k, b_kp1 in cells k-1, k, k+1, is the wall cell:
 * a partial cell with fluid above it, or a fluid cell on top of a solid cell.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool is_wall_cell (amrex::Real b_km1, amrex::Real b_k, amrex::Real b_kp1, amrex::Real small) noexcept
{
    const bool partial = !is_fluid(b_k, small) && !is_solid(b_k, small);
    return (partial && is_fluid(b_kp1, small)) || (is_fluid(b_k, small) && is_solid(b_km1, small));
}

/** Solid fraction of a wall cell as the wall law uses it (zero for a fluid wall cell) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real wall_beta (amrex::Real b_k, amrex::Real small) noexcept
{
    return is_fluid(b_k, small) ? amrex::Real(0.0) : amrex::min(b_k, amrex::Real(1.0));
}

/**
 * Friction velocity of the fraction-stress law.
 * @param ut     horizontal speed in the wall cell
 * @param ut_ref horizontal speed in the cell above
 * @param beta   wall_beta of the wall cell
 * @param dz     cell height
 * @param z0     roughness length
 * The distances are floored at 2 z0 so that ln(d/z0) stays above ln 2 in a thin sliver,
 * where the (1 - beta) weight of the own estimate is small anyway.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real friction_velocity (amrex::Real ut, amrex::Real ut_ref, amrex::Real beta,
                               amrex::Real dz, amrex::Real z0) noexcept
{
    const amrex::Real d_min = amrex::Real(2.0) * z0;
    const amrex::Real d1 = amrex::max(amrex::Real(0.5) * (amrex::Real(1.0) - beta) * dz, d_min);
    const amrex::Real d2 = amrex::max((amrex::Real(1.5) - beta) * dz, d1);
    const amrex::Real us_own = KAPPA * ut     / std::log(d1 / z0);
    const amrex::Real us_ref = KAPPA * ut_ref / std::log(d2 / z0);
    return std::sqrt((amrex::Real(1.0) - beta) * us_own * us_own + beta * us_ref * us_ref);
}


/**
 * Factor on the diffusive coupling through the face above the wall cell (the face between the
 * wall cell, whose fluid centroid sits d1 = (1-beta) dz/2 above the surface, and the cell above,
 * at d2 = (3/2-beta) dz). The solver differences across that face over dz; the factor dz/dz_eff
 * replaces dz by
 *   spacing = 0 (grid):     dz                                    (factor 1)
 *   spacing = 1 (physical): d2 - d1 = (1 - beta/2) dz
 *   spacing = 2 (log):      (d1 + d2) ln(d2/d1) / (2 ln 3), the spacing that gives the face the
 *                           same relative error for a log profile as the first interior face of
 *                           flat ground (d1 = dz/2, d2 = 3 dz/2), where it equals dz
 * All three are 1 for beta = 0, so a surface on a cell face stays flat ground.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real face_above_factor (amrex::Real beta, amrex::Real dz, amrex::Real z0, int spacing) noexcept
{
    const amrex::Real d1 = amrex::max(amrex::Real(0.5) * (amrex::Real(1.0) - beta) * dz, amrex::Real(2.0) * z0);
    const amrex::Real d2 = amrex::max((amrex::Real(1.5) - beta) * dz, d1);
    const amrex::Real dz_phys = d2 - d1;
    const amrex::Real dz_log  = (d1 + d2) * std::log(d2 / d1) / (amrex::Real(2.0) * amrex::Real(1.0986122886681098));
    const amrex::Real dz_eff  = (spacing == 1) ? dz_phys : ((spacing == 2) ? dz_log : dz);
    return dz / amrex::max(dz_eff, std::numeric_limits<amrex::Real>::min());
}

/** Surface condition of the fraction-stress law (erf.if_* inputs of the terrain immersed forcing) */
struct WallCond
{
    static constexpr int neutral     = 0;   //!< no heat exchange, neutral stability functions
    static constexpr int heat_flux   = 1;   //!< given kinematic surface heat flux q [K m/s]
    static constexpr int surface_temp = 2;  //!< given surface potential temperature theta_s [K]
    static constexpr int obukhov     = 3;   //!< given Obukhov length (inverse Linv) [1/m]
    int type{neutral};
    amrex::Real q{0.0};
    amrex::Real theta_s{0.0};
    amrex::Real Linv{0.0};
};

/** Friction velocity, temperature scale and inverse Obukhov length of a wall cell */
struct WallState
{
    amrex::Real ustar{0.0};
    amrex::Real thetastar{0.0};
    amrex::Real Linv{0.0};
};

/** ln(d/z0) - psi, floored at ln 2 (the value the 2 z0 distance floor gives in neutral air) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real phi (amrex::Real d, amrex::Real z0, amrex::Real psi) noexcept
{
    return amrex::max(std::log(d / z0) - psi, amrex::Real(0.6931471805599453));
}

/**
 * Wall state of the fraction-stress law with MOST stability functions (the Dyer functions of
 * the flat surface layer, similarity_funs::calc_psi_m / calc_psi_h, zeta clamped to +-100):
 *     u*^2     = (1-b) (kappa ut     / phi_m(d1))^2 + b (kappa ut_ref / phi_m(d2))^2
 *     theta*   = -q / u*                                   (heat_flux)
 *              = (1-b) kappa (theta - theta_s) / phi_h(d1)
 *              +   b   kappa (theta_ref - theta_s) / phi_h(d2)  (surface_temp)
 *              = theta u*^2 / (kappa g L)                  (obukhov)
 *     1/L      = kappa g theta* / (u*^2 theta)             (heat_flux, surface_temp: iterated)
 * with d1 = (1-b) dz/2, d2 = (3/2-b) dz floored at 2 z0 and phi_{m,h}(d) = ln(d/z0) - psi(d/L).
 * The iteration starts from neutral and runs a fixed number of passes, so it is deterministic
 * and has no data-dependent trip count.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
WallState wall_state (amrex::Real ut, amrex::Real ut_ref, amrex::Real theta, amrex::Real theta_ref,
                      amrex::Real beta, amrex::Real dz, amrex::Real z0, const WallCond& c) noexcept
{
    similarity_funs sfuns;
    const amrex::Real one_r  = amrex::Real(1.0);
    const amrex::Real d_min  = amrex::Real(2.0) * z0;
    const amrex::Real d1 = amrex::max(amrex::Real(0.5) * (one_r - beta) * dz, d_min);
    const amrex::Real d2 = amrex::max((amrex::Real(1.5) - beta) * dz, d1);
    const amrex::Real zmax   = amrex::Real(100.0);
    const amrex::Real us_min = amrex::Real(1.0e-3);
    const amrex::Real kg     = KAPPA * CONST_GRAV;
    const int n_iter = (c.type == WallCond::heat_flux || c.type == WallCond::surface_temp) ? 12 : 1;

    WallState w;
    w.Linv = (c.type == WallCond::obukhov) ? c.Linv : amrex::Real(0.0);
    for (int it = 0; it < n_iter; ++it) {
        const amrex::Real z1 = amrex::max(amrex::min(d1 * w.Linv, zmax), -zmax);
        const amrex::Real z2 = amrex::max(amrex::min(d2 * w.Linv, zmax), -zmax);
        const amrex::Real us_own = KAPPA * ut     / phi(d1, z0, sfuns.calc_psi_m(z1));
        const amrex::Real us_ref = KAPPA * ut_ref / phi(d2, z0, sfuns.calc_psi_m(z2));
        w.ustar = std::sqrt((one_r - beta) * us_own * us_own + beta * us_ref * us_ref);
        const amrex::Real us_f = amrex::max(w.ustar, us_min);
        if (c.type == WallCond::heat_flux) {
            w.thetastar = -c.q / us_f;
        } else if (c.type == WallCond::surface_temp) {
            w.thetastar = (one_r - beta) * KAPPA * (theta     - c.theta_s) / phi(d1, z0, sfuns.calc_psi_h(z1))
                        +          beta  * KAPPA * (theta_ref - c.theta_s) / phi(d2, z0, sfuns.calc_psi_h(z2));
        } else if (c.type == WallCond::obukhov) {
            w.thetastar = theta * us_f * us_f * c.Linv / kg;
        } else {
            w.thetastar = amrex::Real(0.0);
        }
        if (c.type == WallCond::heat_flux || c.type == WallCond::surface_temp) {
            w.Linv = kg * w.thetastar / (us_f * us_f * theta);
        }
    }
    return w;
}

/**
 * Source of the wall heat flux in the wall cell, per unit mass (d theta / dt): the kinematic
 * flux -u* theta* spread over the cell height. With a given surface temperature theta* is
 * proportional to theta - theta_s and the term is a decay, integrated exactly over dt so it
 * cannot overshoot the surface temperature; otherwise the flux is applied as it stands.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real heat_source (const WallState& w, amrex::Real theta, const WallCond& c,
                         amrex::Real dz, amrex::Real dt) noexcept
{
    const amrex::Real flux = -w.ustar * w.thetastar / dz;
    if (c.type != WallCond::surface_temp) { return flux; }
    const amrex::Real dth  = theta - c.theta_s;
    const amrex::Real eps  = std::numeric_limits<amrex::Real>::epsilon();
    const amrex::Real rate = -flux / (dth + std::copysign(eps, dth));
    const amrex::Real dt_safe = amrex::max(dt, std::numeric_limits<amrex::Real>::min());
    // a decay (rate > 0) is integrated exactly; anything else is applied as it stands
    const amrex::Real exact = std::expm1(-amrex::max(rate, amrex::Real(0.0)) * dt_safe) / dt_safe * dth;
    return (rate > amrex::Real(0.0) && std::abs(dth) > eps) ? exact : flux;
}

/**
 * Effective rate C of the wall stress u*^2 spread over the cell height dz, integrated exactly
 * over dt: the momentum source is -rho C u_t.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real stress_rate (amrex::Real ustar, amrex::Real ut, amrex::Real dz, amrex::Real dt) noexcept
{
    // dt is floored rather than branched on, so no guarded 0/0 can be speculated (FPE traps);
    // (1 - exp(-r dt)) / dt tends to r as dt -> 0
    const amrex::Real r = ustar * ustar / (dz * (ut + std::numeric_limits<amrex::Real>::epsilon()));
    const amrex::Real dt_safe = amrex::max(dt, std::numeric_limits<amrex::Real>::min());
    return -std::expm1(-r * dt_safe) / dt_safe;
}

/**
 * Wall value of the k-eqn TKE in the wall cell (erf.if_wall_tke), the value erf.dirichlet_k holds in
 * the first cell over flat ground (Axell & Liungman 2001, Eq. 16):
 *   k = (u*^3 + kappa B d1)^(2/3) / Cmu0^2,  B = -g u* theta* / theta (unstable only),
 * with u*, theta* from the wall state and d1 the fluid-centroid height of the wall cell (floored
 * at 2 z0, as the wall state does). kynema-sgf relaxes its immersed wall cells toward the same value.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real wall_tke (const WallState& w, amrex::Real beta, amrex::Real dz, amrex::Real z0,
                      amrex::Real theta, amrex::Real Cmu0) noexcept
{
    const amrex::Real d1 = amrex::max(amrex::Real(0.5) * (amrex::Real(1.0) - beta) * dz,
                                      amrex::Real(2.0) * z0);
    const amrex::Real B  = (w.thetastar < amrex::Real(-1.e-8))
                         ? -CONST_GRAV * w.ustar * w.thetastar / theta : amrex::Real(0.0);
    return std::pow(w.ustar * w.ustar * w.ustar + KAPPA * B * d1, amrex::Real(2.0) / amrex::Real(3.0))
         / (Cmu0 * Cmu0);
}

} // namespace ib_wall

#endif
