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
 *   spacing = 2 (log):      (d1 + d2) ln(d2/d1) / (2 ln 3), the spacing that gives the face flux,
 *                           with the face viscosity the mean of the two cells' kappa u* z, the
 *                           same relative error for a log profile as the first interior face of
 *                           flat ground (d1 = dz/2, d2 = 3 dz/2), where it equals dz: ln 3 u*^2
 *                           for every beta
 * All three are 1 for beta = 0, so a surface on a cell face stays flat ground. The factor is
 * bounded only while z0 < dz/4 (d2 > d1 for every beta): the inputs require that.
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
    similarity_funs sfuns{};  //!< the similarity functions (erf.most.stable_functions)
    // Surface moisture condition (a moisture model with water vapour): none, a given kinematic
    // moisture flux, or a given surface mixing ratio
    static constexpr int moist_none    = 0;   //!< no moisture exchange, q* = 0
    static constexpr int moist_flux    = 1;   //!< given kinematic surface moisture flux [kg/kg m/s]
    static constexpr int surface_moist = 2;   //!< given surface mixing ratio q_s [kg/kg]
    int moist_type{moist_none};
    amrex::Real qflux{0.0};
    amrex::Real q_s{0.0};
    //! a moisture model is active: the Obukhov length is formed from the virtual potential
    //! temperature and its scale theta_v*, as the flat surface layer forms it
    bool moist{false};
};

/** Friction velocity, temperature and moisture scales and inverse Obukhov length of a wall cell */
struct WallState
{
    amrex::Real ustar{0.0};
    amrex::Real thetastar{0.0};
    amrex::Real qstar{0.0};
    amrex::Real Linv{0.0};
};

/**
 * Virtual potential temperature of the wall law, theta (1 + epsv q), water vapour only (the flat
 * surface layer's plane average subtracts rain water too; the two agree without rain)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real virtual_theta (amrex::Real theta, amrex::Real q) noexcept
{
    return theta * (amrex::Real(1.0) + epsv * q);
}

/**
 * Scale of the virtual heat flux, theta_v* = theta* (1 + epsv q) + epsv theta q*, so that
 * -u* theta_v* = (1 + epsv q) <w'theta'> + epsv theta <w'q'> is the buoyancy flux (the flat
 * surface layer's tstv)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real virtual_thetastar (const WallState& w, amrex::Real theta, amrex::Real q) noexcept
{
    return w.thetastar * (amrex::Real(1.0) + epsv * q) + epsv * theta * w.qstar;
}

/** ln(d/z0) - psi, floored at ln 2 (the value the 2 z0 distance floor gives in neutral air) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real phi (amrex::Real d, amrex::Real z0, amrex::Real psi) noexcept
{
    return amrex::max(std::log(d / z0) - psi, amrex::Real(0.6931471805599453));
}

/**
 * Wall state of the fraction-stress law with MOST stability functions (the functions of the
 * flat surface layer, similarity_funs::calc_psi_m / calc_psi_h with Dyer's or QNSE's stable
 * side as erf.most.stable_functions says, zeta clamped to +-100):
 *     u*^2     = (1-b) (kappa ut     / phi_m(d1))^2 + b (kappa ut_ref / phi_m(d2))^2
 *     theta*   = -q / u*                                   (heat_flux)
 *              = (1-b) kappa (theta - theta_s) / phi_h(d1)
 *              +   b   kappa (theta_ref - theta_s) / phi_h(d2)  (surface_temp)
 *              = theta u*^2 / (kappa g L)                  (obukhov)
 *     q*       = -F_q / u*                                 (moist_flux)
 *              = (1-b) kappa (q - q_s) / phi_h(d1)
 *              +   b   kappa (q_ref - q_s) / phi_h(d2)      (surface_moist)
 *     1/L      = kappa g theta* / (u*^2 theta)             (heat_flux, surface_temp: iterated)
 *              = kappa g theta_v* / (u*^2 theta_v)          (the same with a moisture model, c.moist:
 *                                                           theta_v = theta (1 + epsv q),
 *                                                           theta_v* = theta* (1 + epsv q) + epsv theta q*,
 *                                                           iterated also for a moisture condition alone)
 * with d1 = (1-b) dz/2, d2 = (3/2-b) dz floored at 2 z0 and phi_{m,h}(d) = ln(d/z0) - psi(d/L).
 * With c.moist a given L (obukhov) sets theta_v*, and theta* = (theta_v* - epsv theta q*) / (1 + epsv q).
 * q and q_ref are the water vapour mixing ratios of the wall cell and the cell above (ignored
 * without c.moist). The iteration starts from neutral and runs a fixed number of passes, so it is
 * deterministic and has no data-dependent trip count.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
WallState wall_state (amrex::Real ut, amrex::Real ut_ref, amrex::Real theta, amrex::Real theta_ref,
                      amrex::Real q, amrex::Real q_ref,
                      amrex::Real beta, amrex::Real dz, amrex::Real z0, const WallCond& c) noexcept
{
    const similarity_funs& sfuns = c.sfuns;
    const amrex::Real one_r  = amrex::Real(1.0);
    const amrex::Real d_min  = amrex::Real(2.0) * z0;
    const amrex::Real d1 = amrex::max(amrex::Real(0.5) * (one_r - beta) * dz, d_min);
    const amrex::Real d2 = amrex::max((amrex::Real(1.5) - beta) * dz, d1);
    const amrex::Real zmax   = amrex::Real(100.0);
    const amrex::Real us_min = amrex::Real(1.0e-3);
    const amrex::Real kg     = KAPPA * CONST_GRAV;
    // 1/L is solved for (iterated) with a heat condition that is not a given length, and with a
    // moisture condition alone when the moisture enters the buoyancy
    const bool moist_cond = c.moist && (c.moist_type != WallCond::moist_none);
    const bool solve_L = (c.type == WallCond::heat_flux || c.type == WallCond::surface_temp ||
                          (moist_cond && c.type != WallCond::obukhov));
    const int n_iter = solve_L ? 12 : 1;

    WallState w;
    w.Linv = (c.type == WallCond::obukhov) ? c.Linv : amrex::Real(0.0);
    // 1/L is kept where d1/L stays within the zeta clamp: the stable given-flux problem has no
    // solution below a wind threshold (about 3 m/s at q = -0.02 K m/s, d1 = 10 m, z0 = 0.1 m)
    // and the iteration would otherwise run 1/L away while u* collapses to the floor
    const amrex::Real Linv_max = zmax / d1;
    for (int it = 0; it < n_iter; ++it) {
        w.Linv = amrex::max(amrex::min(w.Linv, Linv_max), -Linv_max);
        const amrex::Real z1 = amrex::max(amrex::min(d1 * w.Linv, zmax), -zmax);
        const amrex::Real z2 = amrex::max(amrex::min(d2 * w.Linv, zmax), -zmax);
        const amrex::Real us_own = KAPPA * ut     / phi(d1, z0, sfuns.calc_psi_m(z1));
        const amrex::Real us_ref = KAPPA * ut_ref / phi(d2, z0, sfuns.calc_psi_m(z2));
        w.ustar = std::sqrt((one_r - beta) * us_own * us_own + beta * us_ref * us_ref);
        const amrex::Real us_f = amrex::max(w.ustar, us_min);
        if (c.moist_type == WallCond::moist_flux) {
            w.qstar = -c.qflux / us_f;
        } else if (c.moist_type == WallCond::surface_moist) {
            w.qstar = (one_r - beta) * KAPPA * (q     - c.q_s) / phi(d1, z0, sfuns.calc_psi_h(z1))
                    +          beta  * KAPPA * (q_ref - c.q_s) / phi(d2, z0, sfuns.calc_psi_h(z2));
        } else {
            w.qstar = amrex::Real(0.0);
        }
        if (c.type == WallCond::heat_flux) {
            w.thetastar = -c.q / us_f;
        } else if (c.type == WallCond::surface_temp) {
            w.thetastar = (one_r - beta) * KAPPA * (theta     - c.theta_s) / phi(d1, z0, sfuns.calc_psi_h(z1))
                        +          beta  * KAPPA * (theta_ref - c.theta_s) / phi(d2, z0, sfuns.calc_psi_h(z2));
        } else if (c.type == WallCond::obukhov) {
            if (c.moist) {
                // the given length sets theta_v*; theta* is what is left after the moisture part
                const amrex::Real thv_star = virtual_theta(theta, q) * us_f * us_f * c.Linv / kg;
                w.thetastar = (thv_star - epsv * theta * w.qstar) / (one_r + epsv * q);
            } else {
                w.thetastar = theta * us_f * us_f * c.Linv / kg;
            }
        } else {
            w.thetastar = amrex::Real(0.0);
        }
        if (solve_L) {
            w.Linv = c.moist ? kg * virtual_thetastar(w, theta, q) / (us_f * us_f * virtual_theta(theta, q))
                             : kg * w.thetastar / (us_f * us_f * theta);
        }
    }
    // the stress uses the same floored u* as theta* and 1/L
    w.ustar = amrex::max(w.ustar, us_min);
    w.Linv  = amrex::max(amrex::min(w.Linv, Linv_max), -Linv_max);
    return w;
}

/** The wall state without moisture (q = q_ref = 0; c.moist and c.moist_type are expected unset) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
WallState wall_state (amrex::Real ut, amrex::Real ut_ref, amrex::Real theta, amrex::Real theta_ref,
                      amrex::Real beta, amrex::Real dz, amrex::Real z0, const WallCond& c) noexcept
{
    return wall_state(ut, ut_ref, theta, theta_ref, amrex::Real(0.0), amrex::Real(0.0), beta, dz, z0, c);
}

/**
 * A wall flux (per unit mass, already divided by the cell height) that is proportional to
 * delta = value - surface value, integrated exactly over dt as a decay so the cell cannot
 * overshoot the surface value; a flux that is not a decay (rate <= 0) is applied as it stands.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real decay_toward (amrex::Real flux, amrex::Real delta, amrex::Real dt) noexcept
{
    const amrex::Real eps  = std::numeric_limits<amrex::Real>::epsilon();
    const amrex::Real rate = -flux / (delta + std::copysign(eps, delta));
    const amrex::Real dt_safe = amrex::max(dt, std::numeric_limits<amrex::Real>::min());
    // a decay (rate > 0) is integrated exactly; anything else is applied as it stands
    const amrex::Real exact = std::expm1(-amrex::max(rate, amrex::Real(0.0)) * dt_safe) / dt_safe * delta;
    return (rate > amrex::Real(0.0) && std::abs(delta) > eps) ? exact : flux;
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
    return decay_toward(flux, theta - c.theta_s, dt);
}

/**
 * Source of the wall moisture flux in the wall cell, per unit mass (d q / dt): the kinematic
 * flux -u* q* spread over the cell height; with a given surface mixing ratio the term is a
 * decay toward q_s, integrated exactly as heat_source does toward theta_s (as there, the
 * reference-cell part of a partial cell's q* is folded into the decay rate, so the wall cell
 * never passes q_s). Zero without a moisture condition (q* = 0).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real moist_source (const WallState& w, amrex::Real q, const WallCond& c,
                          amrex::Real dz, amrex::Real dt) noexcept
{
    const amrex::Real flux = -w.ustar * w.qstar / dz;
    if (c.moist_type != WallCond::surface_moist) { return flux; }
    return decay_toward(flux, q - c.q_s, dt);
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
 *   k = (u*^3 + kappa B d1)^(2/3) / Cmu0^2,  B = -g u* bstar / btheta (unstable only),
 * with u* from the wall state, the buoyancy scale bstar and its temperature btheta (theta* and
 * theta; with a moisture model theta_v* and theta_v, see virtual_thetastar) and d1 the
 * fluid-centroid height of the wall cell (floored at 2 z0, as the wall state does). kynema-sgf
 * relaxes its immersed wall cells toward the same value.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real wall_tke (const WallState& w, amrex::Real beta, amrex::Real dz, amrex::Real z0,
                      amrex::Real Cmu0, amrex::Real bstar, amrex::Real btheta) noexcept
{
    const amrex::Real d1 = amrex::max(amrex::Real(0.5) * (amrex::Real(1.0) - beta) * dz,
                                      amrex::Real(2.0) * z0);
    const amrex::Real B  = (bstar < amrex::Real(-1.e-8))
                         ? -CONST_GRAV * w.ustar * bstar / btheta : amrex::Real(0.0);
    return std::pow(w.ustar * w.ustar * w.ustar + KAPPA * B * d1, amrex::Real(2.0) / amrex::Real(3.0))
         / (Cmu0 * Cmu0);
}

/** The dry wall value: B from theta* and theta */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real wall_tke (const WallState& w, amrex::Real beta, amrex::Real dz, amrex::Real z0,
                      amrex::Real theta, amrex::Real Cmu0) noexcept
{
    return wall_tke(w, beta, dz, z0, Cmu0, w.thetastar, theta);
}

} // namespace ib_wall

#endif
