#ifndef ERF_MOSTUtils_H
#define ERF_MOSTUtils_H

#include <AMReX_Algorithm.H>
#include <AMReX_BLassert.H>
#include <AMReX_Math.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <string>
#include "ERF_Constants.H"

/**
 * Structure of plain old data relevant to MOST BCs
 */
struct most_data
{
public:
    amrex::Real z0_const{amrex::Real(0.1)};    ///< Roughness height -- default constant value(m)
    amrex::Real kappa{KAPPA};                  ///< von Karman constant
    amrex::Real gravity{CONST_GRAV};           ///< Acceleration due to gravity (m/s^2)
    amrex::Real surf_temp_flux{zero};          ///< Heat flux TODO: decide whether this is <θ'w'> or <θv'w'> under moist conditions
    amrex::Real surf_moist_flux{zero};         ///< Moisture flux

    amrex::Real Cnk_a{amrex::Real(0.0185)};    ///< Standard Charnock constant https://doi.org/amrex::Real(10.1175)/JAMC-D-17-amrex::Real(0137.1)
    amrex::Real Cnk_b1{one/amrex::Real(30.0)}; ///< Modified Charnock Eq (4) https://doi.org/amrex::Real(10.1175)/JAMC-D-17-amrex::Real(0137.1)
    amrex::Real Cnk_b2{amrex::Real(1260.0)};   ///< Modified Charnock Eq (4) https://doi.org/amrex::Real(10.1175)/JAMC-D-17-amrex::Real(0137.1)
    amrex::Real Cnk_d{amrex::Real(30.0)};      ///< Modified Charnock Eq (4) https://doi.org/amrex::Real(10.1175)/JAMC-D-17-amrex::Real(0137.1)
    amrex::Real Cnk_b = Cnk_b1 * std::log(Cnk_b2 / Cnk_d); ///< Modified Charnock derived quantity
    bool visc{true};                           ///< Use viscous Charnock formulation

    amrex::Real Cd{amrex::Real(0.001)};
    amrex::Real Ch{amrex::Real(0.001)};
    amrex::Real Cq{amrex::Real(0.001)};

    const amrex::Real Bjr_beta = amrex::Real(1.2); // Empirical parameter from Beljaars 1995 QJRMS
};


/**
 * Lower bound of the log-law bracket ln(zref/z0) - psi of the MOST flux kernels: C/2 held
 * between 1/2 and 1.
 *
 * The bound is 1 for zref >= e^2 z0 and C/2 down to zref = e z0, so it does not exceed the
 * neutral bracket C wherever zref >= e^(1/2) z0. The lower limit of 1/2 keeps the bracket,
 * and so u*, positive and finite when z0 approaches or exceeds zref (C <= 0).
 *
 * The bound is set by the under-relaxed u* iteration, not by the physics. At an unstable
 * fixed point the relaxed map has slope (1 - alpha) - alpha s / bracket, with
 * s = 3 (z/L) dpsi_m/d(z/L) = 3 (1 - phi_m) between 1.3 (z/L = -0.5) and 2.4 (z/L = -50) for
 * Businger-Dyer.
 * With alpha = 0.5 the iteration converges only for bracket > s/3 (0.4 to 0.8), so a bound of
 * 1 keeps the switch between the bounded and free branches inside the convergent range, while
 * 0.1 C (0.3 at C = 3) does not and the iteration oscillates between two values. If alpha
 * changes, the bound needs revisiting; see PR #4866.
 *
 * @param[in] C ln(zref/z0)
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
most_bracket_floor (amrex::Real C)
{
    return amrex::max(amrex::min(amrex::Real(0.5) * C, amrex::Real(1.0)), amrex::Real(0.5));
}


/**
 * Bound on |z/L| of the MOST kernels (see most_rib_zeta_bracketed, most_clamp_zeta and
 * similarity_funs::calc_zeta_fold).
 */
inline constexpr amrex::Real most_zeta_lim = amrex::Real(100.0);

/**
 * Structure of similarity functions for Moeng formulation.
 *
 * The unstable side is Businger-Dyer. The stable side is Businger-Dyer (the default,
 * phi = 1 + beta zeta, psi = -beta zeta) or, after set_stable_functions, the QNSE form
 * phi = (1 + a zeta)/(1 + b zeta), psi = -((a - b)/b) ln(1 + b zeta) of Sukoriansky, Galperin
 * and Perov (2005), with its own coefficients for momentum and heat (most.stable_functions,
 * read by read_most_similarity_funs). The Jimenez (2012) pair (calc_psi_m2, calc_psi_h2) is
 * not affected by the choice.
 */
struct similarity_funs
{
    //
    // Similarity functions from Jimenez (2012)
    //
    /**
     * Compute the Jimenez momentum stability correction.
     *
     * @param[in] zeta stability coordinate z/L
     */
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    calc_psi_m2 (amrex::Real zeta) const
    {
        if (zeta > 0) {
            amrex::Real x = std::pow(one + std::pow(zeta, amrex::Real(2.5)), one/amrex::Real(2.5));
            return ( -amrex::Real(6.1)*std::log(zeta + x) );
        } else {
            amrex::Real x = std::pow(one - amrex::Real(16.0)*zeta, fourth);
            amrex::Real psi_k_m =  two * std::log(myhalf * (one + x)) + std::log(myhalf * (one + x * x)) -
                                   two * std::atan(x) + PIoTwo;
            amrex::Real y = std::pow(one - amrex::Real(10.0)*zeta, one/three);
            amrex::Real psi_c_m = (three/two)*std::log((y*y + y + one)/three)
                                - std::sqrt(three)*std::atan((two*y + one)/std::sqrt(three))
                                + PI/std::sqrt(three);
            return ( (psi_k_m + zeta*zeta*psi_c_m) / (one + zeta*zeta) );
        }
    }

    /**
     * Compute the Jimenez heat stability correction.
     *
     * @param[in] zeta stability coordinate z/L
     */
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    calc_psi_h2 (amrex::Real zeta) const
    {
        if (zeta > 0) {
            amrex::Real x = std::pow(one + std::pow(zeta, amrex::Real(1.1)), one/amrex::Real(1.1));
            return ( -amrex::Real(5.3)*std::log(zeta + x) );
        } else {
            amrex::Real x = std::sqrt(one - amrex::Real(16.0)*zeta);
            amrex::Real psi_k_h = two * std::log(myhalf * (one + x));
            amrex::Real y = std::pow(one - amrex::Real(34.0)*zeta, one/three);
            amrex::Real psi_c_h = (three/two)*std::log((y*y + y + one)/three)
                                - std::sqrt(three)*std::atan((two*y + one)/std::sqrt(three))
                                + PI/std::sqrt(three);
            return ( (psi_k_h + zeta*zeta*psi_c_h) / (one + zeta*zeta) );
        }
    }


    //
    // Similarity functions from Businger & Dyer (1966); stable side Businger-Dyer or QNSE
    //
    /**
     * Compute the momentum stability correction (Businger-Dyer, or QNSE on the stable side).
     *
     * @param[in] zeta stability coordinate z/L
     */
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    calc_psi_m (amrex::Real zeta) const
    {
        if (zeta > 0) {
            return m_qnse ? qnse_psi(zeta, m_qnse_am, m_qnse_bm) : -beta_m * zeta;
        } else {
            amrex::Real x = std::sqrt(std::sqrt(one - gamma_m * zeta));
            return two * std::log(myhalf * (one + x)) + std::log(myhalf * (one + x * x)) -
                   two * std::atan(x) + PIoTwo;
        }
    }

    /**
     * Compute the heat stability correction (Businger-Dyer, or QNSE on the stable side).
     *
     * @param[in] zeta stability coordinate z/L
     */
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    calc_psi_h (amrex::Real zeta) const
    {
        if (zeta > 0) {
            return m_qnse ? qnse_psi(zeta, m_qnse_ah, m_qnse_bh) : -beta_h * zeta;
        } else {
            amrex::Real x = std::sqrt(one - gamma_h * zeta);
            return two * std::log(myhalf * (one + x));
        }
    }

    //
    // Gradient (phi) functions
    //
    /**
     * Compute the momentum gradient function (Businger-Dyer, or QNSE on the stable side).
     *
     * This is the exact companion of calc_psi_m: the two are tied by
     * psi_m(zeta) = int_0^zeta (1 - phi_m(x))/x dx, so differentiating the MOST
     * wind profile u(z) = (u_star/kappa)[ln(z/z0) - psi_m(z/L)] gives
     * du/dz = u_* phi_m(zeta) / (kappa z).
     *
     * NOTE: the surface_temp* flux calculators use the Jimenez psi pair
     *       (calc_psi_m2/calc_psi_h2) rather than this one; their companions are
     *       not implemented yet, so a prescribed-surface-temperature run gets the
     *       Businger-Dyer gradient here.
     *
     * @param[in] zeta stability coordinate z/L
     */
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    calc_phi_m (amrex::Real zeta) const
    {
        if (zeta > 0) {
            return m_qnse ? qnse_phi(zeta, m_qnse_am, m_qnse_bm) : one + beta_m * zeta;
        } else {
            return one / std::sqrt(std::sqrt(one - gamma_m * zeta));
        }
    }

    /**
     * Compute the heat gradient function (Businger-Dyer, or QNSE on the stable side), the
     * exact companion of calc_psi_h: dtheta/dz = theta_* phi_h(zeta) / (kappa z).
     *
     * See the note in calc_phi_m about the Jimenez pair.
     *
     * @param[in] zeta stability coordinate z/L
     */
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    calc_phi_h (amrex::Real zeta) const
    {
        if (zeta > 0) {
            return m_qnse ? qnse_phi(zeta, m_qnse_ah, m_qnse_bh) : one + beta_h * zeta;
        } else {
            return one / std::sqrt(one - gamma_h * zeta);
        }
    }

    /**
     * Stable z/L at which the prescribed-flux problem folds.
     *
     * With L fixed by u* and the flux, the wind a given u* supports is
     * u* max(C - psi_m, b) / kappa, with b = most_bracket_floor(C). Below the smallest such
     * wind the prescribed stable flux has no solution with z/L below most_zeta_lim (Basu et
     * al., BLM 2008; with QNSE the supported wind rises again past its first minimum to a
     * local maximum at the second zero of f below, at z/L of order 1e3 to 1e4, and only then
     * tends to zero with u*; the kernels never reach that range). The physical
     * (large u*) branch has z/L below this value, so bounding z/L by it leaves every point that
     * has a solution unchanged. The result is positive for every C.
     *
     * Businger-Dyer (psi_m = -beta_m z/L): for C > 0 the minimum is at z/L = C / (2 beta_m);
     * for small or negative C (z0 near or above zref) the floor binds first and the minimum is
     * at the kink z/L = (b - C) / beta_m.
     *
     * QNSE (psi_m = -c ln(1 + b_m z/L), c = (a_m - b_m)/b_m): with u* proportional to
     * (z/L)^(-1/3), the supported wind falls with z/L while
     * f(z/L) = C + c ln(1 + b_m z/L) - 3 c b_m (z/L) / (1 + b_m z/L) is positive, and f falls
     * from C at z/L = 0 to its minimum at z/L = 2 / b_m and rises without bound after it. The
     * fold is the first zero of f beyond the kink where the floor stops binding
     * (C - psi_m = b), found by Newton's method on the bracket [kink, 2 / b_m] where f is
     * convex and decreasing (f' = c b_m (b_m z/L - 2) / (1 + b_m z/L)^2 < 0, f'' > 0 for
     * b_m z/L < 5): after the first step the iterates approach the root from the left, a step
     * that would leave the bracket is replaced by its midpoint, and the iteration stops at a
     * relative change of 1e-10 (a handful of steps against the 60 passes of a bisection to
     * the same precision; the roughness kernels call this once per u* iteration). It is the
     * kink itself when f is already negative there, and most_zeta_lim when f has no zero
     * (C > c (2 - ln 3), 10.4 with the default coefficients): the wind then falls
     * monotonically with z/L up to the bound, where the kernels hold z/L and reduce the flux
     * as at a fold. Every return value is at most most_zeta_lim.
     *
     * @param[in] C ln(zref/z0)
     */
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    calc_zeta_fold (amrex::Real C) const
    {
        if (!m_qnse) {
            return amrex::max(C / (two * beta_m), (most_bracket_floor(C) - C) / beta_m);
        }
        // A non-finite C (bad zref or z0) is returned as it is, as most_rib_zeta_bracketed does,
        // so that it shows up in the kernel's results instead of a bisection artefact
        if (!amrex::Math::isfinite(C)) { return C; }
        const amrex::Real b = m_qnse_bm;
        const amrex::Real c = (m_qnse_am - m_qnse_bm) / m_qnse_bm;
        if (c <= zero) { return most_zeta_lim; } // a_m <= b_m is refused by the reader; psi_m = 0 would have no fold
        const amrex::Real bf = most_bracket_floor(C);
        amrex::Real lo = (C < bf) ? (std::exp((bf - C) / c) - one) / b : zero;
        auto f = [&] (amrex::Real z) {
            const amrex::Real s = one + b * z;
            return C + c * std::log(s) - three * c * b * z / s;
        };
        auto df = [&] (amrex::Real z) {
            const amrex::Real s = one + b * z;
            return c * b * (b * z - two) / (s * s);
        };
        const amrex::Real zlim = most_zeta_lim; // by value: a namespace constant must not bind to min's reference in device code
        if (f(lo) <= zero) { return amrex::min(lo, zlim); }
        amrex::Real hi = two / b;
        if (lo >= hi || f(hi) >= zero) { return zlim; }
        // Safeguarded Newton on [lo, hi]: f(lo) > 0 > f(hi), f convex and decreasing, so from
        // the second step on the iterates rise monotonically to the root; a step outside the
        // bracket (the first one can be, from the midpoint) falls back to the midpoint
        amrex::Real z = myhalf * (lo + hi);
        constexpr amrex::Real rtol = amrex::Real(1.0e-10);
        for (int n = 0; n < 30; ++n) {
            const amrex::Real fz = f(z);
            if (fz > zero) { lo = z; } else if (fz < zero) { hi = z; } else { break; }
            const amrex::Real d = df(z);
            amrex::Real zn = (d < zero) ? z - fz / d : myhalf * (lo + hi);
            if (!(zn > lo && zn < hi)) { zn = myhalf * (lo + hi); }
            const bool done = std::abs(zn - z) <= rtol * zn;
            z = zn;
            if (done) { break; }
        }
        return amrex::min(z, zlim);
    }

    /**
     * Select the QNSE stable-side functions phi = (1 + a zeta)/(1 + b zeta) for momentum
     * (a_m, b_m) and heat (a_h, b_h), or return to Businger-Dyer (use_qnse false, the
     * coefficients then unused). read_most_similarity_funs checks the inputs: b_m, b_h > 0
     * and a_m > b_m, a_h > b_h.
     *
     * @param[in] use_qnse whether the stable side is QNSE
     * @param[in] am       a_m
     * @param[in] bm       b_m
     * @param[in] ah       a_h
     * @param[in] bh       b_h
     */
    AMREX_GPU_HOST_DEVICE
    void
    set_stable_functions (bool use_qnse, amrex::Real am, amrex::Real bm,
                          amrex::Real ah, amrex::Real bh)
    {
        m_qnse    = use_qnse;
        m_qnse_am = am;
        m_qnse_bm = bm;
        m_qnse_ah = ah;
        m_qnse_bh = bh;
    }

    /** Whether the stable side is QNSE. */
    AMREX_GPU_HOST_DEVICE
    bool uses_qnse () const { return m_qnse; }

private:
    /** QNSE integrated function, psi = -((a - b)/b) ln(1 + b zeta), zeta > 0. */
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    static amrex::Real
    qnse_psi (amrex::Real zeta, amrex::Real a, amrex::Real b)
    {
        return -((a - b) / b) * std::log(one + b * zeta);
    }

    /** QNSE gradient function, phi = (1 + a zeta)/(1 + b zeta), zeta > 0. */
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    static amrex::Real
    qnse_phi (amrex::Real zeta, amrex::Real a, amrex::Real b)
    {
        return (one + a * zeta) / (one + b * zeta);
    }

    amrex::Real beta_m{amrex::Real(5.0)};         ///< Constants from Dyer, BLM, 1974
    amrex::Real beta_h{amrex::Real(5.0)};         ///< https://doi.org/amrex::Real(10.1007)/BF00240838
    amrex::Real gamma_m{amrex::Real(16.0)};
    amrex::Real gamma_h{amrex::Real(16.0)};
    bool        m_qnse{false};                    ///< stable side QNSE instead of Businger-Dyer
    amrex::Real m_qnse_am{amrex::Real(2.5)};      ///< QNSE a_m (Sukoriansky et al. 2005)
    amrex::Real m_qnse_bm{amrex::Real(0.2)};      ///< QNSE b_m
    amrex::Real m_qnse_ah{amrex::Real(2.5)};      ///< QNSE a_h
    amrex::Real m_qnse_bh{amrex::Real(0.2)};      ///< QNSE b_h
};

/**
 * Stable-side similarity functions of a surface layer from its inputs:
 * <prefix>.most.stable_functions = businger_dyer (the default) | qnse, and with qnse the
 * coefficients <prefix>.most.qnse_am, qnse_bm, qnse_ah, qnse_bh (2.5, 0.2, 2.5, 0.2 by
 * default, the values of the PBL schemes' erf.qnse_*). The run stops when a coefficient is
 * given without qnse, is not finite, b_m or b_h is not positive, or a_m <= b_m or a_h <= b_h
 * (phi would not exceed 1 on the stable side: no stability correction, or one of the wrong sign).
 *
 * @param[in] pp  parameter parser of the reader (prefix erf, or erf.<face> for a side-face
 *                surface layer; the immersed wall laws read the unqualified erf.most.* keys)
 * @param[in] who the reader, for the start-up line
 */
inline similarity_funs
read_most_similarity_funs (amrex::ParmParse& pp, const std::string& who = "Surface layer")
{
    similarity_funs sfuns;
    const std::string prefix = pp.getPrefix() + ".most.";
    std::string name = "businger_dyer";
    pp.queryAdd("most.stable_functions", name);
    const char* keys[4] = {"qnse_am", "qnse_bm", "qnse_ah", "qnse_bh"};
    amrex::Real coef[4] = {amrex::Real(2.5), amrex::Real(0.2), amrex::Real(2.5), amrex::Real(0.2)};
    bool given = false;
    for (int n = 0; n < 4; ++n) {
        const std::string key = std::string("most.") + keys[n];
        given = pp.query(key.c_str(), coef[n]) || given;
    }
    if (name == "qnse") {
        for (int n = 0; n < 4; ++n) {
            if (!amrex::Math::isfinite(coef[n])) {
                amrex::Abort(prefix + keys[n] + " must be finite");
            }
        }
        for (int n = 1; n < 4; n += 2) {
            if (!(coef[n] > zero)) {
                amrex::Abort(prefix + keys[n] + " must be positive (phi = (1 + a zeta)/(1 + b zeta))");
            }
            if (coef[n-1] <= coef[n]) {
                amrex::Abort(prefix + keys[n-1] + " must exceed " + prefix + keys[n] +
                             " (phi > 1 on the stable side)");
            }
        }
        sfuns.set_stable_functions(true, coef[0], coef[1], coef[2], coef[3]);
        amrex::Print() << who << " stable-side similarity functions: qnse (a_m " << coef[0]
                       << ", b_m " << coef[1] << ", a_h " << coef[2] << ", b_h " << coef[3] << ")\n";
    } else if (name == "businger_dyer") {
        if (given) {
            amrex::Abort(prefix + "qnse_am, qnse_bm, qnse_ah and qnse_bh belong to " + prefix +
                         "stable_functions = qnse");
        }
    } else {
        amrex::Abort(prefix + "stable_functions must be businger_dyer or qnse, not " + name);
    }
    return sfuns;
}


/**
 * Log-law bracket ln(zref/z0) - psi of the MOST flux kernels, held at or above
 * most_bracket_floor(C).
 *
 * In strongly unstable air psi exceeds C = ln(zref/z0) and the bracket changes sign, which
 * turns u* negative and puts the diagnosed surface temperature on the wrong side of the air
 * for a prescribed flux. With the floor the bracket is at least 1/2, so u* stays positive and
 * finite for every C.
 *
 * @param[in] C   ln(zref/z0)
 * @param[in] psi stability correction (momentum or heat)
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
most_bracket (amrex::Real C, amrex::Real psi)
{
    return amrex::max(C - psi, most_bracket_floor(C));
}

/**
 * Residual zeta - Rib (C - psi_m)^2 / (C - psi_h) of the prescribed-temperature kernels'
 * z/L equation, with the Jimenez functions and the brackets of most_bracket.
 *
 * @param[in] sfuns similarity functions
 * @param[in] C     ln(zref/z0)
 * @param[in] Rib   bulk Richardson number
 * @param[in] zeta  z/L
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
most_rib_resid (const similarity_funs& sfuns, amrex::Real C, amrex::Real Rib, amrex::Real zeta)
{
    const amrex::Real CmPsim = most_bracket(C, sfuns.calc_psi_m2(zeta));
    const amrex::Real CmPsih = most_bracket(C, sfuns.calc_psi_h2(zeta));
    return zeta - Rib * CmPsim * CmPsim / CmPsih;
}

/**
 * Whether a root of most_rib_resid lies within delta of zeta.
 *
 * The relaxed fixed-point iteration stops on a small step, which a near-touch of the
 * residual also gives where no root is near; the residual must change sign across
 * [zeta - delta, zeta + delta] for zeta to be taken as a solution.
 *
 * @param[in] sfuns similarity functions
 * @param[in] C     ln(zref/z0)
 * @param[in] Rib   bulk Richardson number
 * @param[in] zeta  z/L the relaxed iteration stopped at
 * @param[in] delta half-width of the interval checked
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool
most_rib_zeta_is_root (const similarity_funs& sfuns, amrex::Real C, amrex::Real Rib,
                       amrex::Real zeta, amrex::Real delta)
{
    const amrex::Real rlo = most_rib_resid(sfuns, C, Rib, zeta - delta);
    const amrex::Real rhi = most_rib_resid(sfuns, C, Rib, zeta + delta);
    return (rlo > zero) != (rhi > zero) || rlo == zero || rhi == zero;
}

/**
 * Solve zeta = Rib (C - psi_m)^2 / (C - psi_h) of the prescribed-temperature kernels by
 * bracketing and bisection, with the Jimenez functions and the brackets of most_bracket.
 *
 * The relaxed fixed-point iteration stalls where zeta - Rib F(zeta) nearly touches zero
 * without crossing it, which the Jimenez stable functions produce when z0 is comparable to or
 * above zref (C below about 1): each pass then moves zeta by about the tolerance. From zeta0
 * this steps, with doubling steps, in the direction the iteration moves until the residual
 * changes sign, and bisects that interval, so it returns the root the iteration was heading
 * to, or the next one past a near-touch. The first step is limited so that two nearby roots
 * are not stepped over together. The search stays within +-most_zeta_lim, and returns that
 * bound when no root lies inside it, as the prescribed-flux kernels bound zeta.
 *
 * @param[in] sfuns similarity functions
 * @param[in] C     ln(zref/z0)
 * @param[in] Rib   bulk Richardson number
 * @param[in] zeta0 z/L to start from
 * @param[in] tol   width of the final interval
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
most_rib_zeta_bracketed (const similarity_funs& sfuns, amrex::Real C, amrex::Real Rib,
                         amrex::Real zeta0, amrex::Real tol)
{
    auto resid = [&] (amrex::Real zeta) { return most_rib_resid(sfuns, C, Rib, zeta); };

    amrex::Real a  = zeta0;
    amrex::Real ra = resid(a);
    if (ra == zero) { return a; }

    // A non-finite residual (a NaN or infinite C, Rib or zeta0) cannot be bracketed: with it the
    // step below is NaN and the search never reaches zlim. Return it, so the bad input shows up
    // in the caller's results rather than as a hang.
    if (!amrex::Math::isfinite(ra)) { return ra; }

    // Step against the residual, within +-most_zeta_lim, until it changes sign
    const amrex::Real dir  = (ra > zero) ? amrex::Real(-1.0) : amrex::Real(1.0);
    const amrex::Real zlim = dir * most_zeta_lim;
    auto step = [&] (amrex::Real z, amrex::Real dz) {
        return (dir > zero) ? amrex::min(z + dz, zlim) : amrex::max(z - dz, zlim);
    };
    amrex::Real h  = amrex::max(amrex::min(std::abs(ra),
                                           amrex::Real(0.1) * amrex::max(amrex::Real(one), std::abs(a))),
                                tol);
    amrex::Real b  = step(a, h);
    amrex::Real rb = resid(b);
    // With finite values the doubling step reaches zlim in well under 200 steps (from tol to
    // 2 most_zeta_lim); the cap only guarantees that this loop ends
    int nstep = 0;
    while ((rb > zero) == (ra > zero) && b != zlim && nstep < 200) {
        a  = b;
        ra = rb;
        h *= two;
        b  = step(a, h);
        rb = resid(b);
        ++nstep;
    }
    // No root within the bound: hold zeta there
    if ((rb > zero) == (ra > zero)) { return zlim; }

    // Bisect [a, b]
    int n = 0;
    while (std::abs(b - a) > tol && n < 60) {
        const amrex::Real m  = myhalf * (a + b);
        const amrex::Real rm = resid(m);
        if ((rm > zero) == (ra > zero)) {
            a  = m;
            ra = rm;
        } else {
            b = m;
        }
        ++n;
    }
    return myhalf * (a + b);
}

/**
 * z/L of the prescribed-flux kernels, bounded below by -most_zeta_lim and above by the fold (see
 * similarity_funs::calc_zeta_fold), where the prescribed stable flux stops having a solution.
 *
 * @param[in]  zeta      z/L from the current u* and flux
 * @param[in]  zeta_fold z/L at the fold
 * @param[out] at_fold   whether the fold bound was applied
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
most_clamp_zeta (amrex::Real zeta, amrex::Real zeta_fold, bool& at_fold)
{
    at_fold = (zeta >= zeta_fold);
    return amrex::max(amrex::min(zeta, zeta_fold), -most_zeta_lim);
}

/**
 * Factors applied to the prescribed sensible and moisture fluxes at the stable fold.
 */
struct most_flux_scale
{
    amrex::Real temp{amrex::Real(1.0)};  ///< factor on surf_temp_flux
    amrex::Real moist{amrex::Real(1.0)}; ///< factor on surf_moist_flux
};

/**
 * Part of a prescribed stable flux the wind can carry at the fold.
 *
 * The virtual heat flux is tflux_s + tflux_q, with tflux_s = surf_temp_flux (1 + epsv qv) and
 * tflux_q = qflux epsv theta. When z/L is held at the fold (see
 * similarity_funs::calc_zeta_fold), the wind is below the minimum the prescribed flux needs and
 * carries only -u*^3 thv zeta_fold / (kappa g zref). Only the prescribed components that are
 * downward drive the fold, so only they are reduced, by a common factor, until the net virtual
 * flux equals that value. An upward component, or a moisture flux diagnosed from q_surf, is
 * left unchanged; if no prescribed component is downward, nothing is reduced. The factor is
 * not taken below 0, so when a diagnosed condensation flux alone is stronger than the fold
 * value, the downward prescribed components are removed and the net flux remains stronger
 * than the fold value. Returns 1 for both factors away from the fold.
 *
 * @param[in] at_fold    whether z/L was held at the fold
 * @param[in] zeta_fold  z/L at the fold
 * @param[in] ustar      friction velocity
 * @param[in] thv        virtual potential temperature at zref
 * @param[in] tflux_s    sensible part of the virtual heat flux
 * @param[in] tflux_q    moisture part of the virtual heat flux
 * @param[in] spec_qflux whether the moisture flux is prescribed
 * @param[in] kappa      von Karman constant
 * @param[in] gravity    gravitational acceleration
 * @param[in] zref       reference height
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
most_flux_scale
most_fold_flux_scale (bool at_fold, amrex::Real zeta_fold, amrex::Real ustar,
                      amrex::Real thv, amrex::Real tflux_s, amrex::Real tflux_q,
                      bool spec_qflux, amrex::Real kappa, amrex::Real gravity,
                      amrex::Real zref)
{
    most_flux_scale scale;
    const amrex::Real tflux = tflux_s + tflux_q;
    if (!at_fold || tflux >= zero) { return scale; }
    const amrex::Real tflux_fold = -ustar * ustar * ustar * thv * zeta_fold / (kappa * gravity * zref);
    if (tflux_fold <= tflux) { return scale; }
    const bool down_s = (tflux_s < zero);
    const bool down_q = (spec_qflux && tflux_q < zero);
    const amrex::Real down = (down_s ? tflux_s : zero) + (down_q ? tflux_q : zero);
    if (down >= zero) { return scale; }
    const amrex::Real r = amrex::max((down + tflux_fold - tflux) / down, amrex::Real(0.0));
    if (down_s) { scale.temp  = r; }
    if (down_q) { scale.moist = r; }
    return scale;
}

/**
 * Minimum friction velocity [m/s] of the prescribed-flux kernels' outputs (see most_flux_outputs).
 */
inline constexpr amrex::Real most_min_ustar = amrex::Real(0.001);

/**
 * Surface state of a prescribed-flux kernel once its u* iteration has converged.
 *
 * Sets t_surf, t_star and q_star (and q_surf when the moisture flux is prescribed). At the
 * fold the prescribed fluxes are reduced by most_fold_flux_scale and the Obukhov length is set
 * to the fold value. u* is floored at most_min_ustar here, so every kernel applies the same floor.
 *
 * @param[in]    mdata      MOST data (prescribed fluxes, kappa, gravity)
 * @param[in]    spec_qflux whether the moisture flux is prescribed
 * @param[in]    at_fold    whether z/L was held at the fold on the last pass
 * @param[in]    zeta_fold  z/L at the fold
 * @param[in]    C          ln(zref/z0)
 * @param[in]    psi_h      heat stability correction of the last pass
 * @param[in]    ustar      friction velocity (floored at most_min_ustar)
 * @param[in]    zref       reference height
 * @param[in]    tm         potential temperature at zref
 * @param[in]    tvm        virtual potential temperature at zref
 * @param[in]    qvm        water vapor at zref
 * @param[in]    qflux      moisture flux of the last pass
 * @param[inout] olen       Obukhov length
 * @param[out]   t_surf     surface potential temperature
 * @param[out]   t_star     temperature scale
 * @param[inout] q_surf     surface water vapor (read when the moisture flux is diagnosed)
 * @param[out]   q_star     moisture scale
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
most_flux_outputs (const most_data& mdata, bool spec_qflux, bool at_fold,
                   amrex::Real zeta_fold, amrex::Real C, amrex::Real psi_h,
                   amrex::Real ustar, amrex::Real zref, amrex::Real tm,
                   amrex::Real tvm, amrex::Real qvm, amrex::Real qflux,
                   amrex::Real& olen, amrex::Real& t_surf, amrex::Real& t_star,
                   amrex::Real& q_surf, amrex::Real& q_star)
{
    ustar = amrex::max(ustar, amrex::Real(most_min_ustar));
    const most_flux_scale scale = most_fold_flux_scale(at_fold, zeta_fold, ustar, tvm,
                                                       mdata.surf_temp_flux * (one + epsv * qvm),
                                                       qflux * epsv * tm, spec_qflux,
                                                       mdata.kappa, mdata.gravity, zref);
    if (at_fold) { olen = zref / zeta_fold; }
    const amrex::Real CmPsih = most_bracket(C, psi_h);
    const amrex::Real hflux  = scale.temp * mdata.surf_temp_flux;
    t_surf = hflux * CmPsih / (ustar * mdata.kappa) + tm;
    t_star = -hflux / ustar;
    if (spec_qflux) {
        const amrex::Real mflux = scale.moist * mdata.surf_moist_flux;
        q_surf = mflux * CmPsih / (ustar * mdata.kappa) + qvm;
        q_star = -mflux / ustar;
    } else {
        q_star = mdata.kappa * (qvm - q_surf) / CmPsih;
    }
}

/**
 * Aitken delta-squared acceleration of the relaxed u* iteration in the prescribed-flux kernels.
 *
 * Near the stable fold the slope of the fixed-point map approaches 1 and the plain iteration
 * needs hundreds of passes. From pass `start` on, every third iterate is replaced by the
 * Aitken extrapolation of the last three (kept only within a factor of 2 of them), so points
 * that converge before `start` are unchanged. The map is not stationary when w* or a
 * u*-dependent z0 also changes from pass to pass, so the kernels do not stop on an
 * extrapolated pass: the plain pass after it, with zeta and psi from the same u*, decides
 * convergence.
 */
struct most_aitken
{
    /**
     * Return the iterate to continue from.
     *
     * @param[in] u_new iterate after pass `iter`
     * @param[in] iter  number of passes done
     */
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    update (amrex::Real u_new, int iter)
    {
        extrapolated = false;
        if (iter < start) { return u_new; }
        u0 = u1;
        u1 = u2;
        u2 = u_new;
        if (++n < 3) { return u_new; }
        const amrex::Real d = u2 - two * u1 + u0;
        const amrex::Real ua = (d != zero) ? u2 - (u2 - u1) * (u2 - u1) / d : u_new;
        const amrex::Real lo = myhalf * amrex::min(u0, amrex::min(u1, u2));
        const amrex::Real hi = two * amrex::max(u0, amrex::max(u1, u2));
        if (d == zero || !(ua > lo && ua < hi)) {
            n = 0;
            return u_new;
        }
        // The extrapolated value starts the next triple
        u2 = ua;
        n  = 1;
        extrapolated = true;
        return ua;
    }

    int start{25};
    bool extrapolated{false}; ///< whether the last update returned an extrapolated value
    int n{0};
    amrex::Real u0{zero};
    amrex::Real u1{zero};
    amrex::Real u2{zero};
};


/**
 * Empirical kinematic viscosity [m2/s] formula from Andreas (1989) CRREL Rep.
 * 89-11, valid between -173 and 277 deg C.
 *
 * @param[in] T_degK air temperature in Kelvin
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
air_viscosity (amrex::Real T_degK)
{
    amrex::Real TC = T_degK - amrex::Real(273.15);
    return amrex::Real(1.326e-5)*(one + amrex::Real(6.542e-3)*TC + amrex::Real(8.301e-6)*TC*TC - amrex::Real(4.84e-9)*TC*TC*TC);
}

/**
 * Charnock roughness with optional aerodynamically-smooth contribution
 *
 * Smooth-flow limit follows z0 = 0.11 nu / u* (Smith 1988; ECMWF; WRF)
 *
 * @param[in] ustar friction velocity
 * @param[in] cnk_a charnock constant
 * @param[in] nu    intrinsic fluid viscosity
 * @param[in] visc  flag to apply smooth flow limit
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Charnock_roughness (amrex::Real ustar,
                    amrex::Real Cnk_a,
                    amrex::Real nu,
                    bool visc)
{
    amrex::Real z0 = (Cnk_a / CONST_GRAV) * ustar * ustar;
    if (visc) {
        z0 += amrex::Real(0.11) * nu / std::max(ustar, amrex::Real(0.01));
    }
    // From Davis et al 2008 MWR
    return std::min(std::max(z0, amrex::Real(1.27e-7)), amrex::Real(2.85e-3));
}

/**
 * Modified Charnock roughness
 *
 * Smooth-flow limit is intrinsic.
 *
 * @param[in] ustar friction velocity
 * @param[in] cnk_b modified charnock constant
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Mod_Charnock_roughness (amrex::Real ustar,
                        amrex::Real Cnk_b)
{
    amrex::Real z0 = std::exp( (amrex::Real(2.7)*ustar - amrex::Real(1.8)/Cnk_b) / (ustar + amrex::Real(0.17)/Cnk_b) );
    return std::min(std::max(z0, amrex::Real(1.27e-7)), amrex::Real(1.0e-2));
}

/**
 * COARE3.0, following WRF phys/module_sf_mynn.F
 *
 * Smooth-flow limit follows z0 = 0.11 nu / u* (Smith 1988; ECMWF; WRF)
 *
 * @param[in] zref  reference height for the wind speed
 * @param[in] umm   horizontal wind speed magnitude at zref
 * @param[in] ustar friction velocity
 * @param[in] nu    intrinsic fluid viscosity
 * @param[in] visc  flag to apply smooth flow limit
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
COARE3_roughness (amrex::Real zref,
                  amrex::Real umm,
                  amrex::Real ustar,
                  amrex::Real nu,
                  bool visc)
{
    amrex::Real ws10m = umm*std::log(amrex::Real(10.)/amrex::Real(1e-4)) / std::log(zref/amrex::Real(1e-4));
    amrex::Real Czc = amrex::Real(0.011) + amrex::Real(0.007)*std::min(std::max((ws10m-amrex::Real(10.))/amrex::Real(8.), amrex::Real(0.0)), amrex::Real(1.0));
    amrex::Real z0 = (Czc / CONST_GRAV) * ustar * ustar;
    if (visc) {
        z0 += amrex::Real(0.11) * nu / std::max(ustar, amrex::Real(0.01));
    }
    // From Davis et al 2008 MWR
    return std::min(std::max(z0, amrex::Real(1.27e-7)), amrex::Real(2.85e-3));
}

/**
 * Updated Donelan 2004 as found in AHW
 *
 * Matches implementations in WRF phys/physics_mmm/sf_sfclayrev.F90 and
 * phys/module_sf_mynn.F
 *
 * Smooth flow limit is intrinsic in zn2.
 *
 * @param[in] ustar friction velocity
 * @param[in] nu    intrinsic fluid viscosity
 * @param[in] visc  flag to apply smooth flow limit
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Donelan_roughness (amrex::Real ustar,
                   amrex::Real nu,
                   bool visc)
{
    constexpr amrex::Real ozo{amrex::Real(1.59e-5)};
    amrex::Real zw  = std::min(std::pow(ustar/amrex::Real(1.06), amrex::Real(0.3)),amrex::Real(1.0));
    amrex::Real zn1 = amrex::Real(0.011)*ustar*ustar/CONST_GRAV + ozo;
    amrex::Real zn2 = amrex::Real(10.)*std::exp(-amrex::Real(9.5)*std::pow(ustar, -amrex::Real(.3333)));
    if (visc) {
        zn2 += amrex::Real(0.11) * nu / std::max(ustar, amrex::Real(0.01));
    }
    amrex::Real z0 = (amrex::Real(1.0)-zw) * zn1 + zw * zn2;
    // From Davis et al 2008 MWR
    return std::min(std::max(z0, amrex::Real(1.27e-7)), amrex::Real(2.85e-3));
}


/**
 * Wave coupled roughness
 *
 * Smooth-flow limit follows z0 = 0.11 nu / u* (Smith 1988; ECMWF; WRF)
 *
 * @param[in] Hwave wave height
 * @param[in] Lwave wave length
 * @param[in] eps   avoid divide by zero
 * @param[in] ustar friction velocity
 * @param[in] nu    intrinsic fluid viscosity
 * @param[in] visc  flag to apply smooth flow limit
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
WaveCoupled_roughness (amrex::Real Hwave,
                       amrex::Real Lwave,
                       amrex::Real eps,
                       amrex::Real ustar,
                       amrex::Real nu,
                       bool visc)
{

    amrex::Real z0 = amrex::Real(1200.0) * Hwave * std::pow( Hwave/(Lwave+eps), amrex::Real(4.5) );
    if (visc) {
        z0 += amrex::Real(0.11) * nu / std::max(ustar, amrex::Real(0.01));
    }
    return std::min(std::max(z0, amrex::Real(1.27e-7)), amrex::Real(1.0e-2));
}


/// Lower bound of the roughness recovered by Compute_roughness; shared with
/// the guard in Compute_roughness so the two cannot drift apart.
constexpr amrex::Real inverted_roughness_lo = amrex::Real(1.0e-7);

/**
 * Bounds on the roughness recovered by Compute_roughness.  The inversion is
 * unbounded above: in weak wind under strong stability (u* from an LSM or bulk
 * flux large relative to kappa*U, psi_m << 0) it returns z0 >> zref.  The
 * upper bound is above every tabulated land-cover roughness in NoahmpTable.TBL
 * (Z0MVT <= 1.10 m) and below the 2 m diagnostic height; z0 is also kept below
 * half of zref so that ln(zref/z0) stays positive.
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
bound_inverted_roughness (amrex::Real z0, amrex::Real zref) noexcept
{
    const amrex::Real z0_lo = inverted_roughness_lo;
    const amrex::Real z0_hi = amrex::min(amrex::Real(1.5), amrex::Real(0.5) * zref);
    return amrex::min(amrex::max(z0, z0_lo), z0_hi);
}

/**
 * Invert log-law to get roughness
 *
 * @param[in] zref  reference height for the wind speed
 * @param[in] Olen  Obukhov length
 * @param[in] umm   horizontal wind speed magnitude at zref
 * @param[in] ustar friction velocity
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_roughness (amrex::Real zref,
                   amrex::Real Olen,
                   amrex::Real umm,
                   amrex::Real ustar)
{
    similarity_funs sfuns;
    amrex::Real zeta  = zref / Olen;
    amrex::Real psi_m = sfuns.calc_psi_m2(zeta);
    // z0 = zref / exp(arg) is below the lower bound once arg > log(zref / z0_lo),
    // so the bound is returned without evaluating exp(arg), which overflows when
    // u* sits at its floor of machine epsilon (zero stress: bulk_coeff with
    // Cd = 0, or an LSM that returns none). At or below the threshold exp(arg)
    // is finite and the inversion runs as before, so no other input changes.
    const amrex::Real arg = KAPPA*umm/ustar + psi_m;
    if (arg > std::log(zref / inverted_roughness_lo)) {
        return bound_inverted_roughness(inverted_roughness_lo, zref);
    }
    amrex::Real z0    = zref / std::exp(arg);
    return bound_inverted_roughness(z0, zref);
}
#endif
