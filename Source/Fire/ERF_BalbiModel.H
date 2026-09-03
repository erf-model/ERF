#ifndef ERF_BALBI_MODEL_H
#define ERF_BALBI_MODEL_H

#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>
#include <cmath>
#include <vector>
#include <ERF_Rothermel.H>
#include <ERF_FireParams.H>
#include <ERF_MoistureExtinction.H>
#include <ERF_LevelSetAdvection.H>

/**
 * @file ERF_BalbiModel.H
 * @brief Balbi physical fire spread models (2009 and 2020 formulations)
 *
 * Two formulations are available, selected with
 * erf.fire.balbi.formulation = "2009" (default) or "2020".
 *
 * Balbi (2009) — steady explicit form. The rate of spread follows from the
 * radiant energy flux received by unburned fuel from a tilted flame front:
 *
 *   R       = A × (1 + sin α − cos α)                 [m/s]
 *   A       = χ σ_m δ_m / (2 τ₀ B*)                   [m/s]
 *   χ       = r₀₀ σ_m / (1 + r₀₀ σ_m)                 [−]
 *   B*      = (C_pf (T_i − T_a) + M_f Λ) / h          [−]
 *   tan α   = U / v_b + tan θ
 *   v_b     = sqrt(g δ_m (T_f − T_a) / T_a)           [m/s]
 *
 * This form has no no-wind base spread (R → 0 as wind and slope vanish) and
 * saturates at 2A for large flame tilt.
 *
 * Balbi (2020) — convective-radiative implicit form. R is the root of an
 * equation combining a radiative base term, a convective term and the flame
 * radiation term, which restores a nonzero no-wind spread and an unbounded
 * wind response:
 *
 *   R    = R_b + R_c + R_r
 *   R_b  = min(S/π, 1) σ_B T⁴ / (β ρ_v q)                       (eq. 13)
 *   R_c  = s ΔH/(q τ₀) min(δ_m, 2π/(s β))
 *          × [ δ_m/(2δ_m + H) tan α + U e^(−K₁ √β R) / u₀ ]     (eq. 27)
 *   R_r  = A R (1 + sin γ − cos γ) / (1 + R cos γ / (s r₀₀))    (eq. 15)
 *   q    = C_pf (T_i − T_a) + M_f (Λ + C_pw (T_vap − T_a))      (eq. 9)
 *   A    = min(S/2π, 1) χ₀ ΔH / (4 q)                           (eq. 16, 17)
 *   χ    = χ₀ / (1 + R cos γ / (s r₀₀))                         (eq. C7)
 *   T    = T_a + ΔH (1 − χ) / (C_pa (s_t + 1))                  (eq. B11)
 *   u₀   = 2 (s_t + 1)/τ₀ × T/T_a × ρ_v/ρ_a × min(S, 2π)        (eq. B9)
 *   γ    = atan(tan α + U/u₀)                                   (eq. 2)
 *   H    = u₀² / (g (T/T_a − 1))                                (eq. 23)
 *   β    = w / (δ_m ρ_v),  S = s β δ_m
 *
 * Inputs:
 *   U:   wind speed along the spread direction [m/s]
 *   δ_m: fuel bed depth [m]
 *   σ_m (s): surface-area-to-volume ratio of the 1-hr dead fuel [1/m]
 *   w:   1-hr dead fuel load [kg/m²]
 *   ρ_v: fuel particle density [kg/m³]
 *   h (ΔH): heat of combustion [J/kg]
 *   M_f: fuel moisture content [fraction]
 *   All US/Rothermel unit inputs (sigma_d1 [1/ft], delta [ft], w_d1 [lb/ft²],
 *   rho_p [lb/ft³], heat_content [BTU/lb]) are converted to SI internally.
 *
 * Optional couplings, all off by default:
 *   balbi.directional       ROS evaluated along the local front normal on the
 *                           level-set path instead of using wind and slope
 *                           magnitudes.
 *   balbi.use_surface_temp  ambient temperature taken per cell from
 *                           fire_surface_temp.
 *   balbi.heat_flux_coupling  buoyancy velocity augmented by the fire-induced
 *                           component built from the surface heat flux.
 *   balbi.use_cell_moisture fuel moisture taken per cell from the moisture ODE
 *                           state instead of the domain average.
 *   balbi.use_moisture_extinction  zero ROS at and above the fuel model's
 *                           moisture of extinction.
 *   balbi.herb_curing       fraction of the live herbaceous load carried as
 *                           cured dead fine fuel (2020 form only).
 *   balbi.wind_source       "midflame" (default) or "reference": which wind
 *                           the model consumes.
 *
 * Limitations:
 *   Not appropriate for fuels with very low fuel bed depth (δ_m < 0.01 m)
 *   — returns zero ROS for such cases.
 *
 * The 2009 form is a parity port of wildfire_levelset/src/balbi_model.H and
 * wildfire_levelset/src/compute_balbi_R.H
 * (https://github.com/hgopalan/wildfire_levelset).
 *
 * References:
 *   Balbi, J.-H., Rossi, J.-L., Marcelli, T. & Santoni, P.-A. (2009).
 *     A physical model for wildland fires. Combustion and Flame, 156(12),
 *     2217–2230.
 *   Balbi, J.-H., Chatelon, F.-J., Morvan, D., Rossi, J.-L., Marcelli, T. &
 *     Morandini, F. (2020). A convective-radiative propagation model for
 *     wildland fires. International Journal of Wildland Fire, 29(8), 723–738.
 *   Andrews, P.L. (1986). BEHAVE: Fire Behavior Prediction and Fuel Modeling
 *     System. USDA Forest Service GTR INT-194.
 */

/**
 * @struct BalbiComputed
 * @brief Precomputed Balbi coefficients (fuel + moisture + ambient dependent)
 *
 * Everything in the model that does not depend on the rate of spread itself,
 * on the wind or on the slope. Trivially copyable so it can live in a device
 * lookup table indexed by fuel code.
 */
struct BalbiComputed
{
    // --- Balbi (2009) ---
    amrex::Real A_coeff;       ///< ROS amplitude coefficient A [m/s]
    amrex::Real v_b;           ///< Buoyancy velocity [m/s]

    // --- Balbi (2020) ---
    amrex::Real A_rad;         ///< Radiant factor A (eq. 16) [-]
    amrex::Real Rb_coef;       ///< R_b = Rb_coef T⁴ (eq. 13) [m/(s·K⁴)]
    amrex::Real Rc_coef;       ///< Leading factor of R_c (eq. 27) [m/s]
    amrex::Real u0_coef;       ///< u₀ = u0_coef T/T_a (eq. B9) [m/s]
    amrex::Real chi_coef;      ///< χ = χ₀/(1 + chi_coef R cos γ) (eq. C7) [s/m]
    amrex::Real s_r00;         ///< s r₀₀ (eq. 15 denominator) [-]
    amrex::Real K1_sqrt_beta;  ///< K₁ √β (eq. 27) [s/m]
    amrex::Real h_depth;       ///< Fuel bed depth δ_m [m]
    amrex::Real T_a;           ///< Ambient temperature used to build these [K]
    amrex::Real T_flame_coef;  ///< ΔH/(C_pa (s_t + 1)) (eq. B11) [K]
    amrex::Real chi_0;         ///< Radiant fraction χ₀ [-]
    amrex::Real tol;           ///< Root-solve bracket width [m/s]
    int         max_iter;      ///< Root-solve iteration cap
    int         formulation;   ///< 0 = Balbi (2009), 1 = Balbi (2020)
};

// Unit conversions, US/Rothermel → SI
inline constexpr amrex::Real BALBI_FT_TO_M         = 0.3048;   ///< [ft] → [m]
inline constexpr amrex::Real BALBI_FT_INV_TO_M_INV = 3.28084;  ///< [1/ft] → [1/m]
inline constexpr amrex::Real BALBI_BTU_LB_TO_J_KG  = 2326.0;   ///< [BTU/lb] → [J/kg]
inline constexpr amrex::Real BALBI_LB_FT2_TO_KG_M2 = 4.8824;   ///< [lb/ft²] → [kg/m²]
inline constexpr amrex::Real BALBI_LB_FT3_TO_KG_M3 = 16.0185;  ///< [lb/ft³] → [kg/m³]

/// Gravitational acceleration [m/s²]
inline constexpr amrex::Real BALBI_GRAVITY = 9.81;

/// Minimum safe buoyancy velocity [m/s]. Prevents division by zero in the
/// tilt-angle kernel while staying orders of magnitude below physically
/// meaningful v_b values (≥ 1 m/s).
inline constexpr amrex::Real BALBI_V_B_MIN = 1.0e-3;

/// Minimum dimensionless ignition energy B* before a fuel is treated as
/// non-burnable (A_coeff = 0). Physical values are 0.01–0.10.
inline constexpr amrex::Real BALBI_B_STAR_MIN = 1.0e-6;

/// Minimum fuel bed depth [m] below which the model is not applicable.
inline constexpr amrex::Real BALBI_MIN_DEPTH_M = 0.01;

/// Sanity cap on the rate of spread [m/s]. Only reachable if the Balbi (2020)
/// fixed point runs away; it keeps a diverging iterate finite.
inline constexpr amrex::Real BALBI_MAX_ROS = 30.0;

/// Extreme-wind guard [m/s]: cells above this speed (or holding a non-finite
/// wind) fall back to the domain-mean wind before the ROS is evaluated.
inline constexpr amrex::Real BALBI_MAX_WIND_SPEED_MS = 25.0;
inline constexpr amrex::Real BALBI_MAX_WIND_SPEED_SQ =
    BALBI_MAX_WIND_SPEED_MS * BALBI_MAX_WIND_SPEED_MS;

/// Balbi table size: Anderson FBFM13 codes 1–13 plus index 0 (non-burnable).
inline constexpr int BALBI_TABLE_SIZE = 14;

/// |∇φ| below which the front normal is taken as undefined and the isotropic
/// wind and slope magnitudes are used instead.
inline constexpr amrex::Real BALBI_GRAD_PHI_MIN = 1.0e-12;

/**
 * @brief Compute Balbi fuel/thermal coefficients
 *
 * Precomputes everything that depends on fuel type, fuel moisture and ambient
 * temperature but not on wind, slope or the rate of spread. Callable from the
 * device so that per-cell fuel and per-cell ambient temperature can be used
 * inside a kernel.
 *
 * @param[in] fp    FuelModelParams for the selected fuel model (US units)
 * @param[in] bp    BalbiParams with thermal properties and constants
 * @param[in] M_f   Fuel moisture content [fraction]
 * @param[in] T_amb Ambient temperature [K]; values ≤ 0 select bp.T_a
 * @return BalbiComputed with precomputed coefficients
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
BalbiComputed compute_balbi_params(
    const FuelModelParams& fp,
    const FireParams::BalbiParams& bp,
    amrex::Real M_f,
    amrex::Real T_amb = -1.0) noexcept
{
    BalbiComputed bc{};
    bc.v_b          = 1.0;  // avoid division by zero in the ROS kernel
    bc.h_depth      = 0.0;
    bc.chi_0        = bp.chi_0;
    bc.tol          = bp.tol;
    bc.max_iter     = bp.max_iter;
    bc.formulation  = bp.formulation;

    const amrex::Real T_a = (T_amb > 0.0) ? T_amb : bp.T_a;
    bc.T_a = T_a;

    // Convert fuel parameters from US/Rothermel units to SI
    const amrex::Real sigma_m = fp.sigma_d1     * BALBI_FT_INV_TO_M_INV;  // [1/m]
    const amrex::Real delta_m = fp.delta        * BALBI_FT_TO_M;          // [m]
    const amrex::Real h_si    = fp.heat_content * BALBI_BTU_LB_TO_J_KG;   // [J/kg]
    // Dead fine fuel load, plus the cured fraction of the live herbaceous load,
    // which behaves as dead fine fuel once cured.
    const amrex::Real w_fine  = fp.w_d1 + bp.herb_curing * fp.w_lh;        // [lb/ft²]
    const amrex::Real w_si    = w_fine          * BALBI_LB_FT2_TO_KG_M2;  // [kg/m²]
    const amrex::Real rho_v   = fp.rho_p        * BALBI_LB_FT3_TO_KG_M3;  // [kg/m³]

    // Guard against degenerate fuel beds. The depth threshold is in metres,
    // matching the documented limit of the model.
    if (delta_m < BALBI_MIN_DEPTH_M || sigma_m <= 0.0 || h_si <= 0.0) {
        return bc;
    }

    // Moisture of extinction. Neither Balbi formulation carries an extinction
    // limit of its own, so a fuel bed wetter than its published M_x keeps
    // spreading unless this cutoff is applied.
    if (bp.use_moisture_extinction && M_f >= fp.Mx) {
        return bc;
    }

    bc.h_depth = delta_m;

    // ------------------------------------------------------------------
    // Balbi (2009)
    // ------------------------------------------------------------------
    // Radiation fraction χ = χ₀ / (1 + χ₀), with χ₀ = r₀₀ σ_m dimensionless
    const amrex::Real chi_rad_0 = bp.r_00 * sigma_m;
    const amrex::Real chi_rad   = chi_rad_0 / (1.0 + chi_rad_0);

    // Dimensionless ignition energy B* = (C_pf (T_i − T_a) + M_f Λ) / h.
    // Λ (bp.delta_H) is the latent heat of vaporisation of the fuel water;
    // h (fp.heat_content) is the heat of combustion. These are distinct.
    const amrex::Real E_ig   = bp.C_pf * (bp.T_i - T_a) + M_f * bp.delta_H;
    const amrex::Real B_star = E_ig / h_si;

    // Mean flame temperature: either the fixed input, or derived from the
    // combustion energy the same way the 2020 form does it (eq. B11 with
    // χ = χ₀).
    const amrex::Real T_f = bp.flame_temp_from_combustion
                          ? T_a + h_si * (1.0 - bp.chi_0) / (bp.C_pa * (bp.st + 1.0))
                          : bp.T_f;

    // Buoyancy velocity v_b = sqrt(g δ_m (T_f − T_a) / T_a)
    const amrex::Real dT_f = T_f - T_a;
    const amrex::Real v_b  = (dT_f > 0.0)
                           ? std::sqrt(BALBI_GRAVITY * delta_m * dT_f / T_a)
                           : BALBI_V_B_MIN;

    // ROS amplitude coefficient. σ_m δ_m is dimensionless, so dividing by
    // τ₀ [s/m] leaves [m/s].
    const amrex::Real A_coeff = (B_star > BALBI_B_STAR_MIN)
                              ? chi_rad * sigma_m * delta_m / (2.0 * bp.tau_0 * B_star)
                              : 0.0;

    bc.A_coeff = amrex::max(A_coeff, amrex::Real(0.0));
    bc.v_b     = amrex::max(v_b, BALBI_V_B_MIN);

    // ------------------------------------------------------------------
    // Balbi (2020)
    // ------------------------------------------------------------------
    if (w_si > 0.0 && rho_v > 0.0) {
        const amrex::Real beta = w_si / (delta_m * rho_v);            // packing ratio
        const amrex::Real S    = sigma_m * beta * delta_m;            // = 2 × LAI

        // Ignition energy (eq. 9). Unlike B* above this one carries the
        // sensible heating of the fuel water up to T_vap, through the specific
        // heat of water rather than that of the fuel.
        const amrex::Real q = bp.C_pf * (bp.T_i - T_a)
                            + M_f * (bp.delta_H + bp.C_pw * (bp.T_vap - T_a));

        if (beta > 0.0 && q > 0.0) {
            constexpr amrex::Real two_pi = 6.283185307179586;
            const amrex::Real a_r = amrex::min(S / two_pi, amrex::Real(1.0));   // eq. 17

            bc.A_rad   = a_r * bp.chi_0 * h_si / (4.0 * q);                     // eq. 16
            bc.Rb_coef = amrex::min(S / (0.5 * two_pi), amrex::Real(1.0))
                       * bp.sigma_B / (beta * rho_v * q);                       // eq. 13
            bc.Rc_coef = sigma_m * h_si / (q * bp.tau_0)
                       * amrex::min(delta_m, two_pi / (sigma_m * beta));        // eq. 27
            bc.u0_coef = 2.0 * (bp.st + 1.0) / bp.tau_0 * (rho_v / bp.rho_a)
                       * amrex::min(S, two_pi);                                 // eq. B9
            bc.s_r00   = sigma_m * bp.r_00;

            // eq. C7: χ = χ₀/(1 + p R τ₀ cos γ/(2 s)), and with the paper's
            // definition r₀₀ = 2/(p τ₀) that denominator is 1 + R cos γ/(s r₀₀)
            // — the same form as eq. 15. Note it is the flame tilt γ here, not
            // the slope angle; ForeFire's Balbi2020 uses cos(alpha) and groups
            // the denominator as R s cos α / r₀₀, neither of which follows the
            // paper.
            bc.chi_coef = (bc.s_r00 > 0.0) ? 1.0 / bc.s_r00 : 0.0;

            bc.K1_sqrt_beta  = bp.K1 * std::sqrt(beta);                         // eq. 27
            bc.T_flame_coef  = h_si / (bp.C_pa * (bp.st + 1.0));                // eq. B11
        }
    }

    return bc;
}

/**
 * @brief Build a per-fuel-code lookup table of Balbi coefficients
 *
 * Mirrors build_fcwh_table(): the table is indexed by integer fuel code, with
 * index 0 reserved for the non-burnable code (zero spread) and 1–13 holding
 * the Anderson FBFM13 models.
 *
 * @param[in] bp    BalbiParams with thermal properties and constants
 * @param[in] M_f   Fuel moisture content [fraction]
 * @param[in] T_amb Ambient temperature [K]; values ≤ 0 select bp.T_a
 * @return host-side vector of size BALBI_TABLE_SIZE
 */
inline std::vector<BalbiComputed> build_fuel_balbi_table(
    const FireParams::BalbiParams& bp,
    amrex::Real M_f,
    amrex::Real T_amb = -1.0)
{
    std::vector<BalbiComputed> table(BALBI_TABLE_SIZE);

    // Code 0 is non-burnable: no spread.
    table[0] = BalbiComputed{};
    table[0].v_b = 1.0;

    for (int code = 1; code < BALBI_TABLE_SIZE; ++code) {
        table[code] = compute_balbi_params(get_anderson_fuel_params(code), bp, M_f, T_amb);
    }

    return table;
}

/**
 * @brief Fire-induced vertical velocity from the surface heat flux
 *
 * v_b_Q = k_upward sqrt(g Q H_ref / (ρ_a C_pa T_a)), combined in quadrature
 * with the model's own vertical velocity scale by the callers.
 *
 * @param[in] Q   Surface heat flux [W/m²]
 * @param[in] bp  BalbiParams (supplies k_upward, hf_ref_height, ρ_a, C_pa, T_a)
 * @return induced vertical velocity [m/s], zero for non-positive flux
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real balbi_induced_velocity(amrex::Real Q,
                                   const FireParams::BalbiParams& bp) noexcept
{
    if (!(Q > 0.0)) return 0.0;
    const amrex::Real denom = bp.rho_a * bp.C_pa * bp.T_a;
    if (!(denom > 0.0)) return 0.0;
    return bp.k_upward * std::sqrt(BALBI_GRAVITY * Q * bp.hf_ref_height / denom);
}

/**
 * @brief Right-hand side of the Balbi (2020) ROS equation (eq. 28)
 *
 * Evaluates R_b + R_c + R_r for a trial rate of spread. The radiant fraction
 * (eq. C7) needs the flame tilt, which needs the flame temperature that the
 * radiant fraction sets, so γ is refined once from the slope-only seed before
 * the contributions are formed. The coupling is weak and this converges the
 * inner pair to well below the outer solve's tolerance.
 *
 * @param[in] bc        Precomputed coefficients for this cell
 * @param[in] R         Trial rate of spread [m/s]
 * @param[in] U         Wind speed along the spread direction, clamped ≥ 0 [m/s]
 * @param[in] tan_slope Terrain slope along the spread direction [-]
 * @param[in] vz_extra  Additional vertical velocity, added in quadrature [m/s]
 * @return R_b + R_c + R_r [m/s]
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real balbi2020_rhs(const BalbiComputed& bc,
                          amrex::Real R,
                          amrex::Real U,
                          amrex::Real tan_slope,
                          amrex::Real vz_extra) noexcept
{
    amrex::Real cos_gamma = 1.0 / std::sqrt(1.0 + tan_slope * tan_slope);
    amrex::Real sin_gamma = tan_slope * cos_gamma;
    amrex::Real T         = bc.T_a;
    amrex::Real u0        = BALBI_V_B_MIN;

    for (int pass = 0; pass < 2; ++pass) {
        const amrex::Real chi = bc.chi_0 / (1.0 + bc.chi_coef * R * cos_gamma);  // eq. C7
        T = bc.T_a + bc.T_flame_coef * (1.0 - chi);                              // eq. B11

        u0 = bc.u0_coef * (T / bc.T_a);                                          // eq. B9
        u0 = std::sqrt(u0 * u0 + vz_extra * vz_extra);
        u0 = amrex::max(u0, BALBI_V_B_MIN);

        const amrex::Real tan_gamma = tan_slope + U / u0;                        // eq. 2
        const amrex::Real inv_sec   = 1.0 / std::sqrt(1.0 + tan_gamma * tan_gamma);
        sin_gamma = tan_gamma * inv_sec;
        cos_gamma = inv_sec;
    }

    const amrex::Real T_ratio = T / bc.T_a;

    // Flame height (eq. 23); a flame no hotter than ambient has none.
    const amrex::Real H = (T_ratio > 1.0)
                        ? u0 * u0 / (BALBI_GRAVITY * (T_ratio - 1.0))
                        : 0.0;

    const amrex::Real T2 = T * T;
    const amrex::Real Rb = bc.Rb_coef * T2 * T2;                                 // eq. 13

    const amrex::Real Rc = bc.Rc_coef
                         * (bc.h_depth / (2.0 * bc.h_depth + H) * tan_slope
                            + U * std::exp(-bc.K1_sqrt_beta * R) / u0);          // eq. 27

    const amrex::Real Rr = bc.A_rad * R * (1.0 + sin_gamma - cos_gamma)
                         / (1.0 + R * cos_gamma / bc.s_r00);                     // eq. 15

    return Rb + Rc + Rr;
}

/**
 * @brief Balbi (2020) rate of spread for one cell
 *
 * Solves R = R_b + R_c + R_r (eq. 28). The reference implementations iterate
 * by plain substitution, which oscillates rather than converging whenever the
 * radiative coefficient A is well below one — exactly the low-fuel-load case
 * the paper plots in Fig. 3. The residual is bracketed instead: R_c decays
 * with R and R_r saturates, so f(R) = R_b + R_c + R_r − R is positive at R = 0
 * and negative at the cap, and bisection always returns a root in a fixed
 * number of steps.
 *
 * @param[in] bc        Precomputed coefficients for this cell
 * @param[in] U_n       Wind speed along the spread direction [m/s]
 * @param[in] tan_slope Terrain slope along the spread direction [-]
 * @param[in] vz_extra  Additional vertical velocity, added in quadrature [m/s]
 * @return ROS [m/s]
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real balbi2020_ros(const BalbiComputed& bc,
                          amrex::Real U_n,
                          amrex::Real tan_slope,
                          amrex::Real vz_extra) noexcept
{
    if (bc.Rb_coef <= 0.0 && bc.A_rad <= 0.0) return 0.0;

    // Only the component of the wind pushing into the unburned fuel drives
    // the convective term; a backing wind leaves the base and slope terms.
    const amrex::Real U = amrex::max(U_n, amrex::Real(0.0));

    amrex::Real lo = 0.0;
    amrex::Real hi = BALBI_MAX_ROS;

    // A non-positive residual at R = 0 means the fuel bed radiates less than it
    // takes to ignite: no spread.
    if (!(balbi2020_rhs(bc, lo, U, tan_slope, vz_extra) - lo > 0.0)) return 0.0;

    // A still-positive residual at the cap means the fixed point has run away;
    // return the cap rather than an arbitrary bisection result.
    if (balbi2020_rhs(bc, hi, U, tan_slope, vz_extra) - hi > 0.0) return hi;

    for (int it = 0; it < bc.max_iter; ++it) {
        const amrex::Real mid = 0.5 * (lo + hi);
        if (balbi2020_rhs(bc, mid, U, tan_slope, vz_extra) - mid > 0.0) {
            lo = mid;
        } else {
            hi = mid;
        }
        if (hi - lo <= bc.tol) break;
    }

    return 0.5 * (lo + hi);
}

/**
 * @brief Balbi rate of spread for one cell, in the configured formulation
 *
 * @param[in] bc        Precomputed coefficients for this cell
 * @param[in] U_n       Wind speed along the spread direction [m/s]
 * @param[in] tan_slope Terrain slope along the spread direction [-]
 * @param[in] vz_extra  Additional vertical velocity, added in quadrature [m/s]
 * @return ROS [m/s]
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real balbi_ros_cell(const BalbiComputed& bc,
                           amrex::Real U_n,
                           amrex::Real tan_slope,
                           amrex::Real vz_extra = 0.0) noexcept
{
    if (bc.formulation == 1) {
        return balbi2020_ros(bc, U_n, tan_slope, vz_extra);
    }

    // Balbi (2009): R = A (1 + sin α − cos α), tan α = U/v_b + tan θ
    amrex::Real v_b = std::sqrt(bc.v_b * bc.v_b + vz_extra * vz_extra);
    v_b = amrex::max(v_b, BALBI_V_B_MIN);

    const amrex::Real tan_alpha = U_n / v_b + tan_slope;

    // sin α = tan α / sqrt(1 + tan²α), cos α = 1 / sqrt(1 + tan²α)
    const amrex::Real inv_sec   = 1.0 / std::sqrt(1.0 + tan_alpha * tan_alpha);
    const amrex::Real sin_alpha = tan_alpha * inv_sec;
    const amrex::Real cos_alpha = inv_sec;

    return amrex::max(bc.A_coeff * (1.0 + sin_alpha - cos_alpha), amrex::Real(0.0));
}

/**
 * @brief Compute the flame tilt angle for the Balbi model
 *
 * Given wind speed U and terrain slope θ, compute the resultant angle α
 * via tan α = U/v_b + tan θ. Provided as a diagnostic helper; the ROS kernels
 * obtain sin α and cos α from tan α directly.
 *
 * @param[in] U        Wind speed [m/s]
 * @param[in] theta    Terrain slope [radians]
 * @param[in] v_b      Buoyancy velocity [m/s]
 * @return angle α [radians]
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_balbi_angle(
    amrex::Real U,
    amrex::Real theta,
    amrex::Real v_b) noexcept
{
    // Local constexpr copies so the constants are usable from device code
    constexpr amrex::Real v_b_min = BALBI_V_B_MIN;
    if (v_b < v_b_min) v_b = v_b_min;  // Avoid division by zero

    amrex::Real tan_alpha = U / v_b + std::tan(theta);
    return std::atan(tan_alpha);
}

// NOTE: macarthur_ros() is defined in ERF_MacArthurModel.H — do not duplicate here.

/**
 * @struct BalbiFieldInputs
 * @brief Optional per-cell fields consumed by the Balbi ROS kernels
 *
 * Every member is optional. With none of them set the kernels evaluate the
 * default coefficients with the domain-uniform moisture and ambient
 * temperature, which is the behaviour of the plain 2009 model.
 */
struct BalbiFieldInputs
{
    const amrex::MultiFab* fuel_model   = nullptr;  ///< Per-cell fuel code
    const BalbiComputed*   table        = nullptr;  ///< Device table by fuel code
    int                    table_size   = 0;        ///< Entries in table
    const amrex::MultiFab* heat_flux    = nullptr;  ///< Surface heat flux [W/m²]
    const amrex::MultiFab* surface_temp = nullptr;  ///< Surface temperature [K]
    const amrex::MultiFab* fuel_mc      = nullptr;  ///< Fuel moisture state [fraction]
    FuelModelParams        fp{};                    ///< Uniform fuel, for per-cell rebuilds
    amrex::Real            M_f          = 0.0;      ///< Fuel moisture [fraction]
};

/**
 * @brief Select the Balbi coefficients that apply at one cell
 *
 * With a per-cell ambient temperature the coefficients have to be rebuilt in
 * the kernel, since the ambient temperature enters the ignition energy, the
 * flame temperature and the buoyancy velocity. Otherwise the per-fuel table
 * (or the single default entry) is used.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
BalbiComputed balbi_state_at_cell(int i, int j, int k,
                                  const BalbiComputed& bc_default,
                                  const FireParams::BalbiParams& bp,
                                  const FuelModelParams& fp_uniform,
                                  amrex::Real M_f,
                                  const amrex::Array4<const amrex::Real>& fuel_model,
                                  bool use_table,
                                  const BalbiComputed* tbl,
                                  int tbl_sz,
                                  const amrex::Array4<const amrex::Real>& surf_temp,
                                  bool use_surf_temp,
                                  const amrex::Array4<const amrex::Real>& fuel_mc,
                                  bool use_cell_moisture) noexcept
{
    int fuel_code = -1;
    if (use_table) {
        const int code = static_cast<int>(fuel_model(i, j, k));
        if (code >= 0 && code < tbl_sz) fuel_code = code;
    }

    if (use_surf_temp || use_cell_moisture) {
        if (fuel_code == 0) {
            BalbiComputed bc{};
            bc.v_b = 1.0;
            return bc;  // non-burnable
        }
        const amrex::Real T_amb = use_surf_temp ? surf_temp(i, j, k) : amrex::Real(-1.0);
        amrex::Real M_cell = M_f;
        if (use_cell_moisture) {
            // Component 0 is the 1-hr dead moisture, the class that carries the
            // flaming front. Clamped the same way the domain average is.
            M_cell = amrex::max(amrex::Real(0.01),
                                amrex::min(fuel_mc(i, j, k, 0), amrex::Real(0.40)));
        }
        const FuelModelParams fp = (fuel_code > 0)
                                 ? get_anderson_fuel_params(fuel_code)
                                 : fp_uniform;
        return compute_balbi_params(fp, bp, M_cell, T_amb);
    }

    return (fuel_code >= 0) ? tbl[fuel_code] : bc_default;
}

/**
 * @brief Fill fire_ros [m/s] using the Balbi model
 *
 * Mirrors the GPU kernel structure of compute_ros_field() in ERF_Rothermel.H.
 * The ROS is isotropic here: the flame tilt is built from the wind magnitude
 * and the slope magnitude, so this is the head-fire rate of spread. The
 * FARSITE path turns it into a directional front through the ellipse; the
 * level-set path can instead use fill_balbi_ros_directional() below.
 *
 * @param[out] fire_ros       Rate of spread [m/s]
 * @param[in]  fire_wind_eff  Effective wind after WAF + terrain [m/s]
 * @param[in]  fire_slopes    Terrain slopes [dz/dx, dz/dy]
 * @param[in]  bc_default     Default Balbi coefficients (uniform fuel)
 * @param[in]  bp             Balbi parameters
 * @param[in]  in             Optional per-cell fields
 */
inline void fill_balbi_ros(
    amrex::MultiFab& fire_ros,
    const amrex::MultiFab& fire_wind_eff,
    const amrex::MultiFab& fire_slopes,
    const BalbiComputed& bc_default,
    const FireParams::BalbiParams& bp,
    const BalbiFieldInputs& in = BalbiFieldInputs{})
{
    const bool use_table = (in.fuel_model != nullptr &&
                            in.table      != nullptr &&
                            in.table_size > 0);
    const bool use_surf_temp = (bp.use_surface_temp && in.surface_temp != nullptr);
    const bool use_heat_flux = (bp.heat_flux_coupling && in.heat_flux != nullptr);
    const bool use_cell_mc   = (bp.use_cell_moisture && in.fuel_mc != nullptr);

    // Domain-mean wind, used as the fallback for cells whose wind is
    // unphysically large or non-finite.
    const amrex::Long ncells = fire_wind_eff.boxArray().numPts();
    const amrex::Real ux_mean = (ncells > 0) ? fire_wind_eff.sum(0) / amrex::Real(ncells)
                                             : amrex::Real(0.0);
    const amrex::Real uy_mean = (ncells > 0) ? fire_wind_eff.sum(1) / amrex::Real(ncells)
                                             : amrex::Real(0.0);

    const BalbiComputed*  tbl        = in.table;
    const int             tbl_sz     = in.table_size;
    const FuelModelParams fp_uniform = in.fp;
    const amrex::Real     M_f        = in.M_f;

    for (amrex::MFIter mfi(fire_ros, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto ros_arr    = fire_ros.array(mfi);
        auto wind_arr   = fire_wind_eff.const_array(mfi);
        auto slopes_arr = fire_slopes.const_array(mfi);

        amrex::Array4<const amrex::Real> fuel_model;
        if (use_table) fuel_model = in.fuel_model->const_array(mfi);

        amrex::Array4<const amrex::Real> surf_temp;
        if (use_surf_temp) surf_temp = in.surface_temp->const_array(mfi);

        amrex::Array4<const amrex::Real> heat_flux;
        if (use_heat_flux) heat_flux = in.heat_flux->const_array(mfi);

        amrex::Array4<const amrex::Real> fuel_mc;
        if (use_cell_mc) fuel_mc = in.fuel_mc->const_array(mfi);

        // Local constexpr copy so the constant is usable from device code
        constexpr amrex::Real max_wind_sq = BALBI_MAX_WIND_SPEED_SQ;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (const amrex::IntVect& iv) {
            int i = iv[0];
            int j = iv[1];
            int k = 0;

            // Wind components. The comparison is written so that a non-finite
            // wind also takes the domain-mean branch.
            amrex::Real u = wind_arr(i, j, k, 0);
            amrex::Real v = wind_arr(i, j, k, 1);
            if (!(u*u + v*v <= max_wind_sq)) {
                u = ux_mean;
                v = uy_mean;
            }
            const amrex::Real wind_speed = std::sqrt(u*u + v*v);

            const BalbiComputed bc = balbi_state_at_cell(i, j, k, bc_default, bp,
                                                         fp_uniform, M_f,
                                                         fuel_model, use_table,
                                                         tbl, tbl_sz,
                                                         surf_temp, use_surf_temp,
                                                         fuel_mc, use_cell_mc);

            const amrex::Real vz_extra = use_heat_flux
                                       ? balbi_induced_velocity(heat_flux(i, j, k), bp)
                                       : amrex::Real(0.0);

            const amrex::Real dzdx = slopes_arr(i, j, k, 0);
            const amrex::Real dzdy = slopes_arr(i, j, k, 1);
            const amrex::Real tan_slope = std::sqrt(dzdx*dzdx + dzdy*dzdy);

            ros_arr(i, j, k, 0) = balbi_ros_cell(bc, wind_speed, tan_slope, vz_extra);
        });
    }
}

/**
 * @brief Fill a direction-dependent Balbi ROS field from the front normal
 *
 * The level-set front advances along n̂ = ∇φ/|∇φ|, which points from burned
 * into unburned fuel. Evaluating the Balbi flame tilt with the wind and slope
 * components along n̂,
 *
 *   tan γ = (U · n̂)/v_b + (∇z · n̂),
 *
 * gives head, flank and backing spread from the model itself instead of an
 * imposed ellipse. Cells away from the front, where |∇φ| underflows, fall back
 * to the isotropic magnitudes.
 *
 * @param[out] R_mf        Directional ROS [m/s]
 * @param[in]  phi         Level-set field, needs one ghost cell
 * @param[in]  wind        Effective wind [m/s]
 * @param[in]  slopes      Terrain slopes [dz/dx, dz/dy]
 * @param[in]  geom        Fire-grid geometry
 * @param[in]  bc_default  Default Balbi coefficients
 * @param[in]  bp          Balbi parameters
 * @param[in]  in          Optional per-cell fields
 */
inline void fill_balbi_ros_directional(
    amrex::MultiFab& R_mf,
    const amrex::MultiFab& phi,
    const amrex::MultiFab& wind,
    const amrex::MultiFab& slopes,
    const amrex::Geometry& geom,
    const BalbiComputed& bc_default,
    const FireParams::BalbiParams& bp,
    const BalbiFieldInputs& in = BalbiFieldInputs{})
{
    const bool use_table = (in.fuel_model != nullptr &&
                            in.table      != nullptr &&
                            in.table_size > 0);
    const bool use_surf_temp = (bp.use_surface_temp && in.surface_temp != nullptr);
    const bool use_heat_flux = (bp.heat_flux_coupling && in.heat_flux != nullptr);
    const bool use_cell_mc   = (bp.use_cell_moisture && in.fuel_mc != nullptr);

    const amrex::Long ncells = wind.boxArray().numPts();
    const amrex::Real ux_mean = (ncells > 0) ? wind.sum(0) / amrex::Real(ncells)
                                             : amrex::Real(0.0);
    const amrex::Real uy_mean = (ncells > 0) ? wind.sum(1) / amrex::Real(ncells)
                                             : amrex::Real(0.0);

    const amrex::Real dx = geom.CellSize(0);
    const amrex::Real dy = geom.CellSize(1);

    const BalbiComputed*  tbl        = in.table;
    const int             tbl_sz     = in.table_size;
    const FuelModelParams fp_uniform = in.fp;
    const amrex::Real     M_f        = in.M_f;

    for (amrex::MFIter mfi(R_mf); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto R_arr      = R_mf.array(mfi);
        auto phi_arr    = phi.const_array(mfi);
        auto wind_arr   = wind.const_array(mfi);
        auto slopes_arr = slopes.const_array(mfi);

        amrex::Array4<const amrex::Real> fuel_model;
        if (use_table) fuel_model = in.fuel_model->const_array(mfi);

        amrex::Array4<const amrex::Real> surf_temp;
        if (use_surf_temp) surf_temp = in.surface_temp->const_array(mfi);

        amrex::Array4<const amrex::Real> heat_flux;
        if (use_heat_flux) heat_flux = in.heat_flux->const_array(mfi);

        amrex::Array4<const amrex::Real> fuel_mc;
        if (use_cell_mc) fuel_mc = in.fuel_mc->const_array(mfi);

        constexpr amrex::Real max_wind_sq = BALBI_MAX_WIND_SPEED_SQ;
        constexpr amrex::Real grad_min    = BALBI_GRAD_PHI_MIN;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (const amrex::IntVect& iv) {
            const int i = iv[0];
            const int j = iv[1];
            const int k = iv[2];

            amrex::Real u = wind_arr(i, j, k, 0);
            amrex::Real v = wind_arr(i, j, k, 1);
            if (!(u*u + v*v <= max_wind_sq)) {
                u = ux_mean;
                v = uy_mean;
            }

            const amrex::Real dzdx = slopes_arr(i, j, k, 0);
            const amrex::Real dzdy = slopes_arr(i, j, k, 1);

            // Front normal from central differences; it points toward
            // increasing phi, i.e. into the unburned fuel.
            const amrex::Real dpdx = (phi_arr(i+1, j, k) - phi_arr(i-1, j, k)) / (2.0 * dx);
            const amrex::Real dpdy = (phi_arr(i, j+1, k) - phi_arr(i, j-1, k)) / (2.0 * dy);
            const amrex::Real grad = std::sqrt(dpdx*dpdx + dpdy*dpdy);

            amrex::Real U_n, tan_slope;
            if (grad > grad_min) {
                const amrex::Real nx = dpdx / grad;
                const amrex::Real ny = dpdy / grad;
                U_n       = u*nx + v*ny;
                tan_slope = dzdx*nx + dzdy*ny;
            } else {
                U_n       = std::sqrt(u*u + v*v);
                tan_slope = std::sqrt(dzdx*dzdx + dzdy*dzdy);
            }

            const BalbiComputed bc = balbi_state_at_cell(i, j, k, bc_default, bp,
                                                         fp_uniform, M_f,
                                                         fuel_model, use_table,
                                                         tbl, tbl_sz,
                                                         surf_temp, use_surf_temp,
                                                         fuel_mc, use_cell_mc);

            const amrex::Real vz_extra = use_heat_flux
                                       ? balbi_induced_velocity(heat_flux(i, j, k), bp)
                                       : amrex::Real(0.0);

            R_arr(i, j, k) = balbi_ros_cell(bc, U_n, tan_slope, vz_extra);
        });
    }
}

/**
 * @brief Advect the level set with a direction-dependent Balbi ROS
 *
 * Same WENO5-Z / SSP-RK3 scheme as advect_levelset_weno5z_rk3(), except that
 * the ROS field is rebuilt from the front normal at the start of every RK
 * stage, so R follows the front as it turns.
 *
 * @param[in,out] phi        Level-set field
 * @param[in]     wind       Effective wind [m/s]
 * @param[in]     slopes     Terrain slopes [dz/dx, dz/dy]
 * @param[in]     geom       Fire-grid geometry
 * @param[in]     dt         Timestep [s]
 * @param[in]     eps_visc   Artificial viscosity coefficient
 * @param[in]     bc_default Default Balbi coefficients
 * @param[in]     bp         Balbi parameters
 * @param[in]     in         Optional per-cell fields
 */
inline void advect_levelset_balbi_rk3(amrex::MultiFab& phi,
                                      const amrex::MultiFab& wind,
                                      const amrex::MultiFab& slopes,
                                      const amrex::Geometry& geom,
                                      amrex::Real dt,
                                      amrex::Real eps_visc,
                                      const BalbiComputed& bc_default,
                                      const FireParams::BalbiParams& bp,
                                      const BalbiFieldInputs& in = BalbiFieldInputs{})
{
    using namespace amrex;

    const auto* dx = geom.CellSize();

    const int nghost = 3;
    MultiFab phi_n(phi.boxArray(), phi.DistributionMap(), 1, nghost);
    MultiFab phi_star(phi.boxArray(), phi.DistributionMap(), 1, nghost);
    MultiFab phi_hat(phi.boxArray(), phi.DistributionMap(), 1, nghost);
    MultiFab rhs(phi.boxArray(), phi.DistributionMap(), 1, 0);
    MultiFab R_dir(phi.boxArray(), phi.DistributionMap(), 1, 0);

    phi_n.ParallelCopy(phi);
    phi_n.FillBoundary(geom.periodicity());

    // Stage 1: φ* = φ^n + Δt L(φ^n)
    phi_star.ParallelCopy(phi_n);
    phi_star.FillBoundary(geom.periodicity());
    fill_balbi_ros_directional(R_dir, phi_star, wind, slopes, geom, bc_default, bp, in);
    fire_levelset::compute_levelset_rhs(rhs, phi_star, R_dir, dx[0], dx[1], eps_visc, &slopes);
    for (MFIter mfi(phi_star); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto ps = phi_star.array(mfi);
        auto r  = rhs.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) noexcept {
            ps(iv) = ps(iv) + dt * r(iv);
        });
    }
    phi_star.FillBoundary(geom.periodicity());

    // Stage 2: φ** = (3/4)φ^n + (1/4)(φ* + Δt L(φ*))
    fill_balbi_ros_directional(R_dir, phi_star, wind, slopes, geom, bc_default, bp, in);
    fire_levelset::compute_levelset_rhs(rhs, phi_star, R_dir, dx[0], dx[1], eps_visc, &slopes);
    phi_hat.ParallelCopy(phi_n);
    for (MFIter mfi(phi_hat); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto ph = phi_hat.array(mfi);
        auto ps = phi_star.const_array(mfi);
        auto r  = rhs.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) noexcept {
            ph(iv) = 0.75_rt * ph(iv) + 0.25_rt * (ps(iv) + dt * r(iv));
        });
    }
    phi_hat.FillBoundary(geom.periodicity());

    // Stage 3: φ^{n+1} = (1/3)φ^n + (2/3)(φ** + Δt L(φ**))
    fill_balbi_ros_directional(R_dir, phi_hat, wind, slopes, geom, bc_default, bp, in);
    fire_levelset::compute_levelset_rhs(rhs, phi_hat, R_dir, dx[0], dx[1], eps_visc, &slopes);
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto pnew = phi.array(mfi);
        auto ph   = phi_hat.const_array(mfi);
        auto r    = rhs.const_array(mfi);
        auto pn   = phi_n.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) noexcept {
            pnew(iv) = (1.0_rt/3.0_rt) * pn(iv) + (2.0_rt/3.0_rt) * (ph(iv) + dt * r(iv));
        });
    }
}

#endif
