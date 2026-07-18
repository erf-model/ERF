#ifndef ERF_BEHAVE_MODEL_H
#define ERF_BEHAVE_MODEL_H

/**
 * @file ERF_BehaveModel.H
 * @brief BEHAVE multi-class Rothermel rate-of-spread model (Phase 15)
 *
 * Implements the full multi-class Rothermel fire spread model following
 * Andrews (2018). This extends the single-class model in ERF_Rothermel.H
 * to account for six fuel size classes:
 *   1. Dead 1-hour
 *   2. Dead 10-hour
 *   3. Dead 100-hour
 *   4. Live herbaceous
 *   5. Live woody
 *   6. Dead herbaceous (transferred from live when moisture < 30%)
 *
 * Key features:
 * - Dynamic load transfer from live herbaceous to dead category based on moisture
 * - Per-class reaction intensity with separate dead and live damping coefficients
 * - Live moisture of extinction computed from dead/live fuel ratio
 * - Separate heat sink for dead and live fuels
 *
 * This model is selected by ros_model = "behave" in ParmParse configuration.
 * Selected references:
 * - Andrews, P.L. (2018). The Rothermel surface fire spread model and associated
 *   developments: A comprehensive explanation. USDA Forest Service GTR-RMRS-371.
 * - WRF-SFIRE openwfm/WRF-SFIRE/phys/module_fr_sfire_phys.F, subroutine fire_ros_behave().
 * - Anderson, H.E. (1982). Aids to determining fuel models for estimating fire
 *   behavior. USDA Forest Service INT-122.
 */

#include <AMReX_MultiFab.H>
#include <ERF_Rothermel.H>
#include <cmath>

/**
 * @struct BehaveState
 * @brief Pre-computed BEHAVE multi-class Rothermel state for one cell.
 *
 * Contains all pre-computed quantities that depend on fuel model and
 * per-cell moisture. Computed once per timestep on CPU, used in GPU ROS kernel.
 */
struct BehaveState
{
    // Net fuel loads [lb/ft²] after mineral content correction
    amrex::Real wn_dead;         ///< Net dead fuel load
    amrex::Real wn_live;         ///< Net live fuel load

    // Heat content factors [BTU/lb × dimensionless]
    amrex::Real fuelheat_dead;   ///< Weighted dead heat content factor
    amrex::Real fuelheat_live;   ///< Weighted live heat content factor

    // Damping coefficients (Rothermel equations)
    amrex::Real etam_dead;       ///< Dead moisture damping (Rothermel Eq. 29)
    amrex::Real etam_live;       ///< Live moisture damping (Rothermel Eq. 29 with live M_x)
    amrex::Real etas_dead;       ///< Dead mineral damping (Rothermel Eq. 30)
    amrex::Real etas_live;       ///< Live mineral damping (Rothermel Eq. 30)

    // Heat sink [BTU/lb × lb/ft³]
    amrex::Real heat_sink;       ///< Combined dead+live heat sink

    // Reaction velocity and propagating flux ratio
    amrex::Real gamma;           ///< Optimum reaction velocity Γ [1/min]
    amrex::Real xifr;            ///< Propagating flux ratio ξ

    // No-wind ROS [m/s]
    amrex::Real r_0;

    // Wind/slope coefficients (identical to single-class path)
    amrex::Real C;               ///< Wind factor coefficient
    amrex::Real B;               ///< Wind factor exponent
    amrex::Real beta_ratio_E;    ///< (β/β_opt)^E
    amrex::Real phi_s_const;     ///< Slope factor constant
    amrex::Real U_max_ftmin;     ///< MEWS cap [ft/min]
    amrex::Real wind_conv;       ///< m/s → ft/min conversion factor
};

/**
 * @brief Compute BEHAVE multi-class Rothermel state for a single cell.
 *
 * Ports fire_ros_behave() from WRF-SFIRE module_fr_sfire_phys.F, implementing:
 *   - Dynamic load transfer from live herbaceous to dead at low moisture
 *   - Area-weighted SAV for dead and live fuel classes
 *   - Live moisture of extinction from fuel load ratio
 *   - Per-class moisture and mineral damping coefficients
 *   - Multi-class reaction intensity
 *   - Multi-class heat sink
 *
 * @param[in] fp         FuelModelParams — all fuel properties
 * @param[in] mc_1hr     1-hour dead fuel moisture [fraction]
 * @param[in] mc_10hr    10-hour dead fuel moisture [fraction]
 * @param[in] mc_100hr   100-hour dead fuel moisture [fraction]
 * @param[in] mc_herba   Live herbaceous fuel moisture [fraction]
 * @param[in] mc_woody   Live woody fuel moisture [fraction]
 *
 * @return BehaveState with all precomputed coefficients
 *
 * SAV (Surface Area to Volume) values:
 *   Dead 1-hr:      fp.sigma_d1       [1/ft]  (from FuelModelParams)
 *   Dead 10-hr:     109.0             [1/ft]  (fixed standard value)
 *   Dead 100-hr:    30.0              [1/ft]  (fixed standard value)
 *   Live herba:     fp.sigma_lh       [1/ft]  (from FuelModelParams)
 *   Live woody:     fp.sigma_lw       [1/ft]  (from FuelModelParams)
 *   Dead herba (transferred): fp.sigma_lh [1/ft] (same as live herba)
 *
 * Reference: Andrews (2018) GTR-RMRS-371 pp. 16–38.
 */
inline BehaveState compute_behave_state(
    const FuelModelParams& fp,
    amrex::Real mc_1hr,
    amrex::Real mc_10hr,
    amrex::Real mc_100hr,
    amrex::Real mc_herba,
    amrex::Real mc_woody)
{
    BehaveState bs;

    using amrex::Real;

    // Standard constants (Rothermel 1972, Andrews 2018)
    const Real S_T = 0.0555;   // Total mineral content
    const Real S_e = 0.010;    // Effective mineral content
    const Real FT_MIN_TO_M_S = 0.00508;

    // Fixed SAV values for standard fuel classes [1/ft]
    const Real sigma_d10 = 109.0;
    const Real sigma_d100 = 30.0;

    // =====================================================================
    // Step 1: Dynamic load transfer (live herbaceous to dead)
    // =====================================================================
    // When fmc_herba < 0.30: entire live herba load transfers to dead
    // When 0.30 <= fmc_herba < 1.20: partial transfer
    // Reference: WRF-SFIRE fire_ros_behave, lines ~384–396

    Real fuelload_herba_dead = 0.0;
    Real fuelload_herba_live = fp.w_lh;  // [lb/ft²]

    if (mc_herba < 0.30) {
        fuelload_herba_dead = fuelload_herba_live;
        fuelload_herba_live = 0.0;
    } else if (mc_herba < 1.20) {
        // Linear interpolation: transfer decreases from 100% to 0%
        // as moisture goes from 0.30 to 1.20
        Real transferred = fuelload_herba_live * (1.333 - 1.11 * mc_herba);
        fuelload_herba_dead = transferred;
        fuelload_herba_live -= transferred;
    }

    // =====================================================================
    // Step 2: Area-weighted SAV for dead and live
    // =====================================================================
    // A_k = (sigma_k / rho_p) * fuelload_k  (surface area per cell per class)
    // f_k = A_k / A_total                    (weighting factor within category)
    // savr = f_dead * savr_dead + f_live * savr_live

    Real A_d1 = (fp.sigma_d1 / fp.rho_p) * fp.w_d1;
    Real A_d10 = (sigma_d10 / fp.rho_p) * fp.w_d10;
    Real A_d100 = (sigma_d100 / fp.rho_p) * fp.w_d100;
    Real A_dh = (fp.sigma_lh / fp.rho_p) * fuelload_herba_dead;
    Real A_d_total = A_d1 + A_d10 + A_d100 + A_dh;

    Real A_lh = (fp.sigma_lh / fp.rho_p) * fuelload_herba_live;
    Real A_lw = (fp.sigma_lw / fp.rho_p) * fp.w_lw;
    Real A_l_total = A_lh + A_lw;

    // Weighting factors
    Real f_d1 = (A_d_total > 1.0e-12) ? A_d1 / A_d_total : 0.0;
    Real f_d10 = (A_d_total > 1.0e-12) ? A_d10 / A_d_total : 0.0;
    Real f_d100 = (A_d_total > 1.0e-12) ? A_d100 / A_d_total : 0.0;
    Real f_dh = (A_d_total > 1.0e-12) ? A_dh / A_d_total : 0.0;
    Real f_lh = (A_l_total > 1.0e-12) ? A_lh / A_l_total : 0.0;
    Real f_lw = (A_l_total > 1.0e-12) ? A_lw / A_l_total : 0.0;

    // Area-weighted SAV
    Real sigma_dead = f_d1 * fp.sigma_d1 + f_d10 * sigma_d10 + f_d100 * sigma_d100 + f_dh * fp.sigma_lh;
    Real sigma_live = f_lh * fp.sigma_lh + f_lw * fp.sigma_lw;

    // Total fuels
    Real w_dead = fp.w_d1 + fp.w_d10 + fp.w_d100 + fuelload_herba_dead;
    Real w_live = fuelload_herba_live + fp.w_lw;
    Real w_0 = w_dead + w_live;

    Real f_dead = (w_0 > 1.0e-12) ? w_dead / w_0 : 0.0;
    Real f_live = (w_0 > 1.0e-12) ? w_live / w_0 : 0.0;

    // =====================================================================
    // Step 3: Live moisture of extinction
    // =====================================================================
    // Computed from dead/live fuel load ratio following Rothermel (1972) Eq. 41
    // Reference: WRF-SFIRE lines ~421–430

    Real fuelmce = fp.Mx;  // Dead fuel moisture of extinction (from FuelModelParams)

    // Effective heating number weighting: exp(-138/sigma)
    Real exp_dead = w_dead * std::exp(-138.0 / amrex::max(sigma_dead, 1.0));
    Real exp_live = w_live * std::exp(-138.0 / amrex::max(sigma_live, 1.0));

    Real W_prima = (exp_live > 1.0e-12) ? (exp_dead / exp_live) : 0.0;

    // Weighted dead fuel moisture
    Real fmc_dead = f_d1 * mc_1hr + f_d10 * mc_10hr + f_d100 * mc_100hr + f_dh * 0.15;
    Real M_prima_f = fmc_dead;

    // Live moisture of extinction from Rothermel live fuel path
    Real fuelmce_live = 2.9 * W_prima * (1.0 - M_prima_f / fuelmce) - 0.226;
    fuelmce_live = amrex::max(fuelmce_live, fuelmce);  // Ensure ≥ dead fuel Mx

    // =====================================================================
    // Step 4: Per-class moisture damping
    // =====================================================================
    // Rothermel Eq. 29: η_M = max(0, 1 - 2.59*r + 5.11*r² - 3.52*r³)
    //   where r = min(1, M / M_x)
    // Reference: WRF-SFIRE lines ~431–445

    fmc_dead = f_d1 * mc_1hr + f_d10 * mc_10hr + f_d100 * mc_100hr + f_dh * 0.15;
    Real fmc_live = f_lh * mc_herba + f_lw * mc_woody;

    Real r_dead = amrex::min(1.0, fmc_dead / fuelmce);
    Real r_live = amrex::min(1.0, fmc_live / fuelmce_live);

    // Moisture damping coefficients
    auto etam_poly = [](Real r) -> Real {
        return amrex::max(0.0, 1.0 - 2.59*r + 5.11*r*r - 3.52*r*r*r);
    };
    bs.etam_dead = etam_poly(r_dead);
    bs.etam_live = etam_poly(r_live);

    // Mineral damping (Rothermel Eq. 30)
    bs.etas_dead = 0.174 * std::pow(S_e, -0.19);
    bs.etas_live = 0.174 * std::pow(S_e, -0.19);

    // =====================================================================
    // Step 5: Net fuel loads and heat content
    // =====================================================================

    bs.wn_dead = w_dead * (1.0 - S_T);
    bs.wn_live = w_live * (1.0 - S_T);

    bs.fuelheat_dead = fp.heat_content;  // Dead heat content per Anderson
    bs.fuelheat_live = fp.heat_content;  // Live heat content assumed same

    // =====================================================================
    // Step 6: Bulk density and packing ratio
    // =====================================================================

    Real rho_b = w_0 / fp.delta;  // Bulk density [lb/ft³]
    Real beta = rho_b / fp.rho_p;  // Packing ratio

    // =====================================================================
    // Step 7: Reaction velocity parameters
    // =====================================================================
    // Use dead fuel SAV as characteristic value (weighted across dead classes)
    Real sigma = sigma_dead;
    sigma = amrex::max(sigma, 100.0);  // Minimum SAV guard

    Real beta_op = 3.348 * std::pow(sigma, -0.8189);
    Real sigma_1p5 = std::pow(sigma, 1.5);
    Real Gamma_max = sigma_1p5 / (495.0 + 0.0594 * sigma_1p5);
    Real A = 133.0 * std::pow(sigma, -0.7913);
    Real beta_ratio = beta / beta_op;
    Real Gamma_prime = Gamma_max * std::pow(beta_ratio, A) * std::exp(A * (1.0 - beta_ratio));

    bs.gamma = Gamma_prime;

    // =====================================================================
    // Step 8: Propagating flux ratio and heat of preignition
    // =====================================================================
    // Rothermel Eq. 42 and Eq. 12

    bs.xifr = std::exp((0.792 + 0.681 * std::sqrt(sigma)) * (beta + 0.1))
              / (192.0 + 0.2595 * sigma);

    // Weighted moisture for heat of preignition
    Real M_weighted = f_dead * fmc_dead + f_live * fmc_live;
    Real Q_ig = 250.0 + 1116.0 * M_weighted;

    // =====================================================================
    // Step 9: Heat sink
    // =====================================================================
    // Combined dead + live heat sink per Rothermel
    // heat_sink = ρ_b * ε_h * (Σ(ε_k * q_ig_k))
    // Simplified: use constant ε_h and q_ig

    Real eps_h = std::exp(-138.0 / sigma);
    bs.heat_sink = rho_b * eps_h * Q_ig;

    // =====================================================================
    // Step 10: Reaction intensity and no-wind ROS
    // =====================================================================
    // I_R = Γ * (w_n,dead * heat_dead * η_M,dead * η_s,dead +
    //            w_n,live * heat_live * η_M,live * η_s,live)
    // R_0 = (I_R * ξ) / heat_sink

    Real I_R = bs.gamma * (bs.wn_dead * bs.fuelheat_dead * bs.etam_dead * bs.etas_dead
                          + bs.wn_live * bs.fuelheat_live * bs.etam_live * bs.etas_live);
    I_R = amrex::max(I_R, 0.01);

    Real R0_ft_min = (I_R * bs.xifr) / amrex::max(bs.heat_sink, 1.0);

    // =====================================================================
    // Step 11: Wind and slope factor coefficients
    // =====================================================================
    // Identical to single-class Rothermel path (Eqs. 47–51)

    bs.r_0 = R0_ft_min * FT_MIN_TO_M_S;  // Convert ft/min → m/s

    bs.C = 7.47 * std::exp(-0.133 * std::pow(sigma, 0.55));
    bs.B = 0.02526 * std::pow(sigma, 0.54);
    Real E = 0.715 * std::exp(-3.59e-4 * sigma);

    bs.beta_ratio_E = std::pow(beta_ratio, -E);
    bs.phi_s_const = 5.275 * std::pow(beta, -0.3);

    // MEWS wind speed cap (Andrews 2018)
    bs.U_max_ftmin = (sigma > 1000.0) ? 300.0 : 500.0;
    bs.wind_conv = 196.85;  // m/s → ft/min

    return bs;
}

/**
 * @brief Compute BEHAVE ROS for a single cell (GPU device function).
 *
 * Applies wind and slope factors to the precomputed no-wind ROS.
 * The wind and slope factor formulae are identical to the single-class
 * Rothermel path — only the precomputed coefficients and R_0 differ.
 *
 * @param[in] ux_eff  Effective wind u-component [m/s]
 * @param[in] uy_eff  Effective wind v-component [m/s]
 * @param[in] sx      Slope dz/dx [dimensionless]
 * @param[in] sy      Slope dz/dy [dimensionless]
 * @param[in] bs      BehaveState with precomputed coefficients
 *
 * @return Rate of spread [m/s]
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real behave_ros_cell(
    amrex::Real ux_eff, amrex::Real uy_eff,
    amrex::Real sx, amrex::Real sy,
    const BehaveState& bs) noexcept
{
    using amrex::Real;

    // Effective wind speed [m/s]
    Real U_eff = std::sqrt(ux_eff*ux_eff + uy_eff*uy_eff);

    // Wind factor: φ_w = C * U_eff^B * (β/β_opt)^E
    // Convert to ft/min for Rothermel calculation
    Real U_eff_ftmin = U_eff * bs.wind_conv;

    // Apply MEWS cap
    U_eff_ftmin = amrex::min(U_eff_ftmin, bs.U_max_ftmin);

    // Slope magnitude and squared tangent
    Real tan_slope_sq = sx*sx + sy*sy;

    // Cap at tan²(80°) ≈ 32.16
    tan_slope_sq = amrex::min(tan_slope_sq, 32.0);

    // Slope factor (Rothermel Eq. 51)
    Real phi_s = bs.phi_s_const * tan_slope_sq;

    // Wind factor
    Real phi_w = bs.C * std::pow(U_eff_ftmin, bs.B) * bs.beta_ratio_E;

    // Rate of spread: R = R0 * (1 + φ_w + φ_s)
    Real ROS = bs.r_0 * (1.0 + phi_w + phi_s);

    return amrex::max(ROS, 0.0);
}

/**
 * @brief Fill fire_ros MultiFab using BEHAVE multi-class Rothermel model.
 *
 * @param[out] fire_ros        ROS field [m/s] on fire grid
 * @param[in]  fire_wind_eff   Effective wind (u,v) [m/s]
 * @param[in]  fire_slopes     Terrain slopes (dz/dx, dz/dy)
 * @param[in]  bs_default      Default BehaveState for uniform moisture
 * @param[in]  fire_fuel_mc    Per-cell fuel moisture (5 components, optional)
 * @param[in]  moisture_dynamic Whether per-cell moisture should be used
 */
inline void fill_behave_ros(
    amrex::MultiFab& fire_ros,
    const amrex::MultiFab& fire_wind_eff,
    const amrex::MultiFab& fire_slopes,
    const FuelModelParams& fp,
    const BehaveState& bs_default,
    const amrex::MultiFab* fire_fuel_mc = nullptr,
    bool moisture_dynamic = false)
{
    // When moisture_dynamic and fire_fuel_mc is provided, compute per-cell
    // BehaveState using per-cell moisture. Otherwise use bs_default for all cells.

    if (!moisture_dynamic || fire_fuel_mc == nullptr) {
        // Uniform moisture path: use bs_default for all cells
        for (amrex::MFIter mfi(fire_ros, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();
            auto ros = fire_ros.array(mfi);
            auto wind = fire_wind_eff.const_array(mfi);
            auto slopes = fire_slopes.const_array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (const amrex::IntVect& iv) {
                int i = iv[0], j = iv[1];
                ros(i, j, 0) = behave_ros_cell(wind(i,j,0,0), wind(i,j,0,1),
                                               slopes(i,j,0,0), slopes(i,j,0,1),
                                               bs_default);
            });
        }
    } else {
        // Per-cell moisture path: compute BehaveState for each cell on CPU,
        // store in temporary array, then run GPU kernel

        // Allocate temporary per-cell BehaveState storage
        // Pack into flat array: nRealPerBehaveState = sizeof(BehaveState)/sizeof(Real)
        constexpr int nRealPerBehaveState = sizeof(BehaveState) / sizeof(amrex::Real);
        auto bs_temp = std::make_unique<amrex::MultiFab>(
            fire_ros.boxArray(), fire_ros.DistributionMap(),
            nRealPerBehaveState, 0);

        // CPU pass: Compute per-cell BehaveState
        for (amrex::MFIter mfi(*fire_fuel_mc, false); mfi.isValid(); ++mfi) {
            auto mc = fire_fuel_mc->const_array(mfi);
            auto bs_arr = bs_temp->array(mfi);

            // LoopOnCpu to compute BehaveState (not GPU-compatible)
            amrex::LoopOnCpu(mfi.tilebox(), [&](int i, int j, int k) {
                BehaveState bs = compute_behave_state(fp,
                    mc(i,j,k,0), mc(i,j,k,1), mc(i,j,k,2),
                    mc(i,j,k,3), mc(i,j,k,4));

                // Pack BehaveState into flat Real array
                BehaveState* bs_ptr = reinterpret_cast<BehaveState*>(&bs_arr(i,j,k,0));
                *bs_ptr = bs;
            });
        }

        // GPU pass: Compute ROS using per-cell BehaveState
        for (amrex::MFIter mfi(fire_ros, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();
            auto ros = fire_ros.array(mfi);
            auto wind = fire_wind_eff.const_array(mfi);
            auto slopes = fire_slopes.const_array(mfi);
            auto bs_arr = bs_temp->const_array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (const amrex::IntVect& iv) {
                int i = iv[0], j = iv[1];

                // Unpack BehaveState from flat Real array
                const BehaveState* bs_ptr = reinterpret_cast<const BehaveState*>(&bs_arr(i,j,0,0));
                const BehaveState& bs = *bs_ptr;

                ros(i, j, 0) = behave_ros_cell(wind(i,j,0,0), wind(i,j,0,1),
                                               slopes(i,j,0,0), slopes(i,j,0,1),
                                               bs);
            });
        }
    }
}

#endif  // ERF_BEHAVE_MODEL_H
